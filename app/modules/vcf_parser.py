"""
VCF v4.2 Parser Module
Parses VCF files and extracts pharmacogenomic variants (GENE, STAR, RS INFO tags).
"""
from __future__ import annotations

import re
import io
import logging
from dataclasses import dataclass, field
from typing import Optional

from app.models import ParsedVariant

logger = logging.getLogger(__name__)

# The 6 pharmacogenomically critical genes we care about
PGX_GENES = {"CYP2D6", "CYP2C19", "CYP2C9", "TPMT", "DPYD", "UGT1A1"}

# Known PGx rsIDs for fallback gene identification
RS_TO_GENE: dict[str, str] = {
    "rs3892097": "CYP2D6",  "rs1065852": "CYP2D6", "rs5030655": "CYP2D6",
    "rs16947":   "CYP2D6",  "rs28371706": "CYP2D6","rs28371725": "CYP2D6",
    "rs4244285": "CYP2C19", "rs4986893": "CYP2C19","rs28399504": "CYP2C19",
    "rs56337013": "CYP2C19","rs72552267": "CYP2C19","rs12248560": "CYP2C19",
    "rs1799853": "CYP2C9",  "rs1057910": "CYP2C9", "rs28371686": "CYP2C9",
    "rs7900194": "CYP2C9",  "rs2256871": "CYP2C9",
    "rs1800460": "TPMT",    "rs1142345": "TPMT",   "rs1800462": "TPMT",
    "rs1051334": "TPMT",
    "rs3918290": "DPYD",    "rs55886062": "DPYD",  "rs67376798": "DPYD",
    "rs75017182": "DPYD",
    "rs4148323": "UGT1A1",  "rs35350960": "UGT1A1","rs887829": "UGT1A1",
    "rs8175347": "UGT1A1",
}


class VCFParseError(Exception):
    """Raised when a VCF file is malformed or cannot be parsed."""
    pass


def _parse_info(info_str: str) -> dict[str, str]:
    """Parse the VCF INFO field into a dict. Handles FLAG and KEY=VALUE entries."""
    result: dict[str, str] = {}
    if info_str in (".", ""):
        return result
    for token in info_str.split(";"):
        if "=" in token:
            k, _, v = token.partition("=")
            result[k.strip()] = v.strip()
        else:
            result[token.strip()] = "true"
    return result


def _determine_zygosity(gt_field: Optional[str]) -> str:
    """Determine zygosity from GT sample field."""
    if not gt_field:
        return "unknown"
    # Extract only the GT part (first field in sample column)
    gt = gt_field.split(":")[0]
    alleles = re.split(r"[/|]", gt)
    if len(alleles) < 2:
        return "unknown"
    if alleles[0] == alleles[1]:
        return "homozygous_alt" if alleles[0] != "0" else "homozygous_ref"
    return "heterozygous"


def _resolve_gene(info: dict, rsid: str) -> Optional[str]:
    """Resolve gene name from INFO tags or fallback rsID lookup."""
    # Explicit GENE tag (highest priority)
    gene = info.get("GENE") or info.get("gene")
    if gene:
        return gene.upper()
    # PharmVar / PharmGKB-style tags
    for tag in ("PHARMVAR_GENE", "PGX_GENE", "ANN"):
        if tag in info:
            val = info[tag]
            for g in PGX_GENES:
                if g in val.upper():
                    return g
    # rsID lookup (fallback)
    if rsid:
        return RS_TO_GENE.get(rsid.lower())
    return None


def parse_vcf(file_content: bytes | str, filename: str = "uploaded.vcf") -> tuple[list[ParsedVariant], dict]:
    """
    Parse a VCF v4.2 file and return a list of pharmacogenomic variants
    and metadata dict.

    Args:
        file_content: Raw bytes or string content of the VCF file.
        filename: Original filename (used in metadata).

    Returns:
        (list of ParsedVariant, metadata dict)

    Raises:
        VCFParseError: If the file is malformed or not a valid VCF.
    """
    if isinstance(file_content, bytes):
        try:
            text = file_content.decode("utf-8")
        except UnicodeDecodeError:
            text = file_content.decode("latin-1")
    else:
        text = file_content

    lines = text.splitlines()
    if not lines:
        raise VCFParseError("Empty VCF file provided.")

    # ── Header Parsing ─────────────────────────────────────────────────────
    vcf_version = "unknown"
    column_header: Optional[list[str]] = None
    has_sample = False
    data_lines: list[str] = []

    for i, line in enumerate(lines):
        if line.startswith("##fileformat=VCFv"):
            vcf_version = line.split("VCFv")[-1].strip()
        elif line.startswith("#CHROM"):
            cols = line.lstrip("#").split("\t")
            column_header = [c.strip() for c in cols]
            has_sample = len(column_header) > 9
        elif not line.startswith("#"):
            data_lines.append(line)

    if column_header is None:
        raise VCFParseError("Invalid VCF: missing #CHROM header line.")

    if not vcf_version.startswith("4"):
        logger.warning("VCF version %s detected — designed for v4.x.", vcf_version)

    # ── Variant Parsing ────────────────────────────────────────────────────
    total_variants = 0
    pgx_variants: list[ParsedVariant] = []
    parse_errors = 0

    for lineno, line in enumerate(data_lines, start=1):
        line = line.strip()
        if not line:
            continue
        total_variants += 1
        fields = line.split("\t")

        if len(fields) < 8:
            parse_errors += 1
            logger.warning("Line %d: insufficient columns (%d), skipping.", lineno, len(fields))
            continue

        try:
            chrom = fields[0].lstrip("chr")
            pos   = int(fields[1])
            rsid  = fields[2] if fields[2] != "." else None
            ref   = fields[3]
            alt   = fields[4].split(",")[0]  # take first ALT allele
            qual  = fields[5]
            filt  = fields[6]
            info  = _parse_info(fields[7])

            # Genotype
            if has_sample and len(fields) > 9:
                zygosity = _determine_zygosity(fields[9])
            else:
                zygosity = "unknown"

            # Extract PGx-specific tags from INFO
            star_allele = info.get("STAR") or info.get("star") or info.get("HAPLOTYPE")
            gene = _resolve_gene(info, rsid or "")

            # Only keep variants that belong to a PGx gene
            if gene and gene in PGX_GENES:
                pgx_variants.append(ParsedVariant(
                    gene=gene,
                    rsid=rsid,
                    star_allele=star_allele,
                    chromosome=chrom,
                    position=pos,
                    ref_allele=ref,
                    alt_allele=alt,
                    zygosity=zygosity,
                    raw_info={k: v for k, v in info.items() if k in ("GENE", "STAR", "RS", "AF", "DP", "HAPLOTYPE")},
                ))

        except (ValueError, IndexError) as exc:
            parse_errors += 1
            logger.warning("Line %d: parse error — %s", lineno, exc)
            continue

    if parse_errors > total_variants * 0.5 and total_variants > 0:
        raise VCFParseError(
            f"Too many malformed lines ({parse_errors}/{total_variants}). "
            "The file may not be a valid VCF."
        )

    metadata = {
        "file_name": filename,
        "vcf_version": vcf_version,
        "total_variants": total_variants,
        "pgx_variants_found": len(pgx_variants),
    }

    logger.info(
        "VCF parse complete: %d total variants, %d PGx variants identified.",
        total_variants, len(pgx_variants),
    )
    return pgx_variants, metadata
