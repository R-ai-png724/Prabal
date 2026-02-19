"""
Pharmacogenomic Analysis Module
Maps detected variants → diplotypes → phenotypes → drug risk assessments.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

from app.models import ParsedVariant, PhenotypePrediction, DrugRiskAssessment, OverallRisk

logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent.parent / "data"

_SEVERITY_ORDER = ["none", "low", "moderate", "high", "critical"]


def _load_json(filename: str) -> dict:
    path = DATA_DIR / filename
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# ── Module-level knowledge base (loaded once at startup) ─────────────────────
try:
    DIPLOTYPE_MAP: dict = _load_json("diplotype_phenotype_map.json")
    DRUG_INTERACTIONS: dict = _load_json("drug_gene_interactions.json")
    # Remove meta key
    DRUG_INTERACTIONS.pop("_meta", None)
except FileNotFoundError as e:
    logger.error("Knowledge base file not found: %s", e)
    DIPLOTYPE_MAP = {}
    DRUG_INTERACTIONS = {}


# ── Phenotype Prediction ──────────────────────────────────────────────────────

def _activity_score_to_phenotype(gene: str, score: float) -> tuple[str, str]:
    """Map activity score to phenotype code and full name."""
    gene_data = DIPLOTYPE_MAP.get(gene, {})
    for bucket in gene_data.get("activity_score_to_phenotype", []):
        if bucket["min"] <= score <= bucket["max"]:
            return bucket["phenotype"], bucket["full_name"]
    return "Unknown", "Unknown Metabolizer"


def _star_allele_activity(gene: str, star: str) -> float:
    """Return activity value for a star allele. Defaults to 1.0 (normal) if unknown."""
    gene_data = DIPLOTYPE_MAP.get(gene, {})
    alleles = gene_data.get("alleles", {})
    return alleles.get(star, {}).get("activity_value", 1.0)


def _build_diplotype(variants: list[ParsedVariant], gene: str) -> tuple[str, float, list[str]]:
    """
    Build the diplotype string and total activity score from variants for a gene.
    Returns (diplotype_str, activity_score, list_of_star_alleles).
    """
    gene_variants = [v for v in variants if v.gene == gene]
    if not gene_variants:
        return "*1/*1", 2.0, ["*1", "*1"]  # assume wildtype if no variants

    star_alleles: list[str] = []
    for v in gene_variants:
        if v.star_allele:
            star_alleles.append(v.star_allele)
        elif v.rsid and v.rsid in _build_additional_rsid_star_map():
            star_alleles.append(_build_additional_rsid_star_map()[v.rsid])

    # If we could not find any star alleles, treat as unknown *1 (wildtype)
    if not star_alleles:
        star_alleles = ["*1"]

    # Ensure we have exactly two haplotypes (diploid)
    if len(star_alleles) == 1:
        # Single variant found — assume the other allele is wildtype (*1)
        star_alleles = ["*1", star_alleles[0]]

    # Sort for canonical representation
    a1, a2 = star_alleles[0], star_alleles[1]
    activity = _star_allele_activity(gene, a1) + _star_allele_activity(gene, a2)
    diplotype = f"{a1}/{a2}"
    return diplotype, activity, star_alleles


def _build_additional_rsid_star_map() -> dict[str, str]:
    """Common rsID → star allele mapping for important variants."""
    return {
        "rs3892097":  "*4",
        "rs1065852":  "*10",
        "rs5030655":  "*6",
        "rs16947":    "*2",
        "rs28371706": "*41",
        "rs28371725": "*9",
        "rs4244285":  "*2",
        "rs4986893":  "*3",
        "rs1799853":  "*2",
        "rs1057910":  "*3",
        "rs1800460":  "*3B",
        "rs1142345":  "*3C",
        "rs1800462":  "*2",
        "rs3918290":  "*2A",
        "rs55886062": "*13",
        "rs4148323":  "*6",
        "rs8175347":  "*28",
    }


def identify_phenotypes(variants: list[ParsedVariant]) -> list[PhenotypePrediction]:
    """
    For each unique gene found in variants, compute the diplotype and phenotype.
    Returns a list of PhenotypePrediction objects.
    """
    predictions: list[PhenotypePrediction] = []
    genes_seen: set[str] = {v.gene for v in variants}

    for gene in genes_seen:
        diplotype, activity_score, _ = _build_diplotype(variants, gene)
        phenotype_code, phenotype_full = _activity_score_to_phenotype(gene, activity_score)

        # Confidence based on whether we had explicit star allele annotations
        gene_variants = [v for v in variants if v.gene == gene]
        has_star = any(v.star_allele for v in gene_variants)
        confidence = 0.92 if has_star else 0.72

        predictions.append(PhenotypePrediction(
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype_code,
            phenotype_full=phenotype_full,
            activity_score=activity_score,
            confidence=confidence,
        ))
        logger.debug("Gene %s: diplotype=%s, phenotype=%s, score=%.1f",
                     gene, diplotype, phenotype_code, activity_score)

    return predictions


# ── Drug Risk Prediction ──────────────────────────────────────────────────────

def predict_drug_risks(
    phenotypes: list[PhenotypePrediction],
    drugs: list[str],
) -> list[DrugRiskAssessment]:
    """
    For each (drug, gene) interaction pair, predict the risk category
    based on the patient's phenotype and CPIC rules.
    Returns a list of DrugRiskAssessment objects.
    """
    assessments: list[DrugRiskAssessment] = []
    phenotype_map: dict[str, PhenotypePrediction] = {p.gene: p for p in phenotypes}

    for drug_raw in drugs:
        drug = drug_raw.strip().lower()
        if drug not in DRUG_INTERACTIONS:
            # Drug not in knowledge base → Unknown risk
            assessments.append(DrugRiskAssessment(
                drug=drug,
                gene="Unknown",
                risk_category="Unknown",
                severity="low",
                confidence_score=0.0,
                cpic_guideline="Not in CPIC database",
                recommendation_brief=(
                    f"No pharmacogenomic interaction data found for '{drug}'. "
                    "Consult prescribing information and clinical pharmacist."
                ),
                mechanism=None,
            ))
            continue

        interaction = DRUG_INTERACTIONS[drug]
        gene = interaction["gene"]
        cpic_level = interaction.get("cpic_level", "Unknown")
        phenotype_rules: dict = interaction.get("phenotype_rules", {})

        # Get patient's phenotype for this gene
        patient_phenotype = phenotype_map.get(gene)
        if patient_phenotype is None:
            phenotype_code = "NM"  # assume Normal if no variants found in this gene
            pheno_confidence = 0.65
        else:
            phenotype_code = patient_phenotype.phenotype
            pheno_confidence = patient_phenotype.confidence

        rule = phenotype_rules.get(phenotype_code, phenotype_rules.get("NM"))
        if rule:
            risk_category     = rule["risk_category"]
            severity          = rule["severity"]
            confidence_base   = rule["confidence_base"]
            recommendation    = rule["recommendation"]
            mechanism         = rule.get("mechanism", "")
        else:
            risk_category   = "Unknown"
            severity        = "low"
            confidence_base = 0.5
            recommendation  = "Consult clinical pharmacist for dosing guidance."
            mechanism       = None

        # Combine confidence: phenotype detection × rule base confidence
        final_confidence = round(pheno_confidence * confidence_base, 3)

        assessments.append(DrugRiskAssessment(
            drug=drug,
            gene=gene,
            risk_category=risk_category,
            severity=severity,
            confidence_score=final_confidence,
            cpic_guideline=f"CPIC Level {cpic_level}",
            recommendation_brief=recommendation,
            mechanism=mechanism,
        ))

    return assessments


# ── Overall Risk Aggregation ──────────────────────────────────────────────────

def compute_overall_risk(assessments: list[DrugRiskAssessment]) -> OverallRisk:
    """Aggregate all drug risk assessments into an overall patient risk level."""
    if not assessments:
        return OverallRisk(level="none", flags=[])

    max_severity = "none"
    flags: list[str] = []

    for a in assessments:
        sev = a.severity
        if _SEVERITY_ORDER.index(sev) > _SEVERITY_ORDER.index(max_severity):
            max_severity = sev

        if sev in ("high", "critical"):
            flags.append(f"{a.risk_category} required for {a.drug} ({a.gene}: {a.severity} risk)")
        elif a.risk_category in ("Toxic", "Ineffective"):
            flags.append(f"{a.drug} may be {a.risk_category.lower()} — {a.gene} variant detected")

    return OverallRisk(level=max_severity, flags=list(dict.fromkeys(flags)))  # deduplicate
