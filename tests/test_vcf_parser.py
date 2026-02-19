"""
Unit tests for VCF Parser module.
Run with: python -m pytest tests/test_vcf_parser.py -v
"""
import pytest
from app.modules.vcf_parser import parse_vcf, VCFParseError

SAMPLE_VCF = """\
##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=STAR,Number=1,Type=String,Description="Star allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr22\t42526694\trs3892097\tC\tT\t99\tPASS\tGENE=CYP2D6;STAR=*4\tGT\t0/1
chr10\t96741053\trs4244285\tG\tA\t99\tPASS\tGENE=CYP2C19;STAR=*2\tGT\t0/0
chr1\t12345\t.\tA\tG\t70\tPASS\tDP=30\tGT\t0/1
"""


def test_parse_basic_vcf():
    variants, meta = parse_vcf(SAMPLE_VCF.encode())
    assert meta["vcf_version"] == "4.2"
    assert meta["total_variants"] == 3
    assert meta["pgx_variants_found"] == 2
    assert len(variants) == 2


def test_variant_fields():
    variants, _ = parse_vcf(SAMPLE_VCF.encode())
    v = variants[0]
    assert v.gene == "CYP2D6"
    assert v.rsid == "rs3892097"
    assert v.star_allele == "*4"
    assert v.chromosome == "22"
    assert v.position == 42526694
    assert v.ref_allele == "C"
    assert v.alt_allele == "T"
    assert v.zygosity == "heterozygous"


def test_homozygous_ref():
    variants, _ = parse_vcf(SAMPLE_VCF.encode())
    # Second variant is 0/0 → homozygous_ref
    cyp2c19_variant = next(v for v in variants if v.gene == "CYP2C19")
    assert cyp2c19_variant.zygosity == "homozygous_ref"


def test_rsid_fallback_gene_resolution():
    vcf = """\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr22\t42526694\trs3892097\tC\tT\t99\tPASS\t.\tGT\t0/1
"""
    variants, _ = parse_vcf(vcf.encode())
    # rs3892097 is known to be CYP2D6 via fallback map
    assert len(variants) == 1
    assert variants[0].gene == "CYP2D6"


def test_empty_file_raises():
    with pytest.raises(VCFParseError, match="Empty"):
        parse_vcf(b"")


def test_missing_header_raises():
    vcf = "chr1\t100\t.\tA\tG\t99\tPASS\t.\tGT\t0/1\n"
    with pytest.raises(VCFParseError, match="#CHROM"):
        parse_vcf(vcf.encode())


def test_bytes_and_string_input_equivalent():
    v1, _ = parse_vcf(SAMPLE_VCF.encode())
    v2, _ = parse_vcf(SAMPLE_VCF)
    assert len(v1) == len(v2)
    assert v1[0].gene == v2[0].gene


def test_filename_in_metadata():
    _, meta = parse_vcf(SAMPLE_VCF.encode(), filename="patient_001.vcf")
    assert meta["file_name"] == "patient_001.vcf"
