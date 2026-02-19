"""
Unit tests for PGx Analyzer module.
Run with: python -m pytest tests/test_pgx_analyzer.py -v
"""
import pytest
from app.models import ParsedVariant
from app.modules.pgx_analyzer import (
    identify_phenotypes, predict_drug_risks, compute_overall_risk
)


def make_variant(gene, star_allele=None, rsid=None, zygosity="heterozygous"):
    return ParsedVariant(
        gene=gene,
        rsid=rsid,
        star_allele=star_allele,
        chromosome="22",
        position=1000,
        ref_allele="C",
        alt_allele="T",
        zygosity=zygosity,
    )


class TestPhenotypeIdentification:

    def test_cyp2d6_poor_metabolizer(self):
        """*4/*4 homozygous → PM (activity 0.0+0.0=0.0)"""
        variants = [
            make_variant("CYP2D6", "*4", zygosity="homozygous_alt"),
            make_variant("CYP2D6", "*4"),
        ]
        phenos = identify_phenotypes(variants)
        assert len(phenos) == 1
        assert phenos[0].gene == "CYP2D6"
        assert phenos[0].phenotype == "PM"
        assert phenos[0].activity_score == 0.0

    def test_cyp2d6_intermediate_metabolizer(self):
        """*1/*4 → IM (activity 1.0+0.0=1.0)"""
        variants = [make_variant("CYP2D6", "*4", rsid="rs3892097")]
        phenos = identify_phenotypes(variants)
        assert phenos[0].phenotype == "IM"

    def test_cyp2c19_normal(self):
        """*1/*1 (assumed wildtype when no variants) → NM"""
        # Empty variants for CYP2C19 → wildtype assumed
        phenos = identify_phenotypes([])
        assert phenos == []

    def test_multiple_genes(self):
        variants = [
            make_variant("CYP2D6", "*4"),
            make_variant("CYP2C19", "*2"),
        ]
        phenos = identify_phenotypes(variants)
        genes = {p.gene for p in phenos}
        assert "CYP2D6" in genes
        assert "CYP2C19" in genes


class TestDrugRiskPrediction:

    def test_codeine_pm_ineffective(self):
        variants = [make_variant("CYP2D6", "*4"), make_variant("CYP2D6", "*4")]
        phenos = identify_phenotypes(variants)
        risks = predict_drug_risks(phenos, ["codeine"])
        assert len(risks) == 1
        r = risks[0]
        assert r.drug == "codeine"
        assert r.gene == "CYP2D6"
        assert r.risk_category == "Ineffective"
        assert r.severity == "high"

    def test_codeine_um_toxic(self):
        variants = [make_variant("CYP2D6", "*1xN")]
        phenos = identify_phenotypes(variants)
        risks = predict_drug_risks(phenos, ["codeine"])
        r = risks[0]
        assert r.risk_category == "Toxic"
        assert r.severity == "critical"

    def test_unknown_drug(self):
        risks = predict_drug_risks([], ["unknowndrug123"])
        assert risks[0].risk_category == "Unknown"
        assert risks[0].gene == "Unknown"

    def test_warfarin_normal(self):
        """No CYP2C9 variants → assume NM → Safe"""
        risks = predict_drug_risks([], ["warfarin"])
        r = risks[0]
        assert r.drug == "warfarin"
        assert r.risk_category == "Safe"

    def test_fluorouracil_pm_toxic(self):
        variants = [make_variant("DPYD", "*2A"), make_variant("DPYD", "*2A")]
        phenos = identify_phenotypes(variants)
        risks = predict_drug_risks(phenos, ["fluorouracil"])
        r = risks[0]
        assert r.risk_category == "Toxic"
        assert r.severity == "critical"


class TestOverallRisk:

    def test_overall_critical_from_toxic(self):
        risks = predict_drug_risks(
            identify_phenotypes([make_variant("DPYD", "*2A"), make_variant("DPYD", "*2A")]),
            ["fluorouracil"]
        )
        overall = compute_overall_risk(risks)
        assert overall.level == "critical"

    def test_overall_none_when_safe(self):
        risks = predict_drug_risks([], ["warfarin"])
        overall = compute_overall_risk(risks)
        # Warfarin NM → Safe → severity none
        assert overall.level == "none"

    def test_flags_populated_for_high_severity(self):
        risks = predict_drug_risks(
            identify_phenotypes([make_variant("TPMT", "*2"), make_variant("TPMT", "*2")]),
            ["azathioprine"]
        )
        overall = compute_overall_risk(risks)
        assert len(overall.flags) > 0
