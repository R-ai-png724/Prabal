"""
Pydantic models for request/response schemas.
"""
from pydantic import BaseModel, Field
from typing import Optional
from datetime import datetime


# ── VCF Parsing ──────────────────────────────────────────────────────────────

class ParsedVariant(BaseModel):
    gene: str
    rsid: Optional[str] = None
    star_allele: Optional[str] = None
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    zygosity: str  # homozygous_ref | heterozygous | homozygous_alt
    raw_info: dict = Field(default_factory=dict)


class VCFMetadata(BaseModel):
    file_name: str
    vcf_version: str
    total_variants: int
    pgx_variants_found: int


# ── PGx Analysis ─────────────────────────────────────────────────────────────

class PhenotypePrediction(BaseModel):
    gene: str
    diplotype: str
    phenotype: str        # PM | IM | NM | UM | RM
    phenotype_full: str
    activity_score: float
    confidence: float


class DrugRiskAssessment(BaseModel):
    drug: str
    gene: str
    risk_category: str    # Safe | Adjust Dosage | Toxic | Ineffective | Unknown
    severity: str         # none | low | moderate | high | critical
    confidence_score: float
    cpic_guideline: str
    recommendation_brief: str
    mechanism: Optional[str] = None


class OverallRisk(BaseModel):
    level: str            # none | low | moderate | high | critical
    flags: list[str]


# ── LLM Analysis ─────────────────────────────────────────────────────────────

class VariantCitation(BaseModel):
    variant: str
    pmid: Optional[str] = None
    note: str


class LLMAnalysis(BaseModel):
    clinical_summary: str
    mechanism_explanation: str
    dosing_recommendations: list[str]
    variant_citations: list[VariantCitation]
    llm_model_used: str
    llm_confidence: float


# ── Full Analysis Response ────────────────────────────────────────────────────

class AnalysisResponse(BaseModel):
    status: str
    patient_id: Optional[str] = None
    analysis_timestamp: str
    vcf_metadata: VCFMetadata
    detected_variants: list[ParsedVariant]
    phenotype_predictions: list[PhenotypePrediction]
    drug_risk_assessments: list[DrugRiskAssessment]
    llm_analysis: Optional[LLMAnalysis] = None
    overall_risk: OverallRisk
    processing_time_ms: int


class HealthResponse(BaseModel):
    status: str
    version: str
    llm_provider: str
    llm_model: str
    genes_supported: list[str]
    drugs_supported: list[str]
