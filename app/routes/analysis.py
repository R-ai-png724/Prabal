"""
Analysis API Route — POST /api/v1/analyze
"""
from __future__ import annotations

import time
import logging
from datetime import datetime, timezone

from fastapi import APIRouter, File, Form, UploadFile, HTTPException, Depends
from fastapi.responses import JSONResponse

from app.config import Settings, get_settings
from app.models import AnalysisResponse, VCFMetadata, OverallRisk
from app.modules.vcf_parser import parse_vcf, VCFParseError
from app.modules.pgx_analyzer import identify_phenotypes, predict_drug_risks, compute_overall_risk
from app.modules.llm_service import generate_llm_analysis

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/v1", tags=["Analysis"])


@router.post(
    "/analyze",
    response_model=AnalysisResponse,
    summary="Pharmacogenomic Risk Analysis",
    description=(
        "Upload a VCF v4.2 file and a comma-separated list of drug names. "
        "Returns a comprehensive pharmacogenomic risk assessment with LLM clinical narrative."
    ),
)
async def analyze(
    vcf_file: UploadFile = File(..., description="VCF v4.2 file containing patient genomic variants"),
    drugs: str = Form(..., description="Comma-separated drug names, e.g. 'codeine,warfarin'"),
    patient_id: str | None = Form(None, description="Optional patient identifier (not stored)"),
    skip_llm: bool = Form(False, description="Set to true to skip LLM analysis (faster)"),
    settings: Settings = Depends(get_settings),
):
    start_time = time.monotonic()

    # ── 1. Validate upload size ───────────────────────────────────────────
    content = await vcf_file.read()
    if len(content) > settings.max_vcf_size_bytes:
        raise HTTPException(
            status_code=413,
            detail=f"VCF file exceeds maximum size of {settings.max_vcf_size_mb} MB.",
        )

    # ── 2. Parse VCF ─────────────────────────────────────────────────────
    try:
        variants, vcf_meta_raw = parse_vcf(content, filename=vcf_file.filename or "uploaded.vcf")
    except VCFParseError as e:
        raise HTTPException(status_code=422, detail=f"VCF Parse Error: {str(e)}")
    except Exception as e:
        logger.exception("Unexpected VCF parse error.")
        raise HTTPException(status_code=500, detail=f"Internal error during VCF parsing: {str(e)}")

    vcf_metadata = VCFMetadata(**vcf_meta_raw)

    # ── 3. Parse drug list ────────────────────────────────────────────────
    drug_list = [d.strip().lower() for d in drugs.split(",") if d.strip()]
    if not drug_list:
        raise HTTPException(status_code=422, detail="At least one drug name must be provided.")

    # ── 4. PGx Analysis ───────────────────────────────────────────────────
    phenotypes = identify_phenotypes(variants)
    drug_risks = predict_drug_risks(phenotypes, drug_list)
    overall_risk = compute_overall_risk(drug_risks)

    # ── 5. LLM Analysis (optional) ────────────────────────────────────────
    llm_analysis = None
    if not skip_llm:
        llm_analysis = await generate_llm_analysis(variants, phenotypes, drug_risks, drug_list)

    # ── 6. Assemble Response ──────────────────────────────────────────────
    elapsed_ms = int((time.monotonic() - start_time) * 1000)

    response = AnalysisResponse(
        status="success",
        patient_id=patient_id,
        analysis_timestamp=datetime.now(timezone.utc).isoformat(),
        vcf_metadata=vcf_metadata,
        detected_variants=variants,
        phenotype_predictions=phenotypes,
        drug_risk_assessments=drug_risks,
        llm_analysis=llm_analysis,
        overall_risk=overall_risk,
        processing_time_ms=elapsed_ms,
    )

    logger.info(
        "Analysis complete | patient=%s | variants=%d | drugs=%s | overall_risk=%s | time=%dms",
        patient_id or "anon",
        len(variants),
        drug_list,
        overall_risk.level,
        elapsed_ms,
    )

    return response


@router.post(
    "/analyze/batch",
    summary="Batch Pharmacogenomic Analysis (no LLM)",
    description="Faster endpoint — runs PGx analysis without LLM for batch/screening use cases.",
)
async def analyze_batch(
    vcf_file: UploadFile = File(...),
    drugs: str = Form(...),
    patient_id: str | None = Form(None),
    settings: Settings = Depends(get_settings),
):
    """Same as /analyze but always skips LLM for speed."""
    # Delegate to main endpoint with skip_llm=True
    return await analyze(
        vcf_file=vcf_file,
        drugs=drugs,
        patient_id=patient_id,
        skip_llm=True,
        settings=settings,
    )
