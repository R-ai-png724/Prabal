"""
Health & Info Routes
"""
from fastapi import APIRouter, Depends
from app.config import Settings, get_settings
from app.models import HealthResponse
from app.modules.pgx_analyzer import DRUG_INTERACTIONS, DIPLOTYPE_MAP

router = APIRouter(prefix="/api/v1", tags=["Health"])


@router.get(
    "/health",
    response_model=HealthResponse,
    summary="Health Check",
    description="Returns system health status and supported genes/drugs.",
)
async def health(settings: Settings = Depends(get_settings)):
    return HealthResponse(
        status="healthy",
        version=settings.app_version,
        llm_provider=settings.llm_provider,
        llm_model=settings.active_llm_model,
        genes_supported=sorted(DIPLOTYPE_MAP.keys()),
        drugs_supported=sorted(DRUG_INTERACTIONS.keys()),
    )
