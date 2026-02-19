"""
FastAPI Application Factory
"""
from __future__ import annotations

import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from app.config import get_settings
from app.routes.analysis import router as analysis_router
from app.routes.health import router as health_router

settings = get_settings()

logging.basicConfig(
    level=logging.DEBUG if settings.debug else logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
)
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("🧬 PGx Risk Prediction API starting up...")
    logger.info("LLM Provider: %s | Model: %s", settings.llm_provider, settings.active_llm_model)
    yield
    logger.info("🧬 PGx Risk Prediction API shutting down.")


def create_app() -> FastAPI:
    app = FastAPI(
        title=settings.app_name,
        version=settings.app_version,
        description=(
            "End-to-End Pharmacogenomic Risk Prediction System. "
            "Upload a VCF file and drug names to receive CPIC-aligned risk assessments "
            "and LLM-generated clinical narratives."
        ),
        docs_url="/docs",
        redoc_url="/redoc",
        lifespan=lifespan,
    )

    # ── CORS ──────────────────────────────────────────────────────────────────
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # Tighten in production
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # ── Global Exception Handler ──────────────────────────────────────────────
    @app.exception_handler(Exception)
    async def global_exception_handler(request: Request, exc: Exception):
        logger.exception("Unhandled exception: %s", exc)
        return JSONResponse(
            status_code=500,
            content={"status": "error", "detail": "An unexpected internal error occurred."},
        )

    # ── Routers ───────────────────────────────────────────────────────────────
    app.include_router(analysis_router)
    app.include_router(health_router)

    @app.get("/", tags=["Root"])
    async def root():
        return {
            "name": settings.app_name,
            "version": settings.app_version,
            "docs": "/docs",
            "health": "/api/v1/health",
            "analyze": "POST /api/v1/analyze",
        }

    return app


app = create_app()
