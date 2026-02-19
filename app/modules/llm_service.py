"""
LLM Integration Module
Supports Google Gemini and OpenAI GPT for pharmacogenomic narrative generation.
"""
from __future__ import annotations

import json
import logging
import re
from typing import Optional

from app.config import get_settings
from app.models import (
    ParsedVariant, PhenotypePrediction, DrugRiskAssessment, LLMAnalysis, VariantCitation
)

logger = logging.getLogger(__name__)
settings = get_settings()


# ── Prompt Builder ────────────────────────────────────────────────────────────

def _build_prompt(
    variants: list[ParsedVariant],
    phenotypes: list[PhenotypePrediction],
    drug_risks: list[DrugRiskAssessment],
    drugs: list[str],
) -> str:
    """Construct the pharmacogenomics analysis prompt for the LLM."""

    # Summarize variants
    variant_lines = []
    for v in variants:
        parts = [f"Gene: {v.gene}"]
        if v.rsid:
            parts.append(f"rsID: {v.rsid}")
        if v.star_allele:
            parts.append(f"Star Allele: {v.star_allele}")
        parts.append(f"Position: chr{v.chromosome}:{v.position}")
        parts.append(f"REF/ALT: {v.ref_allele}/{v.alt_allele} ({v.zygosity})")
        variant_lines.append(" | ".join(parts))

    # Summarize phenotypes
    phenotype_lines = [
        f"- {p.gene}: {p.diplotype} → {p.phenotype} ({p.phenotype_full}), Activity Score: {p.activity_score}"
        for p in phenotypes
    ]

    # Summarize drug risks
    risk_lines = []
    for r in drug_risks:
        risk_lines.append(
            f"- Drug: {r.drug.upper()} | Gene: {r.gene} | Risk: {r.risk_category} "
            f"| Severity: {r.severity} | CPIC: {r.cpic_guideline}\n"
            f"  Recommendation: {r.recommendation_brief}"
        )

    prompt = f"""You are a clinical pharmacogenomics expert. A patient's VCF file has been analyzed. 
Provide a comprehensive pharmacogenomic clinical report based on the following findings.

=== DETECTED GENETIC VARIANTS ===
{chr(10).join(variant_lines) if variant_lines else "No pharmacogenomic variants detected."}

=== PHENOTYPE PREDICTIONS ===
{chr(10).join(phenotype_lines) if phenotype_lines else "No phenotype predictions."}

=== DRUG RISK ASSESSMENTS (CPIC-based) ===
{chr(10).join(risk_lines) if risk_lines else "No drug interactions assessed."}

=== DRUGS ANALYZED ===
{", ".join(drugs)}

Based on the above, provide a structured JSON response with EXACTLY this schema:
{{
  "clinical_summary": "<2-3 sentence patient-friendly summary of key pharmacogenomic findings>",
  "mechanism_explanation": "<detailed paragraph explaining the biological mechanisms for the identified variants and how they affect drug metabolism>",
  "dosing_recommendations": [
    "<specific actionable dosing recommendation for each drug>"
  ],
  "variant_citations": [
    {{
      "variant": "<rsID or star allele>",
      "pmid": "<PubMed ID if known, else null>",
      "note": "<brief clinical significance note>"
    }}
  ],
  "llm_confidence": <float 0-1 representing your confidence in these recommendations>
}}

Rules:
- Be specific and clinically actionable.
- For each drug, reference the relevant gene variant and phenotype.
- CPIC guideline levels should be mentioned where applicable.
- If a drug is CONTRAINDICATED or requires urgent dose adjustment, state it clearly.
- Return ONLY valid JSON, no markdown fences, no extra text.
"""
    return prompt


# ── Gemini Client ─────────────────────────────────────────────────────────────

async def _call_gemini(prompt: str) -> str:
    """Call Google Gemini API asynchronously using the new google-genai SDK."""
    try:
        from google import genai  # type: ignore
        from google.genai import types  # type: ignore

        client = genai.Client(api_key=settings.gemini_api_key)
        response = await client.aio.models.generate_content(
            model=settings.gemini_model,
            contents=prompt,
            config=types.GenerateContentConfig(
                temperature=0.2,
                max_output_tokens=2048,
                response_mime_type="application/json",
            ),
        )
        return response.text
    except Exception as e:
        logger.error("Gemini API error: %s", e)
        raise



# ── OpenAI Client ─────────────────────────────────────────────────────────────

async def _call_openai(prompt: str) -> str:
    """Call OpenAI-compatible API asynchronously (supports Groq, Ollama, etc.)."""
    try:
        from openai import AsyncOpenAI  # type: ignore
        # Initialize client with optional base_url for non-OpenAI providers
        client = AsyncOpenAI(
            api_key=settings.openai_api_key,
            base_url=settings.openai_base_url if settings.openai_base_url else None
        )
        response = await client.chat.completions.create(
            model=settings.openai_model,
            messages=[
                {
                    "role": "system",
                    "content": "You are a clinical pharmacogenomics expert. Always respond with valid JSON only.",
                },
                {"role": "user", "content": prompt},
            ],
            temperature=0.2,
            max_tokens=2048,
            response_format={"type": "json_object"},
        )
        return response.choices[0].message.content or "{}"
    except Exception as e:
        logger.error("OpenAI API error: %s", e)
        raise


# ── Response Parser ───────────────────────────────────────────────────────────

def _parse_llm_response(raw: str, model_name: str) -> LLMAnalysis:
    """Parse the raw LLM JSON string into an LLMAnalysis model."""
    # Strip markdown fences if present
    raw = raw.strip()
    raw = re.sub(r"^```[a-z]*\n?", "", raw)
    raw = re.sub(r"\n?```$", "", raw)
    raw = raw.strip()

    try:
        data = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.warning("LLM JSON parse failed: %s. Using fallback.", e)
        data = {}

    citations = []
    for c in data.get("variant_citations", []):
        if isinstance(c, dict):
            citations.append(VariantCitation(
                variant=c.get("variant", "unknown"),
                pmid=c.get("pmid"),
                note=c.get("note", ""),
            ))

    return LLMAnalysis(
        clinical_summary=data.get("clinical_summary", "LLM analysis unavailable."),
        mechanism_explanation=data.get("mechanism_explanation", ""),
        dosing_recommendations=data.get("dosing_recommendations", []),
        variant_citations=citations,
        llm_model_used=model_name,
        llm_confidence=float(data.get("llm_confidence", 0.0)),
    )


# ── Main Entry Point ──────────────────────────────────────────────────────────

async def generate_llm_analysis(
    variants: list[ParsedVariant],
    phenotypes: list[PhenotypePrediction],
    drug_risks: list[DrugRiskAssessment],
    drugs: list[str],
) -> LLMAnalysis:
    """
    Generate a pharmacogenomic clinical narrative using the configured LLM provider.
    Returns an LLMAnalysis object. Falls back gracefully on API errors.
    """
    if not settings.active_llm_key:
        logger.warning("No LLM API key configured — skipping LLM analysis.")
        return _fallback_analysis(settings.active_llm_model)

    prompt = _build_prompt(variants, phenotypes, drug_risks, drugs)

    try:
        if settings.llm_provider == "gemini":
            raw_response = await _call_gemini(prompt)
        else:
            raw_response = await _call_openai(prompt)
        return _parse_llm_response(raw_response, settings.active_llm_model)
    except Exception as e:
        logger.error("LLM call failed: %s — using fallback.", e)
        return _fallback_analysis(settings.active_llm_model)


def _fallback_analysis(model_name: str) -> LLMAnalysis:
    """Return a graceful fallback when LLM is unavailable."""
    return LLMAnalysis(
        clinical_summary=(
            "LLM analysis is currently unavailable (API key not configured or service error). "
            "Please review the drug risk assessments and phenotype predictions for clinical guidance."
        ),
        mechanism_explanation=(
            "Manual review required. Refer to CPIC guidelines at https://cpicpgx.org "
            "for detailed pharmacogenomic interaction information."
        ),
        dosing_recommendations=["Consult clinical pharmacist for personalized dosing guidance."],
        variant_citations=[],
        llm_model_used=model_name,
        llm_confidence=0.0,
    )
