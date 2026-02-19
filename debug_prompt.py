import sys
import asyncio
from app.modules.vcf_parser import parse_vcf
from app.modules.pgx_analyzer import identify_phenotypes, predict_drug_risks
from app.modules.llm_service import _build_prompt
from pathlib import Path

# Mock data path
DATA_DIR = Path("data")
SAMPLE_VCF = DATA_DIR / "sample_test.vcf"

async def main():
    print("--- Generating LLM Prompt for Testing ---\n")
    
    # 1. Read VCF
    if not SAMPLE_VCF.exists():
        print(f"Error: {SAMPLE_VCF} not found.")
        return

    with open(SAMPLE_VCF, "rb") as f:
        content = f.read()

    variants, _ = parse_vcf(content, "sample.vcf")
    
    # 2. Analyze
    phenotypes = identify_phenotypes(variants)
    drugs = ["warfarin", "codeine", "fluorouracil"]
    drug_risks = predict_drug_risks(phenotypes, drugs)

    # 3. Build Prompt
    prompt = _build_prompt(variants, phenotypes, drug_risks, drugs)
    
    print("\nCopy the text below into Google AI Studio (https://aistudio.google.com/):\n")
    print("="*80)
    print(prompt)
    print("="*80)

if __name__ == "__main__":
    asyncio.run(main())
