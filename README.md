# 🧬 Pharmacogenomic Risk Prediction System — Backend

End-to-end backend API for pharmacogenomic (PGx) risk prediction. Upload a patient VCF file and a list of drugs to receive CPIC-aligned risk assessments powered by an LLM clinical narrative.

---

## ⚡ Quick Start

```bash
# 1. Clone and navigate to project
cd /path/to/PRABAL

# 2. Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Configure environment
cp .env.example .env
# Edit .env and add your GEMINI_API_KEY (or OPENAI_API_KEY)

# 5. Run the server
python main.py
```

API will be live at: **http://localhost:8000**  
Interactive docs: **http://localhost:8000/docs**

---

## 🔌 API Endpoints

### `POST /api/v1/analyze`
Main analysis endpoint. Accepts a VCF file and drug list.

**Form fields:**
| Field | Type | Required | Description |
|---|---|---|---|
| `vcf_file` | file | ✅ | VCF v4.2 file |
| `drugs` | string | ✅ | Comma-separated drug names |
| `patient_id` | string | ❌ | Optional ID (not stored) |
| `skip_llm` | bool | ❌ | Skip LLM for faster response |

**Example (curl):**
```bash
curl -X POST http://localhost:8000/api/v1/analyze \
  -F "vcf_file=@data/sample_test.vcf" \
  -F "drugs=codeine,warfarin,fluorouracil" \
  -F "patient_id=PT-001"
```

### `POST /api/v1/analyze/batch`
Same as `/analyze` but always skips LLM (fastest mode).

### `GET /api/v1/health`
Returns system status, supported genes and drugs.

---

## 🧬 Supported Genes & Drug Coverage

| Gene | Function | Example Drugs |
|---|---|---|
| **CYP2D6** | Opioid & antidepressant metabolism | Codeine, Tramadol, Amitriptyline |
| **CYP2C19** | Antiplatelet & PPI metabolism | Clopidogrel, Omeprazole |
| **CYP2C9** | Anticoagulant metabolism | Warfarin |
| **TPMT** | Thiopurine metabolism | Azathioprine, Mercaptopurine |
| **DPYD** | Fluoropyrimidine metabolism | Fluorouracil, Capecitabine |
| **UGT1A1** | SN-38 glucuronidation | Irinotecan |

---

## 📋 Response Schema

```json
{
  "status": "success",
  "patient_id": "PT-001",
  "analysis_timestamp": "2024-01-15T10:30:00Z",
  "vcf_metadata": { "file_name": "...", "vcf_version": "4.2", "total_variants": 10, "pgx_variants_found": 5 },
  "detected_variants": [{ "gene": "CYP2D6", "rsid": "rs3892097", "star_allele": "*4", ... }],
  "phenotype_predictions": [{ "gene": "CYP2D6", "diplotype": "*1/*4", "phenotype": "IM", ... }],
  "drug_risk_assessments": [{ "drug": "codeine", "risk_category": "Adjust Dosage", "severity": "moderate", ... }],
  "llm_analysis": { "clinical_summary": "...", "dosing_recommendations": [...], ... },
  "overall_risk": { "level": "moderate", "flags": ["..."] },
  "processing_time_ms": 1240
}
```

Risk categories: `Safe` | `Adjust Dosage` | `Toxic` | `Ineffective` | `Unknown`  
Severity levels: `none` | `low` | `moderate` | `high` | `critical`

---

## 🔧 LLM Configuration

Set `LLM_PROVIDER` in `.env`:

| Provider | Env Var | Model |
|---|---|---|
| Google Gemini (default) | `GEMINI_API_KEY` | `gemini-1.5-flash` |
| OpenAI | `OPENAI_API_KEY` | `gpt-4o-mini` |

The system **gracefully degrades** — if no API key is set, you still get full PGx analysis without the LLM narrative.

---

## 🧪 Running Tests

```bash
# All tests
python -m pytest tests/ -v

# Specific modules
python -m pytest tests/test_vcf_parser.py -v
python -m pytest tests/test_pgx_analyzer.py -v
```

---

## 📁 Project Structure

```
PRABAL/
├── main.py                         # Entry point (uvicorn)
├── requirements.txt
├── .env.example                    # Environment template
├── app/
│   ├── app.py                      # FastAPI factory
│   ├── config.py                   # Settings (pydantic-settings)
│   ├── models.py                   # Pydantic request/response models
│   ├── modules/
│   │   ├── vcf_parser.py           # VCF v4.2 parser
│   │   ├── pgx_analyzer.py         # Diplotype → phenotype → risk engine
│   │   └── llm_service.py          # Gemini / OpenAI integration
│   └── routes/
│       ├── analysis.py             # POST /api/v1/analyze
│       └── health.py               # GET /api/v1/health
├── data/
│   ├── diplotype_phenotype_map.json
│   ├── drug_gene_interactions.json # CPIC-based rules
│   └── sample_test.vcf             # Test VCF file
└── tests/
    ├── test_vcf_parser.py
    └── test_pgx_analyzer.py
```

---

## 📖 References

- [CPIC Guidelines](https://cpicpgx.org)
- [PharmGKB](https://www.pharmgkb.org)
- [PharmVar (Star Alleles)](https://www.pharmvar.org)
- [FastAPI Documentation](https://fastapi.tiangolo.com)
