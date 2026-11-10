<p align="center">
  <img src="src/img/magnamined_icontext.png" alt="Magna-Mined: AI-curated literature miner for Magnaporthe oryzae research" width="600">
</p>

An automated literature mining, experimental evidence extraction, and gene annotation tool developed for the automated curation of ***Magnaporthe oryzae*** literature and ongoing maintenance of [**MagnaGenes v2.0**](https://github.com/nehasahu486/magnagenes_v2) ([web portal](https://nehasahu486.github.io/magnagenes_v2/)).

---

## 1. Purpose & Core Capabilities

[**MagnaGenes v.1.0**](https://zenodo.org/records/4647766) (Foster et al., 2021) catalogs **12,522 genes** in *M. oryzae* (the rice blast fungus, also known as *Pyricularia oryzae*):
* **~1,636 genes** already hold phenotype annotations in MagnaGenes v.1.0 (1,274 with mutant phenotypes + 362 tested as `none`).
* **10,886 genes** are marked **`Not reported`** (unclassified).

Magna-Mined keeps MagnaGenes v2.0 current automatically (without manual literature searches or paid AI API keys) through two core functions:

1. **Primary Function: Curating Unclassified Genes (Default)**:
   Every quarter or six months, it fetches newly published *M. oryzae* experimental genetics papers, skips the ~1,636 genes already classified in `MagnaGenes.csv`, and extracts phenotype and knockout records strictly for the **10,886 unclassified (`Not reported`) genes**.
2. **Secondary Function: Double-Checking Classified Genes for Additional Functions**:
   When a new paper demonstrates an **additional biological role** for an already-classified gene (for example, a gene previously listed only under `virulence` is shown to also regulate `conidiation`), `double-check` preserves the original MagnaGenes v.1.0 classification and logs the new finding and citation into **`additional_functions`** and **`additional_function_references`**.

> **Reference (MagnaGenes v.1.0):** Foster AJ, Were VM, Yan X, Win J, Harant A, Langner T, Kamoun S, Talbot NJ (2021). *MagnaGenes (v.1.0): an open science database of gene function studies in the blast fungus Magnaporthe oryzae*. Zenodo. [10.5281/zenodo.4647766](https://doi.org/10.5281/zenodo.4647766)

---

## 2. How Magna-Mined Works

Instead of querying Europe PMC 12,522 separate times (once per gene), Magna-Mined follows a **Paper-First Workflow** that scans newly published *M. oryzae* papers directly and extracts findings for every gene studied in each paper.

### Does It Scan Just the Abstract or the Entire Paper?
It scans **the full paper whenever full text is accessible**, and only falls back to the abstract for paywalled articles:

1. **Priority 1: Europe PMC Full-Text XML (`jats_xml`)**:
   When a paper has an open-access PMC ID (`pmcid`), Magna-Mined downloads the structured `fullTextXML` from Europe PMC. It uses all sections (title, abstract, introduction/other, results, discussion, methods) to detect gene names and new inline aliases, and extracts phenotype and mutant evidence specifically from the **Results** and **Methods** sections.
2. **Priority 2: Full-Text PDF Extraction (`pdf`)**:
   When XML is unavailable, it checks for an open-access PDF link via Europe PMC or Unpaywall (when `MAGNAMINED_CONTACT_EMAIL` is configured), or reads a locally supplied PDF via `ingest-pdf`. It converts the PDF into reading-order text with `pymupdf4llm`, splits out **Results** and **Methods** headings (or falls back to the full PDF body text if section headings cannot be split cleanly), and runs extraction on those sections.
3. **Priority 3: Abstract-Only Fallback (`abstract_only`)**:
   If a paper is behind a publisher paywall and neither XML nor a PDF is available, Magna-Mined scans the **title and abstract** only and automatically caps its evidence rating at **`tier_4_low`** so curators know the full text was not inspected.

### Step-by-Step Pipeline
```
Europe PMC Paper Stream (year range or quarterly/6-month publication window)
   │
   ▼ 1. Processed-Paper Tracking (`processed_papers`)
   │    Records evaluated PMIDs/PMCIDs in SQLite and skips already-checked papers
   ▼ 2. Abstract Screening
   │    Reads title & abstract first; skips reviews and non-experimental surveys
   ▼ 3. Full-Text Retrieval (XML -> PDF -> Abstract fallback) & Alias Matching
   │    Maps MGG_#####, MoGene, Δgene/∆gene, Broad v6 IDs, and newly coined symbols
   │    defined inside the paper (e.g., "MoPct1 (MGG_01003)")
   ▼ 4. Focused Snippet Windowing (1,500 characters)
   │    Extracts 1,500-character windows around each gene mention in Results/Methods
   ▼ 5. Sentence-Level Text-Based Search & Exact-Quote Check (`extractor.py`)
   │    Extracts phenotype, knockout status, and checks a supporting quote from the text
   ▼ 6. Multi-Paper Ranking, Conflict Flagging & MagnaGenes v2.0 Sync (`paper_miner.py`, `reconciler.py`, `explorer.py`, `magnagenes_sync.py`)
        Combines evidence across papers, marks `knockout_reported` / `additional_functions`,
        flags conflicts, and exports HTML charts + dated `fin_data_YYYY-MM-DD.tsv` updates
```

### How It Prevents Hallucinations & AI-Slop
1. **Text-Based Pattern Search by Default (No Generative AI Required)**:
   Without `GEMINI_API_KEY`, Magna-Mined does not call any LLM. It uses deterministic **text-based pattern search** directly over sentences from the downloaded paper, so no claims or biological details can be invented.
2. **Mandatory Exact-Quote Check**:
   Every extraction (whether from the built-in text-based search or optional Gemini LLM) must carry a supporting quote copied from the paper. Before anything is stored in SQLite, an **exact-quote check** verifies that the quote matches text inside the actual paper (exact or whitespace/punctuation-normalized match). Any record whose quote cannot be found in the paper is automatically rejected.
3. **Strict Sentence Scoping & Irrelevant-Mention Filtering**:
   To prevent phenotypes from one gene bleeding onto another gene mentioned nearby, extraction looks only at sentences that explicitly name the target gene or its mutant (plus an immediately following complementation/rescue sentence only if no other gene's mutant is named). Sentences where the gene is merely a `qRT-PCR`/`RNA-seq` expression target, fluorescent marker (`-GFP`), housekeeping/reference gene, protein domain reference, part of a multi-gene heatmap/list (≥3 `MGG_` IDs), or a sentence describing a *different* `Δgene` mutant are filtered out automatically.
4. **Confidence Tiers & Curator Audit Trail**:
   Every record is assigned a clear evidence tier based on experimental proof and text availability:
   * **`tier_1_high`**: Targeted gene knockout, deletion, or CRISPR/Cas9 mutation **with genetic complementation rescue** verified in full text (XML or PDF).
   * **`tier_2_high`**: Targeted knockout, deletion, or CRISPR/Cas9 mutant phenotype **without complementation** reported in full text.
   * **`tier_3_medium`**: Overexpression, RNAi silencing, site-directed/point mutagenesis, or biochemical assay without a full knockout.
   * **`tier_4_low`**: Evidence extracted **only from the abstract** (for paywalled papers without XML/PDF), expression correlation, or computational predictions.

### How Multiple Papers for the Same Gene Are Handled
When multiple papers report findings on the **same gene**, Magna-Mined preserves every paper's record in SQLite (`extractions`) and resolves the gene's summary row in two ways:
* **Unclassified genes (`Not reported` in baseline)**:
  `reconcile` ranks all extractions for that gene by precedence:
  `human-reviewed` > **confidence tier** (`tier_1_high` > `tier_2_high` > `tier_3_medium` > `tier_4_low`) > **text source** (`XML` > `PDF` > `abstract_only`) > **method** (`CRISPR/knockout/complementation` > `overexpression/RNAi` > `biochemical/computational`) > **complementation tested** > **newer publication year**.
  * If multiple papers at the same top confidence tier report complementary non-conflicting roles (for example, one reports `virulence` and another reports `conidiation`), both processes are merged cleanly (e.g., `virulence, conidiation`).
  * If one paper reports that deleting the gene has **no defect** (`dispensable / wild-type-like`) while another reports a mutant defect, or if two high-tier papers report completely disjoint phenotypes, Magna-Mined sets **`conflict_status = conflict`** and queues the competing PMIDs and quotes for human curator review (`magna-mined reconcile --review-id <ID>`).
* **Already-classified baseline genes (`double-check` mode)**:
  Original baseline annotations in `baseline_process` are never overwritten; newly reported phenotypes from additional papers are appended to **`additional_functions`** alongside their PMIDs and years in **`additional_function_references`**.

### How Non-`MGG_xxxx` Gene Naming Styles Are Handled
* **Supported automatically**:
  * Formatting variants of `MGG` IDs (`MGG_12345`, `MGG-12345`, `MGG12345`, transcript suffixes such as `MGG_12345T0`, and known legacy Broad v6 IDs such as `MGG_00692` -> `MGG_09565` `PMK1`).
  * Standard gene symbols indexed from `MagnaGenes.csv` (e.g., `PMK1`, `ENG1`, `ZIP8/CREB4/ZIP6`), automatically matched with or without `Mo`/`Po` organism prefixes and mutant notation (`Δpmk1`, `delta-pmk1`, `pmk1-null`, `pmk1-KO`).
  * Newly coined gene names defined alongside an `MGG_` ID anywhere inside a paper (e.g., `"MoPct1 (MGG_01003)"`), which allows subsequent sentences using only `MoPct1` or `ΔMoPct1` to be linked accurately to `MGG_01003`.
* **Not mapped automatically (current limitation)**:
  * Completely different strain or genome locus conventions (for example `Mg01234`, `FR13_...`, `GUY11_...`, or GenBank/RefSeq accessions) **cannot** be resolved automatically unless the paper also mentions the corresponding `MGG_#####` ID or a known gene symbol from the MagnaGenes baseline. Additional external ID mappings can be added to SQLite `gene_aliases` or `LEGACY_MGG_MAP` in `entity_linker.py` when needed.

### Strengths (Pros) & Limitations

**Pros**
* **Zero API cost & offline extraction**: Runs out of the box with text-based pattern search, avoiding LLM API costs and generative hallucinations.
* **Full-text aware**: Reads open-access Europe PMC XML and PDF Results/Methods sections rather than relying solely on abstracts.
* **Resumable incremental updates**: Tracks every checked paper in `processed_papers`, allowing quarterly or semi-annual updates without re-reading old papers.
* **Multi-paper conflict detection & baseline protection**: Keeps baseline `MagnaGenes.csv` classifications intact, flags contradicting studies, and separates new discoveries (`newly_classified`) from secondary roles (`additional_functions`).
* **Visual exploration & direct MagnaGenes v2.0 export**: Generates self-contained interactive HTML dashboards (comparing database coverage prior to search vs. newly classified genes, phenotypes, tiers, and years) and dated `fin_data_YYYY-MM-DD.tsv` files ready for [`magnagenes_v2`](https://github.com/nehasahu486/magnagenes_v2) with SHA-256 deduplicated archiving.

**Limitations**
* **Paywalled papers**: Closed-access papers without open-access XML or PDF fall back to title/abstract only (`tier_4_low`), unless a local PDF is provided via `ingest-pdf`.
* **Figures & complex tables**: Phenotype data embedded only inside image figures or supplementary spreadsheets without supporting sentences in Results/Methods will not be captured by text extraction.
* **Fixed phenotype categories in text-based mode**: Text-based pattern search classifies across 5 main phenotype categories (`virulence`, `appressorium formation`, `conidiation`, `stress response`, `growth/development`, plus `dispensable / wild-type-like`) using curated keywords; atypical wording or specialized biochemical phenotypes may be missed or marked `not reported`.
* **Alternate locus schemes without `MGG_` IDs**: Papers using only non-MGG locus names (such as `Mg01234` or strain-specific IDs) without referencing an `MGG_#####` ID or known baseline symbol cannot be matched automatically.
* **Curator oversight on conflicts**: Contradictory findings across papers (`conflict_status = conflict`) and novel curations still benefit from brief curator inspection of the stored supporting quote before publishing updates to MagnaGenes v2.0.

---

## 3. Setup & Database Initialization

Creates a Python virtual environment and installs the `magna-mined` package and CLI entrypoint from `pyproject.toml`:

```bash
cd /path/to/magna-mined
python3 -m venv venv
source venv/bin/activate
pip install -e .
```

Seeds the SQLite database (`magnaporthe.db`) from the MagnaGenes v.1.0 dataset (Foster et al., 2021; required once):

```bash
# Seed from primary Excel dataset (recommended)
magna-mined init-db --excel data/Magnagenes.xlsx

# Or seed from CSV fallback
magna-mined init-db --csv data/MagnaGenes.csv

# Verify baseline counts (12,522 total genes, 10,886 unclassified)
magna-mined stats
```

---

## 4. Production Literature Mining Workflows

> All commands below assume the virtual environment is active (`source venv/bin/activate`, or invoked directly via `./venv/bin/magna-mined`).

### A. Quarterly or 6-Month Incremental Sync (Main Workflow)
Fetches newly published *M. oryzae* experimental papers from the specified window and curates **only previously unclassified (`Not reported`) genes**:

```bash
# Quarterly update (last 90 days)
magna-mined mine-papers --days 90 --max-papers 100

# Semi-annual update (last 180 days)
magna-mined mine-papers --days 180 --max-papers 200

# Incremental update since a specific date (YYYY-MM-DD)
magna-mined mine-papers --since 2026-01-01 --max-papers 200
```

### B. Year-by-Year Historical Literature Sync
Processes a multi-year publication window systematically year by year:

```bash
# Process up to 50 papers per year across 2022 and 2023
magna-mined mine-papers --year-from 2022 --year-to 2023 --max-papers 50
```

### C. Double-Checking Classified Genes for Additional Functions
Checks whether newly published papers report **additional biological roles** or knockout phenotypes for genes that are *already classified* in `MagnaGenes.csv`:

```bash
# Double-check papers from the last 90 days (or --days 180 / --year-from 2022 --year-to 2023)
magna-mined double-check --days 90 --max-papers 100
```
* Combined single-pass mode (unclassified curation + additional-function check):
  ```bash
  magna-mined mine-papers --days 180 --max-papers 150 --double-check
  ```
* When a baseline `virulence` gene is reported in a new paper to also regulate `conidiation`, the original `baseline_process` remains intact while `additional_functions` records `conidiation` and `additional_function_references` logs the PMID, publication year, and phenotype summary.

### D. Safe Re-Runs, Custom Queries & Paywalled PDFs
* **Resume without duplicates (`processed_papers`)**: Every evaluated `pmid` / `pmcid` and `checked_additional` state is stored in SQLite. Re-running a sync automatically skips papers already checked, avoiding redundant downloads and duplicate extractions.
* **Custom topic queries**: Overrides the default *M. oryzae* search query via `--query`:
  ```bash
  magna-mined mine-papers --query '("Magnaporthe oryzae") AND (appressorium OR effector) AND PUB_YEAR:2024' --max-papers 50
  ```
* **Single-gene lookup & local paywalled PDF reading**:
  ```bash
  # Curate one gene directly (supports MGG ID, symbol, Mo-prefix, Δmutant, or Broad v6 ID)
  magna-mined curate MoPMK1

  # Extract text from a local paywalled PDF using pymupdf4llm
  magna-mined ingest-pdf /path/to/paper.pdf --pmid 12345678
  ```

---

## 5. Ranking Evidence, Searching & Exporting Results

After a sync run, `reconcile` ranks competing evidence and flags conflicts, `search` / `show` inspect records, and `export` generates curated reports:

```bash
# 1. Rank evidence strength & flag conflicting studies for curator review
magna-mined reconcile

# 2. Filter newly classified, knockout-reported, or additional-function genes
magna-mined search --newly-classified
magna-mined search --knockout-only --has-extractions
magna-mined search --with-additional

# 3. Display a single gene's dossier, knockout summary & supporting quotes
magna-mined show MGG_01311

# 4. Export newly classified genes to Interactive HTML, Excel, CSV, TSV, or JSON
magna-mined export exports/newly_classified.html --newly-classified
magna-mined export exports/newly_classified.xlsx --newly-classified
magna-mined export exports/newly_classified.csv --newly-classified
```

### Key Export Columns
| Column | Description |
| :--- | :--- |
| `gene` / `gene_name` | Canonical MGG ID (`MGG_#####`) and symbol (`ENG1`, `MoPah1`, etc.) |
| `baseline_process` | Original classification from `MagnaGenes.csv` (`Not reported` if unclassified) |
| `process_affected` | Curated affected biological process(es) (`virulence`, `conidiation`, etc.) |
| `newly_classified` | `1` if previously `Not reported` in baseline and newly curated by Magna-Mined |
| **`knockout_reported`** | `1` / `Yes` if a targeted knockout / deletion / CRISPR mutant was reported |
| **`knockout_phenotype`** | Confidence tier, affected process, phenotype summary, and PMID for the knockout |
| **`additional_functions`** | New biological processes discovered beyond `baseline_process` |
| **`additional_function_references`** | Citations (`PMID`, year, summary) supporting each additional function |
| `consensus_tier` | `tier_1_high` (deletion + complementation) down to `tier_4_low` |
| `conflict_status` | `none`, `conflict` (queued for curator review), or `resolved` |

### Conflict Review & Ground-Truth Benchmarking
```bash
# Mark a specific extraction ID as verified by the curator to resolve a conflict
magna-mined reconcile --review-id 42 --notes "Verified in Guy11 background"

# Evaluate Precision, Recall, F1 & quote accuracy against MagnaGenes v.1.0 (Foster et al., 2021) ground truth
magna-mined eval-benchmark --mode db
magna-mined eval-benchmark --mode live --sample-size 6
```

---

## 6. Syncing with MagnaGenes v2.0

[MagnaGenes v2.0](https://github.com/nehasahu486/magnagenes_v2) ([live app](https://nehasahu486.github.io/magnagenes_v2/)) builds its dataset from a long-format `(Gene, variable, value)` TSV file (`data-raw/source/annotation/fin_data.tsv`) via `data-raw/prepare_data.R`.

`magna-mined export-magnagenes` generates a dated drop-in file (`fin_data_YYYY-MM-DD.tsv`), is used to merge newly curated annotations into the existing `fin_data*.tsv`, and automatically moves older files into `archive/` with **SHA-256 deduplication** (preventing duplicate copies from accumulating in `archive/` and skipping new writes if the file content has not changed):

```bash
# Option A: Export directly into a local magnagenes_v2 folder and merge into the latest existing file
magna-mined export-magnagenes \
  --output-dir /path/to/magnagenes_v2/data-raw/source/annotation \
  --merge-existing

# Option B: Export to exports/magnagenes/fin_data_YYYY-MM-DD.tsv merged with a base fin_data.tsv
magna-mined export-magnagenes \
  --output-dir exports/magnagenes \
  --base-fin-data /path/to/fin_data.tsv
```

When [`scripts/resolve_latest_fin_data.R`](scripts/resolve_latest_fin_data.R) is sourced inside `magnagenes_v2/data-raw/prepare_data.R`, rebuilding `app/data/genes.rds` automatically selects the newest `fin_data_YYYY-MM-DD.tsv` and archives older files into `archive/` without duplicates.

---

## 7. Project Structure & Tests

```
magna-mined/
├── README.md                      # Documentation & production CLI workflows
├── pyproject.toml                 # Package metadata, dependencies & CLI entrypoint
├── data/                          # MagnaGenes v.1.0 dataset (Foster et al., 2021)
│   ├── Magnagenes.xlsx
│   └── MagnaGenes.csv
├── exports/                       # Generated reports & dated MagnaGenes sync outputs
│   └── magnagenes/
│       └── archive/               # Deduplicated archive of prior fin_data files
├── scripts/
│   └── resolve_latest_fin_data.R  # R helper for magnagenes_v2/data-raw/prepare_data.R
├── src/
│   ├── img/                       # Project logo & icon assets
│   └── magnamined/
│       ├── __init__.py
│       ├── schema.py              # Data models & controlled phenotype vocabularies
│       ├── database.py            # SQLite schema, migrations, KO & additional-function tracking
│       ├── entity_linker.py       # Gene symbol/mutant matching & inline alias detection
│       ├── literature.py          # Europe PMC search, XML parsing, Unpaywall & PDF extraction
│       ├── extractor.py           # Abstract screening, text-based rule search & optional Gemini LLM
│       ├── paper_miner.py         # Bulk paper screening, multi-gene extraction & double-check engine
│       ├── reconciler.py          # Evidence ranking, conflict detection & curator review queue
│       ├── benchmark.py           # MagnaGenes v.1.0 ground-truth Precision/Recall/F1 evaluator
│       ├── explorer.py            # Search, dossier & CSV/TSV/XLSX/JSON/HTML exports
│       ├── magnagenes_sync.py     # Dated fin_data_YYYY-MM-DD.tsv exporter & deduplicated archiver
│       └── pipeline.py            # Click CLI entrypoint (`magna-mined`)
└── tests/                         # Pytest suite (54 unit & integration tests)
```

Executes the test suite:
```bash
./venv/bin/pytest -v
```

### Optional Environment Variables
* `GEMINI_API_KEY`: Optional. Enables Google Gemini structured extraction instead of the default zero-dependency text-based rule search.
* `MAGNAMINED_CONTACT_EMAIL`: Optional. Enables DOI-based open-access PDF lookup via the Unpaywall API when Europe PMC lacks full-text XML.

Example terminal setup (macOS / Linux):
```bash
# 1. Enable Unpaywall open-access PDF lookup for papers without PMC XML
export MAGNAMINED_CONTACT_EMAIL="curator@example.org"

# 2. (Optional) Enable Gemini LLM extraction instead of built-in text-based search
export GEMINI_API_KEY="AIzaSy...example_api_key..."

# 3. Run quarterly mining with those environment variables active in the current shell
magna-mined mine-papers --days 90 --max-papers 100
```

Or pass them for a single command without storing them in the shell session:
```bash
MAGNAMINED_CONTACT_EMAIL="curator@example.org" magna-mined mine-papers --days 90 --max-papers 100
```
