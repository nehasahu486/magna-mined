# Project Context & System Specification: Magna-Mined

## 1. Executive Summary & Biological Mission

**Magna-Mined** is an automated literature mining and gene curation tool for ***Magnaporthe oryzae*** (the rice blast fungus, also known as *Pyricularia oryzae*). Rice blast is one of the most destructive fungal diseases of cultivated rice worldwide. Maintaining an accurate catalog of experimentally studied genes supports research across all major infection and developmental processes:

- **Appressorium formation**: Specialized dome-shaped infection structures that generate turgor pressure to penetrate the rice leaf cuticle.
- **Virulence / Pathogenicity**: Colonization of host tissues, evasion of plant immunity, and necrotic lesion development.
- **Conidiation / Sporulation**: Spore development, morphology, and dispersal driving secondary infection cycles.
- **Vegetative growth & Stress response**: Osmotic/oxidative stress adaptation, cell wall integrity, and melanin biosynthesis.

### Baseline Dataset (Zenodo Record 4647766)
- **Source**: [Zenodo Record 4647766: MagnaGenes Dataset](https://zenodo.org/records/4647766)
- **Primary Data Files**: `data/Magnagenes.xlsx` and `data/MagnaGenes.csv`
- **Coverage**: **12,522 genes** across the *M. oryzae* proteome:
  - **~1,636 baseline-classified genes** (1,274 with documented mutant phenotypes + 362 tested as `none`)
  - **10,886 unclassified genes** marked `Not reported`

Magna-Mined automates continuous discovery, quote-verified evidence extraction, knockout tracking, additional-function detection, and synchronization with [MagnaGenes v2.0](https://nehasahu486.github.io/magnagenes_v2/).

---

## 2. Genomic Identifiers & Nomenclature

Literature mining accounts for standard fungal genetics nomenclature:
1. **Systematic Identifiers (MGG IDs)**: `MGG_\d{5}` (e.g., `MGG_01001`, `MGG_00692`), including transcript suffixes (`MGG_01001T0`) and Broad v6 aliases (`MGG_#####.6`).
2. **Gene Names & Organism Prefixes**: 3–4 letter symbols (`ENG1`, `PMK1`, `MAC1`) and *M. oryzae*-prefixed variants (`MoENG1`, `MoPMK1`), plus inline paper definitions such as `MoPct1 (MGG_01003)`.
3. **Mutant Designations**: `Δeng1`, `∆mgg_01001`, `eng1-null`, `pmk1Δ`, CRISPR knockouts.
4. **Reference Strains**: Primarily **70-15** and **Guy11**.

---

## 3. Paper-First Workflow & Full-Text Retrieval

Instead of issuing 12,522 separate gene-by-gene queries, Magna-Mined streams *M. oryzae* publications chronologically (quarterly, semi-annually, or by year range) through:
1. **Processed-Paper Tracking (`processed_papers`)**: Skips already-evaluated PMIDs/PMCIDs on subsequent runs.
2. **Abstract Screening**: Filters out reviews and non-experimental surveys before downloading full text.
3. **Gene Name & Alias Matching (`entity_linker.py`)**: Resolves MGG IDs, symbols, mutant notation, and inline aliases.
4. **Full-Text Retrieval Priority (`literature.py`)**:
   - **Priority 1**: Europe PMC structured full-text XML (`abstract`, `results`, and `methods` sections).
   - **Priority 2**: PDF-to-Markdown extraction via **`pymupdf4llm`** (Unpaywall open-access PDFs or local paywalled PDFs).
   - **Priority 3**: Abstract-only fallback (capped at `tier_4_low`).
5. **Quote-Verified Extraction (`extractor.py`)**: Zero-dependency sentence-level rule engine by default (or optional Gemini API), requiring an exact `verbatim_quote` verified in the source document and filtering out control/qPCR marker mentions.
6. **Evidence Ranking & MagnaGenes Sync (`reconciler.py`, `magnagenes_sync.py`)**: Ranks tiers, flags conflicts, tracks `knockout_reported` / `knockout_phenotype` and `additional_functions`, and exports dated `fin_data_YYYY-MM-DD.tsv` updates with SHA-256 deduplicated archiving.

---

## 4. Evidence Confidence Tiers & Core Tables

- **Tier 1 (`tier_1_high`)**: Targeted gene deletion / CRISPR knockout + genetic complementation rescue.
- **Tier 2 (`tier_2_high`)**: Targeted knockout / deletion mutant phenotype without complementation.
- **Tier 3 (`tier_3_medium`)**: Overexpression, RNAi silencing, or site-directed/point mutagenesis.
- **Tier 4 (`tier_4_low`)**: Abstract-only or expression/contextual evidence.

### SQLite Tables (`magnaporthe.db`)
- **`genes`**: Master catalog (`gene`, `gene_name`, `baseline_process`, `process_affected`, `detailed_mutant_phenotype`, `knockout_reported`, `knockout_phenotype`, `additional_functions`, `additional_function_references`, `consensus_tier`, `conflict_status`, `pmid_or_link_to_publication`, `blast2go_annotation`, `pfam`, etc.).
- **`gene_aliases`**: Canonical MGG-to-symbol/Broad-v6 mappings.
- **`extractions`**: Verified per-paper evidence (`gene_id`, `pmid`, `pmcid`, `doi`, `publication_year`, `process_affected`, `experimental_method`, `complementation_tested`, `evidence_source`, `confidence_tier`, `verbatim_quote`, `reviewed_by_human`, `review_notes`).
- **`processed_papers`**: Processed-paper ledger (`pmid`, `pmcid`, `checked_additional`, `processed_at`).
- **`clusters` / `cluster_assignments`**: Baseline Zenodo functional cluster mappings.
