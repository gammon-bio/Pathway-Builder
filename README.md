**Pathway Builder**

[![PyPI](https://img.shields.io/pypi/v/PACKAGE_NAME.svg)](https://pypi.org/project/PACKAGE_NAME/)
[![Python Versions](https://img.shields.io/pypi/pyversions/PACKAGE_NAME.svg)](https://pypi.org/project/PACKAGE_NAME/)
[![Codecov](https://img.shields.io/codecov/c/github/gammon-bio/Pathway-Builder?logo=codecov)](https://app.codecov.io/gh/gammon-bio/Pathway-Builder)
[![CI](https://github.com/gammon-bio/Pathway-Builder/actions/workflows/ci.yml/badge.svg)](https://github.com/gammon-bio/Pathway-Builder/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

- Purpose: Compute weighted pathway scores for bulk RNA-seq or single-cell/nucleus RNA-seq using a simple, file-driven workflow or a Python API.
- Inputs: One or more pathway CSV files with columns `gene` and `weight` (case-insensitive). For bulk, default scoring is R-style (per-gene z-score across samples, then a weighted mean with optional KEGG/Reactome evidence boost).
- Outputs: CSVs with weighted scores (per sample for bulk; per celltype/condition for single-cell). Optional figures when using the reference scripts.

> **Quick Start**
>
> ```bash
> pip install -e '.[dev]'   # quote extras on zsh to avoid globbing
> make toy-bulk
> ```
>
> Inspect `out/toy-bulk/` for the demo scores and PDF report.

**Getting Started**

- Install development extras: `pip install -e '.[dev]'` (quote extras in zsh)
- One‑liner (toy bulk demo and validation):
  - `make toy-bulk && make check OUT=out/toy-bulk`
- Examples and datasets: see `docs/examples.md`.

**Setup**

- Use the Makefile from inside `pathway-builder-starter/`:
  - Bulk/report extras: `make setup`
  - Single‑cell + report extras: `make setup-sc`
- Optional (advanced): direct install without Makefile
  - Bulk/report only: `python3 -m pip install -e '.[report]'`
  - Single‑cell + bulk/report: `python3 -m pip install -e '.[singlecell,report]'`
- Python version: 3.9–3.12

**Quick Start**

- Bulk (genes x samples table)
  - Install extras: `pip install -e '.[dev]'`
  - Toy example: `make toy-bulk` (writes to `out/toy-bulk`; validates required deps)
  - Validate outputs: `make check OUT=out/toy-bulk`
  - Quick lint/type sweep: `black --check . && isort --check . && flake8 . && mypy .`

**Cheat Sheet**

- Setup bulk/report extras: `make setup`
- Setup single‑cell extras: `make setup-sc`
- Run toy bulk: `make toy-bulk`
- Validate any bulk run: `make check OUT=out/my_run`
- Run your bulk: `make bulk COUNTS=/path/to/counts.tsv PATHWAYS='genes/pw_a.csv genes/pw_b.csv' LABELS='A B' OUT=out/bulk`
- Run single-cell: `make sc SN_DATA='data/sc/*.h5' PATHWAYS='genes/pw_a.csv genes/pw_b.csv' LABELS='A B' OUT=out/sc`
- Single‑cell with sample info: `make sc SN_DATA='data/sc/*.h5' PATHWAYS='genes/pw_a.csv' SAMPLE_INFO=data/zhang_sample_info.csv OUT=out/sc`
- Zhang (download + run): `make sc-zhang`
- Zhang with sample info: `make sc-zhang-with-info`
- Validate any SC run: `make sc-check OUT=out/sc`

`make toy-bulk` now runs in "simple" scoring mode with `--no_pdf` so the smoke test avoids heavy plotting dependencies. Use `make bulk ... SCORING_STYLE=r REPORT_STYLE=box` for the full PDF workflow.

**Toy Bulk Example**

- Inputs included:
  - `data/toy_bulk_counts.tsv` (3 genes x 2 samples)
  - `genes/toy_pathway.csv` (weights for GeneA and GeneC)
  - `data/toy_sample_info.csv` (maps Sample1→Ctrl, Sample2→Case)
- Run via Makefile: `make toy-bulk`
- Expected results:
  - Output `out/toy-bulk/scores_bulk.csv` contains two rows (Sample1=20.0, Sample2=5.0), from 0.5*GeneA + 0.5*GeneC per sample.
  - If sample info used: `scores_bulk_with_treatment.csv`, `scores_bulk_by_treatment.csv`, `scores_bulk_by_treatment_wide.csv`, `stats_bulk.csv`, and `report_bulk.pdf`.

 

**Sample Info Sheet**

- CSV with at least two columns: `sample` and `treatment` (you can rename with `--sample_col/--treatment_col`).
- Example:

sample,treatment
S1,Control
S2,Control
S3,Case
S4,Case

- The CLI will join treatments onto sample-level scores, aggregate by treatment, run Welch's t-test when there are exactly two treatments, and write a PDF report.

**Statistical Tests (Bulk)**

- Exactly two treatments: Welch's t-test (two-sided) on pathway scores (case minus control). Control is auto-detected when labels include any of `control|ctrl|vehicle|veh` (case-insensitive). Override with `--control_label`.
- More than two treatments: One-way ANOVA across treatments. Report includes F and p; BH q-values via `statsmodels` (installed with the `report` extra).
- Post-hoc for ANOVA: Tukey HSD pairwise comparisons are written to `stats_bulk_tukey.csv` with columns: `group1, group2, meandiff, p_adj, lower, upper, reject`.
- The PDF annotates each pathway plot with the relevant p (and q where available). A `stats_bulk.csv` file is also written alongside the PDF.

**Run Your Own Bulk Dataset**

- Place your files anywhere; the Makefile/CLI takes explicit paths.
- Required inputs:
  - Counts table (genes x samples; TPM/CPM/VST): `COUNTS=/path/to/counts.tsv`
  - One or more pathway CSVs with `gene,weight` (+ optional `evidence_source`): `PATHWAYS='genes/pw1.csv genes/pw2.csv'`
  - Optional sample info to name groups for stats: `SAMPLE_INFO=/path/to/sample_info.csv` with columns `sample,treatment` (rename via `--sample_col/--treatment_col`)
- Use the Makefile:
  - `make bulk \
      COUNTS=/path/to/counts.tsv \
      PATHWAYS='genes/pw_trkb.csv genes/pw_p75.csv' \
      LABELS='TRKB P75' \
      SAMPLE_INFO=/path/to/sample_info.csv \
      CONTROL=Control \
      SCORING_STYLE=r COLLAPSE=mean REPORT_STYLE=box \
      OUT=/path/to/out_dir`
- Outputs: PDF report with group plots + CSVs (`scores_bulk*.csv`, `stats_bulk.csv`, and Tukey post-hoc if >2 groups).
  - Validate outputs: `make check OUT=/path/to/out_dir`

**Run Single-Cell**

- First, install single‑cell dependencies: `make setup-sc`
- Your dataset (point to `.h5ad`, 10x `.h5`, a 10x matrix directory, or a glob — multiple inputs are concatenated):
  - `make sc \
      SN_DATA='data/sc/*.h5' \
      PATHWAYS='genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv' \
      LABELS='TRKB' \
      CELLTYPE=celltype_guess \
      CONDITION=condition \
      OUT=out/sc_run`
  - Add a sample info sheet (to set/override conditions) with columns `sample,condition` and pass it:
    - `make sc \
        SN_DATA='data/sc/*.h5' \
        PATHWAYS='genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv' \
        LABELS='TRKB' \
        SAMPLE_INFO=data/zhang_sample_info.csv \
        OUT=out/sc_run`
- Using your BDNF panels (example):
  - `make sc \
      SN_DATA='data/GSE272085/*.h5' \
      PATHWAYS='genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv' \
      LABELS='TRKB' \
      SAMPLE_INFO=data/zhang_sample_info.csv \
      OUT=out/zhang2024_both`
- Outputs:
  - `scores_sc_by_celltype.csv` and `deltas.csv` (case-control deltas when two conditions)
  - `scores_sc_cells.csv` (cell-level scores for each pathway)
  - `scores_by_celltype_condition_tidy.csv` (pathway, celltype, condition, mean, n)
  - `deltas_by_celltype.csv` (delta per celltype for each pathway)
  - Figures in `figs/` and `report.pdf` (line plots by celltype and boxplots by condition)
  - Additional PDFs: `comparisons_by_celltype_condition.pdf` (grouped bars) and `comparisons_deltas.pdf` (delta bars)
- Validate outputs:
  - `make sc-check OUT=out/sc_run`
- Zhang 2024 dataset (auto-download via script):
  - `make sc-zhang` (requires network + conda as per `scripts/run_bdnf_zhang2024.sh`)
  - Use bundled sample info directly: `make sc-zhang-with-info` (writes to `out/zhang2024_with_info`)

**Pathway CSV Format**

- Required columns: `gene`, `weight` (case-insensitive). If `weight` is missing, all weights default to 1.0.
- Duplicate genes are summed, then re-normalized across genes found in your dataset.
- Gene matching is case-insensitive. Missing genes are ignored (weights renormalize over those present).

Example `genes/my_pathway.csv`:

gene,weight
NTRK2,1.0
NGFR,0.5
SORT1,0.5

**What the score means**

- Bulk (R-style, default): For each panel (e.g., TrkB, p75), z-score each gene across samples, then compute a weighted mean per sample.
  - `Z_i(sample) = (expr_i(sample) - mean(expr_i)) / sd(expr_i)`
  - `score(sample) = Σ_i (w_i' / Σ_i w_i') * Z_i(sample)`
  - Evidence boost (optional): `w_i' = min( w_i + 0.05·[has KEGG] + 0.05·[has Reactome], 1.2 )`
  - Duplicate symbols are collapsed before scoring (default: mean across duplicates).
- Bulk (simple): `score(sample) = Σ_i (w_i / Σ_i w_i) * expr_i(sample)` (no z-score, no boost).
- Single-cell: `score(cell) = Σ_i w_i * expr_i`, where `expr_i` is library-size normalized (`normalize_total(target_sum=1e4)`) and `log1p` transformed.

**Python API**

- Bulk:
  - `from pathway_builder.core import score_bulk_from_table`
  - `df_scores = score_bulk_from_table(df_counts, gene_col=None, pathway_csvs=["genes/my_pathway.csv"], labels=["MY_PW"])`

```python
import pandas as pd
from pathway_builder.core import score_bulk_from_table
counts = pd.DataFrame({"gene": ["NTRK2"], "Sample1": [1.0], "Sample2": [2.0]})
scores = score_bulk_from_table(counts, gene_col="gene", pathway_csvs=["genes/toy_pathway.csv"])
print(scores.head())
```

- Single-cell:
  - `from pathway_builder.core import load_singlecell, score_singlecell_adata`
  - `adata = load_singlecell("data/sc/*.h5")`
  - `summary = score_singlecell_adata(adata, pathway_csvs=["genes/my_pathway.csv"], labels=["MY_PW"], celltype_col="celltype", condition_col="condition")`

**Single-Cell: Pointing to Any Dataset**

- Inputs: pass a single path or a quoted glob (recommended to avoid shell expansion):
  - Single `.h5ad`, 10x `.h5`, or a 10x matrix directory
  - Multiple files via glob: `SN_DATA='data/*.h5'` (quotes ensure the CLI expands it)
- What happens:
  - Inputs are concatenated (outer join on genes). A `sample` column is derived from the base filename (or directory name) and preserved.
  - Duplicate gene symbols are collapsed by summing columns (fast sparse mapping) before filtering/scoring.
  - Genes are filtered (`min_cells>=5`), library-size normalized, and log1p transformed before scoring.
  - If `condition` is missing, all cells are labeled `All`. To set/override conditions, provide a sample info sheet (CSV with `sample,condition`).
    - Makefile: add `SAMPLE_INFO=/path/to/sample_info.csv` to the `sc` target.
    - CLI: pass `--sample_info`, with `--sample_col` and `--info_condition_col` if names differ from defaults.
- Outputs (key files):
  - `scores_sc_by_celltype.csv`, `scores_sc_cells.csv`
  - `deltas.csv` (case − control across conditions when exactly two conditions)
  - Comparisons: `scores_by_celltype_condition_tidy.csv`, `deltas_by_celltype.csv`
  - Figures: `report.pdf` (includes a delta page when two conditions), `comparisons_by_celltype_condition.pdf`, `comparisons_deltas.pdf`, and PNGs in `figs/`
- Use the Makefile `sc` target as shown in “Run Single‑Cell” to pass `SN_DATA`, `PATHWAYS`, `LABELS`, `OUT`, and optional `SAMPLE_INFO`.

**Reference Scripts**

- Generic CLI (existing script): `python run_bdnf_enrichment.py --help`
- Zhang 2024 headless runner: `bash scripts/run_bdnf_zhang2024.sh` (see `docs/README_RUN.md`).
  - To point that pipeline at a different dataset, prefer the Makefile `sc` target with your `SN_DATA` glob as shown above, or run `python run_bdnf_enrichment.py --input_dir <dir> --out_dir <out>` directly.

**Reproducible LLM Prompt (to generate a weighted gene list)**

- Use the prompt below verbatim, then provide context about your pathway and species. The model must output a CSV with the same columns used in our curated lists (e.g., `genes/mouse_probdnf_p75_sortilin_gene_list_full_annotated_sources.csv`).

---
You are assisting in building a weighted gene panel representing a signaling pathway for downstream scoring on RNA‑seq (bulk and single‑cell). Produce a CSV with exactly these six columns and a header (no extra columns, no prose):

gene_symbol,gene_name,module,annotation,evidence_source,weight

Instructions:
- gene_symbol: Official gene symbol for the specified species (mouse or human), one per row.
- gene_name: Descriptive full name for the gene.
- module: High‑level role within the pathway (e.g., ligand, receptor, co‑receptor, adaptor, kinase, transcription factor).
- annotation: One short phrase capturing the specific role (e.g., “core: p75NTR receptor”, “extended adaptor: NF‑κB arm”).
- evidence_source: Curated sources that justify inclusion (e.g., “KEGG:mmu04722; Reactome: p75NTR signaling; peer‑reviewed reviews”).
- weight: Real number > 0 encoding relative importance; higher means more influence. Suggested scale: 0.3–1.0 (core ≈1.0; extended ≈0.3–0.7). Keep weights within [0, 1.2].

Coverage guidance:
- Include core ligand(s), receptor(s)/co‑receptors, critical adaptors/kinases, and key transcriptional effectors.
- Prefer 10–60 genes total. If uncertain, include with a lower weight rather than omitting.

Output only CSV rows following the header above.

Provide the panel for: [describe pathway], species: [human|mouse].
---

Tips:
- Start with a “core” set at weight 1.0; add extended modulators at 0.3–0.7.
- If using a mouse CSV on human data (or vice versa), map symbols via orthologs before use.

**CI**

- GitHub Actions workflow installs the package and runs tests (`.github/workflows/ci.yml`).
- Tests cover bulk and a minimal single-cell case.

**Environment**

- For single-cell workflows, install extras: `pip install -e .[singlecell]` or create the provided conda env: `env/scanpy_env.yml`.

**Outputs**

- Bulk: `scores_bulk.csv` with columns: `pathway, sample, score_w`.
- Single-cell: `scores_sc_by_celltype.csv` with columns: `celltype, condition, score_w__<LABEL>`.

**Common Errors**

- Gene column not found: The tool auto-detects `gene`, `symbol`, `gene_symbol`, etc. If your counts file uses a different name (or the first column isn’t genes), pass `--gene_col <name>`.
- Unmatched `--sample_info` samples: If some samples in `scores_bulk.csv` aren’t present in your sample info sheet (or vice versa), you’ll see a warning. Ensure the `sample` IDs in `--sample_info` exactly match the column names in your counts table (after gene column removal).
- No files matched for `--sn_data`: Quote globs to avoid shell expansion (`SN_DATA='data/sc/*.h5'`). Supported inputs: 10x `.h5`, 10x matrix dir, or `.h5ad`.
- Pathway CSV schema: Requires `gene` (or `gene_symbol`) and `weight`. If `weight` is missing, it defaults to `1.0`. Gene matching is case-insensitive; missing genes are ignored and weights re-normalize.
- No genes matched: If zero pathway genes match your dataset, scores are `0.0`. Check species/symbols and case. For bulk, ensure your counts table gene symbols match the pathway list.
