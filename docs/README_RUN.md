# Zhang 2024 (Cell Reports) — BDNF-axis enrichment runner

This adds a **headless** script so you can run the analysis without Jupyter, using your existing `scanpy_env.yml`.

## Quick start
```bash
# 1) Ensure conda is on your PATH (miniconda/anaconda)
conda --version

# 2) Put these files (in the project root)
#    - scanpy_env.yml
#    - run_bdnf_zhang2024.sh
#    - run_bdnf_enrichment.py

# 3) Run the driver (downloads GEO, extracts, runs analysis)
bash run_bdnf_zhang2024.sh
```

Outputs will appear in `./out/zhang2024/`:
- `bdnf_axis_scores_by_celltype.csv`
- `bdnf_axis_scores_by_condition_celltype.csv`
- `axis_genes_mean_expression.csv`
- `zhang2024_bdnf_axis.h5ad`
- `figures/umap_zhang2024_meta.png` and `figures/umap_zhang2024_scores.png` (under Scanpy's `figures/` dir)

## Notes
- The script infers **Control** vs **KIC** from file names (case-insensitive: `control`, `KIC`, `cachexia`, etc.).
- If you previously hit `invalid non-printable character U+00A0`, that is a non-breaking space copied into code.
  Replace it with a normal space or retype the line. The headless script here avoids that copy/paste pitfall.
- Gene panels used:
  - proBDNF axis: `Ngfr`, `Sort1`, `Sorcs2`, `Sorcs3`
  - TrkB axis: `Ntrk2`
  - Catabolic readouts: `Trim63` (MuRF1), `Fbxo32` (Atrogin-1)

## New features (Aug 2025)
- **Weighted pathway scores**: If you place this CSV next to the scripts, it will be auto‑detected and used:
    - `mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv`
  It should follow the six‑column schema `gene_symbol,gene_name,module,annotation,evidence_source,weight`. We match symbols to the dataset, keep those present, normalize weights to sum=1 over detected genes, and compute a weighted mean of log1p expression per cell.
  Outputs include `score_TrkB_w` in the summary CSVs and PDF.

- **Selected-gene violin plots**: By default we plot `Ngfr,Bdnf,Ntrk2,Rela` grouped by `celltype_guess` to
  `figures/umap_zhang2024_violin_selected_by_celltype.png`. Override with `--violin_genes`.

### CLI flags
You can pass custom CSVs or genes:
```bash
python run_bdnf_enrichment.py --input_dir ./data/GSE272085 --out_dir ./out/zhang2024 \
  --trkb_csv ./mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv \
  --violin_genes "Ngfr,Bdnf,Ntrk2,Rela"
```
The `run_bdnf_zhang2024.sh` runner will pass these automatically if the files exist in the working directory.
