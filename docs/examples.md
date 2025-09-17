# Examples

## Toy bulk example

- Inputs:
  - `data/toy_bulk_counts.tsv` (3 genes × 2 samples)
  - `genes/toy_pathway.csv` (schema: `gene_symbol,gene_name,module,annotation,evidence_source,weight`)
- Run:
  - `make toy-bulk`
- Outputs (in `out/toy-bulk/`):
  - `scores_bulk.csv`, `scores_bulk_wide.csv`
  - `scores_bulk_by_treatment.csv`, `scores_bulk_by_treatment_wide.csv`
  - `stats_bulk.csv`, `report_bulk.pdf`

## Single-cell (Zhang 2024) with bundled sample info

- Inputs: `.h5` files under `data/GSE272085/` (downloaded by `make sc-zhang`)
- Sample info: `data/zhang_sample_info.csv` (maps sample → Control/Cachexia)
- Run:
  - `make sc-zhang-with-info`
- Outputs (in `out/zhang2024_with_info/`):
  - `scores_sc_by_celltype.csv`, `scores_sc_cells.csv`, `deltas.csv`
  - `scores_by_celltype_condition_tidy.csv`, `deltas_by_celltype.csv`
  - `report.pdf`, `comparisons_by_celltype_condition.pdf`, `comparisons_deltas.pdf`

## Pathway CSV schema

Exactly these columns:

```
gene_symbol,gene_name,module,annotation,evidence_source,weight
```

See `genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv` for an example.
