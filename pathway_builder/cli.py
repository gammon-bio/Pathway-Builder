from __future__ import annotations

import argparse
import os
from typing import Optional, Sequence

from .core import (
    _is_tsv_or_csv,
    _read_delim,
    load_singlecell,
    score_bulk_from_table,
    score_singlecell_adata,
)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Score weighted gene pathways for bulk RNA-seq or single-cell data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = ap.add_argument_group("Inputs")
    src.add_argument("--bulk_counts", help="Bulk normalized counts table (genes x samples; CSV/TSV)")
    src.add_argument("--sn_data", help="Single-cell input: .h5ad, 10x .h5, 10x mtx dir, or a glob pattern")
    src.add_argument("--pathway_csv", action="append", help="Pathway CSV(s) with columns [gene,weight]; can be repeated")
    src.add_argument("--label", action="append", help="Label(s) for each pathway CSV; can be repeated; defaults to filename")

    scg = ap.add_argument_group("Single-cell options")
    scg.add_argument("--celltype_col", help="obs column for cell-type")
    scg.add_argument("--condition_col", help="obs column for condition/case-control grouping")
    scg.add_argument("--no_pdf", action="store_true", help="Do not produce a PDF report (bulk mode)")

    bg = ap.add_argument_group("Bulk options")
    bg.add_argument("--gene_col", help="Gene symbol column name in bulk table (auto-detected)")
    bg.add_argument("--sample_info", help="CSV mapping samples to treatments (columns: sample,treatment)")
    bg.add_argument("--sample_col", default="sample", help="Column name for sample ID in --sample_info")
    bg.add_argument("--treatment_col", default="treatment", help="Column name for treatment/group in --sample_info")
    bg.add_argument("--pdf_report", help="Optional path for PDF report (default: <output_dir>/report_bulk.pdf)")
    bg.add_argument("--control_label", help="Optional explicit control label (used for Welch t-test when two treatments)")
    bg.add_argument("--scoring_style", choices=["simple", "r"], default="r", help="Scoring formula: simple weighted mean vs R-style (z-score then weighted mean with evidence boost)")
    bg.add_argument("--collapse_duplicates", choices=["mean", "median", "sum", "maxvar"], default="mean", help="How to collapse duplicate genes in bulk counts")
    bg.add_argument("--no_boost", action="store_true", help="Disable evidence-based weight boosting (R-style only)")
    bg.add_argument("--report_style", choices=["bar", "box"], default="box", help="Plot style for PDF (bulk)")

    out = ap.add_argument_group("Outputs")
    out.add_argument("--output_dir", required=True, help="Directory to write outputs")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = build_parser()
    args = ap.parse_args(argv)

    has_bulk = bool(args.bulk_counts)
    has_sc = bool(args.sn_data)
    if has_bulk == has_sc:
        raise SystemExit("Specify exactly one of --bulk_counts or --sn_data")
    if not args.pathway_csv:
        raise SystemExit("At least one --pathway_csv is required")

    os.makedirs(args.output_dir, exist_ok=True)

    if has_bulk:
        if not _is_tsv_or_csv(args.bulk_counts):
            raise SystemExit("--bulk_counts must be .csv, .tsv, or .txt (genes x samples)")
        df = _read_delim(args.bulk_counts)
        if args.scoring_style == "r":
            from .core import score_bulk_r_style_from_table
            out = score_bulk_r_style_from_table(
                df,
                gene_col=args.gene_col,
                pathway_csvs=list(args.pathway_csv),
                labels=list(args.label) if args.label else None,
                collapse=args.collapse_duplicates,
                boost_by_evidence=not args.no_boost,
            )
        else:
            out = score_bulk_from_table(df, gene_col=args.gene_col, pathway_csvs=list(args.pathway_csv), labels=list(args.label) if args.label else None)
        # If sample info provided, join and produce summaries + PDF
        scores = out.copy()
        if args.sample_info:
            import pandas as pd
            si = pd.read_csv(args.sample_info)
            if args.sample_col not in si.columns or args.treatment_col not in si.columns:
                raise SystemExit(f"--sample_info must contain columns: {args.sample_col},{args.treatment_col}")
            si2 = si[[args.sample_col, args.treatment_col]].rename(columns={args.sample_col: "sample", args.treatment_col: "treatment"})
            scores = scores.merge(si2, on="sample", how="left")
            scores.to_csv(os.path.join(args.output_dir, "scores_bulk_with_treatment.csv"), index=False)
            # Generate PDF unless disabled
            if not args.no_pdf:
                from .report import make_bulk_pdf_report, _infer_control_label
                pdf_path = args.pdf_report or os.path.join(args.output_dir, "report_bulk.pdf")
                methods_note = (
                    ("R-style scoring: per-gene z-score across samples, then weighted mean per sample; weights optionally boosted by evidence (KEGG/Reactome).\n" if args.scoring_style == "r" else
                     "Bulk weighted scores: score = sum_i w_i * expr_i. Weights normalized over detected genes.\n") +
                    "Samples matched to treatments via the provided sample info sheet. Welch's t-test when two treatments; one-way ANOVA (plus Tukey) when more than two."
                )
                # Control detection when needed
                ctrl = args.control_label
                if ctrl is None:
                    levels = list(scores["treatment"].dropna().unique())
                    ctrl = _infer_control_label(levels)
                make_bulk_pdf_report(scores_with_treatment=scores, output_pdf_path=pdf_path, methods_note=methods_note, control_label=ctrl, style=args.report_style)
        out.to_csv(os.path.join(args.output_dir, "scores_bulk.csv"), index=False)
    else:
        adata = load_singlecell(args.sn_data)
        summary = score_singlecell_adata(
            adata=adata,
            pathway_csvs=list(args.pathway_csv),
            labels=list(args.label) if args.label else None,
            celltype_col=args.celltype_col,
            condition_col=args.condition_col,
        )
        summary.to_csv(os.path.join(args.output_dir, "scores_sc_by_celltype.csv"), index=False)


if __name__ == "__main__":
    main()
