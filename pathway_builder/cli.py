from __future__ import annotations

import argparse
import os
from typing import Optional, Sequence
import sys

from .core import (
    _is_tsv_or_csv,
    _read_delim,
    load_singlecell,
    score_bulk_from_table,
    score_singlecell_adata,
    read_pathway_csv,
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


def _resolve_bulk_gene_col(counts_table, requested: Optional[str]) -> str:
    if requested:
        return requested
    cols = {c.lower(): c for c in counts_table.columns}
    return cols.get("gene") or cols.get("symbol") or list(counts_table.columns)[0]


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
        pathway_csvs = list(args.pathway_csv)
        labels = list(args.label) if args.label else None
        resolved_gene_col = _resolve_bulk_gene_col(df, args.gene_col)
        gene_col_for_match = resolved_gene_col if resolved_gene_col in df.columns else None
        detected_gene_col = gene_col_for_match or _resolve_bulk_gene_col(df, None)
        unmatched_pathways = []
        if gene_col_for_match:
            present_genes = (
                df[gene_col_for_match]
                .astype(str)
                .str.strip()
            )
            present_keys = {g.lower() for g in present_genes if g}
            for p in pathway_csvs:
                pw = read_pathway_csv(p)
                pw_genes = pw["gene"].astype(str).str.strip()
                has_match = any(g and g.lower() in present_keys for g in pw_genes)
                if not has_match:
                    unmatched_pathways.append(os.path.basename(p))
        if unmatched_pathways:
            msg = "No pathway genes matched input counts \u2014 check gene column (--gene_col) or provide a gene_map.csv"
            print(msg, file=sys.stderr)
            print(
                f"Detected gene column: '{detected_gene_col}'. Override with --gene_col <column> if needed.",
                file=sys.stderr,
            )
            print(
                "Pathways with no matches: " + ", ".join(unmatched_pathways),
                file=sys.stderr,
            )
        if args.scoring_style == "r":
            from .core import score_bulk_r_style_from_table
            out = score_bulk_r_style_from_table(
                df,
                gene_col=resolved_gene_col,
                pathway_csvs=pathway_csvs,
                labels=labels,
                collapse=args.collapse_duplicates,
                boost_by_evidence=not args.no_boost,
            )
        else:
            out = score_bulk_from_table(
                df,
                gene_col=resolved_gene_col,
                pathway_csvs=pathway_csvs,
                labels=labels,
            )
        # If sample info provided, join and produce summaries + PDF
        scores = out.copy()
        if args.sample_info:
            import pandas as pd
            si = pd.read_csv(args.sample_info)
            if args.sample_col not in si.columns or args.treatment_col not in si.columns:
                raise SystemExit(f"--sample_info must contain columns: {args.sample_col},{args.treatment_col}")
            si2 = si[[args.sample_col, args.treatment_col]].rename(columns={args.sample_col: "sample", args.treatment_col: "treatment"})
            # Warn on unmatched samples between scores and sample_info
            score_samples = set(scores["sample"].astype(str).unique())
            info_samples = set(si2["sample"].astype(str).unique())
            missing_in_info = sorted(score_samples - info_samples)
            extra_in_info = sorted(info_samples - score_samples)
            if missing_in_info:
                msg = (
                    f"Warning: {len(missing_in_info)} sample(s) in results are missing from --sample_info: "
                    + ", ".join(missing_in_info[:10])
                    + (" …" if len(missing_in_info) > 10 else "")
                )
                print(msg, file=sys.stderr)
            if extra_in_info:
                msg = (
                    f"Warning: {len(extra_in_info)} sample(s) in --sample_info not found in results: "
                    + ", ".join(extra_in_info[:10])
                    + (" …" if len(extra_in_info) > 10 else "")
                )
                print(msg, file=sys.stderr)

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
                else:
                    levels = set(scores["treatment"].dropna().astype(str).unique())
                    if ctrl not in levels and len(levels) > 0:
                        print(
                            f"Warning: --control_label '{ctrl}' not found in treatments: " + ", ".join(sorted(levels)),
                            file=sys.stderr,
                        )
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
