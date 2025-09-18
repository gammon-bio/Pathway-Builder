#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unified CLI to score weighted signaling pathways on bulk RNA-seq and single-nucleus/cell data.

Two modes (choose exactly one input):
  - Bulk:      --bulk_counts <csv/tsv> (genes x samples; normalized counts)
  - Single-cell: --sn_data <path|glob> (.h5ad, 10x .h5, or 10x mtx dir)

Multiple --pathway_csv are supported. Each file must have columns [gene, weight] (case-insensitive).
Weights are auto-normalized to sum=1 over detected genes before scoring.

Weighted score definition:
  score_w(sample or cell) = sum_i w_i * expr_i
  - Single-cell: expr_i uses library-size normalization (normalize_total target_sum=1e4) + log1p.
  - Bulk: expr_i uses provided normalized values (TPM/CPM/VST); no additional transform.

Outputs (deterministic layout under --output_dir):
  - Bulk:  scores_bulk.csv
  - SC:    scores_sc_by_celltype.csv, deltas.csv (if two conditions), figs/*.png, report.pdf (optional)

Notes:
  - Gene symbol matching is case-insensitive. Duplicates are collapsed by sum (bulk and single-cell).
  - For single-cell, if --celltype_col is not present, a simple marker-based guess is computed.
  - For single-cell, if --condition_col is not provided, the column is inferred if present; otherwise all cells are one condition.
  - For deltas (Case - Control), if labels are not specified, control is inferred by matching 'control|ctrl|vehicle|veh'.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import warnings
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# ----------------------------- Utilities ------------------------------------


def _is_tsv_or_csv(path: str) -> bool:
    return any(path.lower().endswith(ext) for ext in (".tsv", ".csv", ".txt"))


def _read_delim(path: str) -> pd.DataFrame:
    if path.lower().endswith(".tsv") or path.lower().endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _collapse_duplicates_sum(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    # Collapse duplicate gene symbols by summing expression across duplicates
    # Keep the first occurrence of gene name casing
    df2 = df.copy()
    df2["__key"] = df2[gene_col].astype(str).str.strip().str.lower()
    num_cols = [c for c in df2.columns if c not in (gene_col, "__key")]
    agg = df2.groupby("__key", as_index=False, sort=False)[num_cols].sum()
    # map back representative gene symbol (first seen)
    first_map = (
        df2.drop_duplicates("__key")["__key"]
        .reset_index(drop=True)
        .to_frame()
        .assign(symbol=df2.drop_duplicates("__key")[gene_col].values)
    )
    rep = {k: v for k, v in zip(first_map["__key"], first_map["symbol"])}
    agg[gene_col] = agg["__key"].map(rep)
    cols = [gene_col] + num_cols
    return agg[cols]


def _normalize_weights_present(
    pw_df: pd.DataFrame, present_genes: Iterable[str]
) -> pd.DataFrame:
    if pw_df.empty:
        return pw_df
    gset = {g.lower() for g in present_genes}
    df = pw_df.copy()
    df["__key"] = df["gene"].astype(str).str.strip().str.lower()
    df = df[df["__key"].isin(gset)]
    if df.empty:
        return df.drop(columns=["__key"])  # empty
    s = df["weight"].sum()
    if s > 0:
        df["weight"] = df["weight"] / s
    return df.drop(columns=["__key"])


def _read_pathway_csv(path: str) -> pd.DataFrame:
    df = _read_delim(path)
    colmap = {c.lower(): c for c in df.columns}
    gcol = colmap.get("gene") or colmap.get("symbol") or list(df.columns)[0]
    if "weight" in colmap:
        wcol = colmap["weight"]
    else:
        df["weight"] = 1.0
        wcol = "weight"
    out = (
        df[[gcol, wcol]]
        .rename(columns={gcol: "gene", wcol: "weight"})
        .assign(gene=lambda d: d["gene"].astype(str).str.strip())
    )
    # Coerce weights
    out["weight"] = pd.to_numeric(out["weight"], errors="coerce").fillna(1.0)
    # Collapse duplicate genes by summing weights, then normalize later
    out = out.groupby(out["gene"].str.lower(), as_index=False).agg(
        gene=("gene", "first"), weight=("weight", "sum")
    )[["gene", "weight"]]
    return out


def _infer_control_label(levels: Sequence[str]) -> Optional[str]:
    for lv in levels:
        s = str(lv).lower()
        if any(k in s for k in ("control", "ctrl", "vehicle", "veh")):
            return lv
    return None


def _safe_log(msg: str) -> None:
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def _sanitize_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in label)


# ----------------------------- Single-cell ----------------------------------


def _import_scanpy_if_needed():
    import importlib

    sc = importlib.import_module("scanpy")
    ad = importlib.import_module("anndata")
    return sc, ad


def _read_10x_h5(path: str):
    sc, _ = _import_scanpy_if_needed()
    try:
        return sc.read_10x_h5(path, gex_only=True)
    except TypeError:
        return sc.read_10x_h5(path)


def _read_10x_mtx(path: str):
    sc, _ = _import_scanpy_if_needed()
    # Accept either the matrix dir, or a dir containing filtered_feature_bc_matrix
    mdir = path
    if os.path.isdir(os.path.join(path, "filtered_feature_bc_matrix")):
        mdir = os.path.join(path, "filtered_feature_bc_matrix")
    A = sc.read_10x_mtx(mdir, var_names="gene_symbols", make_unique=False)
    # Enforce RNA only if feature_types present
    if "feature_types" in A.var.columns:
        keep = (
            A.var["feature_types"]
            .astype(str)
            .str.contains("Gene Expression", case=False, regex=False)
        )
        try:
            A = A[:, keep].copy()
        except Exception:
            pass
    return A


def _collapse_adata_duplicates_sum(adata) -> any:
    # Collapse duplicate var_names by summing columns using a fast sparse mapping
    import numpy as _np
    from scipy import sparse as _sp

    keys = pd.Index(adata.var_names).str.lower().to_numpy()
    # Factorize keys preserving first appearance order
    cats = pd.unique(keys)
    code_map = {k: i for i, k in enumerate(cats)}
    codes = _np.fromiter((code_map[k] for k in keys), dtype=_np.int64)
    n_groups = len(cats)
    if n_groups == adata.n_vars:
        return adata
    # Build mapping matrix (n_vars x n_groups)
    M = _sp.csr_matrix(
        (_np.ones(adata.n_vars, dtype=_np.float32), (_np.arange(adata.n_vars), codes)),
        shape=(adata.n_vars, n_groups),
    )
    # Multiply X (n_obs x n_vars) by M to sum duplicate columns
    if _sp.issparse(adata.X):
        X_new = adata.X.dot(M).tocsr()
    else:
        X_new = _sp.csr_matrix(adata.X).dot(M).tocsr()
    # Representative gene names (first seen)
    rep_names = []
    seen = set()
    for name, code in zip(adata.var_names, codes):
        if code not in seen:
            rep_names.append(name)
            seen.add(code)
            if len(seen) == n_groups:
                break
    # Build new AnnData
    _, ad = _import_scanpy_if_needed()
    var_new = pd.DataFrame(index=pd.Index(rep_names))
    return ad.AnnData(X=X_new, obs=adata.obs.copy(), var=var_new)


def _annotate_celltypes_simple(adata) -> None:
    sc, _ = _import_scanpy_if_needed()
    markers = {
        "Myofiber": ["Acta1", "Ttn", "Myh1", "Myh2", "Myh4", "Myh7", "Tnnt3"],
        "Satellite": ["Pax7", "Myf5"],
        "FAP": ["Pdgfra", "Dcn", "Col1a1", "Col1a2"],
        "Endothelial": ["Pecam1", "Kdr", "Cdh5"],
        "Pericyte": ["Rgs5", "Pdgfrb", "Kcnj8"],
        "Schwann": ["Mpz", "Plp1", "S100b", "Ngfr"],
        "Immune": ["Ptprc", "Lyz2", "Cd68", "Itgam"],
    }
    scores = {}
    # log-normalized required for scanpy score_genes; compute lightweight if missing
    if "__norm_done__" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns["__norm_done__"] = True
    for ct, genes in markers.items():
        g = [g for g in genes if g in adata.var_names]
        if len(g) >= 2:
            sc.tl.score_genes(
                adata, gene_list=g, score_name=f"__score_{ct}", use_raw=False
            )
            scores[ct] = f"__score_{ct}"
    score_cols = [v for v in scores.values() if v in adata.obs.columns]
    if score_cols:
        adata.obs["celltype_guess"] = (
            adata.obs[score_cols]
            .idxmax(axis=1)
            .str.replace("__score_", "", regex=False)
        )
    else:
        adata.obs["celltype_guess"] = "Unknown"


def load_singlecell(sn_data: str) -> any:
    sc, ad = _import_scanpy_if_needed()
    paths: List[str] = []
    if os.path.isdir(sn_data):
        paths = [sn_data]
    elif any(ch in sn_data for ch in "*?[]"):
        paths = sorted(glob.glob(sn_data))
    else:
        paths = [sn_data]
    if not paths:
        raise FileNotFoundError(f"No files matched for --sn_data={sn_data}")

    adatas = []
    for p in paths:
        base = os.path.basename(p)
        sample = os.path.splitext(base)[0]
        if os.path.isdir(p):
            A = _read_10x_mtx(p)
        elif p.lower().endswith(".h5"):
            A = _read_10x_h5(p)
        elif p.lower().endswith(".h5ad"):
            A = sc.read_h5ad(p)
        elif p.lower().endswith(".rds"):
            raise RuntimeError(
                "Seurat .rds import not supported in this CLI. Please convert to .h5ad (SeuratDisk) and retry."
            )
        else:
            raise RuntimeError(f"Unsupported file: {p}")
        # Enforce RNA-only features if annotated
        if "feature_types" in getattr(A, "var", pd.DataFrame()).columns:
            keep = (
                A.var["feature_types"]
                .astype(str)
                .str.contains("Gene Expression", case=False, regex=False)
            )
            try:
                A = A[:, keep].copy()
            except Exception:
                pass
        # Prefer gene symbols if available
        if "gene_symbols" in A.var.columns:
            A.var_names = A.var["gene_symbols"].astype(str)
        # Collapse duplicates by summing
        A = _collapse_adata_duplicates_sum(A)
        # Prefix obs with sample and add sample column
        A.obs_names = pd.Index([f"{sample}-{i}" for i in range(A.n_obs)])
        A.obs["sample"] = sample
        adatas.append(A)

    if len(adatas) == 1:
        adata = adatas[0]
    else:
        # Preserve per-file sample labels already stored in A.obs['sample']
        adata = ad.concat(adatas, join="outer", index_unique=None)

    # Basic filter: drop genes not seen in >= 5 cells
    sc.pp.filter_genes(adata, min_cells=5)
    return adata


def score_singlecell(
    adata,
    pathway_csvs: List[str],
    labels: Optional[List[str]],
    celltype_col: Optional[str],
    condition_col: Optional[str],
    output_dir: str,
    make_pdf: bool = True,
    min_cells_per_group: int = 20,
) -> None:
    sc, _ = _import_scanpy_if_needed()
    warnings.filterwarnings("ignore", category=UserWarning)

    # Normalize and log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.uns["__norm_done__"] = True

    # Celltype column
    if celltype_col and celltype_col in adata.obs.columns:
        ctc = celltype_col
    else:
        _safe_log(
            "[info] --celltype_col not provided or missing; computing simple marker-based celltype_guess"
        )
        _annotate_celltypes_simple(adata)
        ctc = "celltype_guess"

    # Condition column
    if condition_col and condition_col in adata.obs.columns:
        cond = condition_col
    else:
        cond = None
        if "condition" in adata.obs.columns:
            cond = "condition"
        else:
            adata.obs["condition"] = "All"
            cond = "condition"

    # Read pathways
    if not pathway_csvs:
        raise SystemExit("--pathway_csv required (one or more)")
    if labels and len(labels) != len(pathway_csvs):
        _safe_log(
            "[warn] Number of --label entries does not match --pathway_csv; using filenames as labels where missing."
        )

    label_list: List[str] = []
    pw_tables: List[pd.DataFrame] = []
    for i, p in enumerate(pathway_csvs):
        lab = None
        if labels and i < len(labels):
            lab = labels[i]
        if not lab:
            lab = os.path.splitext(os.path.basename(p))[0]
        label_list.append(lab)
        pw_tables.append(_read_pathway_csv(p))

    # Build case-insensitive map for var_names
    var_map = {v.lower(): v for v in adata.var_names}

    # Compute cell-level scores per pathway
    from scipy import sparse as sp

    for lab, pw in zip(label_list, pw_tables):
        present = [var_map[g.lower()] for g in pw["gene"] if g.lower() in var_map]
        pw_norm = _normalize_weights_present(pw, present_genes=present)
        lab_s = _sanitize_label(lab)
        out_col = f"score_w__{lab_s}"
        if pw_norm.empty:
            _safe_log(f"[warn] No genes detected for pathway '{lab}' — setting zeros")
            adata.obs[out_col] = 0.0
            continue
        genes = [var_map[g.lower()] for g in pw_norm["gene"]]
        w = pw_norm["weight"].to_numpy().astype(float)
        Xsub = adata[:, genes].X
        if sp.issparse(Xsub):
            scv = Xsub.dot(w)
            if hasattr(scv, "A1"):
                scv = scv.A1
        else:
            scv = (Xsub * w).sum(axis=1)
        adata.obs[out_col] = np.asarray(scv).ravel()

    # Summaries per celltype x condition
    score_cols = [c for c in adata.obs.columns if c.startswith("score_w__")]
    summary = (
        adata.obs[[ctc, cond] + score_cols]
        .groupby([ctc, cond], observed=False)
        .agg(["mean", "count"])
    )
    # Flatten columns
    summary.columns = [
        "__".join([a] + ([b] if isinstance(b, str) else list(b)))
        for a, b in summary.columns.values
    ]
    summary = summary.reset_index().rename(columns={ctc: "celltype", cond: "condition"})
    # Extract per-pathway mean and count
    rows = []
    for _, r in summary.iterrows():
        for lab in label_list:
            lab_s = _sanitize_label(lab)
            rows.append(
                dict(
                    pathway=lab,
                    celltype=r["celltype"],
                    condition=r["condition"],
                    score_w=r.get(f"score_w__{lab_s}__mean", np.nan),
                    n_cells=int(r.get(f"score_w__{lab_s}__count", 0)),
                )
            )
    df_ct = pd.DataFrame(rows)
    # Guard tiny groups
    df_ct["flag_small_group"] = df_ct["n_cells"] < int(min_cells_per_group)
    _ensure_dir(output_dir)
    df_ct.to_csv(os.path.join(output_dir, "scores_sc_by_celltype.csv"), index=False)

    # Write cell-level scores for downstream plotting
    cell_df = adata.obs[["sample", ctc, cond] + score_cols].copy()
    cell_df = cell_df.rename(columns={ctc: "celltype", cond: "condition"})
    cell_df.to_csv(os.path.join(output_dir, "scores_sc_cells.csv"), index=False)

    # Deltas (Case - Control) if two conditions
    deltas_path = os.path.join(output_dir, "deltas.csv")
    cond_levels = list(df_ct["condition"].dropna().unique())
    if len(cond_levels) == 2:
        ctrl = _infer_control_label(cond_levels) or sorted(cond_levels)[0]
        case = [x for x in cond_levels if x != ctrl][0]
        piv = df_ct.pivot_table(
            index=["pathway", "celltype"], columns="condition", values="score_w"
        )
        piv["delta_case_minus_control"] = piv.get(case) - piv.get(ctrl)
        piv = piv.reset_index()[["pathway", "celltype", "delta_case_minus_control"]]
        piv.to_csv(deltas_path, index=False)
    else:
        # write empty placeholder for deterministic outputs
        pd.DataFrame(
            columns=["pathway", "celltype", "delta_case_minus_control"]
        ).to_csv(deltas_path, index=False)

    # Figures + optional PDF
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        figdir = os.path.join(output_dir, "figs")
        _ensure_dir(figdir)
        pdf_path = os.path.join(output_dir, "report.pdf")
        pdf = PdfPages(pdf_path) if make_pdf else None
        for lab in label_list:
            plt.figure(figsize=(8, 4))
            sub = df_ct[df_ct["pathway"] == lab].copy()
            # Sort by condition then celltype for determinism
            sub = sub.sort_values(["condition", "celltype"]) if not sub.empty else sub
            for cond_val in sub["condition"].unique():
                ss = sub[sub["condition"] == cond_val]
                plt.plot(
                    ss["celltype"],
                    ss["score_w"],
                    marker="o",
                    linestyle="-",
                    label=str(cond_val),
                )
            plt.xticks(rotation=45, ha="right")
            plt.ylabel("Weighted score (mean)")
            plt.title(f"{lab} signaling (weighted)")
            plt.legend()
            plt.tight_layout()
            png_path = os.path.join(figdir, f"{_sanitize_label(lab)}_grouped.png")
            plt.savefig(png_path, dpi=150)
            if pdf is not None:
                pdf.savefig()
                plt.close()
            else:
                plt.close()
        # Additional PDF pages: deltas (case - control) per pathway if exactly two conditions
        cond_levels_pdf = list(df_ct["condition"].dropna().unique())
        if len(cond_levels_pdf) == 2 and pdf is not None:
            ctrl_pdf = (
                _infer_control_label(cond_levels_pdf) or sorted(cond_levels_pdf)[0]
            )
            case_pdf = [x for x in cond_levels_pdf if x != ctrl_pdf][0]
            for lab in label_list:
                sub = df_ct[df_ct["pathway"] == lab]
                if sub.empty:
                    continue
                piv = sub.pivot_table(
                    index="celltype", columns="condition", values="score_w"
                )
                if ctrl_pdf in piv.columns and case_pdf in piv.columns:
                    delta = (piv[case_pdf] - piv[ctrl_pdf]).rename("delta")
                    plt.figure(figsize=(max(6, len(delta) * 0.5), 3.5))
                    plt.bar(delta.index, delta.values, color="#4C78A8")
                    plt.axhline(0, color="#333", lw=1)
                    plt.xticks(rotation=45, ha="right")
                    plt.ylabel(f"Delta ({case_pdf} - {ctrl_pdf})")
                    plt.title(f"{lab} — delta by celltype")
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()
        # Additional pages: per-pathway boxplots of cell-level scores by condition, faceted by celltype
        try:
            import seaborn as sns  # improves boxplots
        except Exception:
            sns = None
        for lab in label_list:
            lab_s = _sanitize_label(lab)
            col = f"score_w__{lab_s}"
            if col not in cell_df.columns:
                continue
            plt.figure(figsize=(10, 6))
            if sns is not None:
                sns.boxplot(data=cell_df, x="celltype", y=col, hue="condition")
                sns.stripplot(
                    data=cell_df,
                    x="celltype",
                    y=col,
                    hue="condition",
                    dodge=True,
                    color="#222",
                    size=2,
                    alpha=0.4,
                )
                # remove duplicate legends
                handles, labels2 = plt.gca().get_legend_handles_labels()
                if handles:
                    plt.legend(handles[:2], labels2[:2])
            else:
                # Fallback: basic matplotlib boxplots per condition for each celltype
                ct_list = list(cell_df["celltype"].dropna().unique())
                cond_list = list(cell_df["condition"].dropna().unique())
                data = []
                xticks = []
                for i, ct in enumerate(ct_list):
                    for j, cv in enumerate(cond_list):
                        vals = cell_df[
                            (cell_df["celltype"] == ct) & (cell_df["condition"] == cv)
                        ][col].to_numpy()
                        data.append(vals)
                        xticks.append(f"{ct}\n{cv}")
                plt.boxplot(data, labels=xticks, notch=False)
                plt.xticks(rotation=45, ha="right")
            plt.ylabel("Cell score")
            plt.title(f"{lab} — cell-level scores by celltype/condition")
            plt.tight_layout()
            if pdf is not None:
                pdf.savefig()
                plt.close()
            else:
                plt.close()
        if pdf is not None:
            pdf.close()
    except Exception as e:
        _safe_log(f"[warn] plotting/report failed: {e}")
    # Additional comparisons: tidy mean per celltype x condition and delta bars
    try:
        # Tidy table from df_ct
        tidy = df_ct.rename(columns={"score_w": "mean", "n_cells": "n"})[
            ["pathway", "celltype", "condition", "mean", "n"]
        ]
        tidy.to_csv(
            os.path.join(output_dir, "scores_by_celltype_condition_tidy.csv"),
            index=False,
        )

        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        # Grouped bars by condition
        pdf_cmp = PdfPages(
            os.path.join(output_dir, "comparisons_by_celltype_condition.pdf")
        )
        for lab in df_ct["pathway"].unique():
            sub = df_ct[df_ct["pathway"] == lab]
            piv = sub.pivot_table(
                index="celltype", columns="condition", values="score_w"
            )
            piv = piv.sort_index()
            fig, ax = plt.subplots(figsize=(max(6, len(piv) * 0.5), 4))
            conds = list(piv.columns)
            width = 0.8 / max(1, len(conds))
            x = range(len(piv))
            for i, c in enumerate(conds):
                ax.bar(
                    [xi + i * width for xi in x],
                    piv[c].values,
                    width=width,
                    label=str(c),
                )
            ax.set_xticks([xi + (len(conds) - 1) * width / 2 for xi in x])
            ax.set_xticklabels(piv.index, rotation=45, ha="right")
            ax.set_ylabel("Score (mean)")
            ax.set_title(f"{lab}: mean by celltype and condition")
            ax.legend()
            fig.tight_layout()
            pdf_cmp.savefig(fig)
            plt.close(fig)
        pdf_cmp.close()

        # Delta table and PDF
        cond_levels = list(df_ct["condition"].dropna().unique())
        if len(cond_levels) >= 2:
            ctrl = _infer_control_label(cond_levels) or cond_levels[0]
            case = [x for x in cond_levels if x != ctrl][0]
            dels = []
            for lab in df_ct["pathway"].unique():
                sub = df_ct[df_ct["pathway"] == lab]
                piv = sub.pivot_table(
                    index="celltype", columns="condition", values="score_w"
                )
                d = (piv.get(case) - piv.get(ctrl)).rename("delta").reset_index()
                d.insert(0, "pathway", lab)
                dels.append(d)
            dtab = pd.concat(dels, ignore_index=True)
            dtab.to_csv(os.path.join(output_dir, "deltas_by_celltype.csv"), index=False)
            pdf_delta = PdfPages(os.path.join(output_dir, "comparisons_deltas.pdf"))
            for lab in dtab["pathway"].unique():
                sub = dtab[dtab["pathway"] == lab].sort_values("celltype")
                fig, ax = plt.subplots(figsize=(max(6, len(sub) * 0.5), 3.5))
                ax.bar(sub["celltype"], sub["delta"].values, color="#4C78A8")
                ax.axhline(0, color="#333", lw=1)
                ax.set_ylabel(f"Delta ({case} - {ctrl})")
                ax.set_title(f"{lab}: delta by celltype")
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
                fig.tight_layout()
                pdf_delta.savefig(fig)
                plt.close(fig)
            pdf_delta.close()
    except Exception as e:
        _safe_log(f"[warn] additional comparisons failed: {e}")


# ----------------------------- Bulk -----------------------------------------


def score_bulk(
    bulk_counts: str,
    pathway_csvs: List[str],
    labels: Optional[List[str]],
    output_dir: str,
    gene_col: Optional[str] = None,
) -> None:
    if not _is_tsv_or_csv(bulk_counts):
        raise SystemExit("--bulk_counts must be .csv, .tsv, or .txt (genes x samples)")
    df = _read_delim(bulk_counts)
    if df.empty or df.shape[1] < 2:
        raise SystemExit(
            "Bulk counts table appears malformed: need gene column + at least 1 sample column"
        )

    # Auto-detect gene column if not provided
    if gene_col is None:
        cols = {c.lower(): c for c in df.columns}
        gene_col = cols.get("gene") or cols.get("symbol") or list(df.columns)[0]
    if gene_col not in df.columns:
        raise SystemExit(f"Gene column '{gene_col}' not found in bulk table")

    # Coerce numeric samples and drop non-numeric columns except gene_col
    sample_cols = [c for c in df.columns if c != gene_col]
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    df = df[[gene_col] + sample_cols]

    # Collapse duplicate genes by sum
    df = _collapse_duplicates_sum(df, gene_col=gene_col)

    # Case-insensitive matching helper
    gene_map = {g.lower(): g for g in df[gene_col].astype(str)}

    # Read pathways
    if not pathway_csvs:
        raise SystemExit("--pathway_csv required (one or more)")
    if labels and len(labels) != len(pathway_csvs):
        _safe_log(
            "[warn] Number of --label entries does not match --pathway_csv; using filenames as labels where missing."
        )
    label_list: List[str] = []
    pw_tables: List[pd.DataFrame] = []
    for i, p in enumerate(pathway_csvs):
        lab = None
        if labels and i < len(labels):
            lab = labels[i]
        if not lab:
            lab = os.path.splitext(os.path.basename(p))[0]
        label_list.append(lab)
        pw_tables.append(_read_pathway_csv(p))

    # Compute weighted scores per sample
    rows = []
    for lab, pw in zip(label_list, pw_tables):
        present = [gene_map[g.lower()] for g in pw["gene"] if g.lower() in gene_map]
        pw_norm = _normalize_weights_present(pw, present_genes=present)
        if pw_norm.empty:
            _safe_log(f"[warn] No genes detected for pathway '{lab}' — emitting zeros")
            for s in sample_cols:
                rows.append(dict(pathway=lab, sample=s, score_w=0.0))
            continue
        genes = [gene_map[g.lower()] for g in pw_norm["gene"]]
        w = pw_norm["weight"].to_numpy().astype(float)
        # Extract submatrix (genes x samples) and compute weighted sum across genes
        sub = df.set_index(gene_col).loc[genes, sample_cols].to_numpy(dtype=float)
        scv = (w.reshape(-1, 1) * sub).sum(axis=0)
        for s, val in zip(sample_cols, scv):
            rows.append(dict(pathway=lab, sample=s, score_w=float(val)))

    out = pd.DataFrame(rows)
    _ensure_dir(output_dir)
    out.to_csv(os.path.join(output_dir, "scores_bulk.csv"), index=False)


# ----------------------------- CLI ------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Score weighted signaling pathways on bulk RNA-seq or single-cell data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = ap.add_argument_group("Inputs")
    src.add_argument(
        "--bulk_counts", help="Bulk normalized counts table (genes x samples; CSV/TSV)"
    )
    src.add_argument(
        "--sn_data",
        help="Single-cell input: .h5ad, 10x .h5, 10x mtx dir, or a glob pattern",
    )
    src.add_argument(
        "--pathway_csv",
        action="append",
        help="Pathway CSV(s) with columns [gene,weight]; can be repeated",
    )
    src.add_argument(
        "--label",
        action="append",
        help="Label(s) for each pathway CSV; can be repeated; defaults to filename",
    )

    scg = ap.add_argument_group("Single-cell options")
    scg.add_argument(
        "--celltype_col",
        help="obs column for cell-type; default uses a simple marker-based guess if missing",
    )
    scg.add_argument(
        "--condition_col", help="obs column for condition/case-control grouping"
    )
    scg.add_argument(
        "--min_cells_per_group",
        type=int,
        default=20,
        help="Flag groups with fewer cells as small",
    )
    scg.add_argument(
        "--no_pdf",
        action="store_true",
        help="Do not produce a PDF report; still writes PNGs",
    )
    scg.add_argument(
        "--sample_info",
        help="CSV with sample-to-condition mapping (columns: sample, condition)",
    )
    scg.add_argument(
        "--sample_col",
        default="sample",
        help="Column name in --sample_info for sample IDs",
    )
    scg.add_argument(
        "--info_condition_col",
        default="condition",
        help="Column in --sample_info to use as condition label",
    )

    bg = ap.add_argument_group("Bulk options")
    bg.add_argument(
        "--gene_col", help="Gene symbol column name in bulk table (auto-detected)"
    )

    out = ap.add_argument_group("Outputs")
    out.add_argument("--output_dir", required=True, help="Directory to write outputs")

    ex = ap.add_argument_group("Examples")
    ex.add_argument(
        "--example_bulk",
        action="store_true",
        help="Show bulk example usage and exit",
    )
    ex.add_argument(
        "--example_sc",
        action="store_true",
        help="Show single-cell example usage and exit",
    )
    return ap


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = build_parser()
    args = ap.parse_args(argv)

    if args.example_bulk:
        print(
            "python run_bdnf_enrichment.py --bulk_counts data/bulk_counts.tsv --pathway_csv genes/my_pathway.csv --label MY_PATHWAY --output_dir out/bulk_demo"
        )
        return
    if args.example_sc:
        print(
            "python run_bdnf_enrichment.py --sn_data 'data/GSE272085/*.h5' --pathway_csv genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv --label TRKB --condition_col condition --celltype_col celltype_guess --output_dir out/sc_demo"
        )
        return

    # Validate mode
    has_bulk = bool(args.bulk_counts)
    has_sc = bool(args.sn_data)
    if has_bulk == has_sc:
        raise SystemExit("Specify exactly one of --bulk_counts or --sn_data")
    if not args.pathway_csv:
        raise SystemExit("At least one --pathway_csv is required")

    if has_bulk:
        score_bulk(
            bulk_counts=args.bulk_counts,
            pathway_csvs=list(args.pathway_csv),
            labels=list(args.label) if args.label else None,
            output_dir=args.output_dir,
            gene_col=args.gene_col,
        )
        _safe_log(f"[ok] Wrote: {os.path.join(args.output_dir, 'scores_bulk.csv')}")
    else:
        adata = load_singlecell(args.sn_data)
        # Optional: merge sample info to set/override condition labels
        if getattr(args, "sample_info", None):
            try:
                sinfo = pd.read_csv(args.sample_info)
                scol = args.sample_col or "sample"
                ccol = args.info_condition_col or "condition"
                if scol in sinfo.columns and ccol in sinfo.columns:
                    m = sinfo[[scol, ccol]].rename(
                        columns={scol: "sample", ccol: "__info_condition"}
                    )
                    adata.obs = adata.obs.merge(m, on="sample", how="left")
                    if "__info_condition" in adata.obs.columns:
                        adata.obs["condition"] = adata.obs["__info_condition"].fillna(
                            adata.obs.get("condition", "All")
                        )
                        adata.obs.drop(columns=["__info_condition"], inplace=True)
                else:
                    _safe_log(
                        f"[warn] --sample_info missing required columns: {scol} and/or {ccol}"
                    )
            except Exception as e:
                _safe_log(f"[warn] failed to merge --sample_info: {e}")
        score_singlecell(
            adata=adata,
            pathway_csvs=list(args.pathway_csv),
            labels=list(args.label) if args.label else None,
            celltype_col=args.celltype_col,
            condition_col=args.condition_col,
            output_dir=args.output_dir,
            make_pdf=not bool(args.no_pdf),
            min_cells_per_group=int(args.min_cells_per_group),
        )
        _safe_log(
            f"[ok] Wrote: {os.path.join(args.output_dir, 'scores_sc_by_celltype.csv')} and deltas.csv; figs/* and report.pdf if enabled"
        )


if __name__ == "__main__":
    main()
