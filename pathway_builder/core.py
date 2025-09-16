from __future__ import annotations

import glob
import os
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _is_tsv_or_csv(path: str) -> bool:
    return any(path.lower().endswith(ext) for ext in (".tsv", ".csv", ".txt"))


def _read_delim(path: str) -> pd.DataFrame:
    if path.lower().endswith(".tsv") or path.lower().endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _collapse_duplicates_sum(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    df2 = df.copy()
    df2["__key"] = df2[gene_col].astype(str).str.strip().str.lower()
    num_cols = [c for c in df2.columns if c not in (gene_col, "__key")]
    agg = df2.groupby("__key", as_index=False, sort=False)[num_cols].sum()
    # representative symbol (first seen)
    first = df2.drop_duplicates("__key")[["__key", gene_col]]
    rep = dict(zip(first["__key"], first[gene_col]))
    agg[gene_col] = agg["__key"].map(rep)
    cols = [gene_col] + num_cols
    return agg[cols]


def read_pathway_csv(path: str) -> pd.DataFrame:
    df = _read_delim(path)
    if df.empty:
        return pd.DataFrame(columns=["gene", "weight"])  # type: ignore
    colmap = {c.lower(): c for c in df.columns}
    gcol = (
        colmap.get("gene")
        or colmap.get("gene_symbol")
        or colmap.get("symbol")
        or list(df.columns)[0]
    )
    if "weight" in colmap:
        wcol = colmap["weight"]
    else:
        df["weight"] = 1.0
        wcol = "weight"
    # carry evidence_source if present (for R-style boost)
    evcol = colmap.get("evidence_source") or colmap.get("evidence")
    cols = [gcol, wcol] + ([evcol] if evcol else [])
    out = (
        df[cols]
        .rename(columns={gcol: "gene", wcol: "weight", (evcol or "evidence_source"): "evidence_source"})
        .assign(gene=lambda d: d["gene"].astype(str).str.strip())
    )
    out["weight"] = pd.to_numeric(out["weight"], errors="coerce").fillna(1.0)
    grp = {"gene": ("gene", "first"), "weight": ("weight", "sum")}
    if "evidence_source" in out.columns:
        grp["evidence_source"] = ("evidence_source", "first")
    out = out.groupby(out["gene"].str.lower(), as_index=False).agg(**grp)
    return out[[c for c in ("gene", "weight", "evidence_source") if c in out.columns]]


def normalize_weights_present(pw_df: pd.DataFrame, present_genes: Iterable[str]) -> pd.DataFrame:
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


def load_vst_counts_table(
    counts_table: pd.DataFrame,
    gene_col: Optional[str] = None,
    collapse: str = "mean",
) -> pd.DataFrame:
    """Like load_vst_csv.R: return numeric DataFrame with index=gene symbols, columns=samples.
    collapse: one of mean, maxvar, median, sum
    """
    if counts_table.empty or counts_table.shape[1] < 2:
        raise ValueError("counts_table malformed: need gene column + >=1 sample columns")
    df = counts_table.copy()
    if gene_col is None:
        candidates = [
            "gene",
            "Gene",
            "symbol",
            "Symbol",
            "gene_symbol",
            "GeneSymbol",
            "Gene.Name",
            "GeneID",
            "Geneid",
        ]
        for c in candidates:
            if c in df.columns:
                gene_col = c; break
        if gene_col is None:
            gene_col = df.columns[0]
    genes = df[gene_col].astype(str)
    mat = df.drop(columns=[gene_col])
    mat = mat.apply(pd.to_numeric, errors="coerce")
    mat.index = genes
    if not mat.index.has_duplicates:
        return mat
    collapse = collapse.lower()
    if collapse == "mean":
        grouped = mat.groupby(mat.index).mean()
        return grouped
    elif collapse == "median":
        return mat.groupby(mat.index).median()
    elif collapse == "sum":
        return mat.groupby(mat.index).sum()
    elif collapse == "maxvar":
        # keep row with highest across-sample variance
        def pick_maxvar(group: pd.DataFrame) -> pd.Series:
            if group.shape[0] == 1:
                return group.iloc[0]
            v = group.var(axis=1, ddof=1)
            return group.loc[v.idxmax()]
        out = mat.groupby(mat.index, group_keys=False).apply(pick_maxvar)
        return out
    else:
        raise ValueError("collapse must be one of mean, median, sum, maxvar")


def _boost_weights_by_evidence(pw_df: pd.DataFrame, boost: bool = True) -> pd.DataFrame:
    if not boost or "evidence_source" not in pw_df.columns:
        return pw_df
    df = pw_df.copy()
    es = df["evidence_source"].astype(str).fillna("")
    has_kegg = es.str.contains("KEGG", case=False, na=False)
    has_reac = es.str.contains("Reactome", case=False, na=False)
    df["weight"] = df["weight"].astype(float) + 0.05 * has_kegg + 0.05 * has_reac
    df["weight"] = df["weight"].clip(upper=1.2)
    return df


def score_bulk_r_style_from_table(
    counts_table: pd.DataFrame,
    gene_col: Optional[str],
    pathway_csvs: List[str],
    labels: Optional[List[str]] = None,
    collapse: str = "mean",
    boost_by_evidence: bool = True,
) -> pd.DataFrame:
    """R-style scoring: per-gene z-score across samples, then weighted mean per sample.
    weights optionally boosted by evidence_source.
    """
    X = load_vst_counts_table(counts_table, gene_col=gene_col, collapse=collapse)
    # Precompute row-wise mean and std for z-score
    mu = X.mean(axis=1)
    sd = X.std(axis=1, ddof=1).replace(0.0, 1.0)

    if labels and len(labels) != len(pathway_csvs):
        labels = labels + [None] * (len(pathway_csvs) - len(labels))
    rows = []
    for i, p in enumerate(pathway_csvs):
        lab = labels[i] if labels and i < len(labels) and labels[i] else os.path.splitext(os.path.basename(p))[0]
        pw = read_pathway_csv(p)
        pw = _boost_weights_by_evidence(pw, boost=boost_by_evidence)
        # intersect with X rows
        present = [g for g in pw["gene"] if g in X.index]
        if len(present) == 0:
            for s in X.columns:
                rows.append(dict(pathway=lab, sample=s, score_w=0.0))
            continue
        # z-score subset rows
        sub = X.loc[present]
        subZ = (sub.sub(mu.loc[present], axis=0)).div(sd.loc[present], axis=0)
        w = pw.set_index("gene").loc[present, "weight"].astype(float)
        sumw = float(w.sum()) if float(w.sum()) != 0 else 1.0
        wnorm = (w / sumw).to_numpy()
        # Weighted mean: 1 x samples
        vals = wnorm @ subZ.to_numpy()
        for s, v in zip(X.columns, vals):
            rows.append(dict(pathway=lab, sample=s, score_w=float(v)))
    return pd.DataFrame(rows)


def score_bulk_from_table(
    counts_table: pd.DataFrame,
    gene_col: Optional[str],
    pathway_csvs: List[str],
    labels: Optional[List[str]] = None,
) -> pd.DataFrame:
    if counts_table.empty or counts_table.shape[1] < 2:
        raise ValueError("counts_table malformed: require gene column + >=1 sample column")
    # detect gene column
    if gene_col is None:
        cols = {c.lower(): c for c in counts_table.columns}
        gene_col = cols.get("gene") or cols.get("symbol") or list(counts_table.columns)[0]
    if gene_col not in counts_table.columns:
        raise ValueError(f"Gene column '{gene_col}' not found")
    df = counts_table.copy()
    sample_cols = [c for c in df.columns if c != gene_col]
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    df = df[[gene_col] + sample_cols]
    df = _collapse_duplicates_sum(df, gene_col=gene_col)
    gene_map = {g.lower(): g for g in df[gene_col].astype(str)}
    if labels and len(labels) != len(pathway_csvs):
        # align shorter labels with csvs
        labels = labels + [None] * (len(pathway_csvs) - len(labels))
    rows = []
    for i, p in enumerate(pathway_csvs):
        lab = labels[i] if labels and i < len(labels) and labels[i] else os.path.splitext(os.path.basename(p))[0]
        pw = read_pathway_csv(p)
        present = [gene_map[g.lower()] for g in pw["gene"] if g.lower() in gene_map]
        pw_norm = normalize_weights_present(pw, present_genes=present)
        if pw_norm.empty:
            for s in sample_cols:
                rows.append(dict(pathway=lab, sample=s, score_w=0.0))
            continue
        genes = [gene_map[g.lower()] for g in pw_norm["gene"]]
        w = pw_norm["weight"].to_numpy().astype(float)
        sub = df.set_index(gene_col).loc[genes, sample_cols].to_numpy(dtype=float)
        scv = (w.reshape(-1, 1) * sub).sum(axis=0)
        for s, val in zip(sample_cols, scv):
            rows.append(dict(pathway=lab, sample=s, score_w=float(val)))
    return pd.DataFrame(rows)


def welch_t_test(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float, float, float, float]:
    """Return t, df, p, mean_diff, ci_low, ci_high for Welch's t-test (two-sided)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float('nan'), float('nan'), float('nan'), float(np.nan), float(np.nan), float(np.nan)
    mx, my = float(np.mean(x)), float(np.mean(y))
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    se = np.sqrt(vx / nx + vy / ny)
    if se == 0:
        return float('inf'), float('nan'), 0.0, mx - my, float('nan'), float('nan')
    t = (mx - my) / se
    df = (vx / nx + vy / ny) ** 2 / ((vx**2) / ((nx**2) * (nx - 1)) + (vy**2) / ((ny**2) * (ny - 1)))
    # two-sided p-value
    try:
        from scipy.stats import t as tdist
        p = 2 * tdist.sf(abs(t), df)
        # 95% CI for mean difference
        ql, qh = tdist.ppf([0.025, 0.975], df)
        ci_low = (mx - my) + ql * se
        ci_high = (mx - my) + qh * se
    except Exception:
        p = float('nan'); ci_low = float('nan'); ci_high = float('nan')
    return float(t), float(df), float(p), float(mx - my), float(ci_low), float(ci_high)


def hedges_g(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float('nan')
    mx, my = float(np.mean(x)), float(np.mean(y))
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    s_pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    if s_pooled == 0:
        return float('inf')
    d = (mx - my) / s_pooled
    J = 1 - (3 / (4 * (nx + ny) - 9))
    return float(d * J)


def load_singlecell(sn_data: str):  # -> anndata.AnnData
    import importlib

    sc = importlib.import_module("scanpy")
    ad = importlib.import_module("anndata")

    paths: List[str]
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
            # 10x mtx directory
            mdir = p
            if os.path.isdir(os.path.join(p, "filtered_feature_bc_matrix")):
                mdir = os.path.join(p, "filtered_feature_bc_matrix")
            A = sc.read_10x_mtx(mdir, var_names="gene_symbols", make_unique=False)
        elif p.lower().endswith(".h5"):
            try:
                A = sc.read_10x_h5(p, gex_only=True)
            except TypeError:
                A = sc.read_10x_h5(p)
        elif p.lower().endswith(".h5ad"):
            A = sc.read_h5ad(p)
        else:
            raise RuntimeError(f"Unsupported file: {p}")
        if "gene_symbols" in A.var.columns:
            A.var_names = A.var["gene_symbols"].astype(str)
        if "feature_types" in getattr(A, "var", pd.DataFrame()).columns:
            keep = A.var["feature_types"].astype(str).str.contains("Gene Expression", case=False, regex=False)
            try:
                A = A[:, keep].copy()
            except Exception:
                pass
        A.obs_names = pd.Index([f"{sample}-{i}" for i in range(A.n_obs)])
        A.obs["sample"] = sample
        adatas.append(A)

    if len(adatas) == 1:
        adata = adatas[0]
    else:
        adata = ad.concat(adatas, join="outer", label="sample", index_unique=None)
    sc.pp.filter_genes(adata, min_cells=5)
    return adata


def score_singlecell_adata(
    adata,
    pathway_csvs: List[str],
    labels: Optional[List[str]],
    celltype_col: Optional[str],
    condition_col: Optional[str],
) -> pd.DataFrame:
    import importlib
    from typing import Sequence

    sc = importlib.import_module("scanpy")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # choose grouping columns
    if celltype_col and celltype_col in adata.obs.columns:
        ctc = celltype_col
    else:
        ctc = "celltype_guess"
        if ctc not in adata.obs.columns:
            adata.obs[ctc] = "Unknown"
    if condition_col and condition_col in adata.obs.columns:
        cond = condition_col
    else:
        if "condition" not in adata.obs.columns:
            adata.obs["condition"] = "All"
        cond = "condition"

    if not pathway_csvs:
        raise ValueError("pathway_csvs is required")

    if labels and len(labels) != len(pathway_csvs):
        labels = labels + [None] * (len(pathway_csvs) - len(labels))
    label_list: List[str] = []
    pw_tables: List[pd.DataFrame] = []
    for i, p in enumerate(pathway_csvs):
        lab = labels[i] if labels and i < len(labels) and labels[i] else os.path.splitext(os.path.basename(p))[0]
        label_list.append(lab)
        pw_tables.append(read_pathway_csv(p))

    var_map = {v.lower(): v for v in adata.var_names}
    from scipy import sparse as sp
    for lab, pw in zip(label_list, pw_tables):
        present = [var_map[g.lower()] for g in pw["gene"] if g.lower() in var_map]
        pw_norm = normalize_weights_present(pw, present_genes=present)
        out_col = f"score_w__{lab}"
        if pw_norm.empty:
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

    score_cols = [c for c in adata.obs.columns if c.startswith("score_w__")]
    summary = (
        adata.obs[[ctc, cond] + score_cols].groupby([ctc, cond], observed=False).mean().reset_index()
    )
    return summary
