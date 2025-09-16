from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd


def _as_wide_table(scores: pd.DataFrame) -> pd.DataFrame:
    # Expect columns: pathway, sample, score_w, [treatment]
    idx_cols = [c for c in ["treatment", "sample"] if c in scores.columns]
    if not idx_cols:
        idx_cols = ["sample"]
    wide = scores.pivot_table(index=idx_cols, columns="pathway", values="score_w")
    return wide.reset_index()


def _compute_group_stats(scores: pd.DataFrame, treatment_col: str) -> pd.DataFrame:
    # Per pathway x treatment mean, sd, n, sem
    rows: List[Dict] = []
    for (pathway, trt), grp in scores.groupby(["pathway", treatment_col], observed=False):
        vals = grp["score_w"].astype(float).to_numpy()
        n = int(np.isfinite(vals).sum())
        mean = float(np.nanmean(vals)) if n else np.nan
        sd = float(np.nanstd(vals, ddof=1)) if n > 1 else np.nan
        sem = float(sd / np.sqrt(n)) if (n and np.isfinite(sd)) else np.nan
        rows.append(dict(pathway=pathway, treatment=trt, mean=mean, sd=sd, n=n, sem=sem))
    return pd.DataFrame(rows)


def _welch_stats(scores: pd.DataFrame, treatment_col: str, control_label: Optional[str] = None) -> pd.DataFrame:
    from .core import welch_t_test, hedges_g

    rows: List[Dict] = []
    for pw, grp in scores.groupby("pathway", observed=False):
        levels = list(grp[treatment_col].dropna().unique())
        if len(levels) != 2:
            continue
        if control_label and control_label in levels:
            a, b = control_label, [x for x in levels if x != control_label][0]
        else:
            a, b = levels[0], levels[1]
        xa = grp.loc[grp[treatment_col] == a, "score_w"].to_numpy()
        xb = grp.loc[grp[treatment_col] == b, "score_w"].to_numpy()
        t, df, p, mean_diff, ci_low, ci_high = welch_t_test(xb, xa)  # case minus control
        g = hedges_g(xb, xa)
        rows.append(dict(
            pathway=pw,
            group_a=a,
            group_b=b,
            n_a=int(np.isfinite(xa).sum()),
            n_b=int(np.isfinite(xb).sum()),
            mean_diff=mean_diff,
            t=t,
            df=df,
            p=p,
            ci_low=ci_low,
            ci_high=ci_high,
            hedges_g=g,
            test="welch",
        ))
    out = pd.DataFrame(rows)
    if not out.empty:
        try:
            from statsmodels.stats.multitest import multipletests
            out["q_BH"] = multipletests(out["p"].to_numpy(), method="fdr_bh")[1]
        except Exception:
            out["q_BH"] = np.nan
    else:
        out["q_BH"] = []
    return out


def _anova_stats(scores: pd.DataFrame, treatment_col: str) -> pd.DataFrame:
    """One-way ANOVA across all treatments for each pathway."""
    rows: List[Dict] = []
    try:
        from scipy.stats import f_oneway
    except Exception:
        # SciPy missing: return empty
        return pd.DataFrame(columns=["pathway", "k", "F", "p", "test", "q_BH"])  # type: ignore

    for pw, grp in scores.groupby("pathway", observed=False):
        groups = []
        labels = []
        for lv, g in grp.groupby(treatment_col, observed=False):
            arr = g["score_w"].astype(float).to_numpy()
            arr = arr[np.isfinite(arr)]
            if len(arr) >= 2:
                groups.append(arr)
                labels.append(lv)
        if len(groups) >= 2:
            F, p = f_oneway(*groups)
            rows.append(dict(pathway=pw, k=len(groups), F=float(F), p=float(p), test="anova"))
    out = pd.DataFrame(rows)
    if not out.empty:
        try:
            from statsmodels.stats.multitest import multipletests
            out["q_BH"] = multipletests(out["p"].to_numpy(), method="fdr_bh")[1]
        except Exception:
            out["q_BH"] = np.nan
    else:
        out["q_BH"] = []
    return out


def _infer_control_label(levels: Sequence[str]) -> Optional[str]:
    for lv in levels:
        s = str(lv).lower()
        if any(k in s for k in ("control", "ctrl", "vehicle", "veh")):
            return lv
    return None


def _tukey_posthoc(scores: pd.DataFrame, treatment_col: str) -> pd.DataFrame:
    """Tukey HSD for all treatment pairs per pathway. Returns long DataFrame with comparisons."""
    try:
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
    except Exception:
        return pd.DataFrame(columns=["pathway", "group1", "group2", "meandiff", "p_adj", "lower", "upper", "reject"])  # type: ignore

    rows: List[Dict] = []
    for pw, grp in scores.groupby("pathway", observed=False):
        vals = grp["score_w"].astype(float).to_numpy()
        groups = grp[treatment_col].astype(str).to_numpy()
        # Need at least 2 obs per group for stable intervals; let statsmodels handle edge cases
        try:
            tukey = pairwise_tukeyhsd(endog=vals, groups=groups, alpha=0.05)
        except Exception:
            continue
        res = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
        # Expected columns: group1, group2, meandiff, p-adj, lower, upper, reject
        res = res.rename(columns={"p-adj": "p_adj"})
        res.insert(0, "pathway", pw)
        rows.append(res)
    if not rows:
        return pd.DataFrame(columns=["pathway", "group1", "group2", "meandiff", "p_adj", "lower", "upper", "reject"])  # type: ignore
    out = pd.concat(rows, ignore_index=True)
    # Ensure types
    for c in ["meandiff", "p_adj", "lower", "upper"]:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def make_bulk_pdf_report(
    scores_with_treatment: pd.DataFrame,
    output_pdf_path: str,
    methods_note: Optional[str] = None,
    control_label: Optional[str] = None,
    style: str = "bar",
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    os.makedirs(os.path.dirname(output_pdf_path) or ".", exist_ok=True)
    tbl = _as_wide_table(scores_with_treatment)
    # Compute stats: Welch for 2 groups; ANOVA for >2
    levels = list(scores_with_treatment["treatment"].dropna().unique())
    tukey_df = pd.DataFrame()
    if len(levels) == 2:
        if not control_label:
            control_label = _infer_control_label(levels)
        stats = _welch_stats(scores_with_treatment, treatment_col="treatment", control_label=control_label)
    elif len(levels) > 2:
        stats = _anova_stats(scores_with_treatment, treatment_col="treatment")
        # Tukey post-hoc for multi-group comparisons
        tukey_df = _tukey_posthoc(scores_with_treatment, treatment_col="treatment")
    else:
        stats = pd.DataFrame()
    group_stats = _compute_group_stats(scores_with_treatment, treatment_col="treatment")

    with PdfPages(output_pdf_path) as pdf:
        # Page 1: Methods
        fig = plt.figure(figsize=(8.5, 11))
        plt.axis("off")
        title = "Weighted Pathway Scoring — Bulk RNA-seq"
        default_note = (
            "Score = sum_i w_i * expr_i. Weights are normalized over detected genes.\n"
            "Input: normalized gene×sample matrix (TPM/CPM/VST). Gene matching is case-insensitive; duplicates summed.\n"
            "If provided, treatment groups are used for summaries and Welch's t-test (case minus control)."
        )
        text = (methods_note or default_note)
        plt.text(0.05, 0.95, title, fontsize=16, va="top", ha="left", weight="bold")
        plt.text(0.05, 0.90, text, fontsize=11, va="top", ha="left")
        if not stats.empty:
            # Small stats preview
            if "test" in stats.columns and (stats["test"] == "anova").all():
                cols = [c for c in ["pathway", "k", "F", "p", "q_BH"] if c in stats.columns]
            else:
                cols = [c for c in ["pathway", "group_a", "group_b", "n_a", "n_b", "p", "q_BH", "hedges_g"] if c in stats.columns]
            prev = stats[cols].sort_values([cols[0]]).head(12)
            plt.text(0.05, 0.60, "Stats preview (top 12):", fontsize=12, weight="bold")
            ax = plt.axes([0.05, 0.15, 0.9, 0.40])
            ax.axis("off")
            ax.table(cellText=prev.values, colLabels=list(prev.columns), loc="center")
        pdf.savefig(fig); plt.close(fig)

        # Page 2: Scores table (wide)
        fig = plt.figure(figsize=(11, 8.5))
        ax = plt.gca(); ax.axis("off")
        ax.set_title("Scores (per sample) — wide table", fontsize=14, loc="left")
        # Limit the table rows to fit; if too many, show head and note
        max_rows = 30
        show = tbl if len(tbl) <= max_rows else pd.concat([tbl.head(max_rows - 1), pd.DataFrame([["…"] * len(tbl.columns)], columns=tbl.columns)])
        the_table = ax.table(cellText=show.values, colLabels=list(show.columns), loc="center")
        the_table.auto_set_font_size(False); the_table.set_fontsize(8)
        the_table.scale(1, 1.2)
        pdf.savefig(fig); plt.close(fig)

        # Pages: Bar plots per pathway
        pws = list(scores_with_treatment["pathway"].unique())
        for pw in pws:
            sub = scores_with_treatment[scores_with_treatment["pathway"] == pw]
            gs = group_stats[group_stats["pathway"] == pw]
            fig = plt.figure(figsize=(8.5, 6))
            ax = plt.gca()
            if style == "box":
                # Boxplot with jitter
                import seaborn as sns  # seaborn improves boxplots if available
                try:
                    sns.boxplot(data=sub, x="treatment", y="score_w", ax=ax, color="#DDDDDD")
                    sns.stripplot(data=sub, x="treatment", y="score_w", ax=ax, color="#4C78A8", size=5, jitter=0.2)
                except Exception:
                    # Fallback to matplotlib
                    labels = list(gs["treatment"]) if not gs.empty else sorted(sub["treatment"].dropna().unique())
                    data = [sub.loc[sub["treatment"] == t, "score_w"].to_numpy() for t in labels]
                    ax.boxplot(data, labels=labels)
            else:
                # bars with sem error bars
                x = np.arange(len(gs))
                ax.bar(x, gs["mean"], yerr=gs["sem"], capsize=4, color="#4C78A8")
                ax.set_xticks(x)
                ax.set_xticklabels(gs["treatment"], rotation=45, ha="right")
            ax.set_ylabel("Score")
            ax.set_title(f"{pw} by treatment")
            # annotate p if Welch available
            st = pd.DataFrame()
            if isinstance(stats, pd.DataFrame) and ("pathway" in stats.columns):
                st = stats[stats["pathway"] == pw]
            if isinstance(st, pd.DataFrame) and not st.empty:
                p = float(st["p"].values[0])
                q = float(st["q_BH"].values[0]) if "q_BH" in st.columns else np.nan
                test = str(st["test"].values[0]) if "test" in st.columns else ""
                if test == "anova":
                    txt = f"ANOVA p={p:.3g}"
                else:
                    txt = f"Welch t-test p={p:.3g}"
                if np.isfinite(q):
                    txt += f", q={q:.3g}"
                ax.text(0.98, 0.92, txt, transform=ax.transAxes, ha="right", va="top", fontsize=10,
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.6))
            pdf.savefig(fig); plt.close(fig)

    # Save auxiliary CSVs next to PDF for programmatic access
    base = os.path.dirname(output_pdf_path) or "."
    _as_wide_table(scores_with_treatment).to_csv(os.path.join(base, "scores_bulk_wide.csv"), index=False)
    gs = _compute_group_stats(scores_with_treatment, "treatment")
    gs.to_csv(os.path.join(base, "scores_bulk_by_treatment.csv"), index=False)
    # Wide by treatment (means only)
    wide_trt = gs.pivot_table(index="treatment", columns="pathway", values="mean")
    wide_trt.reset_index().to_csv(os.path.join(base, "scores_bulk_by_treatment_wide.csv"), index=False)
    # Persist whichever stats were computed
    st = stats.copy()
    st.to_csv(os.path.join(base, "stats_bulk.csv"), index=False)
    if not tukey_df.empty:
        tukey_df.to_csv(os.path.join(base, "stats_bulk_tukey.csv"), index=False)
