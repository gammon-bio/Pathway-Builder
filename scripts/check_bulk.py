#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Output directory from bulk run")
    args = ap.parse_args()
    out = args.out

    def load(name: str) -> pd.DataFrame:
        p = os.path.join(out, name)
        return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()

    scores = load("scores_bulk.csv")
    if scores.empty:
        print(f"[error] {os.path.join(out, 'scores_bulk.csv')} not found or empty", file=sys.stderr)
        return 2
    n_samp = scores["sample"].nunique()
    n_pw = scores["pathway"].nunique()
    print(f"[ok] Found: {n_samp} samples, {n_pw} pathways")

    bytr = load("scores_bulk_by_treatment.csv")
    if not bytr.empty and {"treatment", "pathway", "mean"}.issubset(bytr.columns):
        piv = bytr.pivot_table(index="treatment", columns="pathway", values="mean")
        print("\nGroup means (treatment x pathway):")
        print(piv)

    stats = load("stats_bulk.csv")
    if not stats.empty:
        cols = [c for c in ["pathway", "test", "p", "q_BH", "hedges_g", "F"] if c in stats.columns]
        if cols:
            print("\nStats (first 6 rows):")
            print(stats[cols].head(6))

    tk = os.path.join(out, "stats_bulk_tukey.csv")
    if os.path.exists(tk):
        tdf = pd.read_csv(tk)
        print(f"\nTukey pairs: {len(tdf)} rows (showing first 6)")
        cols = [c for c in ["pathway", "group1", "group2", "meandiff", "p_adj"] if c in tdf.columns]
        print(tdf[cols].head(6))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

