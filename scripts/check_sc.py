#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = args.out
    ct_path = os.path.join(out, "scores_sc_by_celltype.csv")
    if not os.path.exists(ct_path):
        print(f"[error] missing {ct_path}", file=sys.stderr)
        return 2
    df = pd.read_csv(ct_path)
    n_ct = df["celltype"].nunique() if "celltype" in df.columns else 0
    n_pw = df["pathway"].nunique() if "pathway" in df.columns else 0
    print(f"[ok] Cell-type summary: {n_ct} celltypes, {n_pw} pathways")
    # Show small pivot of means if plausible
    try:
        piv = df.pivot_table(index=["celltype","condition"], columns="pathway", values="score_w")
        print("\nCelltype x condition means (head):")
        print(piv.head(10))
    except Exception:
        pass
    dpath = os.path.join(out, "deltas.csv")
    if os.path.exists(dpath):
        d = pd.read_csv(dpath)
        print(f"\nDeltas rows: {len(d)}")
        print(d.head(10))
    cells_path = os.path.join(out, "scores_sc_cells.csv")
    if os.path.exists(cells_path):
        cells = pd.read_csv(cells_path)
        print(f"\nCell-level scores rows: {len(cells)} (showing columns)")
        print(list(cells.columns))
    figs = os.path.join(out, "figs")
    if os.path.isdir(figs):
        print(f"\nFigures: {len(os.listdir(figs))} files in figs/")
    pdf = os.path.join(out, "report.pdf")
    if os.path.exists(pdf):
        print("PDF report found.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

