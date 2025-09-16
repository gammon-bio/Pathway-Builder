#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad


def main() -> None:
    genes = [
        "Ngfr",
        "Sort1",
        "Sorcs2",
        "Ntrk2",
        "Acta1",
        "Gapdh",
        "GeneX",
        "GeneY",
        "GeneZ",
        "Rela",
    ]
    n = 40
    X = np.random.poisson(lam=1.5, size=(n, len(genes)))
    cond = np.array(["Control"] * 20 + ["Case"] * 20)
    celltype = np.array(["Myofiber"] * 20 + ["Satellite"] * 20)
    obs = pd.DataFrame({"condition": cond, "celltype_guess": celltype})
    var = pd.DataFrame(index=pd.Index(genes, name=None))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    import os

    os.makedirs("data/sc_toy", exist_ok=True)
    out = "data/sc_toy/sc.h5ad"
    adata.write(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

