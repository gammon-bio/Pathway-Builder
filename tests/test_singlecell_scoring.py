import importlib
import numpy as np
import pandas as pd
import pytest


@pytest.mark.filterwarnings("ignore:.*")
def test_singlecell_weighted_scoring_small(tmp_path):
    scanpy = pytest.importorskip("scanpy")
    anndata = importlib.import_module("anndata")
    from pathway_builder.core import score_singlecell_adata, read_pathway_csv

    # 3 cells x 3 genes, simple counts
    X = np.array(
        [
            [10, 0, 0],  # cell0 expresses A
            [0, 10, 0],  # cell1 expresses B
            [0, 0, 10],  # cell2 expresses C
        ],
        dtype=float,
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])  # gene names
    obs = pd.DataFrame(index=["c0", "c1", "c2"]).assign(celltype=["T1", "T1", "T2"])  # two groups
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Pathway favors GeneA only
    pw = pd.DataFrame({"gene": ["genea"], "weight": [1.0]})
    pw_path = tmp_path / "pw.csv"
    pw.to_csv(pw_path, index=False)

    try:
        summary = score_singlecell_adata(
            adata,
            pathway_csvs=[str(pw_path)],
            labels=["PW"],
            celltype_col="celltype",
            condition_col=None,
        )
    finally:
        pass

    # The T1 group (cells 0 and 1) includes a cell with GeneA expression; after normalize/log1p,
    # mean score should be greater than T2 (only cell with GeneC)
    t1 = float(summary.loc[summary["celltype"] == "T1", "score_w__PW"].values[0])
    t2 = float(summary.loc[summary["celltype"] == "T2", "score_w__PW"].values[0])
    assert t1 > t2
