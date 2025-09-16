import os
import pandas as pd

from pathway_builder.core import score_bulk_from_table


def test_bulk_weighted_scoring_simple(tmp_path):
    # Gene names with mixed case to test case-insensitive matching
    df = pd.DataFrame(
        {
            "Gene": ["GeneA", "GENEB", "genec"],
            "Sample1": [10, 20, 30],
            "Sample2": [5, 0, 5],
        }
    )

    # Weight pathway with A and C equally weighted
    pw = pd.DataFrame({"gene": ["genea", "GENEC"], "weight": [0.5, 0.5]})
    pw_path = tmp_path / "pw.csv"
    pw.to_csv(pw_path, index=False)

    out = score_bulk_from_table(
        counts_table=df,
        gene_col=None,  # auto-detect first column
        pathway_csvs=[str(pw_path)],
        labels=["MY_PW"],
    )
    # Expected: S1 = 0.5*10 + 0.5*30 = 20; S2 = 0.5*5 + 0.5*5 = 5
    out = out.sort_values(["pathway", "sample"]).reset_index(drop=True)
    assert list(out["pathway"]) == ["MY_PW", "MY_PW"]
    assert list(out["sample"]) == ["Sample1", "Sample2"]
    assert out.loc[0, "score_w"] == 20.0
    assert out.loc[1, "score_w"] == 5.0

