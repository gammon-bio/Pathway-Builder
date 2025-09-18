from __future__ import annotations

import math
import numpy as np
import pandas as pd
import pytest

from pathway_builder import core


def test_read_pathway_csv_merges_duplicates_and_preserves_evidence(tmp_path):
    csv_path = tmp_path / "pathway.csv"
    csv_path.write_text(
        "gene,weight,evidence_source\n"
        "GeneA,0.5,KEGG\n"
        "genea,0.5,Reactome\n"
        "GeneB,,\n"
    )

    df = core.read_pathway_csv(str(csv_path))
    assert list(df["gene"]) == ["GeneA", "GeneB"]
    # Duplicate weights summed, GeneB defaults to weight 1.0
    assert math.isclose(df.loc[df["gene"] == "GeneA", "weight"].iloc[0], 1.0)
    assert math.isclose(df.loc[df["gene"] == "GeneB", "weight"].iloc[0], 1.0)
    assert df.loc[df["gene"] == "GeneA", "evidence_source"].iloc[0] in {"KEGG", "Reactome"}


def test_normalize_weights_present_filters_and_normalizes():
    pw = pd.DataFrame({"gene": ["A", "B", "C"], "weight": [1.0, 2.0, 3.0]})
    present = ["A", "C"]

    norm = core.normalize_weights_present(pw, present_genes=present)
    assert list(norm["gene"]) == ["A", "C"]
    assert math.isclose(norm["weight"].sum(), 1.0)
    assert math.isclose(norm.loc[norm["gene"] == "C", "weight"].iloc[0], 0.75)


def test_collapse_duplicates_sum():
    df = pd.DataFrame(
        {
            "gene": ["X", "x", "Y"],
            "s1": [1, 2, 3],
            "s2": [4, 6, 8],
        }
    )
    collapsed = core._collapse_duplicates_sum(df, gene_col="gene")
    assert list(collapsed["gene"]) == ["X", "Y"]
    assert collapsed.loc[0, "s1"] == 3
    assert collapsed.loc[0, "s2"] == 10


def test_load_vst_counts_table_duplicate_handling():
    df = pd.DataFrame(
        {
            "gene": ["A", "A", "B"],
            "sample1": [1.0, 5.0, 2.0],
            "sample2": [2.0, 6.0, 4.0],
        }
    )

    mean = core.load_vst_counts_table(df, collapse="mean")
    assert list(mean.index) == ["A", "B"]
    assert math.isclose(mean.loc["A", "sample1"], 3.0)

    median = core.load_vst_counts_table(df, collapse="median")
    assert math.isclose(median.loc["A", "sample1"], 3.0)

    summed = core.load_vst_counts_table(df, collapse="sum")
    assert math.isclose(summed.loc["A", "sample1"], 6.0)

    maxvar = core.load_vst_counts_table(df, collapse="maxvar")
    # The second row has higher variance, so ensure its values are present
    expected = df.iloc[1][["sample1", "sample2"]]
    selected = maxvar.loc["A"]
    if isinstance(selected, pd.Series):
        selected_rows = selected.to_frame().T
    else:
        selected_rows = selected
    match_mask = (selected_rows[["sample1", "sample2"]] == expected.values).all(axis=1)
    assert match_mask.any()


def test_score_bulk_r_style_handles_label_padding(tmp_path):
    counts = pd.DataFrame(
        {
            "gene": ["A", "B"],
            "S1": [1.0, 2.0],
            "S2": [2.0, 4.0],
        }
    )
    pw1 = tmp_path / "first.csv"
    pw1.write_text("gene,weight\nA,1.0\n")
    pw2 = tmp_path / "second.csv"
    pw2.write_text("gene,weight\nC,1.0\n")

    scores = core.score_bulk_r_style_from_table(
        counts,
        gene_col="gene",
        pathway_csvs=[str(pw1), str(pw2)],
        labels=["FIRST"],
    )

    first = scores[scores["pathway"] == "FIRST"].sort_values("sample")
    assert list(first["sample"]) == ["S1", "S2"]
    assert first["score_w"].tolist() == pytest.approx([-0.7071, 0.7071], rel=1e-3)

    second = scores[scores["pathway"] == "second"].sort_values("sample")
    assert list(second["sample"]) == ["S1", "S2"]
    assert second["score_w"].tolist() == [0.0, 0.0]


def test_score_bulk_from_table_requires_gene_column(tmp_path):
    counts = pd.DataFrame(
        {
            "Symbol": ["G1", "G2"],
            "S1": [1, 2],
            "S2": [3, 4],
        }
    )
    pw = tmp_path / "pw.csv"
    pw.write_text("gene,weight\nG1,1\n")

    with pytest.raises(ValueError):
        core.score_bulk_from_table(
            counts_table=counts,
            gene_col="gene",
            pathway_csvs=[str(pw)],
        )
    try:
        core.load_vst_counts_table(df, collapse="unknown")
    except ValueError as exc:
        assert "collapse" in str(exc)
    else:
        raise AssertionError("Expected ValueError for unknown collapse option")


def test_boost_weights_by_evidence_clips_values():
    pw = pd.DataFrame(
        {
            "gene": ["A", "B"],
            "weight": [1.1, 1.1],
            "evidence_source": ["KEGG", "Reactome"],
        }
    )
    boosted = core._boost_weights_by_evidence(pw, boost=True)
    assert boosted["weight"].max() <= 1.2
    unboosted = core._boost_weights_by_evidence(pw, boost=False)
    assert np.allclose(unboosted["weight"], pw["weight"])


def test_welch_and_hedges_metrics():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    y = np.array([1.5, 2.5, 3.5, 4.5])

    t, df, p, diff, ci_low, ci_high = core.welch_t_test(x, y)
    assert math.isfinite(t)
    assert math.isfinite(df)
    assert math.isfinite(p)
    assert math.isfinite(diff)
    assert math.isfinite(ci_low)
    assert math.isfinite(ci_high)

    g = core.hedges_g(x, y)
    assert math.isfinite(g)
    assert not math.isclose(g, 0.0)


def test_is_tsv_or_csv_variants(tmp_path):
    base = tmp_path / "file"
    for ext in [".csv", ".tsv", ".txt"]:
        path = base.with_suffix(ext)
        path.write_text("dummy")
        assert core._is_tsv_or_csv(str(path))

    assert not core._is_tsv_or_csv(str(base.with_suffix(".xlsx")))


def test_read_delim_respects_extension(tmp_path):
    csv_path = tmp_path / "table.csv"
    tsv_path = tmp_path / "table.tsv"
    csv_path.write_text("a,b\n1,2\n")
    tsv_path.write_text("a\tb\n3\t4\n")

    csv_df = core._read_delim(str(csv_path))
    tsv_df = core._read_delim(str(tsv_path))

    assert csv_df.iloc[0, 0] == 1
    assert tsv_df.iloc[0, 0] == 3
