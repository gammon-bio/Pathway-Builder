import os
import pandas as pd

from pathway_builder.core import score_bulk_from_table
from pathway_builder.report import make_bulk_pdf_report


def test_bulk_pdf_report(tmp_path):
    # Minimal dataset with two treatments and two samples each
    df = pd.DataFrame(
        {
            "Gene": ["A", "B", "C"],
            "S1": [10, 0, 0],
            "S2": [8, 0, 0],
            "S3": [0, 5, 0],
            "S4": [0, 6, 0],
        }
    )
    pw = pd.DataFrame({"gene": ["A", "B"], "weight": [1.0, 1.0]})
    pw_path = tmp_path / "pw.csv"
    pw.to_csv(pw_path, index=False)

    scores = score_bulk_from_table(df, gene_col=None, pathway_csvs=[str(pw_path)], labels=["PW"])
    # Attach sample info
    si = pd.DataFrame({"sample": ["S1", "S2", "S3", "S4"], "treatment": ["Ctrl", "Ctrl", "Case", "Case"]})
    scores = scores.merge(si, on="sample", how="left")

    pdf_path = tmp_path / "rep.pdf"
    make_bulk_pdf_report(scores_with_treatment=scores, output_pdf_path=str(pdf_path))
    assert pdf_path.exists() and pdf_path.stat().st_size > 0
    # Also check that helper CSVs were written alongside (by report function)
    assert (tmp_path / "scores_bulk_wide.csv").exists()
    assert (tmp_path / "scores_bulk_by_treatment.csv").exists()
    assert (tmp_path / "stats_bulk.csv").exists()


def test_bulk_anova_stats(tmp_path):
    import numpy as np
    # 3 groups with different means
    df = pd.DataFrame(
        {
            "Gene": ["A", "B"],
            "C1": [10, 0],
            "C2": [11, 0],
            "T1": [0, 5],
            "T2": [0, 6],
            "T3": [0, 8],
        }
    )
    pw = pd.DataFrame({"gene": ["A", "B"], "weight": [1.0, 1.0]})
    pw_path = tmp_path / "pw.csv"
    pw.to_csv(pw_path, index=False)
    scores = score_bulk_from_table(df, gene_col=None, pathway_csvs=[str(pw_path)], labels=["PW"])
    si = pd.DataFrame({
        "sample": ["C1", "C2", "T1", "T2", "T3"],
        "treatment": ["Control", "Control", "TreatA", "TreatA", "TreatB"],
    })
    scores = scores.merge(si, on="sample", how="left")
    pdf_path = tmp_path / "rep.pdf"
    make_bulk_pdf_report(scores_with_treatment=scores, output_pdf_path=str(pdf_path))
    stats = pd.read_csv(tmp_path / "stats_bulk.csv")
    assert not stats.empty
    assert "test" in stats.columns and set(stats["test"]) == {"anova"}
    assert float(stats.loc[0, "p"]) >= 0.0
    tuk = pd.read_csv(tmp_path / "stats_bulk_tukey.csv")
    assert not tuk.empty
    assert set(["pathway", "group1", "group2", "meandiff", "p_adj"]).issubset(tuk.columns)
