from __future__ import annotations

from pathlib import Path

import pandas as pd

from pathway_builder.cli import main


def test_cli_bulk_runs_with_toy_data(tmp_path):
    out_dir = tmp_path / "toy-out"
    args = [
        "--bulk_counts",
        str(Path("data/toy_bulk_counts.tsv")),
        "--pathway_csv",
        str(Path("genes/toy_pathway.csv")),
        "--label",
        "TOY",
        "--sample_info",
        str(Path("data/toy_sample_info.csv")),
        "--no_pdf",
        "--scoring_style",
        "simple",
        "--output_dir",
        str(out_dir),
    ]

    main(args)

    scores_path = out_dir / "scores_bulk.csv"
    assert scores_path.exists()
    scores = pd.read_csv(scores_path)
    assert not scores.empty


def test_cli_emits_warning_for_missing_genes(tmp_path, capsys):
    counts = tmp_path / "counts.csv"
    counts.write_text("gene,sample1\nFOO,1\n")
    pathway = tmp_path / "pathway.csv"
    pathway.write_text("gene,weight\nBAR,1\n")
    out_dir = tmp_path / "out"

    args = [
        "--bulk_counts",
        str(counts),
        "--pathway_csv",
        str(pathway),
        "--scoring_style",
        "simple",
        "--output_dir",
        str(out_dir),
    ]

    main(args)

    captured = capsys.readouterr()
    assert "No pathway genes matched input counts" in captured.err
    assert (out_dir / "scores_bulk.csv").exists()
