# Contributing

Thanks for your interest in improving Pathway Builder! A few guidelines:

- Issues: Please include steps to reproduce and environment details (OS, Python, package versions).
- Pull requests: Keep changes focused and include a short description of the rationale.
- Style: Follow existing code style; keep functions small and testable.
- Tests: Add tests for new features where applicable (see `pathway-builder-starter/tests/`).
- Large data: Do not commit large datasets or outputs. Use the Makefile to generate outputs locally.

## Development setup

```
make setup          # bulk/report extras
make setup-sc       # single-cell + report extras
pytest pathway-builder-starter/tests -q
```

## Pathway CSV schema

Pathway CSVs should have exactly these columns:

```
gene_symbol,gene_name,module,annotation,evidence_source,weight
```

- `gene_symbol`: official symbol for the target species
- `gene_name`: descriptive full name
- `module`: role (ligand, receptor, co-receptor, adaptor, kinase, TF, etc.)
- `annotation`: short phrase describing function/role in the pathway
- `evidence_source`: citations (e.g., KEGG/Reactome/reviews)
- `weight`: relative importance (>0), typically 0.3–1.0 (core nodes ≈1.0)

## Code of Conduct

Be respectful and constructive. We welcome feedback and contributions from researchers and developers.
