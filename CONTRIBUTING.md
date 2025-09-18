# Contributing

Thanks for your interest in improving Pathway Builder! A few guidelines:

- Issues: Please include steps to reproduce and environment details (OS, Python, package versions).
- Pull requests: Keep changes focused and include a short description of the rationale.
- Style: Follow existing code style; keep functions small and testable.
- Tests: Add tests for new features where applicable (see `pathway-builder-starter/tests/`).
- Large data: Do not commit large datasets or outputs. Use the Makefile to generate outputs locally.

## Development setup

- Install dev dependencies: `pip install -e '.[dev]'` (quote extras in zsh to avoid glob expansion)
- Optional: add single-cell extras with `make setup-sc`
- Run the bundled smoke test once: `make toy-bulk` (after installing `.[dev]` extras)

## Local testing

- Unit tests with coverage (writes `coverage.xml`):
  - `pytest --cov=pathway_builder --cov-report=xml:coverage.xml --cov-report=term-missing`
- Regenerate the toy outputs after changes that touch scoring: `make toy-bulk`
- Lint/type checks: `black --check .`, `isort --check .`, `flake8 .`, `mypy .`

## Security & static analysis

- Dependency audit: `pip-audit`
- Security lint: `bandit -r pathway_builder`
- Complexity / dead code: `radon cc pathway_builder -s`, `vulture pathway_builder`

## Branch naming

- Feature work: `feature/<short-desc>`
- Bug fixes: `bugfix/<short-desc>`
- Automation chores: `automation/<short-desc>`

## Pull request checklist

- ✅ Tests updated or added
- ✅ CHANGELOG entry when behavior changes
- ✅ Docs/tutorials refreshed if user-facing behavior shifts
- ✅ CI passes (`coverage.xml` uploaded via GitHub Actions)

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
