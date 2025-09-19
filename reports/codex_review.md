# Pathway-Builder Repository Health Review

_Date: 2025-09-18 12:06 MDT_

## Summary
- **Tests & coverage**: ✅ Host run (`pytest --cov=pathway_builder --cov-report=xml:coverage.xml --cov-report=term-missing`) passes all 17 tests in ~3s with total coverage 79.7%; pytest gate raised to 75%. Sandbox still terminates early, so CI should rely on host/runner results.
- **Formatting / lint / types**: ⚠️ Tooling still unavailable in the sandbox, but `pip install -e '.[dev]'` now installs `black`, `isort`, `flake8`, `mypy`, plus security/static tools for local or CI use.
- **Security scanning**: ⚠️ `pip-audit` and `bandit` missing in sandbox; available via the updated dev extras.
- **Static analysis**: ⚠️ `radon`/`vulture` not present in sandbox; installable through the dev extras.
- **Smoke test**: ⚠️ Sandbox `make toy-bulk` still segfaults; the target now provisions a `.venv`, installs `.[dev]`, runs simple scoring with `--no_pdf`, and checks numpy/pandas first. Confirm on CI runners.

## Test & Coverage Notes
```
# sandbox attempt
$ pytest --maxfail=1 --disable-warnings -q --cov=pathway_builder --cov-report=xml:coverage.xml
command was killed by a signal
```

```
# host run (Caleb)
$ pytest --cov=pathway_builder --cov-report=xml:coverage.xml --cov-report=term-missing
17 passed in 2.98s; total coverage 79.70%
```

```
# sandbox smoke test
$ make toy-bulk
make[1]: *** [bulk] Segmentation fault: 11
```

## Linting / Typing / Formatting
- `black --check .` → command not found
- `isort --check .` → command not found
- `flake8 .` → command not found
- `mypy .` → command not found

## Security & Dependency Scans
- `pip-audit` → command not found
- `bandit -r pathway_builder` → command not found
- GitHub Actions CodeQL status: not checked (no API access from sandbox)

## Coverage Gaps (<80%)
Latest successful host run (79.70% total):
- `pathway_builder.core` – 76% lines covered
- `pathway_builder.cli` – 77% lines covered

## Complexity Hotspots
- Not evaluated (missing `radon` / `vulture`).

## Dependency Vulnerabilities
- Not evaluated (`pip-audit` unavailable).

## Prioritized Recommendations

### P0 (Must fix before merge)
- None identified.

### P1 (Should fix soon)
1. **Ensure toy bulk smoke test runs in CI-quality environments**  \
   *Rationale*: Segfault in the sandbox still blocks automated validation; the target now provisions `.venv`, installs `.[dev]`, enforces simple scoring/`--no_pdf`, and checks dependencies, but CI runners need verification.  \
   *Suggested action*: Confirm on GitHub Actions; if failures persist, add dependency guards or a fallback smoke target.  \
   *Follow-up*: `make toy-bulk` then `pytest tests/test_cli_entrypoint.py -q`.

2. **Provide formatter / linter / type-checker dependencies**  \
   *Status*: ✅ `.[dev]` now installs `black`, `isort`, `flake8`, `mypy`, `pip-audit`, `bandit`, `radon`, `vulture`; docs list the lint/type/security commands.  \
   *Follow-up*: `pip install -e '.[dev]'` then run the lint/type/security commands locally or in CI.

3. **Continue pushing `core.py` coverage toward ≥80%**  \
   *Rationale*: New tests raised coverage to 76%, but critical helpers (`load_singlecell`, evidence boost branches) remain lightly exercised.  \
   *Suggested action*: Add paramized tests hitting single-cell paths with fakes/mocks.  \
   *Follow-up*: `pytest tests/test_core_utils.py --cov=pathway_builder.core`.

### P2 (Nice-to-have improvements)
1. **Integrate security and static-analysis tooling into CI**  \
   *Rationale*: Automating `pip-audit`, `bandit`, `radon`, and `vulture` will surface issues early.  \
   *Suggested action*: Add a dedicated GitHub Actions job installing these tools and archiving results.  \
   *Follow-up*: Run each tool locally or in CI.

2. **Document troubleshooting for LibreSSL/sandbox setups**  \
   *Rationale*: Persistent `NotOpenSSLWarning` and sandbox failures suggest contributors may encounter similar issues.  \
   *Suggested action*: Update README/CONTRIBUTING with OpenSSL upgrade guidance and sandbox workarounds.  \
   *Follow-up*: N/A (documentation only).

## Items Not Verified
- Linting/formatting/type checks
- Security and dependency audit
- Complexity/dead-code analysis
- CodeQL status
