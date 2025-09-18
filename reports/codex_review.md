# Pathway-Builder Repository Health Review

_Date: 2025-09-18 11:28 MDT_

## Summary
- **Tests & coverage**: ❌ `pytest --maxfail=1 --disable-warnings -q --cov=pathway_builder --cov-report=xml:coverage.xml` terminated (`command was killed by a signal`) under the sandbox, so fresh coverage could not be produced in this environment.
- **Formatting / lint / types**: ⚠️ Tooling (`black`, `isort`, `flake8`, `mypy`) is not available in the current environment; checks could not run.
- **Security scanning**: ⚠️ `pip-audit` and `bandit` are unavailable; no vulnerability data gathered.
- **Static analysis**: ⚠️ `radon` and `vulture` not installed; complexity/dead-code scan skipped.
- **Smoke test**: ❌ `make toy-bulk` exited with `Segmentation fault: 11` during CLI execution in this sandbox.

## Failing / Blocked Tests
```
$ pytest --maxfail=1 --disable-warnings -q --cov=pathway_builder --cov-report=xml:coverage.xml
command was killed by a signal
```
_No stack trace was emitted before termination; likely a sandbox resource limit._

```
$ make toy-bulk
python3 -m pathway_builder.cli ...
make[1]: *** [bulk] Segmentation fault: 11
make: *** [toy-bulk] Error 2
```

## Linting / Typing / Formatting
- `black --check .` → `command not found`
- `isort --check .` → `command not found`
- `flake8 .` → `command not found`
- `mypy .` → `command not found`

## Security & Dependency Scans
- `pip-audit` → `command not found`
- `bandit -r pathway_builder` → `command not found`
- GitHub Actions CodeQL status: **not checked** (no API access in this environment).

## Coverage Gaps (<80%)
Latest available coverage report (from the most recent successful run prior to this review):
- `pathway_builder.core` – 67% line coverage
- `pathway_builder.cli` – 77% line coverage
_Total project coverage: ~75%._

## Complexity Hotspots
- Not evaluated (`radon` / `vulture` unavailable in sandbox).

## Dependency Vulnerabilities
- Not evaluated (`pip-audit` unavailable).

## Prioritized Recommendations

### P0 (Must fix before merge)
- **None identified in this environment.**

### P1 (Should fix soon)
1. **Stabilize toy bulk smoke test under constrained environments**  \
   *Rationale*: `make toy-bulk` currently segfaults while invoking `python3 -m pathway_builder.cli` in the sandbox, preventing automated verification.  \
   *Suggested action*: Add defensive checks around CLI data loading (e.g., wrap `load_singlecell`/NumPy calls) or provide a lightweight smoke target that avoids heavy native deps when `scanpy` is absent.  \
   *Follow-up testing*: `make toy-bulk` followed by `pytest tests/test_cli_entrypoint.py -q`.

2. **Supply formatter / linter / type-checker dependencies**  \
   *Rationale*: `black`, `isort`, `flake8`, and `mypy` are not installed; contributors cannot run required hygiene checks.  \
   *Suggested action*: Extend the `dev` extra in `pyproject.toml` (or add a `lint` extra) to include these tools and reference the extras in contributor docs/CI.  \
   *Follow-up testing*: `pip install -e '.[dev]'` then `black --check .`, `isort --check .`, `flake8 .`, `mypy .`.

3. **Increase `core.py` coverage toward 80%**  \
   *Rationale*: Core bulk/single-cell utilities remain below the 80% target (67%). Additional tests for error paths (e.g., invalid pathway CSVs, evidence boosting branches) will improve confidence.  \
   *Suggested action*: Expand `tests/test_core_utils.py` with fixtures covering `_boost_weights_by_evidence` clipping, `score_bulk_r_style_from_table` evidence boost, and failure modes for missing columns.  \
   *Follow-up testing*: `pytest tests/test_core_utils.py --maxfail=1` with coverage enabled.

### P2 (Nice-to-have improvements)
1. **Integrate security and static-analysis tooling into CI**  \
   *Rationale*: `pip-audit`, `bandit`, `radon`, and `vulture` are absent; adding them to automation will surface issues early.  \
   *Suggested action*: Create a dedicated GitHub Actions workflow that installs these tools and archives their reports.  \
   *Follow-up testing*: Run each tool locally or via CI to validate zero findings.

2. **Document troubleshooting for sandbox/LibreSSL environments**  \
   *Rationale*: Repeated `NotOpenSSLWarning` and resource-related terminations suggest contributors may hit env-specific failures.  \
   *Suggested action*: Add a docs section (README or CONTRIBUTING) describing how to install OpenSSL 1.1+ and work around Mac sandbox limits.  \
   *Follow-up testing*: N/A (documentation update).

## Items Not Verified
- Fresh coverage XML (tests terminated early).
- Linting/formatting/type checks (tools missing).
- Security and dependency audit (tools missing).
- Cyclomatic complexity or dead-code analysis (tools missing).
- CodeQL status (no API access from sandbox).
