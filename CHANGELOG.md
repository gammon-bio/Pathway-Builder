# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] - 2025-09-16
### Added
- Bulk and single-cell weighted pathway scoring.
- Makefile workflows for bulk (`toy-bulk`, `bulk`, `check`) and single-cell (`sc`, `sc-check`, `sc-zhang`, `sc-zhang-with-info`).
- Sample info support for single-cell condition labels.
- Comparison outputs: tidy per-celltype means and deltas; PDFs for means/deltas; delta page in main report.
- Curated TrkB example list and toy bulk dataset.
- Unified LLM prompt (six-column schema) and documentation.
- CI matrix for Python 3.9–3.12 running tests and toy-bulk smoke.

### Changed
- Removed legacy user data and redundant scripts.
- Preserved string sample IDs during AnnData concatenation.

### Fixed
- Faster duplicate-gene collapse in single-cell loader; avoided sparse slicing hangs.
