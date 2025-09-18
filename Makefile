SHELL := /bin/bash
.ONESHELL:

# Defaults (override on invocation)
PYTHON ?= python3
BULK_PY ?= $(PYTHON)
SMOKE_VENV ?= .venv
SMOKE_PY := $(SMOKE_VENV)/bin/python3
CLI ?= $(BULK_PY) -m pathway_builder.cli
SCORING_STYLE ?= r
COLLAPSE ?= mean
REPORT_STYLE ?= box
CONTROL ?=

.PHONY: help setup setup-sc bulk toy-bulk check sc sc-zhang sc-check sc-toy

help:
	@echo "Targets:"
	@echo "  setup             Install package with report extras"
	@echo "  bulk              Run bulk scoring (set variables below)"
	@echo "  toy-bulk          Run toy bulk example and write to out/toy-bulk"
	@echo "  check             Validate outputs and print a short summary (OUT=...)"
	@echo "  sc                Run single-cell scoring (SN_DATA=..., PATHWAYS=..., OUT=...)"
	@echo "  sc-zhang          Run Zhang 2024 pipeline (downloads if needed)"
	@echo "  sc-check          Validate single-cell outputs (OUT=...)"
	@echo "  setup-sc          Install single-cell + report extras"
	@echo "  sc-toy            Generate a tiny synthetic .h5ad and run sc (OUT=out/sc-toy)"
	@echo "  sc-zhang-with-info Run sc on Zhang .h5s using bundled sample info"
	@echo ""
	@echo "Variables for 'make bulk':"
	@echo "  COUNTS=path/to/counts.tsv (required)"
	@echo "  PATHWAYS='path/a.csv path/b.csv' (one or more)"
	@echo "  LABELS='LABEL_A LABEL_B'      (optional; match PATHWAYS order)"
	@echo "  SAMPLE_INFO=path/to/sample_info.csv (optional)"
	@echo "  OUT=out/my_run (required)"
	@echo "  CONTROL=Control (optional; control label)"
	@echo "  SCORING_STYLE=r|simple (default: r)"
	@echo "  COLLAPSE=mean|median|sum|maxvar (default: mean)"
	@echo "  REPORT_STYLE=box|bar (default: box)"

setup:
	$(PYTHON) -m pip install -e '.[report]'

setup-sc:
	$(PYTHON) -m pip install -e '.[singlecell,report]'

# Helper to expand repeated flags for PATHWAYS and LABELS
define map_pathways
$(foreach p,$(PATHWAYS),--pathway_csv $(p))
endef
define map_labels
$(foreach l,$(LABELS),--label $(l))
endef

bulk:
	@if [ -z "$(COUNTS)" ] || [ -z "$(OUT)" ] || [ -z "$(PATHWAYS)" ]; then \
		echo "Usage: make bulk COUNTS=... PATHWAYS='a.csv b.csv' OUT=out/my_run [LABELS='A B'] [SAMPLE_INFO=...]"; \
		exit 2; \
	fi
	mkdir -p "$(OUT)"
	@$(BULK_PY) scripts/check_smoke_deps.py
	$(CLI) \
	  --bulk_counts "$(COUNTS)" \
	  $(call map_pathways) \
	  $(call map_labels) \
	  $(if $(SAMPLE_INFO),--sample_info "$(SAMPLE_INFO)",) \
	  $(if $(CONTROL),--control_label "$(CONTROL)",) \
	  --scoring_style "$(SCORING_STYLE)" \
	  --collapse_duplicates "$(COLLAPSE)" \
	  --report_style "$(REPORT_STYLE)" \
	  $(if $(NO_PDF),--no_pdf,) \
	  --output_dir "$(OUT)"

$(SMOKE_PY):
	$(PYTHON) -m venv $(SMOKE_VENV)
	. $(SMOKE_VENV)/bin/activate
	python3 -m pip install --upgrade pip
	python3 -m pip install -e '.[dev]'

toy-bulk: $(SMOKE_PY)
	$(MAKE) bulk \
	  PYTHON=$(SMOKE_PY) \
	  BULK_PY=$(SMOKE_PY) \
	  CLI='$(SMOKE_PY) -m pathway_builder.cli' \
	  COUNTS=data/toy_bulk_counts.tsv \
	  PATHWAYS='genes/toy_pathway.csv' \
	  LABELS='TOY' \
	  SAMPLE_INFO=data/toy_sample_info.csv \
	  OUT=out/toy-bulk \
	  SCORING_STYLE=simple \
	  COLLAPSE=mean \
	  REPORT_STYLE=box \
	  NO_PDF=1

check:
	@if [ -z "$(OUT)" ]; then echo "Usage: make check OUT=out/my_run"; exit 2; fi
	@if [ ! -d "$(OUT)" ]; then echo "Output dir not found: $(OUT)"; exit 1; fi
	@test -f "$(OUT)/scores_bulk.csv" || { echo "[error] missing $(OUT)/scores_bulk.csv"; exit 1; }
	@if [ ! -f "$(OUT)/report_bulk.pdf" ]; then echo "[warn] $(OUT)/report_bulk.pdf not found (PDF may be disabled)"; fi
	python3 scripts/check_bulk.py --out "$(OUT)"
sc:
	@if [ -z "$(SN_DATA)" ] || [ -z "$(OUT)" ] || [ -z "$(PATHWAYS)" ]; then \
		echo "Usage: make sc SN_DATA='data/*.h5' PATHWAYS='genes/pw.csv ...' [LABELS='A B'] OUT=out/sc_run [CELLTYPE=celltype] [CONDITION=condition]"; \
		exit 2; \
	fi
	mkdir -p "$(OUT)"
	python3 run_bdnf_enrichment.py \
	  --sn_data "$(SN_DATA)" \
	  $(call map_pathways) \
	  $(call map_labels) \
	  $(if $(CELLTYPE),--celltype_col "$(CELLTYPE)",) \
	  $(if $(CONDITION),--condition_col "$(CONDITION)",) \
	  $(if $(SAMPLE_INFO),--sample_info "$(SAMPLE_INFO)" $(if $(SAMPLE_COL),--sample_col "$(SAMPLE_COL)",) $(if $(INFO_CONDITION_COL),--info_condition_col "$(INFO_CONDITION_COL)",),) \
	  --output_dir "$(OUT)"

sc-zhang:
	bash scripts/run_bdnf_zhang2024.sh

sc-check:
	@if [ -z "$(OUT)" ]; then echo "Usage: make sc-check OUT=out/sc_run"; exit 2; fi
	@test -f "$(OUT)/scores_sc_by_celltype.csv" || { echo "[error] missing $(OUT)/scores_sc_by_celltype.csv"; exit 1; }
	python3 scripts/check_sc.py --out "$(OUT)"

sc-zhang-with-info:
	$(MAKE) sc \
	  SN_DATA='data/GSE272085/*.h5' \
	  PATHWAYS='genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv' \
	  LABELS='TRKB' \
	  SAMPLE_INFO=data/zhang_sample_info.csv \
	  OUT=out/zhang2024_with_info

sc-toy:
	mkdir -p data/sc_toy out/sc-toy
	python3 scripts/gen_sc_toy.py
	$(MAKE) sc \
	  SN_DATA='data/sc_toy/sc.h5ad' \
	  PATHWAYS='genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv' \
	  LABELS='TRKB' \
	  CELLTYPE=celltype_guess \
	  CONDITION=condition \
	  OUT=out/sc-toy
	$(MAKE) sc-check OUT=out/sc-toy
