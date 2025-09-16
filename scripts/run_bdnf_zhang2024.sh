#!/usr/bin/env bash
set -euo pipefail

# --- Config ---
ENV_NAME="bdnf-snrna"
ENV_FILE="env/scanpy_env.yml"
# Manage env and installs with conda to avoid prefix mismatch issues
PKG=conda
RUN_TOOL="conda run"
DATA_DIR="./data/GSE272085"
OUT_DIR="./out/zhang2024"
GEO_TAR_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272085/suppl/GSE272085_RAW.tar"

mkdir -p "$DATA_DIR" "$OUT_DIR"

echo "[1/4] Ensuring conda env exists and up-to-date: $ENV_NAME"
if [ ! -f "$ENV_FILE" ]; then
  echo "[error] Could not find $ENV_FILE. Ensure you are running from the pathway-builder-starter folder." >&2
  exit 1
fi
if ! conda env list | grep -q "^$ENV_NAME\s"; then
  conda env create -f "$ENV_FILE"
else
  echo "Updating environment $ENV_NAME from $ENV_FILE (use --prune)"
  conda env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
fi

echo "[2/4] Verifying env with conda run"
# Prefer conda run to avoid activation issues in non-interactive shells
$RUN_TOOL -n "$ENV_NAME" python -c "import sys; print('Python', sys.version)" >/dev/null 2>&1 || {
  echo "[warn] conda run check failed; attempting activation fallback"
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
}

echo "[2.1] Resolving env python path"
# Parse conda env list to find the prefix for $ENV_NAME
ENV_PREFIX="$(conda env list | awk -v n="$ENV_NAME" '$1==n{print $NF}' | tail -n1)"
if [ -z "$ENV_PREFIX" ]; then
  # Best-effort fallback
  ENV_PREFIX="$(conda info --base)/envs/$ENV_NAME"
fi
PYBIN="$ENV_PREFIX/bin/python"
if [ ! -x "$PYBIN" ]; then
  echo "[warn] Could not find $PYBIN; attempting 'conda run -n $ENV_NAME which python'"
  PYBIN="$($RUN_TOOL -n "$ENV_NAME" which python 2>/dev/null || true)"
fi
if [ ! -x "$PYBIN" ]; then
  echo "[error] Unable to resolve python interpreter for env $ENV_NAME" >&2
  exit 2
fi
echo "Using python: $PYBIN"

echo "[2.2] Checking core Python packages in $ENV_NAME"
if ! "$PYBIN" - <<'PY'
import sys, pkgutil
missing=[m for m in ("numpy","scanpy","anndata","pandas","scipy","matplotlib") if pkgutil.find_loader(m) is None]
sys.exit(1 if missing else 0)
PY
then
  echo "[fix] Installing missing core packages into $ENV_NAME via pip (env-local)"
  "$PYBIN" -m pip install --upgrade pip >/dev/null 2>&1 || true
  "$PYBIN" -m pip install numpy pandas scipy scanpy anndata matplotlib
fi

echo "[3/4] Downloading GEO tar if missing"
RAW_TAR="$DATA_DIR/GSE272085_RAW.tar"
if [ ! -f "$RAW_TAR" ]; then
  echo "Fetching $GEO_TAR_URL"
  curl -L -o "$RAW_TAR" "$GEO_TAR_URL"
else
  echo "Found existing $RAW_TAR"
fi

echo "[4/4] Extracting tar to $DATA_DIR"
tar -xf "$RAW_TAR" -C "$DATA_DIR"

echo "[RUN] Unified CLI analysis -> $OUT_DIR"
# Avoid thread oversubscription / hangs on macOS BLAS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# Run with the env's python directly (bypasses conda run issues)
"$PYBIN" run_bdnf_enrichment.py \
  --sn_data "$DATA_DIR/*.h5" \
  --pathway_csv genes/mouse_bdnf_mature_trkb_gene_list_full_annotated_sources.csv \
  --label TRKB \
  --condition_col condition --celltype_col celltype_guess \
  --output_dir "$OUT_DIR"
echo "Done. Outputs in $OUT_DIR"
