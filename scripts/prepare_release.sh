#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"

cd "$ROOT"
rm -rf dist
python3 -m build --outdir "$ROOT/dist"

echo "Release artifacts written to $ROOT/dist"
