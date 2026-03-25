#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="$SCRIPT_DIR/out"
OUT_PDF="$OUT_DIR/NEB_interpo_Manual_zh.pdf"
OUT_TEX="$OUT_DIR/tmp.tex"

mkdir -p "$OUT_DIR"

pandoc "$SCRIPT_DIR/NEB_interpo_Manual_zh.md" \
  --defaults "$SCRIPT_DIR/pandoc.yaml" \
  -t latex \
  --output "$OUT_TEX"

python "$SCRIPT_DIR/patch_tables.py"

cd "$OUT_DIR"
xelatex -interaction=nonstopmode -halt-on-error tmp.tex >/tmp/mgt_xelatex_1.log 2>&1
xelatex -interaction=nonstopmode -halt-on-error tmp.tex >/tmp/mgt_xelatex_2.log 2>&1
mv -f tmp.pdf "$OUT_PDF"

echo "[OK] wrote $OUT_PDF"
