#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=/Users/dr/Source/rust/umol
CRATE_DIR="$ROOT_DIR/umol-models-graph"
INPUT_DEFAULT="$ROOT_DIR/materials/formats/opensmiles/examples/ZINC.FL.smi"

INPUT_PATH="${1:-$INPUT_DEFAULT}"

echo "Building parser binary..." >&2
cargo build --release -p umol-models-graph --bin smiles_parse >/dev/null

BIN="$CRATE_DIR/target/release/smiles_parse"
hyperfine --warmup 2 --min-runs 5 "${BIN} ${INPUT_PATH}"


