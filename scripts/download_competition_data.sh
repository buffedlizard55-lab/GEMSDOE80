#!/usr/bin/env bash
set -euo pipefail
# Autonomous zero-manual-input competition, external federal geoscience, and
# reference scored raster restoration for GEMSDOE80 (Lane H80-DPVC).
# Every file is SHA-256 verified against registry/data_manifest.json.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python3 scripts/restore_data.py
