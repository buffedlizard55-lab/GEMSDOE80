#!/usr/bin/env bash
# Place competition inputs in data/.
#
# Official source (login required):
#   https://www.drivendata.org/competitions/306/competition-doe-gems/data/
#
# This sandbox cannot log in to DrivenData. A previous session restored
# sha256-pinned owner mirrors from buffedlizard55-lab/GEMSDOE24 and GEMSDOE
# (NOT organizer-authenticated). Pins:
#   sample_submission.tif  2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc
#   labels.tif             7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093
#   training_features.tif  4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$REPO_DIR/data"
cat <<EOF
Put these files in $REPO_DIR/data :
  training_features.tif
  labels.tif
  sample_submission.tif
  1m_DEM_links.csv   (optional)

Then:
  python "$REPO_DIR/scripts/prepare_data.py"
  python "$REPO_DIR/scripts/validate_submission.py" \\
      "$REPO_DIR/docs/downloads/gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif" \\
      --template "$REPO_DIR/data/sample_submission.tif"
EOF
exit 0
