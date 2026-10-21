#!/usr/bin/env bash
# Download competition inputs into data/.
#
# NOTICE: The DrivenData data page requires a login. This sandbox has no
# DrivenData credentials, so the files CANNOT be downloaded automatically
# from this machine. Run this on any unrestricted, logged-in machine:
#
#   1. Log in to https://www.drivendata.org/competitions/306/competition-doe-gems/data/
#   2. Download into data/:
#        training_features.tif
#        labels.tif
#        sample_submission.tif
#        1m_DEM_links.csv
#   3. Then run:
#        python scripts/prepare_data.py
#        python scripts/generate_submission.py --template data/sample_submission.tif --features data/training_features.tif
#        python scripts/validate_submission.py docs/submissions/<name>.tif --template data/sample_submission.tif
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$REPO_DIR/data"
mkdir -p "$DATA_DIR"

cat <<EOF

Official competition data page (login required):
  https://www.drivendata.org/competitions/306/competition-doe-gems/data/

Files to download into $DATA_DIR:
  training_features.tif   (multi-band feature raster, EPSG:32611, 100 m)
  labels.tif              (fault label raster)
  sample_submission.tif   (official submission template - match CRS/shape/geotransform exactly)
  1m_DEM_links.csv        (links to 1 m DEM downloads)

After placing the files:
  python "$REPO_DIR/scripts/prepare_data.py"
  python "$REPO_DIR/scripts/generate_submission.py" \\
      --template "$DATA_DIR/sample_submission.tif" \\
      --features "$DATA_DIR/training_features.tif"
  python "$REPO_DIR/scripts/validate_submission.py" "$REPO_DIR/docs/submissions/<name>.tif" --template "$DATA_DIR/sample_submission.tif"

EOF
exit 0
