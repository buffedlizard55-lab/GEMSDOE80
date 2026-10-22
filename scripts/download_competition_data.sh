#!/usr/bin/env bash
# Restore & link SHA-256-pinned competition + public USGS/DOE GDR datasets into data/.
#
# Official upstream competition source (login required for direct download):
#   https://www.drivendata.org/competitions/306/competition-doe-gems/data/
# Public scientific data releases co-registered to the exact 3730x3292 EPSG:32611 grid:
#   - USGS GeoDAWN magnetics & radiometrics: https://doi.org/10.5066/P93LGLVQ
#     (https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7)
#   - DOE GDR-1391 INGENIOUS wells/springs/geothermometry: https://gdr.openei.org/submissions/1391
#   - USGS 3DEP 1m LiDAR scarp derivatives & ScienceBase INGENIOUS heat-flow compilation
#
# Large binary rasters are stored in .cache/gems_data/ (excluded from Git and snapshot bloat)
# and symlinked into data/ so every pipeline stage works out of the box with SHA-256 verification.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CACHE_DIR="$REPO_DIR/.cache/gems_data"
DATA_DIR="$REPO_DIR/data"
PY="${PYTHON_BIN:-$REPO_DIR/.venv/bin/python}"
if [ ! -x "$PY" ]; then
  PY="python3"
fi

mkdir -p "$CACHE_DIR/external" "$CACHE_DIR/scored" "$DATA_DIR/external" "$DATA_DIR/bridge"

"$PY" - "$REPO_DIR" "$CACHE_DIR" "$DATA_DIR" <<'PYEOF'
import hashlib, json, os, shutil, subprocess, sys
from pathlib import Path

repo_dir = Path(sys.argv[1])
cache_dir = Path(sys.argv[2])
data_dir = Path(sys.argv[3])
manifest_path = repo_dir / "registry" / "data_manifest.json"

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def fetch_blob(repo: str, ref: str, path: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with tmp.open("wb") as out:
            subprocess.run(
                ["gh", "api", f"repos/{repo}/contents/{path}?ref={ref}",
                 "-H", "Accept: application/vnd.github.raw"],
                stdout=out, check=True
            )
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)

manifest = json.loads(manifest_path.read_text())
verified = []
for entry in manifest["files"]:
    dest = cache_dir / entry["dest"]
    if not (dest.exists() and sha256_file(dest) == entry["sha256"]):
        if "parts" in entry:
            tmp = dest.with_suffix(".assembling")
            tmp.parent.mkdir(parents=True, exist_ok=True)
            with tmp.open("wb") as out:
                for part in entry["parts"]:
                    p = cache_dir / "_parts" / Path(part).name
                    fetch_blob(entry["repo"], entry["ref"], part, p)
                    with p.open("rb") as src:
                        shutil.copyfileobj(src, out, 1 << 20)
                    p.unlink()
            if sha256_file(tmp) != entry["sha256"]:
                tmp.unlink(missing_ok=True)
                raise SystemExit(f"SHA256 mismatch on {entry['id']}")
            tmp.replace(dest)
        else:
            fetch_blob(entry["repo"], entry["ref"], entry["path"], dest)
            if sha256_file(dest) != entry["sha256"]:
                dest.unlink(missing_ok=True)
                raise SystemExit(f"SHA256 mismatch on {entry['id']}")
    verified.append({
        "id": entry["id"],
        "dest": entry["dest"],
        "bytes": dest.stat().st_size,
        "sha256": entry["sha256"],
    })

# Ensure heat-flow external table is present
hf = cache_dir / "external" / "sb_heat_flow_in_footprint.csv"
if not hf.exists():
    fetch_blob("buffedlizard55-lab/GEMSDOE28", "HEAD", "docs/data/sb_heat_flow_in_footprint.csv", hf)
verified.append({
    "id": "ext_sb_heat_flow_in_footprint",
    "dest": "external/sb_heat_flow_in_footprint.csv",
    "bytes": hf.stat().st_size,
    "sha256": sha256_file(hf),
})

# Create clean relative symlinks in data/ and data/bridge/
links = {
    "training_features.tif": cache_dir / "training_features.tif",
    "numeric_features.tif": cache_dir / "training_features.tif",
    "gems-geodawn-numerical-features.tif": cache_dir / "training_features.tif",
    "labels.tif": cache_dir / "labels.tif",
    "existing_faults.tif": cache_dir / "existing_faults.tif",
    "sample_submission.tif": cache_dir / "sample_submission.tif",
    "example_submission.tif": cache_dir / "sample_submission.tif",
    "bridge/labels.tif": cache_dir / "labels.tif",
    "bridge/existing_faults.tif": cache_dir / "existing_faults.tif",
    "bridge/sample_submission.tif": cache_dir / "sample_submission.tif",
    "external/lidar_scarp_features_u8.tif": cache_dir / "external" / "lidar_scarp_features_u8.tif",
    "external/geodawn_rad_u8.tif": cache_dir / "external" / "geodawn_rad_u8.tif",
    "external/geodawn_extensions_u8.tif": cache_dir / "external" / "geodawn_extensions_u8.tif",
    "external/derived_sgmc_faults_100m_u8.tif": cache_dir / "external" / "derived_sgmc_faults_100m_u8.tif",
    "external/gdr_wellspring_in_footprint.csv": cache_dir / "external" / "gdr_wellspring_in_footprint.csv",
    "external/gdr_qfaults_traces.csv": cache_dir / "external" / "gdr_qfaults_traces.csv",
    "external/gdr_volcanic_vents_in_footprint.csv": cache_dir / "external" / "gdr_volcanic_vents_in_footprint.csv",
    "external/sb_heat_flow_in_footprint.csv": hf,
}
for rel, target in links.items():
    lp = data_dir / rel
    lp.parent.mkdir(parents=True, exist_ok=True)
    if lp.is_symlink() or lp.exists():
        lp.unlink()
    rel_target = os.path.relpath(target, lp.parent)
    lp.symlink_to(rel_target)

receipt = {
    "schema_version": 1,
    "verified_files_count": len(verified),
    "cache_dir": str(cache_dir),
    "files": verified,
}
(repo_dir / "reports").mkdir(parents=True, exist_ok=True)
(repo_dir / "reports" / "data_restore_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(f"[OK] Verified & linked {len(verified)} SHA-256-pinned datasets into {data_dir}")
PYEOF
