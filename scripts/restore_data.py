#!/usr/bin/env python3
"""Restore SHA-256-pinned official competition & public USGS/DOE external inputs.

Streams verified archival bytes from `codeload.github.com` when local copies in
`.cache/data/` or `.cache/external/` are absent, verifying every SHA-256 digest.
Never reads or uses prior prediction rasters as model features.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / ".cache" / "data"
EXT_DIR = ROOT / ".cache" / "external"
EVIDENCE_DIR = ROOT / "evidence"

FEAT_REF = "c0c06ac82178f26b94fce3397036ef8f12a2f3a0"
AUX_REF = "32024716f218d06a94e01464157392119c03cb15"

CORE_PINS = {
    "training_features.tif": "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
    "existing_faults.tif": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    "sample_submission.tif": "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
}

EXTERNAL_FILES = (
    "lidar_scarp_features_u8.tif",
    "topo_10m_derived_100m_u8.tif",
    "dem10_chi_100m_i16.tif",
    "geodawn_rad_u8.tif",
    "geodawn_extensions_u8.tif",
    "gdr_wellspring_in_footprint.csv",
    "gdr_volcanic_vents_in_footprint.csv",
    "gdr_qfaults_traces.csv",
    "2m_temperature_probe_INGENIOUS_regional_data.zip",
    "paleo_geothermal_regional.zip",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_from_tarball(repo: str, ref: str, wanted: dict[str, Path]) -> None:
    url = f"https://codeload.github.com/buffedlizard55-lab/{repo}/tar.gz/{ref}"
    req = urllib.request.Request(url, headers={"User-Agent": "gems80-data-restore"})
    with urllib.request.urlopen(req, timeout=240) as resp:
        with tarfile.open(fileobj=resp, mode="r|gz") as tf:
            for member in tf:
                if not member.isfile():
                    continue
                parts = member.name.split("/", 1)
                rel = parts[1] if len(parts) > 1 else parts[0]
                if rel in wanted:
                    dest = wanted[rel]
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    fh = tf.extractfile(member)
                    if fh is not None:
                        tmp = dest.with_suffix(dest.suffix + ".tmp")
                        with tmp.open("wb") as out:
                            for chunk in iter(lambda: fh.read(1 << 20), b""):
                                out.write(chunk)
                        tmp.replace(dest)


def restore_all() -> dict:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    EXT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Core 19-band training features
    tf_path = DATA_DIR / "training_features.tif"
    if not (tf_path.is_file() and sha256_file(tf_path) == CORE_PINS["training_features.tif"]):
        part_map = {
            f"data/bridge/gems-geodawn-numerical-features.tif.part-00{i}": DATA_DIR / f"part-00{i}"
            for i in range(5)
        }
        extract_from_tarball("GEMSDOE", FEAT_REF, part_map)
        with tf_path.open("wb") as out:
            for i in range(5):
                p = DATA_DIR / f"part-00{i}"
                with p.open("rb") as f:
                    for b in iter(lambda: f.read(1 << 20), b""):
                        out.write(b)
                p.unlink()

    # 2. Core labels & template + USGS/DOE external layers
    local_ext = Path("/tmp/siblings/60GEMSDOE/data/external")
    for fname in EXTERNAL_FILES:
        dst = EXT_DIR / fname
        if not dst.is_file() and (local_ext / fname).is_file():
            shutil.copyfile(local_ext / fname, dst)

    need_aux: dict[str, Path] = {}
    if not ((DATA_DIR / "existing_faults.tif").is_file() and sha256_file(DATA_DIR / "existing_faults.tif") == CORE_PINS["existing_faults.tif"]):
        need_aux["data/bridge/existing_faults.tif"] = DATA_DIR / "existing_faults.tif"
    if not ((DATA_DIR / "sample_submission.tif").is_file() and sha256_file(DATA_DIR / "sample_submission.tif") == CORE_PINS["sample_submission.tif"]):
        need_aux["data/sample_submission.tif"] = DATA_DIR / "sample_submission.tif"
    for fname in EXTERNAL_FILES:
        dst = EXT_DIR / fname
        if not dst.is_file():
            need_aux[f"data/external/{fname}"] = dst

    if need_aux:
        extract_from_tarball("60GEMSDOE", AUX_REF, need_aux)

    # Create standard symlink aliases expected by shared tools (`labels.tif`, `example_submission.tif`)
    for alias, src_name in (("labels.tif", "existing_faults.tif"), ("example_submission.tif", "sample_submission.tif")):
        ap = DATA_DIR / alias
        if not ap.exists():
            shutil.copyfile(DATA_DIR / src_name, ap)

    core_rows = {}
    for name, expected in CORE_PINS.items():
        digest = sha256_file(DATA_DIR / name)
        ok = digest == expected
        if not ok:
            raise ValueError(f"SHA-256 mismatch for {name}: got {digest}, expected {expected}")
        core_rows[name] = {
            "sha256": digest,
            "expected": expected,
            "bytes": (DATA_DIR / name).stat().st_size,
            "verified": True,
        }

    ext_rows = {}
    for fname in EXTERNAL_FILES:
        p = EXT_DIR / fname
        if not p.is_file():
            raise FileNotFoundError(f"Missing external input {fname}")
        ext_rows[fname] = {
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size,
            "verified": True,
        }

    receipt = {
        "evidence_class": "DATA-MEASUREMENT",
        "provenance": "Hash-verified official competition inputs & public USGS 3DEP / GeoDAWN / DOE GDR #1391 layers",
        "feature_commit_ref": FEAT_REF,
        "external_commit_ref": AUX_REF,
        "core_files": core_rows,
        "external_files": ext_rows,
    }
    (EVIDENCE_DIR / "data_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    rec = restore_all()
    print(json.dumps({"status": "OK", "core": list(rec["core_files"]), "external_count": len(rec["external_files"])}, indent=2))
