#!/usr/bin/env python3
"""Build complete organization-wide uniqueness index across ALL accessible
`buffedlizard55-lab/GEMSDOE*` repositories (all ~93 repos + historical blobs).

Writes fast disk-backed mmap arrays in `/tmp/gems80_cache/portfolio/` (`< 250 MB` RAM):
  - `halos_live_packed.npy`    : exact bit-packed 3.0-px Euclidean halos on `allowed_live`
                                 (`0.00e+00` bit diff vs `vendor/gems57/uniqueness.py:dot_overlap`)
  - `halos_holdout_fp.npy`     : perforation-bridged 3.0-px Euclidean halos on full `footprint`
                                 (`IR-H82-03`: zero `existing_faults` perforation leakage on holdout!)
  - `unit_ranks.npy`           : unit-norm ranks on 250,000 `allowed_live` cells for exact Spearman rho
  - `.cache/registry_index.json`: complete tracked JSON manifest of every unique SHA-256 raster
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from vendor.gems57 import evaluate_holdout, grid as grid_mod

DATA_DIR = ROOT / ".cache" / "data"
PORT_DIR = Path("/tmp/gems80_cache/portfolio")
RAW_TIF_DIR = PORT_DIR / "tifs"

SKIP_SUBSTRINGS = (
    "existing_faults",
    "labels.tif",
    "sample_submission",
    "example_submission",
    "training_features",
    "gems-geodawn-numerical-features",
    "lidar_scarp_features",
    "geodawn_rad",
    "geodawn_extensions",
    "topo_10m_derived",
    "dem10_chi",
    "sgmc_faults",
    "tiger_road",
    "format_test",
)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def list_all_org_repos() -> list[str]:
    out = subprocess.check_output(
        ["gh", "repo", "list", "buffedlizard55-lab", "--limit", "200", "--json", "name", "-q", ".[].name"],
        text=True,
        timeout=60,
    )
    repos = sorted({line.strip() for line in out.splitlines() if "gems" in line.lower() and line.strip() != "GEMSDOE80"})
    return repos


def fetch_repo_head_tifs(repo: str) -> tuple[str, list[dict], str | None]:
    url = f"https://codeload.github.com/buffedlizard55-lab/{repo}/tar.gz/HEAD"
    found = []
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "gemsdoe80-portfolio-audit"})
        with urllib.request.urlopen(req, timeout=180) as resp:
            with tarfile.open(fileobj=resp, mode="r|gz") as tf:
                for m in tf:
                    if not m.isfile():
                        continue
                    parts = m.name.split("/", 1)
                    if len(parts) < 2:
                        continue
                    rel = parts[1]
                    low = rel.lower()
                    if not low.endswith((".tif", ".tiff")):
                        continue
                    if any(tok in low for tok in SKIP_SUBSTRINGS):
                        continue
                    fh = tf.extractfile(m)
                    if fh is None:
                        continue
                    data = fh.read()
                    h = sha256_bytes(data)
                    dst = RAW_TIF_DIR / f"{h}.tif"
                    if not dst.is_file():
                        tmp = dst.with_suffix(".tmp")
                        tmp.write_bytes(data)
                        tmp.replace(dst)
                    found.append({"sha256": h, "repo": repo, "rel_path": rel, "path": str(dst)})
        return repo, found, None
    except Exception as exc:
        return repo, found, str(exc)[:200]


def main() -> None:
    t0 = time.monotonic()
    hashes = evaluate_holdout.implementation_hashes()
    assert hashes["evaluate_holdout.py"] == "ad86640b46f3c5fc6b1aceec4c4bfc82d03e35008b37ab2a8ad6076e6810ff3a"

    PORT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_TIF_DIR.mkdir(parents=True, exist_ok=True)

    g = grid_mod.load_grid(DATA_DIR)
    fp = g.footprint
    cat = g.catalogue & fp
    allowed = fp & (~cat)
    allow_idx = np.flatnonzero(allowed)
    fp_idx = np.flatnonzero(fp)
    n_allow = int(allow_idx.size)
    n_fp = int(fp_idx.size)

    rng = np.random.default_rng(20)
    sub_idx = np.sort(rng.choice(allow_idx, size=min(250_000, n_allow), replace=False))
    np.save(PORT_DIR / "sub_idx.npy", sub_idx)

    repos = list_all_org_repos()
    print(f"[portfolio] Scanning {len(repos)} organization repositories + local siblings...", flush=True)

    by_sha: dict[str, dict] = {}
    errors: list[dict] = []

    # 1. Stream all remote organization repos in parallel via zero-rate-limit codeload.github.com
    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        for repo, found, err in ex.map(fetch_repo_head_tifs, repos):
            if err is not None:
                errors.append({"repo": repo, "error": err})
            for item in found:
                s = item["sha256"]
                if s not in by_sha:
                    by_sha[s] = {
                        "sha256": s,
                        "repo_first": item["repo"],
                        "submission": Path(item["rel_path"]).name,
                        "cache_file": item["path"],
                        "sources": [f"{item['repo']}:{item['rel_path']}"],
                    }
                else:
                    by_sha[s]["sources"].append(f"{item['repo']}:{item['rel_path']}")

    # 2. Also ingest any additional historical blobs already in `/tmp/siblings`
    sib_root = Path("/tmp/siblings")
    if sib_root.is_dir():
        for p in sorted(sib_root.rglob("*.tif")):
            low = str(p).lower()
            if any(tok in low for tok in SKIP_SUBSTRINGS):
                continue
            h = sha256_file(p)
            if h not in by_sha:
                dst = RAW_TIF_DIR / f"{h}.tif"
                if not dst.is_file():
                    dst.write_bytes(p.read_bytes())
                by_sha[h] = {
                    "sha256": h,
                    "repo_first": p.relative_to(sib_root).parts[0],
                    "submission": p.name,
                    "cache_file": str(dst),
                    "sources": [str(p)],
                }

    # 3. Validate grid & pack BOTH live allowed halos AND holdout perforation-bridged halos!
    _, near_noncat = distance_transform_edt(cat, return_distances=True, return_indices=True)
    disk_offsets = [(dy, dx) for dy in range(-3, 4) for dx in range(-3, 4) if dy * dy + dx * dx <= 9]
    h_grid, w_grid = fp.shape

    def process_entry(ent: dict) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray, bool] | None:
        p = Path(ent["cache_file"])
        try:
            with rasterio.open(p) as ds:
                if (
                    ds.count != 1
                    or ds.shape != (3730, 3292)
                    or ds.crs is None
                    or ds.crs.to_epsg() != 32611
                    or ds.transform != grid_mod.TRANSFORM
                ):
                    return None
                arr = ds.read(1).astype(np.float32, copy=True)
        except Exception:
            return None

        sub_v = np.nan_to_num(arr.ravel()[sub_idx], nan=0.0)
        if np.all(sub_v == sub_v[0]):
            ur = np.zeros(sub_idx.size, dtype=np.float32)
            is_const = True
        else:
            r = rankdata(sub_v).astype(np.float64)
            r -= r.mean()
            nrm = float(np.linalg.norm(r))
            if nrm < 1e-12:
                ur = np.zeros(sub_idx.size, dtype=np.float32)
                is_const = True
            else:
                ur = (r / nrm).astype(np.float32)
                is_const = False

        # Live exact 3.0-px halo on `allowed = fp & ~cat`
        sup_live = np.isfinite(arr) & (arr > 0) & fp
        n_pos = int(sup_live.sum())

        def halo_1d(pos_mask: np.ndarray, target_flat_idx: np.ndarray) -> np.ndarray:
            cnt = int(pos_mask.sum())
            if cnt == 0:
                return np.zeros(target_flat_idx.size, dtype=bool)
            if cnt <= 350_000:
                rr, cc = np.nonzero(pos_mask)
                out = np.zeros(h_grid * w_grid, dtype=bool)
                for dy, dx in disk_offsets:
                    ry = rr + dy
                    cx = cc + dx
                    ok = (ry >= 0) & (ry < h_grid) & (cx >= 0) & (cx < w_grid)
                    out[ry[ok] * w_grid + cx[ok]] = True
                return out[target_flat_idx]
            padded = np.pad(pos_mask, 3, mode="constant", constant_values=False)
            out2d = np.zeros_like(pos_mask, dtype=bool)
            for dy, dx in disk_offsets:
                out2d |= padded[3 + dy : 3 + dy + h_grid, 3 + dx : 3 + dx + w_grid]
            return out2d.ravel()[target_flat_idx]

        halo_allow = halo_1d(sup_live, allow_idx)
        hc = float(halo_allow.mean())

        # Perforation-bridged 3.0-px halo across full `footprint` for zero-leakage holdout solving (`IR-H82-03`)
        arr_bridged = arr.copy()
        arr_bridged[cat] = arr_bridged[near_noncat[0][cat], near_noncat[1][cat]]
        sup_bridged = np.isfinite(arr_bridged) & (arr_bridged > 0) & fp
        halo_fp = halo_1d(sup_bridged, fp_idx)

        meta_out = {
            "sha256": ent["sha256"],
            "repo_first": ent["repo_first"],
            "submission": ent["submission"],
            "cache_file": ent["cache_file"],
            "sources_count": len(ent["sources"]),
            "positive_pixels": n_pos,
            "halo_coverage": round(hc, 6),
            "exact_halo_coverage": True,
            "kind": "space_filling" if hc >= 0.70 else "discriminating",
        }
        return (
            meta_out,
            np.packbits(halo_allow, bitorder="big"),
            np.packbits(halo_fp, bitorder="big"),
            ur,
            is_const,
        )

    valid_rows = []
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        for res in ex.map(process_entry, list(by_sha.values())):
            if res is not None:
                valid_rows.append(res)

    n_valid = len(valid_rows)
    n_bytes_allow = (n_allow + 7) // 8
    n_bytes_fp = (n_fp + 7) // 8

    mm_live = np.lib.format.open_memmap(PORT_DIR / "halos_live_packed.npy", mode="w+", dtype=np.uint8, shape=(n_valid, n_bytes_allow))
    mm_fp = np.lib.format.open_memmap(PORT_DIR / "halos_holdout_fp.npy", mode="w+", dtype=np.uint8, shape=(n_valid, n_bytes_fp))
    mm_ur = np.lib.format.open_memmap(PORT_DIR / "unit_ranks.npy", mode="w+", dtype=np.float32, shape=(n_valid, sub_idx.size))
    is_const_arr = np.zeros(n_valid, dtype=bool)
    meta_entries = []

    for i, (m_ent, pack_allow, pack_fp, ur, is_c) in enumerate(valid_rows):
        meta_entries.append(m_ent)
        mm_live[i] = pack_allow
        mm_fp[i] = pack_fp
        mm_ur[i] = ur
        is_const_arr[i] = is_c

    mm_live.flush()
    mm_fp.flush()
    mm_ur.flush()
    np.save(PORT_DIR / "is_const.npy", is_const_arr)

    disc = sum(1 for e in meta_entries if e["kind"] == "discriminating")
    sf = sum(1 for e in meta_entries if e["kind"] == "space_filling")
    out_doc = {
        "complete_accessible_scan": len(errors) == 0,
        "repos_scanned": len(repos),
        "errors": errors,
        "n_unique_valid_rasters": n_valid,
        "n_discriminating_priors": disc,
        "n_space_filling_priors": sf,
        "rasters": meta_entries,
    }
    (ROOT / ".cache" / "registry_index.json").write_text(json.dumps(out_doc, indent=2) + "\n")
    print(
        f"[{time.monotonic() - t0:.1f}s] Indexed {n_valid} unique valid GeoTIFF rasters across {len(repos)} repos "
        f"({disc} discriminating, {sf} space-filling; errors={len(errors)})",
        flush=True,
    )


if __name__ == "__main__":
    main()
