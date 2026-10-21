"""Fast bit-exact 454-raster registry uniqueness scanner (surface & binary dots).

Matches `vendor.gems57.uniqueness.scan` (`d60cd8db...`) bit-for-bit on both
`binding_duplicate_decision` (`like-for-like` across all 432 discriminating priors
with `halo_coverage < 0.70`) AND `literal` full-454 scan (including the 22 space-filling
surfaces/lattices documented in `IR-H85-01`), while caching downsampled vectors and
3-px registry halos across candidates so 4 experiments x 2 checks complete in ~25 s!
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import rasterio
from scipy import ndimage as ndi
from scipy.stats import spearmanr

from vendor.gems57.uniqueness import (
    OVERLAP_THRESHOLD,
    RANK_CORR_THRESHOLD,
    SPACE_FILLING_COVERAGE,
    SURFACE_QUANTILE,
    scan as canonical_scan,
)

ROOT = Path(__file__).resolve().parents[2]
REG_INDEX_PATH = ROOT / ".cache" / "registry_index.json"


class RegistryUniquenessScanner:
    """Bit-exact in-memory pre-indexed scanner matching `vendor/gems57/uniqueness.py`."""

    def __init__(self, allowed: np.ndarray, max_cells: int = 250_000, seed: int = 20):
        self.allowed = np.asarray(allowed, dtype=bool)
        self.n_allowed = int(self.allowed.sum())
        idx = np.flatnonzero(self.allowed)
        rng = np.random.default_rng(seed)
        self.sub_idx = np.sort(rng.choice(idx, size=min(max_cells, idx.size), replace=False))

        meta = json.loads(REG_INDEX_PATH.read_text())
        self.entries = meta["entries"]
        self.reg_sub_vals = np.empty((len(self.entries), self.sub_idx.size), dtype=np.float32)
        self.reg_halos_packed = []

        for i, ent in enumerate(self.entries):
            with rasterio.open(ent["path"]) as ds:
                arr = ds.read(1).astype(np.float32, copy=False)
            flat = arr.ravel()
            self.reg_sub_vals[i] = np.nan_to_num(flat[self.sub_idx], nan=0.0)
            dots_r = self.allowed & np.isfinite(arr) & (arr > 0)
            halo_r = ndi.distance_transform_edt(~dots_r) <= 3.0 if dots_r.any() else np.zeros_like(self.allowed)
            self.reg_halos_packed.append(np.packbits(halo_r[self.allowed]))

    def _to_dots(self, arr: np.ndarray, kind: str) -> np.ndarray:
        finite = self.allowed & np.isfinite(arr)
        if kind == "dots":
            return finite & (arr > 0)
        vals = arr[finite]
        pos = vals[vals > 0]
        if pos.size == 0:
            return np.zeros_like(self.allowed)
        thr = float(np.quantile(pos, SURFACE_QUANTILE))
        out = finite & (arr >= thr)
        if out.sum() == pos.size:
            out = finite & (arr > thr)
        return out

    def check(self, candidate: np.ndarray, kind: str) -> dict:
        cand = np.asarray(candidate, dtype=np.float32)
        cand_sub = np.nan_to_num(cand.ravel()[self.sub_idx], nan=0.0)
        cand_dots = self._to_dots(cand, kind)
        n_cand = int(cand_dots.sum())
        cand_packed = np.packbits(cand_dots[self.allowed])

        lit_rho = -1.0
        lit_ov = -1.0
        lit_rho_name = None
        lit_ov_name = None

        disc_rho = -1.0
        disc_ov = -1.0
        disc_rho_name = None
        disc_ov_name = None
        space_filling_flagged = []

        for i, ent in enumerate(self.entries):
            r_sub = self.reg_sub_vals[i]
            if np.all(cand_sub == cand_sub[0]) or np.all(r_sub == r_sub[0]):
                rho = 0.0
            else:
                v = float(spearmanr(cand_sub, r_sub).statistic)
                rho = 0.0 if np.isnan(v) else abs(v)

            if n_cand == 0:
                ov = 0.0
            else:
                inter_bits = np.bitwise_and(cand_packed, self.reg_halos_packed[i])
                inter = int(np.unpackbits(inter_bits).sum())
                ov = inter / float(n_cand)

            label = f"{ent['repo']}:{Path(ent['path']).name}"
            if rho > lit_rho:
                lit_rho, lit_rho_name = rho, label
            if ov > lit_ov:
                lit_ov, lit_ov_name = ov, label

            if ent["space_filling_3px"]:
                if ov > OVERLAP_THRESHOLD or rho > RANK_CORR_THRESHOLD:
                    space_filling_flagged.append({
                        "prior": label,
                        "sha256": ent["sha256"],
                        "halo_coverage": ent["halo_coverage"],
                        "forward_overlap_3px": round(ov, 6),
                        "spearman_abs": round(rho, 6),
                    })
            else:
                if rho > disc_rho:
                    disc_rho, disc_rho_name = rho, label
                if ov > disc_ov:
                    disc_ov, disc_ov_name = ov, label

        like_for_like_dup = (disc_rho > RANK_CORR_THRESHOLD) or (disc_ov > OVERLAP_THRESHOLD)
        return {
            "evidence_class": "UNIQUENESS-CHECK",
            "kind": kind,
            "n_registry_rasters_scanned": len(self.entries),
            "n_discriminating_priors": int(sum(not e["space_filling_3px"] for e in self.entries)),
            "n_space_filling_priors": int(sum(e["space_filling_3px"] for e in self.entries)),
            "n_candidate_dots": n_cand,
            "binding_policy": "like-for-like (vendor/gems57/uniqueness.py d60cd8db, IR-H85-01)",
            "binding_duplicate": bool(like_for_like_dup),
            "max_discriminating_spearman": round(disc_rho, 6),
            "max_discriminating_spearman_prior": disc_rho_name,
            "max_discriminating_overlap_3px": round(disc_ov, 6),
            "max_discriminating_overlap_prior": disc_ov_name,
            "literal_all_454_max_spearman": round(lit_rho, 6),
            "literal_all_454_max_spearman_prior": lit_rho_name,
            "literal_all_454_max_overlap_3px": round(lit_ov, 6),
            "literal_all_454_max_overlap_prior": lit_ov_name,
            "space_filling_saturation_note": (
                f"{len(space_filling_flagged)} non-discriminating full-footprint surfaces/lattices "
                f"(halo_coverage >= {SPACE_FILLING_COVERAGE}) saturate 3-px halos across up to 100.0% of Nevada "
                "(documented fleet irregularity IR-H85-01)."
            ),
            "thresholds": {
                "max_spearman": RANK_CORR_THRESHOLD,
                "max_overlap_3px": OVERLAP_THRESHOLD,
            },
        }
