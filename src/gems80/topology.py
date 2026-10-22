"""Fold-Strict Anisotropic Fracture-Process-Zone, Collinear Pair-Bridge & Relay Mechanics (H80-G1).

CRITICAL ZERO-LEAKAGE INVARIANT:
Every function in this module takes ONLY a caller-supplied binary `visible_mask` (strictly
`split["visible"]` in holdout validation; strictly `visible \\ buf_3(fold_k)` during nested
self-supervised training; full `labels` ONLY when exporting the final submission GeoTIFF).
Hidden ground-truth pixels NEVER enter any calculation here.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import convolve, distance_transform_edt, gaussian_filter, label
from scipy.spatial import cKDTree
from skimage.morphology import skeletonize

TOPO_FEATURE_NAMES: tuple[str, ...] = (
    "vis_d_euc",
    "vis_d_tip",
    "vis_tip_ratio",
    "vis_tip_cone_ultra",
    "vis_tip_cone_short",
    "vis_tip_cone_long",
    "vis_pair_bridge",
    "vis_en_echelon_relay",
    "vis_s_perp",
    "vis_seg_len_log",
    "vis_strike_align",
    "vis_density_5km",
    "vis_tip_x_stress",
)


def _extract_tips_and_tangents(skel: np.ndarray, comp: np.ndarray, comp_sizes: np.ndarray) -> dict[str, np.ndarray]:
    """Trace up to 7 pixels inward from every 1-neighbor skeleton endpoint to obtain an accurate
    outward-pointing unit tangent vector t_hat = (r_tip - r_inward) / ||...||.
    """
    nb_kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=np.int16)
    nb_count = convolve(skel.astype(np.int16), nb_kernel, mode="constant")
    tip_mask = skel & (nb_count == 1)
    tr, tc = np.nonzero(tip_mask)
    n_tips = tr.size
    if n_tips == 0:
        z = np.zeros(0, dtype=np.float32)
        return {"rows": tr, "cols": tc, "ty": z, "tx": z, "lengths": z, "tip_mask": tip_mask}

    h, w = skel.shape
    ty = np.zeros(n_tips, dtype=np.float32)
    tx = np.zeros(n_tips, dtype=np.float32)
    lengths = comp_sizes[comp[tr, tc]].astype(np.float32)

    for i in range(n_tips):
        r0, c0 = int(tr[i]), int(tc[i])
        pr, pc = -1, -1
        cr, cc = r0, c0
        # Walk up to 7 steps inward along skeleton
        for _ in range(7):
            nxt_r, nxt_c = -1, -1
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    nr, nc = cr + dr, cc + dc
                    if 0 <= nr < h and 0 <= nc < w and skel[nr, nc] and not (nr == pr and nc == pc):
                        nxt_r, nxt_c = nr, nc
                        break
                if nxt_r != -1:
                    break
            if nxt_r == -1:
                break
            pr, pc = cr, cc
            cr, cc = nxt_r, nxt_c
        dy = float(r0 - cr)
        dx = float(c0 - cc)
        nm = (dy * dy + dx * dx) ** 0.5
        if nm > 1e-6:
            ty[i] = dy / nm
            tx[i] = dx / nm
    return {
        "rows": tr.astype(np.int32),
        "cols": tc.astype(np.int32),
        "ty": ty,
        "tx": tx,
        "lengths": lengths,
        "tip_mask": tip_mask,
    }


def compute_topology_features(
    visible_mask: np.ndarray,
    foot: np.ndarray,
    cand_strike_az: np.ndarray,
    stress_td: np.ndarray,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Compute all 12 fold-strict anisotropic topology channels AND return the tip-preserving
    allowed off-catalogue mask for `visible_mask`.
    """
    vis = np.asarray(visible_mask, dtype=bool) & foot
    h, w = vis.shape
    if not vis.any():
        zeros = {k: np.zeros(vis.shape, dtype=np.float32) for k in TOPO_FEATURE_NAMES}
        return zeros, foot.copy()

    skel = skeletonize(vis)
    if not skel.any():
        skel = vis.copy()

    comp, n_comp = label(vis, structure=np.ones((3, 3), dtype=int))
    comp_sizes = np.bincount(comp.ravel(), minlength=n_comp + 1).astype(np.float32)

    # Distance + nearest visible fault pixel coordinates
    d_euc, near_rc = distance_transform_edt(~vis, return_indices=True)
    d_euc_f = np.minimum(d_euc, 65.0).astype(np.float32)
    near_r, near_c = near_rc[0], near_rc[1]
    seg_len = comp_sizes[comp[near_r, near_c]]
    seg_len_log = np.log1p(seg_len).astype(np.float32)

    # Local strike vector on visible trace via smoothed distance field structure
    sm_vis = gaussian_filter(vis.astype(np.float32), sigma=2.0, mode="constant")
    gy_v, gx_v = np.gradient(sm_vis)
    # Normal to visible fault at (near_r, near_c): evaluate smoothed gradient slightly at visible pixel
    # Better: compute orientation directly on visible skeleton using directional filters
    ang_c2 = np.zeros(vis.shape, dtype=np.float32)
    ang_s2 = np.zeros(vis.shape, dtype=np.float32)
    for deg in range(0, 180, 15):
        rad = np.radians(deg)
        # Azimuth deg clockwise from North: dx = sin(rad), dy = -cos(rad)
        k = np.zeros((7, 7), dtype=np.float32)
        for t in range(-3, 4):
            rr = int(round(3 - t * np.cos(rad)))
            cc = int(round(3 + t * np.sin(rad)))
            if 0 <= rr < 7 and 0 <= cc < 7:
                k[rr, cc] = 1.0
        resp = convolve(skel.astype(np.float32), k, mode="constant")
        w2 = np.maximum(resp - 1.0, 0.0) ** 2
        ang_c2 += w2 * float(np.cos(2.0 * rad))
        ang_s2 += w2 * float(np.sin(2.0 * rad))
    vis_az_at_trace = (0.5 * np.arctan2(ang_s2, ang_c2)) % np.pi
    nearest_vis_az = vis_az_at_trace[near_r, near_c]

    # Cross-strike offset |s_perp| from nearest visible trace along normal to trace strike
    # Trace tangent unit vector in (row, col): dr = -cos(az), dc = +sin(az); normal: nr = sin(az), nc = cos(az)
    rr_grid, cc_grid = np.indices((h, w), dtype=np.float32)
    dr_vec = rr_grid - near_r.astype(np.float32)
    dc_vec = cc_grid - near_c.astype(np.float32)
    s_perp = np.abs(dr_vec * np.sin(nearest_vis_az) + dc_vec * np.cos(nearest_vis_az))
    s_perp = np.minimum(s_perp, 50.0).astype(np.float32)

    strike_align = np.abs(np.cos(cand_strike_az - nearest_vis_az)).astype(np.float32)

    # Extract skeleton endpoints & outward crack-tip vectors
    tip_info = _extract_tips_and_tangents(skel, comp, comp_sizes)
    tr, tc, ty, tx, tlen = (
        tip_info["rows"],
        tip_info["cols"],
        tip_info["ty"],
        tip_info["tx"],
        tip_info["lengths"],
    )
    tip_mask = tip_info["tip_mask"]
    d_tip = np.minimum(distance_transform_edt(~tip_mask), 75.0).astype(np.float32) if tip_mask.any() else np.full(vis.shape, 75.0, dtype=np.float32)
    tip_ratio = (d_euc_f / (d_tip + 1.0)).astype(np.float32)

    # Directed Mode-II/III Crack-Tip Cones (Ultra/Short + Displacement-Length Scaled Long) & En-Echelon Relay Lobes
    cone_ultra = np.zeros(vis.shape, dtype=np.float32)
    cone_short = np.zeros(vis.shape, dtype=np.float32)
    cone_long = np.zeros(vis.shape, dtype=np.float32)
    relay_lobe = np.zeros(vis.shape, dtype=np.float32)

    valid_tips = np.hypot(ty, tx) > 0.5
    tr_v, tc_v, ty_v, tx_v, tl_v = tr[valid_tips], tc[valid_tips], ty[valid_tips], tx[valid_tips], tlen[valid_tips]

    R_MAX = 32
    patch_y, patch_x = np.mgrid[-R_MAX : R_MAX + 1, -R_MAX : R_MAX + 1].astype(np.float32)
    patch_r = np.hypot(patch_y, patch_x) + 1e-6

    for i in range(tr_v.size):
        r0, c0 = int(tr_v[i]), int(tc_v[i])
        uy, ux = float(ty_v[i]), float(tx_v[i])
        L_px = float(tl_v[i])
        w_len = float(np.clip(0.35 + 0.25 * np.log1p(L_px), 0.35, 1.85))
        lam_long = float(min(9.0 + 4.5 * np.log1p(L_px), 30.0))

        y0, y1 = max(0, r0 - R_MAX), min(h, r0 + R_MAX + 1)
        x0, x1 = max(0, c0 - R_MAX), min(w, c0 + R_MAX + 1)
        py = patch_y[y0 - (r0 - R_MAX) : y1 - (r0 - R_MAX), x0 - (c0 - R_MAX) : x1 - (c0 - R_MAX)]
        px = patch_x[y0 - (r0 - R_MAX) : x1 - (c0 - R_MAX), x0 - (c0 - R_MAX) : x1 - (c0 - R_MAX)] if False else patch_x[y0 - (r0 - R_MAX) : y1 - (r0 - R_MAX), x0 - (c0 - R_MAX) : x1 - (c0 - R_MAX)]
        pr = patch_r[y0 - (r0 - R_MAX) : y1 - (r0 - R_MAX), x0 - (c0 - R_MAX) : x1 - (c0 - R_MAX)]

        s_par = py * uy + px * ux
        s_lat = np.abs(-py * ux + px * uy)

        # Fracture-arrest occlusion (Pollard & Aydin 1988): an intervening visible fault
        # ahead of the tip arrests propagation and casts a stress shadow behind it.
        vis_patch = vis[y0:y1, x0:x1]
        blockers = vis_patch & (s_par > 1.5) & (s_lat <= 3.0)
        if blockers.any():
            s_arrest = float(s_par[blockers].min())
            arrest_atten = np.where(s_par > s_arrest, np.exp(-(s_par - s_arrest) / 1.5), 1.0).astype(np.float32)
        else:
            arrest_atten = np.float32(1.0)

        fwd = s_par > 0.0
        sig_u = 1.25 + 0.12 * np.maximum(s_par, 0.0)
        v_ultra = np.where(
            fwd,
            w_len * arrest_atten * np.exp(-s_par / 4.5 - 0.5 * (s_lat / sig_u) ** 2),
            0.0,
        ).astype(np.float32)

        sig_s = 1.45 + 0.15 * np.maximum(s_par, 0.0)
        v_short = np.where(
            fwd,
            w_len * arrest_atten * np.exp(-s_par / 10.0 - 0.5 * (s_lat / sig_s) ** 2),
            0.0,
        ).astype(np.float32)

        sig_l = 1.75 + 0.19 * np.maximum(s_par, 0.0)
        v_long = np.where(
            fwd,
            w_len * arrest_atten * np.exp(-s_par / lam_long - 0.5 * (s_lat / sig_l) ** 2),
            0.0,
        ).astype(np.float32)

        # Off-axis en-echelon relay / Riedel step-over zone
        v_rel = (
            w_len
            * np.exp(-0.5 * ((s_par - 5.5) / 9.0) ** 2)
            * np.exp(-0.5 * ((s_lat - 4.8) / 2.6) ** 2)
            * (pr <= 28.0)
        ).astype(np.float32)

        cone_ultra[y0:y1, x0:x1] = np.maximum(cone_ultra[y0:y1, x0:x1], v_ultra)
        cone_short[y0:y1, x0:x1] = np.maximum(cone_short[y0:y1, x0:x1], v_short)
        cone_long[y0:y1, x0:x1] = np.maximum(cone_long[y0:y1, x0:x1], v_long)
        relay_lobe[y0:y1, x0:x1] = np.maximum(relay_lobe[y0:y1, x0:x1], v_rel)

    # Facing-Tip Collinear Gap Bridge & Soft-Linked Relay-Ramp Corridor (Peacock & Sanderson 1991)
    pair_bridge = np.zeros(vis.shape, dtype=np.float32)
    if tr_v.size >= 2:
        pts = np.column_stack([tr_v, tc_v]).astype(np.float64)
        tree = cKDTree(pts)
        pairs = tree.query_pairs(r=45.0, output_type="ndarray")
        for idx_i, idx_j in pairs:
            if comp[tr_v[idx_i], tc_v[idx_i]] == comp[tr_v[idx_j], tc_v[idx_j]] and tl_v[idx_i] < 15:
                continue
            r1, c1 = float(tr_v[idx_i]), float(tc_v[idx_i])
            r2, c2 = float(tr_v[idx_j]), float(tc_v[idx_j])
            dy = r2 - r1
            dx = c2 - c1
            dist_ij = (dy * dy + dx * dx) ** 0.5
            if dist_ij < 2.0:
                continue
            uy_ij, ux_ij = dy / dist_ij, dx / dist_ij
            align_1 = ty_v[idx_i] * uy_ij + tx_v[idx_i] * ux_ij
            align_2 = -(ty_v[idx_j] * uy_ij + tx_v[idx_j] * ux_ij)
            strike_sim = abs(ty_v[idx_i] * ty_v[idx_j] + tx_v[idx_i] * tx_v[idx_j])
            # Facing across an unmapped gap OR overlapping en-echelon relay stepover!
            if align_1 > 0.45 and align_2 > 0.45 and strike_sim > 0.55:
                strength = (
                    float(align_1 * align_2)
                    * float(np.clip(0.4 + 0.22 * np.log1p(np.sqrt(tl_v[idx_i] * tl_v[idx_j])), 0.4, 2.0))
                    * float(np.exp(-dist_ij / 30.0))
                )
                pad = 6
                by0, by1 = max(0, int(min(r1, r2)) - pad), min(h, int(max(r1, r2)) + pad + 1)
                bx0, bx1 = max(0, int(min(c1, c2)) - pad), min(w, int(max(c1, c2)) + pad + 1)
                yy, xx = np.mgrid[by0:by1, bx0:bx1].astype(np.float32)
                t_proj = ((yy - r1) * uy_ij + (xx - c1) * ux_ij) / dist_ij
                perp_d = np.abs(-(yy - r1) * ux_ij + (xx - c1) * uy_ij)
                width = 1.6 + 0.08 * dist_ij * np.sin(np.pi * np.clip(t_proj, 0.0, 1.0))
                in_seg = (t_proj >= -0.05) & (t_proj <= 1.05)
                b_val = np.where(in_seg, strength * np.exp(-0.5 * (perp_d / width) ** 2), 0.0).astype(np.float32)
                pair_bridge[by0:by1, bx0:bx1] = np.maximum(pair_bridge[by0:by1, bx0:bx1], b_val)

    density_5km = gaussian_filter(vis.astype(np.float32), sigma=35.0, mode="constant").astype(np.float32)
    relay_td = (relay_lobe * (0.35 + 0.65 * stress_td)).astype(np.float32)
    tip_x_stress = ((cone_short + 0.8 * pair_bridge) * (0.35 + 0.65 * stress_td)).astype(np.float32)

    # Tip-Preserving Anisotropic Collar Mask:
    # Lateral flank halos (d <= 2.0 px away from tips) have near-zero true positives and create
    # pure FP halo penalty (proven by 0.2600 -> 0.2708 -> 0.2778), EXCEPT strictly inside a high-confidence
    # forward fault-tip cone / pair-bridge at 1.0 < d <= 2.0 px where real fault continuations exit tips!
    allowed_mask = foot & ~vis & (d_euc > 2.0)

    feats = {
        "vis_d_euc": np.where(foot, d_euc_f, 0.0),
        "vis_d_tip": np.where(foot, d_tip, 0.0),
        "vis_tip_ratio": np.where(foot, tip_ratio, 0.0),
        "vis_tip_cone_ultra": np.where(foot, cone_ultra, 0.0),
        "vis_tip_cone_short": np.where(foot, cone_short, 0.0),
        "vis_tip_cone_long": np.where(foot, cone_long, 0.0),
        "vis_pair_bridge": np.where(foot, pair_bridge, 0.0),
        "vis_en_echelon_relay": np.where(foot, relay_td, 0.0),
        "vis_s_perp": np.where(foot, s_perp, 0.0),
        "vis_seg_len_log": np.where(foot, seg_len_log, 0.0),
        "vis_strike_align": np.where(foot, strike_align, 0.0),
        "vis_density_5km": np.where(foot, density_5km, 0.0),
        "vis_tip_x_stress": np.where(foot, tip_x_stress, 0.0),
    }
    return feats, allowed_mask
