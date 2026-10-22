#!/usr/bin/env python3
"""Reproduce the existing H80-C holdout ONLY; no tuning, new arm, or submission.

Uses pinned cached features and the shared evaluator. This is a random-segment
proxy audit, not the spatially-blocked validation required for promotion.
"""
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt

from run_h80_experiments import (
    ROOT, FEAT, TEMPLATE, LABELS, build_score, load_band, load_labels,
    load_template, nms_dots, withhold_segments, score_holdout, leakage_canary,
)


def main():
    tpl = load_template(TEMPLATE)
    foot = tpl['foot']
    labels = load_labels(LABELS, foot)
    split = withhold_segments(labels, foot)
    arrays = {}
    valid = foot.copy()
    bands = {}
    for name, index in [('dil', 8), ('inv', 4), ('hg', 3), ('tc', 6), ('cond', 17), ('rtp', 2)]:
        arr, mask, desc, _ = load_band(FEAT, index)
        arrays[name] = arr
        valid &= mask
        bands[name] = (arr, mask)
    from gems80.emission import along_fabric_derivative
    bands['along_fabric_derivative'] = (
        along_fabric_derivative(np.where(valid, arrays['dil'], 0), np.where(valid, arrays['rtp'], 0)), valid)
    dist = distance_transform_edt(~split['visible'])
    emit = foot & valid & ~binary_dilation(split['visible'], iterations=4)
    surface = build_score('C', **arrays, valid=valid, dist=dist, emit_valid=emit)
    bands['catalogue_distance_visible_only'] = (dist, foot)
    bands['combined_surface'] = (surface, emit)
    canary = leakage_canary(bands, split['hidden'], foot, split['visible'])
    dots = nms_dots(surface, emit, sep=3, max_k=34033)
    result = score_holdout(dots, split, foot)
    result.update({
        'audit_kind': 'fixed H80-C reproduction; no new experiment or submission file',
        'leakage_canary': canary,
        'promotion_eligible': False,
        'reason': 'Random segment split, not spatially-blocked; best matched baseline absent; registry not rechecked.',
    })
    path = ROOT / 'reports' / 'h80c-reproduction-20261022.json'
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(path)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
