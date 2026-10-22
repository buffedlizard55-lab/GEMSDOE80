#!/usr/bin/env python3
"""Re-read shipped bytes and evaluate protocol; does not train or spend a slot."""
import hashlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gems80.protocol import strict_registry_gate


def audit():
    card = json.loads((ROOT / 'run_card.json').read_text())
    path = ROOT / card['raster_file']
    with rasterio.open(path) as ds:
        a = ds.read(1)
        checks = dict(single_band=ds.count == 1, float32=ds.dtypes == ('float32',),
                      no_nan_anywhere=bool(np.isfinite(a).all()),
                      values_in_0_1=bool(np.isfinite(a).all() and (a >= 0).all() and (a <= 1).all()),
                      crs_matches_mirror_pin=ds.crs.to_epsg() == 32611,
                      shape_matches_mirror_pin=ds.shape == (3730, 3292),
                      transform_matches_mirror_pin=tuple(ds.transform)[:6] == (100, 0, 243350, 0, -100, 4508550),
                      sha256_matches_card=hashlib.sha256(path.read_bytes()).hexdigest() == card['raster_sha256'])
    with zipfile.ZipFile(ROOT / card['zip_file']) as z:
        checks['zip_single_identical_tif'] = z.namelist() == [path.name] and z.read(path.name) == path.read_bytes()
    gate = strict_registry_gate(json.loads((ROOT / 'evidence/h80_uniqueness_audit.json').read_text()))
    return dict(schema='gems80.release-audit.v1', checks=checks,
                format_consistent_with_mirror_pin=all(checks.values()),
                organizer_template_authenticated=False,
                uniqueness=gate, ok_to_download=True, ok_to_submit=False,
                verdict='negative', new_experiments=0, weekly_slots_used=0,
                reason='Recorded overlap exceeds literal protocol; experiment budget exhausted. Not slot-selected.')


if __name__ == '__main__':
    result = audit()
    (ROOT / 'reports/release-audit-20261022.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
