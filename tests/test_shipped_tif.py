"""The packaged zeros TIF (primary H81-A) must exist, hash-match, and sit on the official grid."""
from __future__ import annotations

import hashlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TIF = ROOT / "docs" / "downloads" / "gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif"
EXPECTED_SHA = "701ed7fed8611861f2fe2ec6437c850aede0a6a792b5fc3dbbca10191953bbbd"
ARCHIVED = ROOT / "docs" / "downloads" / "gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif"


class ShippedTifTests(unittest.TestCase):
    def test_primary_file_exists_and_hashes(self):
        self.assertTrue(TIF.is_file(), f"missing primary {TIF}")
        h = hashlib.sha256(TIF.read_bytes()).hexdigest()
        self.assertEqual(h, EXPECTED_SHA)

    def test_archived_still_present(self):
        self.assertTrue(ARCHIVED.is_file(), f"missing archived {ARCHIVED}")

    def test_grid_when_rasterio_available(self):
        try:
            import numpy as np
            import rasterio
        except ImportError:
            self.skipTest("rasterio not installed")
        with rasterio.open(TIF) as src:
            self.assertEqual(src.count, 1)
            self.assertEqual(src.dtypes[0], "float32")
            self.assertEqual(src.height, 3730)
            self.assertEqual(src.width, 3292)
            self.assertIn("32611", str(src.crs))
            a = src.read(1)
            self.assertTrue(np.isfinite(a).all())
            self.assertGreaterEqual(float(a.min()), 0.0)
            self.assertLessEqual(float(a.max()), 1.0)
            self.assertEqual(int((a > 0.5).sum()), 7415)
            # must be off-catalogue
            import rasterio as rio2
            with rio2.open(ROOT / "data" / "labels.tif") as ls:
                lab = ls.read(1)
            # footprint from template
            with rio2.open(ROOT / "data" / "sample_submission.tif") as tmpl:
                foot = np.isfinite(tmpl.read(1))
            cats = (lab >= 1) & foot
            self.assertEqual(int(((a > 0.5) & cats).sum()), 0)


if __name__ == "__main__":
    unittest.main()
