"""The historical candidate TIF must exist, hash-match, and retain its recorded metadata.

The recorded grid is not organizer-authenticated unless compared with the official template.
"""
from __future__ import annotations

import hashlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TIF = ROOT / "docs" / "downloads" / "gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif"
EXPECTED_SHA = "d9c4ce410af5c413195aa3a463c2e4c6c223e26c6c8cf3a25a8a95c4b7a2bcc5"


class ShippedTifTests(unittest.TestCase):
    def test_file_exists_and_hashes(self):
        self.assertTrue(TIF.is_file(), f"missing {TIF}")
        h = hashlib.sha256(TIF.read_bytes()).hexdigest()
        self.assertEqual(h, EXPECTED_SHA)

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
            self.assertEqual(int((a > 0.5).sum()), 28807)


if __name__ == "__main__":
    unittest.main()
