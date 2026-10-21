"""Integration tests for GeoTIFF format validation (requires numpy and rasterio)."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin
except ImportError:
    np = None


@unittest.skipIf(np is None, "install numpy and rasterio to run GeoTIFF integration tests")
class SubmissionValidatorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.transform = from_origin(100, 200, 10, 10)
        self.template = self.root / "template.tif"
        self.write(self.template, np.zeros((2, 3), dtype="float32"), nodata=-9999)
        self.script = Path(__file__).resolve().parents[1] / "scripts" / "validate_submission.py"

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write(self, path: Path, values, *, transform=None, crs="EPSG:32611", nodata=None) -> None:
        with rasterio.open(
            path, "w", driver="GTiff", height=values.shape[0], width=values.shape[1],
            count=1, dtype=values.dtype, crs=crs,
            transform=transform or self.transform, nodata=nodata,
        ) as dataset:
            dataset.write(values, 1)

    def validate(self, name: str, values, *, transform=None, crs="EPSG:32611") -> dict:
        prediction = self.root / f"{name}.tif"
        self.write(prediction, values, transform=transform, crs=crs)
        result = subprocess.run(
            [sys.executable, str(self.script), str(prediction), "--template", str(self.template)],
            check=False, capture_output=True, text=True,
        )
        self.assertIn(result.returncode, (0, 1), result.stderr)
        report = json.loads(result.stdout)
        self.assertFalse(report["download_ready"])
        return report

    def test_valid_raster_is_format_valid_but_not_auto_promoted(self) -> None:
        report = self.validate("valid", np.full((2, 3), 0.4, dtype="float32"))
        self.assertTrue(report["format_valid"])
        self.assertTrue(report["checks"]["crs_matches"])
        self.assertTrue(report["checks"]["transform_matches"])

    def test_out_of_range_and_nan_are_rejected(self) -> None:
        report = self.validate("range", np.full((2, 3), 1.1, dtype="float32"))
        self.assertFalse(report["format_valid"])
        values = np.full((2, 3), 0.4, dtype="float32")
        values[0, 0] = np.nan
        report = self.validate("nan", values)
        self.assertFalse(report["format_valid"])
        self.assertFalse(report["checks"]["no_nan_or_infinite_inside_footprint"])

    def test_crs_transform_and_shape_mismatch_are_rejected(self) -> None:
        report = self.validate("crs", np.full((2, 3), 0.4, dtype="float32"), crs="EPSG:4326")
        self.assertFalse(report["checks"]["crs_matches"])
        report = self.validate(
            "transform", np.full((2, 3), 0.4, dtype="float32"),
            transform=from_origin(101, 200, 10, 10),
        )
        self.assertFalse(report["checks"]["transform_matches"])
        report = self.validate("shape", np.full((2, 2), 0.4, dtype="float32"))
        self.assertFalse(report["checks"]["shape_matches"])


if __name__ == "__main__":
    unittest.main()
