"""End-to-end integrity checks on packaged GeoTIFFs, ZIP bundle, Run Cards, and GitHub Pages."""
from __future__ import annotations

import hashlib
import json
import unittest
import zipfile
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PRIMARY_TIF = ROOT / "docs" / "downloads" / "gemsdoe80-h80g-choquet-tworegime-20261022T000000Z-zeros.tif"
PRIMARY_ZIP = ROOT / "docs" / "downloads" / "gemsdoe80-h80g-choquet-tworegime-20261022T000000Z-zeros.zip"
PRIMARY_NAN = ROOT / "docs" / "downloads" / "gemsdoe80-h80g-choquet-tworegime-20261022T000000Z-nan.tif"
PRIMARY_SHA = "a3c14d76d1f58802c8bd3f13378cbeb037ad3618d5fbd9dbf5f5a0dae4a35ede"
PRIMARY_ZIP_SHA = "9f7bc446d77ff88ca6d3bd874e0fd82413975d0aea148a8026993bc30b4b81c5"
PRIMARY_DOTS = 51667

LEGACY_TIF = ROOT / "docs" / "downloads" / "gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif"
LEGACY_SHA = "d9c4ce410af5c413195aa3a463c2e4c6c223e26c6c8cf3a25a8a95c4b7a2bcc5"


class _LocalLinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.refs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        d = dict(attrs)
        for key in ("href", "src"):
            v = d.get(key)
            if v and not v.startswith(("http://", "https://", "mailto:", "#")):
                self.refs.append(v.split("#")[0])


class ShippedTifTests(unittest.TestCase):
    def test_primary_h80g_files_exist_and_hash_match(self):
        self.assertTrue(PRIMARY_TIF.is_file(), f"missing {PRIMARY_TIF}")
        self.assertTrue(PRIMARY_ZIP.is_file(), f"missing {PRIMARY_ZIP}")
        self.assertTrue(PRIMARY_NAN.is_file(), f"missing {PRIMARY_NAN}")
        self.assertEqual(hashlib.sha256(PRIMARY_TIF.read_bytes()).hexdigest(), PRIMARY_SHA)
        self.assertEqual(hashlib.sha256(PRIMARY_ZIP.read_bytes()).hexdigest(), PRIMARY_ZIP_SHA)
        with zipfile.ZipFile(PRIMARY_ZIP, "r") as zf:
            names = zf.namelist()
            self.assertEqual(names, [PRIMARY_TIF.name])
            self.assertEqual(hashlib.sha256(zf.read(PRIMARY_TIF.name)).hexdigest(), PRIMARY_SHA)

    def test_primary_h80g_grid_range_and_zero_nans(self):
        import numpy as np
        import rasterio

        with rasterio.open(PRIMARY_TIF) as src:
            self.assertEqual(src.count, 1)
            self.assertEqual(src.dtypes[0], "float32")
            self.assertIsNone(src.nodata)
            self.assertEqual(src.height, 3730)
            self.assertEqual(src.width, 3292)
            self.assertIn("32611", str(src.crs))
            self.assertEqual(
                (src.transform.a, src.transform.b, src.transform.c, src.transform.d, src.transform.e, src.transform.f),
                (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0),
            )
            a = src.read(1)
            self.assertTrue(np.isfinite(a).all())
            self.assertEqual(int((~np.isfinite(a)).sum()), 0)
            self.assertGreaterEqual(float(a.min()), 0.0)
            self.assertLessEqual(float(a.max()), 1.0)
            self.assertEqual(int((a > 0.5).sum()), PRIMARY_DOTS)

    def test_strict_json_run_cards_and_html_links(self):
        for p in (
            ROOT / "research" / "run-card.json",
            ROOT / "docs" / "run-card.json",
            ROOT / "research" / "validator-output.json",
            ROOT / "research" / "h80_experiments_summary.json",
        ):
            raw = p.read_text()
            obj = json.loads(raw, parse_constant=lambda c: (_ for _ in ()).throw(ValueError(f"Non-standard JSON token {c} in {p}")))
            self.assertIsInstance(obj, dict)

        rc = json.loads((ROOT / "docs" / "run-card.json").read_text())
        self.assertEqual(rc["submission_artifacts"]["primary_zeros_sha256"], PRIMARY_SHA)
        self.assertLessEqual(len(rc["submission_note"]), 140)
        self.assertEqual(rc["verdict"], "promote")
        self.assertFalse(rc["leakage_canary"]["leakage_flag"])
        self.assertFalse(rc["registry_uniqueness"]["duplicate_lane"])

        for html_name in ("index.html", "executive-summary.html"):
            page = ROOT / "docs" / html_name
            parser = _LocalLinkParser()
            parser.feed(page.read_text())
            self.assertGreater(len(parser.refs), 0)
            for rel in parser.refs:
                if not rel:
                    continue
                target = (page.parent / rel).resolve()
                self.assertTrue(target.exists(), f"Broken local link {rel!r} in {html_name} -> {target}")

    def test_file_exists_and_hashes(self):
        self.assertTrue(LEGACY_TIF.is_file(), f"missing {LEGACY_TIF}")
        h = hashlib.sha256(LEGACY_TIF.read_bytes()).hexdigest()
        self.assertEqual(h, LEGACY_SHA)

    def test_grid_when_rasterio_available(self):
        import numpy as np
        import rasterio

        with rasterio.open(LEGACY_TIF) as src:
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
