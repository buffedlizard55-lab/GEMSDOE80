#!/usr/bin/env python3
"""Restore integrity-pinned owner mirrors, NOT organizer-authenticated data.
Uses gh authentication without reading or storing credentials. Default: small template only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def restore(entry, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and hashlib.sha256(destination.read_bytes()).hexdigest() == entry['sha256']:
        return
    with tempfile.TemporaryFile() as assembled:
        for part in entry.get('parts', [entry.get('path')]):
            endpoint = f"repos/{entry['repo']}/contents/{part}?ref={entry['ref']}"
            subprocess.run(['gh', 'api', '-H', 'Accept: application/vnd.github.raw+json', endpoint],
                           stdout=assembled, check=True)
        assembled.seek(0)
        digest = hashlib.file_digest(assembled, 'sha256').hexdigest()
        size = assembled.tell()
        if digest != entry['sha256'] or size != entry['bytes']:
            raise ValueError(f"Mirror integrity failure: {entry['id']}")
        assembled.seek(0)
        # Atomic install; never replace a good target with a partial download.
        with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as output:
            staging = Path(output.name)
            try:
                import shutil
                shutil.copyfileobj(assembled, output)
            except BaseException:
                staging.unlink(missing_ok=True)
                raise
        staging.replace(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core', action='store_true', help='Also restore features and labels (~420 MB)')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'research/evidence/gems32-data-manifest.json').read_text())
    for entry in manifest['files']:
        if entry['id'] == 'sample_submission' or (args.core and entry['group'] == 'core'):
            restore(entry, ROOT / 'data' / entry['dest'])
            print(f"Integrity verified: {entry['id']} ({entry['sha256']}); organizer authenticity UNVERIFIED")


if __name__ == '__main__':
    main()
