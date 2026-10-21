#!/usr/bin/env python3
"""Fail-closed entry point pending shared generation infrastructure.

The former private baseline silently generated synthetic geology on an invented
submission grid. It is retired. Do not replace the shared submission_writer.py
with another private writer or infer readiness from a matching grid alone.
Exit 2 means blocked; no raster, ZIP, or sidecar is created.
"""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template', type=Path, default=ROOT / 'data/sample_submission.tif')
    parser.add_argument('--features', type=Path, default=ROOT / 'data/training_features.tif')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args(argv)
    blockers = []
    for label, path in [('template', args.template), ('features', args.features),
                        ('shared writer', ROOT / 'scripts/submission_writer.py'),
                        ('shared evaluator', ROOT / 'scripts/evaluate_holdout.py')]:
        if not path.is_file():
            blockers.append(f'{label} absent: {path}')
    blockers.append('Assigned lane, audited cached stack, registry comparisons and promotion evidence required')
    print('BLOCKED: no output written.\n' + '\n'.join(blockers), file=sys.stderr)
    return 2


if __name__ == '__main__':
    sys.exit(main())
