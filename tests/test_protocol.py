import copy
import json
from pathlib import Path

import pytest
from gems80.protocol import strict_registry_gate

ROOT = Path(__file__).resolve().parents[1]


def sample():
    return json.loads((ROOT / 'evidence/h80_uniqueness_audit.json').read_text())


def test_dense_registry_is_not_exempt():
    result = strict_registry_gate(sample())
    assert result['duplicate'] is True
    assert not result['passed']


@pytest.mark.parametrize('bad', [None, float('nan'), float('inf'), -1, 1.1])
def test_invalid_measurement_fails_closed(bad):
    a = copy.deepcopy(sample())
    a['literal_reading_with_space_filling_priors']['worst_forward_overlap'] = bad
    assert not strict_registry_gate(a)['passed']


def test_exact_threshold_passes_but_above_fails():
    a = sample()
    a['literal_reading_with_space_filling_priors']['worst_forward_overlap'] = .70
    assert strict_registry_gate(a)['passed']
    a['literal_reading_with_space_filling_priors']['worst_forward_overlap'] = .70001
    assert not strict_registry_gate(a)['passed']


def test_absent_evidence_fails_closed():
    assert not strict_registry_gate({})['passed']
