"""Literal, fail-closed parallel-run gate; no density exemptions."""
import math


def strict_registry_gate(audit):
    try:
        measurements = {
            'surface_abs_rho': float(audit['pre_placement_surface_gate']['max_abs_spearman_rho']),
            'dots_abs_rho': float(audit['post_placement_binding_like_for_like_gate']['worst_abs_spearman']),
            'all_registry_forward_overlap': float(audit['literal_reading_with_space_filling_priors']['worst_forward_overlap']),
        }
        complete = int(audit['rasters_compared']) > 0 and all(
            math.isfinite(v) and 0 <= v <= 1 for v in measurements.values())
    except (KeyError, TypeError, ValueError, OverflowError):
        return {'passed': False, 'duplicate': None, 'reason': 'missing or invalid registry evidence'}
    duplicate = (measurements['surface_abs_rho'] > .90
                 or measurements['dots_abs_rho'] > .90
                 or measurements['all_registry_forward_overlap'] > .70)
    return dict(measurements, passed=complete and not duplicate,
                duplicate=duplicate if complete else None,
                reason='literal all-raster thresholds; no space-filling exemption')
