"""Gaussian-smoothing sigma resolution for the spatial map.

Smoothing strength is a live slider (v0.31+); the per-dataset value in
datasets_config.yml is just its default. Kept here -- pure, no Dash/data deps -- so
the backward-compat resolution is unit-testable without importing the Dash app.
"""

import math

DEFAULT_MIN_CELLS_PER_BIN = 5


def smooth_grid(values, valid, sigma, mode='zero_fill'):
    """Smooth valid bins; keep the manuscript's zero-fill behaviour by default."""
    import numpy as np
    from scipy.ndimage import gaussian_filter

    filled = np.where(valid, values, 0.0)
    if sigma and sigma > 0:
        filled = gaussian_filter(filled, sigma=float(sigma))
        if mode == 'mask_normalised':
            weight = gaussian_filter(np.asarray(valid, dtype=float), sigma=float(sigma))
            np.divide(filled, weight, out=filled, where=weight > 0)
    return np.where(valid, filled, np.nan)

# Share-state schema version at which the smoothing slider was introduced. A shared
# link stamped below this (or with no version at all = a pre-slider v2 link) predates
# the slider, so its smoothing was the fixed per-dataset config value of the day --
# reproduce that from LEGACY_SMOOTH_SIGMA rather than the current config default.
SLIDER_SCHEMA_VERSION = 3

# Per-dataset sigma in force BEFORE the slider existed. A pre-slider link encoded only
# smoothing on/off, not the sigma, so to reproduce it faithfully we fall back to these
# OLD config defaults, not the current ones (human_rpc was 0.5 then; its default is now
# 2.0). Unknown id -> historical 1.5.
LEGACY_SMOOTH_SIGMA = {
    'chick_rpc': 2.0,
    'chick_full': 1.5,        # had no configured sigma -> old update_data fallback (1.5)
    'human_rpc': 0.5,         # the old (near-no-op) default; new default is 2.0
    'mouse_rpc': 2.0,
    'mouse_rpc_legacy': 2.0,
}

def _coerce_sigma(value, fallback):
    """Clamp an externally-supplied sigma to the slider's [0, 4] domain.

    A hand-edited / corrupted ``?state=`` could carry a negative sigma (scipy's
    gaussian_filter raises), a huge one, or a non-number; coerce to a safe float,
    falling back to ``fallback`` when it is not numeric at all.
    """
    try:
        if isinstance(value, bool):
            return fallback
        if isinstance(value, int):
            return float(min(4, max(0, value)))
        sigma = float(value)
        return min(4.0, max(0.0, sigma)) if math.isfinite(sigma) else fallback
    except (TypeError, ValueError, OverflowError):
        return fallback


def resolve_smooth_sigma(state, dataset_id, config_default):
    """Pick the smoothing-slider value on dataset load / shared-link restore.

    explicit URL value (clamped to [0, 4]) > legacy per-dataset default (a pre-slider
    link, identified by its schema version) > current config default (plain dataset
    switch, or a current-schema link that carried no sigma). ``state`` is the URL state
    already scoped to this dataset (state_for_dataset), or None when not shared.
    """
    if state is not None and 'smooth_sigma' in state:
        return _coerce_sigma(state['smooth_sigma'], config_default)
    try:
        version = float(state.get('v', 2)) if state is not None else SLIDER_SCHEMA_VERSION
    except (TypeError, ValueError, OverflowError):
        version = SLIDER_SCHEMA_VERSION
    if state is not None and version < SLIDER_SCHEMA_VERSION:
        # Pre-slider link (no sigma key, older schema) -> what it rendered with then.
        return LEGACY_SMOOTH_SIGMA.get(dataset_id, 1.5)
    return config_default
