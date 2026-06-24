"""Gaussian-smoothing sigma resolution for the spatial map.

Smoothing strength is a live slider (v0.31+); the per-dataset value in
datasets_config.yml is just its default. Kept here -- pure, no Dash/data deps -- so
the backward-compat resolution is unit-testable without importing the Dash app.
"""

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

# Slider domain; an externally-supplied (URL) sigma is clamped to this.
SIGMA_MIN, SIGMA_MAX = 0.0, 4.0


def _coerce_sigma(value, fallback):
    """Clamp an externally-supplied sigma to the slider's [0, 4] domain.

    A hand-edited / corrupted ``?state=`` could carry a negative sigma (scipy's
    gaussian_filter raises), a huge one, or a non-number; coerce to a safe float,
    falling back to ``fallback`` when it is not numeric at all.
    """
    try:
        return min(SIGMA_MAX, max(SIGMA_MIN, float(value)))
    except (TypeError, ValueError):
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
    if state is not None and state.get('v', 2) < SLIDER_SCHEMA_VERSION:
        # Pre-slider link (no sigma key, older schema) -> what it rendered with then.
        return LEGACY_SMOOTH_SIGMA.get(dataset_id, 1.5)
    return config_default
