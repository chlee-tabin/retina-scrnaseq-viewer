"""Gaussian-smoothing sigma resolution for the spatial map.

Smoothing strength is a live slider (v0.31+); the per-dataset value in
datasets_config.yml is just its default. Kept here -- pure, no Dash/data deps -- so
the backward-compat resolution is unit-testable without importing the Dash app.
"""

# Per-dataset sigma in force BEFORE the slider existed (v0.31). A share link made then
# encoded only smoothing on/off -- not the sigma -- so to reproduce it faithfully we
# fall back to these OLD config defaults, not the current ones (human_rpc was 0.5 then;
# its default is now 2.0). Unknown id -> historical 1.5.
LEGACY_SMOOTH_SIGMA = {
    'chick_rpc': 2.0,
    'chick_full': 1.5,        # had no configured sigma -> old update_data fallback (1.5)
    'human_rpc': 0.5,         # the old (near-no-op) default; new default is 2.0
    'mouse_rpc': 2.0,
    'mouse_rpc_legacy': 2.0,
}


def resolve_smooth_sigma(state, dataset_id, config_default):
    """Pick the smoothing-slider value on dataset load / shared-link restore.

    explicit URL value (new link) > legacy per-dataset default (old link, no sigma
    key) > current config default (plain dataset switch). ``state`` is the URL state
    already scoped to this dataset (state_for_dataset), or None when not shared.
    """
    if state is not None and 'smooth_sigma' in state:
        return state['smooth_sigma']
    if state is not None:  # shared link predating the slider -> what it rendered with then
        return LEGACY_SMOOTH_SIGMA.get(dataset_id, 1.5)
    return config_default
