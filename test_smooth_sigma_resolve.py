"""Self-checks for the smoothing-slider sigma resolution + backward compat.

The path most likely to regress: a pre-slider share link (schema v<3) must still
render human_rpc at the old 0.5, not the corrected 2.0 default.

Run: ~/repos/_retina_viewer_venv/bin/python test_smooth_sigma_resolve.py
(or via pytest)
"""
from utils.smoothing import (resolve_smooth_sigma, LEGACY_SMOOTH_SIGMA,
                             SLIDER_SCHEMA_VERSION)

CFG = 2.0  # current human_rpc config default (post-fix)


def test_plain_switch_uses_config_default():
    # No shared state (a plain dataset switch) -> current config default.
    assert resolve_smooth_sigma(None, 'human_rpc', CFG) == 2.0


def test_explicit_url_value_wins():
    v3 = SLIDER_SCHEMA_VERSION
    assert resolve_smooth_sigma({'v': v3, 'smooth_sigma': 1.3}, 'human_rpc', CFG) == 1.3
    # 0 is a real value, not "absent".
    assert resolve_smooth_sigma({'v': v3, 'smooth_sigma': 0.0}, 'human_rpc', CFG) == 0.0


def test_legacy_link_uses_old_per_dataset_default():
    # Pre-slider link (older schema, no sigma key) -> the OLD default, NOT current config.
    # Real old links carry v=2; one without any v also defaults to <3.
    assert resolve_smooth_sigma({'v': 2, 'dataset': 'human_rpc'}, 'human_rpc', CFG) == 0.5
    assert resolve_smooth_sigma({'dataset': 'human_rpc'}, 'human_rpc', CFG) == 0.5
    assert resolve_smooth_sigma({'v': 2, 'dataset': 'mouse_rpc'}, 'mouse_rpc', 2.0) == 2.0
    assert resolve_smooth_sigma({'v': 2, 'dataset': 'chick_rpc'}, 'chick_rpc', 2.0) == 2.0
    # Unknown/new dataset in an old link -> historical code fallback (1.5).
    assert resolve_smooth_sigma({'v': 2, 'dataset': 'zzz'}, 'zzz', 9.9) == 1.5


def test_current_schema_link_without_sigma_uses_config_not_legacy():
    # A v3 link that carried no sigma (e.g. shared from a non-spatial view) is NOT a
    # pre-slider link: use the current config default, never the legacy 0.5.
    assert resolve_smooth_sigma({'v': SLIDER_SCHEMA_VERSION, 'dataset': 'human_rpc'},
                                'human_rpc', CFG) == 2.0


def test_url_sigma_is_clamped_and_sanitised():
    v3 = SLIDER_SCHEMA_VERSION
    assert resolve_smooth_sigma({'v': v3, 'smooth_sigma': 40}, 'human_rpc', CFG) == 4.0    # over max
    assert resolve_smooth_sigma({'v': v3, 'smooth_sigma': -3}, 'human_rpc', CFG) == 0.0    # below min
    assert resolve_smooth_sigma({'v': v3, 'smooth_sigma': '1.5'}, 'human_rpc', CFG) == 1.5  # numeric str
    # Non-numeric / null -> fall back to config default, not a silent 0.
    assert resolve_smooth_sigma({'v': v3, 'smooth_sigma': 'abc'}, 'human_rpc', CFG) == 2.0
    assert resolve_smooth_sigma({'v': v3, 'smooth_sigma': None}, 'human_rpc', CFG) == 2.0


def test_legacy_map_preserves_buggy_human_value():
    # So old links look unchanged.
    assert LEGACY_SMOOTH_SIGMA['human_rpc'] == 0.5


if __name__ == '__main__':
    test_plain_switch_uses_config_default()
    test_explicit_url_value_wins()
    test_legacy_link_uses_old_per_dataset_default()
    test_current_schema_link_without_sigma_uses_config_not_legacy()
    test_url_sigma_is_clamped_and_sanitised()
    test_legacy_map_preserves_buggy_human_value()
    print("OK: smoothing sigma resolution + backward compat + clamp + version-gate")
