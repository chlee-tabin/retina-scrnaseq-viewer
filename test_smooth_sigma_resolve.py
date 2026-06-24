"""Self-check for the smoothing-slider sigma resolution + backward compat.

The one path most likely to regress: an OLD share link (no `smooth_sigma` key) must
still render human_rpc at the old 0.5, not the corrected 2.0 default.

Run: ~/repos/_retina_viewer_venv/bin/python test_smooth_sigma_resolve.py
"""
from utils.smoothing import resolve_smooth_sigma, LEGACY_SMOOTH_SIGMA

CFG = 2.0  # current human_rpc config default (post-fix)

# 1. Plain dataset switch (no shared state) -> current config default.
assert resolve_smooth_sigma(None, 'human_rpc', CFG) == 2.0

# 2. New share link with an explicit sigma -> that exact value wins.
assert resolve_smooth_sigma({'dataset': 'human_rpc', 'smooth_sigma': 1.3},
                            'human_rpc', CFG) == 1.3
assert resolve_smooth_sigma({'dataset': 'human_rpc', 'smooth_sigma': 0.0},
                            'human_rpc', CFG) == 0.0  # 0 is a real value, not "absent"

# 3. OLD share link (smoothing on, NO sigma key) -> the legacy per-dataset default,
#    NOT the new config default. This is the whole point of backward compat.
assert resolve_smooth_sigma({'dataset': 'human_rpc'}, 'human_rpc', CFG) == 0.5
assert resolve_smooth_sigma({'dataset': 'mouse_rpc'}, 'mouse_rpc', 2.0) == 2.0
assert resolve_smooth_sigma({'dataset': 'chick_rpc'}, 'chick_rpc', 2.0) == 2.0

# 4. OLD link for an unknown/new dataset id -> historical code fallback (1.5).
assert resolve_smooth_sigma({'dataset': 'zzz'}, 'zzz', 9.9) == 1.5

# 5. The legacy map preserves the buggy human value (so old links look unchanged).
assert LEGACY_SMOOTH_SIGMA['human_rpc'] == 0.5

print("OK: smoothing sigma resolution + backward compat")
