"""Small server-side validators for URL state and forged callback inputs."""
import math

import pandas as pd

from utils.smoothing import _coerce_sigma


# (minimum, maximum, default, integer). Keep in sync with control_panel / deg_panel.
# The LFC input has no UI maximum; require a finite value and honour its minimum.
NUMERIC_CONTROLS = {
    'bins': (2, 120, 50, True),
    'percentile': (0, 1, 0.95, False),
    'wm_rho_nt': (40, 150, 82, False),
    'wm_rho_dv': (40, 150, 64, False),
    'wm_cuts': (0, 8, 4, True),
    'wm_gap': (0, 1.5, 1.0, False),
    'wm_gap_frac': (0, 1, 0.5, False),
    'wm_stretch': (0, 1.5, 1.0, False),
    'wm_pow': (1, 3, 1.6, False),
    'ng_layers': (0, 5, 2, True),
    'ng_frac': (0.04, 0.30, 0.13, False),
    'ng_dorsal': (0, 90, 45, False),
    'ng_ventral': (0, 90, 55, False),
    'volcano_lfc': (0, None, 1.0, False),
    'volcano_padj': (0, 1, 0.05, False),
    'deg_min_cells': (5, 100000, 50, True),
}
CATEGORY_LIMIT = 200


def coerce_number(value, minimum, maximum, fallback, integer=False):
    """Coerce finite numbers/numeric strings, clamp, then optionally truncate to int."""
    try:
        if isinstance(value, bool):
            return fallback
        # Compare Python integers before float conversion (even a huge int can be clamped).
        number = value if isinstance(value, int) else float(value)
        if not isinstance(number, int) and not math.isfinite(number):
            return fallback
        number = max(minimum, number)
        if maximum is not None:
            number = min(maximum, number)
        return int(number) if integer else float(number)
    except (TypeError, ValueError, OverflowError):
        return fallback


def coerce_control(name, value):
    if name == 'smooth_sigma':
        return _coerce_sigma(value, 2.0)
    if name == 'deg_min_cells':
        try:
            if float(value) == 0:
                return 50
        except (TypeError, ValueError, OverflowError):
            pass
    lo, hi, default, integer = NUMERIC_CONTROLS[name]
    return coerce_number(value, lo, hi, default, integer)


def coerce_roi(value):
    """Reject the entire ROI on malformed structure, non-finite points or excess vertices."""
    empty = {'A': [], 'B': []}
    if not isinstance(value, dict) or set(value) - {'A', 'B'}:
        return empty
    result = {}
    for key in ('A', 'B'):
        polygon = value.get(key, [])
        if not isinstance(polygon, list) or len(polygon) > 500:
            return empty
        points = []
        for pair in polygon:
            if not isinstance(pair, list) or len(pair) != 2:
                return empty
            try:
                if any(isinstance(v, bool) or not isinstance(v, (int, float))
                       or not math.isfinite(v) for v in pair):
                    return empty
                points.append([float(pair[0]), float(pair[1])])
            except (ValueError, OverflowError):
                return empty
        result[key] = points
    return result


def coerce_camera(value):
    """Allow only Plotly camera fields. No UI range exists; bound vectors to [-10, 10]."""
    if not isinstance(value, dict) or set(value) - {'eye', 'center', 'up', 'projection'}:
        return None
    result = {}
    for key, vector in value.items():
        if key == 'projection':
            if (not isinstance(vector, dict) or set(vector) != {'type'}
                    or vector['type'] not in ('perspective', 'orthographic')):
                return None
            result[key] = dict(vector)
            continue
        if not isinstance(vector, dict) or set(vector) != {'x', 'y', 'z'}:
            return None
        coords = {axis: coerce_number(v, -10, 10, None) for axis, v in vector.items()}
        if any(v is None for v in coords.values()):
            return None
        result[key] = coords
    return result or None


def sanitize_state(state):
    if not isinstance(state, dict):
        return None
    result = dict(state)
    for name in NUMERIC_CONTROLS:
        if name in result:
            result[name] = coerce_control(name, result[name])
    if 'smooth_sigma' in result:
        # Invalid URL sigma remains an explicit null so resolve_smooth_sigma uses
        # this dataset's config default, rather than a fixed cross-dataset fallback.
        result['smooth_sigma'] = _coerce_sigma(result['smooth_sigma'], None)
    if 'roi' in result:
        result['roi'] = coerce_roi(result['roi'])
    if 'sphere_cam' in result:
        result['sphere_cam'] = coerce_camera(result['sphere_cam'])
    return result


def is_categorical_series(series):
    return (not pd.api.types.is_numeric_dtype(series)
            or pd.api.types.is_bool_dtype(series)
            or (pd.api.types.is_integer_dtype(series) and series.nunique(dropna=False) <= 50))


def category_allowed(series):
    n = len(series.cat.categories) if isinstance(series.dtype, pd.CategoricalDtype) else 0
    return max(n, series.nunique(dropna=False)) <= CATEGORY_LIMIT


def obs_column_types(obs):
    """Exclude high-cardinality categorical columns from both colour and group choices."""
    return {col: ('categorical' if is_categorical_series(s) else 'numeric')
            for col, s in obs.items()
            if not is_categorical_series(s) or category_allowed(s)}
