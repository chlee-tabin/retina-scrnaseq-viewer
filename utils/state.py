import json
import base64
from urllib.parse import parse_qs
from collections import namedtuple
from copy import deepcopy
from functools import partial
import os

from utils.validation import sanitize_state, coerce_control, coerce_roi, coerce_camera
from utils.smoothing import _coerce_sigma, LEGACY_SMOOTH_SIGMA, SLIDER_SCHEMA_VERSION

SCHEMA_VERSION = 4
StateField = namedtuple('StateField', 'component property default coerce introduced owner',
                       defaults=['url'])


def _text(value):
    return value if isinstance(value, str) else None


def _checklist(value):
    return ['enabled'] if isinstance(value, list) and 'enabled' in value else []


def _choice(choices, default, value):
    return value if isinstance(value, str) and value in choices else default


def _projection_default(key, dataset):
    from utils.wholemount import projection_controls
    return projection_controls(dataset)[key]


def _field(component, default=None, coerce=_text, introduced=1, owner='url', prop='value'):
    return StateField(component, prop, default, coerce, introduced, owner)


# One schema for capture, restore and defaults. Owners that also populate options
# retain their callback, but generate their value Outputs and read this same table.
STATE_SCHEMA = {
    'dataset': _field('dataset-select'),
    'mode': _field('viz-mode', 'random', partial(_choice, ('random', 'ordered_asc', 'ordered_desc'), 'random')),
    'bins': _field('bin-number-slider', 50, partial(coerce_control, 'bins')),
    'percentile': _field('percentile-slider', 0.95, partial(coerce_control, 'percentile')),
    'enable_binning': _field('enable-binning', ['enabled'], _checklist),
    'enable_smoothing': _field('enable-smoothing', ['enabled'], _checklist),
    'bin_stat': _field('bin-stat', 'mean', partial(_choice, ('mean', 'frac_pos'), 'mean')),
    'plot_type': _field('plot-type', 'map', partial(_choice, ('map', 'group'), 'map')),
    'compare_genes': _field('compare-genes', [], _checklist),
    'gene2': _field('gene-select-2'),
    'compare_shared_scale': _field('compare-shared-scale', [], _checklist),
    'custom_projection': _field('custom-projection', 'raw', partial(_choice, ('raw', 'flower', 'sphere'), 'raw')),
    'wm_rho_nt': _field('wm-rho-nt', partial(_projection_default, 'wm_rho_nt'), partial(coerce_control, 'wm_rho_nt')),
    'wm_rho_dv': _field('wm-rho-dv', partial(_projection_default, 'wm_rho_dv'), partial(coerce_control, 'wm_rho_dv')),
    'wm_gap': _field('wm-gap', partial(_projection_default, 'wm_gap'), partial(coerce_control, 'wm_gap')),
    'wm_stretch': _field('wm-stretch', partial(_projection_default, 'wm_stretch'), partial(coerce_control, 'wm_stretch')),
    'wm_cuts': _field('wm-cuts', partial(_projection_default, 'wm_cuts'), partial(coerce_control, 'wm_cuts')),
    'wm_symmetric': _field('wm-symmetric', partial(_projection_default, 'wm_symmetric'), _checklist),
    'wm_dewarp': _field('wm-dewarp', partial(_projection_default, 'wm_dewarp'), partial(_choice, ('arcsin', 'pow', 'none'), 'arcsin')),
    'wm_pow': _field('wm-pow', partial(_projection_default, 'wm_pow'), partial(coerce_control, 'wm_pow')),
    'wm_gap_mode': _field('wm-gap-mode', partial(_projection_default, 'wm_gap_mode'), partial(_choice, ('deficit', 'linear'), 'deficit')),
    'wm_gap_frac': _field('wm-gap-frac', partial(_projection_default, 'wm_gap_frac'), partial(coerce_control, 'wm_gap_frac')),
    'wm_pole': _field('wm-pole', partial(_projection_default, 'wm_pole'), partial(_choice, ('origin', 'median'), 'origin')),
    'group_value_source': _field('group-value-source', 'gene', partial(_choice, ('gene', 'meta'), 'gene')),
    'ng_layers': _field('nasal-gap-layers', partial(_projection_default, 'ng_layers'), partial(coerce_control, 'ng_layers'), 2),
    'ng_frac': _field('nasal-gap-frac', partial(_projection_default, 'ng_frac'), partial(coerce_control, 'ng_frac'), 2),
    'ng_dorsal': _field('nasal-gap-dorsal', partial(_projection_default, 'ng_dorsal'), partial(coerce_control, 'ng_dorsal'), 2),
    'ng_ventral': _field('nasal-gap-ventral', partial(_projection_default, 'ng_ventral'), partial(coerce_control, 'ng_ventral'), 2),
    'haa_mode': _field('haa-mode', 'domain', partial(_choice, ('off', 'footprint', 'expression', 'domain', 'peak'), 'domain'), 2),
    'sphere_cam': _field('sphere-camera-store', None, coerce_camera, 2, prop='data'),
    'roi': _field('roi-vertices-store', {'A': [], 'B': []}, coerce_roi, 2, prop='data'),
    'volcano_lfc': _field('volcano-lfc-thresh', 1.0, partial(coerce_control, 'volcano_lfc'), 2),
    'volcano_padj': _field('volcano-padj-thresh', 0.05, partial(coerce_control, 'volcano_padj'), 2),
    'deg_min_cells': _field('deg-min-cells', 50, partial(coerce_control, 'deg_min_cells'), 2),
    'smoothing_mode': _field('smoothing-mode', 'zero_fill', partial(_choice, ('zero_fill', 'mask_normalised'), 'zero_fill'), 4),
    'embedding': _field('embedding-select', lambda ds: ds.get('default_embedding'), owner='data'),
    'color': _field('color-select', 'gene_expression', owner='data'),
    'gene': _field('gene-select', lambda ds: ds.get('default_gene'), owner='data'),
    'smooth_sigma': _field('smooth-sigma-slider', lambda ds: ds.get('smooth_sigma', 1.5),
                           lambda value: _coerce_sigma(value, None), SLIDER_SCHEMA_VERSION, 'data'),
    'custom_x': _field('custom-x-select', 'NT.Score', owner='axes'),
    'custom_y': _field('custom-y-select', 'DV.Score', owner='axes'),
    'group_by': _field('group-by-select', owner='group'),
    'group_split': _field('group-split-select', '', owner='group'),
    'group_gene': _field('group-gene-select', lambda ds: ds.get('default_gene'), owner='group'),
    'group_style': _field('group-style', 'figure', partial(_choice, ('figure', 'violin', 'box', 'strip', 'dotplot'), 'figure'), owner='group'),
    'gene_module': _field('gene-module-select', owner='group'),
    'group_positive_only': _field('group-positive-only', ['enabled'], _checklist, owner='group'),
    'group_replicate': _field('group-replicate-select', owner='group'),
    'group_meta': _field('group-meta-select', owner='group'),
    'data_version': StateField(None, None, None, _text, 4, 'metadata'),
    'v': StateField(None, None, SCHEMA_VERSION, lambda value: _version({'v': value}), 1, 'metadata'),
}


def schema_keys(owner=None):
    return [key for key, field in STATE_SCHEMA.items()
            if field.component and (owner is None or field.owner == owner)]


def share_states():
    from dash import State
    return [State(STATE_SCHEMA[key].component, STATE_SCHEMA[key].property)
            for key in schema_keys()]


def restore_output(key, allow_duplicate=True):
    from dash import Output
    field = STATE_SCHEMA[key]
    return Output(field.component, field.property, allow_duplicate=allow_duplicate)


def restore_outputs(owner):
    return [restore_output(key) for key in schema_keys(owner)]


def restore_state(state, dataset=None):
    """Omissions reset to schema defaults, including documented historical defaults."""
    dataset = dict(dataset or {}, dataset_id=state.get('dataset'))
    result = {}
    for key, field in STATE_SCHEMA.items():
        default = field.default(dataset) if callable(field.default) else deepcopy(field.default)
        result[key] = field.coerce(state[key]) if key in state else default
    # Also used by direct callers, although decoded legacy links already migrate here.
    if 'smooth_sigma' not in state and _version(state) < STATE_SCHEMA['smooth_sigma'].introduced:
        result['smooth_sigma'] = LEGACY_SMOOTH_SIGMA.get(state.get('dataset'), 1.5)
    return result


def restore_values(state, owner, dataset=None):
    values = restore_state(state, dataset)
    return tuple(values[key] for key in schema_keys(owner))


def _version(state):
    try:
        return int(state.get('v', 2))
    except (TypeError, ValueError, OverflowError):
        return SCHEMA_VERSION


def data_version_notice(state, dataset):
    old = state.get('data_version') if state else None
    new = os.path.basename(dataset['file_path'])
    if old and old != new:
        return (f"This link was created with an earlier version of this dataset ({old}); "
                f"it now shows {new}.")
    return ''


def capture_state(values, dataset):
    state = dict(zip(schema_keys(), values))
    state = {key: STATE_SCHEMA[key].coerce(value) for key, value in state.items()}
    state.update(v=SCHEMA_VERSION, data_version=os.path.basename(dataset['file_path']))
    state['roi'] = {key: [[round(x, 4), round(y, 4)] for x, y in points]
                    for key, points in state['roi'].items()}
    return state

def encode_state(state_dict):
    """Encode application state for URL sharing"""
    state_json = json.dumps(state_dict)
    return base64.urlsafe_b64encode(state_json.encode()).decode()

def decode_state(encoded_state):
    """Decode application state from URL"""
    try:
        state_json = base64.urlsafe_b64decode(encoded_state).decode()
        state = sanitize_state(json.loads(state_json))
        if state is None:
            return None
        # v2 -> v3: the old fixed smoothing strength is part of the link's identity.
        if _version(state) < STATE_SCHEMA['smooth_sigma'].introduced:
            state.setdefault('smooth_sigma', LEGACY_SMOOTH_SIGMA.get(state.get('dataset'), 1.5))
        return state
    except Exception:
        return None

def create_share_url(base_url, state_dict):
    """Create a shareable URL carrying the full encoded state (base64 JSON).

    The caller (callbacks.url_callbacks.share_url) decides which keys to include;
    this function no longer strips or derives any of them, so the expression-by-group
    and figure controls survive the round-trip.
    """
    encoded_state = encode_state(state_dict)
    return f"{base_url}?state={encoded_state}"

def parse_url_state(url_search):
    """Parse state from URL search parameters"""
    if not url_search:
        return None
    
    params = parse_qs(url_search.lstrip('?'))
    if 'state' not in params:
        return None

    return decode_state(params['state'][0])


def state_for_dataset(url_search, dataset_id):
    """Parse shared-URL state, but only when it belongs to ``dataset_id`` (else None).

    The single seam for the stale-state guard: a ``?state=`` from a different dataset
    must not re-apply when the user switches datasets (the URL is never cleared). Every
    restore callback that reads state for an ALREADY-LOADED dataset parses through this
    rather than bare ``parse_url_state``; only ``initialize_from_url`` (which SETS the
    dataset, so it has none to compare against) uses the bare parse, by design.
    """
    state = parse_url_state(url_search)
    if state and state.get('dataset') != dataset_id:
        return None
    return state
