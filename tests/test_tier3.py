"""Tier 3 regressions, synthetic only. Also runnable: python tests/test_tier3.py."""
import hashlib
import inspect
import io
import os
from pathlib import Path
import sys
import tarfile
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ['HF_DATA_REPO'] = ''

import anndata as ad
import numpy as np
import pandas as pd
from dash import no_update
from dash.exceptions import PreventUpdate
import yaml

from callbacks import (dataset_callbacks as dataset_cb, main_callbacks as main,
                       url_callbacks as url, deg_callbacks as de_cb)
from scripts.deploy_space import deploy_space
from utils import data_loading as loading, data_provision as provision, deg, plotting, wholemount
from utils.error_handling import handle_callback_error, log_callback_info
from utils.deg_plots import create_volcano_figure
from utils.smoothing import DEFAULT_MIN_CELLS_PER_BIN, smooth_grid, resolve_smooth_sigma
from utils.state import (STATE_SCHEMA, SCHEMA_VERSION, capture_state, create_share_url,
                         data_version_notice, decode_state, encode_state, parse_url_state,
                         restore_outputs, restore_state, restore_values, schema_keys, share_states)
from utils.validation import obs_column_types, is_categorical_series
from utils.version import viewer_revision


def test_every_schema_key_roundtrips():
    dataset = loading.get_dataset_config('human_rpc')
    values = restore_state({'dataset': 'human_rpc', 'v': SCHEMA_VERSION}, dataset)
    values.update(gene='VAX1', gene2='FGF8', mode='ordered_desc', bins=67,
                  percentile=.7, compare_genes=['enabled'], compare_shared_scale=['enabled'],
                  enable_smoothing=[], custom_projection='sphere', smoothing_mode='mask_normalised',
                  smooth_sigma=1.3, roi={'A': [[0., 0.], [1., 0.], [0., 1.]], 'B': []},
                  sphere_cam={'eye': {'x': 1., 'y': 2., 'z': 3.}},
                  group_gene='VAX1', group_by='library', group_split='library',
                  group_style='violin', gene_module='test', group_positive_only=[],
                  group_replicate='library', group_value_source='meta', group_meta='DV.Score')
    state = capture_state([values[key] for key in schema_keys()], dataset)
    decoded = parse_url_state(create_share_url('https://example.test/', state).split('?', 1)[1])
    restored = restore_state(decoded, dataset)
    assert set(decoded) == set(STATE_SCHEMA) | {'v'}
    for key in STATE_SCHEMA:
        assert restored[key] == state[key], key
    assert decoded['v'] == SCHEMA_VERSION
    for key, dash_state in zip(schema_keys(), share_states()):
        field = STATE_SCHEMA[key]
        assert (dash_state.component_id, dash_state.component_property) == (field.component, field.property)
    for owner in ('url', 'data', 'axes', 'group'):
        outputs = restore_outputs(owner)
        vals = restore_values(decoded, owner, dataset)
        assert len(outputs) == len(vals) == len(schema_keys(owner))
        for key, output, value in zip(schema_keys(owner), outputs, vals):
            assert output.component_id == STATE_SCHEMA[key].component
            assert output.component_property == STATE_SCHEMA[key].property
            assert value == state[key]


def test_v2_migration_and_omitted_key_defaults():
    dataset = loading.get_dataset_config('human_rpc')
    for legacy in ({'v': 2, 'dataset': 'human_rpc'}, {'dataset': 'human_rpc'}):
        migrated = decode_state(encode_state(legacy))
        assert migrated['smooth_sigma'] == .5
        assert resolve_smooth_sigma(migrated, 'human_rpc', 2.) == .5
        defaults = restore_state(migrated, dataset)
        assert defaults['smoothing_mode'] == 'zero_fill'
        assert defaults['sphere_cam'] is None
    state = decode_state(encode_state({'v': 3, 'dataset': 'human_rpc', 'gene': 'VAX1'}))
    assert 'smooth_sigma' not in state
    defaults = restore_state(state, dataset)
    assert defaults['compare_genes'] == defaults['compare_shared_scale'] == []
    assert defaults['gene2'] is None
    assert defaults['roi'] == {'A': [], 'B': []}
    assert defaults['smooth_sigma'] == 2.
    with patch.object(url, 'ctx', SimpleNamespace(triggered_id='url.search')):
        restored = url.initialize_from_url('?state=' + encode_state(state), '/')
    by_key = dict(zip(schema_keys('url'), restored))
    for key in ('compare_genes', 'compare_shared_scale', 'gene2', 'roi', 'bins', 'percentile'):
        assert by_key[key] == defaults[key] and by_key[key] is not no_update
    with patch.object(url, 'ctx', SimpleNamespace(triggered_prop_ids={'url.pathname': 'url'})):
        untouched = url.initialize_from_url('?state=' + encode_state(state), '/')
    assert all(value is no_update for value in untouched)


def test_share_callback_resolves_server_data_version():
    dataset = loading.get_dataset_config('chick_rpc')
    values = restore_state({'v': SCHEMA_VERSION, 'dataset': 'chick_rpc'}, dataset)
    _, link = url.share_url(1, *[values[key] for key in schema_keys()], 'https://example.test/?old')
    state = parse_url_state('?' + link.split('?', 1)[1])
    assert state['data_version'] == '20250604_chick_RPC.h5ad'
    assert state['smoothing_mode'] == 'zero_fill'
    assert data_version_notice(state, dataset) == ''
    assert data_version_notice({'dataset': 'chick_rpc'}, dataset) == ''
    assert data_version_notice(dict(state, data_version='older.bin'), dataset) == (
        'This link was created with an earlier version of this dataset (older.bin); '
        'it now shows 20250604_chick_RPC.h5ad.')


def test_legacy_v3_restore_keeps_gene_dropdown_value_and_v4_roundtrips():
    """A v3 whole-mount link must survive each callback that owns gene options."""
    import app

    legacy = {
        'v': 3,
        'dataset': 'chick_rpc',
        'gene': 'CYP26C1',
        'bins': 40,
        'smooth_sigma': 1.0,
        'custom_projection': 'flower',
        'enable_binning': ['enabled'],
        'data_version': '20240101_old_chick.h5ad',
    }
    dataset = loading.get_dataset_config('chick_rpc')
    adata = ad.AnnData(np.ones((2, 3)),
                       obs=pd.DataFrame({'NT.Score': [0., 1.], 'DV.Score': [0., 1.]}),
                       var=pd.DataFrame(index=['FGF8', 'CYP26C1', 'VAX1']))

    def assert_restored(search):
        with patch.object(url, 'ctx', SimpleNamespace(triggered_prop_ids={'url.search': 'url'})):
            url_values = dict(zip(schema_keys('url'), url.initialize_from_url(search, '/')))
        with patch.object(app, 'load_adata', return_value=adata), \
             patch.object(app, 'dataset_column_types', return_value=obs_column_types(adata.obs)), \
             patch.object(app, 'dataset_norm_target', return_value=None):
            data_values = app.update_data('chick_rpc', search, None, None)
        data_store, _, _, color, gene, gene_options = data_values[:6]
        _, search_options = dataset_cb.update_gene_select(
            color, data_store, None, None, search)
        assert color == 'gene_expression'
        assert gene == 'CYP26C1'
        assert url_values['custom_projection'] == 'flower'
        assert {'label': 'CYP26C1', 'value': 'CYP26C1'} in gene_options
        assert {'label': 'CYP26C1', 'value': 'CYP26C1'} in search_options
        return data_store

    legacy_search = '?state=' + encode_state(legacy)
    assert_restored(legacy_search)

    # Capture the restored state exactly as the Share View callback does, then prove the
    # new v4 representation drives the same owners and selected dropdown value.
    restored = restore_state(parse_url_state(legacy_search), dataset)
    v4_state = capture_state([restored[key] for key in schema_keys()], dataset)
    assert v4_state['v'] == 4
    assert_restored('?state=' + encode_state(v4_state))


def test_checksum_cache_and_mismatch_deletion():
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        path = Path(tmp) / 'tiny.bin'
        path.write_bytes(b'verified data')
        expected = hashlib.sha256(path.read_bytes()).hexdigest()
        provision._verified_files.clear()
        assert provision.verify_sha256(path, expected)
        with patch.object(provision.hashlib, 'sha256', side_effect=AssertionError('rehash')):
            assert provision.verify_sha256(path, expected)
            provision._verified_files.clear()  # persistent stamp also survives a process restart
            assert provision.verify_sha256(path, expected)
        path.write_bytes(b'changed data with a different size')
        with patch.object(provision.hashlib, 'sha256', wraps=hashlib.sha256) as digest:
            changed = hashlib.sha256(path.read_bytes()).hexdigest()
            digest.reset_mock()
            assert provision.verify_sha256(path, changed)
            assert digest.call_count == 1
        try:
            provision.verify_sha256(path, '0' * 64)
        except ValueError as error:
            assert 'SHA256 mismatch' in str(error)
        else:
            assert False, 'bad checksum accepted'
        assert not path.exists() and not Path(str(path) + '.sha256.json').exists()


def test_path_rule_nested_provisioning_revision_and_local_verification():
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        root = Path(tmp) / 'datasets'
        dest = root / 'sub' / 'tiny.bin'
        data = b'pinned download'
        config = {'hf_data_revision': 'pinned', 'datasets': {'tiny': {
            'file_path': 'data/sub/tiny.bin', 'sha256': hashlib.sha256(data).hexdigest()}}}

        def download(**kwargs):
            target = Path(kwargs['local_dir']) / kwargs['filename']
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)

        with patch.dict(os.environ, {'DATA_DIR': str(root), 'HF_DATA_REPO': 'owner/data', 'HF_DATA_REVISION': ''}), \
             patch.object(loading, 'load_dataset_config', return_value=config), \
             patch('huggingface_hub.hf_hub_download', side_effect=download) as fetch:
            assert loading.resolve_data_path('data/sub/tiny.bin') == str(dest)
            assert loading.resolve_data_path('sub/tiny.bin') == str(dest)
            assert provision.ensure_one_dataset('data/sub/tiny.bin') == str(dest)
            assert fetch.call_args.kwargs['filename'] == 'sub/tiny.bin'
            assert fetch.call_args.kwargs['revision'] == 'pinned'
            with patch.object(provision, 'verify_sha256', wraps=provision.verify_sha256) as verify:
                assert provision.ensure_one_dataset('data/sub/tiny.bin') == str(dest)
                assert verify.call_count == 1 and fetch.call_count == 1
            dest.unlink()
            with patch.dict(os.environ, {'HF_DATA_REVISION': 'override'}):
                provision.ensure_one_dataset('data/sub/tiny.bin')
                assert fetch.call_args.kwargs['revision'] == 'override'
            # A stale provisioned copy is replaced by the pinned download (self-heal).
            dest.write_bytes(b'stale copy')
            assert provision.ensure_one_dataset('data/sub/tiny.bin') == str(dest)
            assert dest.read_bytes() == data
        # Local development (no HF_DATA_REPO): a developer's own file is never
        # checksum-verified or deleted, even if it differs from the pinned data.
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b'synthetic local copy')
        with patch.dict(os.environ, {'DATA_DIR': str(root), 'HF_DATA_REPO': ''}), \
             patch.object(loading, 'load_dataset_config', return_value=config):
            assert provision.ensure_one_dataset('data/sub/tiny.bin') == str(dest)
        assert dest.read_bytes() == b'synthetic local copy'


def test_bad_download_is_deleted_and_load_refuses_it():
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        config = {'datasets': {'tiny': {'file_path': 'data/tiny.bin', 'sha256': '0' * 64}}}

        def bad_download(**kwargs):
            (Path(kwargs['local_dir']) / kwargs['filename']).write_bytes(b'bad')

        loading.load_adata.cache_clear()
        with patch.dict(os.environ, {'DATA_DIR': tmp, 'HF_DATA_REPO': 'owner/data'}), \
             patch.object(loading, 'load_dataset_config', return_value=config), \
             patch('huggingface_hub.hf_hub_download', side_effect=bad_download), \
             patch.object(loading.ad, 'read_h5ad', side_effect=AssertionError('reader reached')):
            try:
                loading.load_adata('data/tiny.bin')
            except ValueError as error:
                assert 'SHA256 mismatch' in str(error)
            else:
                assert False, 'bad download loaded'
            assert not (Path(tmp) / 'tiny.bin').exists()
        loading.load_adata.cache_clear()


def test_duplicate_basename_rejected_at_config_load():
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        config = Path(tmp) / 'config.yml'
        config.write_text(yaml.safe_dump({'datasets': {'one': {'file_path': 'data/a/x.h5ad'},
                                                      'two': {'file_path': 'data/b/x.h5ad'}}}))
        try:
            loading.load_dataset_config(config)
        except ValueError as error:
            assert "Duplicate dataset basename 'x.h5ad'" in str(error)
        else:
            assert False, 'duplicate basename accepted'


def test_nested_download_route_uses_loader_path():
    import app
    config = {'datasets': {'nested': {'file_path': 'data/sub/synthetic.h5ad'}}}
    with patch.dict(os.environ, {'DATA_DIR': str(ROOT / 'synthetic-data')}), \
         patch.object(app, 'load_dataset_config', return_value=config), \
         patch.object(app.os.path, 'isfile', return_value=True), \
         patch.object(app, 'send_file', return_value='nested download') as send:
        response = app.server.test_client().get('/download/synthetic.h5ad')
        assert response.status_code == 200
        assert send.call_args.args[0] == loading.resolve_data_path('data/sub/synthetic.h5ad')


def test_mask_normalised_preserves_isolated_detected_bin():
    values = np.full((5, 5), np.nan)
    valid = np.zeros((5, 5), bool)
    values[2, 2], valid[2, 2] = 1., True
    normalised = smooth_grid(values, valid, 1.5, 'mask_normalised')
    legacy = smooth_grid(values, valid, 1.5)
    assert normalised[2, 2] == 1. and 0 < legacy[2, 2] < 1.
    assert np.isnan(normalised[~valid]).all()
    x = y = np.array([0., .5, .5, 1.])
    for mode, expected in (('mask_normalised', 1.), ('zero_fill', legacy[2, 2])):
        grid, counts, *_ = plotting._binned_mean(x, y, np.ones(4), bin_size=5,
                                                 min_cells=2, smooth_sigma=1.5,
                                                 stat='frac_pos', smoothing_mode=mode)
        assert np.isclose(grid[2, 2], expected)
        assert np.isnan(grid[counts < 2]).all()
    fig = plotting.create_sphere_binned_figure(y, x, np.ones(4), bin_size=5,
                                               min_cells=2, smooth_sigma=1.5, bin_stat='frac_pos',
                                               smoothing_mode='mask_normalised')
    mesh = next(t for t in fig.data if t.type == 'mesh3d')
    assert np.allclose(mesh.intensity, 1.)


def test_smoothing_mode_reaches_every_binned_dispatch_and_haa():
    obs = pd.DataFrame({'NT.Score': [-1., 0., 1.], 'DV.Score': [-1., 0., 1.]})
    adata = ad.AnnData(np.ones((3, 2)), obs=obs, var=pd.DataFrame(index=['g1', 'g2']))
    settings = dict(loading.get_dataset_config('chick_rpc'), dataset_id='chick_rpc',
                    column_types=obs_column_types(obs), haa_marker='g1')
    args = {name: None for name in inspect.signature(main.update_plot).parameters}
    args.update(data_store={'dataset_id': 'chick_rpc'}, embedding='custom_embedding',
                custom_x='NT.Score', custom_y='DV.Score', color_by='gene_expression', gene='g1',
                gene2='g2', enable_binning=['enabled'], enable_smoothing=['enabled'],
                bin_number=5, percentile=.95, smooth_sigma_ctrl=1.5, haa_mode='domain',
                smoothing_mode='mask_normalised', wm_symmetric=['enabled'])
    for projection, solo, dual in (
            ('raw', 'create_binned_plot', 'create_dual_gene_figure'),
            ('flower', 'create_wholemount_binned_figure', 'create_dual_gene_wholemount_figure'),
            ('sphere', 'create_sphere_binned_figure', 'create_dual_gene_sphere_figure')):
        for compare, renderer in (([], solo), (['enabled'], dual)):
            args.update(custom_projection=projection, compare_genes=compare)
            with patch.object(main, 'load_dataset_state', return_value=(adata, settings)), \
                 patch.object(main, renderer, return_value=plotting.go.Figure()) as render, \
                 patch.object(main, '_haa_center', return_value=None) as haa:
                main.update_plot(**args)
                assert render.call_args.kwargs['smoothing_mode'] == 'mask_normalised'
                if haa.called:
                    assert haa.call_args.kwargs['smoothing_mode'] == 'mask_normalised'


def test_wholemount_dataset_defaults_and_reset():
    for dataset_id, dataset in loading.load_dataset_config()['datasets'].items():
        params = wholemount.dataset_params(dict(dataset, dataset_id=dataset_id))
        if dataset_id.startswith('chick'):
            assert np.array_equal(params['stretch_bumps'], wholemount.LOCKED_PARAMS['stretch_bumps'])
        else:
            assert not params['stretch_bumps']
            live = main._wholemount_params(None, None, None, 1., None,
                                           dataset=dict(dataset, dataset_id=dataset_id))
            assert live['stretch_bumps'] == ()
    custom = {'wholemount_params': {'rho_max_nt_deg': 100, 'gap_gain': .7},
              'nasal_gap': {'layers': 3, 'frac': .2}}
    with patch.object(main, 'get_dataset_config', return_value=custom), \
         patch.object(main, 'ctx', SimpleNamespace(triggered_id='wm-reset')):
        controls = main.reset_wholemount_params(1, {'dataset_id': 'human_rpc'},
                                                '?state=' + encode_state({'dataset': 'human_rpc', 'wm_rho_nt': 50}))
    result = dict(zip(main.PROJECTION_KEYS, controls))
    assert result['wm_rho_nt'] == 100 and result['wm_gap'] == .7 and result['ng_layers'] == 3


def test_figure_region_half_open_edges_and_rounded_restore():
    edges = np.linspace(0, 1, 11)
    x = np.array([0., edges[6], edges[7], edges[8], 1., np.nan, edges[6], edges[7]])
    y = np.array([0., edges[4], edges[5], edges[4], 1., edges[4], edges[6], edges[4]])
    spec = {'n': 10, 'nt': [6, 7], 'dv': [4, 5]}
    assert deg.figure_region_indices(spec, x, y) == [1, 2, 7]
    rectangle = deg.figure_region_rectangle(spec, x, y)
    with patch.object(deg, 'polygon_to_indices', side_effect=AssertionError('polygon used')):
        assert deg.roi_indices(rectangle, x, y, [spec]) == [1, 2, 7]
        rounded = np.round(rectangle, 4).tolist()
        assert deg.roi_indices(rounded, x, y, [spec]) == [1, 2, 7]
    whole = {'n': 10, 'nt': [0, 9], 'dv': [0, 9]}
    assert 4 not in deg.figure_region_indices(whole, x, y)  # axis max excluded


def test_volcano_labels_na_note_and_finite_rest():
    res = pd.DataFrame({'gene': ['a', 'b', 'c'], 'log2FoldChange': [1., -2., 3.],
                        'padj': [.01, .02, np.nan]})
    for mode, solo, contrast in (('A_vs_B', 'A', 'ROI A / ROI B'),
                                 ('A_vs_rest', 'A', 'ROI A / rest'),
                                 ('A_vs_rest', 'B', 'ROI B / rest')):
        fig = create_volcano_figure(res, mode=mode, foreground_label=solo)
        assert fig.layout.xaxis.title.text == f'log2FC ({contrast})'
        assert '1 genes with padj = NA not shown' in fig.layout.title.text
    assert 'padj = NA' not in create_volcano_figure(res.iloc[:2]).layout.title.text
    finite = np.array([True, True, False, False, True])
    for a, b, solo in (([0], [], 'A'), ([], [0], 'B')):
        labels, info = deg.resolve_rois(a, b, 5, finite)
        assert list(labels) == ['A', 'B', '', '', 'B']
        assert info['solo'] == solo and info['n_B'] == 2


def test_gene_search_and_species_match():
    genes = ['Foxd1', 'Vax1', 'Fgf8']
    assert loading.match_gene('VAX1', genes) == 'Vax1'
    assert loading.match_gene('Vax1', genes) == 'Vax1'
    assert loading.match_gene('absent', genes) is None
    for blank in (None, '', ' '):
        assert loading.gene_search_options(genes, blank) == []
    assert loading.gene_search_options(genes, None, 'VAX1') == [
        {'label': 'Vax1', 'value': 'Vax1'}]
    assert loading.gene_search_options(genes, 'v') == [{'label': 'Vax1', 'value': 'Vax1'}]
    matches = loading.gene_search_options([f'gene{i}' for i in range(300)], 'g')
    assert len(matches) == 200
    store = {'dataset_id': 'mouse_rpc', 'genes': genes}
    assert main.restore_comparison_gene(store, 'VAX1', '') == 'Vax1'
    link = '?state=' + encode_state({'v': 4, 'dataset': 'mouse_rpc', 'gene': 'VAX1'})
    assert main.restore_comparison_gene(store, 'Fgf8', link) is None


def test_dataset_switch_and_shared_missing_gene_notices():
    import app
    obs = pd.DataFrame({'NT.Score': [0., 1.], 'DV.Score': [0., 1.]})
    adata = ad.AnnData(np.ones((2, 2)), obs=obs, var=pd.DataFrame(index=['Vax1', 'Fgf8']))
    with patch.object(app, 'load_adata', return_value=adata), \
         patch.object(app, 'dataset_column_types', return_value=obs_column_types(obs)), \
         patch.object(app, 'dataset_norm_target', return_value=None):
        assert app.update_data('mouse_rpc', '', None, 'VAX1')[4] == 'Vax1'
        state = {'v': 4, 'dataset': 'mouse_rpc', 'gene': 'missing', 'data_version': 'old.bin'}
        out = app.update_data('mouse_rpc', '?state=' + encode_state(state), None, 'Vax1')
        assert out[4] == 'Fgf8' and out[5] == [{'label': 'Fgf8', 'value': 'Fgf8'}]
        text = str(out[6])
        assert 'gene missing not found in this dataset; showing Fgf8' in text
        assert 'This link was created with an earlier version' in text
        old = app.update_data('mouse_rpc', '?state=' + encode_state({'dataset': 'mouse_rpc', 'gene': 'VAX1'}), None, None)
        assert old[4] == 'Vax1' and old[6] == []


def test_shared_categorical_helper_and_replicate_fallback():
    obs = pd.DataFrame({'text': pd.Series(['a', 'b'], dtype='string'),
                        'small_int': [0, 1], 'float': [0., 1.]})
    assert is_categorical_series(obs['text']) and is_categorical_series(obs['small_int'])
    assert obs_column_types(obs) == {'text': 'categorical', 'small_int': 'categorical', 'float': 'numeric'}
    adata = SimpleNamespace(obs=pd.DataFrame({'genotype': ['d0', 'd1']}))
    settings = {'column_types': {'genotype': 'categorical'}, 'genes': [], 'dataset_id': 'test'}
    with patch.object(main, 'load_dataset_state', return_value=(adata, settings)):
        controls = main.populate_group_controls({'dataset_id': 'test'}, '', 'gene')
    assert controls[7] is None  # genotype alone is never the sample-hint fallback


def test_selectable_dataset_cache_size_and_unified_floor():
    selectable = loading.validate_datasets(loading.load_dataset_config())
    assert loading.load_adata.cache_info().maxsize == len(selectable)
    assert DEFAULT_MIN_CELLS_PER_BIN == 5
    assert inspect.signature(plotting._binned_mean).parameters['min_cells'].default == DEFAULT_MIN_CELLS_PER_BIN
    assert loading.get_dataset_config('chick_full')['smooth_sigma'] == 1.5


def test_repeated_dataset_derivations_use_one_cache_entry_per_filename():
    """Restore callbacks must not retain a fresh AnnData-derived result per page load."""
    filename = loading.get_dataset_config('chick_rpc')['file_path']
    adata = ad.AnnData(np.ones((8, 3)),
                       obs=pd.DataFrame({'cluster': [0, 1] * 4}),
                       var=pd.DataFrame(index=['A', 'B', 'C']))
    loading.load_adata.cache_clear()
    try:
        with patch.object(loading, 'load_adata', return_value=adata):
            for _ in range(6):
                assert loading.dataset_column_types(filename) == {'cluster': 'categorical'}
                assert isinstance(loading.dataset_norm_target(filename), float)
        assert loading.dataset_column_types.cache_info().currsize == 1
        assert loading.dataset_norm_target.cache_info().currsize == 1
    finally:
        loading.load_adata.cache_clear()


def test_session_ids_are_per_layout_and_polling_is_30_seconds():
    import app

    def component(layout, component_id):
        return next(child for child in layout.children if getattr(child, 'id', None) == component_id)

    layout = app.app._layout
    assert app.app._layout_is_function is False
    assert component(layout, 'session-id').storage_type == 'session'
    assert getattr(component(layout, 'session-id'), 'data', None) is None
    assert component(layout, 'interval-component').interval == 30000
    assert app.app.title == 'Retina scRNA-seq Pattern Viewer'

    # A function-valued Dash layout is evaluated once per /_dash-layout request.
    # Repeated fresh-page layout requests must use this one immutable component tree;
    # the clientside callback is responsible for the per-session UUID.
    client = app.server.test_client()
    for _ in range(6):
        response = client.get('/_dash-layout')
        assert response.status_code == 200
    assert app.app._layout is layout
    dependencies = client.get('/_dash-dependencies').get_json()
    session_callbacks = [item for item in dependencies
                         if item['output'].startswith('session-id.data')]
    assert len(session_callbacks) == 1
    assert session_callbacks[0]['clientside_function'] is not None


def test_callback_errors_preserve_outputs_and_logging_omits_args():
    @handle_callback_error
    def bad():
        raise RuntimeError('/private/server/path')

    with patch('utils.error_handling.logger.exception') as log:
        try:
            bad()
        except PreventUpdate:
            pass
        else:
            assert False, 'callback error did not preserve outputs'
        assert log.called

    @log_callback_info
    def ok(data_store):
        return 1

    with patch('utils.error_handling.logger.debug') as log:
        assert ok({'secret': 'do not log'}) == 1
        assert 'secret' not in str(log.call_args_list)


def test_deploy_refuses_dirty_tree_and_uploads_revision_stamp():
    with patch('scripts.deploy_space.subprocess.check_output', return_value=b' M app.py\n'), \
         patch('huggingface_hub.upload_folder') as upload:
        try:
            deploy_space('owner/space', ROOT)
        except RuntimeError as error:
            assert 'dirty' in str(error)
        else:
            assert False, 'dirty checkout deployed'
        assert not upload.called
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode='w') as tar:
        info = tarfile.TarInfo('README.md')
        payload = b'clean tracked checkout'
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))
    sha = '1234567890abcdef' * 2 + '12345678'

    def git(command, **kwargs):
        assert kwargs['cwd'] == ROOT
        if command[1] == 'status':
            return b''
        if command[1] == 'rev-parse':
            return sha.encode()
        return archive.getvalue()

    def upload(**kwargs):
        root = Path(kwargs['folder_path'])
        assert (root / 'viewer_revision.txt').read_text().strip() == sha
        assert (root / 'README.md').read_bytes() == payload
        assert kwargs['repo_type'] == 'space' and kwargs['repo_id'] == 'owner/space'
        return 'uploaded'

    with patch.dict(os.environ, {'HF_TOKEN': 'synthetic-token'}), \
         patch('scripts.deploy_space.subprocess.check_output', side_effect=git), \
         patch('huggingface_hub.upload_folder', side_effect=upload):
        assert deploy_space('owner/space', ROOT) == 'uploaded'
    with patch.dict(os.environ, {'VIEWER_GIT_SHA': sha}):
        assert viewer_revision() == sha[:7]


def test_gene_vector_never_builds_a_view_that_copies_raw():
    """adata[:, gene] subsets .raw by copying the whole raw matrix (~0.9 GB for chick
    RPC) per call; _gene_vector must read the .X column directly."""
    import scipy.sparse as sp
    X = sp.csr_matrix(np.arange(12, dtype=np.float32).reshape(4, 3))
    a = ad.AnnData(X=X, var=pd.DataFrame(index=['G0', 'G1', 'G2']))
    a.raw = a.copy()
    with patch.object(ad.AnnData, '__getitem__', side_effect=AssertionError('view built')):
        assert np.array_equal(main._gene_vector(a, 'G1'), np.array([1, 4, 7, 10], dtype=np.float32))
    dense = ad.AnnData(X=np.arange(6, dtype=float).reshape(3, 2), var=pd.DataFrame(index=['A', 'B']))
    assert np.array_equal(main._gene_vector(dense, 'B'), np.array([1., 3., 5.]))


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='tier3-cache-', dir=ROOT) as cache, \
         patch.dict(os.environ, {'MPLCONFIGDIR': cache, 'HF_DATA_REPO': ''}):
        for name, test in list(globals().items()):
            if name.startswith('test_') and callable(test):
                test()
                print(f'OK: {name}')
    print('All Tier 3 checks passed.')
