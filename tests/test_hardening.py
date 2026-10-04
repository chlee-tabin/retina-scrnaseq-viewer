"""Public-release regressions. Synthetic data only; run: python tests/test_hardening.py."""
import inspect
import os
from pathlib import Path
import sys
import threading
import tempfile
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from utils import data_loading, data_provision, deg, plotting
from utils.smoothing import _coerce_sigma, resolve_smooth_sigma
from utils.state import encode_state, parse_url_state
from utils.validation import (
    NUMERIC_CONTROLS, category_allowed, coerce_camera, coerce_control, coerce_roi,
    obs_column_types,
)
from callbacks import main_callbacks as main, deg_callbacks as de_cb, url_callbacks


def refuses(call, message=None):
    try:
        call()
    except (deg.DEGError, ValueError) as e:
        if message:
            assert message in str(e), str(e)
    else:
        assert False, 'expected refusal'


def test_controls_and_state():
    for name, (lo, hi, default, integer) in NUMERIC_CONTROLS.items():
        for bad in (None, '', 'garbage', np.nan, np.inf, -np.inf, 'NaN', [], {}, True):
            assert coerce_control(name, bad) == default, (name, bad)
        assert coerce_control(name, -10) == lo
        assert coerce_control(name, str(default)) == default
        if hi is not None:
            assert coerce_control(name, 10**1000) == hi
            assert coerce_control(name, hi + 100) == hi
        if integer:
            assert isinstance(coerce_control(name, default + 0.9), int)
    assert coerce_control('bins', '60.9') == 60
    assert coerce_control('percentile', -1) == 0
    assert coerce_control('deg_min_cells', 0) == 50
    assert coerce_control('deg_min_cells', '0.0') == 50
    assert coerce_control('deg_min_cells', 5) == 5
    assert coerce_control('deg_min_cells', 100) == 100
    for bad in (np.nan, np.inf, 'bad', None):
        assert _coerce_sigma(bad, 1.5) == 1.5
    assert _coerce_sigma(1000, 2) == 4
    assert _coerce_sigma(10**1000, 2) == 4
    assert _coerce_sigma(-1, 2) == 0

    empty = {'A': [], 'B': []}
    for bad in (None, [], 'bad', {'C': []}, {'A': None}, {'A': [[1]]},
                {'A': [[1, 2, 3]]}, {'A': [['1', 2]]}, {'A': [[True, 1]]},
                {'B': [[np.nan, 2]]}, {'A': [[np.inf, 2]]},
                {'A': [[10**1000, 2]]}, {'A': [[1, 2]] * 501}):
        assert coerce_roi(bad) == empty, bad
        assert deg.polygon_to_indices(bad, [0], [0]) == []
    assert len(coerce_roi({'A': [[1, 2]] * 500})['A']) == 500
    assert coerce_roi({'A': [[1, 2], [2, 3], [4, 5]]})['B'] == []

    assert coerce_camera({'eye': {'x': 100, 'y': -100, 'z': '1'}}) == {
        'eye': {'x': 10., 'y': -10., 'z': 1.}}
    for bad in ([], {'eye': []}, {'eye': {'x': 1}}, {'eye': {'x': np.nan, 'y': 0, 'z': 1}},
                {'other': {}}, {'projection': {'type': 'invalid'}}):
        assert coerce_camera(bad) is None

    state = {'dataset': 'chick_rpc', 'bins': 999999, 'percentile': 'NaN',
             'roi': {'A': [[1, 2]] * 501}, 'wm_cuts': 999}
    search = '?state=' + encode_state(state)
    clean = parse_url_state(search)
    assert clean['bins'] == 120 and clean['percentile'] == .95 and clean['wm_cuts'] == 8
    assert clean['roi'] == empty
    assert parse_url_state('?state=' + encode_state(['bad'])) is None
    bad_sigma = parse_url_state('?state=' + encode_state({'v': 3, 'smooth_sigma': 'NaN'}))
    assert resolve_smooth_sigma(bad_sigma, 'chick_full', 1.5) == 1.5
    with patch.object(url_callbacks, 'ctx', SimpleNamespace(triggered_id='url.search')):
        outputs = url_callbacks.initialize_from_url(search, '/')
    assert outputs[2] == 120 and outputs[3] == .95 and outputs[30] == empty

    params = main._wholemount_params('NaN', 9999, -1, 9999, 10**1000,
                                    pow_p=99, gap_frac=-1)
    assert params['rho_max_nt_deg'] == 82 and params['rho_max_dv_deg'] == 150
    assert len(params['cut_angles_deg']) == 8 and params['pow_p'] == 3
    assert params['gap_frac'] == 0 and params['gap_gain'] == 0
    ng = main._nasal_gap_params({'layers': 2}, 9999, 'bad', -1, 999)
    assert ng == {'layers': 5, 'frac': .13, 'dorsal_deg': 0, 'ventral_deg': 90}


def grid():
    return np.array([0., 0., 1., 1.]), np.array([0., 1., 0., 1.])


def test_hover_orientation():
    x, y = grid()
    df = pd.DataFrame({'x': x, 'y': y, 'color': [1., 2., 3., 4.]})
    fig = plotting.create_binned_plot(df, 'custom_embedding', 'gene_expression',
                                     bin_size=2, percentile=1)
    trace = fig.data[0]
    assert np.array_equal(trace.z, [[1, 3], [2, 4]])
    for j in range(2):
        for i in range(2):
            hover = trace.customdata[j][i]
            assert f'Bin: ({i},{j})' in hover
            assert f'x: {trace.x[i]:.2f}' in hover and f'y: {trace.y[j]:.2f}' in hover
            assert f'Value: {trace.z[j][i]:.2f}' in hover
    # Category facet: asymmetric cell counts, rather than expression weights.
    df = pd.DataFrame({'x': np.repeat(x, [1, 2, 3, 4]),
                       'y': np.repeat(y, [1, 2, 3, 4]), 'color': 'group'})
    trace = plotting.create_binned_plot(df, 'custom_embedding', 'group', bin_size=2,
                                       treat_as_categorical=True).data[0]
    for j in range(2):
        for i in range(2):
            assert f'Bin: ({i},{j})' in trace.customdata[j][i]
            assert f'Cells: {int(trace.z[j][i])}' in trace.customdata[j][i]


def test_colour_only_clipping_and_signed_metadata():
    x, y = grid()
    vals = [np.array([1., 1., 1., 10.]), np.array([10., 10., 10., 10.])]
    h, _, _, _, vmax = plotting._binned_mean(x, y, vals[0], bin_size=2, percentile=.5)
    assert vmax == 1 and h[1, 1] == 10
    fig = plotting.create_binned_plot(pd.DataFrame({'x': x, 'y': y, 'color': vals[0]}),
                                     'custom_embedding', 'gene_expression', bin_size=2,
                                     percentile=.5)
    assert fig.layout.coloraxis.cmax == 1 and fig.data[0].z[1][1] == 10
    assert 'Value: 10.00' in fig.data[0].customdata[1][1]
    for shared in (True, False):
        fig = plotting.create_dual_gene_figure(x, y, vals, ['g1', 'g2'], 'custom_embedding',
                                              binned=True, bin_size=2, percentile=.5,
                                              shared_scale=shared)
        a, b = fig.data
        assert a.z[1][1] == b.z[1][1] == 10  # native Heatmap hover reads the unclipped z
        assert a.zmax == (10 if shared else 1) and b.zmax == 10

    # Flower and sphere also use unclipped grid data with one common colour maximum.
    fig = plotting.create_dual_gene_sphere_figure(y, x, vals, ['g1', 'g2'], binned=True,
                                                 bin_size=2, percentile=.5, shared_scale=True)
    meshes = [t for t in fig.data if t.type == 'mesh3d']
    assert len(meshes) == 2 and all(t.cmax == 10 for t in meshes)
    assert max(meshes[0].intensity) == 10
    assert any('10.00<br>' in str(t.text) for t in fig.data if t.type == 'scatter3d')
    fig = plotting.create_dual_gene_wholemount_figure(y, x, vals, ['g1', 'g2'],
                                                     bin_size=2, percentile=.5, shared_scale=True,
                                                     params={'cut_angles_deg': ()})
    bars = [t for t in fig.data if t.type == 'scatter' and t.marker.cmax is not None]
    assert len(bars) == 2 and all(t.marker.cmax == 10 for t in bars)
    assert any('10.00<br>' in str(t.text) for t in fig.data)

    signed = np.array([-4., -2., 1., 3.])
    fig = plotting.create_binned_plot(pd.DataFrame({'x': x, 'y': y, 'color': signed}),
                                     'custom_embedding', 'DV.Score', bin_size=2)
    assert fig.layout.coloraxis.cmin < 0
    assert fig.layout.coloraxis.cmin == -fig.layout.coloraxis.cmax
    sphere = plotting.create_sphere_binned_figure(y, x, signed, bin_size=2)
    mesh = next(t for t in sphere.data if t.type == 'mesh3d')
    assert mesh.cmin < 0 and mesh.cmin == -mesh.cmax
    flower = plotting.create_wholemount_binned_figure(y, x, signed, bin_size=2,
                                                     params={'cut_angles_deg': ()})
    assert any(t.marker.cmin is not None and t.marker.cmin < 0
               and t.marker.cmin == -t.marker.cmax for t in flower.data if t.type == 'scatter')
    all_negative = plotting.create_binned_plot(
        pd.DataFrame({'x': x, 'y': y, 'color': -np.arange(1., 5.)}),
        'custom_embedding', 'DV.Score', bin_size=2)
    assert all_negative.layout.coloraxis.cmin < 0 < all_negative.layout.coloraxis.cmax


def count_object(counts, sparse=False):
    x = sp.csr_matrix(counts) if sparse else np.asarray(counts, dtype=float)
    return ad.AnnData(x, obs=pd.DataFrame(index=[f'c{i}' for i in range(len(counts))]),
                      var=pd.DataFrame(index=[f'g{i}' for i in range(len(counts[0]))]))


def test_count_sources():
    counts = np.array([[1, 2, 0], [2, 0, 3], [0, 4, 2]])
    for sparse in (False, True):
        a = count_object(counts, sparse)
        a.raw = a
        r, genes = deg.get_raw_counts(a)
        assert np.array_equal(r.toarray() if sp.issparse(r) else r, counts)
        assert list(genes) == ['g0', 'g1', 'g2']
        logged = count_object(np.log1p(counts), sparse)
        logged.raw = logged
        refuses(lambda: deg.get_raw_counts(logged), 'raw counts unavailable')
        no_counts = count_object(counts, sparse)  # even integer .X is not a verified source
        refuses(lambda: deg.get_raw_counts(no_counts), 'raw counts unavailable')
        depth = counts.sum(axis=1)
        cp = counts / depth[:, None] * 1e4
        recon = count_object(np.log1p(cp), sparse)
        recon.obs['nCount_RNA'] = depth
        r, _ = deg.get_raw_counts(recon, norm_target=1e4)
        assert np.allclose(r.toarray() if sp.issparse(r) else r, counts)
        recon.raw = recon  # a logged .raw does not prevent a verified .X reconstruction
        assert np.allclose(deg.get_raw_counts(recon, 1e4)[0].toarray()
                           if sparse else deg.get_raw_counts(recon, 1e4)[0], counts)
        recon.obs['nCount_RNA'] = depth * 1.123
        refuses(lambda: deg.get_raw_counts(recon, 1e4), 'raw counts unavailable')
        recon.obs['nCount_RNA'] = -depth
        refuses(lambda: deg.get_raw_counts(recon, 1e4), 'raw counts unavailable')
    a = count_object([[-1., 2.]])
    a.raw = a
    refuses(lambda: deg.get_raw_counts(a), 'raw counts unavailable')


def test_missing_replicates():
    # Five cells per real pseudobulk; malformed labels must not make extra samples.
    obs = pd.DataFrame({'library': ['L1'] * 10 + ['L2'] * 10 + ['L1', np.nan, None, ''],
                        'genotype': ['d0'] * 20 + [np.nan, 'd0', 'd0', 'd0']})
    side = np.array(['A'] * 5 + ['B'] * 5 + ['A'] * 5 + ['B'] * 5 + ['A'] * 4)
    counts = np.ones((24, 2), dtype=int)
    a = ad.AnnData(sp.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=['g0', 'g1']))
    a.raw = a
    reps = deg.replicate_labels(obs, ['library', 'genotype'])
    assert pd.isna(reps[-4:]).all()
    assert pd.isna(deg.replicate_labels(obs.iloc[-4:], ['library', 'genotype'])).all()
    group_reps = main._replicate_series(a, None, ['library', 'genotype'])
    assert pd.isna(group_reps[-4:]).all()
    pb, meta = deg.aggregate_pseudobulk(counts, ['g0', 'g1'], reps, side, min_cells=5)
    assert len(meta) == 4 and meta.n_cells.sum() == 20 and pb.to_numpy().sum() == 40
    with patch.object(deg, 'run_deg', return_value=(pd.DataFrame(), {})):
        _, info = deg.deg_from_labels(a, ['library', 'genotype'], side, min_cells=5)
    assert info['n_excluded_replicate'] == 4
    fig = plotting.create_group_expression_figure(
        pd.DataFrame({'expr': np.ones(24), 'pop': 'RPC', 'replicate': group_reps}),
        'DV.Score', 'pop', value_is_gene=False, min_cells_pseudobulk=5)
    assert any('4 cells excluded' in ann.text for ann in fig.layout.annotations)
    assert all('nan' not in str(t.customdata) for t in fig.data if t.customdata is not None)
    status = de_cb._status(dict(info, mode='A_vs_B', n_A=14, n_B=10, n_pb_A=2,
                               n_pb_B=2, design='~condition', design_reason='test', n_genes_tested=2), 5)
    assert '4 cells excluded' in str(status.children)


def test_design_and_no_refit():
    meta = pd.DataFrame({'side': ['A', 'A', 'B', 'B'],
                         'library': ['L1', 'L2', 'L3', 'L4']})
    design, reason = deg.choose_design(meta)
    assert design == '~condition'
    assert reason == 'library term not estimable: each library is on one side only'
    pb = pd.DataFrame(np.ones((4, 2), dtype=int), columns=['g0', 'g1'])
    designs = []

    class FakeDDS:
        def __init__(self, **kwargs):
            designs.append(kwargs['design'])

        def deseq2(self):
            pass

    class FakeStats:
        def __init__(self, *args, **kwargs):
            self.results_df = pd.DataFrame({'baseMean': [1, 1], 'log2FoldChange': [1, -1],
                                             'pvalue': [.01, .02], 'padj': [.02, .04]},
                                            index=['g0', 'g1'])

        def summary(self):
            pass

    with patch.dict(sys.modules, {
        'pydeseq2.dds': SimpleNamespace(DeseqDataSet=FakeDDS),
        'pydeseq2.ds': SimpleNamespace(DeseqStats=FakeStats),
    }):
        _, info = deg.run_deg(pb, meta, min_cells=5)
        assert designs == ['~condition'] and info['design_reason'] == reason
        paired = meta.copy()
        paired['library'] = ['L1', 'L2', 'L1', 'L2']
        assert deg.choose_design(paired)[0] == '~library + condition'
        # No residual df for the library term -> explicit ~condition, decided before fitting.
        thin = pd.DataFrame({'side': ['A', 'B', 'B'], 'library': ['L1', 'L1', 'L2']})
        assert deg.choose_design(thin) == (
            '~condition', 'library term not estimable: too few pseudobulks to adjust for library')
        with patch.object(FakeDDS, 'deseq2', side_effect=RuntimeError('/secret/path traceback')), \
             patch.object(deg.logger, 'exception') as log:
            refuses(lambda: deg.run_deg(pb, paired, min_cells=5), 'fit failed')
            assert log.called
        assert designs == ['~condition', '~library + condition']  # no retry


def test_de_semaphore():
    entered, release = threading.Event(), threading.Event()

    def fit(*args):
        entered.set()
        assert release.wait(5)
        return 'done'

    with patch.object(deg, '_run_deg', side_effect=fit), ThreadPoolExecutor(max_workers=1) as pool:
        running = pool.submit(deg.run_deg, None, None)
        try:
            assert entered.wait(5)
            refuses(lambda: deg.run_deg(None, None), deg.DE_BUSY_MESSAGE)
        finally:
            release.set()
        assert running.result(timeout=5) == 'done'
    with patch.object(deg, '_run_deg', return_value='retry'):
        assert deg.run_deg(None, None) == 'retry'


def test_load_single_flight_and_identity():
    paths = ['data/synthetic1.h5ad', 'data/synthetic2.h5ad', 'data/synthetic3.h5ad']
    cfg = {'datasets': {str(i): {'file_path': p} for i, p in enumerate(paths)}}
    data_loading.load_adata.cache_clear()
    start = threading.Barrier(8)
    entered, release = threading.Event(), threading.Event()
    fake = SimpleNamespace(n_obs=1)

    def read(path):
        entered.set()
        assert release.wait(5)
        return fake

    def load():
        start.wait(timeout=5)
        return data_loading.load_adata(paths[0])

    with patch.object(data_loading, 'load_dataset_config', return_value=cfg), \
         patch.object(data_provision, 'ensure_one_dataset') as provision, \
         patch.object(data_loading.ad, 'read_h5ad', side_effect=read) as reader, \
         ThreadPoolExecutor(max_workers=8) as pool:
        jobs = [pool.submit(load) for _ in range(8)]
        try:
            assert entered.wait(5)
        finally:
            release.set()
        assert all(f.result(timeout=5) is fake for f in jobs)
        assert reader.call_count == provision.call_count == 1
        assert data_loading.load_adata.cache_info().maxsize == 2
        refuses(lambda: data_loading.load_adata('/secret/unconfigured.h5ad'))
        assert reader.call_count == 1
        assert data_loading.load_dataset_state({'dataset_id': 'unknown', 'filename': paths[0]}) == (None, {})
        for path in paths[1:]:
            data_loading.load_adata(path)
        data_loading.load_adata(paths[0])  # LRU eviction still works
        assert reader.call_count == 4 and data_loading.load_adata.cache_info().currsize == 2
    data_loading.load_adata.cache_clear()
    with patch.object(data_loading, 'load_dataset_config', return_value=cfg), \
         patch.object(data_provision.os.path, 'exists') as exists:
        refuses(lambda: data_provision.ensure_one_dataset('unconfigured.h5ad'))
        assert not exists.called  # refusal precedes any file access


def test_callbacks_and_download():
    # Import layouts/callbacks without prewarming or opening dataset files.
    with patch.dict(os.environ, {'HF_DATA_REPO': ''}):
        import app

    x = np.tile([0., 1., 2.], (201, 1))
    obs = pd.DataFrame({'barcode': [f'cell{i}' for i in range(201)],
                        'group': ['RPC'] * 100 + ['Other'] * 101,
                        'DV.Score': np.linspace(-1, 1, 201), 'NT.Score': np.linspace(-1, 1, 201),
                        'library': 'L1'})
    a = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=['FGF8', 'g1', 'g2']))
    a.raw = a
    assert not category_allowed(obs.barcode)
    assert 'barcode' not in obs_column_types(obs)
    forged = {'dataset_id': 'chick_rpc', 'filename': '/secret/file.h5ad',
              'replicate_columns': ['barcode'], 'x_norm_target': 999, 'min_cells_per_bin': 10**100,
              'figure_regions': [{'name': 'forged', 'n': 0}], 'column_types': {'barcode': 'categorical'}}
    with patch.object(data_loading, 'load_adata', return_value=a) as reader, \
         patch.object(data_loading, 'dataset_norm_target', return_value=None), \
         patch.object(app, 'load_adata', return_value=a), \
         patch.object(app, 'dataset_norm_target', return_value=None):
        loaded, settings = data_loading.load_dataset_state(forged)
        assert loaded is a and settings['replicate_columns'] == ['library', 'genotype']
        assert settings['min_cells_per_bin'] == 5 and settings['x_norm_target'] is None
        assert all(call.args == ('data/20250604_chick_RPC.h5ad',) for call in reader.call_args_list)
        out = app.update_data('chick_rpc', '?state=' + encode_state({'dataset': 'chick_rpc',
                                                                    'color': 'barcode'}), None, None)
        assert out[3] == 'gene_expression'
        assert 'barcode' not in [o['value'] for o in out[2]]
        controls = main.populate_group_controls(forged, '?state=' + encode_state({
            'dataset': 'chick_rpc', 'group_by': 'barcode', 'group_split': 'barcode'}), 'gene')
        assert 'barcode' not in [o['value'] for o in controls[0]] and controls[1] != 'barcode'
        assert controls[3] == ''

        # Supply every callback argument, then forge compute-side bounds independently of URL.
        args = {name: None for name in inspect.signature(main.update_plot).parameters}
        args.update(data_store=forged, embedding='custom_embedding', custom_x='NT.Score',
                    custom_y='DV.Score', color_by='gene_expression', gene='FGF8',
                    bin_number=10**1000, percentile='NaN', smooth_sigma_ctrl=np.inf,
                    enable_binning=['enabled'], enable_smoothing=['enabled'], plot_type='map')
        with patch.object(main, 'create_binned_plot', return_value=plotting.go.Figure()) as plot:
            main.update_plot(**args)
            assert plot.call_args.kwargs['bin_size'] == 120
            assert plot.call_args.kwargs['percentile'] == .95
            assert plot.call_args.kwargs['smooth_sigma'] == 2
            assert plot.call_args.kwargs['min_cells'] == 5
        args.update(color_by='barcode', bin_number=2)
        with patch.object(main, 'create_binned_plot', return_value=plotting.go.Figure()) as plot:
            main.update_plot(**args)
            assert plot.call_args.args[2] == 'gene_expression'
        args.update(plot_type='group', group_by='barcode', group_gene='g1', group_style='box')
        with patch.object(main, 'create_group_expression_plot', return_value=plotting.go.Figure()) as plot:
            main.update_plot(**args)
            assert plot.call_args.args[2] != 'barcode'
        with patch.object(de_cb, 'deg_from_labels', side_effect=deg.DEGError('safe refusal')) as run:
            de_cb.run_deg_cb(1, {'A': [[-2, -2], [2, -2], [2, 2], [-2, 2]]}, forged,
                            '100', 'NT.Score', 'DV.Score', 'custom_embedding', 'raw',
                            'gene_expression', [])
            assert run.call_args.args[1] == ['library', 'genotype']
            assert run.call_args.kwargs['min_cells'] == 100 and run.call_args.kwargs['norm_target'] is None
        assert 'forged' not in [o['value'] for o in de_cb.populate_roi_presets(forged)]
        reader.reset_mock()
        assert data_loading.load_dataset_state({'dataset_id': 'unknown'}) == (None, {})
        assert not reader.called

    # Flask's in-process test client does not start a server or read any .h5ad files.
    client = app.server.test_client()
    with patch.object(app.os.path, 'isfile', return_value=True) as isfile, \
         patch.object(app, 'send_file', return_value='configured download') as send:
        response = client.get('/download/unconfigured.h5ad')
        assert response.status_code == 403 and not isfile.called and not send.called
        assert client.get('/download/subdir/20250604_chick_RPC.h5ad').status_code == 403
        assert client.get('/download/../../secret').status_code == 403
        assert client.get('/download/20250604_chick_RPC.h5ad').status_code == 200
        assert send.call_args.kwargs['download_name'] == '20250604_chick_RPC.h5ad'
    with patch.object(app, 'load_adata', side_effect=RuntimeError('/secret/file.h5ad traceback')), \
         patch.object(app.logger, 'exception') as log:
        message = app.update_data('chick_rpc', '', None, None)[6]
        assert '/secret' not in message and 'traceback' not in message
        assert log.called


def test_replicate_unit_survives_group_cap_and_signature_needs_no_load():
    # 101 libraries x 2 genotypes = 202 composite levels (> 200 cap): still the replicate
    # unit, but not offered as a group-by axis.
    n = 404
    obs = pd.DataFrame({'library': pd.Categorical([f'L{i % 101}' for i in range(n)]),
                        'genotype': pd.Categorical([f'g{(i // 101) % 2}' for i in range(n)])})
    adata = SimpleNamespace(obs=obs)
    settings = {'column_types': {'library': 'categorical', 'genotype': 'categorical'},
                'replicate_columns': ['library', 'genotype']}
    _, reps, composite, group_composite, standalone, _ = main._group_fields(adata, settings)
    assert reps == ['library', 'genotype'] and composite is not None
    assert group_composite is None and composite not in standalone
    # Run-signature / view guards resolve column types from config + cache, never AnnData.
    with patch.object(data_loading, 'load_adata', side_effect=AssertionError('reloaded')), \
         patch.object(de_cb, 'dataset_column_types', return_value={'cluster': 'categorical'}):
        sig = de_cb._run_signature({'A': []}, {'dataset_id': 'chick_rpc'}, 50, 'NT.Score',
                                   'DV.Score', 'custom_embedding', 'raw',
                                   'cluster', ['enabled'])
        assert isinstance(sig, str)


if __name__ == '__main__':
    # Keep Matplotlib's cache (polygon selection) inside the worktree and clean it up.
    with tempfile.TemporaryDirectory(prefix='hardening-cache-', dir=Path(__file__).resolve().parents[1]) as cache, \
         patch.dict(os.environ, {'MPLCONFIGDIR': cache, 'HF_DATA_REPO': ''}):
        for test in (test_controls_and_state, test_hover_orientation,
                     test_colour_only_clipping_and_signed_metadata, test_count_sources,
                     test_missing_replicates, test_design_and_no_refit, test_de_semaphore,
                     test_load_single_flight_and_identity, test_callbacks_and_download,
                     test_replicate_unit_survives_group_cap_and_signature_needs_no_load):
            test()
            print(f'OK: {test.__name__}')
    print('All hardening checks passed.')
