"""Callbacks for the ROI differential-expression view.

Two topographic panels capture ROI A / ROI B by lasso; "Run differential expression"
resolves the (possibly overlapping) ROIs into disjoint sides, builds library x genotype
pseudobulks, and runs a negative-binomial pseudobulk DE (pydeseq2). Results drive the
volcano, a sortable table, the resolved-selection recap, and a CSV download.

ROI panels are rendered with continuous colour and identity plot order, so px.scatter
emits a SINGLE trace whose Plotly pointIndex equals the original adata cell index --
do not switch them to categorical colour or a reordering plot_order, or the captured
indices would no longer map to the right cells (see utils.plotting.create_scatter_plot,
which reorders rows for the map view).
"""
import logging

import numpy as np
import pandas as pd
import scipy.sparse
from dash import Input, Output, State, callback, ctx, dcc, html, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from utils.data_loading import load_adata
from utils.plotting import _message_figure, create_scatter_plot
from utils.deg import (
    resolve_rois, deg_from_labels,
    DEGError, NoReplicate, InsufficientReplicates,
)
from utils.deg_plots import create_volcano_figure, create_recap_figure

logger = logging.getLogger(__name__)

_SCORE_COLS = ('NT.Score', 'DV.Score')


def _gene_vector(adata, gene):
    sub = adata[:, gene].X
    if scipy.sparse.issparse(sub):
        return sub.toarray().flatten()
    return np.asarray(sub).flatten()


def _has_scores(adata):
    return all(c in adata.obs.columns for c in _SCORE_COLS)


# ---- Toggle the main panel: map/group view vs. the ROI-DEG view ----
@callback(
    [Output('map-view-container', 'style'),
     Output('deg-view-container', 'style')],
    Input('plot-type', 'value'),
)
def toggle_deg_view(plot_type):
    if plot_type == 'deg':
        return {'display': 'none'}, {'display': 'block'}
    return {'display': 'block'}, {'display': 'none'}


# ---- Render the two (identical) topographic ROI panels ----
# Coloured by the dataset's default gene for spatial context; NOT re-rendered on every
# gene keystroke, so a lasso selection survives until the user changes it or the dataset.
@callback(
    [Output('roi-plot-a', 'figure'),
     Output('roi-plot-b', 'figure')],
    [Input('data-store', 'data'),
     Input('plot-type', 'value')],
    prevent_initial_call=True,
)
def render_roi_panels(data_store, plot_type):
    if plot_type != 'deg' or not data_store:
        return no_update, no_update
    try:
        adata = load_adata(data_store['filename'])
        if not _has_scores(adata):
            msg = _message_figure(
                "ROI DEG needs DV.Score / NT.Score (the topographic RPC datasets). "
                "This dataset has no spatial scores.")
            return msg, msg

        x = adata.obs['NT.Score'].to_numpy()
        y = adata.obs['DV.Score'].to_numpy()
        gene = data_store.get('default_gene')
        if gene and gene in adata.var_names:
            color, label = _gene_vector(adata, gene), gene
        else:
            color, label = np.zeros(len(x)), ''
        df = pd.DataFrame({'x': x, 'y': y, 'color': color})
        # plot_order='none' -> identity order -> pointIndex == adata cell index (module docstring).
        fig = create_scatter_plot(df, 'NT.Score vs DV.Score', label, plot_order='none')
        fig.update_layout(dragmode='lasso', margin=dict(t=20))
        return fig, fig
    except Exception as e:
        logger.error(f"render_roi_panels failed: {e}", exc_info=True)
        msg = _message_figure("Could not render the ROI panels — see server logs.")
        return msg, msg


# ---- Capture each panel's lasso into the ROI store ----
@callback(
    Output('deg-roi-store', 'data'),
    [Input('roi-plot-a', 'selectedData'),
     Input('roi-plot-b', 'selectedData'),
     Input('clear-roi-btn', 'n_clicks')],
    State('deg-roi-store', 'data'),
    prevent_initial_call=True,
)
def capture_rois(sel_a, sel_b, clear_clicks, store):
    store = store or {'A': [], 'B': []}
    trig = ctx.triggered_id
    if trig == 'clear-roi-btn':
        return {'A': [], 'B': []}

    def idxs(sel):
        return [int(p['pointIndex']) for p in sel.get('points', [])] if sel else []

    if trig == 'roi-plot-a':
        return {**store, 'A': idxs(sel_a)}
    if trig == 'roi-plot-b':
        return {**store, 'B': idxs(sel_b)}
    return no_update


@callback(
    [Output('roi-a-info', 'children'),
     Output('roi-b-info', 'children')],
    Input('deg-roi-store', 'data'),
)
def roi_info(store):
    store = store or {}
    a, b = len(store.get('A') or []), len(store.get('B') or [])
    return (f"{a:,} cells selected",
            f"{b:,} cells selected" if b else "empty → contrast vs. rest")


# ---- Run the pseudobulk DE for the current ROIs ----
@callback(
    [Output('deg-results-store', 'data'),
     Output('deg-recap-plot', 'figure'),
     Output('deg-status', 'children')],
    Input('run-deg-btn', 'n_clicks'),
    [State('deg-roi-store', 'data'),
     State('data-store', 'data'),
     State('deg-min-cells', 'value')],
    prevent_initial_call=True,
)
def run_deg_cb(n_clicks, roi, data_store, min_cells):
    if not n_clicks or not data_store:
        return no_update, no_update, no_update
    try:
        adata = load_adata(data_store['filename'])
        if not _has_scores(adata):
            return None, _message_figure(""), dbc.Alert(
                "This dataset has no DV/NT scores, so ROIs cannot be drawn.", color="warning")

        idx_a = (roi or {}).get('A') or []
        idx_b = (roi or {}).get('B') or []
        labels, sel_info = resolve_rois(idx_a, idx_b, adata.n_obs)
        if sel_info['mode'] is None:
            return None, _message_figure("Lasso at least one region to contrast."), ""

        recap = create_recap_figure(
            adata.obs['NT.Score'].to_numpy(), adata.obs['DV.Score'].to_numpy(),
            labels, sel_info['mode'])
        min_cells = int(min_cells or 50)

        try:
            res, info = deg_from_labels(
                adata, data_store.get('replicate_columns', []), labels,
                min_cells=min_cells, min_frac=0.5,
                norm_target=data_store.get('x_norm_target'))
        except NoReplicate:
            return None, recap, dbc.Alert(
                "This dataset has no configured replicate unit (e.g. library × genotype), "
                "so a pseudobulk test is not possible here.", color="warning")
        except InsufficientReplicates as e:
            return None, recap, dbc.Alert(
                f"Too few pseudobulk replicates ≥ {e.min_cells} cells per side "
                f"(A: {e.n_a}, B: {e.n_b}). Draw a larger ROI, or lower the min-cells floor.",
                color="warning")
        except DEGError as e:
            return None, recap, dbc.Alert(str(e), color="warning")

        info.update(sel_info)
        # NaN padj/LFC (low-count genes) is not valid JSON for a dcc.Store -> None.
        records = res.astype(object).where(pd.notnull(res), None).to_dict('records')
        return records, recap, _status(info, min_cells)
    except Exception as e:
        logger.error(f"run_deg_cb failed: {e}", exc_info=True)
        return (None, _message_figure("Internal error computing DEG — see server logs."),
                dbc.Alert("Internal error computing differential expression.", color="danger"))


def _status(info, min_cells):
    mode = "ROI A vs. rest" if info['mode'] == 'A_vs_rest' else "ROI A vs. ROI B"
    overlap = (f" · {info['n_overlap']:,} overlapping cells → smaller ROI"
               if info.get('n_overlap') else "")
    return dbc.Alert([
        html.B(f"{mode}. "),
        f"{info['n_A']:,} A / {info['n_B']:,} B cells{overlap}. ",
        f"Pseudobulks ≥ {min_cells} cells: {info['n_pb_A']} A / {info['n_pb_B']} B. ",
        f"Design {info['design']} ({info['design_reason']}). ",
        f"{info['n_genes_tested']:,} genes tested.",
    ], color="info", className="py-2 mb-1")


# ---- Render volcano + table from the stored results ----
@callback(
    [Output('volcano-plot', 'figure'),
     Output('deg-table', 'data'),
     Output('deg-table', 'columns')],
    Input('deg-results-store', 'data'),
    prevent_initial_call=True,
)
def render_results(records):
    if not records:
        return _message_figure("Run a contrast to see the volcano."), [], []
    res = pd.DataFrame(records)
    for c in ('baseMean', 'log2FoldChange', 'pvalue', 'padj'):
        if c in res:
            res[c] = pd.to_numeric(res[c], errors='coerce')  # None (from the store) -> NaN
    fig = create_volcano_figure(res)

    show = res.sort_values('padj', na_position='last').copy()
    for c in ('baseMean', 'log2FoldChange', 'pvalue', 'padj'):
        if c in show:
            show[c] = show[c].map(lambda v: f"{v:.3g}" if pd.notnull(v) else "")
    cols_order = [c for c in ('gene', 'log2FoldChange', 'padj', 'pvalue', 'baseMean')
                  if c in show.columns]
    columns = [{'name': c, 'id': c} for c in cols_order]
    return fig, show[cols_order].to_dict('records'), columns


# ---- CSV download of the full results table ----
@callback(
    Output('deg-download', 'data'),
    Input('deg-download-btn', 'n_clicks'),
    State('deg-results-store', 'data'),
    prevent_initial_call=True,
)
def download_deg(n_clicks, records):
    if not n_clicks or not records:
        raise PreventUpdate
    df = pd.DataFrame(records).sort_values('padj', na_position='last')
    return dcc.send_data_frame(df.to_csv, "roi_deg.csv", index=False)
