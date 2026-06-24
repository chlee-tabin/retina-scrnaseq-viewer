"""Callbacks for the ROI differential-expression section.

The user draws one or two polygon ROIs directly on the spatial map (assets/roi_draw.js
captures clicks -> data coords, double-click closes, mirrors vertices into
roi-vertices-store). "Run" selects cells by point-in-polygon on the embedding's axes and
runs a negative-binomial pseudobulk DE (pydeseq2) over the demux replicate unit. Results
drive the volcano, a sortable table, a resolved-selection recap, and a CSV download.

Selection is point-in-polygon in score space, NOT Plotly point selection -- so the ROI
works the same on the per-cell scatter, the binned heatmap, or the smoothed map.
"""
import logging

import numpy as np
import pandas as pd
from dash import (Input, Output, State, callback, clientside_callback, ctx, dcc,
                  html, no_update)
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from utils.data_loading import load_adata
from utils.plotting import _message_figure
from utils.deg import (
    resolve_rois, polygon_to_indices, deg_from_labels,
    DEGError, NoReplicate, InsufficientReplicates,
)
from utils.deg_plots import create_volcano_figure, create_recap_figure

logger = logging.getLogger(__name__)


# ---- Show the ROI-DEG section only in the spatial-map view (not "Expression by group") ----
@callback(
    Output('deg-section', 'style'),
    Input('plot-type', 'value'),
)
def toggle_deg_section(plot_type):
    return {'display': 'block'} if plot_type == 'map' else {'display': 'none'}


# ---- Draw-mode buttons: pick which ROI to draw; Clear resets both ----
@callback(
    [Output('roi-draw-mode-store', 'data', allow_duplicate=True),
     Output('roi-vertices-store', 'data', allow_duplicate=True),
     Output('deg-results-store', 'data', allow_duplicate=True),
     Output('deg-status', 'children', allow_duplicate=True),
     Output('deg-recap-plot', 'figure', allow_duplicate=True)],
    [Input('roi-draw-a-btn', 'n_clicks'),
     Input('roi-draw-b-btn', 'n_clicks'),
     Input('roi-clear-btn', 'n_clicks')],
    prevent_initial_call=True,
)
def set_draw_mode(a_clicks, b_clicks, clear_clicks):
    trig = ctx.triggered_id
    if trig == 'roi-draw-a-btn':
        return 'A', no_update, no_update, no_update, no_update
    if trig == 'roi-draw-b-btn':
        return 'B', no_update, no_update, no_update, no_update
    # Clear: reset the ROIs AND invalidate any stale DEG result -- emptying
    # deg-results-store clears the volcano + table (via render_results) and disables the
    # CSV download; the status + recap are cleared here so nothing dangles for an empty ROI.
    return (None, {'A': [], 'B': []}, None, '',
            _message_figure("Cleared. Draw an ROI (or load a figure region) and Run."))


# ---- Invalidate stale results when the ROI changes (freehand draw, preset, URL restore) ----
# The vertices store is also written by the JS draw tool and the preset loader, neither of
# which goes through set_draw_mode -- so after a run, editing the polygon would leave the old
# volcano/table/CSV on screen describing a region the map no longer shows. Clear them here, the
# one place every vertex change routes through. Guarded so it only fires when there are results
# to clear and the ROI is non-empty (Clear is already handled above; a fresh URL restore has
# nothing stale to drop).
@callback(
    [Output('deg-results-store', 'data', allow_duplicate=True),
     Output('deg-status', 'children', allow_duplicate=True),
     Output('deg-recap-plot', 'figure', allow_duplicate=True)],
    Input('roi-vertices-store', 'data'),
    State('deg-results-store', 'data'),
    prevent_initial_call=True,
)
def invalidate_on_roi_change(verts, results):
    verts = verts or {}
    has_roi = any(verts.get(k) for k in ('A', 'B'))
    if not results or not has_roi:
        raise PreventUpdate
    return (None, '',
            _message_figure("ROI changed — click Run to recompute differential expression."))


# ---- Highlight the active draw button + report vertex counts ----
@callback(
    [Output('roi-draw-a-btn', 'outline'),
     Output('roi-draw-b-btn', 'outline'),
     Output('roi-draw-status', 'children')],
    [Input('roi-draw-mode-store', 'data'),
     Input('roi-vertices-store', 'data')],
)
def draw_state(mode, verts):
    verts = verts or {}
    na, nb = len(verts.get('A') or []), len(verts.get('B') or [])
    if mode:
        msg = (f"Drawing ROI {mode} — click to add vertices, double-click to close. "
               f"(A: {na} / B: {nb} vertices)")
    else:
        msg = f"A: {na} / B: {nb} vertices" if (na or nb) else "No ROI drawn yet."
    return mode != 'A', mode != 'B', msg


# ---- Figure-region presets: reproduce a manuscript area-DEG gate as ROI A ----
@callback(
    Output('roi-preset-select', 'options'),
    Input('data-store', 'data'),
)
def populate_roi_presets(data_store):
    regions = (data_store or {}).get('figure_regions') or []
    return [{'label': r['name'], 'value': r['name']} for r in regions]


@callback(
    [Output('roi-vertices-store', 'data', allow_duplicate=True),
     Output('roi-draw-mode-store', 'data', allow_duplicate=True),
     Output('custom-x-select', 'value', allow_duplicate=True),
     Output('custom-y-select', 'value', allow_duplicate=True)],
    Input('roi-preset-select', 'value'),
    State('data-store', 'data'),
    prevent_initial_call=True,
)
def load_preset_region(region_name, data_store):
    """Set ROI A to the rectangle of a manuscript area-DEG gate. The gate is an axis-aligned
    range in NT.Score × DV.Score fractions of the data's score range (matches the figure's
    `gate()` in area_significant_deg.R); reproduces a Fig 6H/S23 panel after Run. Also forces
    the custom axes back to NT.Score (x) / DV.Score (y) -- the gate is defined in that space,
    and run_deg_cb tests the polygon against the selected axes, so they must match."""
    nope = (no_update, no_update, no_update, no_update)
    if not region_name or not data_store:
        return nope
    spec = {r['name']: r for r in (data_store.get('figure_regions') or [])}.get(region_name)
    if not spec:
        return nope
    adata = load_adata(data_store['filename'])
    if not {'NT.Score', 'DV.Score'} <= set(adata.obs.columns):
        return nope
    nt = pd.to_numeric(adata.obs['NT.Score'], errors='coerce')
    dv = pd.to_numeric(adata.obs['DV.Score'], errors='coerce')
    nlo, nhi, dlo, dhi = float(nt.min()), float(nt.max()), float(dv.min()), float(dv.max())
    n, ni, di = spec['n'], spec['nt'], spec['dv']
    x0 = nlo + (nhi - nlo) * ni[0] / n
    x1 = nlo + (nhi - nlo) * (ni[1] + 1) / n
    y0 = dlo + (dhi - dlo) * di[0] / n
    y1 = dlo + (dhi - dlo) * (di[1] + 1) / n
    rect = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
    # ROI A = the gate; exit draw mode; pin axes to the gate's NT/DV space.
    return {'A': rect, 'B': []}, None, 'NT.Score', 'DV.Score'


# ---- Run the pseudobulk DE for the drawn ROIs ----
@callback(
    [Output('deg-results-store', 'data', allow_duplicate=True),
     Output('deg-recap-plot', 'figure', allow_duplicate=True),
     Output('deg-status', 'children', allow_duplicate=True)],
    Input('run-deg-btn', 'n_clicks'),
    [State('roi-vertices-store', 'data'),
     State('data-store', 'data'),
     State('deg-min-cells', 'value'),
     State('custom-x-select', 'value'),
     State('custom-y-select', 'value'),
     State('embedding-select', 'value'),
     State('custom-projection', 'value')],
    prevent_initial_call=True,
)
def run_deg_cb(n_clicks, verts, data_store, min_cells, custom_x, custom_y,
               embedding, projection):
    if not n_clicks or not data_store:
        return no_update, no_update, no_update
    try:
        # ROI selection is point-in-polygon on the raw DV/NT score axes, so it is only valid
        # on the custom embedding in raw mode. Reject anything else -- UMAP, and every
        # whole-mount projection (flower, sphere, ...), where the on-screen coords are warped
        # away from custom_x/custom_y. Gating on "not raw" (rather than naming each projection)
        # keeps any future projection rejected by default.
        if embedding != 'custom_embedding' or projection not in (None, 'raw'):
            return None, _message_figure(""), dbc.Alert(
                "Draw ROIs on the topographic DV/NT view (custom embedding, raw axes). "
                "UMAP and the whole-mount projections are not supported yet.", color="warning")
        if not custom_x or not custom_y:
            return None, _message_figure(""), dbc.Alert(
                "Select the X / Y axes for the custom embedding first.", color="warning")

        adata = load_adata(data_store['filename'])
        if custom_x not in adata.obs.columns or custom_y not in adata.obs.columns:
            return None, _message_figure(""), dbc.Alert(
                "The selected axes are not present in this dataset.", color="warning")

        xs = pd.to_numeric(adata.obs[custom_x], errors='coerce').to_numpy()
        ys = pd.to_numeric(adata.obs[custom_y], errors='coerce').to_numpy()
        verts = verts or {}
        idx_a = polygon_to_indices(verts.get('A'), xs, ys)
        idx_b = polygon_to_indices(verts.get('B'), xs, ys)
        labels, sel_info = resolve_rois(idx_a, idx_b, adata.n_obs)
        if sel_info['mode'] is None:
            return None, _message_figure(
                "Draw at least one ROI polygon (≥ 3 vertices) on the map."), ""

        recap = create_recap_figure(xs, ys, labels, sel_info['mode'],
                                    x_label=custom_x, y_label=custom_y)
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
    [Input('deg-results-store', 'data'),
     Input('volcano-lfc-thresh', 'value'),
     Input('volcano-padj-thresh', 'value')],
    prevent_initial_call=True,
)
def render_results(records, lfc_thresh, padj_thresh):
    if not records:
        return _message_figure("Draw an ROI and Run to see the volcano."), [], []
    res = pd.DataFrame(records)
    for c in ('baseMean', 'log2FoldChange', 'pvalue', 'padj'):
        if c in res:
            res[c] = pd.to_numeric(res[c], errors='coerce')  # None (from the store) -> NaN
    lfc_t = float(lfc_thresh) if lfc_thresh not in (None, '') else 1.0
    padj_t = float(padj_thresh) if padj_thresh not in (None, '') else 0.05
    fig = create_volcano_figure(res, lfc_thresh=lfc_t, padj_thresh=padj_t)

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


# ---- Clientside: attach/detach the drawing listeners + redraw on vertex changes ----
clientside_callback(
    "function(mode){ window.dash_clientside.roi.setMode(mode); return ''; }",
    Output('roi-mode-dummy', 'children'),
    Input('roi-draw-mode-store', 'data'),
)
clientside_callback(
    "function(verts){ window.dash_clientside.roi.syncVerts(verts); return ''; }",
    Output('roi-verts-dummy', 'children'),
    Input('roi-vertices-store', 'data'),
)
# Redraw the stored polygons whenever the map figure (re)renders -- so a shared-URL ROI
# draws once the plot mounts, and an ROI survives a view change (which replaces the figure
# and would otherwise wipe the overlay). Deferred a tick so the figure settles first.
clientside_callback(
    "function(_fig, verts){ if (verts) setTimeout(function(){ "
    "window.dash_clientside.roi.syncVerts(verts); }, 60); return ''; }",
    Output('roi-verts-dummy2', 'children'),
    Input('main-plot', 'figure'),
    State('roi-vertices-store', 'data'),
    prevent_initial_call=True,
)


# ---- Clientside activity bar: show on Run click, hide when the run finishes ----
clientside_callback(
    "function(n){ return [n ? {display:'block'} : {display:'none'}, !!n]; }",
    Output('deg-progress-wrap', 'style', allow_duplicate=True),
    Output('run-deg-btn', 'disabled', allow_duplicate=True),
    Input('run-deg-btn', 'n_clicks'),
    prevent_initial_call=True,
)
clientside_callback(
    "function(_fig, _data){ return [{display:'none'}, false]; }",
    Output('deg-progress-wrap', 'style', allow_duplicate=True),
    Output('run-deg-btn', 'disabled', allow_duplicate=True),
    Input('deg-recap-plot', 'figure'),
    Input('deg-results-store', 'data'),
    prevent_initial_call=True,
)
