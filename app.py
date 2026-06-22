import dash
from dash import html, dcc, Input, Output, State, callback, ClientsideFunction
import dash_bootstrap_components as dbc
import plotly.express as px
import scanpy as sc
import numpy as np
import pandas as pd
from urllib.parse import urlparse, parse_qs, urlencode
import json
from functools import lru_cache
import logging
import base64
import io
import yaml
from datetime import datetime
from pathlib import Path
import time
import argparse
from flask import send_file
import os
import re
import scipy.sparse

# Import layouts
from layouts.sidebar import create_sidebar
from layouts.main import create_main_panel

# Import callbacks
from callbacks.dataset_callbacks import *
from callbacks.selection_callbacks import *
from callbacks.url_callbacks import *
from callbacks.main_callbacks import *
from callbacks.status_callbacks import *

# Import utilities
from utils.data_loading import load_adata, load_dataset_config, validate_datasets, choose_default_embedding
from utils.config import load_config
from utils.error_handling import handle_callback_error, log_callback_info
from components.status_bar import create_status_bar
from utils.plotting import create_scatter_plot, create_metacell_plot

# Add command line argument parsing. Use parse_known_args (not parse_args) so that
# importing this module under gunicorn -- where sys.argv carries gunicorn's own flags
# (--bind, --workers, ...) -- does not abort with "unrecognized arguments".
parser = argparse.ArgumentParser()
parser.add_argument('-debug', action='store_true', help='Enable debug logging')
args, _ = parser.parse_known_args()

# Configure logging based on command line argument
logging.basicConfig(
    level=logging.DEBUG if args.debug else logging.INFO,
    format='%(asctime)s %(levelname)s %(name)s %(message)s',
    handlers=[
        logging.FileHandler("app_debug.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Directory that holds the .h5ad datasets. Configurable so the data can be
# mounted/uploaded separately (e.g. on Hugging Face Spaces) rather than baked
# into the image. Dataset paths in datasets_config.yml are resolved relative to
# this directory (see utils/data_loading.resolve_data_path).
DATA_DIR = os.getenv("DATA_DIR", "data")

# On platforms where the data is not present locally (e.g. Hugging Face Spaces),
# download the configured datasets into DATA_DIR from a Hugging Face repo. This
# runs in a BACKGROUND thread so it never blocks the web server from binding its
# port -- the platform health check (and Space promotion) must not wait ~30-40 s
# for ~6 GB to download. Each dataset is also fetched on-demand in load_adata if a
# user selects it before the pre-warm reaches it. No-op unless HF_DATA_REPO is set.
import threading


def _prewarm_datasets():
    try:
        from utils.data_provision import ensure_datasets
        ensure_datasets(load_dataset_config())
    except Exception as e:
        logger.error(f"Background dataset provisioning failed (continuing): {e}")


threading.Thread(target=_prewarm_datasets, name="dataset-prewarm", daemon=True).start()

# Initialize the Dash app with bootstrap theme
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True
)
app.title = "Single-cell Data Viewer"  # Set the title for the browser tab
server = app.server

# Load configuration
config = load_config()

# Function to create 2D histogram (metacells) with improved binning and aggregation
def create_metacells(x, y, values=None, n_bins=50):
    H, xedges, yedges = np.histogram2d(x, y, bins=n_bins)
    
    if values is not None:
        H_values = np.zeros_like(H)
        H_counts = np.zeros_like(H)
        
        x_indices = np.digitize(x, xedges) - 1
        y_indices = np.digitize(y, yedges) - 1
        
        # Filter out points outside the bins
        mask = (x_indices >= 0) & (x_indices < H.shape[0]) & \
               (y_indices >= 0) & (y_indices < H.shape[1])
        
        x_indices = x_indices[mask]
        y_indices = y_indices[mask]
        values_masked = values[mask] if values is not None else None
        
        # Use numpy's add.at for efficient binning
        if values_masked is not None:
            np.add.at(H_values, (x_indices, y_indices), values_masked)
            np.add.at(H_counts, (x_indices, y_indices), 1)
            
        # Avoid division by zero
        mask = H_counts > 0
        H_values[mask] = H_values[mask] / H_counts[mask]
        
        return H_values, xedges, yedges
    return H, xedges, yedges


def _detect_cp_target(adata, n_sample=256):
    """Detect the counts-per-X normalisation target of a log1p-normalised ``.X``.

    For CP-normalised data every cell's ``expm1(.X)`` sums to the same target (1e4 for
    CP10K, 1e6 for CPM). Sample cells, sum ``expm1`` across genes, and return the
    median when the sums are tightly clustered; return None otherwise (raw counts,
    z-scored, scran-pooled, or otherwise not a clean log1p(CP)) so the pseudobulk panel
    is omitted rather than reconstructed from an unknown normalisation.
    """
    try:
        n = adata.n_obs
        if n == 0:
            return None
        idx = np.unique(np.linspace(0, n - 1, min(n_sample, n)).astype(int))
        X = adata.X[idx]
        if scipy.sparse.issparse(X):
            sums = np.asarray(np.expm1(X.toarray()).sum(axis=1)).ravel()
        else:
            sums = np.expm1(np.asarray(X)).sum(axis=1).ravel()
        sums = sums[np.isfinite(sums) & (sums > 0)]
        if sums.size < max(5, 0.5 * len(idx)):
            return None
        med = float(np.median(sums))
        # CP normalisation => every cell sums to the same target; require tight spread.
        if med <= 0 or (float(np.max(sums)) - float(np.min(sums))) / med > 0.02:
            return None
        return med
    except Exception as e:  # noqa: BLE001 - best-effort; omit Panel A on any failure
        logger.warning(f"CP-target detection failed ({e}); pseudobulk panel omitted")
        return None

# Main layout with fixed sidebar
app.layout = dbc.Container([
    dcc.Location(id='url', refresh=False),
    dcc.Store(id='data-store'),
    dcc.Store(id='selection-store'),
    dcc.Store(id='url-parameters'),
    dcc.Store(id='initial-load-flag', data=True),
    dcc.Store(id='session-id', data=str(time.time())),
    
    # Add interval component here
    dcc.Interval(
        id='interval-component',
        interval=5000,  # Update every 5 seconds
        n_intervals=0
    ),
    
    # Add status bar at the top
    create_status_bar(),
    
    dbc.Row([
        # Fixed sidebar
        dbc.Col(
            create_sidebar(),
            width=3,
            className="position-fixed",
            style={
                "height": "100vh",
                "overflowY": "auto"
            }
        ),
        # Main content with offset
        dbc.Col(
            create_main_panel(),
            width=9,
            className="offset-3"
        )
    ], className="g-0")  # g-0 removes gutters
], fluid=True)

# Callback to initialize dataset dropdown
@callback(
    Output('dataset-select', 'options'),
    Input('url', 'pathname')  # Remove search as input
)
def initialize_dataset_dropdown(pathname):
    try:
        config = load_dataset_config()
        datasets = validate_datasets(config)
        
        options = [
            {
                'label': data['title'],
                'value': dataset_id
            }
            for dataset_id, data in datasets.items()
        ]
        
        return options
    except Exception as e:
        logger.error(f"Error loading dataset configuration: {str(e)}")
        return []

# GEO accession links: turn "GSE#####" tokens into links to the GEO record page.
_GSE_RE = re.compile(r'(GSE\d+)')

def _linkify_gse(text):
    """Render a description string with any GSE##### accession as a GEO hyperlink."""
    parts = _GSE_RE.split(str(text))
    children = []
    for i, part in enumerate(parts):
        if i % 2 == 1:  # the captured GSE accession
            children.append(html.A(
                part,
                href=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={part}",
                target="_blank",
                rel="noopener noreferrer",
            ))
        elif part:
            children.append(part)
    return children

# Callback to display dataset information
@callback(
    Output('dataset-info', 'children'),
    [Input('dataset-select', 'value'),
     Input('data-store', 'data')]  # Add data-store as an input
)
def update_dataset_info(dataset_id, data_store):
    if not dataset_id:
        return ""
    
    try:
        config = load_dataset_config()
        dataset = config['datasets'][dataset_id]
        
        # Get number of cells from data_store if available
        n_cells = f"{data_store['n_cells']:,}" if data_store and 'n_cells' in data_store else 'N/A'

        # Download link. Use the basename so it matches the /download route,
        # which resolves the request relative to DATA_DIR (the config's leading
        # "data/" is not part of the on-disk filename under DATA_DIR).
        download_name = os.path.basename(dataset['file_path'])
        download_path = f"/download/{download_name}"

        return dbc.Card([
            dbc.CardBody([
                html.H5(dataset['title'], className='card-title'),
                html.P(_linkify_gse(dataset['description']), className='card-text'),
                html.P([
                    html.Strong("Last Updated: "),
                    dataset['last_updated']
                ], className='card-text'),
                html.P([
                    html.Strong("Number of Cells: "),
                    n_cells
                ], className='card-text'),
                html.P([
                    html.Strong("Download: "),
                    html.A(
                        "Download .h5ad file",
                        href=download_path,
                        download=download_name,
                        className="btn btn-outline-primary btn-sm"
                    )
                ], className='card-text')
            ])
        ])
    except Exception as e:
        logger.error(f"Error loading dataset info: {str(e)}")
        return html.Div(f"Error loading dataset information: {str(e)}")

# Simplify the data loading callback to only handle dataset selection
@callback(
    [Output('data-store', 'data', allow_duplicate=True),
     Output('embedding-select', 'options', allow_duplicate=True),
     Output('color-select', 'options', allow_duplicate=True),
     Output('color-select', 'value', allow_duplicate=True),
     Output('gene-select', 'value', allow_duplicate=True),
     Output('gene-select', 'options', allow_duplicate=True),
     Output('loading-output', 'children', allow_duplicate=True),
     Output('embedding-select', 'value', allow_duplicate=True)],
    [Input('dataset-select', 'value'),
     Input('url', 'search')],
    [State('color-select', 'value'),
     State('gene-select', 'value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_data(dataset_id, url_search, current_color, current_gene):
    if not dataset_id:
        return None, [], [], None, None, [], "", 'custom_embedding'
    
    try:
        config = load_dataset_config()
        dataset = config['datasets'][dataset_id]
        adata = load_adata(dataset['file_path'])
        
        # Determine column types. A low-cardinality INTEGER column (cluster ids, Phase
        # codes) is treated as CATEGORICAL -- matching _is_categorical_series in
        # main_callbacks -- so such an annotation is offered in the group-by / colour
        # controls instead of being mis-handled as a continuous axis.
        column_types = {}
        for col in adata.obs.columns:
            s = adata.obs[col]
            is_small_int = pd.api.types.is_integer_dtype(s) and s.nunique() <= 50
            if pd.api.types.is_numeric_dtype(s) and not is_small_int:
                column_types[col] = 'numeric'
            else:
                column_types[col] = 'categorical'
        
        # Create embedding list including custom embedding option
        available_embeddings = list(adata.obsm.keys())
        
        data_store = {
            'filename': dataset['file_path'],
            'dataset_id': dataset_id,
            'n_cells': adata.n_obs,
            'embeddings': available_embeddings + ['custom_embedding'],
            'metadata_cols': list(adata.obs.columns),
            'genes': list(adata.var_names),
            'column_types': column_types,  # Add column types to data store
            # Per-dataset spatial-binning params (mirror the analysis pipeline);
            # used by the binned/smoothed view in main_callbacks.update_plot.
            'smooth_sigma': dataset.get('smooth_sigma', 1.5),
            'min_cells_per_bin': dataset.get('min_cells_per_bin', 1),
            'color_floor': dataset.get('color_floor', 0.05),
            # Optional per-dataset gene sets for the "Expression by group" dot
            # plot, and a path to precomputed DEG results. Both default to
            # absent/empty so datasets without them simply hide the controls.
            'gene_modules': dataset.get('gene_modules', {}),
            'deg_results_path': dataset.get('deg_results_path'),
            # Default gene (for the group view), categorical display (consistent
            # per-type colours + biological order), and the pseudobulk replicate
            # unit -- all optional, all from datasets_config.yml.
            'default_gene': dataset.get('default_gene'),
            'annotation_column': dataset.get('annotation_column'),
            'annotation_order': dataset.get('annotation_order', []),
            'annotation_colors': dataset.get('annotation_colors', {}),
            'replicate_columns': dataset.get('replicate_columns', []),
            # CP-normalisation target of .X (1e4=CP10K, 1e6=CPM, ...) so the group
            # "Figure" reconstructs raw counts for ANY log1p(CP*) normalisation; None
            # when .X is not a clean log1p(CP) (then Panel A is omitted).
            'x_norm_target': _detect_cp_target(adata),
        }
        
        # Create embedding options with custom embedding as first option
        embedding_options = [
            {'label': 'Custom embedding by meta scores', 'value': 'custom_embedding'}
        ] + [
            {'label': emb, 'value': emb} for emb in available_embeddings
        ]
        
        # Check URL state if available. Apply it only when the loaded dataset matches
        # the one in the shared link, so a later manual dataset switch is not silently
        # re-overridden by a stale ?state= (which is never cleared from the URL).
        state = parse_url_state(url_search) if url_search else None
        if state and state.get('dataset') != dataset_id:
            state = None
        url_color = state.get('color') if state else None
        url_gene = state.get('gene') if state else None
        
        # Always include Gene Expression option
        color_options = [
            {'label': 'Gene Expression', 'value': 'gene_expression'}
        ]
        color_options.extend([
            {'label': col, 'value': col} 
            for col in data_store['metadata_cols']
        ])
        
        # Determine color and gene values. Defaults when there is no URL state
        # or prior selection: color by gene expression, and the dataset's
        # configured default_gene (species-correct casing, e.g. FGF8 / Fgf8).
        default_gene = dataset.get('default_gene')
        valid_color_values = ['gene_expression'] + data_store['metadata_cols']
        if url_color in valid_color_values:
            color_value = url_color
        elif current_color in valid_color_values:
            color_value = current_color
        else:
            color_value = 'gene_expression'

        if url_gene in data_store['genes']:
            gene_value = url_gene
        elif current_gene in data_store['genes']:
            gene_value = current_gene
        elif default_gene in data_store['genes']:
            gene_value = default_gene
        else:
            gene_value = None
        
        # Pre-populate gene options so the default gene's value is valid the
        # moment it is set (a Dropdown value absent from its options is dropped
        # client-side; gene options are otherwise filled by update_gene_select).
        gene_options = [{'label': g, 'value': g} for g in sorted(data_store['genes'])]

        # Pick the embedding to open on: a UMAP for datasets without DV/NT spatial
        # scores (e.g. the full-retina object) so they don't open on a blank custom
        # plot; the topographic custom view for the spatial RPC datasets.
        embedding_value = choose_default_embedding(
            available_embeddings, set(data_store['metadata_cols']),
            config_default=dataset.get('default_embedding'),
            url_embedding=(state.get('embedding') if state else None))

        return data_store, embedding_options, color_options, color_value, gene_value, gene_options, "", embedding_value
        
    except Exception as e:
        error_message = f"Error loading data: {str(e)}"
        logger.error(error_message)
        return None, [], [], None, None, [], error_message, 'custom_embedding'

# NOTE: restoring controls from a shared URL is handled solely by
# callbacks/url_callbacks.initialize_from_url (which sets the dataset + the global
# controls) together with update_data, update_custom_embedding_controls and
# populate_group_controls -- each reads the same URL state for the controls it owns.
# The earlier duplicate restore callback that lived here was removed so that no two
# callbacks write the same control values (which raced and dropped the new
# expression-by-group / figure / two-gene state).

# NOTE: cell-selection handling (selection-store / selection-info, driven by
# main-plot.selectedData) lives solely in callbacks/selection_callbacks.py. A
# duplicate copy here registered a second writer for the same two outputs on the
# same input; it was removed so there is one authoritative selection callback.

# NOTE: the main-plot figure callback lives in callbacks/main_callbacks.py
# (update_plot there handles custom embeddings, gene-expression coloring, plot
# ordering, and binning). The earlier copy here was an incomplete duplicate
# (it referenced custom_x/custom_y/gene that were not callback inputs) and has
# been removed so there is one authoritative plot callback.

# Add this after app initialization
@server.route('/download/<path:filepath>')
def download_file(filepath):
    try:
        # Restrict downloads to files that live inside DATA_DIR. Resolve both
        # the data directory and the requested file to absolute, symlink-free
        # paths and confirm the request stays within DATA_DIR. This rejects
        # path-traversal attempts (e.g. "../../etc/passwd", absolute paths, or
        # symlink escapes) instead of relying on os.path.basename alone.
        data_root = os.path.realpath(DATA_DIR)
        requested = os.path.realpath(os.path.join(data_root, filepath))

        # commonpath raises ValueError on mixed drives/relative inputs; treat
        # any such case, or a resolved path outside data_root, as forbidden.
        try:
            within_data_dir = os.path.commonpath([data_root, requested]) == data_root
        except ValueError:
            within_data_dir = False

        if not within_data_dir:
            logger.warning(f"Rejected download outside DATA_DIR: {filepath}")
            return "Forbidden", 403

        if os.path.isfile(requested):
            return send_file(
                requested,
                as_attachment=True,
                download_name=os.path.basename(requested)
            )
        else:
            return "File not found", 404
    except Exception as e:
        logger.error(f"Error serving download '{filepath}': {e}")
        return "Internal server error", 500

if __name__ == '__main__':
    # Host/port/debug are environment-driven so the same entrypoint works
    # locally and on Hugging Face Spaces (which expects 0.0.0.0:7860).
    debug = os.getenv("DASH_DEBUG", "false").lower() == "true"
    host = os.getenv("HOST", "0.0.0.0")
    port = int(os.getenv("PORT", "7860"))
    app.run_server(host=host, port=port, debug=debug)