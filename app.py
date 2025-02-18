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
from utils.data_loading import load_adata, load_dataset_config, validate_datasets
from utils.config import load_config
from utils.error_handling import handle_callback_error, log_callback_info
from components.status_bar import create_status_bar
from utils.plotting import create_scatter_plot, create_metacell_plot

# Add command line argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('-debug', action='store_true', help='Enable debug logging')
args = parser.parse_args()

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
        
        return dbc.Card([
            dbc.CardBody([
                html.H5(dataset['title'], className='card-title'),
                html.P(dataset['description'], className='card-text'),
                html.P([
                    html.Strong("Last Updated: "),
                    dataset['last_updated']
                ], className='card-text'),
                html.P([
                    html.Strong("Number of Cells: "),
                    n_cells
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
     Output('loading-output', 'children', allow_duplicate=True)],
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
        return None, [], [], None, None, ""
    
    try:
        config = load_dataset_config()
        dataset = config['datasets'][dataset_id]
        adata = load_adata(dataset['file_path'])
        
        # Determine column types
        column_types = {}
        for col in adata.obs.columns:
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                column_types[col] = 'numeric'
            else:
                column_types[col] = 'categorical'
        
        # Create embedding list including custom embedding option
        available_embeddings = list(adata.obsm.keys())
        
        data_store = {
            'filename': dataset['file_path'],
            'n_cells': adata.n_obs,
            'embeddings': available_embeddings + ['custom_embedding'],
            'metadata_cols': list(adata.obs.columns),
            'genes': list(adata.var_names),
            'column_types': column_types  # Add column types to data store
        }
        
        # Create embedding options with custom embedding as first option
        embedding_options = [
            {'label': 'Custom embedding by meta scores', 'value': 'custom_embedding'}
        ] + [
            {'label': emb, 'value': emb} for emb in available_embeddings
        ]
        
        # Check URL state if available
        state = parse_url_state(url_search) if url_search else None
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
        
        # Determine color and gene values
        valid_color_values = ['gene_expression'] + data_store['metadata_cols']
        color_value = (url_color if url_color in valid_color_values 
                      else current_color if current_color in valid_color_values 
                      else None)
        
        gene_value = (url_gene if url_gene in data_store['genes']
                     else current_gene if current_gene in data_store['genes']
                     else None)
        
        return data_store, embedding_options, color_options, color_value, gene_value, ""
        
    except Exception as e:
        error_message = f"Error loading data: {str(e)}"
        logger.error(error_message)
        return None, [], [], None, None, error_message

# Update the initialization callback to be more robust
@callback(
    [Output('dataset-select', 'value', allow_duplicate=True),
     Output('embedding-select', 'value', allow_duplicate=True),
     Output('color-select', 'value', allow_duplicate=True),
     Output('viz-mode', 'value', allow_duplicate=True),
     Output('gene-select', 'value', allow_duplicate=True)],
    [Input('url', 'search')],
    [State('dataset-select', 'options'),
     State('dataset-select', 'value')],
    prevent_initial_call='initial_duplicate'
)
def initialize_from_url(search, dataset_options, current_dataset):
    if not search:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update
    
    try:
        state = parse_url_state(search)
        if not state:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
        return (
            state.get('dataset'),
            state.get('embedding'),
            state.get('color'),
            state.get('mode', 'cells'),
            state.get('gene')
        )
    except Exception as e:
        logger.error(f"Error initializing from URL: {str(e)}")
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

# Callback to handle cell selection
@callback(
    [Output('selection-store', 'data'),
     Output('selection-info', 'children')],
    Input('main-plot', 'selectedData'),
    State('data-store', 'data')
)
def update_selection(selected_data, data_store):
    if not selected_data or not data_store:
        return None, "No cells selected"
    
    points = selected_data.get('points', [])
    n_selected = len(points)
    
    selection_data = {
        'indices': [p['pointIndex'] for p in points] if points else []
    }
    
    return selection_data, f"Selected {n_selected} cells"

# Update the main plot callback to handle view state
@callback(
    Output('main-plot', 'figure'),
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('color-select', 'value'),
     Input('viz-mode', 'value'),
     Input('selection-store', 'data'),
     Input('url', 'search')]
)
def update_plot(data_store, embedding, color_by, viz_mode, selection_data, url_search):
    if data_store is None or embedding is None:
        return {}
    
    try:
        adata = load_adata(data_store['filename'])
        
        # Handle custom embedding with better error checking
        if embedding == 'custom_embedding':
            if not custom_x or not custom_y:
                logger.warning("Custom embedding x or y not selected")
                return {}
            x = adata.obs[custom_x].values
            y = adata.obs[custom_y].values
            embedding_name = f'{custom_x} vs {custom_y}'
        else:
            if embedding is None:
                logger.warning("Embedding is None")
                return {}
            x = adata.obsm[embedding][:, 0]
            y = adata.obsm[embedding][:, 1]
            embedding_name = embedding
        
        # Handle color series with proper None checking
        if color_by == 'gene_expression':
            if not gene:
                color_series = None
            else:
                color_series = adata[:, gene].X.toarray().flatten() if scipy.sparse.issparse(adata.X) else adata[:, gene].X
                treat_as_categorical = False
        else:
            color_series = adata.obs[color_by] if color_by else None
            if color_series is not None:
                treat_as_categorical = (
                    color_series.dtype.name in ['category', 'object'] or
                    (pd.api.types.is_integer_dtype(color_series) and len(color_series.unique()) <= 50)
                )
        
        # Create DataFrame with proper handling of None values
        df = pd.DataFrame({
            'x': x,
            'y': y,
            'color': color_series if color_series is not None else pd.Series([None] * len(x))
        })
        
        # Apply visualization mode ordering if specified
        if viz_mode in ['ordered_asc', 'ordered_desc']:
            if color_series is not None and not treat_as_categorical:
                df = df.sort_values('color', 
                                  ascending=(viz_mode == 'ordered_asc'))
        
        # Create the plot using the utility function
        fig = create_scatter_plot(df, embedding, color_by, 
                                treat_as_categorical, selection_data)
        
        # Apply view state from URL if available
        if url_search:
            params = parse_qs(url_search.lstrip('?'))
            try:
                if 'view' in params:
                    view_state = json.loads(base64.urlsafe_b64decode(params['view'][0]).decode('utf-8'))
                    if 'xrange' in view_state and 'yrange' in view_state:
                        fig.update_layout(
                            xaxis_range=view_state['xrange'],
                            yaxis_range=view_state['yrange']
                        )
            except Exception as e:
                logger.warning(f"Failed to apply view state: {e}")
        
        return fig
    
    except Exception as e:
        logger.error(f"Error updating plot: {str(e)}")
        return {}

# Add this after app initialization
@server.route('/download/<path:filepath>')
def download_file(filepath):
    try:
        # Ensure the filepath is safe and within the data directory
        safe_path = os.path.join('data', os.path.basename(filepath))
        if os.path.exists(safe_path):
            return send_file(
                safe_path,
                as_attachment=True,
                download_name=os.path.basename(filepath)
            )
        else:
            return "File not found", 404
    except Exception as e:
        return str(e), 500

if __name__ == '__main__':
    app.run_server(debug=True) 