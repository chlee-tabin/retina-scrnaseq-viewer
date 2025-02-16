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

# Import layouts
from layouts.sidebar import create_sidebar
from layouts.main import create_main_panel

# Import callbacks
from callbacks.dataset_callbacks import *
from callbacks.selection_callbacks import *
from callbacks.url_callbacks import *
from callbacks.main_callbacks import *

# Import utilities
from utils.data_loading import load_adata, load_dataset_config, validate_datasets
from utils.config import load_config
from utils.error_handling import handle_callback_error, log_callback_info

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize the Dash app with bootstrap theme
app = dash.Dash(
    __name__, 
    external_stylesheets=[dbc.themes.BOOTSTRAP], 
    suppress_callback_exceptions=True
)
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
    [Output('data-store', 'data'),
     Output('embedding-select', 'options'),
     Output('color-select', 'options'),
     Output('loading-output', 'children')],
    Input('dataset-select', 'value'),
    prevent_initial_call=True
)
def update_data(dataset_id):
    if not dataset_id:
        return None, [], [], ""
    
    try:
        config = load_dataset_config()
        dataset = config['datasets'][dataset_id]
        file_path = dataset['file_path']
        
        logger.info(f"Loading dataset: {dataset['title']}")
        adata = load_adata(file_path)
        
        # Common processing
        embeddings = list(adata.obsm.keys())
        metadata_cols = list(adata.obs.columns)
        
        data_store = {
            'filename': file_path,
            'n_cells': adata.n_obs,
            'embeddings': embeddings,
            'metadata_cols': metadata_cols
        }
        
        embedding_options = [{'label': emb, 'value': emb} for emb in embeddings]
        color_options = [{'label': col, 'value': col} for col in metadata_cols]
        
        logger.info(f"Successfully loaded dataset with {adata.n_obs} cells")
        return data_store, embedding_options, color_options, ""  # Changed to return empty string
        
    except Exception as e:
        error_message = f"Error loading data: {str(e)}"
        logger.error(error_message)
        return None, [], [], error_message

# Update the initialization callback to be more robust
@callback(
    [Output('dataset-select', 'value', allow_duplicate=True),
     Output('embedding-select', 'value', allow_duplicate=True),
     Output('color-select', 'value', allow_duplicate=True),
     Output('viz-mode', 'value', allow_duplicate=True)],
    [Input('url', 'search')],
    [State('dataset-select', 'options'),
     State('dataset-select', 'value')],
    prevent_initial_call='initial_duplicate'
)
def initialize_from_url(search, dataset_options, current_dataset):
    if not search:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update
    
    try:
        params = parse_qs(search.lstrip('?'))
        
        dataset = params.get('dataset', [None])[0]
        embedding = params.get('embedding', [None])[0]
        color_by = params.get('color', [None])[0]
        viz_mode = params.get('mode', ['cells'])[0]
        
        # Only update if we have valid parameters
        if not dataset or not embedding or not color_by:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
        # Verify dataset exists in options
        if dataset_options:
            valid_datasets = [opt['value'] for opt in dataset_options]
            if dataset not in valid_datasets:
                logger.warning(f"Dataset {dataset} not found in options")
                return dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
        # Only update if values are different from current state
        if dataset == current_dataset:
            logger.info("Dataset already selected, skipping update")
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
        logger.info(f"Initializing from URL with dataset={dataset}, embedding={embedding}, color={color_by}, mode={viz_mode}")
        return dataset, embedding, color_by, viz_mode
        
    except Exception as e:
        logger.error(f"Error initializing from URL: {str(e)}")
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update

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
        
        x = adata.obsm[embedding][:, 0]
        y = adata.obsm[embedding][:, 1]
        
        # Determine if color column should be treated as categorical
        color_series = adata.obs[color_by] if color_by else None
        treat_as_categorical = False
        if color_series is not None:
            if color_series.dtype.name in ['category', 'object']:
                treat_as_categorical = True
            elif np.issubdtype(color_series.dtype, np.integer):
                n_unique = len(color_series.unique())
                if n_unique <= 50:
                    treat_as_categorical = True
                    color_series = color_series.astype('category')
        
        if viz_mode == 'cells':
            df = pd.DataFrame({
                'x': x,
                'y': y,
                'color': color_series if color_by else None
            })
            
            fig = px.scatter(
                df, x='x', y='y', color='color',
                labels={'x': f'{embedding}_1', 'y': f'{embedding}_2'},
                title=f'Single-cell visualization - {embedding}',
                color_discrete_sequence=px.colors.qualitative.Set3 if treat_as_categorical else None,
                color_continuous_scale='viridis' if not treat_as_categorical else None,
                hover_data=None
            )
            
            # Add custom styling for UMAP plots
            if 'umap' in embedding.lower():
                fig.update_layout(
                    plot_bgcolor='white',
                    xaxis=dict(
                        showgrid=False,
                        showticklabels=False,
                        scaleanchor="y",
                        scaleratio=1,
                    ),
                    yaxis=dict(
                        showgrid=False,
                        showticklabels=False,
                        scaleanchor="x",
                        scaleratio=1,
                    )
                )
            else:
                # Maintain original 1:1 aspect ratio for non-UMAP plots
                fig.update_layout(
                    yaxis=dict(
                        scaleanchor="x",
                        scaleratio=1,
                    )
                )
            
            # Update selection styling
            if selection_data and selection_data['indices']:
                selected_indices = selection_data['indices']
                fig.update_traces(
                    selectedpoints=selected_indices,
                    selected=dict(marker=dict(color='red')),
                    unselected=dict(marker=dict(opacity=0.3))
                )
        
        else:
            values = adata.obs[color_by].values if color_by else None
            H, xedges, yedges = create_metacells(x, y, values)
            
            x_centers = (xedges[:-1] + xedges[1:]) / 2
            y_centers = (yedges[:-1] + yedges[1:]) / 2
            
            fig = px.imshow(
                H.T,
                x=x_centers,
                y=y_centers,
                labels={'x': f'{embedding}_1', 'y': f'{embedding}_2'},
                title=f'Metacell visualization - {embedding}',
                aspect='equal'  # This already maintains 1:1 ratio for imshow
            )
        
        # Parse URL to get view state and apply it
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

# Update the share button callback to only handle generating the shareable URL
@callback(
    [Output('share-url', 'style'),
     Output('share-url', 'value')],
    Input('share-button', 'n_clicks'),
    [State('dataset-select', 'value'),
     State('embedding-select', 'value'),
     State('color-select', 'value'),
     State('viz-mode', 'value'),
     State('main-plot', 'relayoutData'),
     State('url', 'href')]
)
def share_url(n_clicks, dataset, embedding, color_by, viz_mode, relay_data, current_url):
    if n_clicks is None:
        return {'display': 'none'}, ''
    
    # Parse the base URL (everything before the query string)
    base_url = current_url.split('?')[0]
    
    # Build query parameters
    params = {}
    if dataset:
        params['dataset'] = dataset
    if embedding:
        params['embedding'] = embedding
    if color_by:
        params['color'] = color_by
    if viz_mode:
        params['mode'] = viz_mode
    
    # Add view state if available
    if relay_data:
        view_state = {}
        if 'xaxis.range[0]' in relay_data and 'xaxis.range[1]' in relay_data:
            view_state['xrange'] = [relay_data['xaxis.range[0]'], relay_data['xaxis.range[1]']]
        if 'yaxis.range[0]' in relay_data and 'yaxis.range[1]' in relay_data:
            view_state['yrange'] = [relay_data['yaxis.range[0]'], relay_data['yaxis.range[1]']]
        
        if view_state:
            params['view'] = base64.urlsafe_b64encode(json.dumps(view_state).encode()).decode()
    
    # Construct the full URL
    full_url = f"{base_url}?{urlencode(params)}"
    
    return {'display': 'block', 'width': '100%', 'marginTop': '10px'}, full_url

if __name__ == '__main__':
    app.run_server(debug=True) 