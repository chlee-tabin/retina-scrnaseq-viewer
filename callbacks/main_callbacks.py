from dash import Input, Output, State, callback
import json
import base64
from urllib.parse import parse_qs, urlencode
from utils.data_loading import load_adata, load_dataset_config
from utils.plotting import create_scatter_plot, create_metacell_plot, create_binned_plot
import pandas as pd
import logging
from utils.processing import create_metacells, calculate_selection_stats
from utils.error_handling import handle_callback_error, log_callback_info
from utils.state import parse_url_state
import scipy.sparse
import traceback
import plotly.graph_objects as go

logger = logging.getLogger(__name__)


def _message_figure(text):
    """Return an empty plot displaying a centered, user-friendly message."""
    fig = go.Figure()
    fig.add_annotation(
        text=text,
        xref="paper", yref="paper",
        x=0.5, y=0.5, showarrow=False,
        font=dict(size=16, color="#666")
    )
    fig.update_layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor="white"
    )
    return fig

@callback(
    Output('main-plot', 'figure', allow_duplicate=True),
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('custom-x-select', 'value'),
     Input('custom-y-select', 'value'),
     Input('color-select', 'value'),
     Input('gene-select', 'value'),
     Input('viz-mode', 'value'),
     Input('bin-number-slider', 'value'),
     Input('percentile-slider', 'value'),
     Input('enable-binning', 'value'),
     Input('enable-smoothing', 'value'),
     Input('selection-store', 'data'),
     Input('url', 'search')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_plot(data_store, embedding, custom_x, custom_y, color_by, gene, viz_mode, bin_number, percentile, enable_binning, enable_smoothing, selection_data, url_search):
    logger.debug("update_plot called with parameters:")
    
    # Initialize treat_as_categorical as False by default
    treat_as_categorical = False
    
    if not data_store or not embedding:
        return {}
        
    try:
        config = load_dataset_config()
        adata = load_adata(data_store['filename'])
        
        # Get coordinates based on embedding type
        if embedding == 'custom_embedding':
            if not custom_x or not custom_y:
                return {}
            x = adata.obs[custom_x]
            y = adata.obs[custom_y]
            x_label = custom_x
            y_label = custom_y
        else:
            coordinates = adata.obsm[embedding]
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            x_label = f"{embedding}_1"
            y_label = f"{embedding}_2"
        
        # Handle gene expression
        if color_by == 'gene_expression':
            logger.debug("Color by gene_expression selected")
            if not gene:
                logger.debug("No gene selected")
                color_series = None
            elif gene not in adata.var_names:
                # Guard against genes absent from this dataset (e.g. a shared
                # URL, a stale selection, or wrong gene-name casing across
                # species). Return a friendly message instead of a 500.
                logger.info(f"Gene '{gene}' not found in dataset '{data_store['filename']}'")
                return _message_figure(
                    f"Gene '{gene}' not found in this dataset. "
                    "Try a different gene or check the name/casing."
                )
            else:
                if scipy.sparse.issparse(adata.X):
                    color_series = adata[:, gene].X.toarray().flatten()
                    logger.debug(f"Gene '{gene}' expression data loaded from sparse matrix")
                else:
                    color_series = adata[:, gene].X
                    logger.debug(f"Gene '{gene}' expression data loaded from dense matrix")
                treat_as_categorical = False
        else:
            logger.debug(f"Color by '{color_by}' selected")
            color_series = adata.obs[color_by] if color_by else None
            if color_series is not None:
                logger.debug(f"Color series dtype: {color_series.dtype}")
                treat_as_categorical = (
                    color_series.dtype.name in ['category', 'object'] or
                    (pd.api.types.is_integer_dtype(color_series) and len(color_series.unique()) <= 50)
                )
            logger.debug(f"treat_as_categorical: {treat_as_categorical}")
        
        # Create DataFrame for plotting
        df = pd.DataFrame({
            'x': x,
            'y': y,
            'color': color_series
        })
        
        # Create the plot
        if enable_binning:
            smooth_on = bool(enable_smoothing) and 'enabled' in enable_smoothing
            smooth_sigma = float(data_store.get('smooth_sigma', 0) or 0) if smooth_on else 0
            min_cells = data_store.get('min_cells_per_bin', 1)
            fig = create_binned_plot(
                df, embedding, color_by,
                bin_size=bin_number,
                percentile=percentile,
                treat_as_categorical=treat_as_categorical,
                smooth_sigma=smooth_sigma,
                min_cells=min_cells
            )
        else:
            fig = create_scatter_plot(
                df, embedding, color_by,
                treat_as_categorical=treat_as_categorical
            )
        
        return fig
        
    except Exception as e:
        logger.error(f"Error updating plot: {str(e)}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        return {}

@callback(
    [Output('custom-x-select', 'options'),
     Output('custom-y-select', 'options'),
     Output('custom-x-select', 'value'),
     Output('custom-y-select', 'value'),
     Output('custom-embedding-container', 'style')],
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('url', 'search')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_custom_embedding_controls(data_store, embedding, url_search):
    if not data_store or not data_store.get('metadata_cols'):
        return [], [], None, None, {'display': 'none'}
    
    # Get numeric columns for custom embedding
    numeric_cols = [
        col for col, type_ in data_store.get('column_types', {}).items()
        if type_ == 'numeric'
    ]
    
    options = [{'label': col, 'value': col} for col in numeric_cols]
    
    # Set container visibility based on embedding type
    container_style = {'display': 'block'} if embedding == 'custom_embedding' else {'display': 'none'}
    
    # Check if we should initialize from URL state
    if url_search:
        state = parse_url_state(url_search)
        if state and state.get('embedding') == 'custom_embedding':
            custom_x = state.get('custom_x')
            custom_y = state.get('custom_y')
            # Only set values if they exist in the options
            valid_values = [opt['value'] for opt in options]
            if custom_x in valid_values and custom_y in valid_values:
                return options, options, custom_x, custom_y, container_style
    
    # Default return if not initializing from URL: pre-select the topographic
    # axes (NT.Score on x, DV.Score on y) when those score columns exist.
    default_x = 'NT.Score' if 'NT.Score' in numeric_cols else None
    default_y = 'DV.Score' if 'DV.Score' in numeric_cols else None
    return options, options, default_x, default_y, container_style

@callback(
    Output('viz-mode', 'style'),
    [Input('embedding-select', 'value'),
     Input('enable-binning', 'value')]
)
def toggle_plot_order(embedding, enable_binning):
    if embedding == 'custom_embedding' and enable_binning and 'enabled' in enable_binning:
        return {'display': 'none'}
    return {'display': 'block'}

@callback(
    Output('binning-controls', 'style'),
    [Input('enable-binning', 'value')]
)
def toggle_binning_controls(enable_binning):
    if enable_binning and 'enabled' in enable_binning:
        return {'display': 'block'}
    return {'display': 'none'} 