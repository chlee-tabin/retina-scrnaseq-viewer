from dash import Input, Output, State, callback
import json
import base64
from urllib.parse import parse_qs, urlencode
from utils.data_loading import load_adata, load_dataset_config
from utils.plotting import create_scatter_plot, create_metacell_plot
import pandas as pd
import logging
from utils.processing import create_metacells, calculate_selection_stats
from utils.error_handling import handle_callback_error, log_callback_info
from utils.state import parse_url_state
import scipy.sparse
import traceback

logger = logging.getLogger(__name__)

@callback(
    Output('main-plot', 'figure', allow_duplicate=True),
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('custom-x-select', 'value'),
     Input('custom-y-select', 'value'),
     Input('color-select', 'value'),
     Input('gene-select', 'value'),
     Input('viz-mode', 'value'),
     Input('selection-store', 'data'),
     Input('url', 'search')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_plot(data_store, embedding, custom_x, custom_y, color_by, gene, viz_mode, selection_data, url_search):
    logger.debug("update_plot called with parameters:")
    logger.debug(f"data_store: {data_store}")
    logger.debug(f"embedding: {embedding}")
    logger.debug(f"custom_x: {custom_x}")
    logger.debug(f"custom_y: {custom_y}")
    logger.debug(f"color_by: {color_by}")
    logger.debug(f"gene: {gene}")
    logger.debug(f"viz_mode: {viz_mode}")
    logger.debug(f"selection_data: {selection_data}")
    logger.debug(f"url_search: {url_search}")

    if data_store is None:
        logger.warning("Data store is None")
        return {}
    
    try:
        adata = load_adata(data_store['filename'])
        logger.debug(f"Loaded adata with {adata.n_obs} cells and {adata.n_vars} genes")
        
        # Handle custom embedding
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
        
        logger.debug(f"Embedding '{embedding_name}' extracted: x[:5]={x[:5]}, y[:5]={y[:5]}")
        
        # Handle gene expression
        if color_by == 'gene_expression':
            logger.debug("Color by gene_expression selected")
            if not gene:
                logger.debug("No gene selected")
                color_series = None
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
        
        logger.debug("Creating DataFrame for plotting")
        df = pd.DataFrame({
            'x': x,
            'y': y,
            'color': color_series
        })
        
        # Order points based on visualization mode
        if viz_mode in ['ordered_asc', 'ordered_desc'] and color_series is not None:
            ascending = viz_mode == 'ordered_asc'
            if treat_as_categorical:
                # For categorical data, group by categories
                df = df.sort_values('color', ascending=ascending)
            else:
                # For continuous data, sort by value
                df = df.sort_values('color', ascending=ascending)
        elif viz_mode == 'random':
            # Randomly shuffle the points
            df = df.sample(frac=1, random_state=42)
        
        logger.debug(f"DataFrame head:\n{df.head()}")
        
        color_label = gene if color_by == 'gene_expression' else color_by
        logger.debug(f"Creating scatter plot with color_by='{color_label}'")
        fig = create_scatter_plot(df, embedding_name, color_label, treat_as_categorical)
        
        if selection_data and selection_data.get('indices'):
            logger.debug(f"Updating plot with selection indices: {selection_data['indices']}")
            fig.update_traces(
                selectedpoints=selection_data['indices'],
                selected=dict(marker=dict(opacity=1)),
                unselected=dict(marker=dict(opacity=0.1))
            )
        
        logger.debug("Plot created successfully")
        return fig
    
    except Exception as e:
        logger.error(f"Error updating plot: {str(e)}")
        logger.error(traceback.format_exc())
        return {} 

@callback(
    [Output('gene-select-container', 'style'),
     Output('gene-select', 'options')],
    [Input('color-select', 'value'),
     Input('data-store', 'data')]
)
@handle_callback_error
@log_callback_info
def update_gene_select(color_value, data_store):
    logger.debug("update_gene_select called with:")
    logger.debug(f"color_value: {color_value}")
    logger.debug(f"data_store: {data_store}")

    if color_value == 'gene_expression' and data_store and 'genes' in data_store:
        gene_options = [{'label': gene, 'value': gene} for gene in data_store['genes']]
        logger.debug(f"Gene options generated: {gene_options[:5]}...")  # Log first 5 for brevity
        return {'display': 'block'}, gene_options
    logger.debug("Hiding gene-select container")
    return {'display': 'none'}, [] 

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
    
    # Default return if not initializing from URL
    return options, options, None, None, container_style 