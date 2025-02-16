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

logger = logging.getLogger(__name__)

@callback(
    Output('main-plot', 'figure', allow_duplicate=True),
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('color-select', 'value'),
     Input('gene-select', 'value'),
     Input('viz-mode', 'value'),
     Input('selection-store', 'data'),
     Input('url', 'search')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_plot(data_store, embedding, color_by, gene, viz_mode, selection_data, url_search):
    if data_store is None or embedding is None:
        return {}
    
    try:
        adata = load_adata(data_store['filename'])
        
        x = adata.obsm[embedding][:, 0]
        y = adata.obsm[embedding][:, 1]
        
        # Handle gene expression
        if color_by == 'gene_expression':
            if not gene:  # If no gene is selected yet
                color_series = None
            else:
                color_series = adata[:, gene].X.toarray().flatten() if scipy.sparse.issparse(adata.X) else adata[:, gene].X
            treat_as_categorical = False
        else:
            color_series = adata.obs[color_by] if color_by else None
            treat_as_categorical = (
                color_series.dtype.name in ['category', 'object'] or
                (pd.api.types.is_integer_dtype(color_series) and len(color_series.unique()) <= 50)
            )
        
        if viz_mode == 'cells':
            df = pd.DataFrame({
                'x': x,
                'y': y,
                'color': color_series
            })
            
            fig = create_scatter_plot(df, embedding, color_by if color_by != 'gene_expression' else gene, 
                                    treat_as_categorical)
            
            if selection_data and selection_data['indices']:
                fig.update_traces(
                    selectedpoints=selection_data['indices'],
                    selected=dict(marker=dict(color='red')),
                    unselected=dict(marker=dict(opacity=0.3))
                )
        else:
            values = adata.obs[color_by].values if color_by else None
            H, xedges, yedges = create_metacells(x, y, values)
            fig = create_metacell_plot(H, xedges, yedges, embedding)
        
        # Apply view state from URL
        view_state = parse_url_state(url_search)
        if view_state and 'view' in view_state:
            fig.update_layout(
                xaxis_range=view_state['view']['xrange'],
                yaxis_range=view_state['view']['yrange']
            )
        
        return fig
    
    except Exception as e:
        logger.error(f"Error updating plot: {str(e)}")
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
    if color_value == 'gene_expression' and data_store and 'genes' in data_store:
        gene_options = [{'label': gene, 'value': gene} for gene in data_store['genes']]
        return {'display': 'block'}, gene_options
    return {'display': 'none'}, [] 