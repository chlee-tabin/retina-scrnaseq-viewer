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

logger = logging.getLogger(__name__)

@callback(
    Output('main-plot', 'figure', allow_duplicate=True),
    [Input('data-store', 'data'),
     Input('embedding-select', 'value'),
     Input('color-select', 'value'),
     Input('viz-mode', 'value'),
     Input('selection-store', 'data'),
     Input('url', 'search')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_plot(data_store, embedding, color_by, viz_mode, selection_data, url_search):
    if data_store is None or embedding is None:
        return {}
    
    try:
        adata = load_adata(data_store['filename'])
        
        x = adata.obsm[embedding][:, 0]
        y = adata.obsm[embedding][:, 1]
        
        color_series = adata.obs[color_by] if color_by else None
        treat_as_categorical = False
        if color_series is not None:
            if color_series.dtype.name in ['category', 'object']:
                treat_as_categorical = True
            elif pd.api.types.is_integer_dtype(color_series):
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
            
            fig = create_scatter_plot(df, embedding, color_by, treat_as_categorical)
            
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