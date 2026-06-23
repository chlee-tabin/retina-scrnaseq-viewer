from dash import Input, Output, State, callback
from utils.error_handling import handle_callback_error, log_callback_info
from utils.data_loading import gene_search_options
import logging

logger = logging.getLogger(__name__)

# NOTE: dataset-select options, data-store loading, and dataset-info rendering
# (including the download link) are handled by the authoritative callbacks in
# app.py (initialize_dataset_dropdown, update_data, update_dataset_info). The
# gene-search callback below is the sole owner of its outputs; its fuzzy search is
# the shared gene_search_options helper (also used by the group / second-gene
# dropdowns in main_callbacks).


@callback(
    [Output('gene-select-container', 'style', allow_duplicate=True),
     Output('gene-select', 'options', allow_duplicate=True)],
    [Input('color-select', 'value'),
     Input('data-store', 'data'),
     Input('gene-select', 'search_value')],
    prevent_initial_call=True
)
@handle_callback_error
@log_callback_info
def update_gene_select(color_value, data_store, search_value):
    # The gene picker is only meaningful when colouring by gene expression.
    if color_value != 'gene_expression' or not data_store or 'genes' not in data_store:
        return {'display': 'none'}, []
    return {'display': 'block'}, gene_search_options(data_store['genes'], search_value)
