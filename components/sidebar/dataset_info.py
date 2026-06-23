from dash import html, dcc
import dash_bootstrap_components as dbc
import logging
from utils.data_loading import load_dataset_config, validate_datasets

logger = logging.getLogger(__name__)

def create_dataset_info():
    # Populate the dataset dropdown at layout-build time and pre-select the
    # first dataset, so the app opens on a dataset instead of a blank state.
    try:
        datasets = validate_datasets(load_dataset_config())
        dataset_options = [{'label': d['title'], 'value': i} for i, d in datasets.items()]
        default_dataset = next(iter(datasets), None)
    except Exception as e:
        logger.error(f"Could not pre-load dataset options: {e}")
        dataset_options, default_dataset = [], None

    return html.Div([
        # Dataset selector dropdown
        html.Label("Dataset:"),
        dcc.Dropdown(
            id='dataset-select',
            options=dataset_options,
            value=default_dataset,
            placeholder="Select dataset",
            className="mb-3"
        ),

        # Dataset info card
        dbc.Card([
            dbc.CardBody([
                html.Div(id='dataset-info'),
            ])
        ], className="mt-3")
    ])
