import yaml
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_CONFIG = {
    'plot_settings': {
        'height': '85vh',
        'marker_size': 5,
        'opacity': 0.7,
        'colorscales': {
            'continuous': 'viridis',
            'categorical': 'Set3'
        },
        'default_viz_mode': 'random'
    },
    'metacell_settings': {
        'bin_size': 50,
        'smoothing': 1.0
    }
}

def load_config(config_path='config.yml'):
    try:
        with open(config_path, 'r') as f:
            user_config = yaml.safe_load(f)
        # Merge user config with defaults
        config = DEFAULT_CONFIG.copy()
        config.update(user_config)
        return config
    except FileNotFoundError:
        logger.warning(f"Config file {config_path} not found, using defaults")
        return DEFAULT_CONFIG
    except Exception as e:
        logger.error(f"Error loading config: {str(e)}")
        return DEFAULT_CONFIG 