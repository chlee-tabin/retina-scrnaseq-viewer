import logging
from functools import wraps
from dash.exceptions import PreventUpdate

logger = logging.getLogger(__name__)

def handle_callback_error(func):
    """Decorator to handle callback errors gracefully"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except PreventUpdate:
            raise
        except Exception:
            logger.exception("Error in %s", func.__name__)
            raise PreventUpdate
    return wrapper

def log_callback_info(func):
    """Decorator to log callback information"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.debug("Callback %s triggered", func.__name__)
        result = func(*args, **kwargs)
        logger.debug(f"Callback {func.__name__} completed")
        return result
    return wrapper
