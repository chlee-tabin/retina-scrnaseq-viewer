import logging
import traceback
from functools import wraps
from dash import html

logger = logging.getLogger(__name__)

def handle_callback_error(func):
    """Decorator to handle callback errors gracefully"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            error_msg = f"Error in {func.__name__}: {str(e)}"
            logger.error(f"{error_msg}\n{traceback.format_exc()}")
            
            # Return appropriate error message based on output type
            if 'figure' in func.__name__:
                return {}
            elif 'children' in func.__name__:
                return html.Div([
                    html.H4("Error", className="text-danger"),
                    html.P(error_msg)
                ])
            else:
                return None
    return wrapper

def log_callback_info(func):
    """Decorator to log callback information"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.debug(f"Callback {func.__name__} triggered with args: {args}, kwargs: {kwargs}")
        result = func(*args, **kwargs)
        logger.debug(f"Callback {func.__name__} completed")
        return result
    return wrapper 