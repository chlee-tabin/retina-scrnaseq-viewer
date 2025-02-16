import json
import base64
from urllib.parse import parse_qs, urlencode

def encode_state(state_dict):
    """Encode application state for URL sharing"""
    state_json = json.dumps(state_dict)
    return base64.urlsafe_b64encode(state_json.encode()).decode()

def decode_state(encoded_state):
    """Decode application state from URL"""
    try:
        state_json = base64.urlsafe_b64decode(encoded_state).decode()
        return json.loads(state_json)
    except Exception:
        return None

def create_share_url(base_url, state_dict):
    """Create shareable URL with encoded state"""
    encoded_state = encode_state(state_dict)
    return f"{base_url}?state={encoded_state}"

def parse_url_state(url_search):
    """Parse state from URL search parameters"""
    if not url_search:
        return None
    
    params = parse_qs(url_search.lstrip('?'))
    if 'state' not in params:
        return None
    
    return decode_state(params['state'][0]) 