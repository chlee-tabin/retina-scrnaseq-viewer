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
    """Create a shareable URL carrying the full encoded state (base64 JSON).

    The caller (callbacks.url_callbacks.share_url) decides which keys to include;
    this function no longer strips or derives any of them, so the expression-by-group
    and figure controls survive the round-trip.
    """
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


def state_for_dataset(url_search, dataset_id):
    """Parse shared-URL state, but only when it belongs to ``dataset_id`` (else None).

    The single seam for the stale-state guard: a ``?state=`` from a different dataset
    must not re-apply when the user switches datasets (the URL is never cleared). Every
    URL-reading restore callback parses through this rather than bare ``parse_url_state``
    so the guard cannot be forgotten by a new reader.
    """
    state = parse_url_state(url_search)
    if state and state.get('dataset') != dataset_id:
        return None
    return state