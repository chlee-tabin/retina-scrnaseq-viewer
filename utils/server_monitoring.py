import psutil
import threading
from collections import defaultdict
import time

# Store active sessions with last activity timestamp
active_sessions = defaultdict(lambda: time.time())
_lock = threading.Lock()

def get_memory_usage():
    """Get current memory usage of the process"""
    process = psutil.Process()
    memory_info = process.memory_info()
    return f"{memory_info.rss / (1024 * 1024):.1f} MB"

def update_session(session_id):
    """Update or add session last activity"""
    with _lock:
        active_sessions[session_id] = time.time()
        # Clean up sessions older than 5 minutes
        current_time = time.time()
        expired = [sid for sid, last_active in active_sessions.items() 
                  if current_time - last_active > 300]
        for sid in expired:
            del active_sessions[sid]
    return len(active_sessions)

def get_active_users():
    """Get number of active users in last 5 minutes"""
    with _lock:
        current_time = time.time()
        return sum(1 for last_active in active_sessions.values() 
                  if current_time - last_active <= 300) 