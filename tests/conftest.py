"""Keep tests synthetic and all temporary/cache files inside this worktree."""
import os
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ['HF_DATA_REPO'] = ''
os.environ['MPLCONFIGDIR'] = str(ROOT / '.test-mpl-cache')


def pytest_configure(config):
    if config.option.basetemp is None:
        config.option.basetemp = str(ROOT / '.test-tmp')


def pytest_sessionfinish(session, exitstatus):
    shutil.rmtree(ROOT / '.test-mpl-cache', ignore_errors=True)
    shutil.rmtree(ROOT / '.test-tmp', ignore_errors=True)
