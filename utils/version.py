"""Viewer release and the source revision stamped by the Space deployment script."""
import os
from pathlib import Path

VIEWER_VERSION = '0.33'


def viewer_revision():
    revision = os.getenv('VIEWER_GIT_SHA')
    if not revision:
        try:
            revision = (Path(__file__).resolve().parents[1] / 'viewer_revision.txt').read_text().strip()
        except OSError:
            revision = 'local'
    return revision[:7]
