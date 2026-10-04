"""Upload a clean, tracked checkout to a Docker Space; never copies local data."""
import argparse
import io
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile


def deploy_space(space_id, checkout=None):
    checkout = Path(checkout or Path(__file__).resolve().parents[1])

    def git(*args):
        return subprocess.check_output(['git', *args], cwd=checkout)

    if git('status', '--porcelain', '--untracked-files=normal').strip():
        raise RuntimeError('Refusing deployment: the git tree is dirty.')
    token = os.getenv('HF_TOKEN')
    if not token:
        raise RuntimeError('Set HF_TOKEN to upload to the Space.')
    sha = git('rev-parse', 'HEAD').decode().strip()
    from huggingface_hub import upload_folder

    with tempfile.TemporaryDirectory(prefix='viewer-deploy-', dir=checkout) as temp:
        # git archive includes tracked files only, excluding ignored data and secrets.
        with tarfile.open(fileobj=io.BytesIO(git('archive', '--format=tar', sha))) as archive:
            archive.extractall(temp, filter='data')
        (Path(temp) / 'viewer_revision.txt').write_text(sha + '\n')
        result = upload_folder(repo_id=space_id, repo_type='space', folder_path=temp,
                               token=token, commit_message=f'Deploy viewer {sha[:7]}',
                               delete_patterns=['*'])
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('space_id', help='Hugging Face Space id (owner/name)')
    args = parser.parse_args()
    print(deploy_space(args.space_id))
