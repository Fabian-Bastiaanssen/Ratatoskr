import getpass
from glob import glob
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

from loguru import logger   

class suppress_stdout:
    def __enter__(self):
        self._original = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args):
        sys.stdout.close()
        sys.stdout = self._original


def src_local(rel_path):
    return os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path))


def get_version():
    with open(src_local("VERSION"), "r") as f:
        version = f.readline().strip()
    return version


def make_dir(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def unzip_file(zip_path, extract_to):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)


def move_and_rename(pattern, target_dir, label):
        for src in glob(str(pattern)):
            src_path = Path(src)
            rename = src_path.parent.name
            dest = target_dir / f"{rename}{src_path.suffix}"
            try:
                shutil.move(str(src_path), str(dest))
            except Exception as e:
                logger.warning(f"Failed to move {label} {src_path} -> {dest}: {e}")


def get_credentials(email=True, password=True, api_key=False, api_being_accessed="LPSN", dev_mode=False):

    if dev_mode:
        if email:
            email = os.getenv("LPSN_DEV_EMAIL")
            try:
                email = email.strip()
            except AttributeError:
                logger.error("LPSN_DEV_EMAIL environment variable is not set or is invalid.")
                sys.exit(1)
        if password:
            password = os.getenv("LPSN_DEV_PASSWORD")
            try:
                password = password.strip()
            except AttributeError:
                logger.error("LPSN_DEV_PASSWORD environment variable is not set or is invalid.")
                sys.exit(1)
        if api_key:
            api_key = os.getenv("GENBANK_DEV_API_KEY")
            try:
                api_key = api_key.strip()
            except AttributeError:
                logger.error("GENBANK_DEV_API_KEY environment variable is not set or is invalid.")
                sys.exit(1)
        return email, password, api_key
    
    else:
        logger.info(f"Alternatively use the -d flag with environment variables: LPSN_EMAIL, LPSN_PASSWORD, GENBANK_API_KEY")
        if email:
            email = input(f"Email for {api_being_accessed} account: ").strip()
        if password:
            password = getpass.getpass(f"{api_being_accessed} Password: ")
        if api_key:
            api_key = getpass.getpass(f"{api_being_accessed} API Key: ")
    return email, password, api_key


def run_supbrocess(command: list, suppress_output: bool = True) -> int:

    if suppress_output:
        with suppress_stdout():
            result = subprocess.run(command, )
    else:
        result = subprocess.run(command)
    
    return result.returncode


def delete_thing(path: str) -> None:
    p = Path(path)
    
    if not path.exists() and not path.is_symlink():
        return

    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()

