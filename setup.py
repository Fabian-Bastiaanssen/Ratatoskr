import os
from distutils.core import Extension

from setuptools import find_packages, setup

def src_local(rel_path):
    return os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path))

def get_version():
    with open(src_local("ratatoskr/VERSION"), "r") as f:
        version = f.readline().strip()
    return version

def get_description():
    with open("README.md", "r") as README:
        description = README.read()

data_files = [(".", ["LICENSE", "README.md", "ratatoskr/VERSION"])]

setup(
    name="ratatoskr",
    # packages=find_packages(where="src"),
    url="https://github.com/Fabian-Bastiaanssen/Ratatoskr",
    python_requires=">=3.12",
    description="Ratatoskr: A tool for collecting and downloading taxonomic type stains.",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Chris Turkington and Fabian Bastiaanssen",
    author_email="chrisjrt1@gmail.com",
    data_files=data_files,
    packages=find_packages(),
    install_requires=[
        "async-dsmz>=2025.0.5",
        "bacdive==1.0.0",
        "biopython>=1.86",
        "loguru>=0.7.2",
        "lpsn>=1.0.0",
        "polars>=1.32.2",
        "pyarrow>=13.0.0",
        "rich-click>=1.8.8",
        "tqdm>=4.66.1",
    ],
    entry_points={"console_scripts": ["ratatoskr=ratatoskr.__main__:main"]},
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
    ],
)
