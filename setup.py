#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The setup script."""
from setuptools import find_packages, setup

NAME = "jsmetrics"
DESCRIPTION = "Library of algorithms and metrics used to characterise and/or identify jet-streams, based on xarray."
URL = "https://github.com/Thomasjkeel/jsmetrics"
AUTHOR = "Thomas Keel"
AUTHOR_EMAIL = "thomas.keel.18@ucl.ac.uk"
REQUIRES_PYTHON = ">=3.7.0"
VERSION = "0.1.7"
LICENSE = "MIT License"

KEYWORDS = "jet-stream climate metrics algorithms xarray"

with open("README.rst") as file:
    readme = file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

# dev_requirements = []
# with open("requirements_dev.txt") as dev:
#     for dependency in dev.readlines():
#         dev_requirements.append(dependency)

requirements = [
    "numpy>=1.21.2",
    "pandas>=0.23",
    "matplotlib>=3.3.2",
    "xarray>=0.19.0",
    "scipy>=1.5.3",
    "dask[array]",
    "netCDF4>=1.5.5.1",
    "bottleneck",
    "cf_xarray",
    "Shapely",
]

dev_requirements = [
    "bump2version",
    "black",
    "pytest",
    "pytest-cov",
    "flake8",
    "parameterized",
    "pre_commit",
    "pylint",
    "pytest",
    "sphinx",
    "numpydoc",
    "sphinx-rtd-theme",
]

setup(
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
    description=DESCRIPTION,
    python_requires=REQUIRES_PYTHON,
    install_requires=requirements,
    license=LICENSE,
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords=KEYWORDS,
    name=NAME,
    packages=find_packages(),
    tests_require=["pytest", "parameterized"],
    extras_require={
        "dev": dev_requirements,
        "data_install": [
            "cdsapi>=0.5.1",
        ],
    },
    url=URL,
    version=VERSION,
    zip_safe=False,
)
