#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The setup script."""
from setuptools import find_packages, setup

NAME = "jet-stream-metrics"
DESCRIPTION = "Implementations of various jet-stream metrics from literature."
URL = "https://github.com/Thomasjkeel/jet-stream-metrics"
AUTHOR = "Thomas Keel"
AUTHOR_EMAIL = "thomas.keel.18@ucl.ac.uk"
REQUIRES_PYTHON = ">=3.8.0"
VERSION = "0.1.0"
LICENSE = "MIT License"

KEYWORDS = "jet-stream metrics climate xarray"

with open("README.md") as file:
    long_description = file.read()

dev_requirements = [
    "black",
    "coverage",
    "flake8",
    "parameterized",
    "pre_commit",
    "pylint",
    "pytest",
]

requirements = [
    "numpy>=1.21.2",
    "pandas>=0.23",
    "xarray>=0.19.0",
    "scipy>=1.5.3",
    "dask[array]>=2.6",
    "netCDF4>=1.5.5.1",
    "bottleneck==1.3.2",
    "cf_xarray",
    "Shapely",
]

setup(
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    classifiers=[
        "Development Status :: 3 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
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
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords=KEYWORDS,
    name=NAME,
    packages=find_packages(include=["metrics"]),
    tests_require=["pytest", "parameterized"],
    extras_require={
        "dev": dev_requirements,
        "plotting": [
            "matplotlib>=3.3.2",
            # "Cartopy>=0.18.0",
        ],
        "data_install": [
            "cdsapi>=0.5.1",
        ],
    },
    url=URL,
    version=VERSION,
    zip_safe=False,
)
