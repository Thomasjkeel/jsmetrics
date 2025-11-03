.PHONY: clean clean-test clean-pyc clean-build docs help install test test-all lint lint/flake8 lint/black coverage dist release servedocs

.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys
from urllib.request import pathname2url
webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys
for line in sys.stdin:
    match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
    if match:
        target, help = match.groups()
        print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help: ## show this help
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

# ----------------------------
# Clean targets
# ----------------------------
clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/ dist/ .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/ htmlcov/ .pytest_cache
	rm -f .coverage

# ----------------------------
# Linting
# ----------------------------
lint/flake8: ## check style with flake8
	flake8 jsmetrics tests

lint/black: ## check style with black
	black --check jsmetrics tests

lint: lint/flake8 lint/black ## check style

# ----------------------------
# Testing
# ----------------------------
test: ## run tests quickly with the default Python
	pytest

test-all: ## run tests on every Python version with tox
	tox

# ----------------------------
# Coverage
# ----------------------------
coverage: ## check code coverage quickly with the default Python
	coverage run --source jsmetrics -m pytest
	coverage report -m
	coverage html
	$(BROWSER) htmlcov/index.html

# ----------------------------
# Documentation
# ----------------------------
docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/jsmetrics.rst docs/modules.rst
	sphinx-apidoc -o docs/ jsmetrics
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

# ----------------------------
# Build / Distribution
# ----------------------------
dist: clean ## builds source and wheel package
	python -m build
	ls -l dist

release: dist ## package and upload a release
	twine upload dist/*

# ----------------------------
# Installation
# ----------------------------
install: clean ## install the package in editable mode with dev deps
	python -m pip install --upgrade pip
	pip install -e .[dev]
