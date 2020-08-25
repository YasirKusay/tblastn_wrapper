PIP = pip3
PYTHON = python3
SPHINX = sphinx-build

.PHONY: setup
setup:
	$(PIP) install -r requirements.txt

.PHONY: setup-dev
setup-dev:
	$(PIP) install -r requirements-dev.txt

install: $(wildcard tblastn_wrapper/*.py)
	$(PIP) install -e .

.PHONY: run
run: install
	tblastn_wrapper

.PHONY: test
test:
	$(PYTHON) -m unittest tests/*.py

.PHONY: docs
docs:
	sphinx-build -b html docs/ docs/_build

.PHONY: format
format:
	black tblastn_wrapper
