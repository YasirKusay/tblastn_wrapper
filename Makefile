PIP = pip3
PYTHON = python3
SPHINX = sphinx-build
POETRY = poetry
BLACK = black

install: pyproject.toml
	$(POETRY) install

build: $(wildcard tblastn_wrapper/*.py) pyproject.toml
	$(POETRY) build

.PHONY: run
run: install
	$(POETRY) run tblastn_wrapper

.PHONY: test
test:
	$(PYTHON) -m unittest tests/*.py

.PHONY: clean
clean:
	rm -rf dist

.PHONY: docs
docs:
	sphinx-build -b html docs/ docs/_build

.PHONY: typecheck
typecheck:
	$(POETRY) run mypy -p tblastn_wrapper --ignore-missing-imports

.PHONY: format
format:
	$(BLACK) tblastn_wrapper
