.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help


clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +
	rm -rf _skbuild
	rm -rf python/goss/_gosscpp*
	rm -rf python/goss/include
	rm -rf python/goss/lib
	rm -rf python/goss/libgoss*

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

install: clean
	python -m pip install -e .

install-no-deps: clean
	python -m pip install -e . --no-deps


docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/source/goss.rst´
	rm -f docs/source/modules.rst
	sphinx-apidoc -o docs/source python/goss
	for file in CONTRIBUTING.md; do \
			cp $$file docs/source/. ;\
	done
	jupytext demo/oscilator_v1.py -o docs/source/oscilator_v1.md
	jupytext demo/oscilator_v2.py -o docs/source/oscilator_v2.md
	jupytext demo/oscilator_v3.py -o docs/source/oscilator_v3.md
	jupytext demo/lorentz.py -o docs/source/lorentz.md
	jupytext demo/tentusscher.py -o docs/source/tentusscher.md
	jupytext demo/tentusscher_field_parameters.py -o docs/source/tentusscher_field_parameters.md
	cd docs && make html

show:
	open docs/build/html/index.html
