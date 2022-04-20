.PHONY: clean clean-test clean-pyc clean-build docs help build-cpp
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
	jupytext demo/bidomain.py -o docs/source/bidomain.md
	jupytext demo/monodomain.py -o docs/source/monodomain.md
	jupytext demo/niederer_benchmark.py -o docs/source/niederer_benchmark.md
	jupytext demo/multicellmodel.py -o docs/source/multicellmodel.md
	cd docs && make html

show:
	open docs/build/html/index.html

delete-build-cpp:
	rm -rf build-cpp

build-cpp-test: delete-build-cpp
	cmake -B build-cpp -S cpp -DBUILD_TESTS=ON
	cmake --build build-cpp

build-cpp: delete-build-cpp
	cmake -B build-cpp -S cpp -DBUILD_TESTS=OFF
	cmake --build build-cpp

test-cpp: build-cpp-test
	cd build-cpp && ctest - V && cd ..


slides:
	rm -f docs/source/presentation.ipynb
	jupytext docs/source/presentation.md -o docs/source/presentation.ipynb
	jupyter nbconvert docs/source/presentation.ipynb --to slides  --post serve --SlidesExporter.reveal_scroll=True
