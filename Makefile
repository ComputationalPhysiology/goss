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

convert-demos:
	mkdir -p docs/demos
	jupytext demo/oscilator_v1.py -o docs/demos/oscilator_v1.md
	jupytext demo/oscilator_v2.py -o docs/demos/oscilator_v2.md
	jupytext demo/oscilator_v3.py -o docs/demos/oscilator_v3.md
	jupytext demo/lorentz.py -o docs/demos/lorentz.md
	jupytext demo/tentusscher.py -o docs/demos/tentusscher.md
	jupytext demo/tentusscher_field_parameters.py -o docs/demos/tentusscher_field_parameters.md
	jupytext demo/odes_on_mesh.py -o docs/demos/odes_on_mesh.md
	jupytext demo/bidomain.py -o docs/demos/bidomain.md
	jupytext demo/monodomain.py -o docs/demos/monodomain.md
	jupytext demo/niederer_benchmark.py -o docs/demos/niederer_benchmark.md
	jupytext demo/multicellmodel.py -o docs/demos/multicellmodel.md


docs: convert-demos ## generate Sphinx HTML documentation, including API doc
	jupyter-book build .

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
