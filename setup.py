import sys
import sysconfig
from pathlib import Path

from skbuild import setup

here = Path(__file__).parent
long_description = (here / "README.md").read_text()

setup(
    name="pygoss",
    python_requires=">=3.7.0",
    version="0.2.1",
    description="Python interface to goss - General ODE System Solver",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ComputationalPhysiology/goss",
    author="Henrik Finsberg, Johan Hake, Cécile Daversin-Catty",
    author_email="henriknf@simula.no",
    maintainer_email="henriknf@simula.no",
    license="LGPLv3+",
    keywords=["ODE", "solver", "system", "equations", "cuda"],
    install_requires=[
        "numpy",
        "modelparameters",
        "gotran",
        "rich-click",
        "pydantic",
        "cppyy",
        "importlib-metadata;python_version<'3.8'",
    ],
    extras_require={
        "test": ["pytest", "pytest-cov"],
        "cbcbeat": ["cbcbeat"],
        "plot": ["matplotlib"],
        "dev": [
            "Sphinx",
            "black",
            "bump2version",
            "flake8",
            "ipython",
            "isort",
            "mypy",
            "pdbpp",
            "pip",
            "pre-commit",
            "pytest",
            "pytest-cov",
            "twine",
            "wheel",
        ],
        "docs": [
            "Sphinx",
            "myst_parser",
            "sphinx_book_theme",
            "sphinxcontrib-bibtex",
            "sphinx-math-dollar",
        ],
    },
    entry_points={"console_scripts": ["goss=goss.cli:app"]},
    cmake_args=[
        "-DPython3_EXECUTABLE=" + sys.executable,
        "-DPython3_LIBRARIES=" + sysconfig.get_config_var("LIBDEST"),
        "-DPython3_INCLUDE_DIRS=" + sysconfig.get_config_var("INCLUDEPY"),
    ],
    packages=["goss"],
    package_dir={"": "python"},
    cmake_install_dir="python/goss/",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
