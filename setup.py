import sys
import sysconfig

from skbuild import setup

setup(
    name="pygoss",
    python_requires=">=3.7.0",
    version="0.1.0",
    description="Basix Python interface",
    url="https://github.com/ComputationalPhysiology/goss",
    author="Henrik Finsberg, Johan Hake, Cécile Daversin-Catty",
    author_email="henriknf@simula.no",
    maintainer_email="henriknf@simula.no",
    license="GNU General Public License v3 (GPLv3)",
    keywords=["ODE", "solver", "system", "equations", "cuda"],
    install_requires=[
        "numpy",
        "modelparameters",
        "gotran",
        "typer",
        "cppyy",
    ],
    extras_require={
        "test": ["pytest", "pytest-cov"],
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
