# Installation


## Installing dependencies
First, you need to install `boost`, i.e

(Debian)
```
apt-get install libboost-all-dev
```
(Fedora)
```
yum install boost-devel
```
Mac OSX
```
brew install boost
```

## Installing goss

Next you can install the python package.
```
python -m pip install pygoss
```

## Install with conda

TODO: We should make it possible to install `goss` with conda. Preferable, we should add it to conda-forge.

## Install from source

If you want the latest version or you want to develop `goss` you can install the code on the `main` branch

```
python -m pip install git+https://github.com/ComputationalPhysiology/goss.git@main
```
or clone the repository and install it from there

```
git clone git@github.com:ComputationalPhysiology/goss.git
cd goss
python -m pip install .
```

(section:docker-install)=
## Docker

`goss` is also available through [Docker](https://docs.docker.com/get-docker/). This is a good choice if you want to use `goss` in an isolated environment.

We provide both a pre-built docker image which you can get by pulling from docker hub
```
docker pull ghcr.io/computationalphysiology/goss:latest
```

### Building your own docker image

An alternative to pulling the image from docker hub, is to build it yourselves.
We provide a Dockerfile in the root of the repo that contain all the instructions for building the docker image. You can do this by executing the following command in the root folder of the project

```
docker build -t goss .
```
This will create a docker image with the name `goss`.


## Development installation

Developers should use editable install and install the development requirements using the following command
```
python -m pip install -e ".[dev]"
```
It is also recommended to install the `pre-commit` hook that comes with the package
```
pre-commit install
```
Note that linters and formatters will run in the CI system.
