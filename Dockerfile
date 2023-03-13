FROM ghcr.io/scientificcomputing/fenics:2023-03-01

COPY . /app
WORKDIR /app

RUN python3 -m pip install .
