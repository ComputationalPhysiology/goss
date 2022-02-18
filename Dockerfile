FROM finsberg/fenics:latest

COPY . /app
WORKDIR /app

RUN python3 -m pip install .
