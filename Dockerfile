# FROM python:3.9.4-buster
FROM continuumio/miniconda3

# make working directory for application
WORKDIR /

# Create the environment and install dependencies
COPY environment.yml .
RUN conda env create -f environment.yml

# Initialize conda in bash config fiiles:
RUN conda init bash

# activate conda environment and check installs
RUN conda activate jsmet
RUN echo "Make sure xarray is installed:"
RUN python -c "import xarray"

# Copy source code
COPY jet_stream_metrics/ .

# run the application
CMD ["python", "main.py"]