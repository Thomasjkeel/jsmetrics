FROM continuumio/miniconda3

# make working directory for application
# WORKDIR /

# Create the environment and install dependencies
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "jsmet", "/bin/bash", "-c"]

# activate conda environment and check installs
RUN echo "Make sure xarray is installed:"
RUN python -c "import xarray"

# Copy source code
COPY . .
# COPY experiments/ .

# The code to run when container is started -> will allow running of experiements
# using docker run image_name experiments/experiments.py
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "jsmet", "python"]