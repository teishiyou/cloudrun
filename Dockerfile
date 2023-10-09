# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM bioconductor/bioconductor_docker:devel
#FROM screen:latest

# Allow statements and log messages to immediately appear in the logs
ENV PYTHONUNBUFFERED True

# Install required Bioconductor package
RUN R -e 'BiocManager::install("readr")'
RUN R -e 'BiocManager::install("readxl")'
RUN R -e 'BiocManager::install("openxlsx")'
RUN R -e 'BiocManager::install("dplyr")'
RUN R -e 'BiocManager::install("tidyr")'
RUN R -e 'BiocManager::install("purrr")'
RUN R -e 'BiocManager::install("stringr")'
RUN R -e 'BiocManager::install("sangeranalyseR")'
RUN R -e 'BiocManager::install("yaml")'

# Copy local code to the container image.
ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

# Install production dependencies.
RUN pip install --no-cache-dir -r requirements.txt

RUN apt-get install git
RUN git clone https://**@github.com/revolka/benchling.git /app/benchling

COPY cred.json /app/

RUN apt-get update && apt-get install -y curl apt-transport-https ca-certificates gnupg

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

RUN apt-get update && apt-get install -y google-cloud-sdk

#gcsguse install pinding
#RUN wget "https://github.com/GoogleCloudPlatform/gcsfuse/releases/download/v1.1.0/gcsfuse_1.1.0_amd64.deb" -P /app/
RUN sudo apt-get update
#RUN sudo apt-get install -fy "/app/gcsfuse_1.1.0_amd64.deb"

#mount gcs folder
#RUN gcloud auth login --cred-file=/app/cred.json
#RUN mkdir /mnt/benchling
#RUN gcsfuse --foreground benchling /mnt/benchling


# Run the web service on container startup. Here we use the gunicorn
# webserver, with one worker process and 8 threads.
# For environments with multiple CPU cores, increase the number of workers
# to be equal to the cores available.
# Timeout is set to 0 to disable the timeouts of the workers to allow Cloud Run to handle instance scaling.
CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 main:app