# Use Python 3.8 base image
FROM python:3.8

# Copy files to container
COPY . /app
WORKDIR /app

# Install build dependencies
RUN apt-get update
RUN apt-get install -y python3-dev build-essential gfortran cmake ibopenblas-* libatlas-* liblapack-* parallel python3-pip ipython3 tk python3-tk

# Install Python dependencies
WORKDIR /app
RUN pip3 install -r requirements.txt

# Install library
WORKDIR /app
RUN python setup.py install
