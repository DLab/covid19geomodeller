# Use Python 3.8 base image
FROM python:3.8

# Copy files to container
COPY . /app
WORKDIR /app

# Install build dependencies
RUN apt-get update
RUN apt-get install -y python3-dev build-essential gfortran cmake ibopenblas-* libatlas-* liblapack-* parallel

# Install Sundials 5.1.0 (as a dependency to install scikits.odes)
ARG SUN_VERS=5.1.0
ARG SUN_DOWN=https://computation.llnl.gov/projects/sundials/download
RUN wget -q ${SUN_DOWN}/sundials-${SUN_VERS}.tar.gz && tar -xzf sundials-${SUN_VERS}.tar.gz && cd sundials-${SUN_VERS}
WORKDIR /app/sundials-${SUN_VERS}
RUN mkdir build && cd build && cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true -DBUILD_STATIC_LIBS:BOOL=true -DBUILD_SHARED_LIBS:BOOL=true -DEXAMPLES_ENABLE:BOOL=true -DEXAMPLES_INSTALL:BOOL=true -DMPI_ENABLE:BOOL=true -DOPENMP_ENABLE:BOOL=true -DPTHREAD_ENABLE:BOOL=true -DCUDA_ENABLE:BOOL=false -DRAJA_ENABLE:BOOL=false -DBLAS_ENABLE:BOOL=false -DLAPACK_ENABLE:BOOL=false && make -j$(nproc) && make install && cd ../../  && rm -rf sundials-${SUN_VERS}*
RUN cd /usr/local/examples && parallel make -C ::: $(find . -type f -name Makefile | xargs dirname)  && cd ../

# Install Python dependencies
WORKDIR /app
RUN pip3 install -r requirements.txt

# Run example to test the dependencies were installed successfully
WORKDIR /app/SEIRHVD
RUN python example.py
