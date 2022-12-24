FROM ubuntu

RUN apt-get update

RUN apt-get install build-essential -y

RUN apt-get install git wget -y

RUN apt-get install gfortran -y

RUN mkdir -p /OpenBlas \
    && cd /OpenBlas && wget https://github.com/xianyi/OpenBLAS/releases/download/v0.3.21/OpenBLAS-0.3.21.tar.gz \
    && tar -xf OpenBLAS-0.3.21.tar.gz && cd OpenBLAS-0.3.21 && make -s FC=gfortran \
    && make PREFIX=/OpenBlas install && rm -rf /OpenBlas/OpenBLAS-0.3.21.tar.gz && rm -rf /OpenBlas/OpenBLAS-0.3.21 

RUN export "LD_LIBRARY_PATH=/OpenBlas/lib:$LD_LIBRARY_PATH" \
    && export "PATH=/OpenBlas/include:$PATH"