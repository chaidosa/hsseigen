FROM ubuntu

RUN apt-get update

RUN apt-get install gcc g++ make -y

RUN apt-get install git wget gdb vim -y

# RUN apt-get install gfortran -y

# RUN mkdir -p /OpenBlas \
#     && cd /OpenBlas && wget https://github.com/xianyi/OpenBLAS/releases/download/v0.3.21/OpenBLAS-0.3.21.tar.gz \
#     && tar -xf OpenBLAS-0.3.21.tar.gz && cd OpenBLAS-0.3.21 && make -s FC=gfortran \
#     && make PREFIX=/OpenBlas install && rm -rf /OpenBlas/OpenBLAS-0.3.21.tar.gz && rm -rf /OpenBlas/OpenBLAS-0.3.21 

# RUN echo "export LD_LIBRARY_PATH=/OpenBlas/lib:$LD_LIBRARY_PATH" >> /etc/profile \
#     && echo  "export PATH=/OpenBlas/include:$PATH" >> /etc/profile

RUN apt-get install libblas-dev liblapack-dev liblapacke-dev libatlas-base-dev libopenblas-dev -y

ARG USERNAME=hsseigen
ARG USER_UID=1001
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    #
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    # && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# ********************************************************
# * Anything else you want to do like clean up goes here *
# ********************************************************

RUN chsh -s /bin/bash $USERNAME

# [Optional] Set the default user. Omit if you want to keep the default as root.
USER $USERNAME