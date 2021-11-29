# Get official docker image fo julia 
FROM julia

# set username as `root`
USER root

# update package list

RUN apt-get update && \
    apt-get install -y software-properties-common apt-transport-https curl && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      build-essential \
      python3-pip \
      python3-h5py \
      python3-numpy \
      python3-scipy \
      python3-dev \
      python3-setuptools \
      && \
    apt-get autoremove --purge -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*

# Install `Plots.jl and precompile it`
#RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Plots;precompile");using Plots'

ARG NB_USER=vscode
ARG NB_UID=1000
RUN useradd -u $NB_UID -m $NB_USER && \
    echo 'triqs ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER $NB_USER
WORKDIR /home/$NB_USER
ENV PATH=/home/$NB_USER/.local/bin:$PATH

#RUN pip3 install -U numpy scipy
RUN pip3 install irbasis irbasis3 xprec

COPY install_julia_pkg.jl /home/$NB_USER/
RUN julia install_julia_pkg.jl
