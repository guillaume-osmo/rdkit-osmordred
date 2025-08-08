FROM continuumio/miniconda3:latest

ENV boost_version=1.82
ENV python="python=3.11"


RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-dev \
    libboost-all-dev \
    libboost-python-dev \
    liblapack-dev \
    liblapacke-dev \
    libgl1-mesa-dev \
    gcc-multilib \
    libglib2.0-dev \
    && rm -rf /var/lib/apt/lists/*

RUN conda create --yes --name rdkit_build $python && \
    conda install --yes --name rdkit_build \
    libcxx cmake \
    libboost=$boost_version libboost-devel=$boost_version \
    libboost-python=$boost_version libboost-python-devel=$boost_version \
    qt \
    numpy matplotlib=3.8 cairo pillow eigen pandas=2.1  \
    jupyter=1.0 ipython=8.20 sphinx myst-parser pytest nbval openblas \
    liblapack liblapacke \
    -c conda-forge && \
    conda clean -afy


    
WORKDIR /work