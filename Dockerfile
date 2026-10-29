FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    git \
    bzip2 \
    jq \
    tar \
    unzip \
    build-essential \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    software-properties-common \
    dirmngr \
    gpg \
    perl \
    jellyfish \
    default-jre \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN mamba create -y -n fungi python=3.9

RUN mamba install -y -n fungi -c conda-forge -c bioconda \
    fastqc=0.12.1 \
    fastp=0.23.4 \
    spades=3.15.5 \
    busco=5.5.0 \
    repeatmasker=4.1.5 \
    augustus=3.5.0 \
    gffread \
    orthofinder \
    mafft \
    trimal \
    raxml-ng \
    antismash \
    eggnog-mapper \
    interproscan \
    r-base \
    r-tidyverse \
    r-fsa \
    r-emmeans \
    r-mass \
    r-car

RUN find /opt/conda/envs/fungi -name "hmmlib_1.75" -exec /opt/conda/envs/fungi/bin/hmmpress -f {} \; || true

RUN /opt/conda/envs/fungi/bin/pip install --no-cache-dir torch torchvision --index-url https://download.pytorch.org/whl/cpu

RUN /opt/conda/envs/fungi/bin/pip install --no-cache-dir skorch numba "biopython<1.80" "numpy<2" tqdm joblib

RUN /opt/conda/envs/fungi/bin/pip install --no-cache-dir tensorflow==2.11.0 keras==2.11.0

RUN /opt/conda/envs/fungi/bin/pip install --no-cache-dir git+https://github.com/ibe-uw/tiara.git --no-deps

RUN git clone https://github.com/BolognaBiocomp/deepsig.git /opt/deepsig && \
    sed -i "s/python_requires='>=3.8,<3.9'/python_requires='>=3.8'/g" /opt/deepsig/setup.py 2>/dev/null || true && \
    /opt/conda/envs/fungi/bin/pip install --no-cache-dir /opt/deepsig && \
    sed -i "s/from keras.preprocessing.sequence import pad_sequences/from tensorflow.keras.preprocessing.sequence import pad_sequences/g" /opt/conda/envs/fungi/lib/python3.9/site-packages/deepsig/helpers.py 2>/dev/null || true

ENV TF_CPP_MIN_LOG_LEVEL=3 \
    DEEPSIG_ROOT=/opt/deepsig \
    MPLCONFIGDIR=/tmp/matplotlib \
    PATH="/opt/conda/envs/fungi/bin:$PATH"

RUN git clone https://github.com/JanaSperschneider/EffectorP-3.0.git /opt/EffectorP-3.0 && \
    cd /opt/EffectorP-3.0 && \
    mkdir -p weka-3-8-4 && \
    wget -O weka-3-8-4/weka.jar "https://repo1.maven.org/maven2/nz/ac/waikato/cms/weka/weka-stable/3.8.4/weka-stable-3.8.4.jar" && \
    chmod +x EffectorP.py

ENV PATH="/opt/EffectorP-3.0:$PATH"

RUN mamba clean -a -y

RUN chmod -R 777 /opt/conda/envs/fungi/lib/python3.9/site-packages/antismash || true

RUN find /opt/conda/envs/fungi -name "hmmlib_1.75" -exec hmmpress -f {} \; || true

RUN /opt/conda/envs/fungi/bin/pip install --no-cache-dir phykit

RUN /opt/conda/envs/fungi/bin/pip install --no-cache-dir "biopython<1.80" "numpy<2" --force-reinstall

WORKDIR /data
