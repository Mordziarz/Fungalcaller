#!/bin/bash
set -e

# ==========================================
# PATHS AND IMAGE CONFIGURATION
# ==========================================
HOST_DB_DIR="./databases"
IMAGE_NAME="mordziarz/fungi_pipeline:latest"

echo "=========================================="
echo " CHECKING DOCKER ENVIRONMENT"
echo "=========================================="

if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo "Docker image $IMAGE_NAME is not available locally. Building from Dockerfile..."
    if [ -f "Dockerfile" ]; then
        docker build -t "$IMAGE_NAME" .
    else
        echo "ERROR: Image not found and no Dockerfile in the current directory!"
        exit 1
    fi
fi

# ==========================================
# CREATING DIRECTORY STRUCTURE FOR DATABASES
# ==========================================
echo "Creating directory structure in $HOST_DB_DIR..."
mkdir -p "$HOST_DB_DIR/eggnog"
mkdir -p "$HOST_DB_DIR/antismash"
mkdir -p "$HOST_DB_DIR/interproscan"

ABS_DB_DIR=$(realpath "$HOST_DB_DIR")

# ==========================================
# 1. EGGNOG-MAPPER DATABASE (~40-50 GB)
# ==========================================
echo "=========================================="
echo " [1/3] Downloading eggNOG-mapper database..."
echo "=========================================="
docker run --rm -v "$ABS_DB_DIR/eggnog:/data_target" "$IMAGE_NAME" bash -c "
    SCRIPT_PATH=\$(which download_eggnog_data.py) && \
    sed -i 's/eggnogdb.embl.de/eggnog5.embl.de/g' \$SCRIPT_PATH && \
    download_eggnog_data.py -y --data_dir /data_target
"


# ==========================================
# 2. ANTISMASH DATABASES (~10-20 GB)
# ==========================================
echo "=========================================="
echo " [2/3] Downloading antiSMASH databases..."
echo "=========================================="
docker run --rm -v "$ABS_DB_DIR/antismash:/data_target" "$IMAGE_NAME" \
    download-antismash-databases --database-dir /data_target


# ==========================================
# 3. INTERPROSCAN DATABASES (~40-50 GB)
# ==========================================
echo "=========================================="
echo " [3/3] Downloading and initializing InterProScan databases..."
echo "=========================================="
docker run --rm -it -v "$ABS_DB_DIR/interproscan:/data_target" "$IMAGE_NAME" bash -c "
    cd /data_target && \
    if [ ! -d "interproscan-5.70-102.0" ]; then
        echo 'Downloading InterProScan archive...' && \
        wget -c https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.70-102.0/interproscan-5.70-102.0-64-bit.tar.gz && \
        tar -pxvzf interproscan-5.70-102.0-64-bit.tar.gz && \
        rm interproscan-5.70-102.0-64-bit.tar.gz
    fi && \
    echo 'Running InterProScan setup and indexing...' && \
    cd interproscan-5.70-102.0 && \
    python3 setup.py -f interproscan.properties && \
    echo 'Running hmmpress on SuperFamily database...' && \
    hmmpress -f data/superfamily/1.75/hmmlib_1.75
"

echo "=========================================="
echo " SUCCESS! All databases have been successfully downloaded and configured in $HOST_DB_DIR"
echo "=========================================="
