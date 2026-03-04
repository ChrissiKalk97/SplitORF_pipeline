#!/bin/bash
# Exit on error
set -e

# Copy your scripts to $PREFIX/bin so they are available in PATH
mkdir -p $PREFIX/bin
cp $SRC_DIR/split-orf-prediction/run_splitorfs_pipeline.sh $PREFIX/bin/split-orf-prediction
chmod +x $PREFIX/bin/split-orf-prediction