#!/bin/bash
# Exit on error
set -e

# Copy your scripts to $PREFIX/bin so they are available in PATH
mkdir -p $PREFIX/bin
cp $SRC_DIR/split-orf-prediction/run_splitorfs_pipeline.sh $PREFIX/bin/split-orf-prediction
cp -r $SRC_DIR/split-orf-prediction/Genomic_scripts_18_10_24 $PREFIX/bin/Genomic_scripts_18_10_24
cp -r $SRC_DIR/split-orf-prediction/SplitOrfs-master $PREFIX/bin/SplitOrfs-master
cp -r $SRC_DIR/split-orf-prediction/Uniqueness_scripts $PREFIX/bin/Uniqueness_scripts
cp $SRC_DIR/split-orf-prediction/*.Rmd $PREFIX/bin/
chmod +x $PREFIX/bin/split-orf-prediction