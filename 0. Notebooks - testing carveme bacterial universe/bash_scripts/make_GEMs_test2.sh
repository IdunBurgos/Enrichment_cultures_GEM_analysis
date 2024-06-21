#!/usr/bin/env bash

# Make GEMs for the three species: 
# 
# GCF_902809985.1 - BL-3
# GCF_902809935.1 - BL-4
# GCF_000022065.1 - R. cellulolyticum

# Had to use model files because the ncbi reference stopped working.


set -eux

if [ -f ./input/curated_models_fasta/GCF_902809985.1.tsv ]; then
    rm ./input/curated_models_fasta/GCF_902809985.1.tsv
fi

if [ -f ./input/curated_models_fasta/GCF_902809935.1.tsv ]; then
    rm ./input/curated_models_fasta/GCF_902809935.1.tsv
fi

if [ -f ./input/curated_models_fasta/GCF_000022065.1.tsv ]; then
    rm ./input/curated_models_fasta/GCF_000022065.1.tsv
fi


carve ./input/curated_models_fasta/GCF_902809985.1.faa --fbc2 -o ./output/GEMs_test/clostridiaBL3_test.xml --verbose --gapfill 'LB[-O2]' -d
carve ./input/curated_models_fasta/GCF_902809935.1.faa --fbc2 -o ./output/GEMs_test/clostridiaBL4_test.xml --verbose --gapfill 'LB[-O2]' -d
carve ./input/curated_models_fasta/GCF_000022065.1.faa --fbc2 -o ./output/GEMs_test/RCell_H10_test.xml --verbose --gapfill 'LB[-O2]' -d

