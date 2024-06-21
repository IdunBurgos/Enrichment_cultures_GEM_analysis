#!/usr/bin/env bash

# To move tsv files created during reconstruction into the the right directory

set -eux

for filename in ./output/MAGs_fasta/*.faa
do
    echo ${filename:20:-4}
    if [ -f ./output/MAGs_fasta/${filename:20:-4}.tsv ]; then
        mv ./output/MAGs_fasta/${filename:20:-4}.tsv ./output/GEMs/GEMs_soft_constraints_score_0_5/${filename:20:-4}.tsv
    fi

done
