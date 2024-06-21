#!/usr/bin/env bash

# Automatically chooses the bacterial model as universe 
# Do not add media at this point

set -eux

for filename in ./output/MAGs_fasta/*.faa
do
    echo ${filename:20:-4}
    if [ -f ./output/MAGs_fasta/${filename:20:-4}.tsv ]; then
        rm ./output/MAGs_fasta/${filename:20:-4}.tsv
    fi
    
    carve ./output/MAGs_fasta/${filename:20:-4}.faa --fbc2 -o ./output/GEMs/${filename:20:-4}.xml --verbose --gapfill 'LB[-O2]'

done
