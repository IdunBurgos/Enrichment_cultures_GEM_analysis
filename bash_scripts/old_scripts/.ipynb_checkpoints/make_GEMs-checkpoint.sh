#!/usr/bin/env bash

# Automatically chooses the bacterial model as universe 
# Do not add media at this point

set -eux

for filename in ./output/MAGs_fasta/*.faa
do
    
    if [ ! -f ./output/GEMs/${filename:20:-4}.xml ]; then
        if [[ ${filename:20:-4} != "CH15-bin.15" ]]; then
            echo ${filename:20:-4}
            carve ./output/MAGs_fasta/${filename:20:-4}.faa --fbc2 -o ./output/GEMs/${filename:20:-4}.xml --verbose --gapfill 'LB[-O2]'
        fi
    fi
done
