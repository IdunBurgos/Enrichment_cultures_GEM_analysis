#!/usr/bin/env bash

# Automatically chooses the bacterial model as universe 
# Do not add media at this point

set -eux

for filename in ./output/MAGs_fasta/*.faa
do
    echo ${filename:20:-4}
    if [ -f ./output/MAGs_fasta/${filename:20:-4}.tsv ]; then
        mv ./output/MAGs_fasta/${filename:20:-4}.tsv ./output/GEMs/GEMs_no_constraints/${filename:20:-4}.tsv
    fi
    
    MAG=${filename:20:-4}"\t"
    sub_sour=$(grep $MAG ./output/soft_constraints/MAG2sour_sub_id.tsv | cut  -f2)  
    
    carve ./output/MAGs_fasta/${filename:20:-4}.faa --fbc2 -o ./output/GEMs/${filename:20:-4}.xml --verbose --gapfill 'LB[-O2]' --soft ./output/soft_constraints/SC_${sub_sour}.tsv


done
