#!/usr/bin/env bash

# Automatically chooses the bacterial model as universe 
# Do not add media at this point

set -eux


for filename in ./output/MAGs_fasta/*.faa
do

    MAG=${filename:20:-4}"\t"
    
    echo $PWD
    if grep -xq $MAG ./output/soft_constraints_v2/relevant_MAGs_0_1.txt
    then
        sub_sour=$(grep $MAG ./output/soft_constraints/MAG2sour_sub_id.tsv | cut  -f2)  
        carve ./output/MAGs_fasta/${filename:20:-4}.faa --fbc2 -o ./output/GEMs/GEMs_soft_constraints_score_0_1/${filename:20:-4}.xml --verbose --gapfill 'LB[-O2]' --soft ./output/soft_constraints_v2/SC_${sub_sour}.tsv --soft-score 0.1
    fi
done
