#!/usr/bin/env bash

# Automatically chooses the bacterial model as universe 
# Do not add media at this point

set -eux

# Run with selected score
soft_score=$1

# For each MAG in directory
for filename in ./output/MAGs_fasta/*.faa
do

    MAG=${filename:20:-4}"\t"
    
    if [ ${filename:20:-4} != "CH15-bin.15" ]; then
       
       # If MAG is among the 99% most abundant
        if grep -xq $MAG ./output/relevant_MAGs_99.txt
        then

            # Run carveme
            sub_sour=$(grep $MAG ./output/MAG2community_id.tsv | cut  -f2)  
            carve ./output/MAGs_fasta/${filename:20:-4}.faa --fbc2 -o ./output/GEMs/GEMs_soft_constraints_score_${soft_score}/${filename:20:-4}.xml --verbose --gapfill 'LB[-O2]' --soft ./output/soft_constraints/SC_${sub_sour}.tsv --soft-score $soft_scoremedia_dfs

            # Move tsv score file
            if [ -f ./output/MAGs_fasta/${filename:20:-4}.tsv ]
            then
                mv ./output/MAGs_fasta/${filename:20:-4}.tsv ./output/GEMs/GEMs_soft_constraints_score_${soft_score}/${filename:20:-4}.tsv
            fi

        fi
    fi
done
