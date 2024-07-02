#!/usr/bin/env bash

# Automatically chooses the bacterial model as universe 
# Do not add media at this point

set -eux

# Run with selected score
soft_score=$1


mkdir ../output/GEMs/GEMs_intermediate/GEMs_soft_constraints_score_${soft_score}/

# For each MAG in directory
for filename in ../output/MAGs_fasta/*.faa
do

    MAG=${filename:21:-4}"\t"
    
    if [ ${filename:21:-4} != "CH15-bin.15" ]; then
       
       # If MAG is among the 99% most abundant
        if grep -xq $MAG ../output/relevant_MAGs_99.txt
        then

            # Run carveme
            sub_sour=$(grep $MAG ../output/MAG2community_id.tsv | cut  -f2)  
            
            
            carve ../output/MAGs_fasta/${filename:21:-4}.faa --fbc2 -o ../output/GEMs/GEMs_intermediate/GEMs_soft_constraints_score_${soft_score}/${filename:21:-4}.xml --verbose --gapfill 'LB_extend' --mediadb ../output/soft_constraints/SC_media_db.tsv --soft ../output/soft_constraints/SC_${sub_sour}.tsv --soft-score $soft_score --solver gurobi

            # Move tsv score file
            if [ -f ../output/MAGs_fasta/${filename:21:-4}.tsv ]
            then
                mv ../output/MAGs_fasta/${filename:21:-4}.tsv ../output/GEMs/GEMs_intermediate/GEMs_soft_constraints_score_${soft_score}/${filename:21:-4}.tsv
            fi

        fi
    fi
done
