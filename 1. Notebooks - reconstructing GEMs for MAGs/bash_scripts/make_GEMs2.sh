#!/usr/bin/env bash

# Automatically chooses the bacterial model as universe 
# Do not add media at this point

set -eux

for filename in ../output/MAGs_fasta/*.faa
do
    echo ${filename:21:-4}
    if [ -f ../output/MAGs_fasta/${filename:21:-4}.tsv ]; then
        rm ../output/MAGs_fasta/${filename:21:-4}.tsv
    fi
    
    if [ -f ../output/GEMs/GEMs_no_constraints/${filename:21:-4}.tsv ]; then
        rm ../output/GEMs/GEMs_no_constraints/${filename:21:-4}.tsv
    fi
    
    MAG=${filename:21:-4}"\t"
       
   # If MAG is among the 99% most abundant
    if grep -xq $MAG ../output/relevant_MAGs_99.txt
    then
    
        carve ../output/MAGs_fasta/${filename:21:-4}.faa --fbc2 -o ../output/GEMs/GEMs_no_constraints/${filename:21:-4}.xml --verbose --gapfill 'LB[-O2]' --solver gurobi
        mv ../output/MAGs_fasta/${filename:21:-4}.tsv ../output/GEMs/GEMs_no_constraints/${filename:21:-4}.tsv
    fi

    
done
