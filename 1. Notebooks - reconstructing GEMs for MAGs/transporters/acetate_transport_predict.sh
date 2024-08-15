#!/usr/bin/env bash

# gapfilling for defined media with SynCon1 and produced compounds from growing community members.

set -eux

for filename in ./output/MAGs_fasta/*.faa
do
    MAG=${filename:20:-4}
    makeblastdb -in ./output/MAGs_fasta/${MAG}.faa -dbtype prot
    
    blastp -query ./transporters/transporters.faa -db ./output/MAGs_fasta/${MAG}.faa -out ./transporters/${MAG}.tsv -outfmt 6
    
done


makeblastdb -in ./output/CH14-bin.0.faa -dbtype prot
blastp -query ./1.\ Notebooks\ -\ reconstructing\ GEMs\ for\ MAGs/transporters/transporters.faa -db ./output/CH14-bin.0.faa -out ./output/CH14-bin.0.tsv -outfmt 6