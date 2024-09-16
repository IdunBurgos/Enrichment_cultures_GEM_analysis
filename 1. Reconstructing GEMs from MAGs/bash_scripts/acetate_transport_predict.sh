#!/opt/homebrew/bin/bash

# gapfilling for defined media with SynCon1 and produced compounds from growing community members.

set -eu

for filename in ../output/MAGs_fasta/*.faa
do
    MAG=${filename:20:-4}
    makeblastdb -in ../output/MAGs_fasta/${MAG}.faa -dbtype prot
    
    blastp -query ./transporters/transporters.faa -db ../output/MAGs_fasta/${MAG}.faa -out ../output/transporters/${MAG}.tsv -outfmt 6
    
done