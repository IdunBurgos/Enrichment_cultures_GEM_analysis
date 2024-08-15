#!/usr/bin/env bash

# gapfilling for defined media with SynCon1 and produced compounds from growing community members.

set -eux

# Run with selected score


# For each MAG in directory
for filename in ../output/GEMs/GEMs_intermediate/GEMs_ACt2r/*.xml
do

    MAG=${filename:44:-4}"\t"
    
    
    # If MAG is among the 99% most abundant
    if grep -xq $MAG ../output/relevant_MAGs_99.txt
    then
    
        # Run carveme
        community_id=$(grep $MAG ../output/MAG2community_id.tsv | cut  -f2)  
        gapfill ../output/GEMs/GEMs_intermediate/GEMs_ACt2r/${filename:44:-4}.xml -o ../output/GEMs/GEMs_final2/${filename:44:-4}.xml -m ${community_id} --mediadb ../output/gapfill_media/gapfill_media.tsv --fbc2 --verbose 
    fi
done
