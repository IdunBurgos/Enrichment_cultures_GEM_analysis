#!/usr/bin/env bash

# Make GEMs for the three species:
# 
# GCF_902809985.1 - BL-3
# GCF_902809935.1 - BL-4
# GCF_000022065.1 - R. cellulolyticum


set -eux

filenames=("GCF_902809985.1" "GCF_902809935.1" "GCF_000022065.1")

for filename in ${filenames[@]}; do
    echo ${filename}
    carve --refseq ${filename} --fbc2 -o ./output/GEMs_test/${filename}.xml --verbose --gapfill 'LB[-O2]' -d
done
