#!/usr/bin/env bash

# Running smetana for communities in different SynCon medias and on different carbon substrates

set -eux

# SC1_C communities
steadiercom ../output/GEMs/GEMs_adapt_media_ACt2r/*.xml -c ../output/steadiercom_sample_0.1.3/communities_99_CM_X.tsv -m SC1_X --mediadb ../output/steadiercom_sample_0.1.3/media_db_CM_X_constrained.tsv -o ../output/steadiercom_sample_0.1.3/test/results_99_CM_X --sample 10 --unlimited ../output/steadiercom_sample_0.1.3/unlimited.txt
