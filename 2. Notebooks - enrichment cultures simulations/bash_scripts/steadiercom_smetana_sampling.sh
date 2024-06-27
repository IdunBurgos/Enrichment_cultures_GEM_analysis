#!/usr/bin/env bash

# Running smetana for communities in different SynCon medias and on different carbon substrates

set -eux

# SC1_X communities
steadiercom ../output/GEMs/GEMs_adapt_media_ACt2r/*.xml -c ../output/steadiercom_sample_0.1.3/communities_99_SC1_X.tsv -m SC1_X --mediadb ../output/steadiercom_sample_0.1.3/media_db_constrained.tsv -o ../output/steadiercom_sample_0.1.3/results/results_99_SC1_X --sample 100 --unlimited ../output/steadiercom_sample_0.1.3/unlimited.txt

# SC1_C communities
steadiercom ../output/GEMs/GEMs_adapt_media_ACt2r/*.xml -c ../output/steadiercom_sample_0.1.3/communities_99_SC1_C.tsv -m SC1_C --mediadb ../output/steadiercom_sample_0.1.3/media_db_constrained.tsv -o ../output/steadiercom_sample_0.1.3/results/results_99_SC1_C --sample 100 --growth 0.0089 --unlimited ../output/steadiercom_sample_0.1.3/unlimited.txt

# SC2_C communities
steadiercom ../output/GEMs/GEMs_adapt_media_ACt2r/*.xml -c ../output/steadiercom_sample_0.1.3/communities_99_SC2_C.tsv -m SC2_C --mediadb ../output/steadiercom_sample_0.1.3/media_db_constrained.tsv -o ../output/steadiercom_sample_0.1.3/results/results_99_SC2_C --sample 100 --growth 0.0089 --unlimited ../output/steadiercom_sample_0.1.3/unlimited.txt

