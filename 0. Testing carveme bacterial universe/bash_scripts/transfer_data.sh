#!/usr/bin/env bash

# To transfer data from universal_model_extension to carveme directory

set -eux

if [ -f ./carveme/carveme/data/generated/bigg_gprs.csv.gz ]; then
    rm ./carveme/carveme/data/generated/bigg_gprs.csv.gz
fi

if [ -f ./carveme/carveme/data/generated/bigg_proteins.faa ]; then
    rm ./carveme/carveme/data/generated/bigg_proteins.faa
fi

if [ -f ./carveme/carveme/data/generated/bigg_proteins.dmnd ]; then
    rm ./carveme/carveme/data/generated/bigg_proteins.dmnd
fi


if [ -f ./carveme/carveme/data/generated/gene_annotations.tsv.gz ]; then
    rm ./carveme/carveme/data/generated/gene_annotations.tsv.gz
fi

if [ -f ./carveme/carveme/data/generated/gene_annotations.tsv.gz ]; then
    rm ./carveme/carveme/data/generated/model_specific_data.csv.gz
fi

if [ -f ./carveme/carveme/data/generated/universe_bacteria.xml.gz ]; then
    rm ./carveme/carveme/data/generated/universe_bacteria.xml.gz
fi



gzip -c ~/universal_model_extension/output/bigg_gprs.csv > ./carveme/carveme/data/generated/bigg_gprs.csv.gz
gzip -c ~/universal_model_extension/output/gene_annotations.tsv > ./carveme/carveme/data/generated/gene_annotations.tsv.gz
gzip -c ~/universal_model_extension/output/model_specific_data.csv > ./carveme/carveme/data/generated/model_specific_data.csv.gz
gzip -c ~/universal_model_extension/output/universe_bacteria.xml > ./carveme/carveme/data/generated/universe_bacteria.xml.gz

cp ~/universal_model_extension/output/bigg_proteins.faa ./carveme/carveme/data/generated
