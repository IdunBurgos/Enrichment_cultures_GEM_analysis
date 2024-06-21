#!/usr/bin/env bash

# Running smetana for communities in different SynCon medias and on different carbon substrates

set -eux


smetana ./output/smetana/GEMs_adapt_95/*.xml -c ./output/smetana/communities.tsv -m SC1_X,SC1_C,SC2_X,SC2_C --mediadb ./output/smetana/media_db.tsv --flavor fbc2 -o ./output/smetana/results/ -v --molweight -d --no-coupling


