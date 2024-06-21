#!/usr/bin/env bash

# Running smetana for communities in different SynCon medias and on different carbon substrates

set -eux


steadiercom ./output/GEMs/GEMs_adapt_media_ACt2r/*.xml -c ./output/smetana/communities_top99.tsv -m SC2_X,SC2_C --mediadb ./output/smetana/media_db.tsv -o ./output/smetana_steadiercom/top99 --sample 50

