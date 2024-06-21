#!/usr/bin/env bash

# Running smetana for communities in different SynCon medias and on different carbon substrates

set -eux


smetana ./output/smetana/GEMs_adapt_99_gapfill_media/*.xml -c ./output/smetana/communities_top95.tsv -m SC1_X,SC1_C,SC2_X,SC2_C --mediadb ./output/smetana/media_db.tsv --flavor fbc2 -o ./output/smetana/smetana_top_95_gapfill_media/coupling -v --molweight -d


smetana ./output/smetana/GEMs_adapt_99_gapfill_media/*.xml -c ./output/smetana/communities_top95.tsv -m SC1_X,SC1_C,SC2_X,SC2_C --mediadb ./output/smetana/media_db.tsv --flavor fbc2 -o ./output/smetana/smetana_top_95_gapfill_media/no_coupling -v --molweight -d --no-coupling


smetana ./output/smetana/GEMs_adapt_99_gapfill_media/*.xml -c ./output/smetana/communities_top99.tsv -m SC2_X,SC2_C --mediadb ./output/smetana/media_db.tsv --flavor fbc2 -o ./output/smetana/smetana_top_99_gapfill_media/coupling -v --molweight -d


smetana ./output/smetana/GEMs_adapt_99_gapfill_media/*.xml -c ./output/smetana/communities_top99.tsv -m SC2_X,SC2_C --mediadb ./output/smetana/media_db.tsv --flavor fbc2 -o ./output/smetana/smetana_top_99_gapfill_media/no_coupling -v --molweight -d --no-coupling


