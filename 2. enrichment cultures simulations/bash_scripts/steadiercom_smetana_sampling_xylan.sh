# Running smetana for communities in different SynCon medias and on different carbon substrates

set -eux

# SC1_X communities
steadiercom ../output/GEMs/GEMs_final/*.xml -c ../output/steadiercom_sample_0.1.3/communities_99_SC1_X.tsv -m SC1_X --mediadb ../output/steadiercom_sample_0.1.3/media_db_constrained.tsv -o ../output/steadiercom_sample_0.1.3/results_1000/results_99_SC1_X --sample 1000 --growth 0.0089 --unlimited ../output/steadiercom_sample_0.1.3/unlimited.txt