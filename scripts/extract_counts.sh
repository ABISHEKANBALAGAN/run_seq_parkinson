#!/bin/bash
set -e

input_file=$1
output_file=$2

mkdir -p ../results

# Extract counts, remove 'Unassigned' etc.
zcat data/Raw_counts.tsv.gz \
  | awk 'NR==1 || $1 !~ /Unassigned|no_feature|ambiguous/' \
  > ../results/counts_matrix.tsv

# Extract coldata (metadata)
paste -d',' \
  <(zcat data/GSE294029_series_matrix.txt.gz | grep '!Sample_geo_accession' | cut -f2- | tr '\t' '\n' | tr -d '"') \
  <(zcat data/GSE294029_series_matrix.txt.gz | grep '!Sample_title' | cut -f2- | tr '\t' '\n' | tr -d '"') \
  | awk -F',' 'BEGIN{OFS=","; print "sample_id,GSM_id,condition,timepoint,replicate"}
  {
    # Titles may need custom parsing based on your dataset
    sample=$2; cond="HC"; tp="d60"; rep="1"; # Placeholder: edit accordingly
    print sample,$1,cond,tp,rep
  }' \
  > ../results/coldata.csv
