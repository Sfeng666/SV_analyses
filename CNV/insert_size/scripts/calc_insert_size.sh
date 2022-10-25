#!/bin/bash
# script to use picard CollectInsertSizeMetrics to calculate metrics and distribution of insert size

# 0. receive arguments
## 0.1. conda environment
path_conda=$1
env_conda=$2

## 0.2. arguments
path_tar=$3
name_sample=$4
name_sample_new=$5
path_out=$6

# 1. local variables
in_bam_original=$path_out/$name_sample\_sort_merge_dedup_indel.bam
in_bam_rename=$path_out/$name_sample_new\_sort_merge_dedup_indel.bam
out_metrics=$path_out/$name_sample_new\_insert_size_metrics.txt
out_hist=$path_out/$name_sample_new\_insert_size_histogram.pdf
log=$path_out/$name_sample_new\_calc_insert_size.log

# 1. configure the conda enviroment
source $path_conda $env_conda

# 1. extract realigned bam files from tar.zip files
tar -zxvf $path_tar/05realign_$name_sample\.tar.gz -C $path_out $name_sample\_sort_merge_dedup_indel.bam > $log 2>&1
mv $in_bam_original $in_bam_rename

# 2. run picard CollectInsertSizeMetrics
picard CollectInsertSizeMetrics \
      I=$in_bam_rename \
      O=$out_metrics \
      H=$out_hist \
      >> $log 2>&1
# rm $in_bam_rename