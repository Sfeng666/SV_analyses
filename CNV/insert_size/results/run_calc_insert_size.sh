#!/bin/bash
# use picard CollectInsertSizeMetrics to calculate metrics and distribution of insert size

# 0. parameters
path_conda=/home/sfeng77/anaconda3/bin/activate
env_conda=WGS_analysis
path_tar=/staging/sfeng77/test_pipeline/out/
path_out=/home/sfeng77/jobs/SV_analyses/CNV/insert_size/results
script=/home/sfeng77/jobs/SV_analyses/CNV/insert_size/scripts/calc_insert_size.sh
rename_dict=/home/sfeng77/jobs/rename_samples/rename_dict.txt

# 1. run calculcation of insert size
# rename_dict=/Users/siyuansmac/bioinfo/project/suzukii_WGS/rename_samples/rename_dict.txt
while IFS= read -r old_name new_name <&3; do
  {
    bash $script \
    $path_conda $env_conda $path_tar $old_name $new_name $path_out &
  } 3<&-
done 3< "$rename_dict"