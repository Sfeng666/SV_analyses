#!/bin/bash

## parameters
wd=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/GEA/data/genotype_data
in_snp_dir=/home/siyuan/jobs/suzukii_WGS/call_snp/amb/updated_poolsnp/paratest/test_maxcov/call_snp_maxcov_99/mincov_12_mincount_1
in_del_af=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/detect_CNV/results/global_CNV_AF/global_clusters_del_af.txt
in_dup_af=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/detect_CNV/results/global_CNV_AF/global_clusters_dup_af.txt
in_sample_size=/home/siyuan/jobs/suzukii_WGS/rename_samples/cloned_input/sample_size.txt
in_chr_assign=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
in_sample_order=/home/siyuan/jobs/suzukii_WGS/EAA/data/env_data/sample_coord_UTM.txt
out_efs_sp=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/GEA/data/genotype_data/median_efs.txt
out_log=log.txt
script=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/GEA/scripts/prep_CNV_ct_split.py

## run the program
python $script \
--wd $wd \
--in_snp_dir $in_snp_dir \
--in_del_af $in_del_af \
--in_dup_af $in_dup_af \
--in_sample_size $in_sample_size \
--in_chr_assign $in_chr_assign \
--in_sample_order $in_sample_order \
--out_efs_sp $out_efs_sp \
--out_log $out_log \
> $out_log 2>&1