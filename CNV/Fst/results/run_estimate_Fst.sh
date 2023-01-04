#!/bin/bash

## parameters
wd=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/Fst/results
in_del_af=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/detect_CNV/results/global_CNV_AF/global_clusters_del_af.txt
in_dup_af=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/detect_CNV/results/global_CNV_AF/global_clusters_dup_af.txt
in_efs=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/GEA/data/genotype_data/median_efs.txt
in_chr_assign=/home/siyuan/jobs/suzukii_WGS/calc_maxcov/amb/assignment_cor_amb.txt
out_log=log.txt
script=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/Fst/scripts/estimate_Fst.py

## run the program
python $script \
--wd $wd \
--in_del_af $in_del_af \
--in_dup_af $in_dup_af \
--in_efs $in_efs \
--in_chr_assign $in_chr_assign \
--out_log $out_log \
> $out_log 2>&1