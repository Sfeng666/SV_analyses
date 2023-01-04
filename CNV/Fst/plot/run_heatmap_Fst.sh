#!/bin/bash

### 1. configure the conda enviroment ###
# source /opt/miniconda3/bin/activate WGS_analysis

### 2. set parameters ###
# rscript=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/Fst/plot/heatmap_Fst.R
# wd=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/Fst/plot
# dd=/home/siyuan/jobs/suzukii_WGS/SV_analyses/CNV/Fst/results
rscript=/Users/siyuansmac/bioinfo/project/suzukii_WGS/SV_analyses/CNV/Fst/plot/heatmap_Fst.R
wd=/Users/siyuansmac/bioinfo/project/suzukii_WGS/SV_analyses/CNV/Fst/plot
dd=/Users/siyuansmac/bioinfo/project/suzukii_WGS/SV_analyses/CNV/Fst/results


### 3. run the plotting R script ###
Rscript $rscript $wd $dd dup >> run_heatmap_Fst.log 2>&1
Rscript $rscript $wd $dd del >> run_heatmap_Fst.log 2>&1