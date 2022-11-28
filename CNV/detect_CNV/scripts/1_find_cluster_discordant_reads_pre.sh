#!/bin/bash
# script to prepare bam file of each contig of each sample for discordant read-pair detection and clustering

# 0. set up arguments and conda environment
## 0.1. conda environment
path_conda=$1
env_conda=$2

## 0.2. arguments
path_tar=$3
name_sample=$4
name_sample_new=$5
path_out_sp=$6
contig=$7
path_in_stage=$8

log=$path_out_sp/log_1_pre_$contig

# 0.3. local variables
in_bam_original=$path_out_sp/$name_sample\_sort_merge_dedup_indel.bam
in_bam_rename=$path_out_sp/$name_sample_new\.bam
in_bam_rename_sort=$path_out_sp/$name_sample_new\_sort.bam
bam_contig=$path_out_sp/$name_sample_new\_$contig\.bam

# 0.4. configure the conda enviroment
source $path_conda $env_conda

echo -e "# Whole job starts at " `date +%F'  '%H:%M` "\n" > $log
# 1. extract realigned bam files from tar.zip files if it has not been extracted
if ! [ -f $in_bam_original ] && ! [ -f $in_bam_rename ] && ! [ -f $in_bam_rename_sort ]; 
then
    tar -zxvf $path_tar/05realign_$name_sample\.tar.gz -C $path_out_sp $(basename $in_bam_original)
    mv $in_bam_original $in_bam_rename
    samtools sort -o $in_bam_rename_sort $in_bam_rename # sort the bam so that an index could be built
    rm $in_bam_rename # remove unsorted version to save space
    samtools index $in_bam_rename_sort # index the bam so that chromosomes could be extracted
fi

# 2. extract the given contig as a bam file if it has not been extracted
if ! [ -f $bam_contig ];
then
    while [ ! -f $in_bam_rename_sort\.bai ]; do sleep 1; done # wait until the index file is generated
    samtools view -b $in_bam_rename_sort $contig > $bam_contig
fi

# 3. further compress bam file by genozip
if ! [ -f $bam_contig\.genozip ];
then
    genozip --no-test $bam_contig >> $log 2>&1
    rm $bam_contig
    # ## remove bam contig file only if a compressed merged bam file has been generated
    # bam_contig_merged=$path_out_sp/../merged/merged_$contig\.bam
    # bam_contig_merged_genozip=$bam_contig_merged\.genozip
    # bam_contig_merged_tar=$path_in_stage/$(basename $bam_contig_merged)
    # if [ -f $bam_contig_merged_genozip] || [ -f $bam_contig_merged_tar];
    # then
    #     rm $bam_contig
    # fi
fi

echo -e "# Whole job ends at " `date +%F'  '%H:%M` "\n" >> $log