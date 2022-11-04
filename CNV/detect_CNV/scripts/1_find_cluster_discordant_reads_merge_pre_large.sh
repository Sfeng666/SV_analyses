#!/bin/bash
# script to prepare bam file of each contig across samples for discordant read-pair detection and clustering

# 0. set up arguments and conda environment
## 0.1. conda environment
path_conda=$1
env_conda=$2

## 0.2. arguments
path_tar=$3
rename_dict=$4
path_out=$5
path_out_sp=$6
contig=$7
path_in_stage=$8

log=$path_out_sp/log_1_pre_$contig
# 0.3. configure the conda enviroment
source $path_conda $env_conda

echo -e "# Whole job starts at " `date +%F'  '%H:%M` "\n" > $log

while read -r name_sample name_sample_new; do
    name_sample=$(echo $name_sample | sed 's/\r$//') # get rid of weird delimiters from operation system
    name_sample_new=$(echo $name_sample_new | sed 's/\r$//') # get rid of weird delimiters from operation system
    in_bam_original=$path_out/$name_sample_new/$name_sample\_sort_merge_dedup_indel.bam
    in_bam_rename=$path_out/$name_sample_new/$name_sample_new\.bam
    in_bam_rename_sort=$path_out/$name_sample_new/$name_sample_new\_sort.bam
    bam_contig=$path_out/$name_sample_new/$name_sample_new\_$contig\.bam
    # 1. extract realigned bam files from tar.zip files if it has not been extracted
    if ! [ -f $in_bam_original ] && ! [ -f $in_bam_rename ] && ! [ -f $in_bam_rename_sort ];  
    then
        tar -zxvf $path_tar/05realign_$name_sample\.tar.gz -C $path_out/$name_sample_new $(basename $in_bam_original) >> $log 2>&1
        mv $in_bam_original $in_bam_rename
        samtools sort -o $in_bam_rename_sort $in_bam_rename >> $log 2>&1 # sort the bam so that an index could be built
        rm $in_bam_rename # remove unsorted version to save space
        samtools index $in_bam_rename_sort >> $log 2>&1 # index the bam so that chromosomes could be extracted

    fi
    # 2. extract the given contig as a bam file if it has not been extracted
    if ! [ -f $bam_contig ];
    then
        while [ ! -f $in_bam_rename_sort\.bai ]; do sleep 1; done # wait until the index file is generated
        samtools view -b $in_bam_rename_sort $contig > $bam_contig >> $log 2>&1
    fi
done < $rename_dict

# 3. merge bam files of each contig across all samples
bam_contig_allsample=$path_out/*-*/*_$contig\.bam
bam_contig_merged=$path_out_sp/merged_$contig\.bam
if ! [ -f $bam_contig_merged ];
then
    samtools merge $bam_contig_merged $bam_contig_allsample
    rm $bam_contig_allsample
fi

# 4. further compress bam file by genozip
bam_contig_merged_genozip=$bam_contig_merged\.genozip
bam_contig_merged_genozip_tar=$path_in_stage/$(basename $bam_contig_merged_genozip)
if ! [ -f $bam_contig_merged_genozip ] && ! [ -f $bam_contig_merged_genozip_tar ];
then
    genozip --no-test -^ --vblock 150 $bam_contig_merged >> $log 2>&1
    rm $bam_contig_merged
    # 5. mv the large zipped file to staging for uploading
    mv $bam_contig_merged_genozip $bam_contig_merged_genozip_tar >> $log 2>&1
fi
if [ -f $bam_contig_merged ];
then
    rm $bam_contig_merged
fi

# 5. mv the large zipped file to staging for uploading
mv $bam_contig_merged_genozip $bam_contig_merged_genozip_tar >> $log 2>&1

echo -e "# Whole job ends at " `date +%F'  '%H:%M` "\n" >> $log