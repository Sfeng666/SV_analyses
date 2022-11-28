#!/bin/bash
# script to generate bam file for each contig for each sample

# 0. variables & environment
sample=$1
contig=$2
conda_pack_squid=$3
path_in_stage=$4

## 0.2. configure the conda enviroment
set -e

workdir_pac=$(basename $conda_pack_squid)
ENVNAME=$(basename $conda_pack_squid .tar.gz)
ENVDIR=$ENVNAME

### these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $workdir_pac -C $ENVDIR
. $ENVDIR/bin/activate
rm $workdir_pac

## 0.3. set local paths
### 0.3.1. input files
in_bam_rename_sort=$sample\_sort.bam
in_bam_rename_sort_stage=$path_in_stage/$in_bam_rename_sort
cp $in_bam_rename_sort_stage ./

### 0.3.2. output & intermediate files
bam_contig=$sample\_$contig\.bam

# 1. extract the given contig as a bam file
samtools view -b $in_bam_rename_sort $contig > $bam_contig

# 2. transfer output contig to /staging
mv $bam_contig $path_in_stage

# 3. clear remained files
rm -rf *

