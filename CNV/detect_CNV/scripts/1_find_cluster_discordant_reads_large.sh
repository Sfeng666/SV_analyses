#!/bin/bash
# script to perform discordant read-pair detection and clustering (step 1) within each sample/across all samples 

# 0. variables & environment
## 0.1. receive arguments from 01qc.sub

sample=$1
contig=$2
conda_pack_squid=$3
insert_size_cutoff=$4
insertSizeDiffCutoff=$5
job=$6
path_in_stage=$7

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
in_bam_compressed=$sample\_$contig\.bam.genozip
in_bam_compressed_stage=$path_in_stage/$in_bam_compressed
in_bam=$sample\_$contig\.bam
in_sam=$sample\_$contig\.sam

### 0.3.2. output & intermediate files
out_everted_inserts=$sample\_$contig\_everted_inserts.txt
out_distant_inserts=$sample\_$contig\_distant_inserts.txt
out_everted_inserts_clustered=$sample\_$contig\_everted_inserts_clustered.txt
out_distant_inserts_clustered=$sample\_$contig\_distant_inserts_clustered.txt

### 0.3.3. software parameters
insert_size_cutoff=$insert_size_cutoff
insertSizeDiffCutoff=$insertSizeDiffCutoff
minNumSupportingInserts=2
minLenCutoff=50

### 0.3.3. software paths
dir_poolDiffCNV=biosoft/poolDiffCNV

### 0.3.4. log file
log=log_$job

# 1. copy bam file from staging and unzip zipped bam files
cp $in_bam_compressed_stage ./
genounzip -^ $in_bam_compressed

# 2. transfer bam to sam
samtools view -h -o $in_sam $in_bam > $log 2>&1
rm $in_bam

# 3. In each population sample, find discordant read pairs suggestive of a CNV, and cluster discordant read pairs into candidate CNVs
## 3.1. for duplications
find_cluster_dup(){
echo -e "### detection of discordant read pairs for duplications starts at" `date +%F'  '%H:%M` "\n" >> $log
cat $in_sam | python $dir_poolDiffCNV/findEvertedInserts.py > $out_everted_inserts # find discordant read pairs for duplications
echo -e "### detection of discordant read pairs for duplications ends at" `date +%F'  '%H:%M` "\n" >> $log

echo -e "### clustering of discordant read pairs for duplications starts at" `date +%F'  '%H:%M` "\n" >> $log
cat $out_everted_inserts | python $dir_poolDiffCNV/clusterEvertedInserts.py $insert_size_cutoff > $out_everted_inserts_clustered 2>>$log # For duplications: cluster discordant read pairs into candidate CNVs
echo -e "### clustering of discordant read pairs for duplications ends at" `date +%F'  '%H:%M` "\n" >> $log
}

## 3.2. for deletions
find_cluster_del(){
echo -e "### detection of discordant read pairs for deletions starts at" `date +%F'  '%H:%M` "\n" >> $log
cat $in_sam | python $dir_poolDiffCNV/findDistantInserts.py $insert_size_cutoff > $out_distant_inserts # find discordant read pairs for deletions
echo -e "### detection of discordant read pairs for deletions ends at" `date +%F'  '%H:%M` "\n" >> $log

echo -e "### clustering of discordant read pairs for deletions starts at" `date +%F'  '%H:%M` "\n" >> $log
cat $out_distant_inserts | python $dir_poolDiffCNV/clusterDistantInserts.py $insertSizeDiffCutoff $minNumSupportingInserts $minLenCutoff > $out_distant_inserts_clustered 2>>$log # For deletions: cluster discordant read pairs into candidate CNVs
echo -e "### clustering of discordant read pairs for deletions ends at" `date +%F'  '%H:%M` "\n" >> $log
}

## 3.3. run duplications and deletions parallely at backstage (one of them)
find_cluster_dup & 
find_cluster_del

## 3.4. before moving on, make sure that the background job is completed by cheking the completion string in log file
completion_str="clustering of discordant read pairs for duplications ends"
while ! grep -q "$completion_str" $log ; do sleep 1; done

## 4. Handling output
shopt -s extglob
rm -rf !($out_everted_inserts|$out_distant_inserts|$out_everted_inserts_clustered|$out_distant_inserts_clustered|$log)
shopt -u extglob
