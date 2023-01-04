#! /bin/bash
# script to compare population CNVs to reference CNVs for each population

# 0. variables & environment
## 0.1. receive arguments from 01qc.sub

contig=$1
sample=$2
insertSizeDiffCutoff=$3
job=$4

## 0.2. set local paths
### 0.2.1. input files
merged_everted_inserts_clustered=merged_$contig\_everted_inserts_clustered.txt
merged_distant_inserts_clustered=merged_$contig\_distant_inserts_clustered.txt
sample_everted_inserts_clustered=$sample\_$contig\_everted_inserts_clustered.txt
sample_distant_inserts_clustered=$sample\_$contig\_distant_inserts_clustered.txt

### 0.2.2. output & intermediate files
compared_everted_clusters=$contig\_everted_clusters_$sample\.txt
compared_distant_clusters=$contig\_distant_clusters_$sample\.txt

### 0.2.3. software parameters
insertSizeDiffCutoff=$insertSizeDiffCutoff
normconst1=1
normconst2=1

### 0.2.4. software paths
dir_poolDiffCNV=biosoft/poolDiffCNV

### 0.2.5. log file
log=log_$job

# 1. compare population CNVs to reference CNVs
## 1.1. for duplications
## using python3 as python2 are no longer built-in at updated CHTC system (centOS8)
compare_cluster_dup(){
echo -e "### comparison of clusters for duplications starts at" `date +%F'  '%H:%M` "\n" >> $log
python3 $dir_poolDiffCNV/combineEvertedClustersAcrossPools_ref.py $merged_everted_inserts_clustered $sample_everted_inserts_clustered $normconst1 $normconst2 $insertSizeDiffCutoff > $compared_everted_clusters 2>>$log
echo -e "### comparison of clusters for duplications ends at" `date +%F'  '%H:%M` "\n" >> $log
}

## 1.2. for deletions
compare_cluster_del(){
echo -e "### comparison of clusters for deletions starts at" `date +%F'  '%H:%M` "\n" >> $log
python3 $dir_poolDiffCNV/combineDistantClustersAcrossPools_ref.py $merged_distant_inserts_clustered $sample_distant_inserts_clustered $normconst1 $normconst2 $insertSizeDiffCutoff > $compared_distant_clusters 2>>$log
echo -e "### comparison of clusters for deletions ends at" `date +%F'  '%H:%M` "\n" >> $log
}

## 1.3. run duplications and deletions parallely at backstage (one of them)
compare_cluster_dup & 
compare_cluster_del

## 1.4. before moving on, make sure that the background job is completed by cheking the completion string in log file
completion_str="comparison of clusters for duplications ends"
while ! grep -q "$completion_str" $log ; do sleep 1; done

## 2. Handling output
shopt -s extglob
rm -rf !($compared_everted_clusters|$compared_distant_clusters|$log)
shopt -u extglob