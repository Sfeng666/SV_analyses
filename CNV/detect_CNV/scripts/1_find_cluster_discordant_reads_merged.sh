#! /bin/bash
# script to perform discordant read-pair detection and clustering (step 1) within each sample/across all samples 

# 0. variables & environment
## 0.1. receive arguments from 01qc.sub

contig=$1
conda_pack_squid=$2
insert_size_cutoff=$3
insertSizeDiffCutoff=$4
job=$5
path_in_stage=$6
split=$7
large=$8

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
bam_contig_allsample_stage=$path_in_stage/*_$contig\.bam
bam_contig_allsample=*_$contig\.bam
bed_split=$contig\_$split.bed
bed_split_distant=$contig\_$split\_distant.bed
bed_split_everted=$contig\_$split\_everted.bed

### 0.3.2. output & intermediate files
sam_contig_merged=merged_$contig\_$split.sam
sam_contig_merged_distant=merged_$contig\_$split\_distant.sam
sam_contig_merged_everted=merged_$contig\_$split\_everted.sam
bam_contig_allsample_merged_sorted=merged_sorted_$contig.bam # intermediate output of sorting merged contig bam by read name. Saving this to disk from pipe will save some runtime for repeated sorting for deletions and duplications.
out_everted_inserts=merged_$contig\_$split\_everted_inserts.txt
out_distant_inserts=merged_$contig\_$split\_distant_inserts.txt
out_everted_inserts_clustered=merged_$contig\_$split\_everted_inserts_clustered.txt
out_distant_inserts_clustered=merged_$contig\_$split\_distant_inserts_clustered.txt

### 0.3.3. software parameters
insert_size_cutoff=$insert_size_cutoff
insertSizeDiffCutoff=$insertSizeDiffCutoff
minNumSupportingInserts=2
minLenCutoff=50

### 0.3.3. software paths
dir_poolDiffCNV=biosoft/poolDiffCNV

### 0.3.4. log file
log=log_$job

# 1. transfer contig bam from staging to working directory
cp $bam_contig_allsample_stage ./

# 2. pipe to generate merged split sam file:
echo -e "### preparing split merged sam input starts at" `date +%F'  '%H:%M` "\n" >> $log
if $large
then
    # 2.1. merge contig bams across samples;
    # 2.2. sort by name (required by bedtools pairtobed); 
    samtools merge -u - $bam_contig_allsample \
    | samtools sort -l 0 -n - -o $bam_contig_allsample_merged_sorted

    # 2.3. extract read pairs that interact with the given genomic region (split)
    # 2.4. transfer bam to sam
    ## for distant inserts
    awk 'NR == 1{print}' $bed_split > $bed_split_distant
    start=$(cut -f2 $bed_split_distant)
    end=$(cut -f3 $bed_split_distant)

    bedtools pairtobed -ubam -abam $bam_contig_allsample_merged_sorted -b $bed_split_distant \
    | samtools view - \
    | awk -v start=$start -v end=$end '$1 !~ /@/ && ($4*2 + $9)/2 > start && ($4*2 + $9)/2 < end{print} $1 ~ /@/{print}' - > $sam_contig_merged_distant \
    2>> $log && echo -e "### preparing split merged sam input for del ends at" `date +%F'  '%H:%M` "\n" >> $log &

    ## for everted inserts
    awk 'NR == 2{print}' $bed_split > $bed_split_everted
    start=$(cut -f2 $bed_split_everted)
    end=$(cut -f3 $bed_split_everted)

    bedtools pairtobed -ubam -abam $bam_contig_allsample_merged_sorted -b $bed_split_everted \
    | samtools view - \
    | awk -v start=$start -v end=$end '$1 !~ /@/ && ($4*2 + $9)/2 > start && ($4*2 + $9)/2 < end{print} $1 ~ /@/{print}' - > $sam_contig_merged_everted \
    2>> $log && echo -e "### preparing split merged sam input for del ends at" `date +%F'  '%H:%M` "\n" >> $log
    
    ## wait until the background run is completed to proceed to next step
    completion_str="preparing split merged sam input for del ends"
    while ! grep -q "$completion_str" $log ; do sleep 1; done
else
    samtools merge -u - $bam_contig_allsample \
    | samtools view -o $sam_contig_merged - >> $log 2>&1
fi
echo -e "### preparing split merged sam input ends at" `date +%F'  '%H:%M` "\n" >> $log

# 3. remove unused large files
rm -f $bam_contig_allsample $bam_contig_allsample_merged_sorted *tmp*.bam

# 4. In merged library, find discordant read pairs suggestive of a CNV, and cluster discordant read pairs into candidate CNVs
## using python3 as python2 are no longer built-in at updated CHTC system (centOS8)
## 3.1. for duplications
find_cluster_dup(){
echo -e "### detection of discordant read pairs for duplications starts at" `date +%F'  '%H:%M` "\n" >> $log
if $large
then
    cat $sam_contig_merged_distant | python3 $dir_poolDiffCNV/findEvertedInserts.py > $out_everted_inserts # find discordant read pairs for duplications
else
    cat $sam_contig_merged | python3 $dir_poolDiffCNV/findEvertedInserts.py > $out_everted_inserts # find discordant read pairs for duplications
fi
echo -e "### detection of discordant read pairs for duplications ends at" `date +%F'  '%H:%M` "\n" >> $log

echo -e "### clustering of discordant read pairs for duplications starts at" `date +%F'  '%H:%M` "\n" >> $log
cat $out_everted_inserts | python3 $dir_poolDiffCNV/clusterEvertedInserts.py $insert_size_cutoff > $out_everted_inserts_clustered 2>>$log # For duplications: cluster discordant read pairs into candidate CNVs
echo -e "### clustering of discordant read pairs for duplications ends at" `date +%F'  '%H:%M` "\n" >> $log
}

## 3.2. for deletions
find_cluster_del(){
echo -e "### detection of discordant read pairs for deletions starts at" `date +%F'  '%H:%M` "\n" >> $log
if $large
then
    cat $sam_contig_merged_everted | python3 $dir_poolDiffCNV/findDistantInserts.py $insert_size_cutoff > $out_distant_inserts # find discordant read pairs for deletions
else
    cat $sam_contig_merged | python3 $dir_poolDiffCNV/findDistantInserts.py $insert_size_cutoff > $out_distant_inserts # find discordant read pairs for deletions
fi
echo -e "### detection of discordant read pairs for deletions ends at" `date +%F'  '%H:%M` "\n" >> $log

echo -e "### clustering of discordant read pairs for deletions starts at" `date +%F'  '%H:%M` "\n" >> $log
cat $out_distant_inserts | python3 $dir_poolDiffCNV/clusterDistantInserts.py $insertSizeDiffCutoff $minNumSupportingInserts $minLenCutoff > $out_distant_inserts_clustered 2>>$log # For deletions: cluster discordant read pairs into candidate CNVs
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
