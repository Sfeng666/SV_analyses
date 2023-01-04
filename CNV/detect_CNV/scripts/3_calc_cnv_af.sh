#! /bin/bash
# script to ccalculate allele frequency of CNVs from each contig of each samoke

# 0. variables & environment
## 0.1. receive arguments from 01qc.sub

contig=$1
sample=$2
conda_pack_squid=$3
path_in_stage=$4
job=$5

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
bed_del_sp_contig=del_$sample\_$contig\.bed
bed_dup_sp_contig=dup_$sample\_$contig\.bed
bam_sample_contig=$path_in_stage/$sample\_$contig\.bam

### 0.3.2. output & intermediate files
del_breakpoint_rp_sorted=del_breakpoint_rp_sorted.bam # uncompressed bam files containing read pairs across breakpoints of deletions
dup_breakpoint_rp_sorted=dup_breakpoint_rp_sorted.bam # uncompressed bam files containing read pairs across breakpoints of duplications
out_del_AF=del_af_$sample\_$contig\.txt # ouput table containing deletion position and allele frequency
out_dup_AF=dup_af_$sample\_$contig\.txt # ouput table containing duplication position and allele frequency

### 0.3.3. log file
log=log_$job

### 0.3.4. transfer input from staging
cp $bam_sample_contig ./

# 1. calculate allele frequency of CNVs for each population
## 1.1. define a function of AF calculation for deletion
calc_af_del(){
echo -e "### calculation of allele frequency for deletions starts at" `date +%F'  '%H:%M` "\n" >> $log
### 1.1.1. generate a trimmed bam file that contain all read pairs across breakpoints (may also contain some unwanted read pairs, e.g. |..<----...|----> )
samtools sort -l 0 -n $(basename $bam_sample_contig) | # sort by name (required by bedtools pairtobed); 
bedtools pairtobed -ubam -abam stdin -b $bed_del_sp_contig -f 1 -type xor > $del_breakpoint_rp_sorted # extract read pairs that have one end within CNV and another end flanking
### 1.1.2. calculate allele frequency for each CNV by further extracting qualified read pairs at breakpoints
while read -r line; do
    current_coord=current_coord_del.bed
    echo "$line" > $current_coord
    ct=$(bedtools pairtobed -bedpe -abam $del_breakpoint_rp_sorted -b $current_coord -f 1 -type xor | # count read pairs that have one end within and another end flanking a CNV
    awk '$9 == "+" && $10 == "-"{print}' | # only include non-everted reads
    wc -l)
    awk -v ct=$ct 'ct != 0{printf "%s\t%s\t%s\t%d\t%d\t%.6f\n", $1, $2, $3, $4, ct, $4/($4 + (ct/2))} ct == 0{printf "%s\t%s\t%s\t%d\t%d\t1\n", $1, $2, $3, $4, ct}' $current_coord >> $out_del_AF # formula for AF calculation: p = Np/(Np + (q/2)), where Np is the number of del-suggestive read pairs, Np+q is the number of breakpoints read pairs
done < $bed_del_sp_contig
echo -e "### calculation of allele frequency for deletions ends at" `date +%F'  '%H:%M` "\n" >> $log
}

## 1.2. define a function of AF calculation for duplications
calc_af_dup(){
echo -e "### calculation of allele frequency for duplications starts at" `date +%F'  '%H:%M` "\n" >> $log
### 1.2.1. generate a trimmed bam file that contain all read pairs across breakpoints (may also contain some unwanted read pairs, e.g. |..<----...|----> )
samtools sort -l 0 -n $(basename $bam_sample_contig) | # sort by name (required by bedtools pairtobed); 
bedtools pairtobed -ubam -abam stdin -b $bed_dup_sp_contig -f 1 -type xor > $dup_breakpoint_rp_sorted # extract read pairs that have one end within CNV and another end flanking
### 1.2.2. calculate allele frequency for each CNV by further extracting qualified read pairs at breakpoints
while read -r line; do
    current_coord=current_coord_dup.bed
    echo "$line" > $current_coord
    ct=$(bedtools pairtobed -bedpe -abam $dup_breakpoint_rp_sorted -b $current_coord -f 1 -type xor | # count read pairs that have one end within and another end flanking a CNV
    awk '$9 == "+" && $10 == "-"{print}' | # only include non-everted reads
    wc -l)
    awk -v ct=$ct 'ct != 0{printf "%s\t%s\t%s\t%d\t%d\t%.6f\n", $1, $2, $3, $4, ct, $4/(ct/2)} ct == 0{printf "%s\t%s\t%s\t%d\t%d\t0\n", $1, $2, $3, $4, ct}' $current_coord >> $out_dup_AF # formula for AF calculation: p = Np/(Np+q/2), where Np is the number of dup-suggestive read pairs, Np+q is the number of breakpoints read pairs
done < $bed_dup_sp_contig
echo -e "### calculation of allele frequency for duplications ends at" `date +%F'  '%H:%M` "\n" >> $log
}

## 1.3. run duplications and deletions parallely at backstage (one of them)
calc_af_del 2>> $log & 
calc_af_dup 2>> $log

## 1.4. before moving on, make sure that the background job is completed by cheking the completion string in log file
completion_str="calculation of allele frequency for deletions ends"
while ! grep -q "$completion_str" $log && ! grep -q fatal $log ; do sleep 1; done

## 2. Handling output
shopt -s extglob
rm -rf !($out_del_AF|$out_dup_AF|$log)
shopt -u extglob
