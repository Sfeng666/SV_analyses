# Perform Copy Number Variation analyses

----
This program detect copy number variations (CNVs) from whole-genome pool-seq data (pair-end reads) and calculate their allele frequency for further population genetic analyses.

## Steps
### Core method: [poolDiffCNV][1];
### By each sample library:
0. generate customized parameters for CNV analyses (*insertSizeCutoff* and *insertSizeDiffCutoff*) from statistics of insert distribution (done in a previous step in R).

### By each contig:
1. Within & across each population sample:
    * find discordant read pairs that support CNVs (tandem duplications and deletions);
    * cluster supportive read pairs into CNV candidates;
2. Using candidate CNVs identified from the merged pool as a reference, get the global identity of candidate CNVs within each population sample by comparing them to reference CNVs;
3. Obtain allele frequency for reference CNVs in each population;
4. Average allele frequency across all populations for each reference CNV;
5. Determine which CNVs exceed upper and lower frequency threshold (e.g. 0.05–0.95);
### By autosome- and X-linked contigs:
6. Convert population CNV frequency to counts based on X- or Autosome-wide median effective sample size (from SNP analzyed sites)

## Notes
1. Since the clustering steps could take ~4 days to complete for even a single sample library, and it is expected to take no less than 4*29 days to complete for the merged library, we must parallel step 1 by each contig to reduce runtime. Step 2-5 would also benefit from paralleling;
2. For best performance of parallelly running tens of thousands of jobs (29*546), we need to run step 2-5 on CHTC.

## Questions
1. determine the value of parameter *distanceCutoff*


----
Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)  
Email: siyuanfeng.bioinfo@gmail.com

----
## References
1. [1]: https://github.com/Sfeng666/poolDiffCNV (poolDiffCNV git repo)