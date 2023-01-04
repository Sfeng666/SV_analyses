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
3. Obtain allele frequency for reference CNVs in each population (using the Schrider method);
4. Average allele frequency across all populations for each reference CNV and filter CNVs that has average MAF below a MAF threshold (e.g. 0.05);

## Notes
1. Since the clustering steps could take ~4 days to complete for even a single sample library, and it is expected to take no less than 4*29 days to complete for the merged library, we must parallel step 1 by each contig to reduce runtime. Step 2-5 would also benefit from paralleling;
2. For best performance of parallelly running tens of thousands of jobs (29*546), we need to run step 2-5 on CHTC. Paralleling is implemented by the DAGMan work flow of HTCondor.
3. For duplications, it is possible that some population duplicaiton frequency is higher than 1. Possible causes are: 1) random measurement errors that may cause AF slightly greater than 1; 2) the detected duplicaiton is due to sequencing artifacts at least in some samples; 3) duplications are present on some haplotypes in more than two copies. Since case 1 is rare and case 2 is not possible to confirm from our data, it is most reasonable for us to divide each population CNV frequency estimate by the maximum population CNV frequency for any duplication that has at least population CNV frequency > 1. After doing that, the max pop now has duplication freq = 1, and the others are scaled down proportionately. The scaled duplication frequencies might not be true ones, but are proportionate to average copy number in each population which will provide a reasonable startnig point for GEA.

## Questions

----
Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)  
Email: siyuanfeng.bioinfo@gmail.com

----
## References
1. [1]: https://github.com/Sfeng666/poolDiffCNV (poolDiffCNV git repo)