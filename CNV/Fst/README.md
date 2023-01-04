# Estimate genome-wide Reynolds Fst for CNVs

----
This program estimates genome-wide Fst (Reynolds, weighted average across CNVs) as a quality check of estimation of CNV allele frequency.

## Steps
Within autosomes and X-chromosome:      
1. From results of previous analyses, retrieve effective sample size of chromosomes for each population as the median effective sample size of analyzed SNPs;
2. Estimate pairwise Reynolds Fst for each CNV;
3. Calculate 'genome-wide' (i.e. autosome-wide/X-chromosome-wide) pairwise Fst as the weighter average of step 2.

## Notes

## Questions
----
Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)  
Email: siyuanfeng.bioinfo@gmail.com

----
## References