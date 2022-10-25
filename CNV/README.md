# Perform Copy Number Variation analyses

----
This program performs environment association analysis between environmental information and genotypes, in order to detect local adaptation linked to a given environmental variable.

## Steps
0. Core method: BayeScEnv;
1. Prepare input files:
    * Environmental data: environmental information in the form of ‘environmental differentiation’, which is usually computed as a contrast to a reference (usually the average environment, but not obligatory), and as a distance. For each environmental variable, the input is one line that has each column as the distance of each sample (cautious about the order!).
    * Genotype data: 


## Notes
1. Contrary to existing environmental association methods, BayeScEnv does not test for an association between the allelic frequencies and the environment. Rather, the test focuses on a relationship between genetic differentiation (as measured by the local Fst in the model) and environmental differentiation.
2. Spatial resolution of variables is set around 5mins (~100 km).
3. SNPs are analyzed separately by chromosome arm (each autosomal arm + X chromosome arms and assigned X contigs + unmapped auto-linked contigs). Ambiguous contigs (cannot be assigned as X- or autosome-linked) are excluded from this analysis, because the haplotype pool size differs between X and autosomes, which is required to compute effective pool size;
4. SNPs are filtered by the across-samples MAF >= 5% (only include top 2 allele across all samples). Across-samples MAF is computed as: sum(minor allele count of each sample)/sum(minor + major allele count of each sample).
5. the input allele count should be based on effective pool size.
6. BayeScEnv only allow integer, so the effective pool size should be rounded (round half up) before input (because BayeScEnv will automatically round down).

## Questions
1. How to narrow down SNPs for analysis:
    * filter SNPs by AF or allele count - average AF across all samples (e.g. 5%)
    * filter by global Fst at each SNP (local adaptation sites tend to have higher Fst)
        * advantages: run faster; avoid multiple-test issues
2. Which environment variables to include:
    * Ratio of built area to vegetation (trees + crops)
    * Ratio of crop to trees
    * Altitude
    * Wind Speed (annual mean)
    * Mean temperature of the coldest quarter
    * Mean temperature of the warmest quarter (correlated with wind)
    * mean Diurnal Range
    * annual precipitation
    * precipitation seasonality
----
Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)  
Email: siyuanfeng.bioinfo@gmail.com

----
## References
Basic concepts and ideas about using NGS to detect