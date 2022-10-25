# Perform Structural Variation analyses

----
This program detect CNV and inversions in natural populations of *Drosophila suzukii* from pool-seq WGS data.

## Steps

## Notes

## Questions

----
Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)  
Email: siyuanfeng.bioinfo@gmail.com

----
## References
### Basic concepts and ideas about using NGS to detect SV
1. [Clarify Paired-end read confusion - library, fragment or insert size?](https://thegenomefactory.blogspot.com/2013/08/paired-end-read-confusion-library.html?lr=1)
2. [paired-end read construction during libary preparation - Illimina](https://support.illumina.com/bulletins/2020/12/how-short-inserts-affect-sequencing-performance.html)
3. [definitions for structural variants as a reference for GATK applications](https://gatk.broadinstitute.org/hc/en-us/articles/9022476791323-Structural-Variants)
### Reference to Picard CollectInsertSizeMetrics
1. [Mannual of Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics)
2. [definitions for metrics of Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics)
3. [extended reading - different read orientations: Forward-Reverse (FR), Reverse-Forward (RF), Tandem](https://gatk.broadinstitute.org/hc/en-us/articles/360035531792-Paired-end-or-mate-pair)
4. [extended reading - difference between pair-end and mate-pair sequencing](https://www.ecseq.com/support/ngs/what-is-mate-pair-sequencing-useful-for)