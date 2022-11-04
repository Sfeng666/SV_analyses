# Calculate insert size of each sample library

----
This program calculate insert size (the length of sequenced DNA fragment between adapters) and generate insert size distribition for each sample library.

## Steps
0. Core method: [picard CollectInsertSizeMetrics][1];
1. Prepare input files:
    * deduplicated and realigned bam file of each sample population
2. Calculate insert size distribution of each population sample using [picard CollectInsertSizeMetrics][1]
3. generate a distribution frequency table across all samples as the input to plotting
4. generate customized parameters for CNV analyses from statistics of insert distribution
5. use ridgeplot (ggplot2) to visualize insert size distributions of all samples

## Notes
1. since the input for plotting is a frequency table instead of raw data (could be very large), ggplot would not be able to calculate stats based on it. Therefore I mannually calculated median and 99th percentile of the distribution, and added these labels to each panel.

## Questions
1. there seems a geographical pattern of insert size distribution, but we will not investigate that.

----
Author: [Siyuan Feng](https://scholar.google.com/citations?user=REHFXSsAAAAJ&hl)  
Email: siyuanfeng.bioinfo@gmail.com

----
## References
1. [1]: https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics (picard CollectInsertSizeMetrics)