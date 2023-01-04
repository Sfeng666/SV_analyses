library(gplots)

# receive cmd arguments
Args <- commandArgs(T)
wd <- Args[1]
dd <- Args[2]
para <- Args[3]

# # # parameters for test only
# wd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/Fst/plot_efs"
# dd <- "/Users/siyuansmac/bioinfo/project/suzukii_WGS/genetic_diff/Fst/ungapped_window_efs"
# para <- "winheter_100_mincov_12_mincount_1"
# chr <- "auto"

# def a function to plot heatmap for Fst
heatmap_Fst <- function(wd, dd, para, chr) {
  
## read the table including Fst
setwd(wd)
df <- read.table(file=paste(dd, "/", para, "/", "Fst_", chr, ".txt", sep = ""), 
                 header=T, row.names = 1, sep="\t", check.names=FALSE)
df_color <- data.frame(sample=colnames(df), 
                       range=factor(c("FR-Run", "European", "European", "Chinese",
                                      "Chinese", "European", "European", "European",
                                      "European", "Japanese", "Japanese", "Chinese", 
                                      "American", "American", "Chinese", "American", 
                                      "American", "American", "US-Haw", "American",
                                      "BR-Pal", "European", "American", "American",
                                      "Japanese", "European", "Japanese", "European", "European"),
                                    levels = c("US-Haw", "American", "BR-Pal", "FR-Run", "European", "Chinese", "Japanese")))
df_color <- df_color[order(df_color$sample),]
df_color <- df_color[order(df_color$range),]
df_color$range <- as.character(df_color$range)
mt <- as.matrix(df[df_color$sample, df_color$sample])

# update 10.09.2021: separate island populations from the mainland population ranges, and assign separate color to them
vc_color <- c(rgb(232, 125,	114, maxColorValue = 255), rgb(83, 182,76, maxColorValue = 255), rgb(109, 157, 248, maxColorValue = 255), rgb(109 - 35, 157 + 35, 248, maxColorValue = 255),
              rgb(232, 125,	114 + 70, maxColorValue = 255), rgb(232 -50, 125,	114 + 20, maxColorValue = 255), rgb(83 + 70, 182,76, maxColorValue = 255))
names(vc_color) <- c("American", "European", "Chinese", "Japanese", "US-Haw", "BR-Pal","FR-Run")
vc_spcolor <- vc_color[df_color$range]
df_color$color <- vc_spcolor

## plot the heatmap without the dendrogram
pdf(file=paste(wd, '/', "Fst_", chr, "_", para, ".pdf", sep = ""))
par(ps = 12, cex = 1, cex.main = 1)
heatmap.2(mt, margins=c(8, 8), 
          ColSideColors = df_color$color,
          col=colorpanel(100,low="#edf8b1",high="#2c7fb8"), 
          dendrogram="none", Rowv=F, Colv=F, adjCol=c(1, 0.5), 
          key=TRUE, key.xlab = 'Pairwise Fst',density.info="none", trace="none")
legend(-0.06, 0.82,legend=names(vc_color),title="Range", cex = 1,
       fill=vc_color, bty="n", xpd = T)
dev.off()

## plot the heatmap and dendrogram based on hierarchical clustering of distances as Fst
pdf(file=paste(wd, '/', "Fst_", chr, "_", para, "_hclust_byFst.pdf", sep = ""))
par(ps = 12, cex = 1, cex.main = 1)
heatmap.2(mt, margins=c(8, 8), 
          ColSideColors = df_color$color,
          col=colorpanel(100,low="#edf8b1",high="#2c7fb8"), dendrogram="column",
          Colv = as.dendrogram(hclust(as.dist(mt))), Rowv = as.dendrogram(hclust(as.dist(mt))), revC = T, adjCol=c(1, 0.5),
          key=TRUE, key.xlab = 'Pairwise Fst',density.info="none", trace="none")
legend(-0.06, 0.82,legend=names(vc_color),title="Range", cex = 1,
       fill=vc_color, bty="n", xpd = T)
dev.off()

## plot the heatmap and dendrogram based on hierarchical clustering of distances calculated from Fst
pdf(file=paste(wd, '/', "Fst_", chr, "_", para, "_hclust.pdf", sep = ""))
par(ps = 12, cex = 1, cex.main = 1)
heatmap.2(mt, margins=c(8, 8), 
          ColSideColors = df_color$color, 
          col=colorpanel(100,low="#edf8b1",high="#2c7fb8"), dendrogram="column", 
          Colv = as.dendrogram(hclust(dist(mt))), Rowv = as.dendrogram(hclust(dist(mt))), revC = T, adjCol=c(1, 0.5),
          key=TRUE, key.xlab = 'Pairwise Fst',density.info="none", trace="none")
legend(-0.06, 0.82,legend=names(vc_color),title="Range", cex = 1,
       fill=vc_color, bty="n", xpd = T)
dev.off()
}

# plot
heatmap_Fst(wd, dd, para, "auto")
heatmap_Fst(wd, dd, para, "X")