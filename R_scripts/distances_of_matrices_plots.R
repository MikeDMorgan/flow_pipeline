#########################
# plot matrix distances #
#########################
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(WGCNA)
library(flashClust)
library(reshape2)

dist.df <- read.table("/ifs/projects/proj052/pipeline_proj052/clustering.dir/P1-mean.distance",
                      h=T, row.names=1, stringsAsFactors=F)
dist.mat <- as.matrix(dist.df)
colnames(dist.mat) <- gsub(colnames(dist.mat), pattern="-", replacement=".")

hmcol <- rev(colorRampPalette(brewer.pal(9, "BuPu"))(100))
heatmap.2(dist.mat, trace="none", col=hmcol,
          labRow=FALSE, labCol=FALSE, density.info="none")

# PCA across cell subsets
pca.dist <- prcomp(dist.df, center=T)
pcs.dist <- data.frame(pca.dist$x)
parent_cell <- lapply(strsplit(rownames(pcs.dist), split="-", fixed=T),
                      FUN=function(x) {paste(x[1], sep="")})
pcs.dist$parent.cell <- unlist(parent_cell)

pca_p <- ggplot(pcs.dist, aes(x=PC1, y=PC2, colour=parent.cell)) + geom_point() + 
  theme_bw() + labs(x="PC1 (60.0%)", y="PC2 (15.8%)")

# use average-linkage hierarchical clustering
# as.dist and dist create empty distance objects
dist.df$cell <- rownames(dist.df)
dist.melt <- melt(dist.df)
colnames(dist.melt) <- c("cellA", "cellB", "distance")

# using solution from:
# http://stackoverflow.com/questions/11343637/convert-a-dataframe-to-an-object-of-class-dist-without-actually-calculating-di
n <- max(table(dist.melt$cellA)) + 1
res <- lapply(with(dist.melt, split(distance, dist.melt$cellA)), 
              function(x) c(rep(NA, n - length(x)), x))
res <- do.call("rbind", res)
res <- rbind(res, rep(0, n))
dist.res <- dist(t(res))

dist.hc <- flashClust(dist.res, method="average")
dist.cut <- cutreeDynamic(dendro=dist.hc, method="tree",
                          deepSplit=T,
                          pamRespectsDendro=F,
                          minClusterSize=78)
dist.cols <- labels2colors(dist.cut)

plotDendroAndColors(dist.hc, colors=dist.cols, dendroLabels=F, hang=0.03,
                    addGuide=T, guideHang=0.05)
