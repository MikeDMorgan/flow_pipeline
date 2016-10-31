##################
# Twins plotting #
##################
library(reshape2)
library(ggplot2)
# load example data sets, generate plots - use CD4+ T cells/naive T cells
fano.df <- read.table("/ifs/projects/proj052/pipeline_proj052/twin_files.dir/CD8_Tmem_fano-CD4.CCR10.CD3-P3",
                 h=T, sep="\t", row.names=1)

mean.df <- read.table("/ifs/projects/proj052/pipeline_proj052/twin_files.dir/CD8_Tmem_mean-CD4.CCR10.CD3-P3",
                      h=T, sep="\t", row.names=1)

mean.melt <- melt(mean.df, id.vars=(c("twin.id", "batch", "panel",
                                      "gate", "twin_num", "flowjo_id",
                                      "family_id", "zygosity",
                                      "age", "replicate", 
                                      "visit")))

fano.melt <- melt(fano.df, id.vars=(c("twin.id", "batch", "panel",
                                      "gate", "twin_num", "flowjo_id",
                                      "family_id", "zygosity",
                                      "age", "replicate", 
                                      "visit")))

merge.melted <- merge(fano.melt, mean.melt, by.x=c("twin.id", "batch", "panel",
                                                   "gate", "twin_num", "flowjo_id",
                                                   "family_id", "zygosity",
                                                   "age", "replicate", 
                                                   "visit", "variable"),
                      by.y=c("twin.id", "batch", "panel",
                             "gate", "twin_num", "flowjo_id",
                             "family_id", "zygosity",
                             "age", "replicate", 
                             "visit", "variable"))

colnames(merge.melted) <- gsub(gsub(gsub(colnames(merge.melted),
                                    pattern="value.x", replacement="fano"),
                               pattern="value.y", replacement="mean"),
                               pattern="variable", replacement="marker")

m_f_p <- ggplot(merge.melted, aes(x=fano, y=mean, colour=marker)) + geom_point() + 
  facet_wrap(~marker) + scale_x_log10() + scale_y_log10() + theme_bw() + 
  labs(x="Expression noise, log10", y="Mean Expression, log10", 
       title="CD8 T memory CCR10+")

fano_hist <- ggplot(merge.melted, aes(fano, fill=marker)) + theme_bw() + 
  geom_histogram() + facet_wrap(~marker, scale="free_x") + 
  labs(x="Unadjusted gene expression noise", y="Counts") +
  theme(axis.title=element_text(size=14, colour="black"))

png("Histogram-fano-CD8_Tmem_CD4.CCR10.CD3", height=720, width=720)
print(fano_hist)
dev.off()

mean_hist <- ggplot(merge.melted, aes(mean, fill=marker)) + theme_bw() + 
  geom_histogram() + facet_wrap(~marker, scale="free_x") + 
  labs(x="Unadjusted mean gene expression", y="Counts") +
  theme(axis.title=element_text(size=14, colour="black"))

png("Histogram-mean-CD8_Tmem_CD4.CCR10.CD3", height=720, width=720)
print(mean_hist)
dev.off()

# calculate correlations
markers <- as.character(unique(merge.melted$marker))
cor_list <- list()
for(x in 1:length(markers)){
  dataset <- merge.melted[merge.melted$marker == markers[x],]
  mcor <- cor(dataset$fano, dataset$mean)
  cor_list[[markers[x]]] <- mcor
}

# age vs fano and mean
f_age_p <- ggplot(merge.melted, aes(x=age, y=fano, colour=marker)) + geom_point() + 
  facet_wrap(~marker) + scale_y_log10() + theme_bw() + 
  labs(x="Age (years)", y="Expression noise, log10", title="CD8 T memory CCR10+")


m_age_p <- ggplot(merge.melted, aes(x=age, y=mean, colour=marker)) + geom_point() + 
  facet_wrap(~marker) + scale_y_log10() + theme_bw() + 
  labs(x="Age (years)", y="Mean Expression, log10", title="CD8 T memory CCR10+")



fano.h2 <- read.table("/ifs/projects/proj052/pipeline_proj052/heritability.dir/CD8_Tmem_fano-CD4.CCR10.CD3-P3.heritability",
                      h=T)
fano.h2$marker <- rownames(fano.h2)

mean.h2 <- read.table("/ifs/projects/proj052/pipeline_proj052/heritability.dir/CD8_Tmem_mean-CD4.CCR10.CD3-P3.heritability",
                      h=T)
mean.h2$marker <- rownames(mean.h2)

h2.merge <- merge(mean.h2, fano.h2, "marker")
colnames(h2.merge) <- c("marker", "h2_mean", "h2_fano")
h2.merge$marker <- gsub(h2.merge$marker, pattern=":", replacement="")

h2.melt <- data.frame(melt(h2.merge, "marker"))
colnames(h2.melt) <- c("marker", "stat", "h2")

h2_bar <- ggplot(h2.melt, aes(y=h2, x=marker)) + geom_bar(stat="identity" , fill="#473C8B") + 
  facet_wrap(~stat) + theme_bw() + theme(axis.text.x=element_text(angle=45, vjust=0.6,
                                                                  size=14, colour="black")) + 
  theme(axis.text.y=element_text(size=14, colour="black")) + 
  theme(axis.title=element_text(size=14, colour="black")) +
  labs(x="Marker", y=expression(paste("Broad-sense heritability ", H^2, sep="")))

png("H2-barchart-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(h2_bar)
dev.off()

# plot noise vs mean h2
h2vsh2 <- ggplot(h2.merge, aes(x=h2_fano, y=h2_mean, colour=marker)) + geom_point(size=4) + 
  theme_bw() + theme(axis.text=element_text(size=14, colour="black")) + 
  labs(x=expression(paste(H^2, " of expression noise", sep="")),
       y=expression(paste(H^2, " of mean expression", sep=""))) + 
  theme(axis.title=element_text(size=14, colour="black")) +
  coord_fixed(ratio=1) + xlim(c(min(h2.merge$h2_mean, h2.merge$h2_fano),
                                max(h2.merge$h2_mean, h2.merge$h2_fano))) + 
  ylim(c(min(h2.merge$h2_mean, h2.merge$h2_fano),
         max(h2.merge$h2_mean, h2.merge$h2_fano)))

png("h2v_vs_h2-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(h2vsh2)
dev.off()

# mean-noise dependence based on zygosity
m_f.byZyg <- ggplot(merge.melted[merge.melted$zygosity != "Rep",],
                    aes(x=fano, y=mean, colour=marker)) + 
  geom_point() + theme_bw() + facet_wrap(~zygosity) +
  scale_x_log10() + scale_y_log10() + labs(x="Expression noise, log10",
                                           y="Mean expression, log10")

png("CD8_Tmem_fano-CD4.CCR10.CD3-mean_vs_fano-by_zygosity.png", height=720, width=720)
print(m_f.byZyg)
dev.off()

# show values after adjustment for age and batch effect
mz_mean_adj <- read.table("/ifs/projects/proj052/pipeline_proj052/residuals.dir/MZ-CD8_Tmem_mean-CD4.CCR10.CD3-P3.adjusted",
                          h=T, sep="\t")
dz_mean_adj <- read.table("/ifs/projects/proj052/pipeline_proj052/residuals.dir/DZ-CD8_Tmem_mean-CD4.CCR10.CD3-P3.adjusted",
                          h=T, sep="\t")
mean_adj <- rbind.data.frame(mz_mean_adj, dz_mean_adj)

mz_fano_adj <- read.table("/ifs/projects/proj052/pipeline_proj052/residuals.dir/MZ-CD8_Tmem_fano-CD4.CCR10.CD3-P3.adjusted",
                          h=T, sep="\t")
dz_fano_adj <- read.table("/ifs/projects/proj052/pipeline_proj052/residuals.dir/DZ-CD8_Tmem_fano-CD4.CCR10.CD3-P3.adjusted",
                          h=T, sep="\t")
fano_adj <- rbind.data.frame(mz_fano_adj, dz_fano_adj)

merge_adj <- merge(fano_adj, mean_adj, 
                   c("indx", "twin.id_x", "twin_num_x", 
                     "family_id", "flowjo_id_x", "age_x", 
                     "gate", "zygosity", "replicate_x",
                     "visit_x", "batch", "panel_x", "marker",
                     "twin.id_y", "twin_num_y", "flowjo_id_y",
                     "age_y", "replicate_y", "visit_y", "panel_y"))

colnames(merge_adj) <- gsub(gsub(colnames(merge_adj), pattern=".x",
                                 replacement=".fano"), pattern=".y",
                            replacement=".mean")

colnames(merge_adj) <- gsub(gsub(colnames(merge_adj),
                                 pattern="twin.id.fano",
                                 replacement="twin.id.twin1"), 
                            pattern="twin.id.mean",
                            replacement="twin.id.twin2")

colnames(merge_adj) <- gsub(colnames(merge_adj),
                            pattern=".meangosi.mean", 
                            replacement="zygosity")

# plot twin concordances using adjusted values
twin_mean_conc <- ggplot(merge_adj, aes(x=twin1_res.mean,
                                        y=twin2_res.mean,
                                        colour=marker)) + geom_point() + 
  theme_bw() + facet_wrap(~marker, scale="free") + 
  stat_smooth(method="lm") + labs(x="Twin1 adjusted mean expression", 
                                  y="Twin2 adjusted mean expression")

png("Twin_vs_twin-mean-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(twin_mean_conc)
dev.off()

twin_fano_conc <- ggplot(merge_adj, aes(x=twin1_res.fano,
                                        y=twin2_res.fano,
                                        colour=marker)) + geom_point() + 
  theme_bw() + facet_wrap(~marker, scale="free") + 
  stat_smooth(method="lm") + labs(x="Twin1 adjusted expression noise", 
                                  y="Twin2 adjusted expression noise")

png("Twin_vs_twin-fano-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(twin_fano_conc)
dev.off()


# plot adjusted vs. unadjusted values
# put together twin1 and twin2 first
twin1.df <- data.frame(merge_adj$twin1_res.mean, merge_adj$twin1.mean, 
                       merge_adj$twin1_res.fano, merge_adj$twin1.fano,
                       merge_adj$marker, merge_adj$zygosity)
colnames(twin1.df) <- c("adj.mean", "raw.mean", "adj.fano", "raw.fano",
                        "marker", "zygosity")

twin2.df <- data.frame(merge_adj$twin2_res.mean, merge_adj$twin2.mean, 
                       merge_adj$twin2_res.fano, merge_adj$twin2.fano,
                       merge_adj$marker, merge_adj$zygosity)
colnames(twin2.df) <- c("adj.mean", "raw.mean", "adj.fano", "raw.fano",
                        "marker", "zygosity")

twin.df <- rbind.data.frame(twin1.df, twin2.df)

mean_adj_raw <- ggplot(twin.df, aes(x=adj.mean,
                                    y=raw.mean,
                                    colour=marker)) + geom_point() + 
  theme_bw() + facet_wrap(~marker, scale="free") + 
  stat_smooth(method="lm") + labs(x="Adjusted mean expression", 
                                  y="Unadjusted mean expression")

png("Adj_vs_unadj-mean-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(mean_adj_raw)
dev.off()


fano_adj_raw <- ggplot(twin.df, aes(x=adj.fano,
                                    y=raw.fano,
                                    colour=marker)) + geom_point() + 
  theme_bw() + facet_wrap(~marker, scale="free") + 
  stat_smooth(method="lm") + labs(x="Adjusted expression noise", 
                                  y="Unadjusted expression noise")

png("Adj_vs_unadj-fano-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(fano_adj_raw)
dev.off()

adj_fano_hist <- ggplot(twin.df, aes(adj.fano, fill=marker)) + theme_bw() + 
  geom_histogram() + facet_wrap(~marker, scale="free_x") + 
  labs(x="Adjusted gene expression noise", y="Counts") +
  theme(axis.title=element_text(size=14, colour="black"))

png("Histogram-adjusted_fano-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(adj_fano_hist)
dev.off()

adj_mean_hist <- ggplot(twin.df, aes(adj.mean, fill=marker)) + theme_bw() + 
  geom_histogram() + facet_wrap(~marker, scale="free_x") + 
  labs(x="Adjusted meand gene expression", y="Counts") +
  theme(axis.title=element_text(size=14, colour="black"))

png("Histogram-adjusted_mean-CD8_Tmem_CD4.CCR10.CD3.png", height=720, width=720)
print(adj_mean_hist)
dev.off()

