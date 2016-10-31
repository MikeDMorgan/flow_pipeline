library(ggplot2)
library(reshape2)
library(ICC)
source('/ifs/devel/projects/proj052/flow_pipeline/R_scripts/multiplot.R')

# one file - test
df1 <- read.table("/ifs/projects/proj052/pipeline_proj052/twin_files.dir/CD4_Tcell_fano-Aq.CD45RA.CD3-P2",
                  h=T, stringsAsFactors=F, sep="\t", row.names=1)

df2 <- read.table("/ifs/projects/proj052/pipeline_proj052/twin_files.dir/CD4_Tcell_mean-Aq.CD45RA.CD3-P2",
                  h=T, stringsAsFactors=F, sep="\t", row.names=1)

# remove Aq data
df1 <- subset(df1, select=-which(colnames(df1) %in% c("Aq")))
df2 <- subset(df2, select=-which(colnames(df2) %in% c("Aq")))

melted1 <- melt(df1, id.vars=c("twin.id", "batch", "panel", "gate", "twin_num", "flowjo_id", 
                               "family_id", "zygosity", "age", "replicate", "visit"))

melted2 <- melt(df2, id.vars=c("twin.id", "batch", "panel", "gate", "twin_num", "flowjo_id", 
                               "family_id", "zygosity", "age", "replicate", "visit"))

# standardise as Z-scores for plotting?
z_scores <- list()
for(x in 1:length(unique(melted1$variable))){
  mark <- unique(melted1$variable)[x]
  mean.mark <- mean(melted1$value[melted1$variable == mark])
  sd.mark <- sd(melted1$value[melted1$variable == mark])
  z.mark <- (mean.mark  - melted1$value[melted1$variable == mark])/sd.mark
  z_scores[[mark]] <- z.mark
}
melted1$z.score <- unlist(z_scores)

mz1 <- subset(melted1, subset=melted1$zygosity == "MZ")
dz1 <- subset(melted1, subset=melted1$zygosity == "DZ")

mz1$family_id <- as.factor(mz1$family_id)
mz1$twin <- rep(c("t1", "t2"), dim(mz1)[1]/2)
mz1_split <- data.frame(split(mz1, f=mz1$twin))

dz1$family_id <- as.factor(dz1$family_id)
dz1$twin <- rep(c("t1", "t2"), dim(dz1)[1]/2)
dz1_split <- data.frame(split(dz1, f=dz1$twin))

mz2 <- subset(melted2, subset=melted2$zygosity == "MZ")
dz2 <- subset(melted2, subset=melted2$zygosity == "DZ")

mz2$family_id <- as.factor(mz2$family_id)
mz2$twin <- rep(c("t1", "t2"), dim(mz2)[1]/2)
mz2_split <- data.frame(split(mz2, f=mz2$twin))

dz2$family_id <- as.factor(dz2$family_id)
dz2$twin <- rep(c("t1", "t2"), dim(dz2)[1]/2)
dz2_split <- data.frame(split(dz2, f=dz2$twin))


p_mz <- ggplot(data.frame(mz1_split), aes(x=t1.z.score, y=t2.z.score, colour=t1.variable)) +
  geom_point() + labs(title="Monozygotic twin correlation") + facet_wrap(~t1.variable)

p_dz <- ggplot(data.frame(dz1_split), aes(x=t1.z.score, y=t2.z.score, colour=t2.variable)) + 
  geom_point() + labs(title="Dizygotic twin correlation") + facet_wrap(~t1.variable)

multiplot(p_mz, p_dz)

# iteratively calculate heritabilities
markers <- unique(as.character(melted1$variable))
h2_list1 = list()
for(i in 1:length(markers)){
  m_mz <- mz1[mz1$variable == markers[i],]
  m_dz <- dz1[dz1$variable == markers[i],]
  
  mz_icc <- ICCest(family_id, value, m_mz)
  dz_icc <- ICCest(family_id, value, m_dz)
  
  h2 <- 2 * (mz_icc$ICC - dz_icc$ICC)
  h2_list1[[markers[i]]] <- h2
  h2.df1 <- data.frame(t(data.frame(h2_list1)))
  h2.df1$marker <- rownames(h2.df1)
  colnames(h2.df1) <- c("H2", "marker")
}

h2_list2 = list()
for(i in 1:length(markers)){
  m_mz <- mz2[mz2$variable == markers[i],]
  m_dz <- dz2[dz2$variable == markers[i],]
  
  mz_icc <- ICCest(family_id, value, m_mz)
  dz_icc <- ICCest(family_id, value, m_dz)
  
  h2 <- 2 * (mz_icc$ICC - dz_icc$ICC)
  h2_list2[[markers[i]]] <- h2
  h2.df2 <- data.frame(t(data.frame(h2_list2)))
  h2.df2$marker <- rownames(h2.df2)
  colnames(h2.df2) <- c("H2", "marker")
}

p1_h2 <- ggplot(h2.df1, aes(x=marker, y=H2)) + geom_bar(stat="identity") + theme_bw() + 
  labs(title="H2 of gene expression noise")
p2_h2 <- ggplot(h2.df2, aes(x=marker, y=H2)) + geom_bar(stat="identity") + theme_bw() + 
  labs(title="H2 of mean gene expression")
multiplot(p1_h2, p2_h2)

h2.merge <- merge(h2.df1, h2.df2, "marker")
colnames(h2.merge) <- c("marker", "H2.fano", "H2.mean")
p_h2merge <- ggplot(h2.merge, aes(x=H2.fano, y=H2.mean, colour=marker)) + geom_point(size=4) + 
  theme_bw() + xlim(-1.0, 1.0) + ylim(-1.0, 1.0) + 
  labs(title="Heritability of gene expression noise vs. mean gene expression CD4 Tcell Aq.CD45RA.CD3 ")
print(p_h2merge)
ggsave(file="h2-fano_vs_mean-CD8_Tcell_AqCD45RACD3.png", p_h2merge)