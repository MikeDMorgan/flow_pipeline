# correlation of matrices
# uses canonical correlations, cancor in R
library(gplots)
library(RColorBrewer)
library(ggplot2)

df1 <- read.table("/ifs/projects/proj052/pipeline_proj052/twin_files.dir/CD4_Tcell_mean-CD127.CCR7.CD244-P1",
                  h=T, sep="\t", row.names=1)

df2 <- read.table("/ifs/projects/proj052/pipeline_proj052/twin_files.dir/CD8_Tcell_mean-AQ.CD3.CD25-P2",
                  h=T, sep="\t", row.names=1)
mat1 <- data.matrix(df1[,1:13])
mat2 <- data.matrix(df2[,1:13])

mat_cor <- cor(mat1, mat2)
mat_cancor <- cancor(mat1, mat2)$cor
hmcol <- colorRampPalette(brewer.pal(9, "PuOr"))(100)

heatmap.2(mat_cor, col=hmcol, trace="none")

summary(mat_cancor)

# calculate distances between matrices using the Frobenius difference
file_list = list.files("/ifs/projects/proj052/pipeline_proj052/twin_files.dir", 
                       full.names=T)
cd4.files <- file_list[grepl(x=file_list, pattern="CD4")]
p1.files <- cd4.files[grepl(cd4.files, pattern="P3")]
my.files <- p1.files[!grepl(p1.files, pattern="(log|tsv)")]

matrix_list = lapply(my.files, FUN=function(x) {read.table(x, sep="\t", h=T)})

for(j in 1:length(matrix_list)){
  rownames(matrix_list[[j]]) <- matrix_list[[j]]$twin.id
}

mats_dist <- lapply(matrix_list, FUN=function(x) {data.matrix(x[,1:13], 
                                                              rownames.force=T)})

# get all possible combinations
all_pairs <- permutations(length(mats_dist), 2, set=T, repeats.allowed=T)

matrixDistance <- function(matrixList, matrixCombinations){
  distMatrix <- matrix(, ncol=length(matrixList),
                         nrow=length(matrixList))
  for(x in 1:length(matrixCombinations)){
    # dimensions of each matrix need to match
    # to be conformable
    gc()
    x1 <- matrixCombinations[x,][1]
    x2 <- matrixCombinations[x,][2]
    mat1 <- matrixList[[x1]]
    mat2 <- matrixList[[x2]]
    names1 <- rownames(mat1)
    names2 <- rownames(mat2)
    setnames1 <- names1[names1 %in% names2]
    setnames2 <- names2[names2 %in% names1]
    intersect.names <- intersect(setnames1, setnames2)

    matDist <- sqrt(sum((mat1[intersect.names,] - mat2[intersect.names,]) ** 2))
    distMatrix[x1, x2] <- matDist
  }
  distMatrix
}