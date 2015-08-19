#######################################################
# Flow cytometry analysis using Bioconductor packages #
#######################################################

library(flowCore)
library(flowQ)
library(flowWorkspace)
library(flowStats)
library(flowViz)
library(openCyto)
library(flowClust)

source("/ifs/devel/projects/proj052/R_scripts/FlowProcess.R")
# need to set system locale to recognise non-UTF-8 encoded character
Sys.setlocale('LC_ALL', 'C')
path.to.files <- "/ifs/projects/proj052/pipeline_proj052/fcs.dir/ZZFA.dir/" # edit as appropriate
panel <- "P3"
flow.data <- read.flowSet(path = path.to.files, pattern=panel, transformation=FALSE)

# reassign sample names as necessary
twin.split <- strsplit(x=sampleNames(flow.data), split=" ", fixed=T)
twin.ids <- unlist(lapply(twin.split, FUN = function(x) {paste("Twin", x[1], x[2], sep="_")}))
sampleNames(flow.data) <- twin.ids

# add marker IDs to flow.data - no assumption about panel number
flow.data <- set_marker_id(flow.data)

comp <- read.table("/ifs/projects/proj052/pipeline_proj052/comp_matrices.dir/P3-compensation_matrix.txt",
                   h=T, row.names=1)
comp_splt <- strsplit(colnames(comp), split=".", fixed=T)
colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1], x[2], sep="-")}))

fs_comp <- compensate(flow.data, spillover=comp)

# remove any samples with < 10,000 events recorded.
fs_filt <- fs_comp[seq(along=fs_comp)]
extrct <- filter_samples(fs_comp)
if(length(extrct)){fs_filt <- fs_comp[-(extrct),]}

biexpTrans <- flowJoTrans(channelRange=256, maxValue=261589, neg=0, pos=4.4176194777,
                          widthBasis=-100)
tf <- transformList(colnames(fs_filt)[c(4:21)], biexpTrans, transformationId="biexp")
#pars <- colnames(fs_filt)[c(1:21)]
#logicleTf <- estimateLogicle(fs_filt, channels=pars)
#fs_filt <- transformList(fs_filt, tf)
fs_filt <- transform(fs_filt, tf)

# use the gatingSet object instead of workFlow
gs <- GatingSet(fs_filt)
#biexpTrans <- flowJoTrans(channelRange=256, maxValue=261589, neg=0, pos=4.4176194777,
#                          widthBasis=-100)
#tf <- transformList(colnames(fs_filt), biexpTrans)
#recompute(gs)

gtFile <- "/ifs/projects/proj052/R_sessions/P3_cyto.tsv"
gt_tcell <- gatingTemplate(gtFile, autostart=1L)
plot(gt_tcell)
gating(gt_tcell, gs, mc.cores=8, parallel_type="multicore")
#gating(gt_tcell, gs)
rm(fs_comp)
gc()


#wf <- workFlow(fs_filt, name="Twin_P5")
rm(list=c("fs_comp"))
gc()

# need to apply the transformation first
#tf <- transformList(colnames(Data(wf[["base view"]]))[c(1:21)], asinh, transformationId="asinh")
#add(wf, tf)
#fs_trans <- transform(fs_comp, tf)
#gs <- GatingSet(fs_trans)

# boundary filter
#boundFilt <- boundaryFilter(filterId="boundFilt", x=c("SSC.A", "FSC.A"))
#add(wf, boundFilt, parent="asinh")
#gc()

# gate on lymphocytes - preselection of CD3 for T-cells
#lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="V800.A",
#                filterId="TCells", eval=F, scale=1)
#add(wf, lg$n2gate, parent="boundFilt+")

# add a boundary filter on forward and side scatter
#boundFilt <- boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))
#add(gs, boundFilt, parent="root")

# select non-NA markers for further analysis
marker_vec = c(1, 2, 3)
for(i in seq_len(length(pData(parameters(fs_filt[[1]]))[["desc"]]))){
  pname <- pData(parameters(fs_filt[[1]]))[["name"]][i]
  if(pname != "NA"){
    marker_vec = c(marker_vec, i)
  }
}

# use warping algorithm to align fluorescence intensity peaks across measured markers
pars <- colnames(Data(wf[["base view"]]))[marker_vec]
norm <- normalization(normFunction= function(x, parameters, ...) warpSet(x, parameters, ...),
                      parameters=pars, normalizationId = "warping")
add(wf, norm, parent="TCells+")
gc()

# generate and add quadrant gates for CD3+CD4+ and CD3+CD8+ subsets
qgate.cd4 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610.A", "V800.A"), filterId="CD3CD4")
add(wf, qgate.cd4, parent="warping")

qgate.cd8 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610.A", "V585.A"), filterId="CD3CD8")
add(wf, qgate.cd8, parent="warping")

tcell.gate <- quadrantGate(Data(wf[["warping"]]), stains=c("V800.A", "V585.A"), filterId="CD4CD8")
add(wf, tcell.gate, parent="warping")

# add separate gates for memory and naive cells based on CCR7 and CD45RA markers
# generate matrix of limits for 2D rectangle gate
rect.mat <- matrix(c(0, Inf, 0, Inf), ncol=2, dimnames=list(c("min", "max"), c("R710.A", "V655.A")))
mem.gate.cd4 <- rectangleGate(filterId="CD4+CCR7CD45RA", .gate=rect.mat)
mem.gate.cd8 <- rectangleGate(filterId="CD8+CCR7CD45RA", .gate=rect.mat)

add(wf, mem.gate.cd4, parent="CD3+CD4+")
add(wf, mem.gate.cd8, parent="CD3+CD8+")

# function to calculate Fano factor, expression noise, across all individuals and markers in a flowSet
# object
get_fano <- function(flowset) {
  fano.vec <- vector(mode="numeric", length=dim(flowset)[2])
  cols <- colnames(flowset)
  for(i in seq_len(dim(flowset)[2])){
    mean.fs = mean(flowset[,i])
    var.fs = var(flowset[,i])
    fano.fs = var.fs/mean.fs
    fano.vec[i] <- fano.fs
  }
  return(fano.vec)
}

# values into a matrix - test on CD3+CD4+ naive T cell subset
fano.mat <- fsApply(Data(wf[["CD8+CCR7CD45RA-"]]), use.exprs=TRUE, FUN=get_fano)
fano.frame <- data.frame(fano.mat)
colnames(fano.frame) <- parameters(Data(wf[["CD8+CCR7CD45RA-"]])[[1]])$desc
fano.frame$twin.id <- rownames(fano.frame)
View(fano.frame)
gc()

fano.frame$subset <- "ZZFSprt2"
write.table(fano.frame, row.names=T,
            file="/ifs/projects/proj052/flow_processing_tables/ZZFSprt2-panel3-Fano_Tcells.tsv",
            sep="\t")
