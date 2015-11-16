#######################################################
# Flow cytometry analysis using Bioconductor packages #
#######################################################

library(flowCore)
library(flowWorkspace)
library(flowStats)
library(flowViz)
library(openCyto)
library(flowClust)

source("/ifs/devel/projects/proj052/flow_pipeline/R_scripts/FlowProcess.R")
# need to set system locale to recognise non-UTF-8 encoded character
Sys.setlocale('LC_ALL', 'C')
path.to.files <- "/ifs/projects/proj052/pipeline_proj052/fcs.dir/ZZFL.dir/" # edit as appropriate
panel <- "P1"
flow.data <- read.flowSet(path = path.to.files, pattern=panel, transformation=FALSE)

# reassign sample names as necessary
twin.split <- strsplit(x=sampleNames(flow.data), split=" ", fixed=T)
twin.ids <- unlist(lapply(twin.split, FUN = function(x) {paste("Twin", x[1], x[2], sep="_")}))
sampleNames(flow.data) <- twin.ids

# add marker IDs to flow.data - no assumption about panel number
flow.data <- set_marker_id(flow.data)
flow.data <- fix_scatter_name(flow.data)

flow.data <- flow.data[1:6]
comp <- read.table("/ifs/projects/proj052/pipeline_proj052/comp_matrices.dir/P1-compensation_matrix.txt",
                   h=T, row.names=1)
comp_splt <- strsplit(colnames(comp), split=".", fixed=T)
colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1], x[2], sep="-")}))

fs_comp <- compensate(flow.data, spillover=comp)

# remove any samples with < 10,000 events recorded.
fs_filt <- fs_comp[seq(along=fs_comp)]
extrct <- filter_samples(fs_comp)
if(length(extrct)){fs_filt <- fs_comp[-(extrct),]}

# apply biexponential transformation to forward and side scatter data, and logicle or
# asinh to flourescence parameters
# biexpTrans <- flowJoTrans(channelRange=256, maxValue=261589, neg=0, pos=4.4176194777,#
#                           widthBasis=-100)
# 
# tf <- transformList(colnames(fs_filt)[c(1:3)], biexpTrans, transformationId="biexp")
# tf <- transformList(colnames(fs_filt)[c(1:21)], asinh, transformationId="asinh")
# 
# pars <- colnames(fs_filt)[c(4:21)]
# logicleTf <- estimateLogicle(fs_filt[[1]], channels=pars)
# fs_filt <- transform(fs_filt, logicleTf)
# fs_filt <- transform(fs_filt, tf)
# 
# # use the gatingSet object to pull out lymphocytes from a gating hierarchy
# lymph_file <- "/ifs/projects/proj052/R_sessions/just_lymphs.txt"
# lymphs <- gatingTemplate(lymph_file, autostart=1L)
# gs_lymphs <- GatingSet(fs_filt)
# gc()
# #gating(lymphs, gs_lymphs, mc.cores=2, parallel_type="multicore")
# gating(lymphs, gs_lymphs)
# proper_lymphs <- getData(gs_lymphs, "/lymphs")
#        
# wf <- workFlow(proper_lymphs, name="workflow")
# pars <- colnames(Data(wf[["base view"]]))[c(4:21)]
# norm <- normalization(normFunction= function(x, parameters, ...) warpSet(x, parameters, ...),
#                       parameters=pars, normalizationId = "warping")
# add(wf, norm, parent="base view")
# fs_norm <- Data(wf[["warping"]])
# 
# gtFile <- "/ifs/projects/proj052/R_sessions/P3_cyto.tsv"
# gt_tcell <- gatingTemplate(gtFile, autostart=1L)
# plot(gt_tcell)
# gs <- GatingSet(fs_norm)
# rm(wf, gs_lymphs, fs_filt, flow.data, proper_lymphs, fs_comp, logicleTf, tf)
# gc()
# #gating(gt_tcell, gs, mc.cores=8, parallel_type="multicore")
# gating(gt_tcell, gs)
# gc()
# 
# #wf <- workFlow(fs_filt, name="Twin_P5")
# rm(list=c("fs_comp"))
# gc()
# 
# 
# 
#biexpTrans <- flowJoTrans(channelRange=256, maxValue=261589, neg=0, pos=4.4176194777,
#                          widthBasis=-100)
#tf_s <- transformList(c("FSC-A", "FSC-H", "SSC-A"), biexpTrans)
#recompute(gs)
#add(wf, tf_s, parent="asinh") 

# need to apply the transformation first
#tf <- transformList(colnames(Data(wf[["base view"]]))[c(1:21)], asinh, transformationId="asinh")
#add(wf, tf)
#fs_trans <- transform(fs_comp, tf)
#gs <- GatingSet(fs_trans)

# boundary filter
#boundFilt <- boundaryFilter(filterId="boundFilt", x=c("SSC.A", "FSC.A"))
#add(wf, boundFilt, parent="asinh")
#gc()

# select non-NA markers for further analysis
wf <- workFlow(fs_filt, name="Twin_P1")
gc()

# apply biexponential transformation to SSC and FSC
biexpTrans <- flowJoTrans(channelRange=256, maxValue=261589, neg=0, pos=4.4176194777,
                          widthBasis=-100)
tf_rem <- transformList(colnames(fs_filt)[c(4:21)], biexpTrans,
                        transformationId="biexp")

tf <- estimateLogicle(fs_filt[[1]], channels=colnames(fs_comp[[1]]))
#tf <- transformList(colnames(fs_filt)[c(4:21)], logicleTransform, transformationId="logicle")

tf_as <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)], asinh, 
                       transformationId="asinh")
add(wf, tf_as)
add(wf, tf, name="logicle")
add(wf, tf_rem)
# add(wf, tf_rem, parent="base view")
# add a boundary filter on forward and side scatter
boundFilt <- boundaryFilter(filterId="boundFilt", x=c("FSC-A", "SSC-A"))
# add(wf, boundFilt, parent="asinh")
add(wf, boundFilt, parent="asinh")

# # gate on lymphocytes - preselection of CD3 for T-cells
# lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC-A", "SSC-A"), preselection="V800-A",
#                 filterId="TCells", eval=F, scale=1)
# add(wf, lg$n2gate, parent="boundFilt+")

mark_vec <- filter_markers((pData(parameters(Data(wf[["asinh"]])[[1]]))[,"desc"])[c(4:21)])
# use warping algorithm to align fluorescence intensity peaks across measured markers
pars <- colnames(Data(wf[["base view"]]))[mark_vec]
norm <- normalization(normFunction= function(x, parameters, ...) warpSet(x, parameters, ...),
                      parameters=pars, normalizationId = "warping")

# norm_log <- normalization(normFunction= function(x, parameters, ...) warpSet(x, parameters, ...),
#                       parameters=pars, normalizationId = "log_warping")

# norm_as <- normalization(normFunction= function(x, parameters, ...) warpSet(x, parameters, ...),
#                          parameters=pars, normalizationId = "asinh_warping")

add(wf, norm, parent="logicle")
add(wf, norm, parent="boundFilt+")
add(wf, norm, parent="biexp")
add(wf, norm, parent="asinh")
# add(wf, norm_log, parent="logicle")
# add(wf, norm_as, parent="asinh")
gc()

# mindensity correctly identifies the top CD3 population
lg_mat <- matrix(c(0, Inf), ncol=1)
colnames(lg_mat) <- "R780-A"
lg <- rectangleGate(.gate=lg_mat, filterId="lymphs")
# lg <- mindensity(Data(wf[["warping"]])[[1]], channel="R780-A", positive=TRUE,
#                  filter_id="lymphs")
add(wf, lg, parent="warping")

# quad gate for CD4+s and CD8+s
tcell_mat <- matrix(c(0, 0), ncol=2)
colnames(tcell_mat) <- c("V585-A", "V605-A")
tcell_g <- quadGate(.gate=tcell_mat, filter_id="CD8CD4")
add(wf, tcell_g, parent="lymphs+")

# get all triboolean gates
params <- pData(parameters(Data(wf[["lymphs+"]])[[1]]))[,"name"]
param_names <- pData(parameters(Data(wf[["lymphs+"]])[[1]]))[,"desc"]
tribools <- make_booleans(Data(wf[["CD8-CD4+"]])[[1]], params, param_names)

# apply all gates
list_of_fanos <- get_frames(wf, "CD8-CD4+", get_fano, tribools)
list_of_means <- get_frames(wf, "CD8-CD4+", get_means, tribools)

# plot the gating strategy
lymph_plot <- xyplot(`CD3` ~ `FSC-A`, Data(wf[["asinh"]]), 
                     xbin=128, smooth=F, filter=lg)

cd4_plot <- xyplot(`CD4` ~ `CD8`, Data(wf[["lymphs+"]]), 
                   smooth=F, xbin=128, filter=tcell_g)

pre_warp <- densityplot(~`CD3`, Data(wf[["asinh"]]))
post_warp <- densityplot(~`CD3`, Data(wf[["warping"]]))

# example triboolean gate
# CD57, CD27, CD31
cd57_mat <- matrix(c(0, Inf), ncol=1)
colnames(cd57_mat) <- c("V705-A")
cd57_g <- rectangleGate(.gate=cd57_mat, filterId="CD57")

cd27_mat <- matrix(c(0, Inf), ncol=1)
colnames(cd27_mat) <- c("B515-A")
cd27_g <- rectangleGate(.gate=cd27_mat, filterId="CD27")

cd31_mat <- matrix(c(0, Inf), ncol=1)
colnames(cd31_mat) <- c("G780-A")
cd31_g <- rectangleGate(.gate=cd31_mat, filterId="CD31")


# make 2-way gates for plotting
cd57cd27_g <- cd57_g * cd27_g
cd57cd31 <- cd57_g * cd31_g
cd27_cd31 <- cd27_g * cd31_g

cd57_plot <- xyplot(`CD57` ~ `CD27`, Data(wf[["CD8-CD4+"]]), smooth=F, xbin=32,
                    filter=cd57cd27_g)
cd27_plot <- xyplot(`CD27` ~ `CD31`, Data(wf[["CD8-CD4+"]]), smooth=F, xbin=32,
                    filter=cd27_cd31)
cd31_plot <- xyplot(`CD31` ~ `CD57`, Data(wf[["CD8-CD4+"]]), smooth=F, xbin=32,
                    filter=cd57cd31)

trigate <- cd57_g * cd27_g * cd31_g
add(wf, trigate, parent="CD8-CD4+")

densityplot(~`CD57`, Data(wf[["defaultRectangleGate+"]]))

for(j in 1:length(list_of_means)){
  df <- list_of_means[[j]]
  df$subset <- "ZZFL"
  gate <- names(list_of_means)[j]
  filename <- paste("ZZFL-P1", sprintf("%s", gate),"mean_CD4_Tmem.tsv", sep="-")
  #write.table(df, row.names=T,
  #            file=filename)
  print(filename)
}

# # select CD3- non Tcells
# nonlg_mat <- matrix(c(0, 2500, 2500, 0, 0,0, 116, 116), ncol=2)
# colnames(nonlg_mat) <- c("SSC-A", "R780-A")
# nonlg_g <- polygonGate(.gate=nonlg_mat, filterId="nonTCells")
# add(wf, nonlg_g, parent="warping")
# 
# # filter out all non-CD16 CD56 cells
# nonnk_mat <- matrix(c(116, 256, 256, 50, 50, 84, 104, 116, 0, 0, 256, 256, 84, 84, 72, 0), ncol=2)
# colnames(nonnk_mat) <- c("V585-A", "V450-A")
# nk_g <- polygonGate(.gate=nonnk_mat, filterId="NKCells")
# add(wf, nk_g, parent="nonTCells+")
# 
# # select NK cell subsets - CD16-CD56hi, CD16hi CD56dim, CD16+CD56+
# cd16neg_mat <- matrix(c(76, 180, 180, 50, 50, 76, 128, 180, 200, 200, 100, 100), ncol=2)
# colnames(cd16neg_mat) <- c("V450-A", "V585-A")
# cd16neg_g <- polygonGate(.gate=cd16neg_mat, filterId="earlyNK")
# add(wf, cd16neg_g, parent="NKCells+")
# 
# cd16hi_mat <- matrix(c(100, 256, 256, 172, 100, 50, 50, 128, 72, 50), ncol=2)
# colnames(cd16hi_mat) <-c ("V450-A", "V585-A")
# cd16hi_g <- polygonGate(.gate=cd16hi_mat, filterId="terminalNK")
# add(wf, cd16hi_g, parent="NKCells+")
# 
# cd56pos_mat <- matrix(c(100, 156, 228, 256, 256, 180, 80, 100, 
#                         54, 76, 128, 150, 172, 172, 128, 50), ncol=2)
# colnames(cd56pos_mat) <- c("V450-A", "V585-A")
# cd56_g <- polygonGate(.gate=cd56pos_mat, filterId="matureNK")
# add(wf, cd56_g, parent="NKCells+")


# # select Vd1+ cells on CD3 and FSC-A
# lg_mat <- matrix(c(0, 250000, 250000, 0, 128, 128, 256, 256), ncol=2)
# lg_log <- matrix(c(3.5, 4.0, 4.5, 4.5, 3.5, 3.5, 2.0, 2.0, 3.0, 5.0, 5.0, 2.0 ), ncol=2)
# colnames(lg_mat) <- c("FSC-A", "V800-A")
# colnames(lg_log) <- c("FSC-A", "V800-A")

# lg <- polygonGate(.gate=lg_mat, filterId="TCells")
# lg_lg <- polygonGate(.gate=lg_log, filterId="logicleTCells")
# 
# # lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, 
# #                num_peaks=2, ref_peak=1, filter_id="TCells")
# 
# add(wf, lg, parent="warping")
# add(wf, lg_lg, parent="log_warping")

# # select all non-CD19 and non-CD20 expressing cells, everything else is a B cell
# non_b <- matrix(c(0, -Inf, 0, -Inf), ncol=2)
# colnames(non_b) <- c("V655-A", "V605-A")
# nonb_g <- rectangleGate(.gate=non_b, filterId="nonBcells")
# add(wf, nonb_g)
# 
# # set a tailgate for CD10+ cells as Immature B cells
# imm_g <- mindensity(fr=Data(wf[["nonBcells-"]])[[1]], channel="R660-A",
#                     filter_id="ImmatureB")
# add(wf, imm_g, parent="nonBcells-")
# 
# # gate on Naive and Memory B cells
# 
# naive_mat <- matrix(c(5, 5, -10, -10, -10, 15, 15, -10), ncol=2)
# colnames(naive_mat) <- c("G660-A", "V450-A")
# naive_g <- polygonGate(.gate=naive_mat, filterId="NaiveBcell")
# add(wf, naive_g, parent="ImmatureB-")
# 
# mem_mat <- matrix(c(6.11, 6.11, 15, 15, -10, 15, 15, -10), ncol=2)
# colnames(mem_mat) <- c("G660-A", "V450-A")
# mem_g <- polygonGate(.gate=mem_mat, filterId="MemBcell")
# add(wf, mem_g, parent="ImmatureB-")
# 
# # generate and add quadrant gates for T cell subsets
# tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")
# 
# # begin combination gating on remaining parameters
# pos_mat <- matrix(c(0, Inf))
# colnames(pos_mat) <- c("V705-A")
# cd57_g <- rectangleGate(.gate=pos_mat, filterId="CD57")
# colnames(pos_mat) <- c("R660-A")
# cd45ra_g <- rectangleGate(.gate=pos_mat, filterId="CD45RA")
# cd57cd45ra <- cd57_g & cd45ra_g
# 
# add(wf, cd57_g, parent="CD4+CD8-")
# add(wf, cd45ra_g, parent="CD4+CD8-")
# add(wf, cd57cd45ra, parent="CD4+CD8-")
# 
# # check > 100 cells per gate before calculating statistics


# function to calculate Fano factor, expression noise, across all individuals and markers in a flowSet
# object
# values into a matrix - test on CD3+CD4+ naive T cell subset
fano.mat <- fsApply(Data(wf[["CD57 and CD45RA+"]]), use.exprs=TRUE, FUN=get_fano)
mean.mat <- fsApply(Data(wf[["CD57 and CD45RA+"]]), use.exprs=TRUE, FUN=get_means)
fano.frame <- data.frame(fano.mat)
mean.frame <- data.frame(mean.mat)
colnames(fano.frame) <- parameters(Data(wf[["CD57 and CD45RA+"]])[[1]])$desc
colnames(mean.frame) <- parameters(Data(wf[["CD57 and CD45RA+"]])[[1]])$desc

fano.frame$twin.id <- rownames(fano.frame)
mean.frame$twin.id <- rownames(mean.frame)
View(fano.frame)
gc()

fano.frame$subset <- "ZZFSprt2"
write.table(fano.frame, row.names=T,
            file="/ifs/projects/proj052/flow_processing_tables/ZZFSprt2-panel3-Fano_Tcells.tsv",
            sep="\t")
