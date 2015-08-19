#######################################################
# Flow cytometry analysis using Bioconductor packages #
#######################################################

set_marker_id <- function(flowset) {
  # add marker IDs to flow.data - no assumption about panel number
  # add FSC and SSC parameter "desc"
  for(i in seq_len(length(flowset))){
    for(x in seq_len(length(colnames(flowset[[i]])))){
      cname <- colnames(flowset[[i]])[x]
      if(cname == "SSC-A" | cname == "FSC-A" | cname == "FSC-H"){
        pData(parameters(flowset[[i]]))[,"desc"][x] <- cname}
    }
    simp.marker = strsplit(pData(parameters(flowset[[i]]))[,"desc"], fixed=T, split=" ")
    new.marker = as.vector(unlist(lapply(simp.marker, FUN = function(x) {paste0(x[1])})))
    pData(parameters(flowset[[i]]))[, "desc"] = new.marker
  }
  return(flowset)
}

get_fano <- function(flowset) {
  fano.vec <- vector(mode="numeric", length=dim(flowset)[2])
  cols <- colnames(flowset)
  for(i in seq_len(dim(flowset)[2])){
    mean.fs = mean(flowset[,i])
    var.fs = var(flowset[,i])
    fano.fs = var.fs/mean.fs
    fano.vec[i] = fano.fs
  }
  return(fano.vec)
}

filter_samples <- function(flowset) {
  extract = numeric(0)
  for(x in seq_len(length(flowset))){
    events <- as.numeric(keyword(flowset[[x]])[["$TOT"]])
    if(events < 10000){
      extract = c(extract, x)
    }
  }
  return(extract)
}

filter_markers <- function(flowset){
  marker_vec = numeric(0)
  for(i in seq_len(length(pData(parameters(flowset[[1]]))[["desc"]]))){
    pname <- pData(parameters(flowset[[1]]))[["desc"]][i]
    if(pname != "NA" || pname != "FSC-A"){
      marker_vec = c(marker_vec, i)
    }
  }
  return(marker_vec)
}
