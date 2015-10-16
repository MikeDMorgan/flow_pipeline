#######################################################
# Flow cytometry analysis using Bioconductor packages #
#######################################################

set_marker_id <- function(flowset) {
  # add marker IDs to flow.data - no assumption about panel number
  # add FSC and SSC parameter "desc"
  for(i in seq_len(length(flowset))){
  simp.marker = strsplit(pData(parameters(flowset[[i]]))[,"desc"], fixed=T, split=" ")
  new.marker = as.vector(unlist(lapply(simp.marker, FUN = function(x) {paste0(x[1])})))
  pData(parameters(flowset[[i]]))[, "desc"] = new.marker
  }
  return(flowset)
}

test_valid_vec <- function(value){
  # this function is used to test the validity
  # of a numerical value, if it is NA, then it will
  # return 0
  tryCatch(
    if(is.na(value)){
      value <- 0.0
      return(value)
    },
    warning = function(w) {
      value = 0.0
      return(value)},
    error = function(e) {
      value = 0.0
      return(value)})
  return(value)}


get_fano <- function(flowset) {
  fano.vec <- vector(mode="numeric", length=dim(flowset)[2])
  cols <- colnames(flowset)
  # mean function may return a NaN, need to
  # handle these nicely
    for(i in seq_len(dim(flowset)[2])){
      # a negative mean expression value indicates a lack of expression
      # expresssion noise is meaningless in these situations
      # gene expression should be normalised to cell size
      # see Metzger et al (2015) Nature for details
      if(i != "FSC-A"){
        norm <- (flowset[,i])^2/log10(flowset[,"FSC-A"])^3
        mean.fs <- mean(norm)}
      else {
        mean.fs <- mean(flowset[,i])
      }
      val <- test_valid_vec(mean.fs)
      if(val <= 0.0){
        fano.vec[i] <- 0.0
      }
      else {
      sd.fs <- sd(flowset[,i])
      fano.fs <- sd.fs/mean.fs
      fano.vec[i] <- fano.fs
      }
    }
  return(fano.vec)
}


check_params <- function(flowset) {
  # go over each parameter and check they is a signal
  cols <- colnames(flowset)
  vector = numeric(0)
  c = 1
  for(i in seq_len(length(cols))){
    # negative or zero mean expression implies a lack of expression
    mean.exprs <- mean(flowset.exprs[,i])
    if(mean.exprs > 0.0){
      vector[c] <- i
      c <- c + 1
    }
  } 
  return(vector)
  }

get_means <- function(flowset) {
  mean.vec <- vector(mode="numeric", length=dim(flowset)[2])
  cols <- colnames(flowset)
  for(i in seq_len(dim(flowset)[2])){
    # a negative mean expression value indicates a lack of expression
    # mean expression is meaningless in these situations
    if(i != "FSC-A"){
      norm <- (flowset[,i])^2/log10(flowset[,"FSC-A"])^3
      mean.fs <- mean(norm)}
    else {
      mean.fs <- mean(flowset[,i])
    }
    val <- test_valid_vec(mean.fs)
    if(val <= 0){
      mean.vec[i] <- 0.0
    }
    else{
      mean.vec[i] <- mean.fs
    }
  }
  return(mean.vec)
}

filter_samples <- function(flowset, min_cells=10000) {
  extract = numeric(0)
  for(x in seq_len(length(flowset))){
    events <-dim(flowset[[x]])[1]
    if(events < min_cells){
      extract = c(extract, x)
    }
  }
  return(extract)
}

filter_markers <- function(markers){
  marker_vec = numeric(0)
  for(i in seq_len(length(markers))){
    pname <- markers[i]
    if(pname != "NA"){
      marker_vec = c(marker_vec, i+3)
    }
  }
  return(marker_vec)
}

fix_scatter_name <- function(flowset){
  for(i in seq_len(length(flowset))){
    for(x in seq_len(length(colnames(flowset[[i]])))){
      cname = colnames(flowset[[i]])[x]
      if(cname == "SSC-A" | cname == "FSC-A" | cname == "FSC-H"){
        pData(parameters(flowset[[i]]))[,"desc"][x] <- cname
      }
    }
  }
  return(flowset)
}


make_booleans <- function(flowset, params, pnames) {
  names(params) <- pnames
  gvec <- list()
  tribools <- list()
  for(i in seq_len(length(params) - 1)){
    # iterate over each non-NA marker
    # create a gate and add it to a vector of gates
    if(names(params)[i] != "NA"){
      prm <- params[i]
      nme <- names(params)[i]
      gmat <- matrix(c(0, Inf), ncol=1)
      colnames(gmat) <- prm
      g1 <- rectangleGate(.gate=gmat, filterId=nme)
      gvec[nme] <- g1
    }
  }
  gate_combs <- combn(names(gvec), m=3)
  for(x in 1:dim(gate_combs)[2]){
    trigates <- gvec[gate_combs[,x]]
    boolgate <- trigates[[1]] * trigates[[2]] * trigates[[3]]
    boolname <- paste(gate_combs[,x][1], gate_combs[,x][2], gate_combs[,x][3], sep=".")
    tribools[boolname] <- boolgate
    }
  return(tribools)
}

get_frames <- function(workspace, data_view, summary_func, gate_list) {
  all_frames <- list()
  for(k in 1:length(gate_list)){
    # name list according to gating markers
    gate_name <- names(gate_list)[[k]]
    trigate <- gate_list[[k]]
    flow_data <- Subset(Data(workspace[[data_view]]), trigate)
    summary_matrix <- fsApply(flow_data, use.exprs=TRUE, FUN=summary_func)

    # set all missing values to 0??
    summary_matrix[is.na(summary_matrix)] <- 0
    summary_matrix[!is.finite(summary_matrix)] <- 0

    summary_frame <- data.frame(summary_matrix)
    colnames(summary_frame) <- pData(parameters(flow_data[[1]]))[,"desc"][c(4:21)]
    gc()
    all_frames[[gate_name]] <- summary_frame
  }
  return(all_frames)
}


iter_filewrite <- function(array_list, path, table) {
  for(i in 1:length(array_list)){
    nme <- names(array_list)[i]
    full_path = paste(path, table, nme, sep="-")
    write.table(file=full_path, array_list[[i]], sep="\t")
    }
}

