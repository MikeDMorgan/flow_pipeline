import CGAT.Experiment as E
import rpy2.robjects as ro
from rpy2.robjects import r as R
from rpy2.robjects import rinterface
from rpy2.robjects import pandas2ri
import re
import pandas as pd
import numpy as np


def get_cd4_Tmem_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3", 
                        statistic="Fano", cell_subset="CD4_Tmem"):
    '''
    Retrieve Fano factors across CD4+ Tmems from panel 3
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')
    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="G610.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''qgate.cd4 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610.A", "V800.A"),'''
      '''filterId = "CD3CD4")''')
    R('''add(wf, qgate.cd4, parent="warping")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''rect.mat <- matrix(c(0, Inf, 0, Inf), ncol=2, '''
      '''dimnames=list(c("min", "max"), c("R710.A", "V655.A")))''')
    R('''mem.gate.cd4 = rectangleGate(filterId="CD4+CCR7CD45RA", .gate=rect.mat)''')

    R('''add(wf, mem.gate.cd4, parent="CD3+CD4+")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ memory T cell subset currently
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD4+CCR7CD45RA+"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD4+CCR7CD45RA+"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def get_cd4_naive_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3",
                         statistic="Fano", cell_subset="CD4_Tnaive"):
    '''
    Retrience Fano factos across CD4+ naive T cells from
    panel 3
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
    #R('''colnames(comp) <- gsub("", "-", colnames(comp))''')
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC-A", "SSC-A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC-A", "SSC-A"), preselection="G610-A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[m_vec]''')
    R('''print(pars)''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''qgate.cd4 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610-A", "V800-A"),'''
      '''filterId = "CD3CD4")''')
    R('''add(wf, qgate.cd4, parent="warping")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''rect.mat <- matrix(c(0, Inf, 0, Inf), ncol=2, '''
      '''dimnames=list(c("min", "max"), c("R710-A", "V655-A")))''')
    R('''mem.gate.cd4 = rectangleGate(filterId="CD4+CCR7CD45RA", .gate=rect.mat)''')

    R('''add(wf, mem.gate.cd4, parent="CD3+CD4+")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ memory T cell subset currently
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD4+CCR7CD45RA-"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD4+CCR7CD45RA-"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def get_cd8_Tmem_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3",
                        statistic="Fano", cell_subset="CD8_Tmem"):
    '''
    Retrieve Fano factors across CD8+ Tmems from panel 3
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="G610.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["TCells+"]]))[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD8+ CD3+ Tcells")
    R('''qgate.cd8 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610.A", "V585.A"),'''
      '''filterId = "CD3CD8")''')
    R('''add(wf, qgate.cd8, parent="warping")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''rect.mat <- matrix(c(0, Inf, 0, Inf), ncol=2, '''
      '''dimnames=list(c("min", "max"), c("R710.A", "V655.A")))''')
    R('''mem.gate.cd8 = rectangleGate(filterId="CD8+CCR7CD45RA", .gate=rect.mat)''')

    R('''add(wf, mem.gate.cd8, parent="CD3+CD8+")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ memory T cell subset currently
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD8+CCR7CD45RA+"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD8+CCR7CD45RA+"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))

def get_cd8_naive_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3",
                         statistic="Fano", cell_subset="CD8_Tnaive"):
    '''
    Retrieve Fano factors across CD8+ naive T cells from
    panel 3
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')
    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="G610.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    # only feature match on non-NA markers
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(fs_filt)[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD8+ CD3+ Tcells")
    R('''qgate.cd8 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610.A", "V585.A"),'''
      '''filterId = "CD3CD8")''')
    R('''add(wf, qgate.cd8, parent="warping")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''rect.mat <- matrix(c(0, Inf, 0, Inf), ncol=2, '''
      '''dimnames=list(c("min", "max"), c("R710.A", "V655.A")))''')
    R('''mem.gate.cd8 = rectangleGate(filterId="CD8+CCR7CD45RA", .gate=rect.mat)''')
    R('''add(wf, mem.gate.cd8, parent="CD3+CD8+")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ memory T cell subset currently
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD8+CCR7CD45RA-"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD8+CCR7CD45RA-"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def get_cd4_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                          statistic="Fano", cell_subset="CD4_Tcells"):
    '''
    Retrieve Fano factors across CD4+ T cells from panel 1
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="R780.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # add a quadrante gate to pull out CD4+, CD8+, DN and DP
    R('''qgate <- quadrantGate(Data(wf[["warping"]]), stains=c("V605.A", "V585.A"),'''
      '''filterId = "CD4CD8")''')
    R('''add(wf, qgate, parent="warping")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ T cell subset
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD4+CD8-"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD4+CD8-"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def get_cd8_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                          statistic="Fano", cell_subset="CD8_Tcells"):
    '''
    Retrieve Fano factors across CD8+ T cells from panel 1
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes R780-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="R780.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # add a quadrante gate to pull out CD4+, CD8+, DN and DP
    R('''qgate <- quadrantGate(Data(wf[["warping"]]), stains=c("V605.A", "V585.A"),'''
      '''filterId = "CD4CD8")''')
    R('''add(wf, qgate, parent="warping")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ T cell subset
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD4-CD8+"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD4-CD8+"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def get_dn_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                         statistic="Fano", cell_subset="DN_Tcells"):
    '''
    Retrieve Fano factors across DN T cells from panel 1
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="R780.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # add a quadrante gate to pull out CD4+, CD8+, DN and DP
    R('''qgate <- quadrantGate(Data(wf[["warping"]]), stains=c("V605.A", "V585.A"),'''
      '''filterId = "CD4CD8")''')
    R('''add(wf, qgate, parent="warping")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ T cell subset
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD4-CD8-"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD4-CD8-"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def get_dp_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                         statistic="Fano", cell_subset="DP_Tcells"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="R780.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # add a quadrante gate to pull out CD4+, CD8+, DN and DP
    R('''qgate <- quadrantGate(Data(wf[["warping"]]), stains=c("V605.A", "V585.A"),'''
      '''filterId = "CD4CD8")''')
    R('''add(wf, qgate, parent="warping")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    # put values into matrix, just the CD3+CD4+ T cell subset
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD4+CD8+"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD4+CD8+"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def get_noise(fcs_dir, out_dir, comp_matrix, panel, setid,
              statistic="Fano", cell_subset="CD4_Tmem"):
    '''
    Pull out the fluorescence intensities based on a specific
    gating strategy and calculate the required measure of noise.
    Write to a tab-delimited file.
    
    fcs_dir - should only contain .fcs files from a single subset.
    panel -  should be in the form Pn, where n=[1, ..., 7]
    setid - four character ID taken from the zip file, used for output file
    statistic - the statistic used to represent expression noise.  Currently
    only implements Fano factor (sigma^2/mu).  TODO: include sigma^2, mu ~ sigma^2,
    sigma, sigma^2 ~ mu.
x
    TODO: 
    * make this generic for different cell types and marker panels
    * split into different functions for different cell types?
    * accept specific gating strategy as a parameter
    * wrap up into proper python objects
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE)''')

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''Sys.setlocale('LC_ALL', 'C')''')
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    #R('''for(i in seq_len(length(flow.data))){''')
    #R('''simp.marker <- strsplit(pData(parameters(flow.data[[i]]))[, "desc"],'''
    #  '''fixed=T, split=" ")''')
    #R('''new.marker <- as.vector(unlist(lapply(simp.marker, '''
    #  '''FUN=function(x) {paste0(x[1])})))''')
    #R('''pData(parameters(flow.data[[i]]))[, "desc"] <- new.marker}''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(flow.data, name="Twins_%s")''' % panel)

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[3:21],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC-A", "SSC-A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC-A", "SSC-A"), preselection="G610-A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[c(4, 6, 9, 10, 11, '''
      '''13, 14, 15, 17, 18, 19, 21)]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ and CD8+ CD3+ Tcells")
    R('''qgate.cd4 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610-A", "V800-A"),'''
      '''filterId = "CD3CD4")''')
    R('''add(wf, qgate.cd4, parent="warping")''')

    R('''qgate.cd8 <- quadrantGate(Data(wf[["warping"]]), stains=c("G610-A", "V585-A"),'''
      '''filterId = "CD3CD8")''')
    R('''add(wf, qgate.cd8, parent="warping")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''rect.mat <- matrix(c(0, Inf, 0, Inf), ncol=2, '''
      '''dimnames=list(c("min", "max"), c("R710-A", "V655-A")))''')
    R('''mem.gate.cd4 = rectangleGate(filterId="CD4+CCR7CD45RA", .gate=rect.mat)''')
    R('''mem.gate.cd8 = rectangleGate(filterId="CD8+CCR7CD45RA", .gate=rect.mat)''')

    R('''add(wf, mem.gate.cd4, parent="CD3+CD4+")''')
    R('''add(wf, mem.gate.cd8, parent="CD3+CD8+")''')

    # function to calculate Fano factor, expression noise, across all individuals
    # and markers in a flowSet object
    # stick this in an .R file and source it

    #R('''get_fano <- function(flowset) {''')
    #R('''fano.vec = vector(mode="numeric", length=dim(flowset)[2])''')
    #R('''cols = colnames(flowset)''')
    #R('''for(i in seq_len(dim(flowset)[2])) {''')
    #R('''mean.fs = mean(flowset[, i])''')
    #R('''var.fs = var(flowset[, i])''')
    #R('''fano.fs = var.fs/mean.fs''')
    #R('''fano.vec[i] <- fano.fs}''')
    #R('''return(fano.vec)}''')

    # put values into matrix, just the CD3+CD4+ memory T cell subset currently
    E.info("calculating Fano factor for all fluorescence parameters")
    R('''fano.mat <- fsApply(Data(wf[["CD4+CCR7CD45RA+"]]), use.exprs=TRUE,'''
      '''FUN=get_fano)''')
    R('''fano.frame <- data.frame(fano.mat)''')
    R('''colnames(fano.frame) <- parameters(Data(wf[["CD4+CCR7CD45RA+"]])[[1]])$desc''')
    R('''fano.frame$twin.id <- rownames(fano.frame)''')
    R('''gc()''')

    R('''fano.frame$subset <- "%s"''' % setid)
    
    R('''sink(file=NULL)''')

    E.info("writing output file")
    R('''write.table(fano.frame, row.names=T,'''
      '''file="%s/%s-%s-%s_%s.tsv", sep="\t")''' % (out_dir, setid, panel,
                                                    statistic, cell_subset))


def merge_flow_tables(file_list, id_column,
                      demo_file, demo_id_column):
    '''
    Take multiple files output from calculating noise measure,
    append together and merge with demographics file that
    contains zygosity and twin pair information

    id_column and demo_id_column can be string or integer.  If integer
    they refer to the respective column numbers in file_list and demo_file
    that correspond to the Twin flow cytometry IDs.  If they are strings
    then they refer explicitly to the column header for each.
    '''

    tab1 = file_list[0]
    df = pd.read_table(tab1, sep="\t", header=0, index_col=None)
    file_list.remove(tab1)

    for tabs in file_list:
        _df = pd.read_table(tabs, sep="\t", header=0, index_col=None)
        df = df.append(_df)

    # get twin IDs in flow data to match those in demographics data
    if type(id_column) == int:
        twin_id_col = df.columns[id_column]
        df["twin_num"] = [int(xi.split("_")[-1]) for xi in df[twin_id_col].values]
    elif type(id_column) == str:
        df["twin_num"] = [int(xs.split("_")[-1]) for xs in df[id_column].values]
    else:
        raise ValueError("Twin ID columns must be either column number or header")

    demo_df = pd.read_table(demo_file, sep="\t", header=0, index_col=None)
    # rename demo_df columns
    demo_cols = ["flowjo_id", "family_id", "zygosity",
                 "age", "replicate", "visit"]
    demo_df.columns = demo_cols

    # get the column header on which to merge with the flow data
    if type(demo_id_column) == int:
        demo_id_col = demo_df.columns[demo_id_column]
    elif type(demo_id_column) == str:
        demo_id_col = "flowjo_id"
    else:
        raise ValueError("Twin ID columns must be a header name "
                         "or a column numnber")
    merge_df = pd.merge(left=df, right=demo_df, left_on="twin_num",
                        right_on=demo_id_col, how="inner")
    merge_indx = [ix for ix, iy in enumerate(merge_df.index)]
    merge_df.index = merge_indx

    return merge_df


def split_zygosity(infile, zygosity_column, id_headers,
                   pair_header):
    '''
    Split the twins file by their zygosity and match up appropriate twin
    pairs based on their family ID
    '''

    twins_df = pd.read_table(infile, sep="\t", header=0, index_col=0)
    MZ_df = twins_df[twins_df[zygosity_column] == "MZ"]
    DZ_df = twins_df[twins_df[zygosity_column] == "DZ"]

    # split each twin pair and match up on family ID
    # TODO: make report/output/log of twins not properly paired
    MZ_1 = MZ_df.loc[MZ_df.duplicated(subset=pair_header)]
    MZ_2 = MZ_df.loc[[not im for im in MZ_df.duplicated(subset=pair_header)]]

    DZ_1 = DZ_df.loc[DZ_df.duplicated(subset=pair_header)]
    DZ_2 = DZ_df.loc[[not di for di in DZ_df.duplicated(subset=pair_header)]]

    # melt the dataframes in order to merge twin pairs
    MZ1_melt = pd.melt(MZ_1, id_vars=id_headers, var_name="marker",
                       value_name="twin1")
    MZ2_melt = pd.melt(MZ_2, id_vars=id_headers, var_name="marker",
                       value_name="twin2")
    matched_MZ = pd.merge(left=MZ1_melt, right=MZ2_melt,
                          on=["marker", zygosity_column, "subset", pair_header],
                          how='inner')

    DZ1_melt = pd.melt(DZ_1, id_vars=id_headers, var_name="marker",
                       value_name="twin1")
    DZ2_melt = pd.melt(DZ_2, id_vars=id_headers, var_name="marker",
                       value_name="twin2")

    matched_DZ = pd.merge(left=DZ1_melt, right=DZ2_melt,
                          on=["marker", zygosity_column, "subset", pair_header],
                          how='inner')

    return {"MZ": matched_MZ, "DZ": matched_DZ}


def regress_out_confounding(infile, confounding_column, group_var):
    '''
    Regress out confounding effects such as process batch,
    age, etc.

    confounding_column - comma-separate list of column headers or 
    column numbers to include in linear model

    Output are regression residuals.  Be careful to check that these
    are heteroscedastic and unstructuted (except for 
    relationship between twins).
    '''

    # need to use pandas2ri.activate() to get conversion between
    # python and rpy - why do they keep changing this?????
    pandas2ri.activate()

    twin_frame = pd.read_table(infile, sep="\t", header=0, index_col=0)
    marker_groups = twin_frame.groupby(group_var)

    # set up rpy function for linear model fit
    r_lm = R["lm"]

    frame_list = []
    if type(confounding_column) == int:
        confound_vars = twin_frame.columns[confounding_column]
    elif type(confounding_column) == str:
        confound_vars = [confounding_column]
    else:
        raise TypeError("type not recognised, enter either a "
                        "column number or header")

    confounding = " + ".join(confound_vars)

    for name, group in marker_groups:
        # convert to an rpy dataframe
        group_df = group.copy()
        # set NAs to 0
        group_df = group_df.fillna(0.0)

        r_group = pandas2ri.py2ri_pandasdataframe(group_df)

        # run linear model, adjusting for confounding variables
        twin1_r = r_lm("twin1 ~ age_x + %s " % confounding, data=r_group)
        twin2_r = r_lm("twin2 ~ age_y + %s " % confounding, data=r_group)

        # pull out residuals
        twin1_residuals = twin1_r.rx("residuals")[0]
        twin2_residuals = twin2_r.rx("residuals")[0]

        # pull out indicies and convert back to python objects,
        # make sure index values are integers
        tw1_indx = [int(ix1) for ix1 in twin1_residuals.dimnames[0]]
        tw2_indx = [int(ix2) for ix2 in twin2_residuals.dimnames[0]]

        tw1_rs = [fi1 for fi1 in twin1_residuals]
        tw2_rs = [fi2 for fi2 in twin2_residuals]

        # transform into Z-scores

        tw1_z = [(fz1 - np.mean(tw1_rs))/np.std(tw1_rs) for fz1 in tw1_rs]
        tw2_z = [(fz2 - np.mean(tw2_rs))/np.std(tw2_rs) for fz2 in tw2_rs]

        tw1_ser = pd.Series(tw1_rs, index=group_df.index, dtype=np.float64)
        tw2_ser = pd.Series(tw2_rs, index=group_df.index, dtype=np.float64)

        res_df = pd.DataFrame([tw1_ser, tw2_ser]).T
        res_df.columns = ["twin1_res", "twin2_res"]

        # merge dataframe of residual and original dataframe
        m_df = pd.merge(left=group_df, right=res_df,
                        left_index=True, right_index=True,
                        how='inner')

        frame_list.append(m_df)

    df = frame_list[0]
    frame_list.remove(df)

    for each in frame_list:
        _df = each
        df = df.append(_df)

    df.index = [xi for xi, yi in enumerate(df.index)]
    return df


def estimate_heritability(mz_file, dz_file):
    '''
    Estimate heritability from twins using Falconers equation:
    2 * intra-twin class difference

    mz_file - text file containing data on MZ twins
    dz_file - text file containing data on DZ twins

    '''

    dz_frame = pd.read_table(dz_file, sep="\t", header=0,
                             index_col=0)
    mz_frame = pd.read_table(mz_file, sep="\t", header=0,
                             index_col=0)    

    DZ_groups = dz_frame.groupby(["marker"])
    MZ_groups = mz_frame.groupby(["marker"])

    # iterate over markers and calculate intra-twin class correlations
    MZ_corr_dict = {}
    DZ_corr_dict = {}

    for dname, dgroup in DZ_groups:
        dztwin_cor = np.corrcoef(dgroup["twin1"], dgroup["twin2"])
        DZ_corr_dict[dname] = abs(dztwin_cor[0, 1])

    for mname, mgroup in MZ_groups:
        mztwin_cor = np.corrcoef(mgroup["twin1"], mgroup["twin2"])
        MZ_corr_dict[mname] = abs(mztwin_cor[0, 1])

    # can only calculate the broad-sense heritability using Falconer's
    # equation as additive and dominance variance components cannot
    # be separated due to close familial relationship
    # these heritability estimates are the upper bounds of heritability
    broad_h = {}

    for key in MZ_corr_dict.keys():
        mz_cor = MZ_corr_dict[key]
        dz_cor = DZ_corr_dict[key]
        H2 = 2 * (mz_cor - dz_cor)
        broad_h[key] = H2

    return broad_h

def get_compensation_matrix(path, infile):
    '''
    Pull out the compensation matrix from an XML workspace file
    '''

    # store dimension: value pairs for each corresponding dimension
    param_dict = {}
    with open(path + "/" + infile, "r") as ofile:
        for line in ofile:
            if re.search("transforms:value", line):
                elements = line.split(" ")
                # get the parameter/dimension
                pars = [par.split("=")[1] for par in elements if re.search("parameter", par)][0]
                pars = pars.strip('"')
                # get the compensation value
                value = [val.split("=")[1] for val in elements if re.search("value", val)][0]
                value = float(value.strip('"'))
                val_dict[pars] = value
                param_dict[param] = val_dict
            elif re.search("transforms:spillover", line):
                if re.search("data-type:parameter", line):
                    # pull out parameter/dimension compensating against
                    l_split = line.split(" ")
                    param = [xv.split("=")[1] for xv in l_split if re.search("parameter", xv)][0]
                    param = param.strip('"')
                    val_dict = {}
                    param_dict[param] = []
                else:
                    pass
            else:
                pass
            
        comp_matrix = pd.DataFrame(param_dict)

    return comp_matrix


def parse_gating_file(infile):
    '''
    Parse a gating design file into multiple dummy files
    '''

    gating_df = pd.read_table(infile, sep="\t", header=0, index_col=None)
    try:
        assert np.any(gating_df["panel"])
        assert np.any(gating_df["cell_type"])
        assert np.any(gating_df["gating"])
    except KeyError:
        raise AttributeError("file header missing, please use the "
                             "columns panel, cell_type and gating")

    file_list = []
    for row in gating_df.index:
        panel = gating_df.loc[row]["panel"]
        if type(panel) != str:
            panel = "P" + str(panel)
        else:
            pass
        cell_type = gating_df.loc[row]["cell_type"]
        out_file = "-".join([panel, cell_type, "gating"]) + ".txt"
        file_list.append(out_file)
    return file_list


def gating_strategy(fcs_dir, out_dir, comp_matrix, panel, setid,
                    gating, statistic="Fano"):
    '''
    Use the gating strategy derived from an XML .wsp file
    and apply with a set of gates to .fcs files
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # set up gating strategy
    E.info("creating gating strategy")
    R('''gs <- GatingSet(fs_comp)''')

    E.info("applying biexponential transformation to data")
    R('''biexpTrans <- flowJoTrans(channelRaneg=256, maxValue=261589,'''
      '''neg=0, pos=4.4176194777, widthBasis=-100)''')
    R('''tf <- transformList(colnames(fs_comp), biexpTrans)''')

    # get polygon vertex co-ordinates from XML file
    # set up n-dimensional matrix and create polygon gate
    # get this from python as an np.array
    E.info("gating around lymphocytes")
    R('''lmat <- matrix(c(915.8591, 1593.8391, 1976.3916, 2204.0872,'''
      '''14.925808, -37.386353, 155.35947, 209920, 145408, 93184,'''
      '''33792, 36864, 151552, 197632), ncol=2)''')
    R('''colnames(lmat) <- c("SSC.A", "FSC.A")''')
    R('''lymphgate <- polygonGate(lmat)''')
    R('''add(gs, lymphgate, parent="root", name="lymphs")''')
    R('''recompute(gs)''')


def get_early_nktcells_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                              statistic="Fano", cell_subset="NKT_early"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/R_scripts/FlowProcess.R")''')
    R('''Sys.setlocale('LC_ALL', 'C')''')

    E.info("reading input .fcs files from %s " % fcs_dir)
    R('''path.to.files = "%s"''' % fcs_dir)
    R('''flowdata = read.flowSet(path = path.to.files,'''
      '''transformation=FALSE, pattern="%s")''' % panel)

    # reassign sample names, retaining original ID as part
    # have to set system locale to recognise non-UTF8 encoded
    # characters (silly file naming from FlowRepository)
    E.info("setting sample names")
    R('''twin.split = strsplit(x=sampleNames(flowdata), split=" ", fixed=T)''')
    R('''twin.ids = unlist(lapply(twin.split, FUN=function(x) {paste("Twin", '''
      '''x[1], x[2], sep="_")}))''')
    R('''sampleNames(flowdata) <- twin.ids''')

    E.info("setting marker IDs")
    # add the marker IDs to flow.data with assumptions of marker panel number
    R('''flow.data <- set_marker_id(flowdata)''')
    R('''colnames(flow.data) <- gsub("-", ".", colnames(flow.data))''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''fs_comp <- compensate(x=flow.data, spillover=comp)''')

    # some samples have very few events which breaks the warping function
    # remove all samples with < 10000 events
    E.info("filtering out samples with insufficient observations")
    R('''extrct <- filter_samples(fs_comp)''')
    R('''fs_filt <- fs_comp[seq(along=fs_comp)]''')
    R('''if(length(extrct)){fs_filt <- fs_filt[-(extrct),]}''')
    R('''rm(list=c("fs_comp"))''')
    R('''gc()''')

    E.info("constructing a workflow object")
    # set up the workflow object
    R('''wf <- workFlow(fs_filt, name="Twins_%s")''' % panel)
    R('''gc()''')

    # transform onto an approx. log-linear scale with asinh
    # TODO: alter this to a user specific transformation
    # exclude SSC and FSC parameters (indices 1-3)
    # assumes 21 parameters in total - different for each panel???

    E.info("transforming data")
    R('''tf <- transformList(colnames(Data(wf[["base view"]]))[c(4:21)],'''
      '''asinh, transformationId="asinh")''')

    # do some garbage collecting inside the R env
    R('''gc()''')

    R('''add(wf, tf)''')

    E.info("setting boundary filter on forward and side scatter")
    # add a boundary filter on forward and side scatter
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC.A", "SSC.A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    R('''lg <- lymphGate(Data(wf[["boundFilt+"]]), channels=c("FSC.A", "SSC.A"), preselection="V800.A",'''
      '''filterId="TCells", eval=F, scale=2.5)''')
    R('''add(wf, lg$n2gate, parent="boundFilt+")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''m_vec <- filter_markers(fs_filt)''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[m_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="TCells+")''')
    R('''gc()''')

    
