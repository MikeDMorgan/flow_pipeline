import CGAT.Experiment as E
import rpy2.robjects as ro
from rpy2.robjects import r as R
from rpy2.robjects import rinterface
from rpy2.robjects import pandas2ri
import re
import pandas as pd
import numpy as np
import math
import itertools
import os
import threading
import sqlite3 as sql
import pandas.io.sql as pdsql
import gc


# MM 14/09/2015 - change to biexponential transformation instead of asinh
# problems encountered with some data, i.e CD4 in panel 5, subset data ZZFY
# try re-downloading?
# why would using a biexponential transformation get around this??
# lots of 0 values??

def get_cd4_Tmem_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3", 
                        cell_subset="CD4_Tmem", db="csvdb"):
    '''
    Retrieve Fano factors across CD4+ Tmems from panel 3
    '''

    pandas2ri.activate()

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')

    # this is to get around annoying file names with stupid non-ASCII characters
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="G610-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ Tcells")
    R('''tcells <- quadGate("V800-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''tmem_mat <- matrix(c(0, Inf, 0, Inf), ncol=2)''')
    R('''colnames(tmem_mat) <- c("V655-A", "R710-A")''')
    R('''tmem_g <- rectangleGate(.gate=tmem_mat, filterId="CCR7CD45RA")''')
    R('''add(wf, tmem_g, parent="CD4+CD8-")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CCR7CD45RA+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CCR7CD45RA+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CCR7CD45RA+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset these to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CCR7CD45RA+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CCR7CD45RA+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    py_con = setupRSqlConnection(db)
    # table_vars = [("gate", "TEXT"), ("batch", "TEXT"), ("rows", "TEXT"),
    #               ("columns", "TEXT"), ("array","TEXT"), ("panel", "TEXT")]
    # key_vars = ["gate", "batch", "panel"]

    # makeSqlRTable(py_con, "%s_fano" % cell_subset,
    #               variable_type=table_vars, key_vars=key_vars)

    # makeSqlRTable(py_con, "%s_mean" % cell_subset,
    #               variable_type=table_vars, key_vars=key_vars)

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    # insert function that puts all arrays into a single table

    E.info("creating SQLite table %s_fano" % cell_subset)
    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    E.info("creating SQLite table %s_mean" % cell_subset)
    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    E.info("serializing arrays to text format")
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''cols <- colnames(fanos_df[,-1])''')
    df_headers = [x for x in R["cols"]]
    variables = zip(df_headers, df_headers)
    fano_query = constructRQuery(table_name="%s_fano" % cell_subset,
                                 variables=variables)
    mean_query = constructRQuery(table_name="%s_mean" % cell_subset,
                                 variables=variables)
    py_fanos = R["fanos_df"]
    py_means = R["means_df"]

    E.info("inserting records into tables")

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())

    fano_query = R["fano_query"]
    loadRSqlDb(r_connect=py_con,
               R_obj=py_fanos,
               R_query=fano_query)

    mean_query = R["mean_query"]
    loadRSqlDb(r_connect=py_con,
               R_obj=py_means,
               R_query=mean_query)

    # try:
    #     # R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
    #     R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
    #       '''(array, gate, batch, rows, columns, panel) values '''
    #       '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    # except rinterface.RRuntimeError:
    #     # R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
    #     R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
    #       '''(array, gate, batch, rows, columns, panel) values '''
    #       '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    # try:
    #     # R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
    #     R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
    #       '''(array, gate, batch, rows, columns, panel) values '''
    #       '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    # except rinterface.RRuntimeError:
    #     # R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
    #     R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
    #       '''(array, gate, batch, rows, columns, panel) values '''
    #       '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    R('''dbDisconnect(con)''')

    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                        cell_subset)
    E.info("output file %s" % outfile)
    os.system("touch %s" % outfile)

    R('''sink(file=NULL)''')
    

def get_cd4_naive_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3",
                         cell_subset="CD4_Tnaive", db="csvdb"):
    '''
    Retrience Fano factos across CD4+ naive T cells from
    panel 3
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="G610-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V800-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''tmem_mat <- matrix(c(0, Inf, 0, Inf), ncol=2)''')
    R('''colnames(tmem_mat) <- c("V655-A", "R710-A")''')
    R('''tmem_g <- rectangleGate(.gate=tmem_mat, filterId="CCR7CD45RA")''')
    R('''add(wf, tmem_g, parent="CD4+CD8-")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CCR7CD45RA-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CCR7CD45RA-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CCR7CD45RA-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CCR7CD45RA-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CCR7CD45RA-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, arrays TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, arrays TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')
    

def get_cd8_Tmem_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3",
                        cell_subset="CD8_Tmem", db="csvdb"):
    '''
    Retrieve Fano factors across CD8+ Tmems from panel 3
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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
    R('''boundFilt = boundaryFilter(filterId="boundFilt", x=c("FSC-A", "SSC-A"))''')
    R('''add(wf, boundFilt, parent="asinh")''')

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="G610-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V800-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''tmem_mat <- matrix(c(0, Inf, 0, Inf), ncol=2)''')
    R('''colnames(tmem_mat) <- c("V655-A", "R710-A")''')
    R('''tmem_g <- rectangleGate(.gate=tmem_mat, filterId="CCR7CD45RA")''')
    R('''add(wf, tmem_g, parent="CD4-CD8+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CCR7CD45RA+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CCR7CD45RA+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CCR7CD45RA+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CCR7CD45RA+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CCR7CD45RA+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # database locking may be an issue with multiple processes trying
    # to write to the database

    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)

    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd8_naive_panel3(fcs_dir, out_dir, setid, comp_matrix, panel="P3",
                         cell_subset="CD8_Tnaive", db="csvdb"):
    '''
    Retrieve Fano factors across CD8+ naive T cells from
    panel 3
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # gate on lymphocytes - preselection on CD3 for T-cells
    # assumes G610-A refers to the CD3 parameter - only applies to panel 3!!
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="G610-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V800-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # add separate gates for memory and naive cells based on CCR7 and CD45RA markers
    # generate a matrix of limits for 2D rectangle gate
    E.info("Selecting memory T cell population")
    R('''tmem_mat <- matrix(c(0, Inf, 0, Inf), ncol=2)''')
    R('''colnames(tmem_mat) <- c("V655-A", "R710-A")''')
    R('''tmem_g <- rectangleGate(.gate=tmem_mat, filterId="CCR7CD45RA")''')
    R('''add(wf, tmem_g, parent="CD4-CD8+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CCR7CD45RA-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CCR7CD45RA-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CCR7CD45RA-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CCR7CD45RA-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CCR7CD45RA-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd4_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                          cell_subset="CD4_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across CD4+ T cells from panel 1
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="R780-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4+CD8-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4+CD8-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4+CD8-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4+CD8-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4+CD8-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd4_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2a",
                          cell_subset="CD4_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4+CD8-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4+CD8-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4+CD8-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4+CD8-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4+CD8-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd4_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2b",
                          cell_subset="CD4_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4+CD8-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4+CD8-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4+CD8-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4+CD8-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4+CD8-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd8_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2b",
                          cell_subset="CD8_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd8_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2a",
                          cell_subset="CD8_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd8_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                          cell_subset="CD8_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across CD8+ T cells from panel 1
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="R780-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_dn_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                         cell_subset="DN_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DN T cells from panel 1
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="R780-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')


    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_dp_tcells_panel1(fcs_dir, out_dir, setid, comp_matrix, panel="P1",
                         cell_subset="DP_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="R780-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4+CD8+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4+CD8+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4+CD8+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4+CD8+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4+CD8+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_dp_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2b",
                         cell_subset="DP_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4+CD8+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4+CD8+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4+CD8+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4+CD8+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4+CD8+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_dp_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2a",
                         cell_subset="DP_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4+CD8+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4+CD8+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4+CD8+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4+CD8+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4+CD8+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_dn_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2b",
                         cell_subset="DN_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')


    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_dn_tcells_panel2(fcs_dir, out_dir, setid, comp_matrix, panel="P2a",
                         cell_subset="DN_Tcells", db="csvdb"):
    '''
    Retrieve Fano factors across DP T cells from panel 
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    R('''lg <- mindensity(Data(wf[["warping"]])[[1]], channel="V800-A",'''
      '''positive=TRUE, filter_id="lymphs")''')
    R('''add(wf, lg, parent="warping")''')

    # generate and add quandrant gates for CD3+CD4+ and CD3+CD8+ subsets
    # this is not generic - need to change this function to accomodate different
    # cell subsets
    E.info("gating around CD4+ CD3+ Tcells")
    R('''tcells <- quadGate("V605-A"=0, "V585-A"=0, filterId="CD4CD8")''')
    R('''add(wf, tcells, parent="lymphs+")''')


    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(arrays, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(arrays, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(arrays, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(arrays, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(arrays, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def synchronized(func):
    '''
    A decorator that mimics the Java synchronized primitive.  It
    adds a lock to the decorated function

    code from http://theorangeduck.com/page/synchronized-python
    '''

    func.__lock__ = threading.Lock()

    def sync_func(*args, **kwargs):
        with func.__lock__:
            return func(*args, **kwargs)

    return sync_func


def setupRSqlConnection(db_path):
    '''
    Setup an SQlite connection in R using RSQLite
    package.

    Arguments
    ---------
    db_path: str
      an absolute path pointing to an SQLite database

    Returns
    -------
    r_con: rpy2.robjects.method.RS4
      an rpy2 object representing an SQLite connection in R
    '''

    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''con <- dbConnect(SQLite(), "%s")''' % db_path)
    py_sqlcon = R["con"]
    E.info("connected to SQLite database %s" % db_path)
    return py_sqlcon


def makeSqlRTable(r_connect, table_name, variable_type,
                  key_vars=None):
    '''
    Create a new SQLite table in a database using RSQLite

    Arguments
    ---------
    r_connect: rpy2.robjects.method.RS4
      an rpy2 object representing an SQLite connection in R

    table_name: str
      the name of the table to create in the SQLite DB

    variable_type: list
      a list of tuples, with (var, TYPE) pairs as elements

    key_vars: list
      a list of variables to act as unique identifiers for records
      and primary key in SQLite table

    Returns
    -------
    None
    '''

    R('''require(RSQLite)''')
    R.assign("con", r_connect)
    variables = ", ".join([" ".join([x, y]) for x, y in variable_type])
    keys = ", ".join(key_vars)

    E.info("creating table %s" % table_name)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(table_name)s '''
      '''(%(variables)s, PRIMARY KEY (%(keys)s))") ''' % locals())


def constructRQuery(table_name, variables):
    '''
    Construct an SQL INSERT query

    Arguments
    ---------
    table_name: str
      the name of the table to insert records int

    variables: list
      a list of tuples with (variable, df header) pairs

    Returns
    -------
    query: str
      the SQL query
    '''

    table_headers = ", ".join([x for x, y in variables])
    df_headers = ", ".join([":" + y for x, y in variables])

    query = '''insert OR ignore into %(table_name)s (%(table_headers)s)
    values (%(df_headers)s)''' % locals()

    return query


@synchronized
def loadRSqlDb(r_connect, R_obj, R_query):
    '''
    Load an R dataframe object into an SQLite database,
    allowing for SQL DB locking using R libraries.

    Arguments
    ---------
    db_path: str
      an absolute path to an SQLite database on the file system

    R_obj: rpy2.robject
      an rpy2.robject that holds an R dataframe to be loaded

    R_query: string
      an SQL query

    Returns
    -------
    value: boolean
      True if successful, False if exception raised
    '''

    R('''require(RSQLite)''')
    R.assign("con", r_connect)
    R.assign("df", R_obj)

    # order of function:
    # need to setup the db connection within R
    # create the table using RSQlite
    # insert records into table with RSQlite
    # deal with database locking issues using "with"?

    try:
        R('''dbGetPreparedQuery(con, "%(R_query)s", bind.data=df)''' % locals())
    except rinterface.RRuntimeError:
        E.info("database locked")
        return False

    return True


def getRSqlTable(r_connect, query):
    '''
    Get records from an SQLite DB using RSQLite

    Arguments
    ---------
    r_connect: rpy2.robjects.method.RS4
      a connection to an SQLite database in R

    query: string
      the query used to recover records from the SQLite DB

    Returns
    -------
    df: rpy2.robjects.vectors.DataFrame
      a dataframe of arrays corresponding to a query
    '''

    R.assign("con", r_connect)
    R('''r.df <- dbGetQuery(con, "%(query)s")''' % locals())
    df = R["r.df"]

    return df


def clearREnvironment():
    '''
    Call this function to clear out all elements of the current
    active R environment.  Retains attached packages
    '''

    E.info("Clearing R environment variables")
    R('''closeAllConnections()''')
    R('''rm(list=ls(all=TRUE))''')
    R('''gc()''')


def mergeGateArrays(db, table_name, filter_zero=True, filter_gate=None):
    '''
    Merge arrays for each triboolean gate across batches,
    filter out zero sum arrays/gates.

    Arguments
    ---------
    db: string
      an absolute path to an SQLite DB containing the arrays

    table_name: string
      the name of the table to pull the arrays out of

    filter_zero: boolean
      if True, drop arrays whose element sum is not greater than 0

    Yields
    -------
    _df: Pandas.core.dataframe
      a pandas dataframe of all subjects across all parameters
      includes sample ID, processing batch ID and antibody panel
    '''

    E.info("connecting to SQLite database")
    # r_connect = setupRSqlConnection(db)
    py_connect = sql.connect(db)
    # get the gate, panel IDs - use pdsql
    query = '''SELECT gate, batch, panel FROM %s''' % table_name
    E.info("retrieving records from DB")
    records = pdsql.read_sql(query, py_connect)
    gates = set(zip(records["gate"], records["panel"]))

    if filter_gate:
        pass
    else:
        filter_gate = ""

    # get the R functions
    E.info("iterating over all %i triboolean gates" % len(gates))
    for gate, panel in gates:
        if re.search(filter_gate, gate):
            E.info("filtering out nonsensical gates %s" % gate)
            pass
        else:
            R('''source("/ifs/devel/projects/proj052/flow_pipeline/R_scripts/FlowProcess.R")''')
            R('''source("/ifs/devel/projects/proj052/flow_pipeline/R_scripts/arrayMunge.R")''')
            # shunt this into an R script and call on each
            # iteration
            # this was causing memory useage to inflate bcos
            # R is shit at memory management
            R('''sink(file="sink.txt")''')
            R('''rm(sub.array)''')
            R('''gc()''')
            E.info("extracting array from gate %s" % gate)
            state = '''SELECT * FROM %(table_name)s''' % locals()
            '''WHERE gate = "%(gate)s" AND panel = "%(panel)s"''' % locals()
            R('''sub.array <- unserializeArray(rdb="%(db)s", '''
              '''statement="%(state)s", gate="%(gate)s", '''
              '''panel="%(panel)s", table="%(table_name)s")''' % locals())
            E.info("unserializing array")
            E.info("appending arrays across processing batches")
            E.info("adding metadata back into dataframes")
            # shunt back into Python for downstream processing
            py_df = R["sub.array"]
            _df = pandas2ri.ri2py_dataframe(py_df)
            element_sum = sum(np.sum(_df.iloc[:,:22], axis=1))
            _df["panel"] = panel
            R('''sink(file=NULL)''')
            # both floats and int zeros evaluate as False in Python!
            if filter_zero and element_sum:
                E.info("retrieving data from gate %s and panel %s" % (gate,
                                                                      panel))
                # finish processing
                # call the Python garbage collector
                gc.collect()
                yield _df
            else:
                E.info("no data in arrays from gate %s and panel %s" % (gate,
                                                                        panel))
                pass


def mergeArrayWithDemographics(flow_arrays, id_column,
                               demo_file, demo_id_column):
    '''
    Merge an array from flow cytometry analysis with demographic
    data on samples

    Arguments
    ---------
    flow_array: generator
      a generators, elements of which are dataframes containing the
      relevant flow cytometry-derived data with additional meta-data
      columns

    id_columns: string
      column header from `flow_array` that identifies samples for merging

    demo_file: string
      absolute path to a file containing the demographics data

    demo_id_column: string
      column header in the `demo_file` that identifies samples

    Yields
    -------
    dataframe: pandas.Core.DataFrame
      a merged dataframe of flow cytometry data and sample demographic data
    '''

    E.info("reading demographics file %s" % demo_file)
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

    E.info("retrieving all arrays")
    for array in flow_arrays:
        # get twin IDs in flow data to match those in demographics data
        if type(id_column) == int:
            twin_id_col = array.columns[id_column]
            array["twin_num"] = [int(xi.split("_")[-1])
                                 for xi in array[twin_id_col].values]
        elif type(id_column) == str:
            array["twin_num"] = [int(xs.split("_")[-1])
                                 for xs in array[id_column].values]
        else:
            raise ValueError("Twin ID columns must be either column number or header")

        merge_df = pd.merge(left=array, right=demo_df, left_on="twin_num",
                            right_on=demo_id_col, how="inner")
        merge_indx = [ix for ix, iy in enumerate(merge_df.index)]
        merge_df.index = merge_indx
        
        # drop duplicates
        merge_df.drop_duplicates(subset=["twin_num", "gate", "age", "family_id",
                                         "batch", "panel"], inplace=True)
        yield merge_df


def split_zygosity(infile, zygosity_column, id_headers,
                   pair_header):
    '''
    Split the twins file by their zygosity and match up appropriate twin
    pairs based on their family ID
    '''

    # catch empty files and fail nicely
    try:
        twins_df = pd.read_table(infile, sep="\t", header=0, index_col=0)
    except ValueError:
        return pd.DataFrame

    MZ_df = twins_df[twins_df[zygosity_column] == "MZ"]
    DZ_df = twins_df[twins_df[zygosity_column] == "DZ"]

    # split each twin pair and match up on family ID
    # TODO: make report/output/log of twins not properly paired
    MZ_1 = MZ_df.loc[MZ_df.duplicated(subset=pair_header)]
    MZ_2 = MZ_df.loc[[not im for im in MZ_df.duplicated(subset=pair_header)]]

    DZ_1 = DZ_df.loc[DZ_df.duplicated(subset=pair_header)]
    DZ_2 = DZ_df.loc[[not di for di in DZ_df.duplicated(subset=pair_header)]]

    # remove cell size (FSC) and keep in all frames as an additional variable
    # melt the dataframes in order to merge twin pairs
    MZ1_melt = pd.melt(MZ_1, id_vars=id_headers, var_name="marker",
                       value_name="twin1")
    MZ2_melt = pd.melt(MZ_2, id_vars=id_headers, var_name="marker",
                       value_name="twin2")
    matched_MZ = pd.merge(left=MZ1_melt, right=MZ2_melt,
                          on=["marker", zygosity_column, "gate", "batch",
                              pair_header],
                          how='inner')

    DZ1_melt = pd.melt(DZ_1, id_vars=id_headers, var_name="marker",
                       value_name="twin1")
    DZ2_melt = pd.melt(DZ_2, id_vars=id_headers, var_name="marker",
                       value_name="twin2")

    matched_DZ = pd.merge(left=DZ1_melt, right=DZ2_melt,
                          on=["marker", zygosity_column, "gate", "batch",
                              pair_header],
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

        E.info("adjusting model for %s" % confounding)
        # run linear model, adjusting for confounding variables
        # hardcode adjustment for age and cell size
        twin1_r = r_lm("twin1 ~ age_x + %s " % confounding, data=r_group)
        twin2_r = r_lm("twin2 ~ age_y + %s " % confounding, data=r_group)

        # pull out residuals
        twin1_residuals = twin1_r.rx("residuals")[0]
        twin2_residuals = twin2_r.rx("residuals")[0]

        # pull out indicies and convert back to python objects,
        # make sure index values are integers
        tw1_indx = [int(ix1) for ix1 in twin1_residuals.names[0]]
        tw2_indx = [int(ix2) for ix2 in twin2_residuals.names[0]]

        tw1_rs = [fi1 for fi1 in twin1_residuals]
        tw2_rs = [fi2 for fi2 in twin2_residuals]

        # transform into Z-scores

        # tw1_z = [(fz1 - np.mean(tw1_rs))/np.std(tw1_rs) for fz1 in tw1_rs]
        # tw2_z = [(fz2 - np.mean(tw2_rs))/np.std(tw2_rs) for fz2 in tw2_rs]

        # tw1_ser = pd.Series(tw1_z, index=group_df.index, dtype=np.float64)
        # tw2_ser = pd.Series(tw2_z, index=group_df.index, dtype=np.float64)

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


def make_kinship_matrix(twins_file, id_column, 
                        family_column, zygosity_column):
    '''
    Using a file of twin data with a zygosity and family header,
    generate and expected genetic relationship matrix.
    '''

    E.info("reading twins file: %s" % twins_file)
    twin_df = pd.read_table(twins_file, sep="\t", index_col=None,
                            header=0)
    # some individuals are repeats, denoted by "rep" in the zygosity
    # column - remove these
    twin_df = twin_df[twin_df[zygosity_column] != "rep"]

    n_twins = len(twin_df)
    # setup a kinship coefficient matrix
    E.info("setting up kinship matrix")
    kinship = pd.DataFrame(columns=twin_df.index,
                           index=twin_df.index)

    # iterate over all twins, kinship = 1 if MZ and same
    # family, kinsip = 0.5 if DZ and same family, else 0
    twin_pairs = itertools.product(twin_df.index, twin_df.index)
    E.info("iterating over all twin pairs")
    for twin1, twin2 in twin_pairs:
        fam1 = twin_df.loc[twin1, family_column]
        fam2 = twin_df.loc[twin2, family_column]
        if twin1 == twin2:
            kin = 1.0
        elif fam1 == fam2:
            E.info("twin pair found %s:%s" % (twin1, twin2))
            if twin_df.loc[twin1, zygosity_column]  == "MZ":
                kin = 1.0
            elif twin_df.loc[twin1, zygosity_column] == "DZ":
                kin = 0.5
            else:
                E.info("These are not twins")
                kin = 0.0
        else:
            kin = 0.0

        kinship.loc[twin1, twin2] = kin

    return kinship
    

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
    # use the R package `ICC` function `ICCest`
    R('''suppressPackageStartupMessages(library(ICC))''')
    py_icc = R["ICCest"]
    pandas2ri.activate()

    MZ_corr_dict = {}
    DZ_corr_dict = {}
    for dname, dgroup in DZ_groups:
        # select family id and measurements
        dsub_df = dgroup[["family_id", "twin1_res", "twin2_res"]]

        # some expression measures are all 0 or numpy.nan, drop
        # these from analysis to prevent downstream errors
        # mask values where one twin is 0 as this may arise
        # due to missing values
        dsub_df["twin1_res"] = np.nan_to_num(dsub_df["twin1_res"])
        dsub_df["twin2_res"] = np.nan_to_num(dsub_df["twin2_res"])
        tw1_zero = dsub_df["twin1_res"] == 0
        tw2_zero = dsub_df["twin2_res"] == 0
        dsub_df = dsub_df[~tw1_zero]
        dsub_df = dsub_df[~tw2_zero]

        dmelt_sub = pd.melt(dsub_df, id_vars="family_id",
                            value_name="value", var_name="twin")
        if len(set(dmelt_sub["value"])) > 1:
            # shunt into R environment
            dmelt_r = pandas2ri.py2ri_pandasdataframe(dmelt_sub)
            # use default parameters for confidence intervals
            # alpha = 0.05
            d_icc_est = py_icc("family_id", "value", dmelt_r)
            d_icc = [fx for fx in d_icc_est.rx("ICC")[0]][0]
            DZ_corr_dict[dname] = d_icc
        else:
            DZ_corr_dict[dname] = 0.0

    for mname, mgroup in MZ_groups:
        # select family id and measurements
        msub_df = mgroup[["family_id", "twin1_res", "twin2_res"]]

        # some expression measures are all 0 or numpy.nan, drop
        # these from analysis to prevent downstream errors
        # mask values where one twin is 0 as this may arise
        # due to missing values
        msub_df["twin1_res"] = np.nan_to_num(msub_df["twin1_res"])
        msub_df["twin2_res"] = np.nan_to_num(msub_df["twin2_res"])
        mtw1_zero = msub_df["twin1_res"] == 0
        mtw2_zero = msub_df["twin2_res"] == 0
        msub_df = msub_df[~mtw1_zero]
        msub_df = msub_df[~mtw2_zero]

        mmelt_sub = pd.melt(msub_df, id_vars="family_id",
                            value_name="value", var_name="twin")
        if len(set(mmelt_sub["value"])) > 1:
            # shunt into R environment
            mmelt_r = pandas2ri.py2ri_pandasdataframe(mmelt_sub)
            # use default parameters for confidence intervals
            # alpha = 0.05
            m_icc_est = py_icc("family_id", "value", mmelt_r)
            m_icc = [mf for mf in m_icc_est.rx("ICC")[0]][0]
            MZ_corr_dict[mname] = m_icc
        else:
            MZ_corr_dict[mname] = 0.0

    # can only calculate the broad-sense heritability using Falconer's
    # equation as additive and dominance variance components cannot
    # be separated due to close familial relationship
    # these heritability estimates are the upper bounds of heritability
    broad_h = {}

    # remove markers with NA name
    for key in MZ_corr_dict.keys():
        if not re.search("NA", key):
            mz_cor = MZ_corr_dict[key]
            dz_cor = DZ_corr_dict[key]
            H2 = 2 * (mz_cor - dz_cor)
            broad_h[key] = H2
        else:
            pass
   
    return broad_h


def merge_heritability(file_list):
    '''
    Merge together heritability estimates
    into a single table?? How feasible is this??
    '''

    # get the cell type and marker panel from the file name
    file1 = file_list[0]
    fname = "_".join(file1.split("/")[-1].split("-")[1:-1])

    df = pd.read_table(file1, header=0, index_col=None, sep=":")
    # assign the heritability calcualtion the cell type
    df.columns = [fname]
    # markers are consistent across each panel
    # merge on index
    file_list.remove(file1)

    for files in file_list:
        finame = "_".join(files.split("/")[-1].split("-")[1:-1])
        _df = pd.read_table(files, header=0, index_col=None, sep=":")
        _df.columns = [finame]

        df = pd.merge(df, _df, left_index=True, right_index=True,
                      how='inner')

    return df


def _matrix_distance(matrix1, matrix2, distance="Euclid"):
    '''
    Calculate the distance between two matrices using
    the difference between matrix norms (Frobenius). Matrices
    must be conformable.

    Arguments
    ---------
    matrix1: numpy.ndarray
      a matrix of values, it must have the same dimension
      as matrix 2

    matrix2: numpy.ndarray
      a matrix of values, it must have the same dimension
      as matrix 1

    distance: string
      the distance measure to use.  Currently on the sum
      of squared differences is used.

    Returns
    -------
    dist: numpy.float64
      a distance metric value of the distances between matrix1
      and matrix 2
    '''

    if len(matrix1) == len(matrix2):
        fnorm1 = np.linalg.norm(matrix1)
        fnorm2 = np.linalg.norm(matrix2)
        dist = abs(fnorm1 - fnorm2)

        return dist

    else:
        raise AttributeError("matrix is not conformable."
                             "Matrices must have the same "
                             "dimension")


def getMatrixDistances(list_of_matrices, distance="Euclid"):
    '''
    Take in a list of matrices and calculate all
    pairwise distances between them.  Uses the distance
    based on matrix norms.  Matrices must be conformable,
    if they are not then the distance will be based on
    the overlap in samples names.  This is assumed to correspond
    to the rownames of each matrix.
    '''

    # dictionaries are unordered, lists are not
    mat_names = list_of_matrices.keys()
    n_pairs = len(mat_names)
    mat_idx = [x for x,y in enumerate(mat_names)]
    idx_dict = dict(zip(mat_names, mat_idx))

    distance_matrix = pd.DataFrame(columns=mat_names,
                                   index=mat_names)
    distance_matrix = distance_matrix.fillna(0.0)

    all_pairs = itertools.combinations_with_replacement(mat_names,
                                                        2)
    # iterate over all matrix pairs, make them conformable
    # then calculate their distances

    E.info("calculating all pair-wise distances "
           "between matrices of size %i by %i" % (n_pairs,
                                                  n_pairs))
    for pair in all_pairs:
        pair1 = list_of_matrices[pair[0]]
        pair2 = list_of_matrices[pair[1]]
        names1 = pair1.index
        names2 = pair2.index
        conform = set(names1).intersection(names2)

        mat1 = pair1.loc[conform].values
        mat2 = pair2.loc[conform].values

        idx1 = idx_dict[pair[0]]
        idx2 = idx_dict[pair[1]]

        distance = _matrix_distance(mat1, mat2)
        distance_matrix.iloc[idx1, idx2] = distance

    distance_df = pd.DataFrame(distance_matrix,
                               columns=mat_names,
                               index=mat_names)
    return distance_df


def makeMatrixList(list_of_files, id_col):
    '''
    Read in a list of files as dataframes, extract the
    important numerical values (1st 13 columns), set the
    rownames as unique IDs and return a dict of pandas
    dataframes
    '''

    E.info("reading matrices")
    count = 0
    matrix_list = {}
    for infile in list_of_files:
        # get the file name
        part1 = "-".join(infile.split("-")[:-1])
        fname = part1.split("/")[-1]

        _df = pd.read_table(infile, sep="\t", header=0,
                            index_col=None)
        _df.index = _df[id_col]
        df = _df.iloc[:, :13]
        matrix_list[fname] = df
        count += 1

    E.info("%i matrices read" % count)

    return matrix_list


def plotHistogram(data, variable, save_path,
                  x_title=None, y_title=None,
                  colour_var=None, scales=None,
                  split_var=None):
    '''
    plot a histogram
    '''
    pandas2ri.activate()
    R('''suppressPackageStartupMessages(library(ggplot2))''')
    
    # convert to R dataframe
    data_r = pandas2ri.py2ri_pandasdataframe(data)
    R.assign("data.df", data_r)

    if colour_var:
        R('''g_hist <- ggplot(data.df, aes(%(variable)s, fill=%(colour_var)s)) + '''
          '''geom_histogram() + theme_bw()''' % locals())
    else:
        R('''g_hist <- ggplot(data.df, aes(%(variable)s)) + '''
          '''geom_histogram()''' % locals())
    if split_var and scales:
        R('''g_hist <- g_hist + facet_wrap(~%(split_var)s, '''
          '''scale="%(scales)s")''' % locals())
    elif split_var and not scales:
        R('''g_hist <- g_hist + facet_wrap(~%(split_var)s)''' % locals())
    else:
        pass

    if x_title:
        R('''g_hist <- g_hist + xlab("%(x_title)s")''' % locals())
    else:
        pass

    if y_title:
        R('''g_hist <- g_hist + ylab("%(y_title)s")''' % locals())
    else:
        pass

    R('''png("%(save_path)s", height=720, width=720)''' % locals())
    R('''print(g_hist)''')
    R('''dev.off()''')


def plotBarchart(data, x_variable, y_variable,
                 save_path,
                 x_title=None, y_title=None,
                 colour_var=None,
                 split_var=None):
    '''
    Make a barchart
    '''

    pandas2ri.activate()
    R('''suppressPackageStartupMessages(library(ggplot2))''')
    
    # convert to R dataframe
    data_r = pandas2ri.py2ri_pandasdataframe(data)
    R.assign("data.df", data_r)

    if colour_var:
        R('''g_bar <- ggplot(data.df, '''
          '''aes(x=%(x_variable)s, y=%(y_variable)s), fill="#473C8B") + '''
          '''theme_bw() + geom_bar(stat="identity", fill="#473C8B") + '''
          '''theme(axis.text.x=element_text(angle=45, vjust=0.6, '''
          '''size=14, colour="black")) '''% locals())
    else:
        R('''g_bar <- ggplot(data.df, '''
          '''aes(x=%(x_variable)s, y=%(y_variable)s)) + '''
          '''theme_bw() + geom_bar(stat="identity", fill="#473C8B") + '''
          '''theme(axis.text.x=element_text(angle=45, vjust=0.6, '''
          '''size=14, colour="black")) '''% locals())

    if split_var:
        R('''g_bar <- g_bar + facet_wrap(~%(split_var)s)''' % locals())
    else:
        pass

    if x_title:
        R('''g_bar <- g_bar + xlab("%(x_title)s")''' % locals())
    else:
        pass

    if y_title:
        if re.search(y_title, "expression(paste"):
            R('''g_bar <- bar + ylab(%(y_title)s)''' % locals())
        else:
            R('''g_bar <- g_bar + ylab("%(y_title)s")''' % locals())
    else:
        pass

    R('''png("%(save_path)s", height=720, width=720)''' % locals())
    R('''print(g_bar)''')
    R('''dev.off()''')


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
    R('''suppressPackageStartupMessages(library(RSQLite))''')
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
                              cell_subset="NKT_early", db="csvdb"):
    '''
    Retrieve summary statistic across early NKT cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak
    R('''lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, '''
      '''num_peaks=3, ref_peak=2, filter_id="TCells")''')
    R('''add(wf, lg, parent="warping")''')

    # gate on the CD1d-multimeric complex positive cells
    R('''nkt <- mindensity(fr=Data(wf[["TCells+"]])[[1]], channel="G560-A", positive=T,'''
      '''filter_id="NKT")''')
    R('''add(wf, nkt, parent="TCells+")''')

    # some NKT subset may be very rare/missing in some individuals
    # check these, if they are missing then remove them
    # NKT subsets are determined on CD4 vs CCR5 and
    # CD4 vs CD8 gating
    R('''nk_sub <- quadGate(.gate=c("R780-A"=0.0, "V605-A"=0.0), filterId="CD4CCR5")''')
    R('''add(wf, nk_sub, parent="NKT+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CCR5-CD4+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CCR5-CD4+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CCR5-CD4+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CCR5-CD4+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CCR5-CD4+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_naive_nktcells_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                              cell_subset="NKT_naive", db="csvdb"):
    '''
    Retrieve summary statistic across early NKT cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak
    R('''lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, '''
      '''num_peaks=3, ref_peak=2, filter_id="TCells")''')
    R('''add(wf, lg, parent="warping")''')

    # gate on the CD1d-multimeric complex positive cells
    R('''nkt <- mindensity(fr=Data(wf[["TCells+"]])[[1]], channel="G560-A", positive=T,'''
      '''filter_id="NKT")''')
    R('''add(wf, nkt, parent="TCells+")''')

    # some NKT subset may be very rare/missing in some individuals
    # check these, if they are missing then remove them
    # NKT subsets are determined on CD4 vs CCR5 and
    # CD4 vs CD8 gating
    R('''nk_sub <- quadGate(.gate=c("R780-A"=0.0, "V605-A"=0.0), filterId="CD4CCR5")''')
    R('''add(wf, nk_sub, parent="NKT+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CCR5+CD4+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CCR5+CD4+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CCR5+CD4+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CCR5+CD4+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CCR5+CD4+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                      cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')

    
def get_terminal_nktcells_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                                 cell_subset="NKT_terminal", db="csvdb"):
    '''
    Retrieve summary statistic across early NKT cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak
    R('''lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, '''
      '''num_peaks=3, ref_peak=2, filter_id="TCells")''')
    R('''add(wf, lg, parent="warping")''')

    # gate on the CD1d-multimeric complex positive cells
    R('''nkt <- mindensity(fr=Data(wf[["TCells+"]])[[1]], channel="G560-A", '''
      '''positive=T, filter_id="NKT")''')
    R('''add(wf, nkt, parent="TCells+")''')

    # some NKT subset may be very rare/missing in some individuals
    # check these, if they are missing then remove them
    # NKT subsets are determined on CD4 vs CCR5 and
    # CD4 vs CD8 gating
    R('''nk_sub <- quadGate(.gate=c("V605-A"=0, "V585-A"=0), filterId="CD4CD8")''')
    R('''add(wf, nk_sub, parent="NKT+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_effector_nktcells_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                                 cell_subset="NKT_effector", db="csvdb"):
    '''
    Retrieve summary statistic across early NKT cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    E.info("Gating around lymphocyte population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak
    R('''lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, '''
      '''num_peaks=3, ref_peak=2, filter_id="TCells")''')
    R('''add(wf, lg, parent="warping")''')

    # gate on the CD1d-multimeric complex positive cells
    R('''nkt <- mindensity(fr=Data(wf[["TCells+"]])[[1]], channel="G560-A", positive=T,'''
      '''filter_id="NKT")''')
    R('''add(wf, nkt, parent="TCells+")''')

    # some NKT subset may be very rare/missing in some individuals
    # check these, if they are missing then remove them
    # NKT subsets are determined on CD4 vs CCR5 and
    # CD4 vs CD8 gating
    R('''nk_sub <- quadGate(.gate=c("V605-A"=0.0, "V585-A"=0.0), filterId="CD4CD8")''')
    R('''add(wf, nk_sub, parent="NKT+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4-CD8+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4-CD8+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4-CD8+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4-CD8+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_Vd1_tcells_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                          cell_subset="Vd1_Tcells", db="csvdb"):
    '''
    Retrieve summary statistic across early Vd1+ invariant T cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    # fill in the missing gating here...
    E.info("Gating around lymphocyte population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak
    R('''lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, '''
      '''num_peaks=3, ref_peak=2, filter_id="TCells")''')
    R('''add(wf, lg, parent="warping")''')

    E.info("gating around Vd1+ T cells")
    R('''vd1_mat <- matrix(c(0, 15), ncol=1)''')
    R('''colnames(vd1_mat) <- "B515-A"''')
    R('''vd_g <- rectangleGate(.gate=vd1_mat, filterId="Vd1")''')
    R('''add(wf, vd_g, parent="TCells+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["Vd1+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["Vd1+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["Vd1+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "Vd1+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "Vd1+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_Vd2_vg9dim_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                          cell_subset="Vd2p_Vg9dim", db="csvdb"):
    '''
    Retrieve summary statistic across Vd2+Vgdim T cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    # fill in the missing gating here...
    E.info("Gating around lymphocyte population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak
    R('''lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, '''
      '''num_peaks=3, ref_peak=2, filter_id="TCells")''')
    R('''add(wf, lg, parent="warping")''')

    E.info("selecting V1-Vg9 expressing T cells")
    R('''vg9_mat <- matrix(c(0, 15, 15, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(vg9_mat) <- c("R660-A", "B515-A")''')
    R('''vg9_g <- polygonGate(.gate=vg9_mat, filterId="Vg9")''')
    R('''add(wf, vg9_g, parent="TCells+")''')

    E.info("gating on Vg9 and Vd2; selecting Vd2+Vg9dim cells")
    R('''vd2_mat <- matrix(c(5, 10, 10, 10, 10, '''
      '''5, 5, 10, 10, 10, -6, -6), ncol=2)''')
    R('''colnames(vd2_mat) <- c("R660-A", "G610-A")''')
    R('''vd2_g <- polygonGate(.gate=vd2_mat, filterId="Vd2Vg9dim")''')
    R('''add(wf, vd2_g, parent="Vg9+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["Vd2Vg9dim-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["Vd2Vg9dim-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["Vd2Vg9dim-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "Vd2Vg9dim-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "Vd2Vg9dim-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_Vd2n_vg9p_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                         cell_subset="Vd2n_Vg9p", db="csvdb"):
    '''
    Retrieve summary statistic across Vd2+Vgdim T cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    # use the warping algorithm to align intensity peaks ( match features)
    # across individuals for each measured marker
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')

    # need to remove parameters without observations, i.e. non-expressed proteins

    E.info("feature matching across fluorescence parameters")
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters,...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId="warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')
    R('''gc()''')

    # fill in the missing gating here...
    E.info("Gating around lymphocyte population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak
    R('''lg <- tailgate(fr=Data(wf[["warping"]])[[1]], channel="V800-A", positive=T, '''
      '''num_peaks=3, ref_peak=2, filter_id="TCells")''')
    R('''add(wf, lg, parent="warping")''')

    E.info("selecting V1-Vg9 expressing T cells")
    R('''vg9_mat <- matrix(c(0, 15, 15, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(vg9_mat) <- c("R660-A", "B515-A")''')
    R('''vg9_g <- polygonGate(.gate=vg9_mat, filterId="Vg9")''')
    R('''add(wf, vg9_g, parent="TCells+")''')

    E.info("gating on Vg9 and Vd2; selecting Vd2-Vgdim cells")
    R('''vd2_mat <- matrix(c(5, 10, 10, 10, 10, '''
      '''5, 5, 10, 10, 10, -6, -6), ncol=2)''')
    R('''colnames(vd2_mat) <- c("R660-A", "G610-A")''')
    R('''vd2_g <- polygonGate(.gate=vd2_mat, filterId="Vd2Vg9dim")''')
    R('''add(wf, vd2_g, parent="Vg9+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["Vd2Vg9dim+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["Vd2Vg9dim+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["Vd2Vg9dim+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "Vd2Vg9dim+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "Vd2Vg9dim+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_hemat_SC_panel5(fcs_dir, out_dir, setid, comp_matrix, panel="P5",
                        cell_subset="hemat_SCs", db="csvdb"):
    '''
    Retrieve summary statistic across Vd2+Vgdim T cells from panel 5
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("Gating around haematopoeitic stem cell population")
    # need to use an openCyto tail gate here, the lymphGate doesn't like this marker
    # panel for some reason - do it after feature matching over parameters
    # this will take a tailgate of the distribution of CD3 using the 1st positive
    # peak as the reference.  Values are taken to the right of this peak

    # haematopoeitic stem cells are identified as small population
    # by CD3 and CD34
    # this population identification may not be the most rigorous
    # define a polygon gate for this population before the warping
    R('''hsc_mat <- matrix(c(8.4,9,11,11,10,9.1,8.4,'''
      '''6.1,8.0,9.0,5.5,5.0,6.1),ncol=2)''')
    R('''colnames(hsc_mat) <- c("V450-A", "V800-A")''')
    R('''hsc_g <- polygonGate(.gate=hsc_mat, filterId="HSC")''')
    R('''add(wf, hsc_g, parent="boundFilt+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["HSC+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["HSC+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["HSC+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "HSC+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "HSC+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_early_nkcells_panel4(fcs_dir, out_dir, setid, comp_matrix, panel="P4",
                             cell_subset="early_NK", db="csvdb"):
    '''
    Retrieve summary statistic across early NK cells (CD56+CD16-) from Panel 4
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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
    R('''add(wf, tf)''')

    E.info("feature matching samples over all parameters")

    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="asinh")''')

    E.info("gating on CD3- cells")
    R('''nonlg_mat <- matrix(c(0, 2500, 2500, 0, 0,0, 116, 116), ncol=2)''')
    R('''colnames(nonlg_mat) <- c("SSC-A", "R780-A")''')
    R('''nonlg_g <- polygonGate(.gate=nonlg_mat, filterId="nonTCells")''')
    R('''add(wf, nonlg_g, parent="warping")''')

    E.info("filtering out CD16 and CD56 double negative cells")
    R('''nonnk_mat <- matrix(c(116, 256, 256, 50, 50, 84, 104, 116, '''
      '''0, 0, 256, 256, 84, 84, 72, 0), ncol=2)''')
    R('''colnames(nonnk_mat) <- c("V585-A", "V450-A")''')
    R('''nk_g <- polygonGate(.gate=nonnk_mat, filterId="NKCells")''')
    R('''add(wf, nk_g, parent="nonTCells+")''')

    R('''cd16neg_mat <- matrix(c(76, 180, 180, 50, 50, 76, 128,'''
    '''180, 200, 200, 100, 100), ncol=2)''')
    R('''colnames(cd16neg_mat) <- c("V450-A", "V585-A")''')
    R('''cd16neg_g <- polygonGate(.gate=cd16neg_mat, filterId="earlyNK")''')
    R('''add(wf, cd16neg_g, parent="NKCells+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["earlyNK+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["earlyNK+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["earlyNK+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "earlyNK+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "earlyNK+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_terminal_nkcells_panel4(fcs_dir, out_dir, setid, comp_matrix, panel="P4",
                             cell_subset="terminal_NK", db="csvdb"):
    '''
    Retrieve summary statistic across terminal NK cells (CD56dim CD16+) from Panel 4
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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
    R('''add(wf, tf)''')

    E.info("feature matching samples over all parameters")

    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="asinh")''')

    E.info("gating on CD3- cells")
    R('''nonlg_mat <- matrix(c(0, 2500, 2500, 0, 0,0, 116, 116), ncol=2)''')
    R('''colnames(nonlg_mat) <- c("SSC-A", "R780-A")''')
    R('''nonlg_g <- polygonGate(.gate=nonlg_mat, filterId="nonTCells")''')
    R('''add(wf, nonlg_g, parent="warping")''')

    E.info("filtering out CD16 and CD56 double negative cells")
    R('''nonnk_mat <- matrix(c(116, 256, 256, 50, 50, 84, 104, 116, '''
      '''0, 0, 256, 256, 84, 84, 72, 0), ncol=2)''')
    R('''colnames(nonnk_mat) <- c("V585-A", "V450-A")''')
    R('''nk_g <- polygonGate(.gate=nonnk_mat, filterId="NKCells")''')
    R('''add(wf, nk_g, parent="nonTCells+")''')

    R('''cd16hi_mat <- matrix(c(100, 256, 256, 172, 100, '''
      '''50, 50, 128, 72, 50), ncol=2)''')
    R('''colnames(cd16hi_mat) <-c ("V450-A", "V585-A")''')
    R('''cd16hi_g <- polygonGate(.gate=cd16hi_mat, filterId="terminalNK")''')
    R('''add(wf, cd16hi_g, parent="NKCells+")''')


    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["terminalNK+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["terminalNK+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["terminalNK+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "terminalNK+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "terminalNK+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())

    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_mature_nkcells_panel4(fcs_dir, out_dir, setid, comp_matrix, panel="P4",
                              cell_subset="mature_NK", db="csvdb"):
    '''
    Retrieve summary statistic across mature NK cells (CD56+ CD16+) from Panel 4
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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
    R('''add(wf, tf)''')

    E.info("feature matching samples over all parameters")

    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="asinh")''')

    E.info("gating on CD3- cells")
    R('''nonlg_mat <- matrix(c(0, 2500, 2500, 0, 0,0, 116, 116), ncol=2)''')
    R('''colnames(nonlg_mat) <- c("SSC-A", "R780-A")''')
    R('''nonlg_g <- polygonGate(.gate=nonlg_mat, filterId="nonTCells")''')
    R('''add(wf, nonlg_g, parent="warping")''')

    E.info("filtering out CD16 and CD56 double negative cells")
    R('''nonnk_mat <- matrix(c(116, 256, 256, 50, 50, 84, 104, 116, '''
      '''0, 0, 256, 256, 84, 84, 72, 0), ncol=2)''')
    R('''colnames(nonnk_mat) <- c("V585-A", "V450-A")''')
    R('''nk_g <- polygonGate(.gate=nonnk_mat, filterId="NKCells")''')
    R('''add(wf, nk_g, parent="nonTCells+")''')

    R('''cd56pos_mat <- matrix(c(100, 156, 228, 256, 256, 180, 80, 100, '''
      '''54, 76, 128, 150, 172, 172, 128, 50), ncol=2)''')
    R('''colnames(cd56pos_mat) <- c("V450-A", "V585-A")''')
    R('''cd56_g <- polygonGate(.gate=cd56pos_mat, filterId="matureNK")''')
    R('''add(wf, cd56_g, parent="NKCells+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["matureNK+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["matureNK+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["matureNK+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "matureNK+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "matureNK+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_immature_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                               cell_subset="immature_Bcell", db="csvdb"):
    '''
    Retrieve summary statistic across immature B cells (CD19CD20+CD21+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting immature B cell subset on CD10+ expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["immature+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["immature+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["immature+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "immature+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "immature+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_mature_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                             cell_subset="mature_Bcell", db="csvdb"):
    '''
    Retrieve summary statistic across mature B cells (CD19CD20+CD21-)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting mature B cell subset on CD10- expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')



    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["immature-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["immature-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["immature-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "immature-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "immature-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_memory_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                             cell_subset="memory_Bcell", db="csvdb"):
    '''
    Retrieve summary statistic across memory B cells (CD19CD20+CD21-CD95+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting mature B cell subset on CD10- expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')

    E.info("gating on memory B cells with CD95+ expression")
    R('''mem_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mem_mat) <- "G660-A"''')
    R('''mem_g <- rectangleGate(.gate=mem_mat, filterId="memory")''')
    R('''add(wf, mem_g, parent="immature-")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["memory+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["memory+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["memory+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "memory+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "memory+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_naive_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                            cell_subset="naive_Bcell", db="csvdb"):
    '''
    Retrieve summary statistic across naive B cells (CD19CD20+CD21-CD95-)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting mature B cell subset on CD10- expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')

    E.info("gating on naive  B cells with CD95- expression")
    R('''mem_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mem_mat) <- "G660-A"''')
    R('''mem_g <- rectangleGate(.gate=mem_mat, filterId="memory")''')
    R('''add(wf, mem_g, parent="immature-")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["memory-"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["memory-"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["memory-"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "memory-", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "memory-", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_IgA_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                          cell_subset="IgA_Bmem", db="csvdb"):
    '''
    Retrieve summary statistic across memory IgA+ B cells (CD19CD20+CD21-CD95+IgA+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting mature B cell subset on CD10- expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')

    E.info("gating on memory B cells with CD95+ expression")
    R('''mem_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mem_mat) <- "G660-A"''')
    R('''mem_g <- rectangleGate(.gate=mem_mat, filterId="memory")''')
    R('''add(wf, mem_g, parent="immature-")''')

    E.info("Gating around IgA+ memory B cells")
    R('''iga_mat <- matrix(c(0, 10, 10, 0 -6, -6, 0 , 0), ncol=2)''')
    R('''colnames(iga_mat) <- c("G560-A", "V450-A")''')
    R('''iga_g <- polygonGate(.gate=iga_mat, filterId="IgA")''')
    R('''add(wf, iga_g, parent="memory+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["IgA+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["IgA+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["IgA+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "IgA+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "IgA+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_IgG_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                          cell_subset="IgG_Bmem", db="csvdb"):
    '''
    Retrieve summary statistic across memory IgG+ B cells (CD19CD20+CD21-CD95+IgG+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting mature B cell subset on CD10- expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')

    E.info("gating on memory B cells with CD95+ expression")
    R('''mem_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mem_mat) <- "G660-A"''')
    R('''mem_g <- rectangleGate(.gate=mem_mat, filterId="memory")''')
    R('''add(wf, mem_g, parent="immature-")''')

    E.info("Gating around IgG+ memory B cells")
    R('''igg_mat <- matrix(c(-6, 0, 0, -6, 0, 0, 10, 10), ncol=2)''')
    R('''colnames(igg_mat) <- c("G560-A", "V450-A")''')
    R('''igg_g <- polygonGate(.gate=igg_mat, filterId="IgG")''')
    R('''add(wf, igg_g, parent="memory+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["IgG+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["IgG+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["IgG+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "IgG+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "IgG+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_IgM_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                          cell_subset="IgM_Bmem", db="csvdb"):
    '''
    Retrieve summary statistic across memory IgM+ B cells (CD19CD20+CD21-CD95+IgG-IgA-IgM+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting mature B cell subset on CD10- expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')

    E.info("gating on memory B cells with CD95+ expression")
    R('''mem_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mem_mat) <- "G660-A"''')
    R('''mem_g <- rectangleGate(.gate=mem_mat, filterId="memory")''')
    R('''add(wf, mem_g, parent="immature-")''')

    E.info("Gating around IgM+ memory B cells")
    R('''noig_mat <- matrix(c(-6, 0, 0, -6, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(noig_mat) <- c("G560-A", "V450-A")''')
    R('''noig_g <- polygonGate(.gate=noig_mat, filterId="nonIg")''')
    R('''add(wf, noig_g, parent="memory+")''')
    R('''igm_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(igm_mat) <- "V585-A"''')
    R('''igm_g <- rectangleGate(.gate=igm_mat, filterId="IgM")''')
    R('''add(wf, igm_g, parent="nonIg+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["IgM+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["IgM+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["IgM+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "IgM+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "IgM+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_IgE_Bcells_panel6(fcs_dir, out_dir, setid, comp_matrix, panel="P6",
                          cell_subset="IgE_Bmem", db="csvdb"):
    '''
    Retrieve summary statistic across memory IgE+ B cells (CD19CD20+CD21-CD95+IgG-IgA-IgM-)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on CD19CD20+ cells")
    R('''cd19_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''cd20_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd19_mat) <- "V655-A"''')
    R('''colnames(cd20_mat) <- "V605-A"''')
    R('''cd19_g <- rectangleGate(.gate=cd19_mat, filterId="CD19")''')
    R('''cd20_g <- rectangleGate(.gate=cd20_mat, filterId="CD20")''')
    R('''bcell_g <- cd19_g | cd20_g''')
    R('''add(wf, bcell_g, parent="warping")''')

    # select on mature/immature subset CD10+/-
    E.info("selecting mature B cell subset on CD10- expression")
    R('''imm_mat <- matrix(c(10, 10, 15, 15, 0, 11.6, 16, 0), ncol=2)''')
    R('''colnames(imm_mat) <- c("R660-A", "G710-A")''')
    R('''imm_g <- polygonGate(.gate=imm_mat, filterId="immature")''')
    R('''add(wf, imm_g, parent="CD19 or CD20+")''')

    E.info("gating on memory B cells with CD95+ expression")
    R('''mem_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mem_mat) <- "G660-A"''')
    R('''mem_g <- rectangleGate(.gate=mem_mat, filterId="memory")''')
    R('''add(wf, mem_g, parent="immature-")''')

    E.info("Gating around IgE+ memory B cells")
    R('''noig_mat <- matrix(c(-6, 0, 0, -6, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(noig_mat) <- c("G560-A", "V450-A")''')
    R('''noig_g <- polygonGate(.gate=noig_mat, filterId="nonIg")''')
    R('''add(wf, noig_g, parent="memory+")''')
    R('''ige_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(ige_mat) <- "V585-A"''')
    R('''ige_g <- rectangleGate(.gate=ige_mat, filterId="IgE")''')
    R('''add(wf, ige_g, parent="nonIg+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["IgE+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["IgE+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["IgE+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "IgE+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "IgE+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_monocytes_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                         cell_subset="monocytes", db="csvdb"):
    '''
    Retrieve summary statistic across monocytes (Lin-HLA-DR+CD14+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14+ monocytes")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD14+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD14+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD14+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD14+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD14+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_cd123cd11c_DC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                             cell_subset="DC_cd123cd11c", db="csvdb"):
    '''
    Retrieve summary statistic across CD123+CD11c double-positive DCs
    (Lin-HLA-DR+CD14-CD11c+CD123+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD123 and CD11c+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD123 and CD11c+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD123 and CD11c+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD123 and CD11c+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD123 and CD11c+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_myeloid_DC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                          cell_subset="M_dendritic", db="csvdb"):
    '''
    Retrieve summary statistic across myeloid DCs  (Lin-HLA-DR+CD14-CD11cCD123+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    E.info("gating around myeloid dendritic cells")
    R('''mdc_mat <- matrix(c(2, 15, 15, 2, -6, -6, 15, 15), ncol=2)''')
    R('''colnames(mdc_mat) <- c("G660-A", "G560-A")''')
    R('''mdc_g <- polygonGate(.gate=mdc_mat, filterId="MDC")''')
    R('''add(wf, mdc_g, parent="CD123 and CD11c-")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["MDC+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["MDC+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["MDC+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "MDC+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "MDC+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_plasmacytoid_DC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                               cell_subset="P_dendritic", db="csvdb"):
    '''
    Retrieve summary statistic across plasmacytoid DCs (Lin-HLA-DR+CD14-CD11c-CD123+)
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    E.info("gating around plasmacytoid dendritic cells")
    R('''pdc_mat <- matrix(c(-6, 0, 0, 15), ncol=2)''')
    R('''colnames(pdc_mat) <- c("G660-A", "G560-A")''')
    R('''pdc_g <- rectangleGate(.gate=pdc_mat, filterId="PDC")''')
    R('''add(wf, pdc_g, parent="CD123 and CD11c-")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["PDC+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["PDC+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["PDC+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "PDC+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "PDC+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_CD8_APC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                       cell_subset="CD8_APCs", db="csvdb"):
    '''
    Retrieve summary statistic across APCs for CD8+ T cells
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    E.info("gating around myeloid dendritic cells")
    R('''mdc_mat <- matrix(c(2, 15, 15, 2, -6, -6, 15, 15), ncol=2)''')
    R('''colnames(mdc_mat) <- c("G660-A", "G560-A")''')
    R('''mdc_g <- polygonGate(.gate=mdc_mat, filterId="MDC")''')
    R('''add(wf, mdc_g, parent="CD123 and CD11c-")''')

    E.info("selecting APCs for CD8 T cells")
    R('''cd8_mat <- matrix(c(-5, 15, 0, 15), ncol=2)''')
    R('''colnames(cd8_mat) <- c("R780-A", "B515-A")''')
    R('''cd8apc_g <- rectangleGate(.gate=cd8_mat, filterId="CD8APC")''')
    R('''add(wf, cd8apc_g, parent="MDC+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD8APC+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD8APC+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD8APC+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD8APC+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD8APC+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_CD4_APC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                       cell_subset="CD4_APCs", db="csvdb"):
    '''
    Retrieve summary statistic across APCs for CD4+ T cells
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    E.info("gating around myeloid dendritic cells")
    R('''mdc_mat <- matrix(c(2, 15, 15, 2, -6, -6, 15, 15), ncol=2)''')
    R('''colnames(mdc_mat) <- c("G660-A", "G560-A")''')
    R('''mdc_g <- polygonGate(.gate=mdc_mat, filterId="MDC")''')
    R('''add(wf, mdc_g, parent="CD123 and CD11c-")''')

    E.info("selecting APCs for CD4 T cells")
    R('''cd4_mat <- matrix(c(0, 15, -6, 0), ncol=2)''')
    R('''colnames(cd4_mat) <- c("R780-A", "B515-A")''')
    R('''cd4apc_g <- rectangleGate(.gate=cd4_mat, filterId="CD4APC")''')
    R('''add(wf, cd4apc_g, parent="MDC+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD4APC+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD4APC+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD4APC+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD4APC+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD4APC+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_CD1cneg_mDC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                           cell_subset="CD1cneg_mDC", db="csvdb"):
    '''
    Retrieve summary statistic across CD1c- myeloid DCs
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    E.info("gating around myeloid dendritic cells")
    R('''mdc_mat <- matrix(c(2, 15, 15, 2, -6, -6, 15, 15), ncol=2)''')
    R('''colnames(mdc_mat) <- c("G660-A", "G560-A")''')
    R('''mdc_g <- polygonGate(.gate=mdc_mat, filterId="MDC")''')
    R('''add(wf, mdc_g, parent="CD123 and CD11c-")''')

    E.info("selecting CD1c- DCs")
    R('''cd1cneg_mat <- matrix(c(-6, 0, -6, 0), ncol=2)''')
    R('''colnames(cd1cneg_mat) <- c("R780-A", "B515-A")''')
    R('''cd1cneg_g <- rectangleGate(.gate=cd1cneg_mat, filterId="CD1cneg")''')
    R('''add(wf, cd1cneg_g, parent="MDC+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD1cneg+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD1cneg+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD1cneg+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD1cneg+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD1cneg+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_inflammatory_DC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                               cell_subset="Inflam_DC", db="csvdb"):
    '''
    Retrieve summary statistic across inflammatory (CD16+) DCs
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    E.info("gating around myeloid dendritic cells")
    R('''mdc_mat <- matrix(c(2, 15, 15, 2, -6, -6, 15, 15), ncol=2)''')
    R('''colnames(mdc_mat) <- c("G660-A", "G560-A")''')
    R('''mdc_g <- polygonGate(.gate=mdc_mat, filterId="MDC")''')
    R('''add(wf, mdc_g, parent="CD123 and CD11c-")''')

    E.info("selecting CD1c- DCs")
    R('''cd1cneg_mat <- matrix(c(-6, 0, -6, 0), ncol=2)''')
    R('''colnames(cd1cneg_mat) <- c("R780-A", "B515-A")''')
    R('''cd1cneg_g <- rectangleGate(.gate=cd1cneg_mat, filterId="CD1cneg")''')
    R('''add(wf, cd1cneg_g, parent="MDC+")''')

    E.info("gating on CD16+ inflammatory DCs")
    R('''cd16_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd16_mat) <- "G710-A"''')
    R('''cd16_g <- rectangleGate(.gate=cd16_mat, filterId="CD16")''')
    R('''add(wf, cd16_g, parent="CD1cneg+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD16+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD16+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD16+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD16+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD16+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')


def get_CD16neg_DC_panel7(fcs_dir, out_dir, setid, comp_matrix, panel="P7",
                          cell_subset="CD16neg_DC", db="csvdb"):
    '''
    Retrieve summary statistic across CD16- DCs
    '''

    R('''sink(file="sink.txt")''')
    R('''suppressPackageStartupMessages(library(flowCore))''')
    R('''suppressPackageStartupMessages(library(flowWorkspace))''')
    R('''suppressPackageStartupMessages(library(flowStats))''')
    R('''suppressPackageStartupMessages(library(openCyto))''')   
    R('''suppressPackageStartupMessages(library(RSQLite))''')
    R('''source("/ifs/projects/proj052/src/flow_pipeline/R_scripts/FlowProcess.R")''')
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
    R('''flow.data <- fix_scatter_name(flow.data)''')

    # apply compensation matrix to all events
    E.info("apply compensations to defined fcs dimensions")
    R('''comp <- read.table(file="%s", sep="\t", h=T, row.names=1)''' % comp_matrix)
    R('''comp_splt <- strsplit(colnames(comp), split=".", fixed=T)''')
    R('''colnames(comp) <- unlist(lapply(comp_splt, FUN=function(x) {paste(x[1],'''
      '''x[2], sep="-")}))''')
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

    E.info("feature matching samples over all parameters")
    R('''mark_vec <- filter_markers((pData(parameters(flow.data[[1]]))[,"desc"])[c(4:21)])''')
    R('''pars <- colnames(Data(wf[["base view"]]))[mark_vec]''')
    R('''norm <- normalization(normFunction=function(x, parameters, ...) warpSet(x, parameters, ...),'''
      '''parameters=pars, normalizationId = "warping")''')
    R('''add(wf, norm, parent="boundFilt+")''')

    E.info("gating on lineage markers")
    R('''nonlin_mat <- matrix(c(0, 16, 16, 0, -6, -6, 0, 0), ncol=2)''')
    R('''colnames(nonlin_mat) <- c("R710-A", "V605-A")''')
    R('''nonlin_g <- polygonGate(.gate=nonlin_mat, filterId="Nonlineage")''')
    R('''add(wf, nonlin_g, parent="warping")''')

    E.info("selecting CD14- dendritic cells")
    R('''mono_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(mono_mat) <- "V800-A"''')
    R('''mono_g <- rectangleGate(.gate=mono_mat, filterId="CD14")''')
    R('''add(wf, mono_g, parent="Nonlineage+")''')
    
    # select out CD123+CD11c DCs
    R('''cd123_mat <- matrix(c(6.3258208265593, 16), ncol=1)''')
    R('''colnames(cd123_mat) <- "V800-A"''')
    R('''cd123_g <- rectangleGate(.gate=cd123_mat, filterId="CD123")''')
    R('''cd11c_mat <- matrix(c(7.93315738174014, 16), ncol=1)''')
    R('''colnames(cd11c_mat) <- "G660-A"''')
    R('''cd11c_g <- rectangleGate(.gate=cd11c_mat, filterId="CD11c")''')
    R('''dp_dc <- cd123_g & cd11c_g''')
    R('''add(wf, dp_dc, parent="CD14-")''')

    E.info("gating around myeloid dendritic cells")
    R('''mdc_mat <- matrix(c(2, 15, 15, 2, -6, -6, 15, 15), ncol=2)''')
    R('''colnames(mdc_mat) <- c("G660-A", "G560-A")''')
    R('''mdc_g <- polygonGate(.gate=mdc_mat, filterId="MDC")''')
    R('''add(wf, mdc_g, parent="CD123 and CD11c-")''')

    E.info("selecting CD1c- DCs")
    R('''cd1cneg_mat <- matrix(c(-6, 0, -6, 0), ncol=2)''')
    R('''colnames(cd1cneg_mat) <- c("R780-A", "B515-A")''')
    R('''cd1cneg_g <- rectangleGate(.gate=cd1cneg_mat, filterId="CD1cneg")''')
    R('''add(wf, cd1cneg_g, parent="MDC+")''')

    E.info("gating on CD16+ inflammatory DCs")
    R('''cd16_mat <- matrix(c(0, Inf), ncol=1)''')
    R('''colnames(cd16_mat) <- "G710-A"''')
    R('''cd16_g <- rectangleGate(.gate=cd16_mat, filterId="CD16")''')
    R('''add(wf, cd16_g, parent="CD1cneg+")''')

    # get all triboolean gates
    R('''params <- pData(parameters(Data(wf[["CD16+"]])[[1]]))[,"name"]''')
    R('''param_names <- pData(parameters(Data(wf[["CD16+"]])[[1]]))[,"desc"]''')
    R('''tribools <- make_booleans(Data(wf[["CD16+"]])[[1]], params, param_names)''')

    # apply the gates to the data
    # how much of a problem are missing and null values?
    # reset theset to 0?
    E.info("calculating summary statistics over all triboolean gates")
    R('''list_of_fanos <- get_frames(wf, "CD16+", get_fano, tribools)''')
    R('''list_of_means <- get_frames(wf, "CD16+", get_means, tribools)''')

    # use the RSQLite interface to generate BLOBS for arrays and insert
    # into RDMS - props to 
    # http://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    # on how to achieve this
    # the unserialize function no longer accepts a character string,
    # will need to get this back out in Python

    R('''con <- dbConnect(SQLite(), "%s")''' % db)
    # create the table if it doesn't already exist

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_fano '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      ''' panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    R('''dbGetQuery(con, "CREATE table if not exists %(cell_subset)s_mean '''
      '''(gate TEXT, batch TEXT, rows TEXT, columns TEXT, array TEXT, '''
      '''panel TEXT, PRIMARY KEY (gate, batch, panel))") ''' % locals())

    # set up column names, row names, gate and batch IDs
    # use '/' as the separator
    R('''markers <- paste(colnames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''twins <- paste(rownames(list_of_means[[1]]), collapse="/", sep="")''')
    R('''batch <- "%s"''' % setid)
    R('''panel <- "%s"''' % panel)

    # serialize the arrays into a dataframe
    # need to convert them to characters for later re-reading
    # to read data back use unserialize(charToRaw)
    R('''means_df <- data.frame(indx=1:length(list_of_means),'''
      '''arrays=I(unlist(lapply(list_of_means, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''fanos_df <- data.frame(indx=1:length(list_of_fanos),'''
      '''arrays=I(unlist(lapply(list_of_fanos, function(x) {rawToChar('''
      '''serialize(x, NULL, ascii=T))}))))''')

    R('''means_df$gate <- rownames(means_df)''')
    R('''fanos_df$gate <- rownames(fanos_df)''')

    R('''means_df$batch <- batch''')
    R('''fanos_df$batch <- batch''')
    R('''means_df$rows <- twins''')
    R('''fanos_df$rows <- twins''')
    R('''means_df$columns <- markers''')
    R('''fanos_df$columns <- markers''')
    R('''means_df$panel <- panel''')
    R('''fanos_df$panel <- panel''')

    # insert values into RDMS
    # these will fail if the combination of gate and batch are not unique
    # is there a way of checking whether a record exists first?
    # SQLite DB locking is a problem, use tryCatch in R to allow N tries??

    R('''fano_query <- "insert OR ignore into %(cell_subset)s_fano '''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    R('''mean_query <- "insert OR ignore into %(cell_subset)s_mean'''
      '''(array, gate, batch, rows, columns, panel) values '''
      '''(:arrays, :gate, :batch, :rows, :columns, :panel)"''' % locals())
    try:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=fano_query, dataframe=fanos_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_fano '''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=fanos_df)''' % locals()) 
    try:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
    except rinterface.RRuntimeError:
        R('''retry_DBwrite(dbh=con, query=mean_query, dataframe=means_df, nTries=3)''')
        # R('''dbGetPreparedQuery(con, "insert OR ignore into %(cell_subset)s_mean'''
        #   '''(array, gate, batch, rows, columns, panel) values '''
        #   '''(:arrays, :gate, :batch, :rows, :columns, :panel)", bind.data=means_df)''' % locals())
        
    outfile = "%s/%s-%s-%s.tsv" % (out_dir, setid, panel,
                                   cell_subset)
    os.system("touch %s" % outfile)
    R('''dbDisconnect(con)''')
    R('''sink(file=NULL)''')
