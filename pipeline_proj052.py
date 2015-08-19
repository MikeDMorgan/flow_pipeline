##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline template
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_proj052.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import re
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))

FLOW_DIR = PARAMS['flow_zip_dir']

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


# ---------------------------------------------------
# parse workspace files first to pull out compensation matrices
@follows(mkdir("comp_matrices.dir"))
@transform("%s/*.wsp"  %PARAMS['flow_workspace_dir'],
           regex("%s/(.+)_(.+).wsp" % PARAMS['flow_workspace_dir']),
           r"comp_matrices.dir/\2-compensation_matrix.txt")
def getCompensationMatrices(infile, outfile):
    '''
    Parse the workspace XML files to retrieve the spectral
    compensation values - put into matrix format
    '''

    job_memory = "2G"

    statement = '''
    python /ifs/devel/projects/proj052/scripts/workspace2tsv.py
    --method=compensation
    --log=%(outfile)s.log
    %(infile)s
    > %(outfile)s
    '''

    P.run()


# Unpack FlowRepository files
@follows(getCompensationMatrices,
         mkdir("fcs.dir"))
@subdivide("%s/Discovery/*.zip" % FLOW_DIR,
           regex("%s/FlowRepository_FR-FCM-(.+)_files.zip" % FLOW_DIR),
           r"fcs.dir/\1.dir")
def unzipFCS(infile, outfiles):
    '''
    Unpack zipped .fcs files into relevant directories
    '''

    # get 4 character file ID
    # unzip into the fcs.dir
    # process in R script and output files
    panel = PARAMS['flow_panel']
    fcs_dir = os.getcwd() + "/" + "fcs.dir"
    fileid = infile.split("/")[-1].split("-")[2].strip("_files.zip")
    command = []

    # make the archive specific directoy
    command.append(''' mkdir %s/%s.dir ''' % (fcs_dir, fileid))

    command.append(''' unzip -qj %s -d %s/%s.dir ''' % (infile,
                                                    fcs_dir,
                                                    fileid))

    statement = " ; ".join(command)

    P.run()


# parse the gating strategy file and generate a set of dummy files
# to dictate further .fcs file processing

@follows(getCompensationMatrices,
         unzipFCS,
         mkdir("gating.dir"))
@subdivide("%s" % PARAMS['flow_gating_file'],
           regex("(.+)/(.+)_design.tsv"),
           r"gating.dir/\2.txt")
def splitGatingFile(infile, outfile):
    '''
    Split the gating design file into separate
    files that contain panel number and cell type info
    '''

    statement = '''
    python /ifs/devel/projects/proj052/scripts/workspace2tsv.py
    --method=parse_gating
    --log=%(outfile)s.log
    --gating-directory=gating.dir
    %(infile)s
    > %(outfile)s
    '''

    P.run()


# create dummy files to match up panel, cell type and .fcs file group
@follows(getCompensationMatrices,
         splitGatingFile,
         mkdir("dummy.dir"))
@subdivide("gating.dir/*-gating.txt",
           regex("gating.dir/(.+)-(.+)-gating.txt"),
           add_inputs(r"fcs.dir/*.dir"),
           r"dummy.dir/\1-\2-*.dummy")
def matchFiles(infiles, outfile):
    '''
    Create dummy files - one per combination of .fcs file group, compensation
    matrix and gating strategy
    '''

    cell_file_ele = (infiles[0].split("/"))[-1].split("-")
    fcs_dirs = infiles[1:]

    for fdir in fcs_dirs:
        panel = cell_file_ele[0]
        cell_type = cell_file_ele[1]
        fcdir = fdir.split("/")[1]
        dummy_file = "-".join([panel, cell_type, fcdir])

        P.touch("dummy.dir/" + dummy_file + ".dummy")

# Process for each cell subset and panel combination
@job_limit(10)
@follows(getCompensationMatrices,
         splitGatingFile,
         matchFiles,
         mkdir("flow_tables.dir"))
@transform("dummy.dir/*.dummy",
           regex("dummy.dir/(.+)-(.+)-(.+).dummy"),
           add_inputs(r"comp_matrices.dir/\1-compensation_matrix.txt"),
           r"flow_tables.dir/\1-\2-fano.tsv")
def processFCS(infiles, outfiles):
    '''
    Sequentially unzip each file selecting for a specific
    marker panel if defined in the .ini file.
    '''

    # get 4 character file ID
    # unzip into the fcs.dir
    # process in R script and output files
    dummy = (infiles[0].strip("dummy.dir/")).split("-")
    comp_matrix = infiles[1]

    panel = dummy[0]
    cell_type = dummy[1]
    fcs_dir = os.getcwd() + "/" + "fcs.dir"
    fileid = dummy[-1].strip(".dummy")
    flow_panel = panel + ".dir"

    statement = '''
    python /ifs/devel/projects/proj052/scripts/fcs2tsv.py
    --fcs-directory=%s/%s.dir
    --summary-stats=fano
    --compensation-matrix=%s
    --cell-type=%s
    --panel-id=%s
    --output-directory=flow_tables.dir
    --fileset-identifier=%s
    --log=%s.log
    ''' % (fcs_dir, fileid, comp_matrix, cell_type,
           panel, fileid, "dummy.dir" + "/" + fileid)

    job_memory = "40G"

    P.run()


@follows(processFCS,
         mkdir("twin_files.dir"))
@collate(processFCS,
         regex("flow_tables.dir/(.+)-(.+)-(.+).tsv"),
         r"twin_files.dir/\2-\3.tsv")
def mergeFlowTables(infiles, outfile):
    '''
    Collate and merge all separate tables into a single
    large table for all MZ and DZ twins
    '''

    input_files = ",".join(infiles)
    demo_file = PARAMS['twins_demographics']
    demo_ids = PARAMS['twins_demo_header']
    twin_id = "twin.id"

    statement = '''
    python /ifs/devel/projects/proj052/scripts/flow2twins.py
    --task=merge_flow
    --twin-id-column=%(twin_id)s
    --demographics-file=%(twins_demographics)s
    --demo-id-column=%(twins_demo_header)s
    --log=%(outfile)s.log
    %(input_files)s
    > %(outfile)s'''

    P.run()


@follows(mergeFlowTables)
@subdivide(mergeFlowTables,
           regex("twin_files.dir/(.+)-(.+).tsv"),
           r"twin_files.dir/\1-\2-MZ.tsv")
def splitByZygosity(infile, outfile):
    '''
    Split the expression noise file into MZ and DZ
    separate files to match up twin pairs
    and process separately
    '''

    zygosity_column = "zygosity"
    id_headers = ",".join(["twin.id", "twin_num", "family_id",
                           "flowjo_id", "age", "subset", "zygosity",
                           "replicate", "visit"])
    pair_header = "family_id"
    output_pattern = infile.rstrip(".tsv")

    statement = '''
    python /ifs/devel/projects/proj052/scripts/flow2twins.py
    --task=split_zygosity
    --zygosity-column=%(zygosity_column)s
    --id-columns=%(id_headers)s
    --family-id-column=%(pair_header)s
    --output-file-pattern=%(output_pattern)s
    %(infile)s
    '''

    P.run()


@follows(splitByZygosity)
@transform("twin_files.dir/*.tsv",
           regex("twin_files.dir/(.+)-(.+)-(MZ|DZ).tsv"),
           r"twin_files.dir/\3-\1-\2-adjusted.tsv")
def regressConfounding(infile, outfile):
    '''
    Regress out confounding variables within zygosity classes
    for each twin.
    '''

    marker_column = "marker"
    confounding = "subset"
    job_memory = "8G"

    statement = '''
    python /ifs/devel/projects/proj052/scripts/flow2twins.py
    --task=regress_confounding
    --confounding-column=%(confounding)s
    --marker-group=%(marker_column)s
    --log=%(outfile)s.log
    %(infile)s
    > %(outfile)s
    '''

    P.run()


@follows(mkdir("heritability.dir"),
         regressConfounding)
@collate(regressConfounding,
         regex("twin_files.dir/(.+)-(.+)-(.+)-adjusted.tsv"),
         r"heritability.dir/\2-\3-heritability.tsv")
def calculateHeritability(infiles, outfile):
    '''
    Use Falconer's equation to calculate heritability
    estimates for all markers
    '''

    mono_file = [mx for mx in infiles if re.search("MZ", mx)][0]
    di_file = [dx for dx in infiles if re.search("DZ", dx)][0]

    statement = '''
    python /ifs/devel/projects/proj052/scripts/heritability.py
    --monozygote-file=%(mono_file)s
    --dizygote-file=%(di_file)s
    --log=%(outfile)s.log
    > %(outfile)s
    '''

    P.run()

# ---------------------------------------------------
# Generic pipeline tasks
@follows(processFCS)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
