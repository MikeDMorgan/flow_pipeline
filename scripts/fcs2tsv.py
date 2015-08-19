'''
fcs2tsv.py - convert fcs file to text files containing fluorescence intensities
===============================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Apply gating strategy to fcs files in :term:``input_directory``, pull out
either intensity information of specfic summary statistics

Usage
-----

.. Example use case

Example::

   python fcs2tsv.py

Type::

   python fcs2tsv.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import rpy2.robjects as ro
from rpy2.robjects import r as R
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr
import PipelineProject052 as P52

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--fcs-directory", dest="fcs_dir", type="string",
                      help="path to directory containing .fcs files"
                      " for processsing")

    parser.add_option("--output-format", dest="out_format", type="choice",
                      choices=("summary", "intensity"), help="output either"
                      " summary of intensities or raw intensity data")

    parser.add_option("--summary-stats", dest="stats", type="choice",
                      choices=("fano", "mean", "std", "var", "median",
                               "geometric", "regress"),
                      help="summary statistics to output if "
                      "out_format == summary")

    parser.add_option("--gating-strategy", dest="gates", type="string",
                      help=".tsv of gating strategy.  See docs for details")

    parser.add_option("--compensation-matrix", dest="comp_matrix", type="string",
                      help="text file containing the compensation/spillover matrix")

    parser.add_option("--cell-type", dest="cell_type", type="string",
                      help="description of cell type gene expression "
                      "is measured on.  Will be added to output file name")

    parser.add_option("--panel-id", dest="panel", type="string",
                      help="ID for marker panel")

    parser.add_option("--output-directory", dest="out_dir", type="string",
                      help="output directory path for files")

    parser.add_option("--fileset-identifier", dest="fileset_id", type="string",
                      help="unique identifier for sets of files/samples processed "
                      "together.  Useful for assigning to batches for processing")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    if options.cell_type == "CD4_Tmem":
        P52.get_cd4_Tmem_panel3(fcs_dir=options.fcs_dir,
                                out_dir=options.out_dir,
                                comp_matrix=options.comp_matrix,
                                panel=options.panel,
                                setid=options.fileset_id,
                                statistic=options.stats,
                                cell_subset=options.cell_type)
    elif options.cell_type == "CD8_Tmem":
        P52.get_cd8_Tmem_panel3(fcs_dir=options.fcs_dir,
                                out_dir=options.out_dir,
                                comp_matrix=options.comp_matrix,
                                panel=options.panel,
                                setid=options.fileset_id,
                                statistic=options.stats,
                                cell_subset=options.cell_type)
    elif options.cell_type == "CD4_Tnaive":
        P52.get_cd4_naive_panel3(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 statistic=options.stats,
                                 cell_subset=options.cell_type)
    elif options.cell_type == "CD8_Tnaive":
        P52.get_cd8_naive_panel3(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 statistic=options.stats,
                                 cell_subset=options.cell_type)
    elif options.cell_type == "CD4_Tcells":
        P52.get_cd4_tcells_panel1(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  statistic=options.stats,
                                  cell_subset=options.cell_type)
    elif options.cell_type == "CD8_Tcells":
        P52.get_cd8_tcells_panel1(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  statistic=options.stats,
                                  cell_subset=options.cell_type)
    elif options.cell_type == "DN_Tcells":
        P52.get_dn_tcells_panel1(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 statistic=options.stats,
                                 cell_subset=options.cell_type)
    elif options.cell_type == "DP_Tcells":
        P52.get_dp_tcells_panel(fcs_dir=options.fcs_dir,
                                out_dir=options.out_dir,
                                comp_matrix=options.comp_matrix,
                                panel=options.panel,
                                setid=options.fileset_id,
                                statistic=options.stats,
                                cell_subset=options.cell_type)
    else:
        outfile = "%s/%s-%s-%s_%s.tsv" % (options.out_dir, options.fileset_id,
                                          options.panel, options.stats,
                                          options.cell_type)
        P.touch(outfile)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
