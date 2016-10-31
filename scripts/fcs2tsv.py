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
import CGATPipelines.Pipeline as P

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

    parser.add_option("--database", dest="database", type="string",
                      help="SQLite database to write results to ")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    E.info("Processing data for %s in panel %s" % (options.cell_type,
                                                   options.panel))

    if options.cell_type == "CD4_Tmem" and options.panel == "P3":
        P52.get_cd4_Tmem_panel3(fcs_dir=options.fcs_dir,
                                out_dir=options.out_dir,
                                comp_matrix=options.comp_matrix,
                                panel=options.panel,
                                setid=options.fileset_id,
                                cell_subset=options.cell_type,
                                db=options.database)
    elif options.cell_type == "CD8_Tmem" and options.panel == "P3":
        P52.get_cd8_Tmem_panel3(fcs_dir=options.fcs_dir,
                                out_dir=options.out_dir,
                                comp_matrix=options.comp_matrix,
                                panel=options.panel,
                                setid=options.fileset_id,
                                cell_subset=options.cell_type,
                                db=options.database)
    elif options.cell_type == "CD4_Tnaive" and options.panel == "P3":
        P52.get_cd4_naive_panel3(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "CD8_Tnaive" and options.panel == "P3":
        P52.get_cd8_naive_panel3(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "CD4_Tcell" and options.panel == "P1":
        P52.get_cd4_tcells_panel1(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "CD8_Tcell" and options.panel == "P1":
        P52.get_cd8_tcells_panel1(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "DN_Tcell" and options.panel == "P1":
        P52.get_dn_tcells_panel1(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "DP_Tcell" and options.panel == "P1":
        P52.get_dp_tcells_panel1(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "CD4_Tcell" and options.panel == "P2a":
        P52.get_cd4_tcells_panel2(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "CD4_Tcell" and options.panel == "P2b":
        P52.get_cd4_tcells_panel2(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "CD8_Tcell" and options.panel == "P2a":
        P52.get_cd8_tcells_panel2(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "CD8_Tcell" and options.panel == "P2b":
        P52.get_cd8_tcells_panel2(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "DN_Tcell" and options.panel == "P2a":
        P52.get_dn_tcells_panel2(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "DN_Tcell" and options.panel == "P2b":
        P52.get_dn_tcells_panel2(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "DP_Tcell" and options.panel == "P2a":
        P52.get_dp_tcells_panel2(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "DP_Tcell" and options.panel == "P2b":
        P52.get_dp_tcells_panel2(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    # elif options.cell_type == "early_NK" and options.panel == "P4":
    #     P52.get_early_nkcells_panel4(fcs_dir=options.fcs_dir,
    #                                  out_dir=options.out_dir,
    #                                  comp_matrix=options.comp_matrix,
    #                                  panel=options.panel,
    #                                  setid=options.fileset_id,
    #                                  cell_subset=options.cell_type,
    #                                  db=options.database)
    # elif options.cell_type == "terminal_NK" and options.panel == "P4":
    #     P52.get_terminal_nkcells_panel4(fcs_dir=options.fcs_dir,
    #                                     out_dir=options.out_dir,
    #                                     comp_matrix=options.comp_matrix,
    #                                     panel=options.panel,
    #                                     setid=options.fileset_id,
    #                                     cell_subset=options.cell_type,
    #                                     db=options.database)
    # elif options.cell_type == "mature_NK" and options.panel == "P4":
    #     P52.get_mature_nkcells_panel4(fcs_dir=options.fcs_dir,
    #                                   out_dir=options.out_dir,
    #                                   comp_matrix=options.comp_matrix,
    #                                   panel=options.panel,
    #                                   setid=options.fileset_id,
    #                                   cell_subset=options.cell_type,
    #                                   db=options.database)
    elif options.cell_type == "NKT_early" and options.panel == "P5":
        P52.get_early_nktcells_panel5(fcs_dir=options.fcs_dir,
                                      out_dir=options.out_dir,
                                      comp_matrix=options.comp_matrix,
                                      panel=options.panel,
                                      setid=options.fileset_id,
                                      cell_subset=options.cell_type,
                                      db=options.database)
    elif options.cell_type == "NKT_naive" and options.panel == "P5":
        P52.get_naive_nktcells_panel5(fcs_dir=options.fcs_dir,
                                      out_dir=options.out_dir,
                                      comp_matrix=options.comp_matrix,
                                      panel=options.panel,
                                      setid=options.fileset_id,
                                      cell_subset=options.cell_type,
                                      db=options.database)
    elif options.cell_type == "NKT_terminal" and options.panel == "P5":
        P52.get_terminal_nktcells_panel5(fcs_dir=options.fcs_dir,
                                         out_dir=options.out_dir,
                                         comp_matrix=options.comp_matrix,
                                         panel=options.panel,
                                         setid=options.fileset_id,
                                         cell_subset=options.cell_type,
                                         db=options.database)
    elif options.cell_type == "NKT_effector" and options.panel == "P5":
        P52.get_effector_nktcells_panel5(fcs_dir=options.fcs_dir,
                                         out_dir=options.out_dir,
                                         comp_matrix=options.comp_matrix,
                                         panel=options.panel,
                                         setid=options.fileset_id,
                                         cell_subset=options.cell_type,
                                         db=options.database)
    elif options.cell_type == "Vd1_Tcells" and options.panel == "P5":
        P52.get_Vd1_tcells_panel5(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "Vd2p_Vg9dim" and options.panel == "P5":
        P52.get_Vd2_vg9dim_panel5(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "Vd2n_Vg9p" and options.panel == "P5":
        P52.get_Vd2n_vg9p_panel5(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "hemat_SCs" and options.panel == "P5":
        P52.get_hemat_SC_panel5(fcs_dir=options.fcs_dir,
                                out_dir=options.out_dir,
                                comp_matrix=options.comp_matrix,
                                panel=options.panel,
                                setid=options.fileset_id,
                                cell_subset=options.cell_type,
                                db=options.database)
    elif options.cell_type == "immature_Bcell" and options.panel == "P6":
        P52.get_immature_Bcells_panel6(fcs_dir=options.fcs_dir,
                                       out_dir=options.out_dir,
                                       comp_matrix=options.comp_matrix,
                                       panel=options.panel,
                                       setid=options.fileset_id,
                                       cell_subset=options.cell_type,
                                       db=options.database)
    elif options.cell_type == "mature_Bcell" and options.panel == "P6":
        P52.get_mature_Bcells_panel6(fcs_dir=options.fcs_dir,
                                     out_dir=options.out_dir,
                                     comp_matrix=options.comp_matrix,
                                     panel=options.panel,
                                     setid=options.fileset_id,
                                     cell_subset=options.cell_type,
                                     db=options.database)
    elif options.cell_type == "memory_Bcell" and options.panel == "P6":
        P52.get_memory_Bcells_panel6(fcs_dir=options.fcs_dir,
                                     out_dir=options.out_dir,
                                     comp_matrix=options.comp_matrix,
                                     panel=options.panel,
                                     setid=options.fileset_id,
                                     cell_subset=options.cell_type,
                                     db=options.database)
    elif options.cell_type == "naive_Bcell" and options.panel == "P6":
        P52.get_naive_Bcells_panel6(fcs_dir=options.fcs_dir,
                                    out_dir=options.out_dir,
                                    comp_matrix=options.comp_matrix,
                                    panel=options.panel,
                                    setid=options.fileset_id,
                                    cell_subset=options.cell_type,
                                    db=options.database)
    elif options.cell_type == "IgA_Bmem" and options.panel == "P6":
        P52.get_IgA_Bcells_panel6(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "IgG_Bmem" and options.panel == "P6":
        P52.get_IgG_Bcells_panel6(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "IgM_Bmem" and options.panel == "P6":
        P52.get_IgM_Bcells_panel6(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "IgE_Bmem" and options.panel == "P6":
        P52.get_IgE_Bcells_panel6(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "monocytes" and options.panel == "P7":
        P52.get_monocytes_panel7(fcs_dir=options.fcs_dir,
                                 out_dir=options.out_dir,
                                 comp_matrix=options.comp_matrix,
                                 panel=options.panel,
                                 setid=options.fileset_id,
                                 cell_subset=options.cell_type,
                                 db=options.database)
    elif options.cell_type == "M_dendritic" and options.panel == "P7":
        P52.get_myeloid_DC_panel7(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    elif options.cell_type == "DC_cd123cd11c" and options.panel == "P7":
        P52.get_cd123cd11c_DC_panel7(fcs_dir=options.fcs_dir,
                                     out_dir=options.out_dir,
                                     comp_matrix=options.comp_matrix,
                                     panel=options.panel,
                                     setid=options.fileset_id,
                                     cell_subset=options.cell_type,
                                     db=options.database)
    elif options.cell_type == "P_dendritic" and options.panel == "P7":
        P52.get_plasmacytoid_DC_panel7(fcs_dir=options.fcs_dir,
                                       out_dir=options.out_dir,
                                       comp_matrix=options.comp_matrix,
                                       panel=options.panel,
                                       setid=options.fileset_id,
                                       cell_subset=options.cell_type,
                                       db=options.database)
    elif options.cell_type == "CD8_APCs" and options.panel == "P7":
        P52.get_CD8_APC_panel7(fcs_dir=options.fcs_dir,
                               out_dir=options.out_dir,
                               comp_matrix=options.comp_matrix,
                               panel=options.panel,
                               setid=options.fileset_id,
                               cell_subset=options.cell_type,
                               db=options.database)
    elif options.cell_type == "CD4_APCs" and options.panel == "P7":
        P52.get_CD4_APC_panel7(fcs_dir=options.fcs_dir,
                               out_dir=options.out_dir,
                               comp_matrix=options.comp_matrix,
                               panel=options.panel,
                               setid=options.fileset_id,
                               cell_subset=options.cell_type,
                               db=options.database)
    elif options.cell_type == "CD1cneg_mDC" and options.panel == "P7":
        P52.get_CD1cneg_mDC_panel7(fcs_dir=options.fcs_dir,
                                   out_dir=options.out_dir,
                                   comp_matrix=options.comp_matrix,
                                   panel=options.panel,
                                   setid=options.fileset_id,
                                   cell_subset=options.cell_type,
                                   db=options.database)
    elif options.cell_type == "Inflam_DC" and options.panel == "P7":
        P52.get_inflammatory_DC_panel7(fcs_dir=options.fcs_dir,
                                       out_dir=options.out_dir,
                                       comp_matrix=options.comp_matrix,
                                       panel=options.panel,
                                       setid=options.fileset_id,
                                       cell_subset=options.cell_type,
                                       db=options.database)
    elif options.cell_type == "CD16neg_DC" and options.panel == "P7":
        P52.get_CD16neg_DC_panel7(fcs_dir=options.fcs_dir,
                                  out_dir=options.out_dir,
                                  comp_matrix=options.comp_matrix,
                                  panel=options.panel,
                                  setid=options.fileset_id,
                                  cell_subset=options.cell_type,
                                  db=options.database)
    else:
        outfile = "%s/%s-%s-%s.tsv" % (options.out_dir, options.fileset_id,
                                          options.panel, options.cell_type)
        P.touch(outfile)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
