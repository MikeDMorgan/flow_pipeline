'''
data2distance.py - calculate pair-wise distances between data matrices
===============================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. distances and clustering of data matrices

Usage
-----

.. Example use case

Example::

   python data2distance.py

Type::

   python data2distance.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import PipelineProject052 as P52
import numpy as np
import pandas as pd


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

    parser.add_option("--id-column", dest="id_col", type="string",
                      help="column header for sample IDs")

    parser.add_option("--matrix-distance", dest="dist", type="choice",
                      choices=["Euclid"], help="distance metric to use "
                      "for distance between matrices")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    parser.set_defaults(filter_zero=True,
                        filt_gate=None)

    infile = argv[-1]
    # the input file is a list of files
    myfunc = lambda x: x.rstrip("\n")

    with open(infile, "r") as ifile:
        all_files = ifile.readlines()

    list_of_files = map(myfunc, all_files)

    matrix_list = P52.makeMatrixList(list_of_files=list_of_files,
                                     id_col=options.id_col)
    distance_matrix = P52.getMatrixDistances(list_of_matrices=matrix_list,
                                             distance=options.dist)

    distance_matrix.to_csv(options.stdout, sep="\t", index_label=None)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
