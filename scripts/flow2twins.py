'''
flow2twins.py - file manipulation, table merging, etc for twins data
===============================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. merging and manipulating files related to twin data

Usage
-----

.. Example use case

Example::

   python flow2twins.py

Type::

   python flow2twins.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
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

    parser.add_option("--twin-id-column", dest="id_column", type="string",
                      help="column number or header for twin IDs in "
                      "input tables")

    parser.add_option("--demographics-file", dest="demo_file", type="string",
                      help="tab-separated text file containing twins "
                      "demographic data")

    parser.add_option("--demo-id-column", dest="demo_id_column", type="string",
                      help="column header or number that indicates the twin IDs "
                      "that match to those from the flow cytometry data")

    parser.add_option("--task", dest="task", type="choice",
                      choices=["merge_flow", "split_zygosity",
                               "regress_confounding"],
                      help="choose a task")

    parser.add_option("--output-file-pattern", dest="out_pattern", type="string",
                      help="output filename pattern to use where output is "
                      "multiple files")

    parser.add_option("--zygosity-column", dest="zygosity_col", type="string",
                      help="column header containing zygosity information")

    parser.add_option("--id-columns", dest="id_headers", type="string",
                      help="comma-separated list of column headers used to "
                      "uniquely identify each sample")

    parser.add_option("--family-id-column", dest="family_id", type="string",
                      help="column header containing family IDs")

    parser.add_option("--confounding-column", dest="confounding", type="string",
                      help="either a comma-separates list or single value, "
                      "column number or column header")

    parser.add_option("--marker-group", dest="marker_col", type="string",
                      help="column header containing marker IDs to group "
                      "regression fits by")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    if options.task == "merge_flow":
        list_of_files = infile.split(",")

        # test using CD4+ Tmem cells, all .fcs files from FlowRepository.org
        out_df = P52.merge_flow_tables(file_list=list_of_files,
                                       id_column=options.id_column,
                                       demo_file=options.demo_file,
                                       demo_id_column=options.demo_id_column)

        out_df.to_csv(options.stdout, sep="\t", index_col="indx")

    elif options.task == "split_zygosity":
        out_frames = P52.split_zygosity(infile=infile,
                                        zygosity_column=options.zygosity_col,
                                        id_headers=(options.id_headers).split(","),
                                        pair_header=options.family_id)
        # expect keys: MZ and DZ
        MZ_frame = out_frames["MZ"]
        DZ_frame = out_frames["DZ"]

        # output filenames using pattern and zygosity as a prefix
        MZ_outfile = "-".join([options.out_pattern, "MZ.tsv"])
        MZ_frame.to_csv(MZ_outfile, sep="\t", index_label="indx")

        DZ_outfile = "-".join([options.out_pattern,"DZ.tsv"])
        DZ_frame.to_csv(DZ_outfile, sep="\t", index_label="indx")

    elif options.task == "regress_confounding":
        out_df = P52.regress_out_confounding(infile=infile,
                                             confounding_column=options.confounding,
                                             group_var=options.marker_col)

        out_df.to_csv(options.stdout, sep="\t", index_label="indx")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
