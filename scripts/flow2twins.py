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
                               "regress_confounding", "kinship"],
                      help="choose a task")

    parser.add_option("--output-file-pattern", dest="out_pattern", type="string",
                      help="output filename pattern to use where output is "
                      "multiple files")

    parser.add_option("--output-directory", dest="out_dir", type="string",
                      help="directory to output files into")

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

    parser.add_option("--filter-zero-arrays", dest="filter_zero", action="store_true",
                      help="Filter out arrays where there are no observations")

    parser.add_option("--database", dest="database", type="string",
                      help="absolute path to SQLite database")

    parser.add_option("--tablename", dest="table", type="string",
                      help="tablename to extract from SQL database")

    parser.add_option("--filter-gates", dest="filt_gate", type="string",
                      help="regex to filter out unwanted triboolean gates")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    parser.set_defaults(filter_zero=True,
                        filt_gate=None)

    infile = argv[-1]

    if options.task == "merge_flow":

        merged_arrays = P52.mergeGateArrays(db=options.database,
                                            table_name=options.table,
                                            filter_zero=options.filter_zero,
                                            filter_gate=options.filt_gate)
        
        all_dfs = P52.mergeArrayWithDemographics(flow_arrays=merged_arrays,
                                                 id_column=options.id_column,
                                                 demo_file=options.demo_file,
                                                 demo_id_column=options.demo_id_column)

        for out_df in all_dfs:
            # construct the table names using cell_type, panel and gate
            # cell type and statistic should be in the out_pattern
            outname = "-".join([options.out_pattern, out_df["gate"][0],
                                out_df["panel"][0]])
            out_file = "/".join([options.out_dir, outname])
            E.info("writing %s data to file" % outname)
            out_df.to_csv(out_file, sep="\t", index_col="indx")

    elif options.task == "split_zygosity":
        out_frames = P52.split_zygosity(infile=infile,
                                        zygosity_column=options.zygosity_col,
                                        id_headers=(options.id_headers).split(","),
                                        pair_header=options.family_id)
        # expect keys: MZ and DZ
        try:
            MZ_frame = out_frames["MZ"]
            DZ_frame = out_frames["DZ"]          

            # output filenames using pattern and zygosity as a prefix
            MZ_outfile = "-".join([options.out_pattern, "MZ.tsv"])
            MZ_frame.to_csv(MZ_outfile, sep="\t", index_label="indx")

            DZ_outfile = "-".join([options.out_pattern,"DZ.tsv"])
            DZ_frame.to_csv(DZ_outfile, sep="\t", index_label="indx")
        except TypeError:
            pass

    elif options.task == "regress_confounding":
        out_df = P52.regress_out_confounding(infile=infile,
                                             confounding_column=options.confounding,
                                             group_var=options.marker_col)

        out_df.to_csv(options.stdout, sep="\t", index_label="indx")
    elif options.task == "kinship":
        out_df = P52.make_kinship_matrix(twins_file=infile,
                                         id_column=options.id_headers,
                                         family_column=options.family_id,
                                         zygosity_column=options.zygosity_col)

        out_df.to_csv(options.stdout, index_label="indx", sep="\t")


    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
