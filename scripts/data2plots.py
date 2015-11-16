'''
data2plots.py - plot generation from mean expression and expression noise analysis
==================================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. distances and clustering of data matrices

Usage
-----
Order of operations:
 * merge dataframes
 * melt dataframe
 * check plotting parameters
 * generate plot

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

    parser.add_option("--plot-type", dest="plot_type", type="choice",
                      choices=["histogram", "scatter", "barchart"],
                      help="type of plot to generate")

    parser.add_option("--x-axis", dest="x_axis", type="string",
                      help="variable to plot on the X-axis."
                      "This is the default axis for plotting.")

    parser.add_option("--y-axis", dest="y_axis", type="string",
                      help="variable to plot on the Y-axis")

    parser.add_option("--split-by", dest="split_by", type="string",
                      help="varible over which to split up plots")

    parser.add_option("--X-title", dest="x_title", type="string",
                      help="label to attach to X-axis")

    parser.add_option("--Y-title", dest="y_title", type="string",
                      help="label to attach to Y-axis")

    parser.add_option("--colour-var", dest="col_var", type="string",
                      help="variable to colour points by")

    parser.add_option("--free-scale", dest="free_scale", type="choice",
                      choices=["free_x", "free_y", "free"], help="whether to "
                      "use free scaling on plot axes")

    parser.add_option("--outfile", dest="outfile", type="string",
                      help="file to save plot to")

    parser.add_option("--melt-data", dest="melt", action="store_true",
                      help="melt the dataframe first, requires ID vars")

    parser.add_option("--melt-id-vars", dest="id_vars", type="string",
                      help="comma separated list of id variables for"
                      " the melted dataframe")

    parser.add_option("--merge-frames", dest="merge", action="store_true",
                      help="merge two input dataframes together")

    parser.add_option("--merge-id-vars", dest="merge_vars", type="string",
                      help="comma separate list of id variables to merge "
                      "two dataframes on")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    parser.set_defaults(free_scale="both",
                        split_by=None,
                        col_var=None,
                        melt=False)

    infile = argv[-1]

    if len(infile.split(",")) == 2:
        infiles = infile.split(",")
        df1 = pd.read_table(infiles[0], sep=":", index_col=0,
                            header=0)
        df2 = pd.read_table(infiles[1], sep=":", index_col=0,
                            header=0)
        ids = options.merge_vars.split(",")
        df1[options.y_axis] = df1.index
        df2[options.y_axis] = df2.index

        df1.columns = [ids[0], options.y_axis]
        df2.columns = [ids[0], options.y_axis]

        df = pd.merge(df1, df2, on=options.y_axis)

        # these need to not be hard-coded!
        df.columns = ["mean_h2", ids[0], "fano_h2"]
    else:
        df = pd.read_table(infile, sep="\t", index_col=0,
                           header=0)

    # assumes the first column is the index
    if options.melt:
        mids = options.id_vars.split(",")
        if options.y_axis:
            _df = pd.melt(df, id_vars=options.x_axis, value_name=options.y_axis,
                          var_name=options.col_var)
        else:
            _df = pd.melt(df, id_vars=mids, value_name=options.x_axis,
                          var_name=options.col_var)
        df = _df
    else:
        pass

    # check variables are present
    try:
        var = df[options.x_axis]
    except ValueError:
        raise ValueError("no plotting variable found")

    if options.col_var:
        try:
            cols = df[options.col_var]
        except ValueError:
            E.warn("Colour variable not found in data frame."
                   "Check the data file is the correct one")
    else:
        pass

    if options.split_by:
        try:
            splits = df[options.split_by]
        except ValueError:
            E.warn("Split-by variable not found in the data "
                   "frame.  Check the data file is the correct"
                   " one.")
    else:
        pass

    try:
        assert options.outfile
    except AssertionError:
        raise IOError("no output file detected")

    if options.plot_type == "histogram":
        P52.plotHistogram(data=df,
                          variable=options.x_axis,
                          save_path=options.outfile,
                          x_title=options.x_title,
                          y_title=options.y_title,
                          colour_var=options.col_var,
                          scales=options.free_scale,
                          split_var=options.split_by)
    elif options.plot_type == "barchart":
        P52.plotBarchart(data=df,
                         x_variable=options.x_axis,
                         y_variable=options.y_axis,
                         save_path=options.outfile,
                         x_title=options.x_title,
                         y_title=options.y_title,
                         colour_var=options.col_var,
                         split_var=options.split_by)
    else:
        pass

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
