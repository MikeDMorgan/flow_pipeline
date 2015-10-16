
'''
heritability.py - heritability calculations
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python heritability.py

Type::

   python heritability.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import PipelineProject052 as P52
from rpy2.robjects import r as R
from rpy2.robjects import rinterface
from rpy2.robjects import pandas2ri


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

    parser.add_option("--task", dest="task", type="choice",
                      choices=["heritability", "merge",
                               "kinship"],
                      help="task")

    parser.add_option("--heritability", dest="equation", type="string",
                      help="equation used to estimate heritability")

    parser.add_option("--monozygote-file", dest="mz_file", type="string",
                      help="file containing monozygotic twin data")

    parser.add_option("--dizygote-file", dest="dz_file", type="string",
                      help="file containing dizygotic twin data")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.task == "heritability":
        h_df = P52.estimate_heritability(mz_file=options.mz_file,
                                         dz_file=options.dz_file)

        # generate scatter plots of each marker
        outdir = "/".join(options.dz_file.split("/")[:-1])
        plot_out = "-".join(options.dz_file.split("/")[-1].split("-")[1:4])
        plot_out = "/".join([outdir, plot_out])
        E.info("plotting correlations to %s" % plot_out)
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        R('''mz.df <- read.table("%s", sep="\t", h=T, row.names=1)''' % options.mz_file)
        R('''dz.df <- read.table("%s", sep="\t", h=T, row.names=1)''' % options.dz_file)
        R('''all.dz <- data.frame(rbind(mz.df, dz.df))''')
        R('''p_cor <- ggplot(all.dz, aes(x=twin1, y=twin2, '''
          '''colour=zygosity)) + '''
          '''geom_point(size=1) + stat_smooth(method=lm) + '''
          '''facet_wrap( ~ marker, scales="free")''')
        R('''png("%s-cors.png", height=720, width=720)''' % plot_out)
        R('''print(p_cor)''')
        R('''dev.off()''')

        options.stdout.write("H^2\n")
        for key in h_df.keys():
            options.stdout.write("%s: %0.3f\n" % (key, h_df[key]))
    elif options.task == "merge":
        infiles = argv[-1]
        infiles = infiles.split(",")
        out_df = P52.merge_heritability(infiles)
        out_df.to_csv(options.stdout, sep="\t")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
