'''
workspace2tsv.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Convert XML workspace files to various text files
by parsing and removing relevant components

Usage
-----

.. Example use case

Example::

   python workspace2tsv.py

Type::

   python workspace2tsv.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
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

    parser.add_option("--method", dest="method", type="choice",
                      choices=["compensation", "parse_gating"],
                      help="select method to perform on workspace "
                      "file.")

    parser.add_option("--gating-directory", dest="gate_dir", type="string",
                      help="directory to store gating dummy files")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # write footer and output benchmark information.
    E.Stop()

    infile = argv[-1]

    if options.method == "compensation":
        split_file = infile.split("/")
        infile = split_file[-1]
        split_file.remove(infile)
        path = "/".join(split_file)
        out_df = P52.get_compensation_matrix(path=path,
                                             infile=infile)
        out_df.to_csv(options.stdout, sep="\t")

    elif options.method == "parse_gating":
        for dfile in P52.parse_gating_file(infile):
            outfile = options.gate_dir + "/" + dfile
            P.touch(outfile)
    else:
        pass

if __name__ == "__main__":
    sys.exit(main(sys.argv))
