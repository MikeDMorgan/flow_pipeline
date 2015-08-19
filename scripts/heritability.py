
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

    parser.add_option("--heritability", dest="equation", type="string",
                      help="equation used to estimate heritability")

    parser.add_option("--monozygote-file", dest="mz_file", type="string",
                      help="file containing monozygotic twin data")

    parser.add_option("--dizygote-file", dest="dz_file", type="string",
                      help="file containing dizygotic twin data")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    h_df = P52.estimate_heritability(mz_file=options.mz_file,
                                     dz_file=options.dz_file)
    options.stdout.write("H^2\n")
    for key in h_df.keys():
        options.stdout.write("%s: %0.3f\n" % (key, h_df[key]))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
