# this is an attempt to create a Python library for handling, processing and analysing flow cytometry data
# to do list:
# - basic fcs representation
# - group of fcs objects
# - compensation
# - gating functions
# - clustering functions
# - meta data annotations
# - data transformations
# - plotting functions, with gates/filters
# - subsetting
# - identification of cell subsets
# - summary statistics
# - automated gating from a template file

# what is the most basic representation of an fcs file?  A single intensity measurement?

# represent measurements on a single cell - events are an ordered list of the parameter measurements on a single cell

# FCS structure:
# HEADER
# TEXT
# DATA
# ANALYSIS
# OPTIONAL SEGMENTS

# HEADER description:
# first six bytes are the version identifier
# bytes 6-9 are space characters (ASCII 32)
# for FCS2.0 bytes 6-9
# at least 3 pairs of ASCII-encoded integers indicating the byte offsets for the start and end (last byte) of the 
# primary TEXT section, the DATA section and the ANALYSIS section, respectively.

# These byte offsets are relative the beginning of the dataset, limited to 8 bytes in fcs 3.1:

# bytes 10-17 offset to the start of TEXT segment
# bytes 18-25 offset to the end of the TEXT segment

# bytes 26-33 offset to the start of the DATA segment
# bytes 34-41 offset to the end of the DATA segment

# bytes 42-49  offset to the start of the ANALYSIS segment
# bytes 50-57 offset to the end of the ANALYSIS segment

# if these segments are not present these values will be 0
# the HEADER and TEXT segments are required


# struc and file.read to access byte representations
