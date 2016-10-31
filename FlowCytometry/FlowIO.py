###########################################################
# Flow cytometry IO classes and functions
###########################################################
'''
Have an object/data structure to represent an FCS file.

Each .fcs file is a single sample, an experiment can be
represented by a collection of .fcs files.

FCS file attributes::
  * header
  * text
  * data
  * analysis
  * optional segments??

Have separate objects that represent these different sections?



'''
import struct
import numpy as np
import sys
import codecs

class FCS(object):
    '''
    The base class representation of an fcs file.  Follows the
    ISAC specification for FCS 3.1 format
    Uses UTF-8 throughout for keywords

    Numerical values are base 10
   
    attributes
    ----------
    version - FCS file version
    header_size - bytes for header, version dependent?
    text_size - bytes for TEXT segment (tuple)
    data_size - bytes for DATA segment (tuple)
    analysis_size - bytes for ANALYSIS segment (tuple)
    key_words - dict of keyword: value pairs
    data - a numpy ndarray of dimension (nEvents, nParams)

    ANALYSIS - ??
    '''

    def __init__(self, file_handle):
        '''
        Instantiate the FCS object
        '''
        
        # read and parse the file upon instantiation
        # accept a file handle or a file path

        if type(file_handle) == str:
            self.ofile = open(file_handle, "r")
        elif type(file_handle) == file:
            self.ofile = file
        else:
            raise IOError("input file must be either a path "
                          "to a valid file, or an open file "
                          "connection.")
        
        # version is the first 6 bytes of the file
        self.version = self.ofile.read(6)
        self._get_fcs_segments()
        self._get_text_segment()
        self._get_data_segment()

    def _get_fcs_segments(self):
        '''
        Get byte offsets for file segments 
        and assign version
        '''

        # bytes 6-9 are space
        # read and discard them
        dump = self.ofile.read(4)
        t_start = int(self.ofile.read(8))
        t_end = int(self.ofile.read(8))
        self.text_size = (t_start, t_end)

        d_start = int(self.ofile.read(8))
        d_end = int(self.ofile.read(8))
        self.data_size = (d_start, d_end)
        
        a_start = int(self.ofile.read(8))
        a_end = int(self.ofile.read(8))
        self.analysis_size = (a_start, a_end)

        # use all of the first 57 bytes

    def _get_text_segment(self):
        '''
        Using the attributes `text_size`,
        parse the TEXT segment and assign
        all keyword: value pairs as a dictionary
        '''

        # read all the bytes to the end, then parse
        # the keyword: value pairs
        # the delimiter is at the start and end of
        # each pair, required keywords start with '$'

        self.ofile.seek(self.text_size[0], 0)
        TEXT = self.ofile.read(self.text_size[1] - self.text_size[0] + 1)
        delim = TEXT[0]
        # check delimiter bounds the header text
        if TEXT[-1] == delim:
            pass
        else:
            raise ValueError("The file parser expects the header section to be bounded by "
                             "the same delimiter.  This file may be corrupt")

        # remove the first and last characters (delimiters)
        text = TEXT[1:-1].split(delim)
        kwords, vals =  text[0::2], text[1::2]
        self.key_words = dict(zip(kwords, vals))

    def _get_data_segment(self):
        '''
        Use the data_size attributes to parse
        the DATA segment and assign as a
        collections object
        '''
        byte_map = {"4,3,2,1": ">",
                    "1,2,3,4": "<"}

        data_map = {"I": "I",
                    "F": "f",
                    "D": "d",
                    "A": "b"}

        dump = self.ofile.read(self.data_size[0] - self.text_size[1])

        DATA = self.ofile.read(self.data_size[1] - self.data_size[0])
        # need to unpack the data
        # use the keyword $BYTEORD

        try:
            byte_order = self.key_words["$BYTEORD"]
        except KeyError:
            raise IOError("The $BYTEORD keyword is missing, check "
                          "the input FCS file is complete and valid")

        try:
            dtype = data_map[self.key_words["$DATATYPE"]]
        except KeyError:
            raise IOError("The $DATATYPE keyword is missing, check "
                          "the input FCS is complete and valid")

        try:
            self.key_words["$PAR"]
        except KeyError:
            raise IOError("The $PAR keyword is missing, check "
                          "the input FCS is complete and valid")

        try:
            self.key_words["$TOT"]
        except KeyError:
            raise IOError("The $TOT keyword is missing, check "
                          "the input FCS is complete and valid")

        dsize = int(self.key_words["$PAR"]) * int(self.key_words["$TOT"])

        # need byte order, bits per parameter and bit range if using
        # unsigned integers and number of parameters
        frmt = "%s%i%s" % (byte_map[byte_order],
                           dsize,
                           dtype)

        # DEBUG: need to figure out correcr size to pass to struct.unpack
        print self.key_words["$BEGINDATA"], self.key_words["$ENDDATA"], self.data_size
        data = np.array(struct.unpack(frmt, DATA))
        data.shape = (int(self.key_words["$PAR"]),
                      int(self.key_words["$TOT"]))

        self.data = data
