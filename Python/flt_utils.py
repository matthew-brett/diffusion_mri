#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: flt_utils.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2008/10/27 05:55:55 $
#Version:   $Revision: 1.1 $
#======================================================================

import os
import numpy
import array

from utils import *
from mhd_utils import *


def read_flt_file(filename):
    """
    Read the diffusion signals from a flt file.

    Usage:
        (dsize, data) = read_flt_file(filename)
    please see example below for explanation.

    A flt file starts with a header structure
    +-------------------------------------------------------------------+
    | dim (int) | size (int[dim]) | data_type (int) | data_length (int) |
    +-------------------------------------------------------------------+
    followed by the binary data section.

    Notes:
    (1) The number of dimensions is stored in the first 4 bytes as integer type.
    The range of dimension should be from 1 to 4. It is very likely that the byte
    order used by the file is big-endian rather than the little-endian on your PC.
    So swapping may be required to get the right number.
    (2) The current program only recognizes float data type which is encoded as 4.
    For complex data, the data_type field should be 8.

    Example:
    If 'data.flt' contains 81 diffusion MR images, each of which has size 128x128x64,
    then the header organzied in 4-bytes-unit should look like
         4,  128, 128, 64, 81, 4, 128*128*64*81
    and right after the header are 128*128*64*81 float numbers stored by stacking 81
    float-valued images of size 128*128*64.
    In this case, (dsize,data) = read_flt_file('data.flt') will return
      (1) 'dsize' as a 1-D int array [128,128,64,81] and
      (2) 'data' as a 2-D float array of shape [81,128*128*64]
    """

    fid = open(filename,'rb')
    arr = array.array('i')
    arr.fromfile(fid, 1) # dim
    dim = arr[0]
    #http://www.python.org/search/hypermail/python-1993/0393.html
    if dim>100:
        """print 'Read very high dimension (>100).'
        print 'Endianness may come into play.'
        print 'Try to swap the byte order.'"""
        swap = True;
        arr.byteswap()
        dim = arr[0]
        #print 'dim =',dim
    else:
        swap = False
    assert(dim>=1 and dim<=4) # only accept data up to 4 dimensions.

    arr = array.array('i')
    arr.fromfile(fid,dim+2)
    if swap:
        arr.byteswap()
    volume = reduce(lambda x,y: x*y, arr[0:dim-1], 1)

    binvalues = array.array('f')
    binvalues.read(fid, volume*arr[dim-1])
    if swap:
        binvalues.byteswap()
    fid.close()

    data = numpy.array(binvalues, numpy.float)
    data = numpy.reshape(data, (arr[dim-1], volume))

    return (arr[:dim],data)


def write_flt_file(filename, data, dsize):
    """ Write the data into a flt format file. Big endian is always used. """
    binfile = open(filename,'wb')

    dsize = numpy.array(dsize)
    dsize[-1] = data.shape[0]

    header = [len(dsize)]  # dimension
    header.extend(list(dsize)) # size
    header.append(4)       # data type: float
    header.append(dsize.prod())  # total length of data

    a = array.array('i')
    a.fromlist(header)
    if is_little_endian():
        a.byteswap()

    a.tofile(binfile)

    a = array.array('f')
    for o in data:
        a.fromlist(list(o))
    if is_little_endian():
        a.byteswap()
    a.tofile(binfile)
    binfile.close()


def from_flt_to_mhd(fltfile, mhdfile):
    dsize,data = read_flt_file(fltfile)
    write_mhd_file(mhdfile, data, dsize)


def from_mhd_to_flt(mhdfile, fltfile):
    data, meta_dict = load_raw_data_with_mhd(mhdfile)
    dsize = [int(i) for i in meta_dict['DimSize'].split()]
    write_flt_file(fltfile, data, dsize)
