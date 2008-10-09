#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: mhd_utils.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2008/10/09 06:19:47 $
#Version:   $Revision: 1.2 $
#======================================================================

import os
import numpy
import array

from utils import *

def read_meta_header(filename):
    """Return a dictionary of meta data from meta header file"""
    fileIN = open(filename, "r")
    line = fileIN.readline()

    meta_dict = {}
    tag_set1 = ['ObjectType','NDims','DimSize','ElementType','ElementDataFile']
    tag_set2 = ['BinaryData','BinaryDataByteOrderMSB','CompressedData','CompressedDataSize']
    tag_set3 = ['Offset','CenterOfRotation','AnatomicalOrientation','ElementSpacing','TransformMatrix']
    tag_set4 = ['Comment','SeriesDescription','AcquisitionDate','AcquisitionTime','StudyDate','StudyTime']
    tag_set = []

    tag_set.extend(tag_set1)
    tag_set.extend(tag_set2)
    tag_set.extend(tag_set3)
    tag_set.extend(tag_set4)
    tag_flag = [False]*len(tag_set)
    while line:
        tags = str.split(line,'=')
        #print tags[0]
        for i in range(len(tag_set)):
            tag = tag_set[i]
            if (str.strip(tags[0]) == tag) and (not tag_flag[i]):
                #print tags[1]
                meta_dict[tag] = str.strip(tags[1])
                tag_flag[i] = True
        line = fileIN.readline()
    #print comment
    fileIN.close()
    return meta_dict

def load_raw_data_with_mhd(filename):
    meta_dict = read_meta_header(filename)
    dim = int(meta_dict['NDims'])
    #print dim
    #print meta_dict['ElementType']
    assert(meta_dict['ElementType']=='MET_FLOAT')
    arr = [int(i) for i in meta_dict['DimSize'].split()]
    #print arr
    volume = reduce(lambda x,y: x*y, arr[0:dim-1], 1)
    #print volume
    pwd = os.path.split(filename)[0]
    if pwd:
        data_file = pwd +'/' + meta_dict['ElementDataFile']
    else:
        data_file = meta_dict['ElementDataFile']
    #print data_file
    fid = open(data_file,'rb')
    binvalues = array.array('f')
    binvalues.read(fid, volume*arr[dim-1])
    if is_little_endian(): # assume data in file is always big endian
        binvalues.byteswap()
    fid.close()
    data = numpy.array(binvalues, numpy.float)
    data = numpy.reshape(data, (arr[dim-1], volume))
    return (data, meta_dict)


def write_meta_header(filename, meta_dict):
    header = ''
    # do not use tags = meta_dict.keys() because the order of tags matters
    tags = ['ObjectType','NDims','BinaryData',
       'BinaryDataByteOrderMSB','CompressedData','CompressedDataSize',
       'TransformMatrix','Offset','CenterOfRotation',
       'AnatomicalOrientation',
       'ElementSpacing',
       'DimSize',
       'ElementType',
       'ElementDataFile',
       'Comment','SeriesDescription','AcquisitionDate','AcquisitionTime','StudyDate','StudyTime']
    for tag in tags:
        if tag in meta_dict.keys():
            header += '%s = %s\n'%(tag,meta_dict[tag])
    f = open(filename,'w')
    f.write(header)
    f.close()


def dump_raw_data(filename, data):
    """ Write the data into a raw format file. Big endian is always used. """
    rawfile = open(filename,'wb')
    a = array.array('f')
    for o in data:
        a.fromlist(list(o))
    if is_little_endian():
        a.byteswap()
    a.tofile(rawfile)
    rawfile.close()


def write_mhd_file(mhdfile, data, dsize):
    assert(mhdfile[-4:]=='.mhd')
    meta_dict = {}
    meta_dict['ObjectType'] = 'Image'
    meta_dict['BinaryData'] = 'True'
    meta_dict['BinaryDataByteOrderMSB'] = 'False'
    meta_dict['ElementType'] = 'MET_FLOAT'
    meta_dict['NDims'] = str(len(dsize))
    meta_dict['DimSize'] = ' '.join([str(i) for i in dsize])
    meta_dict['ElementDataFile'] = os.path.split(mhdfile)[1].replace('.mhd','.raw')
    write_meta_header(mhdfile, meta_dict)

    pwd = os.path.split(mhdfile)[0]
    if pwd:
        data_file = pwd +'/' + meta_dict['ElementDataFile']
    else:
        data_file = meta_dict['ElementDataFile']

    dump_raw_data(data_file, data)


