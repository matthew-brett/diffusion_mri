#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: reconstruction.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2008/09/19 22:10:59 $
#Version:   $Revision: 1.17 $
#======================================================================

"""
ABOUT:
    This program is a Python implementation of the diffusion MRI reconstruction
    algorithm presented in the following NeuroImage paper,

     Bing Jian, Baba C. Vemuri, Evren Ozarslan, Paul R. Carney, and Thomas H. Mareci
     A novel tensor distribution model for the diffusion-weighted MR signals,
     NeuroImage 37(1), 2007, pp. 164-176

Usage:
    In ipython, run the following:
     run reconstruction.py example.ini mow (or dot or qbi)
     [dsize, data] = read_flt_file('output.flt')
     x = fmin_l_bfgs_b(real_spherical_harmonics, [1,1], None, args=(data[:,0],8,2))

TODO:
    (1) implement NNLS/FCNNLS in Python
"""

import sys
import math
import numpy
import ConfigParser

from utils import *


class Reconstructor:
    """ The base class for reconstruction. """
    def __init__(self,config):
        print " ***** Initializing the base reconstructor ... *****"
        self.config = config
        parser = ConfigParser.ConfigParser()
        if parser.read(config):
            self.InputFile = parser.get('Files','InputFile')
            print "Reading signals from file", self.InputFile, "..."
            self.dsize, self.Signal = read_flt_file(self.InputFile)
            print "Obtained signals of size", list(self.dsize)

            if int(parser.get('Options','UseBMatrix')) == 1:
                self.BMatrixFile = parser.get('Files','BMatrixFile')
                print "Reading B-matrix from file", self.BMatrixFile, "..."
                self.B = numpy.array(numpy.loadtxt(self.BMatrixFile))
            else:
                self.b = float(parser.get('Constants','bValue'))
                self.DiffusionGradientFile = parser.get('Files','DiffusionGradientFile')
                print "Reading b-value (%f) and gradients from file %s ..."%(self.b, self.DiffusionGradientFile)
                self.g = numpy.array(numpy.loadtxt(self.DiffusionGradientFile))
                self.B = math.sqrt(self.b)*self.g
            print "Obtained B-matrix of size", self.B.shape

            if int(parser.get('Options','UseS0Image')) == 1:
                self.S0ImageFile = parser.get('Files','S0Image')
                print "Reading S0 image from", self.S0ImageFile
                self.S0 = read_flt_file(self.S0ImageFile)[1]
            else:
                self.S0 = float(parser.get('Constants','S0default'))
                print "Use constant S0 value:", self.S0
            self.S0scale = float(parser.get('Constants','S0scale'))
            self.s = self.S0scale*self.Signal/self.S0  # KxM

            self.r = float(parser.get('Constants','r'))
            self.t = float(parser.get('Constants','t'))

            self.OutputPath = parser.get('Files','OutputFile')
        else:
            raise InputError('init','No configuration file')

    def write_output(self):
        # assume we already get self.coeff
        print "The spherical harmonics coefficients approximated, size: ", self.coeff.shape
        OutputPath = self.OutputPath
        dsize = self.dsize
        dsize[-1] = self.coeff.shape[0]
        degree =  int((math.sqrt(8*dsize[-1]+1) - 3)/2)
        write_flt_file(OutputPath+'.flt', self.coeff, dsize)
        print "The real-valued spherical harmonics coefficients written to %s.flt"%OutputPath
        [coeff_real, coeff_imag] = convert_coeff(self.coeff,degree)
        dsize[-1] = coeff_real.shape[0]
        write_flt_file(OutputPath+'_real.flt', coeff_real, dsize)
        write_flt_file(OutputPath+'_imag.flt', coeff_imag, dsize)
        print "Complex-valued coefficients written to %s_(real|imag).flt"%OutputPath
        print "To view the reconstruction result, please launch vis_spharm.sav and open %s_real.flt"%OutputPath


class MOWReconstructor(Reconstructor):
    """ The class for reconstruction using MOW method."""
    def __init__(self,config):
        Reconstructor.__init__(self,config)
        print " ***** Initializing the MOW reconstructor ... ***** "
        parser = ConfigParser.ConfigParser()
        if parser.read(config):
            self.EigenVectorsFile = parser.get('Files','EigenVectorsFile')
            print "Reading eigenvectors for mixture components from", self.EigenVectorsFile
            self.ev = numpy.array(numpy.loadtxt(self.EigenVectorsFile))
            lambda1 = float(parser.get('Constants','Lambda1'))
            lambda2 = float(parser.get('Constants','Lambda2'))
            self.ew = [lambda1,lambda2,lambda2]
            print "Eigenvalues:", self.ew
            self.p = float(parser.get('Constants','p'))
            self.damping_factor = float(parser.get('Constants','lambda'))

            self.TessellationFile = parser.get('Files','TessellationFile')
            print "Reading tessellation vectors from", self.TessellationFile
            self.tessellation = numpy.array(numpy.loadtxt(self.TessellationFile))

            self.SHBasisFile = parser.get('Files','SHBasisFile')
            print "Reading the spherical harmonics basis matrix from", self.SHBasisFile
            self.spharm_basis = numpy.array(numpy.loadtxt(self.SHBasisFile))

        else:
            raise InputError('init','No configuration file')


    def assembleA(self):
        ev = self.ev
        ew = self.ew
        p = self.p
        B = prepareB(self.B)
        self.A = assemble_wishart_matrix(B,ev,ew,p)

        return self.A

    def assembleR(self):
        ev = self.ev
        ew = self.ew
        tessellation = self.tessellation
        self.R = assemble_PDF_matrix(tessellation,ev,ew,self.r,self.t)

        return self.R


    def reconstruct(self):

        A = self.assembleA()  #KxN
        print "The system matrix 'A' assembled, size: ", self.A.shape
        alpha = self.damping_factor
        invA = damped_inverse(A, alpha)
        self.w = numpy.dot(invA, self.s);              #NxM
        print "The weights 'w' solved, size: ", self.w.shape
        R = self.assembleR()  #TxN
        print "The PDF matrix R assembled, size: ", self.R.shape
        self.prob = numpy.dot(R,self.w)        #TxM
        print "The probability profiles computed, size: ", self.prob.shape
        inv_basis = damped_inverse(self.spharm_basis, 0)
        self.coeff = numpy.dot(inv_basis, self.prob)
        self.write_output()


class QBIReconstructor(Reconstructor):
    """ The class for reconstruction using QBI method."""
    def __init__(self,config):
        Reconstructor.__init__(self,config)
        print " ***** Initializing the QBI reconstructor ... ***** "
        parser = ConfigParser.ConfigParser()
        if parser.read(config):
            self.degree = int(parser.get('Constants','degree'))
        else:
            raise InputError('init','No configuration file')


    def reconstruct(self):
        print "Construct spherical harmonics basis from the gradients table with degree = ", self.degree
        B = construct_SH_basis(self.g,self.degree)
        invB = damped_inverse(B, 0)
        coeff_for_s = numpy.dot(invB, self.s)
        self.coeff = QBI_coeff(coeff_for_s, self.degree)
        self.write_output()


class DOTReconstructor(Reconstructor):
    """ The class for reconstruction using QBI method."""
    def __init__(self,config):
        Reconstructor.__init__(self,config)
        print " ***** Initializing the DOT reconstructor ... *****"
        parser = ConfigParser.ConfigParser()
        if parser.read(config):
            self.degree = int(parser.get('Constants','degree'))
            self.TessellationFile = parser.get('Files','TessellationFile')
            print "Reading tessellation vectors from", self.TessellationFile
            self.tessellation = numpy.array(numpy.loadtxt(self.TessellationFile))

            self.SHBasisFile = parser.get('Files','SHBasisFile')
            print "Reading the spherical harmonics basis matrix from", self.SHBasisFile
            self.spharm_basis = numpy.array(numpy.loadtxt(self.SHBasisFile))

        else:
            raise InputError('init','No configuration file')


    def reconstruct(self):
        print "Compute M matrices up to degree", self.degree
        M = [computeM(l, self.g, self.tessellation) for l in range(0,self.degree+1,2)] #TxK
        P = numpy.zeros([self.tessellation.shape[0], self.s.shape[1]])
        for l in range(0, self.degree+1,2):
            Il = [[comp_Il(l,S,self.b,self.r,self.t) for S in image] for image in self.s]
            Il = numpy.array(Il).reshape(self.s.shape) #KxM
            P += numpy.dot(M[l/2],Il) #TxM
        self.prob = P
        print "The probability profiles computed, size: ", self.prob.shape
        inv_basis = damped_inverse(self.spharm_basis, 0)
        self.coeff = numpy.dot(inv_basis, self.prob)
        self.write_output()


reconstructor = {}
reconstructor['MOW'] = MOWReconstructor
reconstructor['DOT'] = DOTReconstructor
reconstructor['QBI'] = QBIReconstructor

def main(config_file, method):
    reconstructor[method.upper()](config_file).reconstruct()

if __name__ == '__main__':
    if len(sys.argv)<3:
        print "Usage: reconstuction.py config_file method(mow,dot,qbi)"
    else:
        main(sys.argv[1],sys.argv[2])

