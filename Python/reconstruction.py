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
Usage:
    In ipython, run the following:
     run reconstruction.py example.ini mow (or dot or qbi)
     [dsize, data] = read_flt_file('output.flt')
     x = scipy.optimize.fmin_l_bfgs_b(real_spherical_harmonics, [1,1], None, args=(data[:,0],8,2))
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
            self.InputFile = parser.get('Input','diffusion_image')
            print "Reading signals from file", self.InputFile, "..."
            self.dsize, self.Signal = read_flt_file(self.InputFile)
            print "Obtained signals of size", list(self.dsize)

            try:
                self.BMatrixFile = parser.get('Input','b_matrices')
                print "Reading B-matrices from file", self.BMatrixFile, "..."
                self.B = numpy.array(numpy.loadtxt(self.BMatrixFile))
            except:
                self.b = float(parser.get('Input','b_value_fixed'))
                self.DiffusionGradientFile = parser.get('Input','diffusion_gradients')
                print "Reading b-value (%f) and gradients from file %s ..."%(self.b, self.DiffusionGradientFile)
                self.g = numpy.array(numpy.loadtxt(self.DiffusionGradientFile))
                self.B = math.sqrt(self.b)*self.g
            print "Obtained B-matrix of size", self.B.shape

            try:
                self.S0ImageFile = parser.get('Input','baseline_image')
                print "Reading S0 image from", self.S0ImageFile
                self.S0 = read_flt_file(self.S0ImageFile)[1]
            except:
                self.S0 = float(parser.get('Input','baseline_value'))
                print "Use constant S0 value:", self.S0
            self.S0scale = float(parser.get('Input','baseline_scale'))
            self.s = self.S0scale*self.Signal/self.S0  # KxM

            self.r = float(parser.get('Constants','r'))
            self.t = float(parser.get('Constants','t'))

            self.OutputPath = parser.get('Output','output_basename')
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
    """ The class for reconstruction using MOW method.

    Reference:
         Bing Jian, Baba C. Vemuri, Evren Ozarslan, Paul R. Carney, and Thomas H. Mareci
         A novel tensor distribution model for the diffusion-weighted MR signals,
         NeuroImage 37(1), 2007, pp. 164-176.

    TODO:
        (1) implement NNLS/FCNNLS in Python

    """

    def __init__(self,config):
        Reconstructor.__init__(self,config)
        print " ***** Initializing the MOW reconstructor ... ***** "
        parser = ConfigParser.ConfigParser()
        if parser.read(config):
            self.EigenVectorsFile = parser.get('MOW','eigenvectors')
            print "Reading eigenvectors for mixture components from", self.EigenVectorsFile
            self.ev = numpy.array(numpy.loadtxt(self.EigenVectorsFile))
            lambda1 = float(parser.get('MOW','eigenvalue1'))
            lambda2 = float(parser.get('MOW','eigenvalue2'))
            self.ew = [lambda1,lambda2,lambda2]
            print "Eigenvalues:", self.ew
            self.p = float(parser.get('MOW','p'))
            self.damping_factor = float(parser.get('MOW','lambda'))

            self.TessellationFile = parser.get('Input','tessellation')
            print "Reading tessellation vectors from", self.TessellationFile
            self.tessellation = numpy.array(numpy.loadtxt(self.TessellationFile))

            self.SHBasisFile = parser.get('Input','sh_basis_matrix')
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
    """ The class for reconstruction using QBI method.

    References:
    [1] David S. Tuch. Q-ball imaging. Magn. Reson. Med., 52(6):1358-1372, 2004.
    [2] Adam W. Anderson.
        Measurement of Fiber Orientation Distributions Using High Angular Resolution Diffusion Imaging.
        Magn. Reson. Med., 54(5):1194-1206, 2005
    [3] Maxime Descoteaux, Elaine Angelino, Shaun Fitzgibbons, and Rachid Deriche.
        Regularized, Fast and Robust Analytical Q-Ball Imaging.
        Magn. Reson. Med., 58:497-510, 2007.

    """
    def __init__(self,config):
        Reconstructor.__init__(self,config)
        print " ***** Initializing the QBI reconstructor ... ***** "
        parser = ConfigParser.ConfigParser()
        if parser.read(config):
            self.degree = int(parser.get('QBI','degree'))
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
    """ The class for reconstruction using QBI method.

    Reference:
    Evren Özarslan, Timothy M. Shepherd, Baba C. Vemuri, Stephen J. Blackband, and Thomas H. Mareci.
    Resolution of complex tissue microarchitecture using the diffusion orientation transform (DOT).
    NeuroImage, 36(3):1086-1103, 2006.
    """
    def __init__(self,config):
        Reconstructor.__init__(self,config)
        print " ***** Initializing the DOT reconstructor ... *****"
        parser = ConfigParser.ConfigParser()
        if parser.read(config):
            self.degree = int(parser.get('DOT','degree'))
            self.TessellationFile = parser.get('Input','tessellation')
            print "Reading tessellation vectors from", self.TessellationFile
            self.tessellation = numpy.array(numpy.loadtxt(self.TessellationFile))

            self.SHBasisFile = parser.get('Input','sh_basis_matrix')
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
            #print Il
            P += numpy.dot(M[l/2],Il) #TxM
        #print P
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

