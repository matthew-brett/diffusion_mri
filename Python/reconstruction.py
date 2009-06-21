#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: reconstruction.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2009/06/21 03:38:44 $
#Version:   $Revision: 1.9 $
#======================================================================

"""
Usage:
    In ipython, go to the data directory and run the following:
     run ../Python/reconstruction.py example.ini mow (or dot or qbi)
"""

import os
import sys
import math
import numpy
import ConfigParser

from utils import *
from mhd_utils import *
from flt_utils import *
from lsqnonneg import lsqnonneg

class Reconstructor:
    """ The base class for reconstruction. """

    def __init__(self):
        print " ***** Initializing the base reconstructor ... *****"

    def load_signal_from_file(self, input_file):
        self.InputFile = input_file
        print "Reading signals from file", self.InputFile, "..."
        if os.path.splitext(self.InputFile)[-1] == '.flt':
            self.dsize, self.Signal = read_flt_file(self.InputFile)
        else:
            self.Signal, meta_dict = load_raw_data_with_mhd(self.InputFile)
            self.dsize = [int(i) for i in meta_dict['DimSize'].split()]
        print "Obtained signals of size", list(self.dsize)

    def load_bmatrix_from_file(self, input_file):
        self.BMatrixFile = input_file
        print "Reading B-matrices from file", self.BMatrixFile, "..."
        self.B = numpy.array(numpy.loadtxt(self.BMatrixFile))

    def load_bmatrix_from_gradient(self, b, gradient_file):
        self.b = b
        self.DiffusionGradientFile = gradient_file
        print "Reading b-value (%f) and gradients from file %s ..."%(self.b, self.DiffusionGradientFile)
        self.g = numpy.array(numpy.loadtxt(self.DiffusionGradientFile))
        self.B = math.sqrt(self.b)*self.g

    def load_S0_from_file(self, input_file):
        self.S0ImageFile = input_file
        print "Reading S0 image from", self.S0ImageFile
        if os.path.splitext(self.S0ImageFile)[-1] == '.flt':
            self.S0 = read_flt_file(self.S0ImageFile)[1]
        else:
            self.S0 = load_raw_data_with_mhd(self.S0ImageFile)[0]

    def init_from_configuration(self, parser):
        self.config = parser
        input_file = parser.get('Input','diffusion_image')
        self.load_signal_from_file(input_file)
        #print self.Signal
        try:
            self.load_bmatrix_from_file(parser.get('Input','b_matrices'))
        except:
            b = float(parser.get('Input','b_value_fixed'))
            gradient_file = parser.get('Input','diffusion_gradients')
            self.load_bmatrix_from_gradient(b, gradient_file)
        print "Obtained B-matrix of size", self.B.shape

        try:
            self.load_S0_from_file(parser.get('Input','baseline_image'))
        except:
            self.S0 = float(parser.get('Input','baseline_value'))
            print "Use constant S0 value:", self.S0
        self.S0scale = float(parser.get('Input','baseline_scale'))
        self.s = self.S0scale*self.Signal/self.S0  # KxM
        self.r = float(parser.get('Constants','r'))
        self.t = float(parser.get('Constants','t'))
        self.OutputPath = parser.get('Output','output_basename')

    def load_tessellation(self, input_file):
        self.TessellationFile = input_file
        print "Reading tessellation vectors from", self.TessellationFile
        self.tessellation = numpy.array(numpy.loadtxt(self.TessellationFile))

    def load_spharm_basis(self, input_file):
        self.SHBasisFile = input_file
        print "Reading the spherical harmonics basis matrix from", self.SHBasisFile
        self.spharm_basis = numpy.array(numpy.loadtxt(self.SHBasisFile))
        self.inv_basis = damped_inverse(self.spharm_basis, 0)

    def write_output(self):
        # assume we already get self.coeff
        print "The spherical harmonics coefficients approximated, size: ", self.coeff.shape
        OutputPath = self.OutputPath
        dsize = self.dsize
        dsize[-1] = self.coeff.shape[0]
        degree =  int((math.sqrt(8*dsize[-1]+1) - 3)/2)
        write_mhd_file(OutputPath+'.mhd', self.coeff, dsize)
        write_flt_file(OutputPath+'.flt', self.coeff, dsize)
        cr,ci = convert_coeff(self.coeff, degree)
        write_flt_file(OutputPath+'_real.flt', cr, dsize)
        write_flt_file(OutputPath+'_imag.flt', ci, dsize)
        print "The real-valued spherical harmonics coefficients written to %s.mhd"%OutputPath



class MOWReconstructor(Reconstructor):
    """ The class for reconstruction using MOW method.

    Reference:
         Bing Jian, Baba C. Vemuri, Evren Ozarslan, Paul R. Carney, and Thomas H. Mareci
         A novel tensor distribution model for the diffusion-weighted MR signals,
         NeuroImage 37(1), 2007, pp. 164-176.

    TODO:
        (1) implement FCNNLS in Python
    """
    def __init__(self):
        Reconstructor.__init__(self)
        print " ***** Initializing the MOW reconstructor ... ***** "

    def init_from_configuration(self, parser):
        Reconstructor.init_from_configuration(self, parser)
        self.load_eigenvectors(parser.get('MOW','eigenvectors'))
        lambda1 = float(parser.get('MOW','eigenvalue1'))
        lambda2 = float(parser.get('MOW','eigenvalue2'))
        self.ew = [lambda1,lambda2,lambda2]
        print "Eigenvalues:", self.ew
        self.p = float(parser.get('MOW','p'))
        self.damping_factor = float(parser.get('MOW','lambda'))
        self.load_tessellation(parser.get('Input','tessellation'))
        self.load_spharm_basis(parser.get('Input','sh_basis_matrix'))

        try:
            self.solver = parser.get('MOW','solver')
        except:
            self.solver = 'nnls'
        self.assembleA()
        self.assembleR()

    def load_eigenvectors(self, input_file):
        self.EigenVectorsFile = input_file
        print "Reading eigenvectors for mixture components from", self.EigenVectorsFile
        self.ev = numpy.array(numpy.loadtxt(self.EigenVectorsFile))


    def assembleA(self):
        ev = self.ev
        ew = self.ew
        p = self.p
        B = prepareB(self.B)
        self.A = assemble_wishart_matrix(B,ev,ew,p)
        print "The system matrix 'A' assembled, size: ", self.A.shape
        #return self.A

    def assembleR(self):
        ev = self.ev
        ew = self.ew
        tessellation = self.tessellation
        self.R = assemble_PDF_matrix(tessellation,ev,ew,self.r,self.t)
        print "The PDF matrix R assembled, size: ", self.R.shape
        #return self.R

    def reconstruct(self):
        A = self.A  #KxN

        if self.solver == 'dls':
            alpha = self.damping_factor
            invA = damped_inverse(A, alpha)
            self.w = numpy.dot(invA, self.s);              #NxM
        else: #if self.solver == 'nnls':
            self.w = numpy.zeros((A.shape[1],self.s.shape[1]))
            for i in range(self.s.shape[1]):
                self.w[:,i] = lsqnonneg(A, self.s[:,i])[0]

        print "The weights 'w' solved, size: ", self.w.shape
        R = self.R  #TxN
        self.prob = numpy.dot(R,self.w)        #TxM
        print "The probability profiles computed, size: ", self.prob.shape
        self.coeff = numpy.dot(self.inv_basis, self.prob)


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
    def __init__(self):
        Reconstructor.__init__(self)
        print " ***** Initializing the QBI reconstructor ... ***** "

    def init_from_configuration(self, parser):
        Reconstructor.init_from_configuration(self, parser)
        self.degree = int(parser.get('QBI','degree'))
        try:
        #if True:
            self.use_spline = parser.get('QBI','use_spline')
            self.NodeFile = parser.get('QBI','nodes')
            self.nodes = numpy.array(numpy.loadtxt(self.NodeFile))

            self.FunkRadonBasisFile = parser.get('QBI','funk_radon_basis')
            self.funk_radon_basis = numpy.array(numpy.loadtxt(self.FunkRadonBasisFile))
            # funk_radon_basis.shape = (#_of_tessellation, #_of_nodes)

            self.TessellationFile = parser.get('Input','tessellation')
            print "Reading tessellation vectors from", self.TessellationFile
            self.tessellation = numpy.array(numpy.loadtxt(self.TessellationFile))

            self.SHBasisFile = parser.get('Input','sh_basis_matrix')
            print "Reading the spherical harmonics basis matrix from", self.SHBasisFile
            self.spharm_basis = numpy.array(numpy.loadtxt(self.SHBasisFile))

        except:
        #else:
            self.use_spline = ''

        if self.use_spline:
            print " ***** Computing the QBI-ODF using spherical splines ... ***** "
            import spherical_splines
            _R = spherical_splines.assemble_kernel_matrix(self.g, self.nodes)
            _one = numpy.ones([81,1])
            _eye = numpy.eye(81)
            _lambda = 0
            self.A = numpy.r_[numpy.c_[_one,_R + _lambda * _eye], numpy.c_[0, _one.T]]
            self.inv_basis = damped_inverse(self.spharm_basis, 0)
        else:
            print "Construct spherical harmonics basis from the gradients table with degree = ", self.degree
            B = construct_SH_basis(self.g,self.degree)
            self.invB = damped_inverse(B, 0)


    def reconstruct(self):
        if self.use_spline:
            A = self.A
            s = self.s
            c = numpy.dot(numpy.linalg.pinv(A),numpy.r_[s,[[0]]])
            self.odf = numpy.dot(self.funk_radon_basis, c[1::])
            self.coeff = numpy.dot(self.inv_basis, self.odf)
        else:
            coeff_for_s = numpy.dot(self.invB, self.s)
            self.coeff = QBI_coeff(coeff_for_s, self.degree)


class DOTReconstructor(Reconstructor):
    """ The class for reconstruction using the DOT method.

    Reference:
    Evren Ozarslan, Timothy M. Shepherd, Baba C. Vemuri, Stephen J. Blackband, and Thomas H. Mareci.
    Resolution of complex tissue microarchitecture using the diffusion orientation transform (DOT).
    NeuroImage, 36(3):1086-1103, 2006.
    """
    def __init__(self):
        Reconstructor.__init__(self)
        print " ***** Initializing the DOT reconstructor ... ***** "

    def init_from_configuration(self, parser):
        Reconstructor.init_from_configuration(self, parser)
        self.degree = int(parser.get('DOT','degree'))
        self.load_tessellation(parser.get('Input','tessellation'))
        self.load_spharm_basis(parser.get('Input','sh_basis_matrix'))

    def comp_M(self):
        print "Compute M matrices up to degree", self.degree
        self.M = [computeM(l, self.g, self.tessellation) for l in range(0,self.degree+1,2)] #TxK

    def reconstruct(self):
        P = numpy.zeros([self.tessellation.shape[0], self.s.shape[1]])
        for l in range(0, self.degree+1,2):
            Il = [[comp_Il(l,S,self.b,self.r,self.t) for S in image] for image in self.s]
            Il = numpy.array(Il).reshape(self.s.shape) #KxM
            #print Il
            P += numpy.dot(self.M[l/2],Il) #TxM
        #print P
        self.prob = P
        print "The probability profiles computed, size: ", self.prob.shape
        self.coeff = numpy.dot(self.inv_basis, self.prob)


reconstructor = {}
reconstructor['MOW'] = MOWReconstructor
reconstructor['DOT'] = DOTReconstructor
reconstructor['QBI'] = QBIReconstructor

def main(config_file, method):
    parser = ConfigParser.ConfigParser()
    if parser.read(config_file):
        recon = reconstructor[method.upper()](parser)
        recon.reconstruct()
        recon.write_output()
    else:
        raise InputError('init','No configuration file')


if __name__ == '__main__':
    if len(sys.argv)<3:
        print "Usage: reconstuction.py config_file method(mow|dot|qbi)"
    else:
        main(sys.argv[1],sys.argv[2])

