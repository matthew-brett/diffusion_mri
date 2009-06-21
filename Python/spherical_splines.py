#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: spherical_splines.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2009/04/09 06:04:19 $
#Version:   $Revision: 1.8 $
#=====================================================================

from numpy import arange, array, dot, ones, zeros, eye, r_, c_, linalg, loadtxt, max
from math import log, sqrt, fabs, pi
import os

def q(z):
    """Evaluate q(z) = \int_0^1 (1-h)^2/sqrt(1-2*h*z+h*h)} dh """
    if z>=1.0: return 0.5

    """
    Grace Wahba, Spline interpolation and smoothing on the sphere,
    SIAM J. SCI. STAT. COMPUT., 2(1) 1981, pp. 5-16. [Equation (3.4)]
    http://www.stat.wisc.edu/%7Ewahba/ftp1/oldie/sphspl.pdf

    W = (1-z)/2.0
    C = 2*sqrt(W)
    A = log(1+1.0/sqrt(W))
    return 0.5*(A*(12*W*W - 4*W) -6*C*W + 6*W + 1)
    """

    """
    H. J. Taijeron, A. G. Gibson, C. Chandler,
    Spline interpolation and smoothing on hyperspheres,
    SIAM J. SCI. COMPUT., 15(5) 1994, pp. 1111-1125. [Table 1]
    """
    try:
        S = sqrt(2-2*z)
        N = (1-z)*log(2-2*z)
        L = (1-z)*log(sqrt(2/(1-z))+1)
        return 0.5*(-L*(3*z-1)+3*S*(z-1)+4-3*z)
    except:
        return 0.5


def sanity_test():
    """A sanity test compared with numerical approximation. """
    dh = .0001
    hh = arange(dh/2.0,1,dh)
    g = lambda z: sum([dh*(1-h)**2/sqrt(1-2*h*z+h*h) for h in hh])
    print max([fabs(q(z)-g(z)) for z in arange(0,1,0.01)]) #8.33332181038e-010

# The reproducing kernel: see Wahba (3.3) and Taijeron et al. (45)
R = lambda z: (0.25*(q(z)+q(-z)) - 1/6.0)/(2*pi)

def assemble_kernel_matrix(x,v):
    """Assemble the kernel matrix from a given set of directions. """
    #use R(fabs(dot(i,j))) in case angles are considered in the range (0,pi/2)
    return array([R(dot(i,j)) for i in x for j in v]).reshape(len(x),len(v))

def test(gradient_file, signal_file, knots_file, _lambda = 0):
    """
    Fit the spherical thin-plate spline model to descrete signal data.
    Reference:
    Ferreira et al. Directional Log-Spline Distributions,
    Bayesian Analysis 3(2) 2008, pp. 297-316 [Eq.(3)]

    In [331]: c = test('81vectors.txt','3fib.mhd')
    2.92138710735e-013

    In [332]: c = test('81vectors.txt','2fib.mhd')
    2.37209920248e-013

    In [333]: c = test('81vectors.txt','1fib.mhd')
    3.90984974495e-013

    """
    import mhd_utils
    import flt_utils
    g = loadtxt(gradient_file) #diffusion-mri/Python/data/81vectors.txt
    v = loadtxt(knots_file)
    #vv = r_[v,-v]
    #gg = r_[g,-g]
    _R = assemble_kernel_matrix(g, v)
    #_one = ones([81,1])
    #_eye = eye(81)
    #_R = assemble_kernel_matrix(gg, vv)
    _one_column = ones([len(g),1])
    _one_row = ones([1,len(v)])
    #_eye = eye(81*2)
    #A = r_[c_[_one,_R + _lambda * _eye], c_[0, _one.T]]
    A = r_[c_[_one_column,_R], c_[0, _one_row]]
    basename = os.path.splitext(signal_file)[1]
    if basename == '.flt':
        [dsize,s] = flt_utils.read_flt_file(signal_file)
    elif basename == '.mhd':
        [s,dsize] = mhd_utils.load_raw_data_with_mhd(signal_file)
    else:
        return

    _zero_row = zeros([1,s.shape[1]])
    c = dot(linalg.pinv(A),r_[s,_zero_row])
    #return A,s,c,_R
    #c = dot(linalg.pinv(A),r_[s,s,[[0]]])
    print max(abs(c[0] + dot(_R,c[1::]) - s))
    return A,c,s,g,_R
