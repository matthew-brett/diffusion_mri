#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: utils.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2008/10/09 06:19:47 $
#Version:   $Revision: 1.5 $
#======================================================================

import math
import numpy
import array
import scipy.special
import spherical_harmonics

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message



def is_little_endian():
    """
    determine the endianess of the local machine
    http://developer.nicta.com.au/pipermail/elefant-dev/2007-May/000005.html
    """
    from struct import pack
    if pack('<h', 1) == pack('=h',1):
        return True
    elif pack('>h', 1) == pack('=h',1):
        return False
    else:
        print "Unknown Endian"
        return True



def convert_coeff(coeff, degree, dl=2):
    """
    Convert real-valued coefficients to complex-valued coefficients in
    the format used by the IDL visualization program vis_spharm.sav
    """

    src_coeff_dim, n = coeff.shape
    assert(2*src_coeff_dim>=degree*(degree-1))
    assert(degree%2==0)

    if dl==2:
        des_coeff_dim = (degree/2 + 1)**2
    elif dl==1:
        des_coeff_dim = degree*(degree+1)/2
    else:
        return None

    coeff_imag = numpy.zeros([des_coeff_dim, n])
    coeff_real = numpy.zeros([des_coeff_dim, n])

    i = 0
    for l in range(0, degree+1, dl):
        if dl==2:
            center = (l+1)*(l+2)/2 - l - 1
        elif dl==1:
            center = (l+1)*(l+1) - l - 1

        coeff_real[i] = coeff[center] #m = 0
        coeff_imag[i] = 0
        i += 1

        for m in range(1,l+1):
            if m%2:
                coeff_imag[i] = coeff[center-m]
                coeff_real[i] = -coeff[center+m]
            else:
                coeff_imag[i] = -coeff[center-m]
                coeff_real[i] = coeff[center+m]
            coeff_imag[i] *= math.sqrt(2.0)/2.0
            coeff_real[i] *= math.sqrt(2.0)/2.0
            i += 1

    return coeff_real, coeff_imag


def damped_inverse(A, alpha):
    [U,W,V] = numpy.linalg.svd(A,False)
    invW = [w/(w*w+alpha*alpha) for w in W]
    invA = numpy.dot(V.T,numpy.dot(numpy.diag(invW),U.T))
    return invA

def prepareB(B):
    K = B.shape[0]  # assert(B.shape[1] = 6 or 3)
    if B.shape[1]==3:
        cellB = [numpy.outer(B[j],B[j]) for j in range(K)]
    elif B.shape[1]==6:
        cellB = []
        for j in range(K):
            cellB.append(numpy.zeros([3,3], float))
            cellB[j][1-1,1-1] = B[j,1-1];
            cellB[j][2-1,2-1] = B[j,2-1];
            cellB[j][3-1,3-1] = B[j,3-1];
            cellB[j][1-1,2-1] = B[j,4-1];
            cellB[j][1-1,3-1] = B[j,5-1];
            cellB[j][2-1,3-1] = B[j,6-1];
            cellB[j][2-1,1-1] = B[j,4-1];
            cellB[j][3-1,1-1] = B[j,5-1];
            cellB[j][3-1,2-1] = B[j,6-1];
    else:
        raise InputError('assembleA','B.shape[1] should be either 3 or 6')
    return cellB;


def outer_product_matrice(ev):
    return [numpy.outer(v,v) for v in ev]

def assemble_prolate_tensors(ev, ew):
    """
    Construct rotationally symmetric tensors from a set of dominant
    eigenvectors (ev) and fixed eigenvalues (ew)
    """
    return [(ew[0]-ew[1])*numpy.outer(v,v)+ew[1]*numpy.eye(3) for v in ev]

def assemble_inverse_matrice(ev, ew):
    return [1.0/ew[1]*numpy.eye(3) - (ew[0]-ew[1])/(ew[0]*ew[1]) *numpy.outer(v,v) for v in ev]


def assemble_wishart_matrix(B,ev,ew,p):
    """
    Assemble a matrix A from B-matrix and a set of SPD matrices with fixed
    eigenvalues (EW) and varying eigenvectors (EV).
    """
    D = assemble_prolate_tensors(ev, ew)
    if p>0: # wishart basis
        A = [math.pow(1+numpy.trace(numpy.dot(b,d))/p, -p) for b in B for d in D]
    else: # reduce to tensor basis
        A = [math.exp(-numpy.trace(numpy.dot(b,d))) for b in B for d in D]
    A = numpy.array(A).reshape(len(B),len(D))
    return A


def assemble_spike_matrix(B,ev,ew,sigma,tessellation):
    """
    Construct the matrix used in Alexander, IPMI05.
    """
    K = len(B)
    T = tessellation.shape[0]
    N = ev.shape[0]
    v = tessellation
    R = numpy.zeros([K,T], float)
    F = numpy.zeros([T,N], float)
    for j in range(K):
        for i in range(T):
            R[j,i] = math.exp(-sum(ew)*numpy.dot(numpy.dot(v[i],B[j]),v[i]))
    for i in range(T):
        for j in range(N):
            dot_prod = math.acos( min(math.fabs(numpy.dot(v[i],ev[j])),1) )
            F[i,j] = math.exp(-(dot_prod*dot_prod)/(sigma*sigma))
    A = numpy.dot(R,F)/T
    return A


def assemble_PDF_matrix(tessellation, ev, ew, r, t):

    factor = r*r/(4*t)
    detD = ew[0]*ew[1]*ew[1]

    B = outer_product_matrice(tessellation)
    D = assemble_inverse_matrice(ev, ew)

    R = [numpy.exp(-factor*numpy.trace(numpy.dot(b,d))) for b in B for d in D]
    R = numpy.array(R).reshape(len(B),len(D))/math.sqrt(detD)
    return R

"""
function [az,elev,r] = cart2sph(x,y,z)
%CART2SPH Transform Cartesian to spherical coordinates.
%   [TH,PHI,R] = CART2SPH(X,Y,Z) transforms corresponding elements of
%   data stored in Cartesian coordinates X,Y,Z to spherical
%   coordinates (azimuth TH, elevation PHI, and radius R).  The arrays
%   X,Y, and Z must be the same size (or any of them can be scalar).
%   TH and PHI are returned in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PHI is the elevation angle from the xy plane.
%
%   Class support for inputs X,Y,Z:
%      float: double, single
%
%   See also CART2POL, SPH2CART, POL2CART.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2008/10/09 06:19:47 $

"""

def cart2sph(x,y,z):
    r = math.sqrt(x**2+y**2+z**2)
    elev = math.atan2(z,math.sqrt(x**2+y**2))
    az = math.atan2(y,x)
    return az, elev, r


def construct_SH_basis(pts, degree):
    sph_coord = numpy.array([cart2sph(x,y,z) for x,y,z in pts])
    B  = [spherical_harmonics.evaluate_SH((math.pi/2 - elev,az),degree,2)[0] for az,elev,r in sph_coord]
    # NOTE: spherical harmonics' convention: PHI - azimuth, THETA - polar angle
    # PI/2 - elev is to convert elevation [-PI/2, PI/2] to polar angle [0, PI]
    B = numpy.array(B).squeeze()
    return B


"""
       for l=2:2:degree
            oddl = prod(3:2:l-1);
            evenl = prod(2:2:l);
            center =(l + 2)*(l + 1)/2 - l - 1
            ll = oddl/evenl;
            if (mod((l/2),2))
                ll = -ll;
            end
            for m=-l:l
                coeff_symm(i, center+m) = coeff_symm(i,center+m)*ll ;
            end %% m-loop
        end  %% l-loop
"""

def QBI_coeff(coeff, degree):
    for l in range(2,degree+1,2):
        oddl = numpy.arange(3,l,2).prod()
        evenl = numpy.arange(2,l+1,2).prod()
        center =(l + 2)*(l + 1)/2 - l - 1
        ll = oddl*1.0/evenl
        if (l/2)%2:
            ll = -ll
        for m in range(-l,l+1):
            coeff[center+m] *= ll
    return coeff




def comp_Il(l,S,b,R0,t):
    D = -math.log(S)/b
    #D(D<0) = eps;  %% in case S>1
    beta = R0 / math.sqrt(D*t)


    expB = math.exp(-0.25*(beta**2))
    erfB = scipy.special.erf(0.5*beta)
    f2 = erfB/(4.0*math.pi*R0**3)
    powD = math.pow(4.0*math.pi*D*t, 1.5)
    f1 = expB/powD

    def g0(beta):
        return 1,0

    def g2(beta):
        return -(1+6/beta**2), 3

    def g4(beta):
        Al = 1.0 + 20.0/(beta**2) + 210.0/(beta**4)
        Bl = (15.0/2.0)*(1.0 - 14.0/(beta**2))
        return Al,Bl

    def g6(beta):
        Al = -(1.0 + 42.0/(beta**2) + (1575.0/2.0)/(beta**4) + 10395.0/(beta**6))
        Bl = (105.0/8.0)*(1.0 - 36.0/(beta**2) + 396.0/(beta**4))
        return Al,Bl

    def g8(beta):
        Al = 1.0 + 72.0/(beta**2) + (10395.0/4.0)/(beta**4) + 45045.0/(beta**6) + 675675.0/(beta**8)
        Bl = (315.0/16.0)*(1.0 - 66.0/(beta**2) + 1716.0/(beta**4) - 17160.0/(beta**6))
        return Al,Bl

    g = [g0,g2,g4,g6,g8]
    Al,Bl = g[l/2](beta)
    Il = Al*f1 + Bl*f2
    return Il


def computeM(l, gradients, tessellation):
    # assert(gradients.shape[1]==3)
    # assert(tessellation.shape[1]==3)
    r_dot_u = numpy.array([numpy.dot(r,u) for r in tessellation for u in gradients])
    r_dot_u[numpy.where(r_dot_u < -1)[0]] = -1.0
    r_dot_u[numpy.where(r_dot_u >  1)[0]] =  1.0
    Pl = [spherical_harmonics.P(l,0,math.acos(cos_theta))[0] for cos_theta in r_dot_u]
    #Pl = [scipy.special.legendre(l)(cos_theta) for cos_theta in r_dot_u] is much slower
    M =  (-1)**(l/2)*(2*l+1)*numpy.array(Pl).reshape(len(tessellation),len(gradients))
    return M



