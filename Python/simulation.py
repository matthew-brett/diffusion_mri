#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: simulation.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2009/06/21 18:11:17 $
#Version:   $Revision: 1.7 $
#======================================================================

import math
import random
import numpy
import scipy.special


BesselJprimeRoots = numpy.array([
   [0.00000000000000,   0.00000000000000,   0.00000000000000,   0.00000000000000,   0.00000000000000,   0.00000000000000,  0.00000000000000,  0.00000000000000,   0.00000000000000,   0.00000000000000],
   [0.00000000000000,   3.83170597020751,   7.01558666981562,   10.17346813506272,  13.32369193631422,  16.47063005087763, 19.61585851046824, 22.76008438059277,  25.90367208761838,  29.04682853491686],
   [1.84118378134066,   5.33144277352503,   8.53631636634629,   11.70600490259206,  14.86358863390903,  18.01552786268180, 21.16436985918879, 27.45705057105924,  30.60192297266909,  33.74618289866739],
   [3.05423692822714,   6.70613319415846,   9.96946782308760,   13.17037085601612,  16.34752231832178,  19.51291278248821, 22.67158177247742, 25.82603714178526,  28.97767277299368,  32.12732702044347],
   [4.20118894121053,   8.01523659837595,   11.34592431074301,  14.58584828616703,  20.97247693653770,  27.31005793020435, 30.47026880629042, 33.62694918279668,  36.78102067546439,  39.93310862365949],
   [5.31755312608399,   9.28239628524161,   12.68190844263889,  15.96410703773155,  19.19602880004891,  22.40103226768900, 25.58975968138674, 28.76783621766650,  31.93853934097278,  35.10391667734677],
   [6.41561637570024,   10.51986087377231,  13.98718863014030,  17.31284248788462,  20.57551452138689,  23.80358147659386, 27.01030789777772, 30.20284907898166,  33.38544390101012,  36.56077768688036],
   [7.50126614468415,   11.73493595304271,  15.26818146109787,  18.63744300966620,  21.93171501780224,  25.18392559949963, 28.40977636251008, 31.61787571610504,  34.81339298429743,  37.99964089771530],
   [8.57783648971408,   12.93238623708958,  16.52936588436694,  19.94185336652734,  23.26805292645757,  26.54503206182358, 29.79074858319661, 33.01517864137514,  36.22438054878716,  39.42227457893926],
   [9.64742165199722,   14.11551890789462,  17.77401236691525,  21.22906262285312,  24.58719748631768,  27.88926942795509, 31.15532655618833, 34.39662855427218,  37.62007804419709,  40.83017868182204],
   [10.71143397069995,  15.28673766733295,  19.00459353794605,  22.50139872677729,  25.89127727683914,  29.21856349993608, 32.50524735237553, 35.76379292880880,  39.00190281151422,  42.22463843075328]
])



def Jprime(n, x):
    _jn = scipy.special.jn
    if n<=0:
        return -_jn(1,x)
    else:
        return 0.5*(_jn(n-1,x) - _jn(n+1,x))


def K(n,m):
    if (n==0 and m==0):
        return 1.0
    elif (n==0 or m==0):
        return 2.0
    else:
        return 4.0

SQR = lambda x: x*x
PI = math.pi

def Alpha(k,m):
    return BesselJprimeRoots[m+1][k-1]

def Theta(q, v):
    q = q/numpy.linalg.norm(q)
    v = v/numpy.linalg.norm(v)
    qv = numpy.fabs(numpy.dot(q,v))
    if qv >= 1.0:
        return 0.0
    else:
        return numpy.arccos(qv)

def Minus1PowN(n):
    if (n%2):
        return -1.0
    else:
        return 1.0


Delta = 17.8E-3    #17.8ms
delta = 2.2E-3     #2.2ms
D = 2.0E-3         #2.0E-3 mm^2/s;
b = 1500.0         #1500.0f s/mm^2;
S0 = 60.0

qMag = math.sqrt(b/(4.0*SQR(PI)*(Delta-delta/3.0)))

def Attenuation(r, l, q, theta):
    eps = 1.0e-5
    if(numpy.fabs(theta) < eps):
        theta = eps
    if(numpy.fabs(theta-PI/2.0) < eps):
        theta = PI/2.0 - eps

    sTheta = numpy.sin(theta)
    s2Theta = numpy.sin(2.0*theta)
    cTheta = numpy.cos(theta)

    nmax = 50    #1000
    kmax = 10     #10
    mmax = 10     #10

    _2PIqr = 2.0*PI*q*r
    _2PIql = 2.0*PI*q*l

    s = 0.0

    kn = 2.0*SQR(r)*SQR(SQR(_2PIqr))*SQR(s2Theta)
    kd = SQR(l)
    #print kn,kd

    c = numpy.zeros(mmax)
    for m in range(mmax):
        c[m] = SQR(Jprime(m, _2PIqr*sTheta))

    an = numpy.zeros([nmax,mmax])
    ad = numpy.zeros(nmax)
    bn = numpy.zeros(nmax)
    for n in range(nmax):
        ad[n] = SQR(SQR(n*PI*r/l)-SQR(_2PIqr*cTheta))
        bn[n] = (1.0-Minus1PowN(n)*numpy.cos(_2PIql*cTheta))
        for m in range(mmax):
            an[n,m] = K(n,m)

    bd = numpy.zeros([kmax,mmax])
    f = numpy.zeros([kmax,mmax])
    for k in range(kmax):
        for m in range(mmax):
            bd[k,m] = SQR( SQR(Alpha(k+1,m))-SQR(_2PIqr*sTheta))
            if(m==0):
                f[k,m] = 1.0
            else:
                f[k,m] = SQR(Alpha(k+1,m))/(SQR(Alpha(k+1,m))-SQR(m))


    for n in range(nmax):
        for k in range(kmax):
            for m in range(mmax):
                d = math.exp(-(SQR(Alpha(k+1,m)/r)+SQR(n*PI/l))*D*Delta)
                s += (an[n,m]/ad[n])*(bn[n]/bd[k,m])*f[k,m]*c[m]*d
    s *= kn/kd
    return s


def gen_1fib(tessellation_file, gradient_file):
    #file1 = '46vectors.txt'
    #file2 = '81vectors.txt'

    r = 0.005
    l = 0.015
    q = qMag

    v = numpy.loadtxt(tessellation_file)
    g = numpy.loadtxt(gradient_file)

    n = len(v)
    if n%3 == 0:
        dsize = [n/3,3,len(g)]
    elif n%2 == 0:
        dsize = [n/2,2,len(g)]
    else:
        dsize = [n,1,len(g)]

    data = numpy.zeros([len(g),n])
    for j in range(len(v)):
        print j
        for i in range(len(g)):
            theta = Theta(v[j],g[i])
            data[i,j] += Attenuation(r, l, q, theta)

    return data,dsize

import flt_utils
def gen_basis(tessellation,gradient):
    file1 = '%dvectors.txt'%tessellation
    file2 = '%dvectors.txt'%gradient
    data, dsize = gen_1fib(file1, file2)
    filename = 's_t%d_g%d_f1.flt'%(tessellation,gradient)
    flt_utils.write_flt_file(filename,data,dsize)
    return data,dsize,filename

def generate_k_subsets(n,k):
    """
    Generate all k-subsets of an n-set sequentially.
    Reference:
        http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
        Chapter 3
    """
    m, h  = 0, k
    a = range(k)
    yield a
    while True:
        if m < n-h: h = 1
        else: h += 1
        m = a[k-h] + 1
        for j in range(h): a[k+j-h] = m + j
        yield a
        if a[0] == n - k:  break


def gen_2fib(basis, file1, threshold):
    v = numpy.loadtxt(file1)
    n = len(v)
    x = generate_k_subsets(n,2)
    vv = []
    for i,j in x:
        if Theta(v[i],v[j])>threshold:
            vv.append((i,j))
    print len(vv)
    s = numpy.zeros([basis.shape[0],len(vv)])
    for k in range(len(vv)):
        print k
        i,j = vv[k]
        s[:,k] = 0.5*(basis[:,i] + basis[:,j])
    return s,vv,v




def add_rician_noise_to_file(filename, stds, shuffle = False):
    #stds = numpy.arange(0.01,0.1,0.01)
    flt_filename = filename + '.flt'
    dsize, data = flt_utils.read_flt_file(flt_filename)
    inds = numpy.arange(data.shape[1])
    for std in stds:
        s_noise = add_rician_noise(data, std)
        k,n = s_noise.shape
        if shuffle:
            numpy.random.shuffle(inds)
        output_filename = '%s_sd%03d.flt'%(filename,std*100)
        txt_filename = '%s_sd%03d.txt'%(filename,std*100)
        flt_utils.write_flt_file(output_filename,s_noise[:,inds],dsize)
        f = open(txt_filename,'w')
        f.write('\n'.join([str(ind) for ind in inds]))
        f.close()




