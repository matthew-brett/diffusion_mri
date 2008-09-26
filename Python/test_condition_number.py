#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: test_condition_number.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2008/09/18 07:00:23 $
#Version:   $Revision: 1.2 $
#======================================================================

import os

import matplotlib
import pylab

from utils import *


def compute_wishart_A(p):
    g = pylab.load('81vectors.txt')
    B = prepareB(math.sqrt(1500.0)*g)
    ew = [0.0015,0.0004,0.0004]
    ev1 = pylab.load('81vectors.txt')
    ev2 = pylab.load('321vectors.txt')
    A1 = assemble_wishart_matrix(B,ev1,ew,p)
    A2 = assemble_wishart_matrix(B,ev2,ew,p)
    return A1,A2


def comp_condA_wishart():
    g = pylab.load('81vectors.txt')
    B = prepareB(math.sqrt(1500.0)*g)
    p_list = [(i+1)*(i+1)+1 for i in range(17)]
    ew = [0.0015,0.0004,0.0004]
    ev1 = pylab.load('81vectors.txt')
    ev2 = pylab.load('321vectors.txt')
    condA_list1 = []
    condA_list2 = []
    for p in p_list:
        A = assemble_wishart_matrix(B,ev1,ew,p)
        U,S,V = numpy.linalg.svd(A)
        condA = S[0]/S[-1]
        condA_list1.append(condA)

        A = assemble_wishart_matrix(B,ev2,ew,p)
        U,S,V = numpy.linalg.svd(A)
        condA = S[0]/S[-1]
        condA_list2.append(condA)

    return p_list,condA_list1,condA_list2

def comp_condA_spike():
    g = pylab.load('81vectors.txt')
    B = prepareB(math.sqrt(1500.0)*g)
    sigma_list = range(10,110,10)
    ew = [1,0,0]
    ev1 = pylab.load('81vectors.txt')
    ev2 = pylab.load('321vectors.txt')
    tessellation = pylab.load('vertices_iso_3.txt')
    condA_list1 = []
    condA_list2 = []
    for sigma in sigma_list:

        A = assemble_spike_matrix(B,ev1,ew,sigma,tessellation)
        U,S,V = numpy.linalg.svd(A)
        condA = S[0]/S[-1]
        condA_list1.append(condA)

        A = assemble_spike_matrix(B,ev2,ew,sigma,tessellation)
        U,S,V = numpy.linalg.svd(A)
        condA = S[0]/S[-1]
        condA_list2.append(condA)

    return sigma_list,condA_list1,condA_list2

def test_cond_wishart():
    #from matplotlib import rc
    #rc('text', usetex=True)
    pylab.rcParams['axes.formatter.limits'] = (-5,5)
    p,c1,c2 = comp_condA_wishart()
    pylab.plot(p,c1,'ro--',label='81x81',linewidth=2)
    pylab.plot(p,c2,'bs-',label='81x321',linewidth=2)
    #pylab.title('condition numbers of assembled matrices',fontsize=20)
    pylab.xlabel('p in Wishart distributions',fontsize=24,fontweight='bold')
    pylab.ylabel('condition numbers, cond(A)',fontsize=24,fontweight='bold')
    fontprop = matplotlib.font_manager.FontProperties(size=24,weight='bold')
    pylab.legend(('size(A):81x81', 'size(A): 81x321'),loc='center right',prop = fontprop)
    for label in pylab.axes().get_xticklabels():
        label.set_fontsize(16)
        label.set_fontweight('bold')
    for label in pylab.axes().get_yticklabels():
        label.set_fontsize(16)
        label.set_fontweight('bold')

    pylab.grid(False)
    pylab.show()
    if not os.path.isdir('figures'):
        os.mkdir('figures')
    pylab.savefig('./figures/cond_wishart.png',dpi=300)
    pylab.savefig('./figures/cond_wishart.eps',dpi=300)


def test_cond_spike():
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', weight='bold')
    pylab.rcParams['axes.formatter.limits'] = (-5,5)
    p,c1,c2 = comp_condA_spike()
    pylab.plot(p,c1,'ro--',label='81x81',linewidth=2)
    pylab.plot(p,c2,'bs-',label='81x321',linewidth=2)
    #pylab.title('condition numbers of assembled matrices',fontsize=20)
    pylab.xlabel(r'$\sigma$ in radial basis functions',fontsize=24,fontweight='bold')
    pylab.ylabel('condition numbers, cond(A)',fontsize=24,fontweight='bold')
    fontprop = matplotlib.font_manager.FontProperties(size=24,weight='bold')
    pylab.legend(('size(A):81x81', 'size(A): 81x321'),loc='upper center',prop = fontprop)

    ## try pylab.setp()
    for label in pylab.axes().get_xticklabels():
        label.set_fontsize(16)
        label.set_fontweight('bold')
    for label in pylab.axes().get_yticklabels():
        label.set_fontsize(16)
        label.set_fontweight('bold')

    pylab.grid(False)
    pylab.show()
    if not os.path.isdir('figures'):
        os.mkdir('figures')
    pylab.savefig('./figures/cond_spike.png',dpi=300)
    pylab.savefig('./figures/cond_spike.eps',dpi=300)


