#!/usr/bin/env python
#coding=utf-8

#======================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: ipmi07_exp.py,v $
#Language:  Python
#Author:    $Author: bjian $
#Date:      $Date: 2009/06/21 18:11:17 $
#Version:   $Revision: 1.5 $
#======================================================================

import ConfigParser
import numpy
from pyx import *
from utils import *
#from flt_utils import *
import reconstruction
#from simulation import *


"""
To repeat the Fig5 in the ipmi07 paper:
   cd Python
   import ipmi07_exp
   cd ../data
   ipmi07_exp.save_result('2fib_result.txt')
   ipmi07_exp.plot_errorbar('2fib_result.txt','ipmi07_fig5')
"""

reconstructor = {}
reconstructor['MOW'] = reconstruction.MOWReconstructor
reconstructor['DOT'] = reconstruction.DOTReconstructor
reconstructor['QBI'] = reconstruction.QBIReconstructor
methods = reconstructor.keys()



def init_reconstructor(method):
    recon = reconstructor[method.upper()]()
    b = 1500.0
    gradient_file = '81vectors.txt'
    recon.load_bmatrix_from_gradient(b, gradient_file)
    recon.r = 0.015 #### radius of the probability surface, measured in mm = 1000um = 1000micron
    recon.t = 0.020 #### diffusion time, measured in seconds
    return recon

def prepare_mow(recon):
    recon.load_eigenvectors('321vectors.txt')
    recon.ew = [0.0015, 0.0004, 0.0004]
    recon.p = 0 #2
    recon.load_tessellation('321vectors.txt')
    recon.load_spharm_basis('basis_321_8.txt')
    recon.assembleA()
    recon.assembleR()
    recon.solver = 'nnls'

def prepare_qbi(recon):
    recon.use_spline = ''
    recon.degree = 8
    B = construct_SH_basis(recon.g, recon.degree)
    recon.invB = damped_inverse(B, 0)

def prepare_dot(recon):
    recon.degree = 8
    recon.load_tessellation('321vectors.txt')
    recon.load_spharm_basis('basis_321_8.txt')
    recon.comp_M()

prepare = {}
prepare['MOW'] = prepare_mow
prepare['DOT'] = prepare_dot
prepare['QBI'] = prepare_qbi


def expt_2fib(method):
    recon = init_reconstructor(method)
    prepare[method.upper()](recon)
    init_angles = (((90,30),), ((90,20), (90,100)))
    init_angles = deg2rad(init_angles)
    sigma = numpy.linspace(0.01,0.10,10)

    s = numpy.loadtxt('data.txt')[:,1] # columns 0,1,2 => 1-,2-,3-fiber
    n = 100
    s_noiseless = numpy.outer(s, numpy.ones(n)) # make n copies
    init_pos = init_angles[1] # ((90,20), (90,100))
    deviations = {}
    for sig in sigma:
        ## option 1: read noisy data from saved files
        #noisyfilename = '%s/%s_sd%03d.flt'%(data_path, root_name, int(sig*100))
        ##parser.set('Input','diffusion_image',noisyfilename)
        #dsize, signal = read_flt_file(self.InputFile)
        #r.s = r.S0scale*signal/r.S0  # KxM
        # option 2: generate noisy data on the fly
        recon.s = add_rician_noise(s_noiseless, sig)
        recon.reconstruct()
        #angles[root_name][int(sig*100)] = fiber_orientation(r.coeff, init_pos)
        angles = fiber_orientation(recon.coeff, init_pos)
        #print angles
        deviation = [deviation_angles(angle, init_pos) for angle in angles]
        deviations[int(sig*100)] = deviation
    return deviations


def save_result(datafile):
    methods = ['dot','mow','qbi']
    dd = [expt_2fib(method) for method in methods]
    mean_2fib = numpy.array([numpy.array(md[i]).mean() for i in range(1,11) for md in dd]).reshape(10,3)
    std_2fib = numpy.array([numpy.array(md[i]).std() for i in range(1,11) for md in dd]).reshape(10,3)
    sigma = numpy.linspace(0.01,0.10,10)
    delta = 0.002
    sigma = numpy.linspace(0.01,0.10,10)
    data = numpy.c_[sigma-delta,sigma,sigma+delta,mean_2fib,std_2fib]
    numpy.savetxt(datafile, data)
    return data


def plot_errorbar(datafile, figname):
    # Initialize graph object
    g = graph.graphxy(width=8, ratio=4./3, key=graph.key.key(pos='tl'),
                      x=graph.axis.linear(min=0.005, max=0.105, title="std. dev. of noise"),
              y=graph.axis.linear(min=0, max=30, title="error angle (degree)"))

    methods = ["DOT","MOW","QBI"]
    errorbar_color = [color.rgb.red,color.rgb.black,color.rgb.blue]
    line_color = [color.rgb.red, color.rgb.black, color.rgb.blue]
    symbol_shape = [graph.style.symbol.square, graph.style.symbol.triangle, graph.style.symbol.circle]
    symbol_color = [color.rgb.red, color.rgb.black, color.rgb.blue]

    for i, method in enumerate(methods):
        # Plot the i-th line
        g.plot(graph.data.file(datafile, x=i+1, y=i+4, dy=i+7, title = method),
               styles=[graph.style.errorbar(errorbarattrs=[errorbar_color[i]]),
                       graph.style.line([line_color[i],
                                         style.linestyle.solid,
                                         style.linewidth.thick]),
                       graph.style.symbol(symbol_shape[i], symbolattrs=[symbol_color[i]])])

    # Now plot the text, horizontally centered
    g.text(g.width/2, g.height + 0.2, "2-fib",
           [text.halign.center, text.valign.bottom, text.size.Large])

    # Write the output
    g.writeEPSfile(figname)
    g.writePDFfile(figname)

