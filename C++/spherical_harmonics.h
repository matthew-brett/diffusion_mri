/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: spherical_harmonics.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:45 $
Version:   $Revision: 1.9 $
=========================================================================*/

#ifndef __SPHERICAL_HARMONICS_H_
#define __SPHERICAL_HARMONICS_H_

#include <vcl_iostream.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_matrix_fixed.h>


#include <cmath>
//#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstring>


void
convert_coeffs(int degree, vnl_vector<float> sym_coeff, vnl_vector<float> & coeff_real, vnl_vector<float> & coeff_imag);

void
compute_qball_ODF(int degree, vnl_vector<float> & sym_coeff, vnl_vector<float> & coeff_real, vnl_vector<float> & coeff_imag);


vnl_matrix <float> 
construct_basis(int degree, vnl_matrix<float> sphere_points);

void
comp_magnitude(int degree, vnl_vector<float> sym_coeff, vnl_vector<float> & coeff_mag);



#endif  //__SPHERICAL_HARMONICS_H_
