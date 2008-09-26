/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: special_functions.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:45 $
Version:   $Revision: 1.1 $
=========================================================================*/

#ifndef __SPECIAL_FUNCTIONS_H_
#define __SPECIAL_FUNCTIONS_H_

#define PI   3.141592653589f
/*
  associated Legendre polynomial
  P_l^m(cos(t)) 
  =(-1)^{l + m} /(2^l l!) sin^m(t)(d/d cos(t))^{l + m} sin^{2l}(t)
*/

double pLegendre(int l, int m, double t);

double logGamma(double xx);

double factorial(int n);

double erf(double x);
double erfc(double x);


#endif  //__SPECIAL_FUNCTIONS_H_