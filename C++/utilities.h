/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: utilities.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:45 $
Version:   $Revision: 1.4 $
=========================================================================*/

#ifndef __AUX_UTILS_H_
#define __AUX_UTILS_H_

/*! \file utilities.h
    \brief Header file including declaration of some auxiliary functions
*/



#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_3x3.h>

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <iomanip> 


#define PI  3.141592653589f
template <class T>
vnl_matrix_fixed < T, 3, 3> projection(vnl_vector_fixed < T, 3> v);

int 
ReadValuesFromFile(const char* filename, std::vector < float > &vectorList);

int 
ReadVectorsFromFile(const char* filename, std::vector < vnl_vector_fixed < float, 3> > &vectorList);

int 
ReadVectorsFromFile(const char* filename, std::vector < vnl_vector_fixed < float, 6> > &vectorList);

int vnl_NNLS(vnl_matrix<float>& A,  vnl_vector<float>& b, vnl_vector<float>& x);


int get_config_fullpath(const char* input_config,char* f_config);


#endif // __AUX_UTILS_H_
