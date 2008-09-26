/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: spherical_harmonics.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.17 $
=========================================================================*/



#include "spherical_harmonics.h"
#include "special_functions.h"



void
convert_coeffs(int degree, vnl_vector<float> sym_coeff, vnl_vector<float> & coeff_real, vnl_vector<float> & coeff_imag)
{

	int i = 0;
	for (int l=0; l<=degree; l+=2)
	{
		int center =(l + 2)*(l + 1)/2 - l -1;
		coeff_real[i] = (sym_coeff[center]); // m = 0
		coeff_imag[i] = (0); 
		++i;
		for (int m=1; m<=l; ++m)
		{
			if (m%2)
			{
				coeff_imag[i] = sym_coeff[center-m]*sqrt(2.0f)/2;
				coeff_real[i] = -sym_coeff[center+m]*sqrt(2.0f)/2; 
			}
			else
			{
				coeff_imag[i] = -sym_coeff[center-m]*sqrt(2.0f)/2;
				coeff_real[i] = sym_coeff[center+m]*sqrt(2.0f)/2;
			}
			++i;
		}
	}
}

void
convert_coeffs_back(int degree, vnl_vector<float>  coeff_real, vnl_vector<float>  coeff_imag, vnl_vector<float> & sym_coeff)
{
	int i = 0;
	for (int l=0; l<=degree; l+=2)
	{
		int center =(l + 2)*(l + 1)/2 - l -1;
		sym_coeff[center] = coeff_real[i]; // m = 0
		++i;
		for (int m=1; m<=l; m++)
		{
			if (m%2)
			{
				//coeff_imag[i] = sym_coeff[center-m]*sqrt(2)/2;
				sym_coeff[center-m] = coeff_imag[i]*sqrt(2.0f);
				//coeff_real[i++] = -sym_coeff[center+m]*sqrt(2)/2;
				sym_coeff[center+m] = -coeff_real[i]*sqrt(2.0f);
			}
			else
			{
				//coeff_imag[i] = -sym_coeff[center-m]*sqrt(2)/2;
				sym_coeff[center-m] = -coeff_imag[i]*sqrt(2.0f);
				//coeff_real[i++] = sym_coeff[center+m]*sqrt(2)/2;
				sym_coeff[center+m] = coeff_real[i]*sqrt(2.0f);
			}
			++i;
		}
	}
}

void
compute_qball_ODF(int degree, vnl_vector<float> & sym_coeff, vnl_vector<float> & coeff_real, vnl_vector<float> & coeff_imag)
{
    convert_coeffs(degree, sym_coeff, coeff_real, coeff_imag);
	int i = 0;
	for (int l=0; l<=degree; l+=2)
	{

		int start = (l/2)*(l/2); 
		int oddl = 1;
		int evenl = 1;
		for (int i=1;i<l;i+=2)
		{
			oddl *= i;
			evenl *= (i+1);
		}
		float ll = oddl/(evenl*1.0f);
		if ((l/2)%2) ll = -ll;
		for (int m=0; m<=l; ++m)
		{
        	   coeff_imag[start+m] *= ll;
               coeff_real[start+m] *= ll;   	  
		}
	}
	convert_coeffs_back(degree, coeff_real, coeff_imag, sym_coeff);
	
}


vnl_matrix <float> 
construct_basis(int degree, vnl_matrix<float> sphere_points)
{
	assert(sphere_points.rows()>0);
	assert(sphere_points.cols()==3);
	int n = sphere_points.rows();
	int k =(degree + 2)*(degree + 1)/2;
	vnl_matrix <float> basis(n, k);
	// basis.set_size(n,k);
	for (int i = 0; i < n; ++i)
	{
		vnl_vector_fixed<float,3> v = sphere_points.get_row(i).normalize();
		float theta = acos(v[2]); //polar
		float phi = atan2(v[1], v[0]);  //azimuth  
		for (int l = 0; l <= degree; l += 2)
		{
			int center =(l + 2)*(l + 1)/2 - l -1;
			float lconstant = sqrt((2*l + 1)/(4*PI));
		
			basis(i, center) = lconstant*pLegendre(l, 0, theta);
			float precoeff;
			for (int m = 1; m <= l; ++m)
			{
				precoeff = lconstant * sqrt(2.0f)*sqrt(factorial(l - m)/factorial(l + m));
				if (m%2) { precoeff = -precoeff;}
				precoeff *= 1.0f*pLegendre(l, m, theta);
				basis(i, center + m) = precoeff*cos(m*phi);
				basis(i, center - m) = precoeff*sin(m*phi);
			}  
		}
	}	
	return basis;		
}

void
comp_magnitude(int degree, vnl_vector<float> sym_coeff, vnl_vector<float> & coeff_mag)
{
	coeff_mag[0] = sym_coeff[0];
	float magnitude = 0;
	for (int l=2; l<=degree; l+=2)
	{
		int start_pos = l*(l-1)/2;
		for (int k=0;k<2*l+1;++k)
		{
			magnitude += sym_coeff[start_pos+k]*sym_coeff[start_pos+k];
		}
		coeff_mag[l/2] = sqrt(magnitude);
	}
}


