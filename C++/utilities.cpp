/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: utilities.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:45 $
Version:   $Revision: 1.7 $
=========================================================================*/

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#define BUFSIZE 4096

#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_trace.h>
#include "utilities.h"
#include "nnls.h"

template <class T>
vnl_matrix_fixed < T, 3, 3> projection(vnl_vector_fixed < T, 3> v)
{                      
    // vnl_vector_fixed<float, 3> v (x);
    v  = v.normalize();
    vnl_matrix_fixed < T, 3, 3> P;
    P.set_identity();
    P = P - outer_product(v, v);
    return P;
}


int 
ReadValuesFromFile(const char* filename, std::vector < float > &vectorList)
{
	std::vector < std::string> file;
	std::string line;
	file.clear();
	std::ifstream infile(filename, std::ios_base::in);
	while (std::getline(infile, line, '\n'))
	{
    		  file.push_back(line);
	}
	std::cout << "Read " << file.size() << " lines.\n";
	char *pEnd;
	for (unsigned int i = 0; i < file.size(); i++)
	{
		char *pch = strtok((char *)file[i].c_str(), " ,");
		int k = 0;
		float v;
		while (pch != NULL)
		{
			// std::cout<<pch<<" ";
			if (k == 0)	
			{ v = static_cast<float>(strtod(pch, &pEnd)); }
			pch = strtok(NULL, " ,");
			k++;
		}
		// std::cout<< "read "<<k<< " tokens" << std::endl;
		vectorList.push_back(v);
	}
	return vectorList.size();
}


int 
ReadVectorsFromFile(const char* filename, std::vector < vnl_vector_fixed < float, 3> > &vectorList)
{
	std::vector < std::string> file;
	std::string line;
	file.clear();
	std::ifstream infile(filename, std::ios_base::in);
	while (std::getline(infile, line, '\n'))
	{
    		  file.push_back(line);
	}
	//!std::cout << "Read " << file.size() << " lines.\n";
	char *pEnd;
	for (unsigned int i = 0; i < file.size(); i++)
	{
		char *pch = strtok((char *)file[i].c_str(), " ,");
		int k = 0;
		vnl_vector_fixed <float, 3> v;
		while (pch != NULL)
		{
			// std::cout<<pch<<" ";
			if (k < 3)	
			{ v[k] = static_cast<float>(strtod(pch, &pEnd)); }
			pch = strtok(NULL, " ,");
			k++;
		}
		// std::cout<< "read "<<k<< " tokens" << std::endl;
		vectorList.push_back(v);
	}
	return vectorList.size();
}


int 
ReadVectorsFromFile(const char* filename, std::vector < vnl_vector_fixed < float, 6> > &vectorList)
{
	std::vector < std::string> file;
	std::string line;
	file.clear();
	std::ifstream infile(filename, std::ios_base::in);
	while (std::getline(infile, line, '\n'))
	{
    		  file.push_back(line);
	}
	//!std::cout << "Read " << file.size() << " lines.\n";
	char *pEnd;
	for (unsigned int i = 0; i < file.size(); i++)
	{
		char *pch = strtok((char *)file[i].c_str(), " ,");
		int k = 0;
		vnl_vector_fixed <float, 6> v;
		while (pch != NULL)
		{
			// std::cout<<pch<<" ";
			if (k < 6)	
			{ v[k] = static_cast<float>(strtod(pch, &pEnd)); }
			pch = strtok(NULL, " ,");
			k++;
		}
		// std::cout<< "read "<<k<< " tokens" << std::endl;
		vectorList.push_back(v);
	}
	return vectorList.size();
}

int vnl_NNLS(vnl_matrix<float>& A,  vnl_vector<float>& b, vnl_vector<float>& sol)
{
	int m = A.rows();
	int n = A.cols();
	// Using column major for results....
	double *a = new double[m*n];

	for (int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			a[i*m+j] = A(j,i);
		}
	}
 
	double *rhs = new double[m];
	for (int j=0;j<m;j++)
	{
		rhs[j] = b[j];
	}

	//float *x = new float[n];
	double rnorm;
	double *w = new double[n];
	double *zz = new double[m];
	int *indx = new int[3*n];
	int mode;
 
	double *x = new double[n];
	nnls(a, m, m, n, rhs, x, &rnorm, w, zz, indx, &mode);
	for (int i = 0; i < n; i++)	sol[i] = x[i];
	delete[] a;
	delete[] rhs;
	delete[] x;
	delete[] w;
	delete[] zz;
	delete[] indx;
	return mode;
	//return sol;
 
}


int get_config_fullpath(const char* input_config,char* f_config)
{
#ifdef WIN32
	char* lpPart[BUFSIZE]={NULL};
	int retval = GetFullPathName(input_config,
			     BUFSIZE, f_config,
			     lpPart);

	if (retval == 0) 
	{
		// Handle an error condition.
		printf ("GetFullPathName failed (%d)\n", GetLastError());
		return -1;
    }
	else {
		printf("The full path name is:  %s\n", f_config);

	}
#else
	strcpy(f_config,input_config);
#endif
	return 0;
}
