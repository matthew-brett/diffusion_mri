/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: QBIReconstructor.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.13 $
=========================================================================*/


#include "Reconstructor.h"
#include "QBIReconstructor.h"

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_trace.h>

void QBIReconstructor::ComputeP(vnl_vector <float> &S, vnl_vector <float> &P)
{
}

void QBIReconstructor::ComputeODF(unsigned i, vnl_vector <float> &S)
{
	//std::cout << "at [" << i << "]:" << S << std::endl;
	vnl_vector <float> SPH_coeff =(inv_SPH_basis_S) * S;
	vnl_vector <float> coeff_real, coeff_imag;
	coeff_real.set_size(nMoreCoeffs);
	coeff_imag.set_size(nMoreCoeffs);
	compute_qball_ODF(degree, SPH_coeff, coeff_real, coeff_imag);
	//convert_coeffs(degree, SPH_coeff, coeff_real, coeff_imag);
	for (int j = 0; j < nCoeffs; j++)
	{
		m_SPH_coeff[j][i] = SPH_coeff[j];
	}
	for (int j = 0; j < nMoreCoeffs; j++)
	{
		m_coeff_real[j][i] = coeff_real[j];
		m_coeff_imag[j][i] = coeff_imag[j];
	}
	
}

void QBIReconstructor::ComputeODF()
{
	std::cout << "start computing ODF ... " << std::endl;
	int nProcessed, nIgnored;
	nProcessed = 0;
	nIgnored = 0;
	
	vnl_vector <float> coeff;
	coeff.set_size(nCoeffs);
	
	for (int i = 0; i < nVoxels; i++)
	{	
		vnl_vector <float> S(nSamples);
		bool bProcess = ReadSignal(i, S);
		if (bProcess)
		{
			SignalAttenuation(i, S);
			ComputeODF(i,S);			
			nProcessed++;
		}
		else
		{
			//std::cout << "at pixel " << i <<" S.min = " << S.min_value() << std::endl;
			for (int j = 0; j < nCoeffs; j++)
			{	
				m_SPH_coeff[j][i] = 0;
			}
			for (int j = 0; j < nMoreCoeffs; j++)
			{
				m_coeff_real[j][i] = 0;
				m_coeff_imag[j][i] = 0;
			}
			nIgnored ++;			
		}
	}
	std::cout << "finish computing ODF ... " << std::endl;
	std::cout << "processed " << nProcessed << " ignored " << nIgnored << std::endl;
}

void
QBIReconstructor::
AssembleSPHBasisFromGradients()
{
	assert(gradients.rows()>0);
	///!!! note it requires gradient directions 
	degree = GetPrivateProfileInt("QBI", "degree", 6, f_config);
	if (degree < 4)
		degree = 6;
	std::cout << "The degree of Spherical Harmonics for QBI is " << degree << std::endl;

	SPH_basis_S = construct_basis(degree, gradients);
	//vcl_cout << SPH_basis_S << vcl_endl;
	inv_SPH_basis_S = vnl_svd < float>(SPH_basis_S).inverse();
	std::cout << "size of SPH_basis_S: "
		<< SPH_basis_S.rows() << "x" 
		<< SPH_basis_S.cols() << std::endl;
#if 0   // in case regulariation is needed
	vnl_svd < float> mySVD(SPH_basis_S);
	// The SVD class holds three matrices U, W, V such that the original 
	// matrix A can be written A = U*W*V^T.  
	//
	// The matrices U and V are orthogonal (i.e. U*U^T = U^T*U = I), and 
	// the matrix W stores the singular values in decreasing order along 
	// its diagonal.
	vnl_matrix < float>      myU = mySVD.U();
	vnl_diag_matrix < float> myW = mySVD.W();  
	vnl_matrix < float>      myV = mySVD.V();
	vnl_diag_matrix < float> myWinv = mySVD.Winverse();
	
	std::cout << "sphbasis: sigma_max = " << mySVD.sigma_max() << " sigma_min = " << mySVD.sigma_min() << std::endl;  
	
	int rnk = mySVD.rank();
	vnl_matrix < float> Winverse(myWinv.rows(), myWinv.columns());
	Winverse.fill(0);
	//float lambda = cond; // mySVD.sigma_max() / cond;
        float lambda = _lambda;
        // see Maxime Descoteaux TR5681
        vnl_vector<float> Ldiag;
        int k =(degree + 2)*(degree + 1)/2;
        Ldiag.set_size(k);
     	for (int l = 0; l <= degree; l += 2)
    	{
		int center =(l + 2)*(l + 1)/2 - l -1;
		Ldiag[center] = l;
		for (int m = 1; m <= l; m++)
		{
			Ldiag[center + m] = l;
			Ldiag[center - m] = l;
		}  
    	}
	for (unsigned int i = 0; i < rnk;++i)
	{
		Winverse(i, i)= myW(i,i) / (myW(i, i)*myW(i,i) + lambda*Ldiag[i]);
	}
	inv_SPH_basis_S = myV * Winverse * myU.conjugate_transpose();
	//inv_SPH_basis_S = vnl_svd < float>(SPH_basis_S).inverse();
#endif
}


void QBIReconstructor::
DoReconstruction() 
{
	// step1: prepare data
	std::cout<< "start QBI_ODF reconstruction... "<<std::endl;
	LoadInputDataFiles();
	PrepareInputData();
	// step2: 
	LoadImagingParams();
	LoadCommonConstants();
	// step3:
	//LoadWishartParameters();
	//step 4:
	//LoadOutputDataFiles();
	std::cout<< "Data and parameters loaded."<<std::endl;
	AssembleSPHBasisFromGradients();
	//step 5: prepare linear system
	//step6: reconstruct
	PrepareOutputData();
	ComputeODF();
	std::cout<< "QBI ODF computed."<<std::endl;
	Write();	
	std::cout<< "QBI ODF saved."<<std::endl;
}



