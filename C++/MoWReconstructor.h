/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: MoWReconstructor.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.9 $
=========================================================================*/

#ifndef _mow_reconstructor_h_
#define _mow_reconstructor_h_

#include "spherical_harmonics.h"
#include "FltFile.h"
#include "Reconstructor.h"

class MoWReconstructor: public Reconstructor
{
	

public:	
	void DoReconstruction();
	
	MoWReconstructor(): p(2)
	{
		ew[0] = 0.0015f;
		ew[1] = 0.0004f;
		ew[2] = 0.0004f;	
	}

private:
	int N;
	

	float ew[3];  		// eigenvalues;	
	vnl_matrix<float> ev; 			// eigenvectors;

	
	float cond;	
	float p;		// p-paremeter of Wishart distribution

	vnl_matrix <float> A, Ainv, R, R_QBI, R_DSI;
	
//private:
	void ComputeP(vnl_vector <float> &S, vnl_vector <float> &P);
	void AssembleA();
	void LoadWishartParameters();

};



#endif //_mow_reconstructor_h_



