/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: DOTReconstructor.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.10 $
=========================================================================*/

#ifndef _dot_reconstructor_h_
#define _dot_reconstructor_h_

#include "spherical_harmonics.h"
#include "FltFile.h"
#include "Reconstructor.h"

class DOTReconstructor: public Reconstructor
{
	
public:	
	void DoReconstruction();

private:
	void ComputeP(vnl_vector <float> &S, vnl_vector <float> &P);
	void AssembleM();
	void LoadDOTConstants();

private:
	vnl_matrix <float> M;
	int degreeDOT;
	vnl_vector <float> wj;

};



#endif //_dot_reconstructor_h_




