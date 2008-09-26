/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: ODFReconstructor.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.7 $
=========================================================================*/

#ifndef _qbi_reconstructor_h_
#define _qbi_reconstructor_h_

#include "spherical_harmonics.h"
#include "FltFile.h"
#include "Reconstructor.h"

class QBIReconstructor: public Reconstructor
{
	
public:	
	void DoReconstruction();

private:
	vnl_matrix < float> SPH_basis_S, inv_SPH_basis_S;
	
	void ComputeP(vnl_vector <float> &S, vnl_vector <float> &P);
	void AssembleSPHBasisFromGradients();	
//private:
	void ComputeODF(unsigned i, vnl_vector <float> &S);
	void ComputeODF();

};



#endif //_qbi_reconstructor_h_



