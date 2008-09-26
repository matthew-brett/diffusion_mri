/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: test_odf.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:45 $
Version:   $Revision: 1.8 $
=========================================================================*/


#include <vnl/vnl_sym_matrix.h>
#include <vcl_iostream.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_cross.h>


#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
using namespace std;

#include "QBIReconstructor.h"

typedef auto_ptr<QBIReconstructor> QBIReconstructorPtr;

int main(int argc, char* argv[])
{
	
	if (argc<2)
	{
		cerr<< "Usage: Reconstruct ConfigFile" << std::endl;
		return -1;
	}
	
	QBIReconstructorPtr reconstructor(new QBIReconstructor);
	reconstructor->SetConfigFile(argv[1]); 
	reconstructor->DoReconstruction();

	return 0;	
}
