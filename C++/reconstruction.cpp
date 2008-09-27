/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: test_mow.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:45 $
Version:   $Revision: 1.7 $
=========================================================================*/

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
using namespace std;

#include <vnl/vnl_sym_matrix.h>
#include <vcl_iostream.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_cross.h>

#include "MoWReconstructor.h"
#include "DOTReconstructor.h"
#include "QBIReconstructor.h"


typedef auto_ptr<MoWReconstructor> MoWReconstructorPtr;
typedef auto_ptr<DOTReconstructor> DOTReconstructorPtr;
typedef auto_ptr<QBIReconstructor> QBIReconstructorPtr;


int main(int argc, char* argv[])
{
    if (argc<3)
    {
        cerr << "Usage: reconstruction config_file method" << std::endl;
        return -1;
    }

    if (!strcmp(argv[2], "mow"))
    {
        MoWReconstructorPtr reconstructor(new MoWReconstructor);
        reconstructor->SetConfigFile(argv[1]);
        reconstructor->DoReconstruction();
    }
    else if (!strcmp(argv[2], "dot"))
    {
        DOTReconstructorPtr reconstructor(new DOTReconstructor);
        reconstructor->SetConfigFile(argv[1]);
        reconstructor->DoReconstruction();
    }
    else if (!strcmp(argv[2], "qbi"))
    {
        QBIReconstructorPtr reconstructor(new QBIReconstructor);
        reconstructor->SetConfigFile(argv[1]);
        reconstructor->DoReconstruction();
    }
    else
    {
        cerr << "Currently only three methods are supported. Try 'mow','dot' or 'qbi'. " << std::endl;
        return -1;
    }
    return 0;
}
