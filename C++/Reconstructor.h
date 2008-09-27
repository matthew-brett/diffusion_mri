/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: Reconstructor.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.42 $
=========================================================================*/

#ifndef _reconstructor_h_
#define _reconstructor_h_

#include <vcl_iostream.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_sym_matrix.h>

#include "spherical_harmonics.h"
#include "FltFile.h"
#include "utilities.h"

// Modified by Guang Cheng 10/27/2007
#define BUFF_MAX 255
// End

typedef vnl_vector_fixed <float, 3>  vector_3;
typedef vnl_vector_fixed <float, 6>  vector_6;

/*! \brief abstract base class for reconstruction task
*   \todo declare a copy ctor and an assignment operator
*/
class Reconstructor
{

public:
    /*! \brief default ctor with initialization list, see explanation of data members.
    *
    */
    Reconstructor():
        fltFile(0), S0Image(0), b2Dslice(0),
        minSig(0.001f), S0Scale(1), S0default(1), bUseS0(0),
        bUseBMatrix(0), bFixedBValue(1), degree(6)
    {
        std::cout << "ctor of Reconstructor is being called ..." << std::endl;
    }

    /*! \brief dtor, make it virtual!
    *
    */
    virtual ~Reconstructor();

    /*! \brief set configuration file
    *
    */
    inline void SetConfigFile(const char* input_config)
    {
        get_config_fullpath(input_config, this->f_config);
    };

    /*! \brief interface for performing reconstruction, see implementation in derived classes
    *
    */
    virtual void DoReconstruction() = 0;

protected:
    virtual void ComputeP(vnl_vector <float> &S, vnl_vector <float> &P) = 0;
    bool ReadSignal(unsigned i, vnl_vector <float> &S);
    void SignalAttenuation(unsigned i, vnl_vector <float> &S);
    void ComputeADC(vnl_vector <float> &S);
    void SaveProbabilityProfile(unsigned i, vnl_vector <float> &P);
    void Reconstruct();
    void Write();
    void Prepare();
    void LoadOutputDataFiles();
    void LoadInputDataFiles();
    void PrepareInputData();
    void PrepareOutputData();
    void LoadImagingParams();
    void LoadCommonConstants();

    void LoadFltFile(const char*);
    void LoadS0File(const char*);
    void WriteFltFile(const char*);

    void AssembleSPHBasis();


protected:
    /*! \brief
    *  diffusion time, measured in seconds; typical value: 0.020s
    */
    float t;
    /*! \brief
    *  radius of the probability surface, measured in mm = 1000um; typical value: 0.015mm
    */
    float r;
    /*! \brief
    *  regularization factor for spherical harmonics coefficients
    */
    float lambda;
    /*! \brief
    *  configuration file
    */
    char f_config[BUFF_MAX];
    // int K, T;
    /*! \brief
    *  K*3 matrix containing gradient directions where K is the # of measurements
    */
    vnl_matrix <float> gradients;
    /*! \brief
    *  K*6 matrix containing b-matrices where K is the # of measurements
    */
    vnl_matrix <float> bmatrix;
    /*! \brief
    *  T*6 matrix containing tessellation vertices where T is the # of vertices
    */
    vnl_matrix <float> tessellation;
    /*! \brief
    *  b-values of diffusion gradients
    */
    vnl_vector <float> b;

    FltFile*  fltFile;
    int nSamples,nVoxels;
    int bUseBMatrix;

    int degree;
    int nCoeffs, nMoreCoeffs;
    std::vector<float *> m_SPH_coeff;
    //std::vector<float *> m_SPH_coeff_mag;
    std::vector<float *> m_coeff_real;
    std::vector<float *> m_coeff_imag;
    std::vector<float *> m_weights;

    std::vector<float *> m_S_coeff;
    std::vector<float *> m_S_coeff_real;
    std::vector<float *> m_S_coeff_imag;
    std::vector<float *> m_qball_coeff_real;
    std::vector<float *> m_qball_coeff_imag;

private:
    float minSig;
    std::vector <float*> data, S0;

    vnl_matrix < float> SPH_basis, inv_SPH_basis;
    FltFile*  S0Image;
    float S0default, S0Scale;
    bool b2Dslice;
    int bUseS0;
    int bFixedBValue;
};



#endif //_reconstructor_h_







