/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: Reconstructor.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.49 $
=========================================================================*/


#include "Reconstructor.h"

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_trace.h>
#include "utilities.h"


/*! \brief Read signal values at #i voxel
 *
 */
bool Reconstructor::ReadSignal(unsigned i, vnl_vector <float> &S)
{
    if (b2Dslice)
    {
        for (int j = 0; j < nSamples; ++j)
        {
            S[j] = data[0][j*nVoxels + i];
            //if (S[j]<=0) S[j] = 1.0e-3;
            // or do try catch later
        }
    }
    else
    {
        for (int j = 0; j < nSamples; ++j)
        {
            S[j] = data[j][i];
            //if (S[j]<=0) S[j] = 1.0e-3;
        }
    }
    return (S.min_value() >= minSig) ? true : false; 
}


/*! \brief Compute signal attenuation S/S_0
 *         
 */
void Reconstructor::SignalAttenuation(unsigned i, vnl_vector <float> &S)
{
    if (bUseS0 && (S0Image != NULL) )
    {
        for (int j = 0; j < nSamples; ++j)
        {
            S[j] = S[j]/(S0[0][i]*S0Scale); 
        }
    }
    else
    {
        for (int j = 0; j < nSamples; ++j)
        {
            S[j] = S[j]/S0default; 
        }
    }
}

/*! \brief Compute apparent diffusion coefficients using
 *              S = S0*exp(-b*ADC)
 */
void Reconstructor::ComputeADC(vnl_vector <float> &S)
{
    for (int j = 0; j < nSamples; ++j)
    {
        S[j] = -log(S[j])/b[j]; 
    }
}

/*! \brief Represent the computed P(r=R0, theta, phi) 
 *         using spherical harmonic coefficients
 */
void Reconstructor::SaveProbabilityProfile(unsigned i, vnl_vector <float> &P)
{
 
    vnl_vector < float> SPH_coeff =(inv_SPH_basis)*P;
    for (int j = 0; j < nCoeffs; ++j)
    {
        m_SPH_coeff[j][i] = SPH_coeff[j];
    }
    
    vnl_vector < float> coeff_real, coeff_imag;
    coeff_real.set_size(nMoreCoeffs);
    coeff_imag.set_size(nMoreCoeffs);
    convert_coeffs(degree, SPH_coeff, coeff_real, coeff_imag);
    for (int j = 0; j < nMoreCoeffs; ++j)
    {
        m_coeff_real[j][i] = coeff_real[j];
        m_coeff_imag[j][i] = coeff_imag[j];
    }
}

/*! \brief Compute P voxel by voxel 
 *
 */
void Reconstructor::Reconstruct()
{
    int n = tessellation.rows();
    assert(n>0);
    std::cout << "start reconstruction ... " << std::endl;
    int nProcessed, nIgnored;
    nProcessed = 0;
    nIgnored = 0;
    vnl_vector < float> P;
    for (int i = 0; i < nVoxels; ++i)
    {   
        vnl_vector <float> S(nSamples);
        bool bProcess = ReadSignal(i, S);
        if (bProcess)
        {
            SignalAttenuation(i, S);
            ComputeP(S, P);
            ++nProcessed;
        }
        else
        {
            //std::cout << "at pixel " << i <<" S.min = " << S.min_value() << std::endl;
            P.set_size(n);
            P.fill(0);
            ++nIgnored;         
        }
        SaveProbabilityProfile(i, P);
    }
    std::cout << "finish reconstruction ... " << std::endl;
    std::cout << "processed " << nProcessed << " ignored " << nIgnored << std::endl;
}


/*! \brief Write the spherical harmonics coefficents into files 
 *
 */
void Reconstructor::Write()
{
    char f_output[BUFF_MAX];
    GetPrivateProfileString("Output", "output_basename", NULL, f_output, BUFF_MAX, f_config);
    if (strlen(f_output)>0 && m_SPH_coeff.size()>0)
    { 
        WriteFltFile(f_output);
        std::cout << "Write output to " << f_output << "_(real|imag).flt" <<std::endl;
    }
}



/*! \brief Read imaging parameters from configuration file
 *
 */
void Reconstructor::LoadImagingParams()
{
    // Part I: get bmatrix information or diffusion gradient directions and bvalue 
    char f_bmatrix[BUFF_MAX] = "\0";
    GetPrivateProfileString("Input", "b_matrices", NULL, f_bmatrix, BUFF_MAX, f_config);
    if (strlen(f_bmatrix) == 0)
    {
        std::cout << "No b_matrices file is provided." <<std::endl;
        char f_gradient[BUFF_MAX] = "\0";
        GetPrivateProfileString("Input", "diffusion_gradients", NULL, f_gradient, BUFF_MAX, f_config);
        if (strlen(f_gradient) == 0)
        {
            throw "No diffusion gradients is provided.";
            return;
        }
        else  
        {
            std::cout << "Reading gradient directions from " << f_gradient << " ... " << std::endl;
            //K = ReadVectorsFromFile(f_gradient, gradients);
            std::ifstream infile(f_gradient);
            gradients.read_ascii(infile);
            assert(gradients.rows()>0);
            assert(gradients.cols()==3);
            //K = gradients.rows();         
            std::cout <<  gradients.rows() << " gradient directions were loaded from " << f_gradient << std::endl;
        }

        char f_bvalues[BUFF_MAX] = "\0"; // measured in s/mm^2
        GetPrivateProfileString("Input", "b_value_infile", NULL, f_bvalues, BUFF_MAX, f_config);
        if (strlen(f_bvalues) == 0)
        {
            std::cout << "No file for b-values is provided. Use single fixed b-value instead." <<std::endl;
            char s_bvalue[BUFF_MAX] ="\0";  // measured in s/mm^2
            GetPrivateProfileString("Input", "b_value_fixed", "1500", s_bvalue, 10, f_config);
            float bv = static_cast<float>(atof(s_bvalue));
            int n = gradients.rows();
            b.set_size(n);
            for (int i = 0; i < n; ++i)
            {   
                b[i] = bv; //push_back(bv);
            }
            std::cout << "fixed b-value: " << b[n - 1] << std::endl;
            return;
        }
        else
        {
            std::cout << "Reading b-values from " << f_bvalues << " ... " << std::endl;
            //K = ReadValuesFromFile(f_bvalues, b);
            std::ifstream infile(f_bvalues);
            b.read_ascii(infile);
            assert(b.size()>0);         
            std::cout <<  b.size() << " b-values were loaded from " << f_bvalues << std::endl;
        }

    }
    else  
    {
        std::cout << "Reading B matrices from " << f_bmatrix << " ... " << std::endl;
        //K = ReadVectorsFromFile(f_bmatrix, bmatrix);
        std::ifstream infile(f_bmatrix);
        bmatrix.read_ascii(infile);
        assert(bmatrix.rows()>0);
        assert(bmatrix.cols()==6);
        //K = bmatrix.rows();           
        std::cout <<  bmatrix.rows() << " B matrices were loaded from " << f_bmatrix << std::endl;
    }   
}

/*! \brief  Get the path to input data by parsing the configuration file
 *
 */
void Reconstructor::
LoadInputDataFiles() 
{
    // Part IV: Input
    char f_input[BUFF_MAX] = "\0";
    GetPrivateProfileString("Input", "diffusion_image", NULL, f_input, BUFF_MAX, f_config);
    if (strlen(f_input) == 0)
    {
        //std::cerr << "No Input Image is provided." <<std::endl;
        throw "No input data file is provided!";
        return;
    }
    else
    {
        //!std::cout << strlen(f_input) << std::endl;
        LoadFltFile(f_input);
        std::cout << "Read input data from " << f_input << std::endl;
    }
    char f_S0Image[BUFF_MAX] = "\0";
    GetPrivateProfileString("Input", "baseline_image", NULL, f_S0Image, BUFF_MAX, f_config);
    if (strlen(f_S0Image) == 0)
    {
        std::cerr << "No baseline image is provided." <<std::endl;
        char s_S0value[BUFF_MAX] ="\0";
        GetPrivateProfileString("Input", "baseline_value", "1", s_S0value, 10, f_config);
        S0default = static_cast<float>(atof(s_S0value));
        std::cout << "Use default baseline value as S0: " << S0default << std::endl;

    }
    else
    {
        //!std::cout << strlen(f_input) << std::endl;
        LoadS0File(f_S0Image);
        std::cout << "Read S0Image Data from " << f_S0Image << std::endl;
    }
    char s_S0scale[BUFF_MAX] ="\0";
    GetPrivateProfileString("Input", "baseline_scale", "1", s_S0scale, 10, f_config);
    /* std::cout <<"scale: " << s_S0scale << std::endl; */
    S0Scale = static_cast<float>(atof(s_S0scale));
    std::cout << "Scaling S0 image by factor: " << S0Scale << std::endl;
}

/*! \brief  Get the path to output files by parsing the configuration file
 *
 */
void Reconstructor::
LoadOutputDataFiles() // const char* f_config
{
    //(3) tessellation for computing spherical_harmonics: degree, 
    char f_tessellation[BUFF_MAX];
    GetPrivateProfileString("Input", "tessellation", NULL, f_tessellation, BUFF_MAX, f_config);

    //T = ReadVectorsFromFile(f_tessellation, tessellation);    
    std::ifstream infile(f_tessellation);
    tessellation.read_ascii(infile);
    assert(tessellation.rows()>0);
    assert(tessellation.cols()==3);
    //T = tessellation.rows();          

    std::cout << "Read " << tessellation.rows() << " tessellation points from " << f_tessellation << std::endl;

/*  degree = GetPrivateProfileInt("Constants", "degree", 6, f_config);
    if (degree < 4)
        degree = 6;
    std::cout << "The degree of Spherical Harmonics for P is " << degree << std::endl;
*/
}


/*! \brief  Read some common parameters from configuration file
 *
 */
void Reconstructor::
LoadCommonConstants() // const char* f_config
{

    char s_minsig[BUFF_MAX] = "0.001";
    GetPrivateProfileString("Constants", "signal_threshold", NULL, s_minsig, 10, f_config);
    minSig = static_cast<float>(atof(s_minsig));
    std::cout << "minSig: " << minSig << std::endl;
    
    char s_rvalue[BUFF_MAX] = "0.015";  // measured in mm
    // millimeter (mm) = 10-3m;  micrometer (um = micron) = 10^-6m
    GetPrivateProfileString("Constants", "r", NULL, s_rvalue, 10, f_config);
    r = static_cast<float>(atof(s_rvalue));
    std::cout << "r: " << r << std::endl;
    
    char s_tvalue[BUFF_MAX] = "0.020";   // measured in seconds
    GetPrivateProfileString("Constants", "t", NULL, s_tvalue, 10, f_config);
    t = static_cast<float>(atof(s_tvalue));
    std::cout << "t: " << t << std::endl;

    char s_lambda[BUFF_MAX] = "0.006";   
    GetPrivateProfileString("Constants", "lambda", NULL, s_lambda, 10, f_config);
    lambda = static_cast<float>(atof(s_lambda));
    std::cout << "lambda: " << lambda << std::endl;
    
}


/*! \brief  Assemble the spherical harmonics basis matrix as well as its pseudoinverse
 *
 */
void
Reconstructor::
AssembleSPHBasis()
{
    //assert(tessellation.rows()>0);
    // call function defined in spherical_harmonics.h
    //SPH_basis = construct_basis(degree, tessellation);    


    char f_sh_basis[BUFF_MAX] = "\0";
    GetPrivateProfileString("Input", "sh_basis_matrix", NULL, f_sh_basis, BUFF_MAX, f_config);
    if (strlen(f_sh_basis) == 0)
    {
            throw "No spherical harmonics basis matrix is provided.";
            return;
    }
    else  
    {
        std::ifstream infile(f_sh_basis);
        SPH_basis.read_ascii(infile);
        std::cout <<  SPH_basis.rows() << "," << SPH_basis.cols() << " SPH basis loaded from " << f_sh_basis << std::endl;
    }
    degree = (sqrt((float)(8*SPH_basis.cols()+1))-3.0)/2;
    std::cout << "The degree of Spherical Harmonics for P is " << degree << std::endl;

    //inv_SPH_basis = vnl_svd < float>(SPH_basis).inverse();
    vnl_svd < float> mySVD(SPH_basis);
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

    if (lambda>0) // with additional regularization, see Maxime Descoteaux TR5681
    {
        vnl_vector<float> Ldiag;
        int k =(degree + 2)*(degree + 1)/2;
        Ldiag.set_size(k);
        for (int l = 0; l <= degree; l += 2)
        {
            int center = (l + 2)*(l + 1)/2 - l -1;
            Ldiag[center] = l*1.0;
            for (int m = 1; m <= l; ++m)
            {
                Ldiag[center + m] = l*1.0;
                Ldiag[center - m] = l*1.0;
            }  
        }
        // check this section
        for (int i = 0; i < rnk;++i)
        {
            Winverse(i, i)= myW(i,i) / (myW(i, i)*myW(i,i) + lambda*Ldiag[i]);
        }
    }
    else
    {
        for (int i = 0; i < rnk;++i)
        {
            Winverse(i, i)= 1 / myW(i, i);
        }

    }
    inv_SPH_basis = myV * Winverse * myU.conjugate_transpose();
}


/*! \brief Write the spherical harmonics coefficents into files 
 *
 */
void Reconstructor::
WriteFltFile(const char* filename)
{
    std::vector < unsigned> dsize = fltFile->getSize();
    
    int nCoeffs =(degree + 2)*(degree + 1)/2;
    int ndim = dsize.size();
    
    
    dsize[ndim - 1] = nCoeffs;
    std::cout << "size: [";
    for (int i = 0; i < ndim; ++i)
    {   std::cout << dsize[i] << " ";}
    std::cout << "]" <<std::endl;
    
    FltFile newFltFile;
    
    char fltname[BUFF_MAX];
    char * str1 = ".flt";
    char * str2 = "_real.flt";
    char * str3 = "_imag.flt";
    char * str4 = "_weight.flt";
    char mag_str[BUFF_MAX];
    
    strcpy(fltname, filename);
    strcat(fltname, str1);
    //!std::cout << fltname << std::endl;
    newFltFile.WriteFile(fltname, m_SPH_coeff, dsize, 'b', 1);
    

    nCoeffs =(degree/2 + 1)*(degree/2 + 1);
    dsize[ndim - 1] = nCoeffs;
    std::cout << "real|imag size: [";
    for (int i = 0; i < ndim; ++i)
    {   std::cout << dsize[i] << " ";}
    std::cout << "]" <<std::endl;
    
    strcpy(fltname, filename);
    strcat(fltname, str2);
    //!std::cout << fltname << std::endl;
    newFltFile.WriteFile(fltname, m_coeff_real, dsize, 'b', 1);
    
    strcpy(fltname, filename);
    strcat(fltname, str3);
    //!std::cout << fltname << std::endl;
    newFltFile.WriteFile(fltname, m_coeff_imag, dsize, 'b', 1);

    /*
    float * coeff_mag = new float[nVoxels];
    std::vector<float *> m_SPH_coeff_mag;
    m_SPH_coeff_mag.push_back(coeff_mag);
    std::cout << " # of mag =" << m_SPH_coeff_mag.size() <<  std::endl;
    dsize[ndim - 1] = 1;
    int pos = 0;
    
    for (int i=0;i<=degree;i+=2)
    {
        sprintf(mag_str, "_mag_%d", i);
        strcpy(fltname, filename);
        strcat(fltname, mag_str);
        strcat(fltname, str1);
        
    
        for (long j=0;j<nVoxels;++j)
        {
            coeff_mag[j] = 0;
            for (int k=0; k<2*i+1; ++k)
            {
                //coeff_mag[j] += (m_SPH_coeff[pos+k])[j]*(m_SPH_coeff[pos+k])[j];
            }
        }
        pos = ((i+2)*(i+1))/2;
        newFltFile.WriteFile(fltname, m_SPH_coeff_mag, dsize, 'b', 1);
      }
    */
    
}

/*! \brief Load attenuated signals saved in .flt file format
 *
 */
void Reconstructor::
LoadFltFile(const char* filename)
{
    this->fltFile = new FltFile(filename, 0);
}

/*! \brief Load baseline signals saved in .flt file format
 *
 */
void Reconstructor::
LoadS0File(const char* filename)
{
    this->S0Image = new FltFile(filename, 0);
    // Note S0 Image must be stored in 3D or 4D format with last dim = 1;
    // that is XxYx1 for 2D slice or XxYxZx1 for 3D volume
}

/*! \brief Allocate memory for input data
 *
 */
void Reconstructor::PrepareInputData()
{
    // std::cout << "Reconstruction ... " << std::endl;
    if (fltFile == NULL)
    {
        throw "Not Ready for Reconstruction. (Failed to read input data)";
        return;
        // todo: check other conditions
    }
    nSamples = fltFile->samples();
    nVoxels = fltFile->volume(); // / nSamples ;
    // bool b2Dslice = 0;
    if (nSamples == 1)  // must be 2D slice
    {
        nSamples =(fltFile->getSize())[2]; 
        nVoxels = nVoxels/nSamples;
        b2Dslice = 1;
    }
    data = fltFile->getData();
    if (bUseS0 &&(S0Image != NULL))
    {
        S0 = S0Image->getData();
    }
}

/*! \brief Allocate memory for output 
 *
 */
void Reconstructor::PrepareOutputData()
{
    // std::cout << "Reconstruction ... " << std::endl;
    nCoeffs =(degree + 2)*(degree + 1)/2;
    nMoreCoeffs =(degree/2 + 1)*(degree/2 + 1);
    std::cout << "nMoreCoeffs =" << nMoreCoeffs <<  std::endl;
    
    for (int j = 0; j < nCoeffs; ++j)
    {
        float * coeff = new float[nVoxels];
        m_SPH_coeff.push_back(coeff);
    }
    for (int j = 0; j < nMoreCoeffs; ++j)
    {
        float * coeff_real = new float[nVoxels];
        float * coeff_imag = new float[nVoxels];
        m_coeff_real.push_back(coeff_real);
        m_coeff_imag.push_back(coeff_imag);
    }
}

/*! \brief dtor: release memory
 *
 */
Reconstructor::~Reconstructor() 
{
    delete fltFile;
    delete S0Image;

    assert(m_coeff_real.size()==m_coeff_imag.size());
    for (unsigned int i=0; i<m_SPH_coeff.size(); ++i){
        delete[] m_SPH_coeff[i];
    }

    
    for (unsigned int i=0; i<m_coeff_real.size(); ++i){
        delete[] m_coeff_real[i];
        delete[] m_coeff_imag[i];
    }

    std::cout << "dtor of Reconstructor has been called." << std::endl;
}
