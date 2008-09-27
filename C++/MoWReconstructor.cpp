/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: MoWReconstructor.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:44 $
Version:   $Revision: 1.23 $
=========================================================================*/

#include "Reconstructor.h"
#include "MoWReconstructor.h"

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_trace.h>
#include "utilities.h"

void MoWReconstructor::DoReconstruction()
{
    try{
        // step1: prepare data
        LoadInputDataFiles();
        PrepareInputData();
        // step2: load imaging parameters and constants
        LoadImagingParams();
        LoadCommonConstants();
        // step3: load mow parameters
        LoadWishartParameters();
        // step 4: prepare output
        LoadOutputDataFiles();
        AssembleSPHBasis();
        // step 5: prepare linear system
        AssembleA();
        // step6: reconstruct
        PrepareOutputData();
        Reconstruct();
        Write();
    }
    catch(char *str){
        std::cout << "Exception raised: " << str << std::endl;
    }
    system("pause");
}


void MoWReconstructor::LoadWishartParameters()
{
    //(2) tessellation for computing weights; eigenvectors, eigenvalues,
    char f_eigenvector[BUFF_MAX];
    GetPrivateProfileString("MOW", "eigenvectors", NULL, f_eigenvector, 80, f_config);
    // N = ReadVectorsFromFile(f_eigenvector, ev);
    std::ifstream infile(f_eigenvector);
    ev.read_ascii(infile);
    assert(ev.cols()==3);
    N = ev.rows();
    std::cout << "Read " << N << " eigenvectors from " << f_eigenvector << std::endl;
    char s_eigenvalues[2][10]; // measured in mm^2/s
    GetPrivateProfileString("MOW", "eigenvalue1", NULL, s_eigenvalues[0], 10, f_config);
    GetPrivateProfileString("MOW", "eigenvalue2", NULL, s_eigenvalues[1], 10, f_config);
    // GetPrivateProfileString("Constants", "Lambda3", NULL, s_eigenvalues[2], 10, f_config);
    ew[0] = static_cast<float>(atof(s_eigenvalues[0]));
    ew[1] = static_cast<float>(atof(s_eigenvalues[1]));
    ew[2] = ew[1];
    std::cout << "Eigenvalues: [ ";
    for (int i = 0; i < 3; ++i)
    {
        std::cout << ew[i] << " ";
    }
    std::cout << "]" <<std::endl;
    char s_pvalue[BUFF_MAX] = "2";
    GetPrivateProfileString("MOW", "p", NULL, s_pvalue, 10, f_config);
    p = static_cast<float>(atof(s_pvalue));
    std::cout << "p: " << p << std::endl;
    //char s_cond[BUFF_MAX] = "22.2";
    //GetPrivateProfileString("MOW", "condition_number", NULL, s_cond, 10, f_config);
    //cond = static_cast<float>(atof(s_cond));
    //std::cout << "cond: " << cond << std::endl;
}


void MoWReconstructor::ComputeP(vnl_vector <float> &S, vnl_vector <float> &P)
{
    vnl_vector <float> weights;
    weights.set_size(ev.rows());
    char s_method[BUFF_MAX] = "dls";
    GetPrivateProfileString("MOW", "solver", NULL, s_method, 10, f_config);
    if (!strcmp(s_method, "dls"))
    {
        weights = Ainv*S;
    }
    else if (!strcmp(s_method, "dls"))
    {
        vnl_NNLS(A,  S, weights);
    }
    else
    {
        throw "The MOW solver has to be either 'dls' (damped least squares) or 'nnls' (non-negative least squares)."
        return;
    }
    int option =  GetPrivateProfileInt("MOW", "PorQBIorDSI", 0, f_config);
    //std::cout << "option: " << option << std::endl;

    // "0", P, "1", QBI, "2", DSI, "3" weights
    if (option==1)
    {P = R_QBI*weights;}
    else if (option==2)
    {P = R_DSI*weights;}
    else if (option == 3)
    {P = weights;  //!!! only when T==N
    }
    else
    { P = R*weights; }

#if 1
            for (unsigned int k = 0; k<R.rows(); ++k)
            {
                if (P[k] < 0)
                    P[k] = 0;
            }
            P = P/P.sum();
#endif

    // save weights
}




void MoWReconstructor:: AssembleA()
{
    int K;
    if (!bUseBMatrix)
        K = gradients.rows();
    else
        K = bmatrix.rows();
    int T = tessellation.rows();
    int N = ev.rows();

    if (N == 0 || K == 0 || T == 0)
    {
        throw " Not Ready for Assembling A ";
        return;
    }
    A.set_size(K, N);
    R.set_size(T, N);
    R_QBI.set_size(T,N);
    R_DSI.set_size(T,N);
    vnl_matrix <float> D(3, 3);
    vnl_matrix <float> Q;
    float factor = r*r/(4*t); // should depend on r and t;
    std::cout << " r^2/(4t) = " << factor << std::endl;
    for (int i = 0; i < N; ++i)
    {
        ev.set_row(i,ev.get_row(i).normalize());
        //assert normalized vector
        D.set_identity();
        D = ew[1]*D +(ew[0] - ew[1]) * outer_product(ev.get_row(i), ev.get_row(i));
        for (int j = 0; j < K; ++j)
        {
            if (bUseBMatrix)
            {
                vnl_vector_fixed <float, 6> v = bmatrix.get_row(j);
                float data[6];  // be careful about the order
                //  vnl_sym_matrix :
                //    [ 0 1 3 ]
                //    [ 1 2 4 ]
                //    [ 3 4 5 ]
                //  while in b_matrix
                //    [ 0 3 4 ]
                //    [ 3 1 5 ]
                //    [ 4 5 2 ]
                //  very annoying !!
                data[0] = v[0];
                data[1] = v[3];
                data[2] = v[1];
                data[3] = v[4];
                data[4] = v[5];
                data[5] = v[2];
                Q = D*vnl_sym_matrix < float>(data, 3).as_matrix();
            }
            else
            {
                vnl_vector_fixed < float, 3> v = gradients.get_row(j);
                Q = b[j]*outer_product(v, v)*D;
            }
            float q = vnl_trace(Q);
            if (p>0)
                A(j, i) = pow(1 + q, -(float)p);
            else
                A(j, i) = exp(-q);
        }
        D.set_identity();
        D =(1/ew[1])*D +(1/(ew[0]) - 1/(ew[1])) * outer_product(ev.get_row(i), ev.get_row(i));
        for (int j = 0; j < T; ++j)
        {
            vnl_vector_fixed < float, 3> v = tessellation.get_row(j);
            Q = outer_product(v, v)*D;
            // float q =(1/ew[0]) * v1.squared_magnitude() +(1/ew[1]) * v2.squared_magnitude();
            float q = vnl_trace(Q);
            R(j, i) = exp(-factor*q) / sqrt(ew[0]*ew[1]*ew[1]);
            R_QBI(j,i) = 1.0f/sqrt(q);
            R_DSI(j,i) = 1.0f/(sqrt(q)*sqrt(q)*sqrt(q));

        }
    }

        /* compute the mutual coherence of A */
    /*vnl_matrix<float> G = A.transpose()*A;
        int m= G.rows();
    std::cout << " size of Gram matrix is " << m << std::endl;
    float mc = -1;
        for (int k=0;k<m;++k)
      {
        for (int j=k+1;j<m;++j)
          {
        float Gkj = fabs(G(k,j)/sqrt(G(k,k)*G(j,j)));

                if (Gkj>mc)
                 { mc = Gkj;
         std::cout << "mc = " << mc << std::endl;
                 }
          }
      }
    std::cout << "mutual coherence = " << mc << std::endl;
    */
    Ainv.set_size(N, K);
    //std::cout<< "K= " <<K << " N= " << N << std::endl;
    //vnl_matrix < float> I(K, K);
    //I.set_identity();
    // INVA = vnl_matrix_inverse<float>(A) * I;
    //std::cout << "R.rows = " << R.rows() << std::endl;
    // Compute the singular value decomposition of A.
    // In VXL this is implemented as an class extending the basic matrix class.
    vnl_svd < float> mySVD(A);
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

    std::cout << "sigma_max = " << mySVD.sigma_max() << " sigma_min = " << mySVD.sigma_min() << std::endl;

    int rnk = mySVD.rank();

    vnl_matrix < float> Winverse(myWinv.rows(), myWinv.columns());
    Winverse.fill(0);

    //float lambda = cond; // mySVD.sigma_max() / cond;
    float lambda = cond;
    for (int i = 0; i < rnk;++i)
    {
        Winverse(i, i)= myW(i,i) / (myW(i, i)*myW(i,i) + lambda*lambda);
    }
    Ainv = myV * Winverse * myU.conjugate_transpose();

    std::cout << "Assemble A OK!" << std::endl;
}




