/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: DOTReconstructor.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/19 07:59:58 $
Version:   $Revision: 1.13 $
=========================================================================*/


#include "Reconstructor.h"
#include "DOTReconstructor.h"

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#include <math.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_trace.h>
#include "utilities.h"

#include "special_functions.h"

//#define PI   3.141592653589f

inline int FLOOR(float x)  {return x>0.0f ? (int)x : -((int)(1.0f-x));}
inline int ROUND(float x)  {return x>0.0f ? (int)(x+0.5f) : -(int)(0.5f-x);}
inline int SQR(int x)    {return x*x;}
inline double SQR(double x)    {return x*x;}
inline float SQR(float x)    {return x*x;}

inline int CUBE(int x)    {return x*x*x;}
inline double CUBE(double x)    {return x*x*x;}
inline float CUBE(float x)    {return x*x*x;}





inline int ipow (int x, int i)
{
	int result = 1;
	while(i){
		if (i & 1) result *= x;
		x *= x;
		i >>= 1;
	}
	return result;
}

inline float ipow (float x, int i)
{
	float result = 1.0f;
	while(i){
		if (i & 1)	result *= x;
		x *= x;
		i >>= 1;
	}
	return result;
}

inline double ipow (double x, int i)
{
	double result = 1.0;
	while(i){
		if (i & 1)    result *= x;
		x *= x;
		i >>= 1;
	}
	return result;
}

double Al(int l, double beta)
{
	switch(l)
	{
		case 0:
			return 1.0;	break;
		case 2:
			return -(1.0+6.0/SQR(beta));	break;
		case 4:
			return 1.0 + 20.0/SQR(beta) + 210.0/ipow(beta, 4);      break;
		case 6:
			return -(1.0 + 42.0/SQR(beta) + (1575.0/2.0)/ipow(beta, 4) + 10395.0/ipow(beta, 6) );
			break;
		case 8:
			return 1.0 + 72.0/SQR(beta) + (10395.0/4.0)/ipow(beta, 4)  + 45045.0/ipow(beta, 6) + 675675.0/ipow(beta, 8);
			break;
		default:
			return 0.0;
	}
}


double Bl(int l, double beta)
{
	switch(l)
	{
		case 0:
			return 0.0;	break;
		case 2:
			return 3.0;	break;
		case 4:
			return (15.0/2.0)*(1.0 - 14.0/SQR(beta));	break;
		case 6:
			return (105.0/8.0)*(1.0 - 36.0/SQR(beta) + 396.0/ipow(beta, 4) );
			break;
		case 8:
			return (315.0/16.0)*(1.0 - 66.0/SQR(beta) + 1716.0/ipow(beta, 4) - 17160.0/ipow(beta, 6));
			break;
		default:
			return 0.0;
	}
}

double Pl(int n, double x)
{
	switch(n){
		case 0:
			return 1.0;	break;
		case 2:
			return 0.5*(3.0*SQR(x) - 1.0);		break;
		case 4:
			return (1.0/8.0)*(35.0*ipow(x,4) - 30.0*SQR(x) + 3.0);      break;
		case 6:
			return (1.0/16.0)*(231.0*ipow(x, 6) -315.0* ipow(x,4) + 105.0*ipow(x,2) - 5.0);
			break;
		case 8:
			return (1.0/128.0)*(6435.0*ipow(x,8) - 12012.0*ipow(x,6) + 6930.0*ipow(x,4) - 1260.0*ipow(x,2) + 35.0);
			break;
	}
	return 0.0;
}


void
DOTReconstructor::
AssembleM()
{
	int K = gradients.rows();
	int T = tessellation.rows();

	if (K == 0 || T == 0)
	{
		throw " Not Ready for Assembling M.";
		return;
	}

	M.set_size(T, K*(degreeDOT/2 + 1));
	vnl_matrix <float> DP;
	DP.set_size(T,K);
	for (int out = 0; out < T; out++)
	{
		for (int in = 0; in < K; in++)
		{
			DP(out, in) = dot_product(gradients.get_row(in), tessellation.get_row(out));
			//printf("DP(out,in)=%.5f\n",DP(out,in));
			
		}	
		//system("pause");
	}	
	
	for (int l = 0; l <=degreeDOT; l+=2)
	{
		float coeff = pow(-1.0, l/2) *(2.0*l + 1); 
		//pow(-1,l/2) = sign( 0.5 - (l/2)%2 )
		for (int in = 0; in < K; in++)
		{
			float coeff_wj = coeff*wj[in];
			for (int out = 0; out < T; ++out)
			{
				// how to compute wj[in] 
				M(out, (l/2)*K+in) = coeff_wj * Pl(l, DP(out, in));
			}
		}
	}
	//vcl_cout << M << vcl_endl;
	//std::ofstream outfile("C:/M.txt",std::ios_base::out);
	//M.print(outfile);
	//system("pause");

}

void DOTReconstructor::
LoadDOTConstants() // const char* f_config
{
	degreeDOT = GetPrivateProfileInt("[DOT]", "degree", 6, f_config);
	if (degreeDOT < 4)
		degreeDOT = 6;
	std::cout << "The degree of Spherical Harmonics for DOT is " << degreeDOT << std::endl;
	
	
	char f_wjvalues[80] = "\0"; // measured in s/mm^2
	GetPrivateProfileString("[DOT]", "weights", NULL, f_wjvalues, 80, f_config);
	if (strlen(f_wjvalues) == 0)
	{
		std::cerr << "No weights file found. Uniform weights are used." <<std::endl;
		int n = gradients.rows();
		wj.set_size(n);
		wj.fill(1);
	}
	else
	{
		std::cout << "Reading weights from " << f_wjvalues << " ... " << std::endl;
		//K = ReadValuesFromFile(f_wjvalues, wj);
		std::ifstream infile(f_wjvalues);
		wj.read_ascii(infile);
		assert(wj.size()>0);	
		std::cout <<  wj.size() << " weights were loaded from " << f_wjvalues << std::endl;
	}
	
}

void DOTReconstructor::
ComputeP(vnl_vector <float> &S, vnl_vector <float> &P)
{
	ComputeADC(S);
	vnl_vector <float> D = S;
	vnl_vector <float> Z;
	Z.set_size(nSamples*(1+degreeDOT/2));
	for (int j = 0; j < nSamples; ++j)
	{
		float beta = r/sqrt(D[j] * t);
		float expB = exp(-0.25*SQR(beta));
		float erfB = erf(0.5*beta);
		float f2 = erfB/(4.0*PI*CUBE(r)); 
		//printf("%.2f %.2f %.2f\n",beta, erfB, f2);

		float powD = pow(4.0*PI*D[j]*t, 1.5);
		float f1 = expB/powD;		
		for (int l = 0; l <= degreeDOT; l += 2)
		{
			//printf("%.2f %.2f %.2f %.2f\n",Al(l, beta),Bl(l, beta),f1,f2);
			//system("pause");
			Z[(l/2)*nSamples+j] = Al(l, beta)*f1 + Bl(l, beta)*f2;
		}
	}
	//std::ofstream outfile("C:/Z.txt",std::ios_base::out);
	//Z.print(outfile);
	//vcl_cout << Z << vcl_endl;
	//system("pause");

	P = M*Z;
	//vcl_cout << P << vcl_endl;
	//system("pause");
}


void DOTReconstructor::
DoReconstruction() // const char* f_config
{
	try{
		// step1: prepare data
		LoadInputDataFiles();
		PrepareInputData();
		// step2: 
		LoadImagingParams();
		LoadCommonConstants();
		// step3:
		LoadDOTConstants();
		//step 4:
		LoadOutputDataFiles();
		AssembleSPHBasis();
		//step 5: prepare linear system
		AssembleM();
		//step6: reconstruct
		PrepareOutputData();
		Reconstruct();
		Write();	
	}
	catch(char *str){
		std::cout << "Exception raised: " << str << std::endl;
	}
	system("pause");
}


