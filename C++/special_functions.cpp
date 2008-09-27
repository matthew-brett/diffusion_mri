/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: special_functions.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2008/09/10 02:12:45 $
Version:   $Revision: 1.1 $
=========================================================================*/

#include <iostream>
#include <cmath>

#include "special_functions.h"

/*
  associated Legendre polynomial
  P_l^m(cos(t))
  =(-1)^{l + m} /(2^l l!) sin^m(t)(d/d cos(t))^{l + m} sin^{2l}(t)
*/
double pLegendre(int l, int m, double t)
{
/*
    compute the associated Legendre polynomials P_l^m(cos(t))
    of degree l and order m
    in double precision by the recurrence relation

(input)
    l(int) degree s.t. l >= 0
    m(int) order s.t. -l <= m <=l
    t(double) variable
 */
    double cs, sn, y, y_1, y_2, yy, c, d;
    int mm, i, k;
    if (l < 0){
        printf("l=%d. l must be non-negative.\n");
        return 0;
    }
    if (m < -l || m>l){
        printf("m=%d. m must be -l <= m <= l.\n");
        return 0;
    }
  /*
    compute P_l^m(x) by the recurrence relation
(l - m)P_l^m(x) = x(2l - 1)P_{l - 1}^m(x) -(l + m - 1)P_{l - 2}^m(x)
    with
    P_m^m(x) =(-1)^m(2m - 1)!!(1 - x)^{m/2},
    P_{m + 1}^m(x) = x(2m + 1) P_m^m(x).
  */
    cs = cos(t);
    sn = sin(t);
    mm = m;                   /*   mm = |m|   */
    if (m < 0)  mm = - mm;
    y_1 = 1.0;
    for (i = 1; i <= mm; ++i)   y_1 *= - 1.0 *(2.0*i - 1) * sn;
    if (l == mm) {
        yy = y_1;
    }
    else  {
        y =(2.0*mm + 1.0) * cs * y_1;
        if (l == mm + 1) {
            yy = y;
        }
        else    {
            c = 2.0 * mm - 1.0;
            for (k = mm + 2; k <= l; ++k){
                y_2 = y_1;
                y_1 = y;
                d = 1.0 /(k - mm);
                y =(2.0 + c * d) * cs * y_1 -(1.0 + c * d) * y_2;
            }
            yy = y;
        }
    }
  /*
    In the case that m < 0,
    compute P_n^{-|m|}(x) by the formula
    P_l^{-|m|}(x) =(-1)^{|m|}((l-|m|)!/(l+|m|)!)^{1/2} P_l^{|m|}(x).
  */
    if (m < 0){
        for (i = l - mm + 1; i <= l + mm; ++i)  yy *= 1.0 / i;
        if (mm%2 == 1)    yy = - yy;
    }

    return yy;
}

double logGamma(double xx)
{
    double x, y, tmp, ser;
    static double cof[6]=
    {76.18009172947146, - 86.50532032941677,
    24.01409824083091, - 1.231739572450155,
    0.1208650973866179e-2, - 0.5395239384953e-5
    };
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -=(x + 0.5)*log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++)
        ser += cof[j]/++y;
    return -tmp + log(2.5066282746310005*ser/x);
}

double factorial(int n)
{
    if (n <= 1)
    {
        return 1.0;
    }
    return exp(logGamma(n + 1.0));
}


double erfc(double x)
{
   double t,z,ans;
   z=fabs(x);
   t=1.0/(1.0+0.5*z);
   ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
   t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
   t*(-0.82215223+t*0.17087277)))))))));
   return x >= 0.0 ? ans : 2.0-ans;
}

double erf(double x)
{
   return 1.0-erfc(x);
}


#if 0
//There is a mistake in the functions below. (found when processing the 1fib.flt using DOT)

static const double rel_error= 1E-12; //calculate 12 significant figures
//you can adjust rel_error to trade off between accuracy and speed
//but don't ask for > 15 figures (assuming usual 52 bit mantissa in a float)

double erfc(double x)
//erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
// = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
// = 1-erf(x)
//expression inside [] is a continued fraction so '+' means add to denominator only
{
    static const double one_sqrtpi= 0.564189583547756287; // 1/sqrt(pi)
    if (fabs(x) < 2.2) {
        return 1.0 - erf(x); //use series when fabs(x) < 2.2
    }
    //if (signbit(x)) { //continued fraction only valid for x>0
    if (x>0){
        return 2.0 - erfc(-x);
    }
    double a=1, b=x; //last two convergent numerators
    double c=x, d=x*x+0.5; //last two convergent denominators
    double q1, q2= b/d; //last two convergents (a/c and b/d)
    double n= 1.0, t;
    do {
        t= a*n+b*x;
        a= b;
        b= t;
        t= c*n+d*x;
        c= d;
        d= t;
        n+= 0.5;
        q1= q2;
        q2= b/d;
    } while (fabs(q1-q2)/q2 > rel_error);
    return one_sqrtpi*exp(-x*x)*q2;
}


double erf(double x)
//erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
// = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
// = 1-erfc(x)
{
    static const double two_sqrtpi= 1.128379167095512574; // 2/sqrt(pi)
    if (fabs(x) > 2.2) {
        return 1.0 - erfc(x); //use continued fraction when fabs(x) > 2.2
    }
    double sum= x, term= x, xsqr= x*x;
    int j= 1;
    do {
        term*= xsqr/j;
        sum-= term/(2*j+1);
        ++j;
        term*= xsqr/j;
        sum+= term/(2*j+1);
        ++j;
    } while (fabs(term/sum) > rel_error); // CORRECTED LINE
    return two_sqrtpi*sum;
}
#endif



