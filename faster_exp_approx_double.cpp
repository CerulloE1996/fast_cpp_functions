
// [[Rcpp::depends(RcppEigen)]]  


#include <RcppEigen.h> 
#include <Rcpp.h> 


#include <unsupported/Eigen/SpecialFunctions>


#define EIGEN_USE_MKL_ALL
#include "Eigen/src/Core/util/MKL_support.h"



// [[Rcpp::plugins(cpp17)]]      

using namespace Rcpp;    
using namespace Eigen;  

 

double __int_as_double (int64_t a) { double r; memcpy (&r, &a, sizeof r); return r;}
int64_t __double_as_int (double a) { int64_t r; memcpy (&r, &a, sizeof r); return r;}




// From:  https://stackoverflow.com/questions/65554112/fast-double-exp2-function-in-c/65562273#65562273s
/* compute 2**p, for p in [-1022, 1024). Maximum relative error: 4.93e-5. RMS error: 9.91e-6 */
// [[Rcpp::export]]
double fast_exp_approx_double (double p)
{
  
  double res;
  
  p  =  1.442695040888963387 * p;
  p = (p < -1022) ? -1022 : p; // clamp below
  
  /* 2**p = 2**(w+z), with w an integer and z in [0, 1) */
  double w = floor(p); // integral part
  double z = p - w;     // fractional part
 
  double c3_recip = 1.0 /  ( 4.84257784485816955566 - z);
  double approx;
  approx   = fma(27.7283337116241455078, c3_recip, -5.72594201564788818359);
  approx   = fma(-0.49013227410614490509, z, approx);
  
  int64_t resi = ((1LL << 52) * (w + 1023 + approx));   /* assemble the exponent and mantissa components into final result */
  
  memcpy (&res, &resi, sizeof res);
  return res;
}

// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, 1  > fast_exp_approx_double_Eigen(  Eigen::Array<double, -1, 1  > x)
{  
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) =  fast_exp_approx_double(x(i));
  }    
  
  return x; 
  
}   
  
// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  >  fast_exp_approx_double_Eigen_mat(  Eigen::Array<double, -1, -1  > x)
{    
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) =  fast_exp_approx_double(x(i, j));
    }  
  }
  
  return x; 
   
}    
