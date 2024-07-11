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


float __int_as_float (int32_t a) { float r; memcpy (&r, &a, sizeof r); return r;}
int32_t __float_as_int (float a) { int32_t r; memcpy (&r, &a, sizeof r); return r;}



/* natural log on [0x1.f7a5ecp-127, 0x1.fffffep127]. Maximum relative error 9.4529e-5 */
// see: https://stackoverflow.com/questions/39821367/very-fast-approximate-logarithm-natural-log-function-in-c
// Note: this function is  the same as the one in the link above, but for * double * instead of * float *.  

// [[Rcpp::export]]
double fast_log_approx_double(double a)
{ 
  double m, r, s, t, i, f;
  int64_t e;
  
  e = (__double_as_int (a) - 0x3fe5555555555555 )    &   0xFFF0000000000000   ;
  // 0x3fe5555555555555
  m = __int_as_double (__double_as_int (a) - e);
  //  i = double(e) * (double)1.19209290e-7; // 0x1.0p-23
  i = (double)e *   0.000000000000000222044604925031308085 ; // ldexp(1.0, -52) ; 
  //  return(i); 
  
  /* m in [2/3, 4/3] */
  f = m - 1.0;
  s = f * f;
  /* Compute log1p(f) for f in [-1/3, 1/3] */
  r = fma (0.230836749076843261719, f, -0.279208570718765258789); // 0x1.d8c0f0p-3, -0x1.1de8dap-2
  t = fma (0.331826031208038330078, f, -0.498910337686538696289); // 0x1.53ca34p-2, -0x1.fee25ap-2
  r = fma (r, s, t);
  r = fma (r, s, f);
  r = fma (i, 0.693147182464599609375, r); // 0x1.62e430p-1 // log(2)
  return r; 
}  



// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, 1  > fast_log_approx_double_Eigen(  Eigen::Array<double, -1, 1  > x)
{  
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) =  fast_log_approx_double(x(i));
  }   
  
  return x; 
  
}   
 
// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  >  fast_log_approx_double_Eigen_mat(  Eigen::Array<double, -1, -1  > x)
{   
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) =  fast_log_approx_double(x(i, j));
    }  
  }
  
  return x; 
  
}    
