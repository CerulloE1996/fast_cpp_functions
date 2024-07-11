
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


 


// see https://stackoverflow.com/questions/39587752/difference-between-ldexp1-x-and-exp2x
/* For a in [0.5, 4), compute a * 2**i, -250 < i < 250 */
// Note: this function is  the same as the one in the link above, but for * double * instead of * float *.  
// [[Rcpp::export]]
double fast_ldexp (double a, int i)
{ 
  int64_t ia = ( (uint64_t)i << 52) + __double_as_int (a); // scale by 2**i
  a = __int_as_double (ia);
  if ((unsigned int)(i + 1021) > 500) { // |i| > 125
    i = (i ^ (1021 << 52)) - i; // ((i < 0) ? -125 : 125) << 52
    a = __int_as_double (ia - i); // scale by 2**(+/-125)
    a = a * __int_as_double ((1023 << 52) + i); // scale by 2**(+/-(i%125))
  }
  return a;
}


// see https://stackoverflow.com/questions/39587752/difference-between-ldexp1-x-and-exp2x
// Note: this function is  the same as the one in the link above, but for * double * instead of * float *. 
// [[Rcpp::export]]
double fast_exp_double(double a)
{  
  
  a  =  1.442695040888963387 * a;
  const double cvt = 8106479329266893.0 ; //  ldexp(1.8, 52) ; //  12582912.0; // 0x1.8p23
  double f, r;
  int i;
  
  // exp2(a) = exp2(i + f); i = rint (a)
  r = (a + cvt) - cvt;
  f = a - r;
  i = (int)r;
  // approximate exp2(f) on interval [-0.5,+0.5]
  r =            0.000153720378875732421875;  // 0x1.426000p-13f
  r = fma (r, f, 0.00133903871756047010422); // 0x1.5f055ep-10f
  r = fma (r, f, 0.00961817800998687744141); // 0x1.3b2b20p-07f
  r = fma (r, f, 0.0555036030709743499756); // 0x1.c6af7ep-05f
  r = fma (r, f, 0.240226522088050842285); // 0x1.ebfbe2p-03f
  r = fma (r, f, 0.693147182464599609375); // 0x1.62e430p-01f
  r = fma (r, f, 1.0); // 0x1.000000p+00f
  // exp2(a) = 2**i * exp2(f);
  
  r = fast_ldexp (r, i);
  
  if (!(abs (a) < 150.0)) {
    const double large  = fast_ldexp(1.0, -1024) ;  // 0x1.0p127 // = 1.70141184e38; 
    r = a + a; // handle NaNs
    if (a < 0.0) r = 0.0;
    if (a > 0.0) r = large * large; // + INF
  }
  return r;
  
}


// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, 1  > fast_exp_double_Eigen(  Eigen::Array<double, -1, 1  > x) {
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) = fast_exp_double(x(i));
  }    
  
  return x; 
  
}





// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  > fast_exp_double_Eigen_mat(  Eigen::Array<double, -1, -1  > x) {
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) = fast_exp_double(x(i, j));
    }
  }
  
  return x; 
  
}


