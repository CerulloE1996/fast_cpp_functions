// [[Rcpp::depends(RcppEigen)]]  


#include <RcppEigen.h> 
#include <Rcpp.h> 


#include <unsupported/Eigen/SpecialFunctions>


#define EIGEN_USE_MKL_ALL
#include "Eigen/src/Core/util/MKL_support.h"



// [[Rcpp::plugins(cpp17)]]      

using namespace Rcpp;    
using namespace Eigen;  


// [[Rcpp::export]]
double  fn_ld_exp_bitshift_1_double(long long b)  { 
  
  return ldexp(1.0, b); 
  
}   




// [[Rcpp::export]]
float  fn_ld_exp_bitshift_1_float(int b)  { 
  
  return ldexpf(1.0f, b); 
  
}  



// [[Rcpp::export]]
double  fn_ld_exp_bitshift_a_double(double a, int b)  { 
  
  return ldexp(a, b); 
  
}    




// [[Rcpp::export]]
float  fn_ld_exp_bitshift_a_float(float a, int b)  { 
  
  return ldexpf(a, b); 
  
}   



double __int_as_double (int64_t a) { double r; memcpy (&r, &a, sizeof r); return r;}
int64_t __double_as_int (double a) { int64_t r; memcpy (&r, &a, sizeof r); return r;}



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





/* compute natural logarithm, maximum error 0.85089 ulps */
double fast_log_double (double a)
{
  double i, m, r, s, t;
  int64_t e;
        

        i = 0.0;
        if (a < 2.225074e-308){ // 0x1.0p-126
          a = a * 8388608.0 ; // 0x1.0p+23
          i = -52.0;
        }
        e = (__double_as_int (a) - 0x3fe5555555555555 )    &   0xFFF0000000000000   ;
        m = __int_as_double (__double_as_int (a) - e);
        i = fma ((double)e, 0.000000000000000222044604925031308085, i); // 0x1.0p-52
        /* m in [2/3, 4/3] */
        m = m - 1.0;
        s = m * m;
        /* Compute log1p(m) for m in [-1/3, 1/3] */
        r =             -0.13031005859375;  // -0x1.0ae000p-3
        t =              0.140869140625;  //  0x1.208000p-3
        r = fma (r, s, -0.121483512222766876221); // -0x1.f198b2p-4
        t = fma (t, s,  0.139814853668212890625); //  0x1.1e5740p-3
        r = fma (r, s, -0.166846126317977905273); // -0x1.55b36cp-3
        t = fma (t, s,  0.200120344758033752441); //  0x1.99d8b2p-3
        r = fma (r, s, -0.249996200203895568848); // -0x1.fffe02p-3
        r = fma (t, m, r);
        r = fma (r, m,  0.333331972360610961914); //  0x1.5554fap-2
        r = fma (r, m, -0.5); // -0x1.000000p-1  
        r = fma (r, s, m);
        r = fma (i,  0.693147182464599609375, r); //  0x1.62e430p-1 // log(2)
        if (!((a > 0.0) && (a < INFINITY))) {
          r = a + a;  // silence NaNs if necessary
          if (a  < 0.0) r = INFINITY - INFINITY; //  NaN
          if (a == 0.0) r = -INFINITY;
        }
        return r;
}



// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, 1  > fast_log_double_Eigen(  Eigen::Array<double, -1, 1  > x)
{ 
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) =  fast_log_double(x(i));
  }   
  
  return x; 
  
}    

// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  >  fast_log_double_Eigen_mat(  Eigen::Array<double, -1, -1  > x)
{   
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) =  fast_log_double(x(i, j));
    }  
  }
  
  return x; 
  
}    


