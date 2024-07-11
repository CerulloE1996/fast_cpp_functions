

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








// From:  https://stackoverflow.com/questions/65554112/fast-double-exp2-function-in-c/65562273#65562273s
/* compute 2**p, for p in [-1022, 1024). Maximum relative error: 4.93e-5. RMS error: 9.91e-6 */  
// Note:  modified -  uses some FMA operations 
// [[Rcpp::export]]
double fast_exp_approx_double_wo_checks (double p)
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
Eigen::Array<double, -1, 1  > fast_exp_approx_double_wo_checks_Eigen(  Eigen::Array<double, -1, 1  > x)
{  
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) =  fast_exp_approx_double_wo_checks(x(i));
  }    
  
  return x; 
  
}   

// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  >  fast_exp_approx_double_wo_checks_Eigen_mat(  Eigen::Array<double, -1, -1  > x)
{    
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) =  fast_exp_approx_double_wo_checks(x(i, j));
    }  
  }
  
  return x; 
  
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



// see https://stackoverflow.com/questions/39587752/difference-between-ldexp1-x-and-exp2x
// Note: this function is  the same as the one in the link above, but for * double * instead of * float *. 
// [[Rcpp::export]]
double fast_exp_double_wo_checks(double a)
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
  
  // if (!(abs (a) < 150.0)) {
  //   const double large  = fast_ldexp(1.0, -1024) ;  // 0x1.0p127 // = 1.70141184e38;
  //   r = a + a; // handle NaNs
  //   if (a < 0.0) r = 0.0;
  //   if (a > 0.0) r = large * large; // + INF
  // }
  return r;
  
}


// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, 1  > fast_exp_double_wo_checks_Eigen(  Eigen::Array<double, -1, 1  > x) {
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) = fast_exp_double_wo_checks(x(i));
  }    
  
  return x; 
  
}





// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  > fast_exp_double_wo_checks_Eigen_mat(  Eigen::Array<double, -1, -1  > x) {
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) = fast_exp_double_wo_checks(x(i, j));
    }
  }
  
  return x; 
  
}



/* natural log on [0x1.f7a5ecp-127, 0x1.fffffep127]. Maximum relative error 9.4529e-5 */
// see: https://stackoverflow.com/questions/39821367/very-fast-approximate-logarithm-natural-log-function-in-c
// Note: this function is  the same as the one in the link above, but for * double * instead of * float *.  

// [[Rcpp::export]]
double fast_log_approx_double_wo_checks(double a)
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
Eigen::Array<double, -1, 1  > fast_log_approx_double_wo_checks_Eigen(  Eigen::Array<double, -1, 1  > x)
{  
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) =  fast_log_approx_double_wo_checks(x(i));
  }   
  
  return x; 
  
}   

// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  >  fast_log_approx_double_wo_checks_Eigen_mat(  Eigen::Array<double, -1, -1  > x)
{   
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) =  fast_log_approx_double_wo_checks(x(i, j));
    }  
  }
  
  return x; 
  
}    



/* compute natural logarithm, maximum error 0.85089 ulps */
double fast_log_double  (double a)
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
  // if (!((a > 0.0) && (a < INFINITY))) {
  //   r = a + a;  // silence NaNs if necessary  
  //   if (a  < 0.0) r = INFINITY - INFINITY; //  NaN
  //   if (a == 0.0) r = -INFINITY;
  // }
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






/* compute natural logarithm, maximum error 0.85089 ulps */
double fast_log_double_wo_checks (double a)
{
  double i, m, r, s, t;
  int64_t e;
  
  i = 0.0;
  // if (a < 2.225074e-308){ // 0x1.0p-126
  //   a = a * 8388608.0 ; // 0x1.0p+23
  //   i = -52.0;
  // }
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
  return r;
}



// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, 1  > fast_log_double_wo_checks_Eigen(  Eigen::Array<double, -1, 1  > x)
{ 
  
  for (int i = 0; i < x.rows(); ++i) {
    x(i) =  fast_log_double_wo_checks(x(i));
  }   
  
  return x; 
  
}    

// Note: Compiler needs  to auto-vectorise for the following to be fast
// [[Rcpp::export]]
Eigen::Array<double, -1, -1  >  fast_log_double_wo_checks_Eigen_mat(  Eigen::Array<double, -1, -1  > x)
{   
  
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x(i, j) =  fast_log_double_wo_checks(x(i, j));
    }  
  }
  
  return x; 
  
}    




