//==============================================================================
// Package:  dnorm.cpp
// Usage:    #include "dnorm.cpp"
// Purpose:  neglogL for normal error
// Includes: admodel.h
// License:  BSD
// Notes:    Martell, S. 2009. ADMB Foundation Newsletter 1(4):3-4.
// History:  2011-02-21 Arni Magnusson added support for
//                      dvector-dvariable-dvariable, renamed std to sigma, moved
//                      pi outside functions, simplified comments
//           2010-06-19 Steve Martell created
//==============================================================================

#include <admodel.h>
#ifdef M_PI
double pi = M_PI;
#else
double pi = 3.14159265358979323844;
#endif



dvariable dnorm(const double& x, const dvariable& mu, const dvariable& sigma)
{
  return 0.5*log(2.0*pi) + log(sigma) + 0.5*square((x-mu)/sigma);
}



dvariable dnorm(const dvector& x, const dvar_vector& mu,
                const dvar_vector& sigma)
{
  int n = size_count(x);
  return 0.5*n*log(2.0*pi) + sum(log(sigma)) +
    0.5*sum(elem_div(square(x-mu),square(sigma)));
}



dvariable dnorm(const dvector& x, const dvariable& mu, const dvariable& sigma)
{
  int n = size_count(x);
  return 0.5*n*log(2.0*pi) + n*log(sigma) +
    sum(square(x-mu))/(2.0*square(sigma));
}
