#ifndef Chi2Repo_H
#define Chi2Repo_H

#include <Rcpp.h>

inline double getChi(const double allocPerArm, const double a, const double c){
  const double f = allocPerArm * (a - c);
  const double chi2 = f * f * 2 * allocPerArm / (allocPerArm * allocPerArm * (a + c) * (2 * allocPerArm - a - c));
  return R::pchisq(chi2, 1, false, false);
}

#endif
