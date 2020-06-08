#include <Rcpp.h>
#include "FisherRepo.h"
#include "ChiRepo.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]]

// [[Rcpp::export]]
NumericMatrix createFisherMat(int allocationsPerArm) {
  size_t outcomesPerArm = allocationsPerArm + 1;
  NumericMatrix returnVar(outcomesPerArm, outcomesPerArm);
  FisherRepo repo(allocationsPerArm);
  for (size_t j = 0; j < outcomesPerArm; ++j) {
    for (size_t i = 0; i < outcomesPerArm; ++i) {
      returnVar(j,i) = repo.getP(j, i);
    }
  }
  return returnVar;
}

// [[Rcpp::export]]
NumericMatrix createChi2Mat(int allocationsPerArm) {
  size_t outcomesPerArm = allocationsPerArm + 1;
  NumericMatrix returnVar(outcomesPerArm, outcomesPerArm);
  for (size_t j = 0; j < outcomesPerArm; ++j) {
    for (size_t i = 0; i < outcomesPerArm; ++i) {
      returnVar(j,i) = getChi(allocationsPerArm, j, i);
    }
  }
  return returnVar;
}
