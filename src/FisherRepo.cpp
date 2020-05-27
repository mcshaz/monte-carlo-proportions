#include <vector>
#include <cmath>
#include <limits>
#include "FisherRepo.h"

using namespace std;

static vector<vector<double>> createFisherGrid(const size_t& allocationsPerArm) {
  size_t n = 2 * allocationsPerArm;
  vector<double> log2Factorials(n + 1);
  log2Factorials[0] = 0;
  for (size_t i = 1; i <= n; ++i) {
    log2Factorials[i] = log2Factorials[i-1] + log2(i);
  }

  vector<vector<double>> ps(allocationsPerArm + 1);
  const double smallestPwr = log2(numeric_limits<double>::denorm_min());
  const double halfway = allocationsPerArm / 2.0;

  for (size_t a = 0; a <= allocationsPerArm; ++a) {
    // casting as int truncates towards 0
    size_t blimit = halfway - abs(halfway - a);
    ps[a] = vector<double>(blimit + 1);
    for (size_t b = 0; b <= blimit; ++b) {
      const double powr = log2Factorials[a + b] + log2Factorials[n - a - b] + 2 * log2Factorials[allocationsPerArm]
        - log2Factorials[a] - log2Factorials[b] - log2Factorials[allocationsPerArm - a]
        - log2Factorials[allocationsPerArm - b] - log2Factorials[n];
      if (powr > smallestPwr) {
        ps[a][b] = powr < -1.0 // 2 ^ -1 = 0.5 which is then multiplied by 2
          ? 2 * exp2(powr)
          : 1.0;
      }
    }
    size_t j = a;
    size_t i = 0;
    double cum = ps[a][0];
    while(j != 0 && ps[--j].size() > ++i) {
      if (cum >= 1.0) {
        ps[j][i] = 1.0;
      } else {
        cum += ps[j][i];
        ps[j][i] = cum >= 1.0 ? 1.0 : cum;
      }
    };
  }
  return ps;
}

FisherRepo::FisherRepo(size_t allocationsPerArm):
  _ps(createFisherGrid(allocationsPerArm)){}
const double FisherRepo::getP(const size_t a, const size_t c){
  size_t i;
  size_t j;
  if (a > c) {
    j = a;
    i = c;
  } else {
    j = c;
    i = a;
  }
  if (i >= _ps[j].size()) {
    size_t tmp = j;
    j = _ps.size() - 1 - i;
    i = _ps.size() - 1 - tmp;
  }
  return _ps[j][i];
}
