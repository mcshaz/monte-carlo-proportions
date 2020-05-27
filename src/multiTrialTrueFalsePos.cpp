#include <Rcpp.h>
#include <fstream>
#include <limits>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <mutex>
#include <algorithm>
#include "FisherRepo.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;

// [[Rcpp::export]]
void multiTrialTrueFalsePos(NumericVector baselineRisks,
                                 IntegerVector participantsPerArm,
                                 double absRRStep,
                                 int monteCarloRuns,
                                 Rcpp::String path) {

#ifdef _OPENMP
  omp_set_num_threads(std::min(omp_get_max_threads(), (int)participantsPerArm.size()));
#endif
  std::mutex outfileMutex;
  std::ofstream outfile(path.get_cstring());
  if (!outfile.is_open()) {
    Rcpp::stop("Unable to write to file\n");
  }
  outfile << "trialN\triskCtrl\triskRx\tnonSig\tRxHarm\tRxBenefit\n";
  const double epsilon = std::numeric_limits<double>::epsilon();
  unsigned int riskSteps = 0;
  for (double& br : baselineRisks) {
    riskSteps += (br + epsilon) / absRRStep;
    if ((br - absRRStep * (double)riskSteps) > epsilon) {
      ++riskSteps;
    }
  }
  Progress p(participantsPerArm.size() * riskSteps, true);

  const size_t runLimit = monteCarloRuns; // micro optimisation
// #pragma omp parallel for schedule(dynamic)
  for(const int& partNo :participantsPerArm)
  {
    FisherRepo repo(partNo);
    for (const double& baseRisk: baselineRisks)
    {
      if (!p.is_aborted()) { // the only way to exit an OpenMP loop
        NumericVector baselineOutcomes = rbinom(monteCarloRuns, partNo, baseRisk);
        for (double intervRisk = baseRisk; intervRisk > epsilon; intervRisk -= absRRStep)
        {
          int nonSig = 0;
          int rxHarm = 0;
          int rxBenefit = 0;
          NumericVector intervOutcomes = rbinom(monteCarloRuns, partNo, intervRisk);
          for (size_t i = 0; i < runLimit; ++i) {
            if (repo.getP(intervOutcomes[i], baselineOutcomes[i]) >= 0.05){
              ++nonSig;
            } else if (intervOutcomes[i] > baselineOutcomes[i]) {
              ++rxHarm;
            } else {
              ++rxBenefit;
            }
          }
          // scoping the lock - i.e. free the lock for waiting threads as soon as file has written row of data
          {
            const std::lock_guard<std::mutex> lock(outfileMutex);
            // "trialN\triskCtrl\triskRx\tnonSig\tRxHarm\tRxBenefit\n"
            outfile << partNo * 2 << "\t" << baseRisk << "\t" << intervRisk << "\t" << nonSig << "\t" << rxHarm << "\t" << rxBenefit << "\n";
          }
          p.increment(); //increment progress bar
          // if no false positives or negatives results, there is no point decreasing risk reduction further
          if ((rxHarm == 0 && nonSig == 0) || Progress::check_abort()) {
            // to do increment progress bar further if skipping for 0
            break;
            // break not allowed in parallelised loops, but this should just break out to the surrounding loop, which is not parallelised
          }
        }
      }
    }
  }
  outfile.close();
}
