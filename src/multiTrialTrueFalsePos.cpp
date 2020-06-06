#include <Rcpp.h>
#include <fstream>
#include <limits>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <mutex>
#include <algorithm>
#include <random>
#include <memory>
#include <cmath>
#include <pcg-cpp/pcg_random.hpp>
#include <randutils/randutils.hpp>
#include "FisherRepo.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;
typedef pcg32 chosenRNG; // pcg32_k64

void multiTrialTrueFalsePos(std::vector<double> &baselineRisks,
                            std::vector<unsigned> &participantsPerArm,
                            double absRRStep,
                            int monteCarloRuns,
                            const char* path) {
#ifdef _OPENMP
  omp_set_num_threads(std::min(omp_get_max_threads(), (int)participantsPerArm.size()));
#endif
  std::mutex outfileMutex;
  std::ofstream outfile(path);
  if (!outfile.is_open()) {
    Rcpp::stop("Unable to write to file\n");
  }
  outfile << "trialN\triskCtrl\triskRx\tnonSig\tRxBenefit\n";
  const double epsilon = std::numeric_limits<double>::epsilon();
  unsigned riskSteps = 0;
  for (double& br : baselineRisks) {
    unsigned stp = (br + epsilon) / absRRStep;
    if ((br - absRRStep * (double)stp) > epsilon) {
      ++riskSteps;
    }
    riskSteps += stp;
  }
  Progress p(participantsPerArm.size() * riskSteps, true);

  const size_t runLimit = monteCarloRuns; // micro optimisation
  #pragma omp parallel for schedule(dynamic)
  for(size_t pni = 0; pni < participantsPerArm.size(); ++pni)
  {
    const unsigned partNo = participantsPerArm[pni];
    FisherRepo repo(partNo);
    randutils::auto_seed_128 seeds;
    chosenRNG mcrng{seeds};
    pcg32::state_type streamNo = 0;
    for (const double& baseRisk: baselineRisks)
    {
      std::binomial_distribution<unsigned> cBinom(partNo, baseRisk);
      bool isPerfectTest = false;
      // NumericVector baselineOutcomes = br(mcrng); // rbinom(monteCarloRuns, partNo, baseRisk);
      for (double intervRisk = baseRisk; intervRisk > epsilon; intervRisk -= absRRStep)
      {
        std::binomial_distribution<unsigned> iBinom(partNo, intervRisk);
        unsigned nonSig = 0;
        unsigned rxBenefit = 0;
        if (isPerfectTest) {
          rxBenefit = runLimit;
        } else {
        // NumericVector intervOutcomes = rbinom(monteCarloRuns, partNo, intervRisk);
          mcrng.set_stream(++streamNo);
          for (size_t i = 0; i < runLimit; ++i) {
            const unsigned ir = iBinom(mcrng);
            const unsigned cr = cBinom(mcrng);
            if (repo.getP(ir, cr) >= 0.05){
              ++nonSig;
            } else if (ir < cr) {
              ++rxBenefit;
            }
          }
          // if no false positives or negatives results, there is no point decreasing risk reduction further
          if (rxBenefit == runLimit) {
            isPerfectTest = true;
          }
        }
        // scoping the lock - i.e. free the lock for waiting threads as soon as file has written row of data
        {
          const std::lock_guard<std::mutex> lock(outfileMutex);
          // "trialN\triskCtrl\triskRx\tnonSig\tRxBenefit\n"
          outfile << partNo * 2 << "\t" << baseRisk << "\t" << intervRisk << "\t" << nonSig << "\t" << rxBenefit << "\n";
        }
        p.increment(); //increment progress bar
        if (Progress::check_abort()) {
          Rcpp::stop("Aborted by User\n");
        }
      }
    }
    if (mcrng.wrapped()) {
      Rcerr << "Insufficient period to chosen random number generator leading to wrapping around 0 state (period of 2^" << mcrng.period_pow2() << ")\n";
    }
  }
  outfile.close();
}
// [[Rcpp::export]]
void multiTrialTrueFalsePos(NumericVector baselineRisks,
                            IntegerVector participantsPerArm,
                            double absRRStep,
                            int monteCarloRuns,
                            Rcpp::String path) {
  std::vector<double> br = Rcpp::as<std::vector<double>>(baselineRisks);
  std::vector<unsigned> ppa = Rcpp::as<std::vector<unsigned>>(participantsPerArm);
  multiTrialTrueFalsePos(br,
                         ppa,
                         absRRStep,
                         monteCarloRuns,
                         path.get_cstring());
  Rcout << "completed\n";
}

