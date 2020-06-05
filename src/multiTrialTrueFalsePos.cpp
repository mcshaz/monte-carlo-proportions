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
#include "FisherRepo.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;
typedef pcg32_k2_fast chosenRNG; // pcg32_k64
inline static chosenRNG getMCRNG() {
  // Seed with a real random value, if available
  pcg_extras::seed_seq_from<std::random_device> seed_source;

  // Make a random number engine - could also forgo seed_source and use memory address with pcg32_unique
  chosenRNG rng(seed_source);
  return rng;
}

void multiTrialTrueFalsePos(std::vector<double> &baselineRisks,
                            std::vector<unsigned int> &participantsPerArm,
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
  unsigned int riskSteps = 0;
  for (double& br : baselineRisks) {
    unsigned int stp = (br + epsilon) / absRRStep;
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
    const unsigned int partNo = participantsPerArm[pni];
    FisherRepo repo(partNo);
    chosenRNG mcrng = getMCRNG();
    std::unique_ptr<chosenRNG> rngChkPtr = nullptr;
    if (pni == participantsPerArm.size() - 1){
      rngChkPtr = std::make_unique<chosenRNG>(mcrng);
    }
    for (const double& baseRisk: baselineRisks)
    {
      std::binomial_distribution<unsigned int> cBinom(partNo, baseRisk);
      bool isPerfectTest = false;
      // NumericVector baselineOutcomes = br(mcrng); // rbinom(monteCarloRuns, partNo, baseRisk);
      for (double intervRisk = baseRisk; intervRisk > epsilon; intervRisk -= absRRStep)
      {
        std::binomial_distribution<unsigned int> iBinom(partNo, intervRisk);
        unsigned int nonSig = 0;
        unsigned int rxBenefit = 0;
        if (isPerfectTest) {
          rxBenefit = runLimit;
        } else {
        // NumericVector intervOutcomes = rbinom(monteCarloRuns, partNo, intervRisk);
          for (size_t i = 0; i < runLimit; ++i) {
            const unsigned int ir = iBinom(mcrng);
            const unsigned int cr = cBinom(mcrng);
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
    if (rngChkPtr) {
      Rcout << "Required 2^" << log2(mcrng - (*rngChkPtr)) << " of 2^" << rngChkPtr->period_pow2() << " random numbers\n";
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
  std::vector<unsigned int> ppa = Rcpp::as<std::vector<unsigned int>>(participantsPerArm);
  multiTrialTrueFalsePos(br,
                         ppa,
                         absRRStep,
                         monteCarloRuns,
                         path.get_cstring());
  Rcout << "completed\n";
}

