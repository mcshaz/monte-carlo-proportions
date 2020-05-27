\name{multiTrialTrueFalsePos}
\alias{multiTrialTrueFalsePos}
\docType{package}
\title{
Monte Carlo Trials with outcome counts random binomial of given probability
}
\description{
Monte Carlo Trials with outcome counts randomly taken from a binomial distribution of a given probability. Fisher's exact p values then used
}
\usage{
multiTrialTrueFalsePos(baselineRisks, participantsPerArm, absRRStep, monteCarloRuns)
}
\arguments{
  \item{baselineRisks}{a vector of risks bounded by (0,1) which will be the underlying risk in the control group}
  \item{participantsPerArm}{a vector of integers denoting how many participants in each arm (N/2) for each set of monte carlo runs}
  \item{absRRStep}{a real value bounded by (0,1) denoting how much to reduce the risk in the intervention arm each run}
  \item{monteCarloRuns}{the number of trials run for each combination of the arguments above}
}
\value{
A DataFrame containing:
$trialN - the number of participants in those trials (with half allocated to each arm)
$riskCtrl - risk for that run in the control group
$riskRx - risk for that run in the intervention/treatment group
$nonSig - count of the number of times p was >= 0.05
$RxHarm - count of the number of time p was < 0.05 indicating control had a lower risk of the outcome
$RxBenefit - count of the number of time p was < 0.05 indicating treatment group had a lower risk of the outcome
}
\examples{
\dontrun{
multiTrialTrueFalsePos(c(0.1, 0.2, 0.3), c(50L, 100L, 200L), 0.05, 1000L)
}
}
