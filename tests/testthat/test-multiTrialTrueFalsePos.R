test_that("multiTrialTrueFalsePos", {
  # note that the first Fisher matrix with any p < 0.05 is 4 allocations per arm
  baseRisk <- 0.3
  mcruns <- 10L
  riskRed <- 0.1
  participantsPerArm  <- 3L
  filenm <- "test.tsv"
  multiTrialTrueFalsePos(c(baseRisk), c(participantsPerArm), riskRed, mcruns, filenm)
  res <- read.table(filenm, header = TRUE)

  testthat::expect_that(res, is_a("data.frame"))
  rowCnt <- (baseRisk + 0.00001) %/% riskRed
  testthat::expect_equal(nrow(res), rowCnt)
  testthat::expect_equal(nrow(na.omit(res)), rowCnt)

  testthat::expect_equal(res$riskCtrl, rep(baseRisk, rowCnt))
  testthat::expect_equivalent(sort(res$riskRx), seq(riskRed, baseRisk, by = riskRed))
  testthat::expect_equal(res$nonSig, rep(mcruns, rowCnt))
  testthat::expect_equal(res$RxHarm, rep(0, rowCnt))
  testthat::expect_equal(res$RxBenefit, rep(0, rowCnt))
  testthat::expect_equal(res$trialN, rep(participantsPerArm * 2L, rowCnt))

#now test multi threading
  mcruns <- 200
  multiTrialTrueFalsePos(c(baseRisk, 0.35), c(participantsPerArm, 30, 70), riskRed, mcruns, filenm)
  res <- read.table(filenm, header = TRUE)
  rowCnt <- 21
  testthat::expect_equal(nrow(res), rowCnt)
  testthat::expect_equal(res$nonSig + res$RxHarm + res$RxBenefit, rep(mcruns, rowCnt))
  # example of real world use:
  # multiTrialTrueFalsePos(seq(0.03, 0.5, by=0.01), c(seq(20,200, by=5),seq(220,400, by=20), seq(450, 1000, by=50), seq(1200, 2000, by=200), seq(2500, 6000, by=500)), 0.01, 1000000, "millionRuns.tsv")

})

