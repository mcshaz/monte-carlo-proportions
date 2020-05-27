test_that("high performance Fisher p values", {
  alloc <- 8L

  fisherRes <- function() {
    nrow <- alloc + 1L
    outcomes <- matrix(0:(nrow^2 - 1),nrow = nrow, ncol =nrow) # all permutations. R also has a t(combn(0:10, 2)) & then add cbind(0:10, 0:10)

    prow <- function(r) {
      i <- r %/% nrow
      j <- r - nrow * i
      f <- fisher.test(matrix(c(j, alloc - j, i, alloc - i), 2L))
      return(unname(f$p.value))
    }
    return(matrix(t(sapply(outcomes,prow)), nrow=nrow, ncol=nrow))
  }

  cfm <- createFisherMat(alloc)
  testthat::expect_that( cfm, is_a("matrix"))
  testthat::expect_equal(dim(cfm), c(alloc + 1L, alloc + 1L))
  testthat::expect_equal(sum(cfm < 0) + sum(cfm > 1), 0)
  testthat::expect_equivalent(cfm, fisherRes())

  alloc <- 9L
  cfm <- createFisherMat(alloc)
  testthat::expect_equivalent(cfm, fisherRes())

  alloc <- 100L
  cfm <- createFisherMat(alloc)
  testthat::expect_equivalent(cfm, fisherRes())
})

# test_that("all permutations equal to R", {
#   alloc <- 10L
#   outcomes <- as.matrix(expand.grid(0L:alloc, 0L:alloc)) # all permutations. R also has a t(combn(0:10, 2)) & then add cbind(0:10, 0:10)
#
#   applyRFisher <- function(){
#     prow <- function(r) {
#       f <- fisher.test(matrix(c(r[1], alloc - r[1], r[2], alloc - r[2]), 2))
#       return(c(f$p.value, unname(f$estimate), f$conf.int[1], f$conf.int[2]))
#     }
#     r <- t(apply(outcomes,1,prow))
#     return(r)
#   }
#
#   applyMCF <- function() {
#     df <- monteCarloFisher(alloc = alloc, outcomes = outcomes)$col0vs1
#     r <- cbind(df$p, df$or, df$ci_lb, df$ci_ub)
#     return(r)
#   }
#
#   fisherRes <- applyRFisher()
#   mcfRes <- applyMCF()
# })
