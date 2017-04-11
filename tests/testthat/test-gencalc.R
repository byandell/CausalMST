context("CMST 4-way tests")

test_that("generalized code works", {
  
  n_ind <- 500
  z <- sample(0:1, n_ind, replace = TRUE, prob = rep(0.5,2))
  y1 <- z + sqrt(0.4) * rnorm(n_ind)
  y2 <- 0.5 * y1 + sqrt(0.1) * rnorm(n_ind) 

  # Get models
  models <- mediationModels(z, y1, y2, qtl2scan::fit1)
  models <- subset(models$models, names(models$models$LR)[1:4])

  # Shat = covariance matrix
  Shat <- calcShat(models)
  oShat <- CausalMST:::oldShat(models$indLR)
  dimnames(oShat) <- dimnames(Shat)
  expect_equal(Shat, oShat)

  # Z scores for 
  Znew <- CausalMST::calcZ(models, flavor = "B")
  n_ind <- nrow(models$indLR)
  BICs <- CausalMST::calcICs(models, "B")
  Zold <- CausalMST:::oldCalcZ(Shat, BICs, n_ind)
  expect_equal(Znew$Z[1:6], t(Zold)[!is.na(t(Zold))])
  
  parNew <- CausalMST::ParametricIUCMST(models)
  parOld <- CausalMST:::oldParCMST(Zold)
  names(parOld) <- names(parNew)
  expect_equal(parOld, parNew)
})