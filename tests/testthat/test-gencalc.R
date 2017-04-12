context("CMST 4-way tests")

test_that("generalized code works", {
  
  n_ind <- 1000
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
  expect_equal(as.vector(Znew), t(Zold)[!is.na(t(Zold))])

  # Parametric test  
  parNew <- normIUCMST(models)
  parOld <- CausalMST:::oldParCMST(Zold)
  expect_equal(parOld, parNew$pv)
  
  # Joint Parametric test  
  jointNew <- normJointIUCMST(models)
  Cor.hat <- CausalMST:::oldCorHat(oShat)
  jointOld <- CausalMST:::oldJointCMST(Zold, Cor.hat)
  expect_equal(jointOld, jointNew$pv)
  
  # Nonparametric (sign) test  
  nonparNew <- binomIUCMST(models)
  nonparOld <- CausalMST:::oldNonparCMST("bic", nrow(models$indLR), models$df, models$indLR)
  names(parOld) <- names(parNew)
  expect_equal(nonparOld, nonparNew$pv)
})