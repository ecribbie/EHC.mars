# testfwd_stepwise.R
library(EHC.MARS)
load("testfwd_stepwise.RData")
testy=marstestdata[,1]
testx=marstestdata[,2:11]
tmp_fwd=fwd_stepwise(testy,testx,testmc)
test_that("fwd_stepwise() returns the correct object", {
  expect_equal(fwd_stepwise(testy,testx,testmc), testfwd)
})
