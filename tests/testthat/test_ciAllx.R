
# test_dir("C:/Users/kanchana/Desktop/02-R-package/proportion/tests/testthat/")
source("./R/103.ConfidenceIntervals_BASE_n_x.R")

context("Testing summary of base CI functions")
test_that("Summary base CI functions", {

  Output.df <- ciAllx(5, 5, 0.05)
  Output10.df <- ciAllx(10, 10, 0.05)
  Output50.df <- ciAllx(50, 50, 0.05)
  h=sort(as.character(unique(Output.df$method)))

  expect_that( Output.df, is_a("data.frame") )
  expect_that( nrow(Output.df),   equals(6) )
  expect_that( nrow(Output10.df), equals(6) )
  expect_that( nrow(Output50.df), equals(6) )

  expect_that( ncol(Output.df),   equals(7) )
  expect_that( ncol(Output10.df), equals(7) )
  expect_that( ncol(Output50.df), equals(7) )

  expect_match(h[1], "ArcSine", ignore.case = TRUE)
  expect_match(h[2], "Likelihood", ignore.case = TRUE)
  expect_match(h[3], "Logit-Wald", ignore.case = TRUE)
  expect_match(h[4], "Score", ignore.case = TRUE)
  expect_match(h[5], "Wald", ignore.case = TRUE)
  expect_match(h[6], "Wald-T", ignore.case = TRUE)

  expect_equal(sum(as.numeric(Output.df$LowerAbb)>1),0)
  expect_equal(sum(as.numeric(Output.df$UpperAbb)>1),1)
  expect_equal(sum(as.numeric(Output.df$ZWI)<2),1)

  expect_equal(sum(as.numeric(Output10.df$LowerAbb)>1),0)
  expect_equal(sum(as.numeric(Output10.df$UpperAbb)>1),1)
  expect_equal(sum(as.numeric(Output10.df$ZWI)<2),1)

  expect_equal(sum(as.numeric(Output50.df$LowerAbb)>1),0)
  expect_equal(sum(as.numeric(Output50.df$UpperAbb)>1),2)
  expect_equal(sum(as.numeric(Output50.df$ZWI)<2),2)

  expect_that( ciAllx(0, 0, 0.5), throws_error("greater") )
  expect_that( ciAllx(-1, 4,0.5), throws_error("between") )
  expect_that( ciAllx(1, 1, 1.5), throws_error("between") )
  expect_that( ciAllx(50, 14, .05), throws_error("between") )
  expect_that( ciAllx("asdf",5, 0.5), throws_error("positive") )
  expect_that( ciAllx(10, 10, c(0.5,5, 0.05)), throws_error("between") )

})
