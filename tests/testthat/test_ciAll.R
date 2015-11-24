
# test_dir("C:/Users/kanchana/Desktop/02-R-package/proportion/tests/testthat/")
source("C:/Users/kanchana/Desktop/02-R-package/proportion/R/101.Confidence_base_n.R")

context("Testing summary of base CI functions")
test_that("Summary base CI functions", {

  Output.df <- ciAll(5, 0.05)
  Output10.df <- ciAll(10, 0.05)
  Output50.df <- ciAll(50, 0.05)
  h=sort(as.character(unique(Output.df$method)))

  expect_that( Output.df, is_a("data.frame") )
  expect_that( nrow(Output.df),   equals(36) )
  expect_that( nrow(Output10.df), equals(66) )
  expect_that( nrow(Output50.df), equals(306) )

  expect_that( ncol(Output.df),   equals(7) )
  expect_that( ncol(Output10.df), equals(7) )
  expect_that( ncol(Output50.df), equals(7) )

  expect_match(h[1], "ArcSine", ignore.case = TRUE)
  expect_match(h[2], "Likelihood", ignore.case = TRUE)
  expect_match(h[3], "Logit-Wald", ignore.case = TRUE)
  expect_match(h[4], "Score", ignore.case = TRUE)
  expect_match(h[5], "Wald", ignore.case = TRUE)
  expect_match(h[6], "Wald-T", ignore.case = TRUE)

  expect_equal(sum(as.numeric(Output.df$LowerAbb)>1),5)
  expect_equal(sum(as.numeric(Output.df$UpperAbb)>1),5)
  expect_equal(sum(as.numeric(Output.df$ZWI)>1),3)

  expect_that( ciAll(0, 0.5), throws_error("greater") )
  expect_that( ciAll(-1, 0.5), throws_error("greater") )
  expect_that( ciAll(1, 1.5), throws_error("between") )
  expect_that( ciAll("asdf", 0.5), throws_error("greater") )
  expect_that( ciAll(10, c(0.5,0.05)), throws_error("between") )

})
