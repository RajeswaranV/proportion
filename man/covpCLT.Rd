% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/221.CoverageProb_CC_All.R
\name{covpCLT}
\alias{covpCLT}
\title{Coverage Probability of Continuity corrected Logit Wald method}
\usage{
covpCLT(n, alp, c, a, b, t1, t2)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{c}{- Continiuty correction}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}

\item{t1}{- Lower tolerance limit to check the spread of coverage Probability}

\item{t2}{- Upper tolerance limit to check the spread of coverage Probability}
}
\value{
A dataframe with
 \item{mcpALT}{  Continuity corrected Logit Wald Coverage Probability}
 \item{micpALT }{   Continuity corrected Logit Wald minimum coverage probability}
 \item{RMSE_N }{   Root Mean Square Error from nominal size}
 \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
 \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
 \item{tol }{   Required tolerance for coverage probability}
}
\description{
Coverage Probability of Continuity corrected Logit Wald method
}
\details{
Evaluation of continuity corrected Wald-type interval based on the logit
transformation of p using coverage probability, root mean square statistic,
and the proportion of proportion lies within the desired level of coverage
}
\examples{
n= 10; alp=0.05;c=1/(2*n); a=1;b=1; t1=0.93;t2=0.97
covpCLT(n,alp,c,a,b,t1,t2)
}
\references{
[1] 1998 Agresti A and Coull BA.
Approximate is better than "Exact" for interval estimation of binomial proportions.
The American Statistician: 52; 119 - 126.

[2] 1998 Newcombe RG.
Two-sided confidence intervals for the single proportion: Comparison of seven methods.
Statistics in Medicine: 17; 857 - 872.

[3] 2008 Pires, A.M., Amado, C.
Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
REVSTAT - Statistical Journal, 6, 165-197.
}
\seealso{
Other Coverage probability for continuity corrected methods: 
\code{\link{PlotcovpCAS}()},
\code{\link{PlotcovpCAll}()},
\code{\link{PlotcovpCLT}()},
\code{\link{PlotcovpCSC}()},
\code{\link{PlotcovpCTW}()},
\code{\link{PlotcovpCWD}()},
\code{\link{covpCAS}()},
\code{\link{covpCAll}()},
\code{\link{covpCSC}()},
\code{\link{covpCTW}()},
\code{\link{covpCWD}()}
}
\concept{Coverage probability for continuity corrected methods}
