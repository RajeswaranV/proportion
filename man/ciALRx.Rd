% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/113.ConfidenceIntervals_ADJ_n_x.R
\name{ciALRx}
\alias{ciALRx}
\title{AdjustedLikelyhood Ratio method of CI estimation}
\usage{
ciALRx(x, n, alp, h)
}
\arguments{
\item{x}{- Number of successes}

\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{h}{- Adding factor}
}
\value{
A dataframe with
 \item{x}{  Number of successes (positive samples)}
 \item{LALRx }{   Likelyhood Ratio Lower limit}
 \item{UALRx }{   Likelyhood Ratio Upper Limit}
 \item{LABB }{   Likelyhood Ratio Lower Abberation}
 \item{UABB }{   Likelyhood Ratio Upper Abberation}
 \item{ZWI }{   Zero Width Interval}
}
\description{
AdjustedLikelyhood Ratio method of CI estimation
}
\details{
Likelihood ratio limits for the data \eqn{x + h} and \eqn{n + (2*h)}
instead of the given \code{x} and \code{n}, where \code{h} is a positive integer
\eqn{(1, 2.)} and for the given \code{x} and \code{n}.
}
\examples{
x=5; n=5; alp=0.05;h=2
ciALRx(x,n,alp,h)
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
\code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
   \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
   \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
   \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
   \code{\link[prevalence]{propCI}} which provides similar functionality.

Other Adjusted methods of CI estimation  given x & n: 
\code{\link{PlotciAAllx}()},
\code{\link{ciAASx}()},
\code{\link{ciAAllx}()},
\code{\link{ciALTx}()},
\code{\link{ciASCx}()},
\code{\link{ciATWx}()},
\code{\link{ciAWDx}()}
}
\concept{Adjusted methods of CI estimation  given x & n}
