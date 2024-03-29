% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/321.Expec_Leng_CC_All.R
\name{lengthCWD}
\alias{lengthCWD}
\title{Expected Length summary of continuity corrected Wald method}
\usage{
lengthCWD(n, alp, c, a, b)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{c}{- Continuity correction}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}
}
\value{
A dataframe with
 \item{sumLen}{  The sum of the expected length}
 \item{explMean}{  The mean of the expected length}
 \item{explSD}{  The Standard Deviation of the expected length}
 \item{explMax}{  The max of the expected length}
 \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
 \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
}
\description{
Expected Length summary of continuity corrected Wald method
}
\details{
Evaluation of Wald-type interval with continuity correction using sum of length of the \eqn{n + 1} intervals
}
\examples{
n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
lengthCWD(n,alp,c,a,b)
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
Other Expected length  of continuity corrected methods: 
\code{\link{PlotexplCAS}()},
\code{\link{PlotexplCAll}()},
\code{\link{PlotexplCLT}()},
\code{\link{PlotexplCSC}()},
\code{\link{PlotexplCTW}()},
\code{\link{PlotexplCWD}()},
\code{\link{PlotlengthCAS}()},
\code{\link{PlotlengthCAll}()},
\code{\link{PlotlengthCLT}()},
\code{\link{PlotlengthCSC}()},
\code{\link{PlotlengthCTW}()},
\code{\link{PlotlengthCWD}()},
\code{\link{lengthCAS}()},
\code{\link{lengthCAll}()},
\code{\link{lengthCLT}()},
\code{\link{lengthCSC}()},
\code{\link{lengthCTW}()}
}
\concept{Expected length  of continuity corrected methods}
