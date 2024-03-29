% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/402.p-Confidence_p-Bias_BASE_All_Graph.R
\name{PlotpCOpBIBA}
\alias{PlotpCOpBIBA}
\title{Plots p-confidence and p-bias for Bayesian method given n and alpha level and priors a & b}
\usage{
PlotpCOpBIBA(n, alp, a1, a2)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{a1}{- Shape parameter 1 for prior Beta distribution in Bayesian model}

\item{a2}{- Shape parameter 2 for prior Beta distribution in Bayesian model}
}
\description{
Plots p-confidence and p-bias for Bayesian method given n and alpha level and priors a & b
}
\details{
Evaluation of Bayesian Highest Probability Density (HPD) and two tailed
intervals using p-confidence and p-bias for the \eqn{n + 1} intervals for the
Beta - Binomial conjugate prior model for the probability of success \code{p}
}
\examples{
n=5; alp=0.05;a1=1;a2=1
PlotpCOpBIBA(n,alp,a1,a2)
}
\references{
[1] 2005 Vos PW and Hudson S.
Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
The American Statistician: 59; 137 - 142.
}
\seealso{
Other p-confidence and p-bias of base methods: 
\code{\link{PlotpCOpBIAS}()},
\code{\link{PlotpCOpBIAll}()},
\code{\link{PlotpCOpBIEX}()},
\code{\link{PlotpCOpBILR}()},
\code{\link{PlotpCOpBILT}()},
\code{\link{PlotpCOpBISC}()},
\code{\link{PlotpCOpBITW}()},
\code{\link{PlotpCOpBIWD}()},
\code{\link{pCOpBIAS}()},
\code{\link{pCOpBIAll}()},
\code{\link{pCOpBIBA}()},
\code{\link{pCOpBIEX}()},
\code{\link{pCOpBILR}()},
\code{\link{pCOpBILT}()},
\code{\link{pCOpBISC}()},
\code{\link{pCOpBITW}()},
\code{\link{pCOpBIWD}()}
}
\concept{p-confidence and p-bias of base methods}
