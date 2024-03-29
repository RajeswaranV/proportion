% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/102.Confidence_base_n_Graph.R
\name{PlotciBA}
\alias{PlotciBA}
\title{Plots the CI estimation of Bayesian method}
\usage{
PlotciBA(n, alp, a, b)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{a}{- Shape parameter 1 for prior Beta distribution in Bayesian model. Can also be a vector of priors.}

\item{b}{- Shape parameter 2 for prior Beta distribution in Bayesian model. Can also be a vector of priors.}
}
\description{
Plots the CI estimation of Bayesian method
}
\details{
The plot of Confidence Interval of \code{n} given \code{alp} using Bayesian method
}
\examples{
n=5; alp=0.05; a=0.5;b=0.5;
PlotciBA(n,alp,a,b)
n=5; alp=0.05; a=c(0.5,2,1,1,2,0.5);b=c(0.5,2,1,1,2,0.5)
PlotciBA(n,alp,a,b)
}
\seealso{
Other Basic methods of CI estimation: 
\code{\link{PlotciAS}()},
\code{\link{PlotciAllg}()},
\code{\link{PlotciAll}()},
\code{\link{PlotciEX}()},
\code{\link{PlotciLR}()},
\code{\link{PlotciLT}()},
\code{\link{PlotciSC}()},
\code{\link{PlotciTW}()},
\code{\link{PlotciWD}()},
\code{\link{ciAS}()},
\code{\link{ciAll}()},
\code{\link{ciBA}()},
\code{\link{ciEX}()},
\code{\link{ciLR}()},
\code{\link{ciLT}()},
\code{\link{ciSC}()},
\code{\link{ciTW}()},
\code{\link{ciWD}()}
}
\concept{Basic methods of CI estimation}
