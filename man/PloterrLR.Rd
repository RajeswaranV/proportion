% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/502.Error-Failure_LimitBased_BASE_All_Graph.R
\name{PloterrLR}
\alias{PloterrLR}
\title{Plots error, long term power and pass/fail criteria using Likelihood Ratio method}
\usage{
PloterrLR(n, alp, phi, f)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{phi}{- Null hypothesis value}

\item{f}{- Failure criterion}
}
\description{
Plots error, long term power and pass/fail criteria using Likelihood Ratio method
}
\details{
Plots of error, long term power and pass/fail
criteria using Likelihood Ratio method
}
\examples{
n=20; alp=0.05; phi=0.05; f=-2
PloterrLR(n,alp,phi,f)
}
\seealso{
Other Error for base methods: 
\code{\link{PloterrAS}()},
\code{\link{PloterrAll}()},
\code{\link{PloterrBA}()},
\code{\link{PloterrEX}()},
\code{\link{PloterrLT}()},
\code{\link{PloterrSC}()},
\code{\link{PloterrTW}()},
\code{\link{PloterrWD}()},
\code{\link{errAS}()},
\code{\link{errAll}()},
\code{\link{errBA}()},
\code{\link{errEX}()},
\code{\link{errLR}()},
\code{\link{errLT}()},
\code{\link{errSC}()},
\code{\link{errTW}()},
\code{\link{errWD}()}
}
\concept{Error for base methods}
