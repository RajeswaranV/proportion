% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/225.CoverageProb_GENERAL_SIMULATEDp.R
\name{PlotcovpSIM}
\alias{PlotcovpSIM}
\title{Plots graph of simulation based Coverage Probability with simulated p in [0,1]}
\usage{
PlotcovpSIM(n, LL, UL, alp, s, a, b, t1, t2)
}
\arguments{
\item{n}{- Number of trials}

\item{LL}{- Lower limit}

\item{UL}{- Upper limit}

\item{alp}{- Alpha value (significance level required)}

\item{s}{- Number of hypothetical "p"}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}

\item{t1}{- Lower tolerance limit to check the spread of coverage Probability}

\item{t2}{- Upper tolerance limit to check the spread of coverage Probability}
}
\description{
Plots graph of simulation based Coverage Probability with simulated p in [0,1]
}
\details{
Graphical evaluation of intervals obtained from any method using
coverage probability, root mean square statistic, and the proportion of
proportion lies within the desired level of coverage for the \eqn{n + 1} intervals
and pre-defined space for the parameter \code{p} using Monte Carle simulation
}
\examples{
\dontrun{
LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
n= 5; alp=0.05; s=5000; a=1; b=1; t1=0.93; t2=0.97
PlotcovpSIM(n,LL,UL,alp,s,a,b,t1,t2)
}
}
\seealso{
Other Simulated methods for coverage probability: 
\code{\link{covpSIM}()}
}
\concept{Simulated methods for coverage probability}
