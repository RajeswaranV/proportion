% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/224.CoverageProb_GENERAL_GIVEN_p.R
\name{PlotcovpGEN}
\alias{PlotcovpGEN}
\title{Plot of simulation based Coverage Probability with discrete values for p}
\usage{
PlotcovpGEN(n, LL, UL, alp, hp, t1, t2)
}
\arguments{
\item{n}{- Number of trials}

\item{LL}{- Lower limit}

\item{UL}{- Upper limit}

\item{alp}{- Alpha value (significance level required)}

\item{hp}{- Hypothetical "p"}

\item{t1}{- Lower tolerance limit to check the spread of coverage Probability}

\item{t2}{- Upper tolerance limit to check the spread of coverage Probability}
}
\description{
Plot of simulation based Coverage Probability with discrete values for p
}
\details{
Graphical evaluation of intervals obtained from any method using
coverage probability, root mean square statistic, and the proportion of proportion
lies within the desired level of coverage for the \eqn{n + 1} intervals and pre-defined
space for the parameter \code{p}
}
\examples{
\dontrun{
LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
hp=seq(0,1,by=0.0001)
n= 5; alp=0.05; t1=0.93; t2=0.97
PlotcovpGEN(n,LL,UL,alp,hp,t1,t2)
}
}
\seealso{
Other General methods for coverage probability: 
\code{\link{covpGEN}()}
}
\concept{General methods for coverage probability}
