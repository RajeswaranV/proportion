% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/213.CoverageProb_ADJ_All_Graph.R
\name{PlotcovpAAll}
\alias{PlotcovpAAll}
\title{Plots the Coverage Probability using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)}
\usage{
PlotcovpAAll(n, alp, h, a, b, t1, t2)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{h}{- Adding factor}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}

\item{t1}{- Lower tolerance limit to check the spread of coverage Probability}

\item{t2}{- Upper tolerance limit to check the spread of coverage Probability}
}
\description{
Plots the Coverage Probability using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
}
\details{
The  plots of the Coverage Probability of
6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
}
\examples{
\dontrun{
n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
PlotcovpAAll(n,alp,h,a,b,t1,t2)
}
}
\seealso{
Other Coverage probability of adjusted methods: 
\code{\link{PlotcovpAAS}()},
\code{\link{PlotcovpALR}()},
\code{\link{PlotcovpALT}()},
\code{\link{PlotcovpASC}()},
\code{\link{PlotcovpATW}()},
\code{\link{PlotcovpAWD}()},
\code{\link{covpAAS}()},
\code{\link{covpAAll}()},
\code{\link{covpALR}()},
\code{\link{covpALT}()},
\code{\link{covpASC}()},
\code{\link{covpATW}()},
\code{\link{covpAWD}()}
}
\concept{Coverage probability of adjusted methods}
