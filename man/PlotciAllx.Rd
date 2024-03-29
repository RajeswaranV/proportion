% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/104.ConfidenceIntervals_BASE_n_x_Graph.R
\name{PlotciAllx}
\alias{PlotciAllx}
\title{Plots the CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)}
\usage{
PlotciAllx(x, n, alp)
}
\arguments{
\item{x}{- Number of success}

\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}
}
\description{
Plots the CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
}
\details{
Plots of the Confidence Intervals of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} and \code{x}
}
\examples{
x=5; n=5; alp=0.05;
PlotciAllx(x,n,alp)
}
\seealso{
Other Base methods of CI estimation given x & n: 
\code{\link{PlotciAllxg}()},
\code{\link{PlotciEXx}()},
\code{\link{ciASx}()},
\code{\link{ciAllx}()},
\code{\link{ciBAx}()},
\code{\link{ciEXx}()},
\code{\link{ciLRx}()},
\code{\link{ciLTx}()},
\code{\link{ciSCx}()},
\code{\link{ciTWx}()},
\code{\link{ciWDx}()}
}
\concept{Base methods of CI estimation given x & n}
