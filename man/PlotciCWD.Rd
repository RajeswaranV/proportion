% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/122.ConfidenceIntervals_CC_n_Graph.R
\name{PlotciCWD}
\alias{PlotciCWD}
\title{Plots the CI estimation of  continuity corrected Wald method given n, alp and c}
\usage{
PlotciCWD(n, alp, c)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{c}{- Continuity correction}
}
\description{
Plots the CI estimation of  continuity corrected Wald method given n, alp and c
}
\details{
Plots the Confidence Interval for continuity corrected Wald method
 for \code{n} given \code{alp} along with Continuity correction \code{c}
}
\examples{
n=5; alp=0.05;c=1/(2*n)
PlotciCWD(n,alp,c)
}
\seealso{
Other Continuity correction methods of CI estimation: 
\code{\link{PlotciCAS}()},
\code{\link{PlotciCAllg}()},
\code{\link{PlotciCAll}()},
\code{\link{PlotciCLT}()},
\code{\link{PlotciCSC}()},
\code{\link{PlotciCTW}()},
\code{\link{ciCAS}()},
\code{\link{ciCAll}()},
\code{\link{ciCLT}()},
\code{\link{ciCSC}()},
\code{\link{ciCTW}()},
\code{\link{ciCWD}()}
}
\concept{Continuity correction methods of CI estimation}
