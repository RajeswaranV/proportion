% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/422.p-Confidence_p-Bias_CC_All_Graph.R
\name{PlotpCOpBICAS}
\alias{PlotpCOpBICAS}
\title{Plots  p-confidence and p-bias for continuity corrected ArcSine method}
\usage{
PlotpCOpBICAS(n, alp, c)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{c}{- Continuity correction}
}
\description{
Plots  p-confidence and p-bias for continuity corrected ArcSine method
}
\details{
p-confidence and p-bias plots for continuity corrected ArcSine method
}
\examples{
n=5; alp=0.05;c=1/(2*n)
PlotpCOpBICAS(n,alp,c)
}
\seealso{
Other p-confidence and p-bias of continuity corrected methods: 
\code{\link{PlotpCOpBICAll}()},
\code{\link{PlotpCOpBICLT}()},
\code{\link{PlotpCOpBICSC}()},
\code{\link{PlotpCOpBICTW}()},
\code{\link{PlotpCOpBICWD}()},
\code{\link{pCOpBICAS}()},
\code{\link{pCOpBICAll}()},
\code{\link{pCOpBICLT}()},
\code{\link{pCOpBICSC}()},
\code{\link{pCOpBICTW}()},
\code{\link{pCOpBICWD}()}
}
\concept{p-confidence and p-bias of continuity corrected methods}
