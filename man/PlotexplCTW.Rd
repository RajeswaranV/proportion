% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/322.Expec_Leng_CC_All_Graph.R
\name{PlotexplCTW}
\alias{PlotexplCTW}
\title{Plots the Expected length using continuity corrected Wald-T method}
\usage{
PlotexplCTW(n, alp, c, a, b)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{c}{- Continuity correction}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}
}
\description{
Plots the Expected length using continuity corrected Wald-T method
}
\details{
The  plots of the expected length of continuity corrected Wald-T method
}
\examples{
\dontrun{
n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
PlotexplCTW(n,alp,c,a,b)
}
}
\seealso{
Other Expected length  of continuity corrected methods: 
\code{\link{PlotexplCAS}()},
\code{\link{PlotexplCAll}()},
\code{\link{PlotexplCLT}()},
\code{\link{PlotexplCSC}()},
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
\code{\link{lengthCTW}()},
\code{\link{lengthCWD}()}
}
\concept{Expected length  of continuity corrected methods}
