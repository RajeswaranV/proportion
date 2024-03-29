% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/323.Sum_Leng_CC_All_Graph.R
\name{PlotlengthCAll}
\alias{PlotlengthCAll}
\title{Plots the sum of length for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)}
\usage{
PlotlengthCAll(n, alp, c, a, b)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{c}{- Continuity correction}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}
}
\description{
Plots the sum of length for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
}
\details{
The  plots for 5 continuity corrected  methods (Wald, Wald-T, Score, Logit-Wald, ArcSine).
Sum of the length is shown as bar graph.
}
\examples{
\dontrun{
n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
PlotlengthCAll(n,alp,c,a,b)
}
}
\seealso{
Other Expected length  of continuity corrected methods: 
\code{\link{PlotexplCAS}()},
\code{\link{PlotexplCAll}()},
\code{\link{PlotexplCLT}()},
\code{\link{PlotexplCSC}()},
\code{\link{PlotexplCTW}()},
\code{\link{PlotexplCWD}()},
\code{\link{PlotlengthCAS}()},
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
