% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/312.Expec_Leng_ADJ_All_Graph.R
\name{PlotexplALT}
\alias{PlotexplALT}
\title{Plots the Expected length using adjusted Logit Wald method}
\usage{
PlotexplALT(n, alp, h, a, b)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{h}{- Adding factor}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}
}
\description{
Plots the Expected length using adjusted Logit Wald method
}
\details{
The  plots of the Expected length of adjusted Wald method
}
\examples{
\dontrun{
n= 10; alp=0.05; h=2;a=1;b=1;
PlotexplALT(n,alp,h,a,b)
}
}
\seealso{
Other Expected length  of adjusted methods: 
\code{\link{PlotexplAAS}()},
\code{\link{PlotexplAAll}()},
\code{\link{PlotexplALR}()},
\code{\link{PlotexplASC}()},
\code{\link{PlotexplATW}()},
\code{\link{PlotexplAWD}()},
\code{\link{PlotlengthAAS}()},
\code{\link{PlotlengthAAll}()},
\code{\link{PlotlengthALR}()},
\code{\link{PlotlengthALT}()},
\code{\link{PlotlengthASC}()},
\code{\link{PlotlengthATW}()},
\code{\link{PlotlengthAWD}()},
\code{\link{lengthAAS}()},
\code{\link{lengthAAll}()},
\code{\link{lengthALR}()},
\code{\link{lengthALT}()},
\code{\link{lengthASC}()},
\code{\link{lengthATW}()},
\code{\link{lengthAWD}()}
}
\concept{Expected length  of adjusted methods}
