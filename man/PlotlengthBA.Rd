% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/303.Sum_Leng_BASE_All_Graph.R
\name{PlotlengthBA}
\alias{PlotlengthBA}
\title{Plot the Bayesian method of length summary}
\usage{
PlotlengthBA(n, alp, a, b, a1, a2)
}
\arguments{
\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{a}{- Beta parameters for hypo "p"}

\item{b}{- Beta parameters for hypo "p"}

\item{a1}{- Beta Prior Parameters for Bayesian estimation}

\item{a2}{- Beta Prior Parameters for Bayesian estimation}
}
\value{
A dataframe with
}
\description{
Plot the Bayesian method of length summary
}
\details{
Plots of Bayesian Highest Probability Density (HPD) and two tailed
intervals using expected length of the \eqn{n + 1}
intervals for the Beta - Binomial conjugate prior model for the probability of success \code{p}
}
\examples{
n=5; alp=0.05;a=1;b=1;a1=1;a2=1
PlotlengthBA(n,alp,a,b,a1,a2)
}
\seealso{
Other Expected length  of base methods: 
\code{\link{PlotexplAS}()},
\code{\link{PlotexplAll}()},
\code{\link{PlotexplBA}()},
\code{\link{PlotexplEX}()},
\code{\link{PlotexplLR}()},
\code{\link{PlotexplLT}()},
\code{\link{PlotexplSC}()},
\code{\link{PlotexplTW}()},
\code{\link{PlotexplWD}()},
\code{\link{PlotlengthAS}()},
\code{\link{PlotlengthAll}()},
\code{\link{PlotlengthEX}()},
\code{\link{PlotlengthLR}()},
\code{\link{PlotlengthLT}()},
\code{\link{PlotlengthSC}()},
\code{\link{PlotlengthTW}()},
\code{\link{PlotlengthWD}()},
\code{\link{lengthAS}()},
\code{\link{lengthAll}()},
\code{\link{lengthBA}()},
\code{\link{lengthEX}()},
\code{\link{lengthLR}()},
\code{\link{lengthLT}()},
\code{\link{lengthSC}()},
\code{\link{lengthTW}()},
\code{\link{lengthWD}()}
}
\concept{Expected length  of base methods}
