% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/612.Empirical_x.R
\name{empericalBAx}
\alias{empericalBAx}
\title{The empirical Bayesian approach for Beta-Binomial model given x}
\usage{
empericalBAx(x, n, alp, sL, sU)
}
\arguments{
\item{x}{- Number of successes}

\item{n}{- Number of trials}

\item{alp}{- Alpha value (significance level required)}

\item{sL}{- Lower support for MLE stats::optimization}

\item{sU}{- Upper support for MLE stats::optimization}
}
\value{
A dataframe with
 \item{x }{- Number of successes (positive samples)}
 \item{pomean }{ - Posterior mean}
 \item{LEBAQ }{ - Lower limits of Quantile based intervals}
 \item{UEBAQ }{ - Upper limits of Quantile based intervals}
 \item{LEBAH }{ - Lower limits of HPD intervals}
 \item{UEBAH }{ - Upper limits of HPD intervals}
}
\description{
The empirical Bayesian approach for Beta-Binomial model given x
}
\details{
Highest Probability Density (HPD) and two tailed intervals are provided for the
required x (any one value from \eqn{ 0, 1, 2 ..n}) based on empirical Bayesian approach for
Beta-Binomial model. Lower and Upper support values are needed to obtain the MLE of
marginal likelihood for prior parameters.
}
\examples{
sL=runif(1,0,2)				#Lower and upper of Support for MLE stats::optimization
sU=runif(1,sL,10)
x=0; n= 5; alp=0.05
empericalBAx(x,n,alp,sL,sU)
}
\references{
[1] 1998 Lehmann EL and Casella G
Theory of Point Estimation, 2nd ed Springer, New York
}
\seealso{
Other Miscellaneous  functions for Bayesian method: 
\code{\link{empericalBA}()},
\code{\link{probPOSx}()},
\code{\link{probPOS}()},
\code{\link{probPREx}()},
\code{\link{probPRE}()}
}
\concept{Miscellaneous  functions for Bayesian method}
