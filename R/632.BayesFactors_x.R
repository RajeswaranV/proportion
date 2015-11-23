#' Bayesain Hypothesis testing given x: Hypothesis 1: Theta = Theta0 Vs Theta <> Theta0
#' @param x - Number of success
#' @param n - Number of trials from data
#' @param th0 - Hypothetical parameter for H0
#' @param a1 - Priors for hypothesis H1
#' @param b1 - Priors for hypothesis H1
#' @details  Computes Bayes factor under Beta-Binomial model for the
#' model: \eqn{p = p0} Vs \eqn{p \ne p0}{p not equal to p0} from the given number of
#' trials \code{n} and and for given number
#' of successes \eqn{x = 0, 1, 2......n }
#' @return A dataframe with
#'  \item{x}{ Number of successes}
#'  \item{BaFa01}{ Bayesian Factor}
#' @family Hypothesis testing
#' @examples
#' x=682; n=925; th0=0.75; a1=3; b1=3
#' hypotestBAF1x(x,n,th0,a1,b1)
#' @references
#' [1] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods.
#' Springer, New York
#' [2] 2014 Sakthivel S, Subbiah M and Ramakrishnan R
#' Default prior approach for Bayesian testing of hypotheses involving single binomial proportion
#' International Journal of Statistics and Analysis, 4 (2), 139 - 153
#' @export
#####Hypothesis 1: Theta = Theta0  Vs  Theta <> Theta0
hypotestBAF1x<-function(x,n,th0,a1,b1)	#n:No of trials, th0:#Hypothetical parameter,a1,b1,Priors for H1
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(th0)) stop("'th0' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(b1)) stop("'b1' is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1|| x>n || x<0 ) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(th0) != "integer") & (class(th0) != "numeric") || length(th0) >1|| th0<0 ) stop("'th0' has to be greater than 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1) >1|| a1<=0 ) stop("'a1' has to be greater than 0")
  if ((class(b1) != "integer") & (class(b1) != "numeric") || length(b1) >1|| b1<=0 ) stop("'b1' has to be greater than 0")


BaFa01=(beta(a1,b1)/beta(x+a1,n-x+b1))*(th0^x)*((1-th0)^(n-x))

return(data.frame(x,BaFa01))
}
######################################################################################################
#' Bayesain Hypothesis testing given x: Hypothesis 2: Theta = Theta0 Vs Theta > Theta0
#' @param x - Number of success
#' @param n - Number of trials from data
#' @param th0 - Hypothetical parameter for H0
#' @param a1 - Priors for hypothesis H1
#' @param b1 - Priors for hypothesis H1
#' @details  Computes Bayes factor under Beta-Binomial model for the
#' model: \eqn{p = p0} Vs \eqn{p > p0} from the given number of
#' trials \code{n} and and for given number
#' of successes \eqn{x = 0, 1, 2......n }
#' @return A dataframe with
#'  \item{x}{ Number of successes}
#'  \item{BaFa01}{ Bayesian Factor}
#' @family Hypothesis testing
#' @examples
#' x=682; n=925; th0=0.75; a1=3; b1=3
#' hypotestBAF2x(x,n,th0,a1,b1)
#' @references
#' [1] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods.
#' Springer, New York
#' [2] 2014 Sakthivel S, Subbiah M and Ramakrishnan R
#' Default prior approach for Bayesian testing of hypotheses involving single binomial proportion
#' International Journal of Statistics and Analysis, 4 (2), 139 - 153
#' @export
#####Hypothesis 1: Theta = Theta0  Vs  Theta <> Theta0
hypotestBAF2x<-function(x,n,th0,a1,b1)	#n:No of trials, th0:#Hypothetical parameter,a1,b1,Priors for H1
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(th0)) stop("'th0' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(b1)) stop("'b1' is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1|| x>n || x<0 ) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(th0) != "integer") & (class(th0) != "numeric") || length(th0) >1|| th0<0 ) stop("'th0' has to be greater than 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1) >1|| a1<=0 ) stop("'a1' has to be greater than 0")
  if ((class(b1) != "integer") & (class(b1) != "numeric") || length(b1) >1|| b1<=0 ) stop("'b1' has to be greater than 0")


bet1=function(p) dbeta(p,shape1=a1,shape2=b1)
bet2=function(p) dbeta(p,shape1=x+a1,shape2=n-x+b1)

t1=integrate(bet1,th0,1)$value
t2=integrate(bet2,th0,1)$value

BaFa01=(t1/t2)*(th0^x)*((1-th0)^(n-x))

return(data.frame(x,BaFa01))
}
######################################################################################################
#' Bayesain Hypothesis testing given x: Hypothesis 3: Theta = Theta0 Vs Theta < Theta0
#' @param x - Number of success
#' @param n - Number of trials from data
#' @param th0 - Hypothetical parameter for H0
#' @param a1 - Priors for hypothesis H1
#' @param b1 - Priors for hypothesis H1
#' @details  Computes Bayes factor under Beta-Binomial model for the
#' model: \eqn{p = p0} Vs \eqn{p < p0} from the given number of
#' trials \code{n} and and for given number
#' of successes \eqn{x = 0, 1, 2......n }
#' @return A dataframe with
#'  \item{x}{ Number of successes}
#'  \item{BaFa01}{ Bayesian Factor}
#' @family Hypothesis testing
#' @examples
#' x=682; n=925; th0=0.75; a1=3; b1=3
#' hypotestBAF3x(x,n,th0,a1,b1)
#' @references
#' [1] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods.
#' Springer, New York
#' [2] 2014 Sakthivel S, Subbiah M and Ramakrishnan R
#' Default prior approach for Bayesian testing of hypotheses involving single binomial proportion
#' International Journal of Statistics and Analysis, 4 (2), 139 - 153
#' @export
#####Hypothesis 3: Theta = Theta0  Vs  Theta < Theta0
hypotestBAF3x<-function(x,n,th0,a1,b1)	#n:No of trials, th0:#Hypothetical parameter,a1,b1,Priors for H1
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(th0)) stop("'th0' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(b1)) stop("'b1' is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1|| x>n || x<0 ) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(th0) != "integer") & (class(th0) != "numeric") || length(th0) >1|| th0<0 ) stop("'th0' has to be greater than 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1) >1|| a1<=0 ) stop("'a1' has to be greater than 0")
  if ((class(b1) != "integer") & (class(b1) != "numeric") || length(b1) >1|| b1<=0 ) stop("'b1' has to be greater than 0")


bet1=function(p) dbeta(p,shape1=a1,shape2=b1)
bet2=function(p) dbeta(p,shape1=x+a1,shape2=n-x+b1)

t1=integrate(bet1,0,th0)$value
t2=integrate(bet2,0,th0)$value

BaFa01=(t1/t2)*(th0^x)*((1-th0)^(n-x))

return(data.frame(x,BaFa01))
}
######################################################################################################
#' Bayesain Hypothesis testing given x: Hypothesis 4: Theta <= Theta0 Vs Theta > Theta0
#' @param x - Number of success
#' @param n - Number of trials from data
#' @param th0 - Hypothetical parameter for H0
#' @param a0 - Priors for hypothesis H0
#' @param b0 - Priors for hypothesis H0
#' @param a1 - Priors for hypothesis H1
#' @param b1 - Priors for hypothesis H1
#' @details  Computes Bayes factor under Beta-Binomial model for the
#' model: \eqn{p <= p0} Vs \eqn{p > p0} from the given number of trials \code{n} and
#'  for given number
#' of successes \eqn{x = 0, 1, 2......n }
#' @return A dataframe with
#'  \item{x}{ Number of successes}
#'  \item{BaFa01}{ Bayesian Factor}
#' @family Hypothesis testing
#' @examples
#' x=682; n=925; th0=0.75;a0=0.5; b0=0.5; a1=3; b1=3
#' hypotestBAF4x(x,n,th0,a0,b0,a1,b1)
#' @references
#' [1] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods.
#' Springer, New York
#' [2] 2014 Sakthivel S, Subbiah M and Ramakrishnan R
#' Default prior approach for Bayesian testing of hypotheses involving single binomial proportion
#' International Journal of Statistics and Analysis, 4 (2), 139 - 153
#' @export
#####Hypothesis 4: Theta <= Theta0  Vs  Theta > Theta0
hypotestBAF4x<-function(x,n,th0,a0,b0,a1,b1)	#n:No of trials, th0:#Hypothetical parameter,a0,b0,Priors for H0,a1,b1,Priors for H1
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(th0)) stop("'th0' is missing")
  if (missing(a0)) stop("'a0' is missing")
  if (missing(b0)) stop("'b0' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(b1)) stop("'b1' is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1|| x>n || x<0) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(th0) != "integer") & (class(th0) != "numeric") || length(th0) >1|| th0<0 ) stop("'th0' has to be greater than 0")
  if ((class(a0) != "integer") & (class(a0) != "numeric") || length(a0) >1|| a0<=0 ) stop("'a0' has to be greater than 0")
  if ((class(b0) != "integer") & (class(b0) != "numeric") || length(b0) >1|| b0<=0 ) stop("'b0' has to be greater than 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1) >1|| a1<=0 ) stop("'a1' has to be greater than 0")
  if ((class(b1) != "integer") & (class(b1) != "numeric") || length(b1) >1|| b1<=0 ) stop("'b1' has to be greater than 0")

bet0=function(p) dbeta(p,shape1=a0,shape2=b0)
bet01=function(p) dbeta(p,shape1=x+a0,shape2=n-x+b0)	#Null Posterior based

bet1=function(p) dbeta(p,shape1=a1,shape2=b1)
bet11=function(p) dbeta(p,shape1=x+a1,shape2=n-x+b1)	#Alternate Posterior based

t0=integrate(bet0,0,th0)$value
t1=integrate(bet1,th0,1)$value
t01=integrate(bet01,0,th0)$value
t11=integrate(bet11,th0,1)$value

BaFa01=t01*t1/(t0*t11)

return(data.frame(x,BaFa01))
}

######################################################################################################
#' Bayesain Hypothesis testing given x: Hypothesis 5: Theta >= Theta0 Vs Theta < Theta0
#' @param x - Number of success
#' @param n - Number of trials from data
#' @param th0 - Hypothetical parameter for H0
#' @param a0 - Priors for hypothesis H0
#' @param b0 - Priors for hypothesis H0
#' @param a1 - Priors for hypothesis H1
#' @param b1 - Priors for hypothesis H1
#' @details  Computes Bayes factor under Beta-Binomial model for the
#' model: \eqn{p = p0} Vs \eqn{p \ne p0}{p not equal to p0} from the given number of
#' trials \code{n} and and for given number
#' of successes \eqn{x = 0, 1, 2......n }
#' @return A dataframe with
#'  \item{x}{ Number of successes}
#'  \item{BaFa01}{ Bayesian Factor}
#' @family Hypothesis testing
#' @examples
#' x=682; n=925; th0=0.75; a0=0.5; b0= 0.5;a1=1; b1=1
#' hypotestBAF5x(x,n,th0,a0,b0,a1,b1)
#' @references
#' [1] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods.
#' Springer, New York
#' [2] 2014 Sakthivel S, Subbiah M and Ramakrishnan R
#' Default prior approach for Bayesian testing of hypotheses involving single binomial proportion
#' International Journal of Statistics and Analysis, 4 (2), 139 - 153
#' @export
#####Hypothesis 1: Theta = Theta0  Vs  Theta <> Theta0
hypotestBAF5x<-function(x,n,th0,a0,b0,a1,b1)	#n:No of trials, th0:#Hypothetical parameter,a0,b0,Priors for H0,a1,b1,Priors for H1
{
  if (missing(x)) stop("'x' is missing")
  if (missing(th0)) stop("'th0' is missing")
  if (missing(a0)) stop("'a0' is missing")
  if (missing(b0)) stop("'b0' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(b1)) stop("'b1' is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1||x>n || x<0) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(th0) != "integer") & (class(th0) != "numeric") || length(th0) >1|| th0<0 ) stop("'th0' has to be greater than 0")
  if ((class(a0) != "integer") & (class(a0) != "numeric") || length(a0) >1|| a0<=0 ) stop("'a0' has to be greater than 0")
  if ((class(b0) != "integer") & (class(b0) != "numeric") || length(b0) >1|| b0<=0 ) stop("'b0' has to be greater than 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1) >1|| a1<=0 ) stop("'a1' has to be greater than 0")
  if ((class(b1) != "integer") & (class(b1) != "numeric") || length(b1) >1|| b1<=0 ) stop("'b1' has to be greater than 0")


bet0=function(p) dbeta(p,shape1=a0,shape2=b0)
bet01=function(p) dbeta(p,shape1=x+a0,shape2=n-x+b0)	#Null Posterior based

bet1=function(p) dbeta(p,shape1=a1,shape2=b1)
bet11=function(p) dbeta(p,shape1=x+a1,shape2=n-x+b1)	#Alternate Posterior based

t0=integrate(bet0,th0,1)$value
t1=integrate(bet1,0,th0)$value
t01=integrate(bet01,th0,1)$value
t11=integrate(bet11,0,th0)$value

BaFa01=t01*t1/(t0*t11)
return(data.frame(x,BaFa01))

}

######################################################################################################
#' Bayesain Hypothesis testing given x: Hypothesis 6: Theta < Theta1 Vs Theta > Theta2
#' @param x - Number of success
#' @param n - Number of trials from data
#' @param th1 - Hypothetical parameter for H1
#' @param a1 - Priors for hypothesis H1
#' @param b1 - Priors for hypothesis H1
#' @param th2 - Hypothetical parameter for H2
#' @param a2 - Priors for hypothesis H2
#' @param b2 - Priors for hypothesis H2
#' @details  Computes Bayes factor under Beta-Binomial model for the
#' model: \eqn{p < p1} Vs \eqn{p > p2} from the given number of
#' trials \code{n} and and for given number
#' of successes \eqn{x = 0, 1, 2......n }
#' @return A dataframe with
#'  \item{x}{ Number of successes}
#'  \item{BaFa12}{ Bayesian Factor}
#' @family Hypothesis testing
#' @examples
#' x=682; n=925; th1=0.5; a1=1; b1=1; th2=0.9; a2=0.5; b2=0.5
#' hypotestBAF6x(x,n,th1,a1,b1,th2,a2,b2)
#' @references
#' [1] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods.
#' Springer, New York
#' [2] 2014 Sakthivel S, Subbiah M and Ramakrishnan R
#' Default prior approach for Bayesian testing of hypotheses involving single binomial proportion
#' International Journal of Statistics and Analysis, 4 (2), 139 - 153
#' @export
#####Hypothesis 6: Theta < Theta1  Vs  Theta > Theta2
hypotestBAF6x<-function(x,n,th1,a1,b1,th2,a2,b2)	#n:No of trials, th1,th2:#Hypothetical parameter,ai,bi,Priors for hypothesis Hi (i:1,2)
{

  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(th1)) stop("'th1' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(b1)) stop("'b1' is missing")
  if (missing(th2)) stop("'th2' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if (missing(b2)) stop("'b2' is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1||x>n || x<0) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(th1) != "integer") & (class(th1) != "numeric") || length(th1) >1|| th1<0 ) stop("'th1' has to be greater than 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1) >1|| a1<=0 ) stop("'a1' has to be greater than 0")
  if ((class(b1) != "integer") & (class(b1) != "numeric") || length(b1) >1|| b1<=0 ) stop("'b1' has to be greater than 0")
  if ((class(th2) != "integer") & (class(th2) != "numeric") || length(th2) >1|| th2<0 ) stop("'th2' has to be greater than 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2) >1|| a2<=0 ) stop("'a2' has to be greater than 0")
  if ((class(b2) != "integer") & (class(b2) != "numeric") || length(b2) >1|| b2<=0 ) stop("'b2' has to be greater than 0")

  #####For Hypothesis 1
bet1=function(p) dbeta(p,shape1=x+a1,shape2=n-x+b1)

	#####For Hypothesis 2
bet2=function(p) dbeta(p,shape1=x+a2,shape2=n-x+b2)

t1=integrate(bet1,0,th1)$value
t2=integrate(bet2,th2,1)$value

BaFa12=t1/t2

return(data.frame(x,BaFa12))

}
######################################################################################################