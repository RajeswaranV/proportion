##### 1.ADJUSTED WALD-Coverage Probability
gcovpAWD<-function(n,alp,h,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pAW=0
qAW=0
seAW=0
LAW=0
UAW=0
s=5000								#Simulation run to generate hypothetical p
cpAW=matrix(0,k,s)
ctAW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAW=0								#Coverage probabilty
ctr=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pAW[i]=y[i]/n1
qAW[i]=1-pAW[i]
seAW[i]=sqrt(pAW[i]*qAW[i]/n1)
LAW[i]=pAW[i]-(cv*seAW[i])
UAW[i]=pAW[i]+(cv*seAW[i])
if(LAW[i]<0) LAW[i]=0
if(UAW[i]>1) UAW[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAW[i] && hp[j] < UAW[i])
{
cpAW[i,j]=dbinom(i-1, n,hp[j])
ctAW[i,j]=1
}
}
cppAW[j]=sum(cpAW[,j])
if(t1<cppAW[j]&&cppAW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPAW=data.frame(hp,cp=cppAW,method="Adj-Wald")
return(CPAW)
}

##### 2.ADJUSTED SCORE - Coverage Probability
gcovpASC<-function(n,alp,h,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pAS=0
qAS=0
seAS=0
LAS=0
UAS=0
s=1000								#Simulation run to generate hypothetical p
cpAS=matrix(0,k,s)
ctAS=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAS=0							#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n1)
cv2=(cv/(2*n1))^2

#SCORE (WILSON) METHOD
for(i in 1:k)
{
pAS[i]=y[i]/n1
qAS[i]=1-pAS[i]
seAS[i]=sqrt((pAS[i]*qAS[i]/n1)+cv2)
LAS[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)-(cv*seAS[i]))
UAS[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)+(cv*seAS[i]))
if(LAS[i]<0) LAS[i]=0
if(UAS[i]>1) UAS[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAS[i] && hp[j] < UAS[i])
{
cpAS[i,j]=dbinom(i-1, n,hp[j])
ctAS[i,j]=1
}
}
cppAS[j]=sum(cpAS[,j])						#Coverage Probability
if(t1<cppAS[j]&&cppAS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPAS=data.frame(hp,cp=cppAS,method="Adj-Wilson")
return(CPAS)
}


##### 3.ADJUSTED ARC SINE - Coverage Probability
gcovpAAS<-function(n,alp,h,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pAA=0
qAA=0
seAA=0
LAA=0
UAA=0
s=5000								#Simulation run to generate hypothetical p
cpAA=matrix(0,k,s)
ctAA=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAA=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ADJUSTED ARC-SINE METHOD
for(i in 1:k)
{
pAA[i]=y[i]/n1
qAA[i]=1-pAA[i]
seAA[i]=cv/sqrt(4*n1)
LAA[i]=(sin(asin(sqrt(pAA[i]))-seAA[i]))^2
UAA[i]=(sin(asin(sqrt(pAA[i]))+seAA[i]))^2
if(LAA[i]<0) LAA[i]=0
if(UAA[i]>1) UAA[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAA[i] && hp[j] < UAA[i])
{
cpAA[i,j]=dbinom(i-1, n,hp[j])
ctAA[i,j]=1
}
}
cppAA[j]=sum(cpAA[,j])						#Coverage Probability
if(t1<cppAA[j]&&cppAA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPAA=data.frame(hp,cp=cppAA,method="Adj-ArcSine")
return(CPAA)
}

##### 4.ADJUSTED LOGIT-WALD - Coverage Probability
gcovpALT<-function(n,alp,h,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)

####INITIALIZATIONS
pALT=0
qALT=0
seALT=0
lgit=0
LALT=0
UALT=0
s=5000								#Simulation run to generate hypothetical p
cpALT=matrix(0,k,s)
ctALT=matrix(0,k,s)							#Cover Pbty quantity in sum
cppALT=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ADJUSTED LOGIT-WALD METHOD
for(i in 1:k)
{
pALT[i]=y[i]/n1
qALT[i]=1-pALT[i]
lgit[i]=log(pALT[i]/qALT[i])
seALT[i]=sqrt(pALT[i]*qALT[i]*n1)
LALT[i]=1/(1+exp(-lgit[i]+(cv/seALT[i])))
UALT[i]=1/(1+exp(-lgit[i]-(cv/seALT[i])))
if(LALT[i]<0) LALT[i]=0
if(UALT[i]>1) UALT[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LALT[i] && hp[j] < UALT[i])
{
cpALT[i,j]=dbinom(i-1, n,hp[j])
ctALT[i,j]=1
}
}
cppALT[j]=sum(cpALT[,j])				#Coverage Probability
if(t1<cppALT[j]&&cppALT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

}
CPALT=data.frame(hp,cp=cppALT,method="Adj-Logit-Wald")
return(CPALT)
}

##### 5. ADJUSTED WALD_t - Coverage Probability
gcovpATW<-function(n,alp,h,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)

####INITIALIZATIONS
pATW=0
qATW=0
seATW=0
LATW=0
UATW=0
DOF=0
cv=0
s=5000								#Simulation run to generate hypothetical p
cpATW=matrix(0,k,s)
ctATW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppATW=0
ctr=0

#MODIFIED_t-WALD METHOD
for(i in 1:k)
{
pATW[i]=y[i]/n1
qATW[i]=1-pATW[i]
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOF[i]=2*((f1(pATW[i],n1))^2)/f2(pATW[i],n1)
cv[i]=qt(1-(alp/2), df=DOF[i])
seATW[i]=cv[i]*sqrt(f1(pATW[i],n1))
LATW[i]=pATW[i]-(seATW[i])
UATW[i]=pATW[i]+(seATW[i])
if(LATW[i]<0) LATW[i]=0
if(UATW[i]>1) UATW[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LATW[i] && hp[j] < UATW[i])
{
cpATW[i,j]=dbinom(i-1, n,hp[j])
ctATW[i,j]=1
}
}
cppATW[j]=sum(cpATW[,j])						#Coverage Probability
if(t1<cppATW[j]&&cppATW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPATW=data.frame(hp,cp=cppATW,method="Adj-Wald-T")
return(CPATW)
}

##### 6.ADJUSTED LIKELIHOOD RATIO - Coverage Probability
gcovpALR<-function(n,alp,h,a,b,t1,t2)
{
####INPUT n
y=0:n
k=n+1
y1=y+h
n1=n+(2*h)
####INITIALIZATIONS
mle=0
cutoff=0
LAL=0
UAL=0
s=5000								#Simulation run to generate hypothetical p
cpAL=matrix(0,k,s)
ctAL=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAL=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ADJUSTED LIKELIHOOD-RATIO METHOD
for(i in 1:k)
{
likelhd = function(p) dbinom(y1[i],n1,p)
loglik = function(p) dbinom(y1[i],n1,p,log=TRUE)
mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
cutoff[i]=loglik(mle[i])-(cv^2/2)
loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
LAL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
UAL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAL[i] && hp[j] < UAL[i])
{
cpAL[i,j]=dbinom(i-1, n,hp[j])
ctAL[i,j]=1
}
}
cppAL[j]=sum(cpAL[,j])						#Coverage Probability
if(t1<cppAL[j]&&cppAL[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPAL=data.frame(hp,cp=cppAL,method="Adj-Likelihood")
return(CPAL)
}

#############################################################################################################
#' Plots the Coverage Probability using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpAAll(n,alp,h,a,b,t1,t2)
#' }
#' @export
##### 9.All methods - Coverage Probability
PlotcovpAAll<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL


  #### Calling functions and creating df
  df.1    = gcovpAWD(n,alp,h,a,b,t1,t2)
  df.2    = gcovpAAS(n,alp,h,a,b,t1,t2)
  df.3    = gcovpALR(n,alp,h,a,b,t1,t2)
  df.4    = gcovpASC(n,alp,h,a,b,t1,t2)
  df.5    = gcovpALT(n,alp,h,a,b,t1,t2)
  df.6    = gcovpATW(n,alp,h,a,b,t1,t2)

  nndf=  rbind(df.1,df.2,df.3,df.4,df.5,df.6)

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of adjusted methods") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Wald method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpAWD(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpAWD<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  Waldcovp.df    = covpAWD(n,alp,h,a,b,t1,t2)
  nndf    = gcovpAWD(n,alp,h,a,b,t1,t2)

  nndf$mcp=Waldcovp.df$mcpAW
  nndf$micp=Waldcovp.df$micpAW

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted ArcSine method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpAAS(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpAAS<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf    = gcovpAAS(n,alp,h,a,b,t1,t2)
  ArcSinecovp.df = covpAAS(n,alp,h,a,b,t1,t2)

  nndf$mcp=ArcSinecovp.df$mcpAA
  nndf$micp=ArcSinecovp.df$micpAA

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted ArcSine method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Likelihood Ratio method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpALR(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpALR<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf   = gcovpALR(n,alp,h,a,b,t1,t2)
  LRcovp.df      = covpALR(n,alp,h,a,b,t1,t2)

  nndf$mcp=LRcovp.df$mcpAL
  nndf$micp=LRcovp.df$micpAL

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Likelihood Ratio method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Score method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpASC(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpASC<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf   = gcovpASC(n,alp,h,a,b,t1,t2)
  Scorecovp.df   = covpASC(n,alp,h,a,b,t1,t2)

  nndf$mcp=Scorecovp.df$mcpAS
  nndf$micp=Scorecovp.df$micpAS

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Score method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Logistic Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Logistic Wald method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpALT(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpALT<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf    = gcovpALT(n,alp,h,a,b,t1,t2)
  WaldLcovp.df   = covpALT(n,alp,h,a,b,t1,t2)

  nndf$mcp=WaldLcovp.df$mcpALT
  nndf$micp=WaldLcovp.df$micpALT

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Logistic Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots the Coverage Probability using adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  plots of the Coverage Probability of
#' adjusted Wald-T method for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpATW(n,alp,h,a,b,t1,t2)
#' }
#' @export
PlotcovpATW<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  #### Calling functions and creating df
  nndf    = gcovpATW(n,alp,h,a,b,t1,t2)
  AdWaldcovp.df  = covpATW(n,alp,h,a,b,t1,t2)

  nndf$mcp=AdWaldcovp.df$mcpATW
  nndf$micp=AdWaldcovp.df$micpATW

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of the adjusted Wald-T method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

