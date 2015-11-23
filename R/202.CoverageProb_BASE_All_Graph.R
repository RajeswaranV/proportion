#1.WALD
##### 1.WALD-Coverage Probability
gcovpW<-function(n,alp,a,b,t1,t2)
{
###INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pW=0
qW=0
seW=0
LW=0
UW=0
s=5000								#Simulation run to generate hypothetical p
cpW=matrix(0,k,s)
ctW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppW=0								#Coverage probabilty
ctr=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pW[i]=x[i]/n
qW[i]=1-(x[i]/n)
seW[i]=sqrt(pW[i]*qW[i]/n)
LW[i]=pW[i]-(cv*seW[i])
UW[i]=pW[i]+(cv*seW[i])
if(LW[i]<0) LW[i]=0
if(UW[i]>1) UW[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LW[i] && hp[j] < UW[i])
{
cpW[i,j]=dbinom(i-1, n,hp[j])
ctW[i,j]=1
}
}
cppW[j]=sum(cpW[,j])
if(t1<cppW[j]&&cppW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPW=data.frame(hp,cp=cppW,method="Wald")

return(CPW)
}

##### 2.SCORE - Coverage Probability
gcovpS<-function(n,alp,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pS=0
qS=0
seS=0
LS=0
US=0
s=5000								#Simulation run to generate hypothetical p
cpS=matrix(0,k,s)
ctS=matrix(0,k,s)							#Cover Pbty quantity in sum
cppS=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n)
cv2=(cv/(2*n))^2

#SCORE (WILSON) METHOD
for(i in 1:k)
{
pS[i]=x[i]/n
qS[i]=1-(x[i]/n)
seS[i]=sqrt((pS[i]*qS[i]/n)+cv2)
LS[i]=(n/(n+(cv)^2))*((pS[i]+cv1)-(cv*seS[i]))
US[i]=(n/(n+(cv)^2))*((pS[i]+cv1)+(cv*seS[i]))
if(LS[i]<0) LS[i]=0
if(US[i]>1) US[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LS[i] && hp[j] < US[i])
{
cpS[i,j]=dbinom(i-1, n,hp[j])
ctS[i,j]=1
}
}
cppS[j]=sum(cpS[,j])						#Coverage Probability
if(t1<cppS[j]&&cppS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPS=data.frame(hp,cp=cppS,method="Score")
return(CPS)
}

##### 3.ARC SINE - Coverage Probability
gcovpA<-function(n,alp,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pA=0
qA=0
seA=0
LA=0
UA=0
s=5000								#Simulation run to generate hypothetical p
cpA=matrix(0,k,s)
ctA=matrix(0,k,s)							#Cover Pbty quantity in sum
cppA=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pA[i]=x[i]/n
qA[i]=1-pA[i]
seA[i]=cv/sqrt(4*n)
LA[i]=(sin(asin(sqrt(pA[i]))-seA[i]))^2
UA[i]=(sin(asin(sqrt(pA[i]))+seA[i]))^2
if(LA[i]<0) LA[i]=0
if(UA[i]>1) UA[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LA[i] && hp[j] < UA[i])
{
cpA[i,j]=dbinom(i-1, n,hp[j])
ctA[i,j]=1
}
}
cppA[j]=sum(cpA[,j])						#Coverage Probability
if(t1<cppA[j]&&cppA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPA=data.frame(hp,cp=cppA,method="ArcSine")
return(CPA)
}

##### 4.LOGIT-WALD - Coverage Probability
gcovpLT<-function(n,alp,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pLT=0
qLT=0
seLT=0
lgit=0
LLT=0
ULT=0
s=5000								#Simulation run to generate hypothetical p
cpLT=matrix(0,k,s)
ctLT=matrix(0,k,s)							#Cover Pbty quantity in sum
cppLT=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
pLT[1]=0
qLT[1]=1
LLT[1] = 0
ULT[1] = 1-((alp/2)^(1/n))

pLT[k]=1
qLT[k]=0
LLT[k]= (alp/2)^(1/n)
ULT[k]=1

for(j in 1:(k-2))
{
pLT[j+1]=x[j+1]/n
qLT[j+1]=1-pLT[j+1]
lgit[j+1]=log(pLT[j+1]/qLT[j+1])
seLT[j+1]=sqrt(pLT[j+1]*qLT[j+1]*n)
LLT[j+1]=1/(1+exp(-lgit[j+1]+(cv/seLT[j+1])))
ULT[j+1]=1/(1+exp(-lgit[j+1]-(cv/seLT[j+1])))
if(LLT[j+1]<0) LLT[j+1]=0
if(ULT[j+1]>1) ULT[j+1]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LLT[i] && hp[j] < ULT[i])
{
cpLT[i,j]=dbinom(i-1, n,hp[j])
ctLT[i,j]=1
}
}
cppLT[j]=sum(cpLT[,j])				#Coverage Probability
if(t1<cppLT[j]&&cppLT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

}
CPLT=data.frame(hp,cp=cppLT,method="Logit-Wald")
return(CPLT)
}

##### 5. WALD_t - Coverage Probability
gcovpTW<-function(n,alp,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pTW=0
qTW=0
seTW=0
LTW=0
UTW=0
DOF=0
cv=0
s=5000								#Simulation run to generate hypothetical p
cpTW=matrix(0,k,s)
ctTW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppTW=0
ctr=0

#MODIFIED_t-WALD METHOD
for(i in 1:k)
{
if(x[i]==0||x[i]==n)
{
pTW[i]=(x[i]+2)/(n+4)
qTW[i]=1-pTW[i]
}else
{
pTW[i]=x[i]/n
qTW[i]=1-pTW[i]
}
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOF[i]=2*((f1(pTW[i],n))^2)/f2(pTW[i],n)
cv[i]=qt(1-(alp/2), df=DOF[i])
seTW[i]=cv[i]*sqrt(f1(pTW[i],n))
LTW[i]=pTW[i]-(seTW[i])
UTW[i]=pTW[i]+(seTW[i])
if(LTW[i]<0) LTW[i]=0
if(UTW[i]>1) UTW[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LTW[i] && hp[j] < UTW[i])
{
cpTW[i,j]=dbinom(i-1, n,hp[j])
ctTW[i,j]=1
}
}
cppTW[j]=sum(cpTW[,j])						#Coverage Probability
if(t1<cppTW[j]&&cppTW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPTW=data.frame(hp,cp=cppTW,method="Wald-T")
return(CPTW)
}

##### 6.LIKELIHOOD RATIO - Coverage Probability
gcovpL<-function(n,alp,a,b,t1,t2)
{
####INPUT n
y=0:n
k=n+1
####INITIALIZATIONS
mle=0
cutoff=0
LL=0
UL=0
s=5000								#Simulation run to generate hypothetical p
cpL=matrix(0,k,s)
ctL=matrix(0,k,s)							#Cover Pbty quantity in sum
cppL=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LIKELIHOOD-RATIO METHOD
for(i in 1:k)
{
likelhd = function(p) dbinom(y[i],n,p)
loglik = function(p) dbinom(y[i],n,p,log=TRUE)
mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
cutoff[i]=loglik(mle[i])-(cv^2/2)
loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
LL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
UL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LL[i] && hp[j] < UL[i])
{
cpL[i,j]=dbinom(i-1, n,hp[j])
ctL[i,j]=1
}
}
cppL[j]=sum(cpL[,j])						#Coverage Probability
if(t1<cppL[j]&&cppL[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPL=data.frame(hp,cp=cppL,method="Likelihood")
return(CPL)
}

##### 7.Exact - Coverage Probability
gcovpEX=function(n,alp,e,a,b,t1,t2)
{
  nvar=length(e)

  res <- data.frame()

  for(i in 1:nvar)
  {
    lu=gintcovpEX202(n,alp,e[i],a,b,t1,t2)
    res <- rbind(res,lu)
  }
  return(res)
}
gintcovpEX202=function(n,alp,e,a,b,t1,t2)
{

  x=0:n
  k=n+1
  LEX=0
  UEX=0
  s=5000					#Simulation run to generate hypothetical p
  cpEX=matrix(0,k,s)
  ctEX=matrix(0,k,s)			#Cover Pbty quantity in sum
  cppEX=0
  ctr=0
  #EXACT METHOD
  LEX[1]=0
  UEX[1]= 1-((alp/(2*e))^(1/n))
  LEX[k]=(alp/(2*e))^(1/n)
  UEX[k]=1

  for(i in 1:(k-2))
  {
    LEX[i+1]=exlim202l(x[i+1],n,alp,e)
    UEX[i+1]=exlim202u(x[i+1],n,alp,e)
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LEX[i] && hp[j] < UEX[i])
      {
        cpEX[i,j]=dbinom(i-1, n,hp[j])
        ctEX[i,j]=1
      }
    }
    cppEX[j]=sum(cpEX[,j])						#Coverage Probability
    if(t1<cppEX[j]&&cppEX[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPEX=data.frame(hp,cpp=cppEX)
  mcpEX=mean(cppEX)							#Mean Cov Prob
  micpEX=min(cppEX)							#Min Cov Prob
  return(data.frame(CPEX,mcpEX,micpEX,e))
}
#####TO FIND LOWER LIMITS
exlim202l=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
  LEX= uniroot(f1,c(0,1))$root
  return(LEX)
}
#####TO FIND UPPER LIMITS
exlim202u=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f2  = function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
  UEX = uniroot(f2,c(0,1))$root
  return(UEX)
}
#############################################################################################################
#' Graphs of Coverage Probability - exact method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e -  Exact method indicator (1:Clop-Pear,0.5:MID-p)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  graphs of basic Coverage Probability methods
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; e=0.5; a=1;b=1; t1=0.93;t2=0.97 # Mid-p
#' PlotcovpEX(n,alp,e,a,b,t1,t2)
#' n= 10; alp=0.05; e=1; a=1;b=1; t1=0.93;t2=0.97 #Clop-Pear
#' PlotcovpEX(n,alp,e,a,b,t1,t2)
#' n=5; alp=0.05;
#' e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpEX(n,alp,e,a,b,t1,t2)
#' }
#' @export
##### 7. EXACT METHOD - Coverage Probability
PlotcovpEX=function(n,alp,e,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10 ) stop("Plot of only 10 intervals of 'e' is possible")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  if(length(e)>1){
    dfex=gcovpEX(n,alp,e,a,b,t1,t2)
    exdf=dfex[,c(1,2,5)]
    exdf$e=as.factor(exdf$e)

ggplot2::ggplot(exdf, ggplot2::aes(x=hp, y=cpp))+
  ggplot2::labs(y = "Coverage Probability") +
  ggplot2::labs(title = "Coverage Probability for exact method for multiple e values") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
  ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
  ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
  ggplot2::geom_line(ggplot2::aes(color=e)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)
  }
  else{
    dfex=gcovpEX(n,alp,e,a,b,t1,t2)
dfex1=data.frame(micp=dfex$micpEX[1]	,mcp=dfex$mcpEX[1]	)

ggplot2::ggplot(dfex, ggplot2::aes(x=hp, y=cpp))+
  ggplot2::labs(title = "Coverage Probability of exact method") +
  ggplot2::labs(y = "Coverage Probability") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_line(ggplot2::aes(color="Coverage Probability"))+
  ggplot2::geom_point(ggplot2::aes(color="CP Values"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp),color="Confidence Level"),linetype = 2)+
  ggplot2::geom_hline(data=dfex1,ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
  ggplot2::geom_hline(data=dfex1,ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
  ggplot2::scale_colour_manual(name='Heading',
                               values=c('Coverage Probability'='red',
                                        'CP Values'='red',
                                        'Minimum Coverage'='black',
                                        'Mean Coverage'='blue',
                                        'Confidence Level'='brown'),
                               guide='legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(2,1,1,1,1),
                                                                     shape=c(NA, NA, 16,NA,NA),
                                                                     linetype=c(1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1))))


}
}

#############################################################################################################
#' Graphs of  Coverage Probability of the Bayesian method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  The  graphs of  Coverage Probability of Bayesian method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97;a1=1;a2=1
#' PlotcovpBA(n,alp,a,b,t1,t2,a1,a2)
#' }
#' @export
#8.BAYESIAN
PlotcovpBA<-function(n,alp,a,b,t1,t2,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1)>1 || a1<0  ) stop("'a1' has to be greater than or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2)>1 || a2<0  ) stop("'a2' has to be greater than or equal to 0")
  ID=method=Value=hp=cp=cpp=mcpBAQ=micpBAQ=mcpBAH=micpBAH=NULL

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LBAQ=0
UBAQ=0
LBAH=0
UBAH=0
s=5000
cpBAQ=matrix(0,k,s)
ctBAQ=matrix(0,k,s)							#Cover Pbty quantity in sum
cppBAQ=0								#Coverage probabilty
ctr=0

cpBAH=matrix(0,k,s)
ctBAH=matrix(0,k,s)							#Cover Pbty quantity in sum
cppBAH=0								#Coverage probabilty
ctrH=0

##############
#library(TeachingDemos)				#To get HPDs
for(i in 1:k)
{
#Quantile Based Intervals
LBAQ[i]=qbeta(alp/2,x[i]+a1,n-x[i]+a2)
UBAQ[i]=qbeta(1-(alp/2),x[i]+a1,n-x[i]+a2)

LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[1]
UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[2]

}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LBAQ[i] && hp[j] < UBAQ[i])
{
cpBAQ[i,j]=dbinom(i-1, n,hp[j])
ctBAQ[i,j]=1
}
if(hp[j] > LBAH[i] && hp[j] < UBAH[i])
{
cpBAH[i,j]=dbinom(i-1, n,hp[j])
ctBAH[i,j]=1
}

}
cppBAQ[j]=sum(cpBAQ[,j])						#Coverage Probability
if(t1<cppBAQ[j]&&cppBAQ[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

cppBAH[j]=sum(cpBAH[,j])						#Coverage Probability
if(t1<cppBAH[j]&&cppBAH[j]<t2) ctrH=ctrH+1		#tolerance for cov prob - user defined

}
CPBAQ=data.frame(hp,cpp=cppBAQ,method="Quantile")
CPBAH=data.frame(hp,cpp=cppBAH,method="HPD")

df.new=rbind(CPBAQ,CPBAH)
df.new$mcpBAQ=mean(cppBAQ)
df.new$micpBAQ=min(cppBAQ)					#Mean Cov Prob

df.new$mcpBAH=mean(cppBAH)
df.new$micpBAH=min(cppBAH)					#Mean Cov Prob


ggplot2::ggplot(df.new, ggplot2::aes(x=hp, y=cpp))+
  ggplot2::labs(title = "Coverage Probability of Bayesian methods") +
  ggplot2::labs(y = "Coverage Probability") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_line(ggplot2::aes(color=method)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp),color="Confidence Level"),linetype = 2)+
  ggplot2::geom_hline(ggplot2::aes(yintercept=micpBAQ,color="Minimum Coverage Quantile"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=mcpBAQ,color="Mean Coverage Quantile"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=micpBAH,color="Minimum Coverage HPD"))+
  ggplot2::geom_hline(ggplot2::aes(yintercept=mcpBAH,color="Mean Coverage HPD"))+
  ggplot2::scale_colour_manual(name='Heading',
                               values=c('Quantile' ='black',
                                        'HPD' = 'red',
                                        'Minimum Coverage Quantile'='red',
                                        'Mean Coverage Quantile'='blue',
                                        'Minimum Coverage HPD'='black',
                                        'Mean Coverage HPD'='cyan',
                                        'Confidence Level'='brown'),
                               guide='legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(2,1,1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1,1,1),
                                                                     linetype=c(1,1,1,1,1,1,1))))

}

#############################################################################################################
#' Graphs of basic Coverage Probability 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  graphs of basic Coverage Probability methods
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpAll(n,alp,a,b,t1,t2)
#' }
#' @export
##### 9.  Coverage Probability - Graph
PlotcovpAll<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  df1=gcovpW(n,alp,a,b,t1,t2)
  df2=gcovpS(n,alp,a,b,t1,t2)
  df3=gcovpA(n,alp,a,b,t1,t2)
  df4=gcovpLT(n,alp,a,b,t1,t2)
  df5=gcovpTW(n,alp,a,b,t1,t2)
  df6=gcovpL(n,alp,a,b,t1,t2)

 nndf=  rbind(df1,df2,df3,df4,df5,df6)

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(title = "Coverage Probability for 6 base methods") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Wald method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpWD(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpWD<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  Waldcovp.df    = covpWD(n,alp,a,b,t1,t2)

  nndf=gcovpW(n,alp,a,b,t1,t2)
  nndf$mcp=Waldcovp.df$mcpW
  nndf$micp=Waldcovp.df$micpW

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base ArcSine method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpAS(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpAS<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ArcSinecovp.df = covpAS(n,alp,a,b,t1,t2)

  nndf=gcovpA(n,alp,a,b,t1,t2)
  nndf$mcp=ArcSinecovp.df$mcpA
  nndf$micp=ArcSinecovp.df$micpA

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for ArcSine method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Likelihood Ratio method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpLR(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpLR<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  LRcovp.df      = covpLR(n,alp,a,b,t1,t2)
 # ss1 = data.frame(method = LRcovp.df$method, MeanCP=LRcovp.df$mcpL, MinCP= LRcovp.df$micpL, RMSE_N=LRcovp.df$RMSE_N,RMSE_M=LRcovp.df$RMSE_M,RMSE_MI=LRcovp.df$RMSE_MI,tol=LRcovp.df$tol)

  nndf=gcovpL(n,alp,a,b,t1,t2)
  nndf$mcp=LRcovp.df$mcpL
  nndf$micp=LRcovp.df$micpL

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Likelihood Ratio method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Score method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpSC(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpSC<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  Scorecovp.df   = covpSC(n,alp,a,b,t1,t2)
  #ss1 = data.frame( MeanCP=Scorecovp.df$mcpS, MinCP= Scorecovp.df$micpS, RMSE_N=Scorecovp.df$RMSE_N,RMSE_M=Scorecovp.df$RMSE_M,RMSE_MI=Scorecovp.df$RMSE_MI,tol=Scorecovp.df$tol)

  nndf=gcovpS(n,alp,a,b,t1,t2)
  nndf$mcp=Scorecovp.df$mcpS
  nndf$micp=Scorecovp.df$micpS

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Score method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Logit Wald method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpLT(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpLT<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  WaldLcovp.df   = covpLT(n,alp,a,b,t1,t2)

  nndf=gcovpLT(n,alp,a,b,t1,t2)
  nndf$mcp=WaldLcovp.df$mcpLT
  nndf$micp=WaldLcovp.df$micpLT

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Logit Wald method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}

#############################################################################################################
#' Plots Coverage Probability for base Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Plots Coverage Probability for base Wald-T method
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1; b=1; t1=0.93; t2=0.97
#' PlotcovpTW(n,alp,a,b,t1,t2)
#' }
#' @export
PlotcovpTW<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
  ID=method=Value=hp=cp=cpp=mcp=micp=NULL

  ####INPUT n
  AdWaldcovp.df  = covpTW(n,alp,a,b,t1,t2)

  nndf=gcovpTW(n,alp,a,b,t1,t2)
  nndf$mcp=AdWaldcovp.df$mcpTW
  nndf$micp=AdWaldcovp.df$micpTW

  ggplot2::ggplot(nndf, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability for Wald-T method") +
    ggplot2::labs(y = "Coverage Probability") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t1), color="red",linetype = 2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=t2), color="blue",linetype = 2) +
    ggplot2::geom_text(ggplot2::aes(y=t1, label="\nLower tolerance(t1)", x=.1), colour="red", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_text(ggplot2::aes(y=t2, label="Higher tolerance(t2)", x=.1), colour="blue", text=ggplot2::element_text(size=11)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=micp,color="Minimum Coverage"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=mcp,color="Mean Coverage"))+
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::guides(colour = ggplot2::guide_legend("Heading"))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=1-(alp)),linetype = 2)

}
