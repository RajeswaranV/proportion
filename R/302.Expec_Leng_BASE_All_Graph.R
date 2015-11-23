##### 1.WALD Expected Length for a given n and alpha level
gexplWD<-function(n,alp,a,b) #n:No of trials,alp:sign level,a&b beta parameters for hypo "p'
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pW=0
qW=0
seW=0
LW=0
UW=0
s=5000
LEW=0 								#LENGTH OF INTERVAL

ewiW=matrix(0,k,s)						#Expected length quantity in sum
ewW=0									#Expected Length
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
LEW[i]=UW[i]-LW[i]
}
sumLEW=sum(LEW)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiW[i,j]=LEW[i]*dbinom(i-1, n,hp[j])
}
ewW[j]=sum(ewiW[,j])						#Expected Length
}
ELW=data.frame(hp,ew=ewW,method="Wald")
return(ELW)
}

##### 2.SCORE - Expected Length for a given n and alpha level
gexplSC<-function(n,alp,a,b)
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
s=5000
LES=0 								#LENGTH OF INTERVAL
ewiS=matrix(0,k,s)						#Expected length quantity in sum
ewS=0									#Expected Length
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
LES[i]=US[i]-LS[i]
}
sumLES=sum(LES)

####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiS[i,j]=LES[i]*dbinom(i-1, n,hp[j])
}
ewS[j]=sum(ewiS[,j])						#Expected Length
}
ELS=data.frame(hp,ew=ewS,method="Wilson")
# windows()
# plot(ELS,xlab="p",ylab="Expected Length",main="Wilson",type="l")
# abline(v=0.5, lty=2)
return(ELS)
}

##### 3. ARC SINE - Expected Length for a given n and alpha level
gexplAS<-function(n,alp,a,b)
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
s=5000
LEA=0 								#LENGTH OF INTERVAL

ewiA=matrix(0,k,s)						#Expected length quantity in sum
ewA=0									#Expected Length
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
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
LEA[i]=UA[i]-LA[i]
}
sumLEA=sum(LEA)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiA[i,j]=LEA[i]*dbinom(i-1, n,hp[j])
}
ewA[j]=sum(ewiA[,j])						#Expected Length
}
ELA=data.frame(hp,ew=ewA,method="ArcSine")
# windows()
# plot(ELA,xlab="p",ylab="Expected Length",main="Arc Sine",type="l")
# abline(v=0.5, lty=2)
return(ELA)
}


##### 4.LOGIT-WALD - Expected Length for a given n and alpha level
gexplLT<-function(n,alp,a,b) #n:No of trials,alp:sign level,a&b beta parameters for hypo "p'
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
s=5000
LELT=0 								#LENGTH OF INTERVAL

ewiLT=matrix(0,k,s)						#Expected length quantity in sum
ewLT=0									#Expected Length
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
for(i in 1:k)
{
LELT[i]=ULT[i]-LLT[i]
}
sumLET=sum(LELT)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiLT[i,j]=LELT[i]*dbinom(i-1, n,hp[j])
}
ewLT[j]=sum(ewiLT[,j])						#Expected Length
}
ELLT=data.frame(hp,ew=ewLT,method="Logit-Wald")
# windows()
# plot(ELLT,xlab="p",ylab="Expected Length",main="Logit-Wald",type="l")
# abline(v=0.5, lty=2)
return(ELLT)
}

##### 5.t-WALD - Expected Length for a given n and alpha level
gexplTW<-function(n,alp,a,b)
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
s=5000
LETW=0 								#LENGTH OF INTERVAL

ewiTW=matrix(0,k,s)						#Expected length quantity in sum
ewTW=0									#Expected Length
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
LETW[i]=UTW[i]-LTW[i]
}
sumLETW=sum(LETW)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiTW[i,j]=LETW[i]*dbinom(i-1, n,hp[j])
}
ewTW[j]=sum(ewiTW[,j])						#Expected Length
}
ELTW=data.frame(hp,ew=ewTW,method="Wald-T")
# windows()
# plot(ELTW,xlab="p",ylab="Expected Length",main="t-Wald",type="l")
# abline(v=0.5, lty=2)
return(ELTW)
}

#####6.LIKELIHOOD RATIO - Expected Length for a given n and alpha level
gexplLR<-function(n,alp,a,b)
{
####INPUT n
y=0:n
k=n+1
####INITIALIZATIONS
mle=0
cutoff=0
LL=0
UL=0
s=5000
LEL=0 								#LENGTH OF INTERVAL

ewiL=matrix(0,k,s)						#Expected length quantity in sum
ewL=0										#Simulation run to generate hypothetical p
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
LEL[i]=UL[i]-LL[i]
}
sumLEL=sum(LEL)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiL[i,j]=LEL[i]*dbinom(i-1, n,hp[j])
}
ewL[j]=sum(ewiL[,j])						#Expected Length
}
ELL=data.frame(hp,ew=ewL,method="likelihood")
return(ELL)
}

###################################################################################################
#' Plot for Exact method of expected length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' The input can also be a range of values between 0 and 1.
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Confidence interval for \code{p} based on inverting equal-tailed
#' binomial tests with null hypothesis \eqn{H0: p = p0} using expected length of the \eqn{n + 1} intervals.
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;e=0.5;a=1;b=1
#' PlotexplEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=1;a=1;b=1 #Clopper-Pearson
#' PlotexplEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=c(0.1,0.5,0.95,1);a=1;b=1 #Range including Mid-p and Clopper-Pearson
#' PlotexplEX(n,alp,e,a,b)
#' @export
##### 1.EXACT EMTHOD Expected Length for a given n and alpha level
PlotexplEX<-function(n,alp,e,a,b) #n:No of trials,alp:sign level,e: Exact method indicator (1:Clop-Pear,0.5:MID-p),a&b beta parameters for hypo "p'
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  ELEX2=gexplEX(n,alp,e,a,b)
  ELEX2$e=as.factor(ELEX2$e)

  ggplot2::ggplot(ELEX2, ggplot2::aes(x=hp, y=ewEX, color=e))+
    ggplot2::labs(title = "Expected length of Exact method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line()

}

gexplEX<-function(n,alp,e,a,b)
{
nvar=length(e)

res <- data.frame()

for(i in 1:nvar)
{
  lu=ncf302(n,alp,e[i],a,b)
  res <- rbind(res,lu)
}
return(res)
}
ncf302<-function(n,alp,e,a,b) #n:No of trials,alp:sign level,e: Exact method indicator (1:Clop-Pear,0.5:MID-p),a&b beta parameters for hypo "p'
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LEX=0
UEX=0
s=5000
LEEX=0 								#LENGTH OF INTERVAL
ewiEX=matrix(0,k,s)						#Expected length quantity in sum
ewEX=0									#Expected Length

#EXACT METHOD
LEX[1]=0
UEX[1]= 1-((alp/(2*e))^(1/n))
LEX[k]=(alp/(2*e))^(1/n)
UEX[k]=1
LEEX[1]=1-((alp/(2*e))^(1/n))
LEEX[k]=1-((alp/(2*e))^(1/n))

for(i in 1:(k-2))
{
LEX[i+1]=exlim302l(x[i+1],n,alp,e)
UEX[i+1]=exlim302u(x[i+1],n,alp,e)
LEEX[i+1]=UEX[i+1]-LEX[i+1]
}
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiEX[i,j]=LEEX[i]*dbinom(i-1, n,hp[j])
}
ewEX[j]=sum(ewiEX[,j])						#Expected Length
}
sumLEEX=sum(LEEX)
ELEX=data.frame(hp,ewEX,e)

return(ELEX)
}
#####TO FIND LOWER LIMITS
exlim302l=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
  LEX= uniroot(f1,c(0,1))$root
  return(LEX)
}
#####TO FIND UPPER LIMITS
exlim302u=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f2  = function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
  UEX = uniroot(f2,c(0,1))$root
  return(UEX)
}

###############################################################################################################
#' Plot the Bayesian method of expected length calculation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  Plots of Bayesian Highest Probability Density (HPD) and two tailed
#' intervals using expected length of the \eqn{n + 1}
#' intervals for the Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1;a1=1;a2=1
#' PlotexplBA(n,alp,a,b,a1,a2)
#' @export
##### 8.BAYESIAN Expected Length for a given n and alpha level
PlotexplBA<-function(n,alp,a,b,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1)>1 || a1<0 ) stop("'a1' has to be greater than or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2)>1 || a2<0 ) stop("'a2' has to be greater than or equal to 0")


####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LBAQ=0
UBAQ=0
LBAH=0
UBAH=0
s=5000
LEBAQ=0 								#LENGTH OF INTERVAL
LEBAH=0
ewiBAQ=matrix(0,k,s)						#Expected length quantity in sum
ewBAQ=0
ewiBAH=matrix(0,k,s)						#Expected length quantity in sum
ewBAH=0									#Expected Length

#library(TeachingDemos)				#To get HPDs
for(i in 1:k)
{
#Quantile Based Intervals
LBAQ[i]=qbeta(alp/2,x[i]+a1,n-x[i]+a2)
UBAQ[i]=qbeta(1-(alp/2),x[i]+a1,n-x[i]+a2)

LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[1]
UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[2]

LEBAQ[i]=UBAQ[i]-LBAQ[i]
LEBAH[i]=UBAH[i]-LBAH[i]
}
sumLEBAQ=sum(LEBAQ)
sumLEBAH=sum(LEBAH)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiBAQ[i,j]=LEBAQ[i]*dbinom(i-1, n,hp[j])
ewiBAH[i,j]=LEBAH[i]*dbinom(i-1, n,hp[j])

}
ewBAQ[j]=sum(ewiBAQ[,j])
ewBAH[j]=sum(ewiBAH[,j])						#Expected Length
}
ELBAQ=data.frame(hp,ew=ewBAQ,method="Quantile")
ELBAH=data.frame(hp,ew=ewBAH,method="HPD")

df.ba=rbind(ELBAQ,ELBAH)

ggplot2::ggplot(df.ba, ggplot2::aes(x=hp, y=ew))+
  ggplot2::labs(title = "Expected length of Bayesian Quantile & HPD based methods") +
  ggplot2::labs(y = "Expected Length") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_line(ggplot2::aes(color=method)) +
  ggplot2::geom_vline(ggplot2::aes(xintercept=0.5,color="Intercept"),linetype = 2)+
  ggplot2::scale_colour_manual(name='Heading',
                               values=c('Quantile' ='black',
                                        'HPD' = 'red',
                                        'Intercept'= 'brown',
                                        'Confidence Level'='blue'),
                               guide='legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(1,2,1))))

}

#############################################################################################################
#' Plots the Expected length using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for the expected length of \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Expected length  of base methods
#' @examples
#' n= 10; alp=0.05; a=1;b=1;
#' PlotexplAll(n,alp,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplAll<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.eall=  explAll(n,alp,a,b)

  ggplot2::ggplot(df.eall, ggplot2::aes(x=hp, y=ew))+
    ggplot2::labs(title = "Expected length of 6 base methods") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0.5),linetype = 1)

}

#############################################################################################################
#' Plots the expected length for Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Wald-type intervals using sum of length of the \eqn{n + 1}
#'  intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotexplWD(n,alp,a,b)
#' @export
##### 1.WALD sum of length for a given n and alpha level
PlotexplWD<-function(n,alp,a,b) #n:No of trials,alp:sign level,a&b beta parameters for hypo "p'
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  df.wd=  gexplWD(n,alp,a,b)
  ddf.wd = lengthWD(n,alp,a,b)
  df.wd$gMean=ddf.wd$explMean
  df.wd$gMax=ddf.wd$explMax
  df.wd$gUL=ddf.wd$explMean+ddf.wd$explSD
  df.wd$gLL=ddf.wd$explMean-ddf.wd$explSD


  ggplot2::ggplot(data=df.wd, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Wald method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE)  +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))


}

#############################################################################################################
#' Plots the expected length for Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of score test approach using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotexplSC(n,alp,a,b)
#' @export
##### 2.SCORE - sum of length for a given n and alpha level
PlotexplSC<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  df.sc=  gexplSC(n,alp,a,b)
  ddf.sc = lengthSC(n,alp,a,b)
  df.sc$gMean=ddf.sc$explMean
  df.sc$gMax=ddf.sc$explMax
  df.sc$gUL=ddf.sc$explMean+ddf.sc$explSD
  df.sc$gLL=ddf.sc$explMean-ddf.sc$explSD

  ggplot2::ggplot(data=df.sc, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Score method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

  }

#############################################################################################################
#' Plots ArcSine method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Wald-type interval for the arcsine transformation of the parameter
#' \code{p} using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotexplAS(n,alp,a,b)
#' @export
##### 3. ARC SINE - sum of length for a given n and alpha level
PlotexplAS<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  df.as=  gexplAS(n,alp,a,b)

  ddf.as = lengthAS(n,alp,a,b)
  df.as$gMean=ddf.as$explMean
  df.as$gMax=ddf.as$explMax
  df.as$gUL=ddf.as$explMean+ddf.as$explSD
  df.as$gLL=ddf.as$explMean-ddf.as$explSD

  ggplot2::ggplot(data=df.as, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of ArcSine method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

}

#############################################################################################################
#' Plots Logit Wald method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Wald-type interval based on the logit
#' transformation of \code{p} using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotexplLT(n,alp,a,b)
#' @export
##### 4.LOGIT-WALD - sum of length for a given n and alpha level
PlotexplLT<-function(n,alp,a,b) #n:No of trials,alp:sign level,a&b beta parameters for hypo "p'
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  df.lt=  gexplLT(n,alp,a,b)
  ddf.lt = lengthLT(n,alp,a,b)
  df.lt$gMean=ddf.lt$explMean
  df.lt$gMax=ddf.lt$explMax
  df.lt$gUL=ddf.lt$explMean+ddf.lt$explSD
  df.lt$gLL=ddf.lt$explMean-ddf.lt$explSD

  ggplot2::ggplot(data=df.lt, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Logit Wald method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

}

#############################################################################################################
#' Plots Wald-T method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of approximate method based on a t_approximation of the
#' standardized point estimator using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotexplTW(n,alp,a,b)
#' @export
##### 5.t-WALD - sum of length for a given n and alpha level
PlotexplTW<-function(n,alp,a,b) #n:No of trials,alp:sign level,a&b beta parameters for hypo "p'

{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  df.tw=  gexplTW(n,alp,a,b)
  ddf.tw = lengthTW(n,alp,a,b)
  df.tw$gMean=ddf.tw$explMean
  df.tw$gMax=ddf.tw$explMax
  df.tw$gUL=ddf.tw$explMean+ddf.tw$explSD
  df.tw$gLL=ddf.tw$explMean-ddf.tw$explSD

  ggplot2::ggplot(data=df.tw, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Wald-T method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))


}

#############################################################################################################
#' Plots likelihood Ratio method of expected length
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Likelihood ratio limits using sum of length of the \eqn{n + 1} intervals
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' PlotexplLR(n,alp,a,b)
#' @export
#####6.LIKELIHOOD RATIO - sum of length for a given n and alpha level
PlotexplLR<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  df.lr=  gexplLR(n,alp,a,b)
  ddf.lr = lengthLR(n,alp,a,b)
  df.lr$gMean=ddf.lr$explMean
  df.lr$gMax=ddf.lr$explMax
  df.lr$gUL=ddf.lr$explMean+ddf.lr$explSD
  df.lr$gLL=ddf.lr$explMean-ddf.lr$explSD

  ggplot2::ggplot(data=df.lr, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of Likelihood Ratio method") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(mapping=ggplot2::aes(colour=method), show_guide = TRUE) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMean, fill="Mean"),color="orange"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gMax, fill="Max"),color="blue"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gLL, fill="Lower Limit"),color="cyan4"  ) +
    ggplot2::geom_hline(mapping=ggplot2::aes(yintercept=gUL, fill="Upper Limit"),color="brown"  ) +
    ggplot2::scale_color_hue("Method") +
    ggplot2::scale_fill_manual(
      "Metric lines", values=c(1,1,1,1),
      guide=ggplot2::guide_legend(override.aes = list(colour=c("orange", "blue", "cyan4","brown"))),
      labels=c("Mean", "Max", "Lower Limit(Mean- 1SD)", "Upper Limit(Mean + 1SD)"))

}
