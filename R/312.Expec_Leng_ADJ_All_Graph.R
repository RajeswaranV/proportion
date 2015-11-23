##### 1.ADJUSTED WALD Expected Length for a given n and alpha level
gexplAWD<-function(n,alp,h,a,b)
{
####INPUT n
x=0:n
k=n+1
y=x+h
m=n+(2*h)
####INITIALIZATIONS
pAW=0
qAW=0
seAW=0
LAW=0
UAW=0
s=2000
LEAW=0 								#LENGTH OF INTERVAL

ewiAW=matrix(0,k,s)						#Expected length quantity in sum
ewAW=0									#Expected Length
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pAW[i]=y[i]/m
qAW[i]=1-pAW[i]
seAW[i]=sqrt(pAW[i]*qAW[i]/m)
LAW[i]=pAW[i]-(cv*seAW[i])
UAW[i]=pAW[i]+(cv*seAW[i])
if(LAW[i]<0) LAW[i]=0
if(UAW[i]>1) UAW[i]=1
LEAW[i]=UAW[i]-LAW[i]
}
sumLEAW=sum(LEAW)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiAW[i,j]=LEAW[i]*dbinom(i-1, n,hp[j])
}
ewAW[j]=sum(ewiAW[,j])						#Expected Length
}
ELAW=data.frame(hp,ew=ewAW,method="Adj-Wald")
return(ELAW)
}
###############################################################################################################
##### 2.ADJUSTED SCORE - Expected Length for a given n and alpha level
gexplASC<-function(n,alp,h,a,b)
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
s=2000
LEAS=0 								#LENGTH OF INTERVAL

ewiAS=matrix(0,k,s)						#Expected length quantity in sum
ewAS=0									#Expected Length
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
LEAS[i]=UAS[i]-LAS[i]
}
sumLEAS=sum(LEAS)

####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiAS[i,j]=LEAS[i]*dbinom(i-1, n,hp[j])
}
ewAS[j]=sum(ewiAS[,j])						#Expected Length
}
ELAS=data.frame(hp,ew=ewAS,method="Adj-Wilson")
return(ELAS)
}
###############################################################################################################
##### 3. ADJUSTED ARC SINE - Expected Length for a given n and alpha level
gexplAAS<-function(n,alp,h,a,b)
{
####INPUT n
x=0:n
k=n+1
y=x+h
m=n+(2*h)
####INITIALIZATIONS
pAA=0
qAA=0
seAA=0
LAA=0
UAA=0
s=2000
LEAA=0 								#LENGTH OF INTERVAL

ewiAA=matrix(0,k,s)						#Expected length quantity in sum
ewAA=0									#Expected Length
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pAA[i]=y[i]/m
qAA[i]=1-pAA[i]
seAA[i]=cv/sqrt(4*m)
LAA[i]=(sin(asin(sqrt(pAA[i]))-seAA[i]))^2
UAA[i]=(sin(asin(sqrt(pAA[i]))+seAA[i]))^2
if(LAA[i]<0) LAA[i]=0
if(UAA[i]>1) UAA[i]=1
LEAA[i]=UAA[i]-LAA[i]
}
sumLEAA=sum(LEAA)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiAA[i,j]=LEAA[i]*dbinom(i-1, n,hp[j])
}
ewAA[j]=sum(ewiAA[,j])						#Expected Length
}
ELAA=data.frame(hp,ew=ewAA,method="Adj-ArcSine")
return(ELAA)
}

###############################################################################################################
##### 4.ADJUSTED LOGIT-WALD - Expected Length for a given n and alpha level
gexplALT<-function(n,alp,h,a,b)
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
s=2000
LEALT=0 								#LENGTH OF INTERVAL

ewiALT=matrix(0,k,s)						#Expected length quantity in sum
ewALT=0									#Expected Length
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD

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
LEALT[i]=UALT[i]-LALT[i]
}
sumLEALT=sum(LEALT)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiALT[i,j]=LEALT[i]*dbinom(i-1, n,hp[j])
}
ewALT[j]=sum(ewiALT[,j])						#Expected Length
}
ELALT=data.frame(hp,ew=ewALT,method="Adj-Logit-Wald")
return(ELALT)
}

###############################################################################################################
##### 5.ADJUSTED t-WALD - Expected Length for a given n and alpha level
gexplATW<-function(n,alp,h,a,b)
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
s=2000
LEATW=0 								#LENGTH OF INTERVAL

ewiATW=matrix(0,k,s)						#Expected length quantity in sum
ewATW=0									#Expected Length
#MODIFIED_t-ADJ_WALD METHOD
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
LEATW[i]=UATW[i]-LATW[i]
}
sumLEATW=sum(LEATW)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiATW[i,j]=LEATW[i]*dbinom(i-1, n,hp[j])
}
ewATW[j]=sum(ewiATW[,j])						#Expected Length
}
ELATW=data.frame(hp,ew=ewATW,method="Adj-Wald-T")
return(ELATW)
}
###############################################################################################################
#####6.ADJUSTED LIKELIHOOD RATIO - Expected Length for a given n and alpha level
gexplALR<-function(n,alp,h,a,b)
{
####INPUT n
y=0:n
y1=y+h
k=n+1
n1=n+(2*h)
####INITIALIZATIONS
mle=0
cutoff=0
LAL=0
UAL=0
s=2000
LEAL=0 								#LENGTH OF INTERVAL

ewiAL=matrix(0,k,s)						#Expected length quantity in sum
ewAL=0									#Expected Length
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LIKELIHOOD-RATIO METHOD
for(i in 1:k)
{
likelhd = function(p) dbinom(y1[i],n1,p)
loglik = function(p) dbinom(y1[i],n1,p,log=TRUE)
mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
cutoff[i]=loglik(mle[i])-(cv^2/2)
loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
LAL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
UAL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
LEAL[i]=UAL[i]-LAL[i]
}
sumLEAL=sum(LEAL)
####Expected Length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiAL[i,j]=LEAL[i]*dbinom(i-1, n,hp[j])
}
ewAL[j]=sum(ewiAL[,j])						#Expected Length
}
ELAL=data.frame(hp,ew=ewAL,method="Adj-Likelyhood")
return(ELAL)
}

#############################################################################################################
#' Plots the Expected length using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)  for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplAAll(n,alp,h,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplAAll<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.new=  explAAll(n,alp,h,a,b)

  ggplot2::ggplot(df.new, ggplot2::aes(x=hp, y=ew))+
    ggplot2::labs(title = "Expected length of 6 adjusted methods") +
    ggplot2::labs(y = "Expected length") +
    ggplot2::labs(x = "p") +
    ggplot2::geom_line(ggplot2::aes(color=method)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0.5),linetype = 1)

}

#############################################################################################################
#' Plots the Expected length using adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplAWD(n,alp,h,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplAWD<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.awd=  gexplAWD(n,alp,h,a,b)
  ddf.awd = lengthAWD(n,alp,h,a,b)
  df.awd$gMean=ddf.awd$explMean
  df.awd$gMax=ddf.awd$explMax
  df.awd$gUL=ddf.awd$explMean+ddf.awd$explSD
  df.awd$gLL=ddf.awd$explMean-ddf.awd$explSD

  ggplot2::ggplot(data=df.awd, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Wald method") +
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
#' Plots the Expected length using adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Score method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplASC(n,alp,h,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplASC<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.asc=  gexplASC(n,alp,h,a,b)
  ddf.asc = lengthASC(n,alp,h,a,b)
  df.asc$gMean=ddf.asc$explMean
  df.asc$gMax=ddf.asc$explMax
  df.asc$gUL=ddf.asc$explMean+ddf.asc$explSD
  df.asc$gLL=ddf.asc$explMean-ddf.asc$explSD

  ggplot2::ggplot(data=df.asc, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Score method") +
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
#' Plots the Expected length using adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplATW(n,alp,h,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplATW<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.atw=  gexplATW(n,alp,h,a,b)
  ddf.atw = lengthATW(n,alp,h,a,b)
  df.atw$gMean=ddf.atw$explMean
  df.atw$gMax=ddf.atw$explMax
  df.atw$gUL=ddf.atw$explMean+ddf.atw$explSD
  df.atw$gLL=ddf.atw$explMean-ddf.atw$explSD

  ggplot2::ggplot(data=df.atw, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Wald-T method") +
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
#' Plots the Expected length using adjusted Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Wald method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplALT(n,alp,h,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplALT<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.alt=  gexplALT(n,alp,h,a,b)
  ddf.alt = lengthALT(n,alp,h,a,b)
  df.alt$gMean=ddf.alt$explMean
  df.alt$gMax=ddf.alt$explMax
  df.alt$gUL=ddf.alt$explMean+ddf.alt$explSD
  df.alt$gLL=ddf.alt$explMean-ddf.alt$explSD

  ggplot2::ggplot(data=df.alt, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Logit Wald method") +
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
#' Plots the Expected length using adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted ArcSine method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplAAS(n,alp,h,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplAAS<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.aas=  gexplAAS(n,alp,h,a,b)
  ddf.aas = lengthAAS(n,alp,h,a,b)
  df.aas$gMean=ddf.aas$explMean
  df.aas$gMax=ddf.aas$explMax
  df.aas$gUL=ddf.aas$explMean+ddf.aas$explSD
  df.aas$gLL=ddf.aas$explMean-ddf.aas$explSD

  ggplot2::ggplot(data=df.aas, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Logit Wald method") +
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
#' Plots the Expected length using adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The  plots of the Expected length of adjusted Likelihood Ratio method
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' PlotexplALR(n,alp,h,a,b)
#' @export
##### 9.All methods - Expected length
PlotexplALR<-function(n,alp,h,a,b)
{

  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.alr=  gexplALR(n,alp,h,a,b)
  ddf.alr = lengthALR(n,alp,h,a,b)
  df.alr$gMean=ddf.alr$explMean
  df.alr$gMax=ddf.alr$explMax
  df.alr$gUL=ddf.alr$explMean+ddf.alr$explSD
  df.alr$gLL=ddf.alr$explMean-ddf.alr$explSD

  ggplot2::ggplot(data=df.alr, mapping=ggplot2::aes(x=hp, y=ew)) +
    ggplot2::labs(title = "Expected length of adjusted Likelihood Ratio method") +
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
