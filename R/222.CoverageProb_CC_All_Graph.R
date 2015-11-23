##### 1.CC-WALD-Coverage Probability
gcovpCWD<-function(n,alp,c,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1

####INITIALIZATIONS
pCW=0
qCW=0
seCW=0
LCW=0
UCW=0
s=5000								#Simulation run to generate hypothetical p
cpCW=matrix(0,k,s)
ctCW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCW=0								#Coverage probabilty
ctr=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pCW[i]=x[i]/n
qCW[i]=1-pCW[i]
seCW[i]=sqrt(pCW[i]*qCW[i]/n)
LCW[i]=pCW[i]-((cv*seCW[i])+c)
UCW[i]=pCW[i]+((cv*seCW[i])+c)
if(LCW[i]<0) LCW[i]=0
if(UCW[i]>1) UCW[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCW[i] && hp[j] < UCW[i])
{
cpCW[i,j]=dbinom(i-1, n,hp[j])
ctCW[i,j]=1
}
}
cppCW[j]=sum(cpCW[,j])
if(t1<cppCW[j]&&cppCW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPCW=data.frame(hp,cp=cppCW,method="Continuity corrected Wald")
# windows()
# plot(hp,cppCW,type="l",xlab="p",ylab="Coverage Probability",main="Wald_CC")
# abline(h=1-(alp), lty=2)
return(CPCW)
}

##### 2.CC-SCORE - Coverage Probability
gcovpCSC<-function(n,alp,c,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCS=0
qCS=0
seCS_L=0
seCS_U=0
LCS=0
UCS=0
s=5000								#Simulation run to generate hypothetical p
cpCS=matrix(0,k,s)
ctCS=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCS=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n)
cv2=cv/(2*n)

#SCORE (WILSON) METHOD
for(i in 1:k)
{
pCS[i]=x[i]/n
qCS[i]=1-pCS[i]
seCS_L[i]=sqrt((cv^2)-(4*n*(c+c^2))+(4*n*pCS[i]*(1-pCS[i]+(2*c))))	#Sq. root term of LL
seCS_U[i]=sqrt((cv^2)+(4*n*(c-c^2))+(4*n*pCS[i]*(1-pCS[i]-(2*c))))	#Sq. root term of LL
LCS[i]=(n/(n+(cv)^2))*((pCS[i]-c+cv1)-(cv2*seCS_L[i]))
UCS[i]=(n/(n+(cv)^2))*((pCS[i]+c+cv1)+(cv2*seCS_U[i]))
if(LCS[i]<0) LCS[i]=0
if(UCS[i]>1) UCS[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCS[i] && hp[j] < UCS[i])
{
cpCS[i,j]=dbinom(i-1, n,hp[j])
ctCS[i,j]=1
}
}
cppCS[j]=sum(cpCS[,j])						#Coverage Probability
if(t1<cppCS[j]&&cppCS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPCS=data.frame(hp,cp=cppCS,method="Continuity corrected Wilson")
# windows()
# plot(hp,cppCS,type="l",xlab="p",ylab="Coverage Probability",main="Wilson_CC")
# abline(h=1-(alp), lty=2)
return(CPCS)
}


##### 3.CC-ARC SINE - Coverage Probability
gcovpCAS<-function(n,alp,c,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCA=0
qCA=0
seCA=0
LCA=0
UCA=0
s=5000								#Simulation run to generate hypothetical p
cpCA=matrix(0,k,s)
ctCA=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCA=0									#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pCA[i]=x[i]/n
qCA[i]=1-pCA[i]
seCA[i]=cv/sqrt(4*n)
LCA[i]=(sin(asin(sqrt(pCA[i]))-seCA[i]-c))^2
UCA[i]=(sin(asin(sqrt(pCA[i]))+seCA[i]+c))^2
if(LCA[i]<0) LCA[i]=0
if(UCA[i]>1) UCA[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCA[i] && hp[j] < UCA[i])
{
cpCA[i,j]=dbinom(i-1, n,hp[j])
ctCA[i,j]=1
}
}
cppCA[j]=sum(cpCA[,j])						#Coverage Probability
if(t1<cppCA[j]&&cppCA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPCA=data.frame(hp,cp=cppCA,method="Continuity corrected ArcSine")
# windows()
# plot(hp,cppCA,type="l",xlab="p",ylab="Coverage Probability",main="Arc Sine_CC")
# abline(h=1-(alp), lty=2)
return(CPCA)
}
##### 4.LOGIT-WALD - Coverage Probability
gcovpCLT<-function(n,alp,c,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCLT=0
qCLT=0
seCLT=0
lgit=0
LCLT=0
UCLT=0
s=5000								#Simulation run to generate hypothetical p
cpCLT=matrix(0,k,s)
ctCLT=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCLT=0								#Coverage probabilty
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
pCLT[1]=0
qCLT[1]=1
LCLT[1] = 0
UCLT[1] = 1-((alp/2)^(1/n))

pCLT[k]=1
qCLT[k]=0
LCLT[k]= (alp/2)^(1/n)
UCLT[k]=1

lgiti=function(t) exp(t)/(1+exp(t))	#LOGIT INVERSE
for(j in 1:(k-2))
{
pCLT[j+1]=x[j+1]/n
qCLT[j+1]=1-pCLT[j+1]
lgit[j+1]=log(pCLT[j+1]/qCLT[j+1])
seCLT[j+1]=sqrt(pCLT[j+1]*qCLT[j+1]*n)
LCLT[j+1]=lgiti(lgit[j+1]-(cv/seCLT[j+1])-c)
UCLT[j+1]=lgiti(lgit[j+1]+(cv/seCLT[j+1])+c)
}
for(i in 1:k)
{
if(LCLT[i]<0) LCLT[i]=0
if(UCLT[i]>1) UCLT[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCLT[i] && hp[j] < UCLT[i])
{
cpCLT[i,j]=dbinom(i-1, n,hp[j])
ctCLT[i,j]=1
}
}
cppCLT[j]=sum(cpCLT[,j])				#Coverage Probability
if(t1<cppCLT[j]&&cppCLT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

}
CPLT=data.frame(hp,cp=cppCLT,method="Continuity corrected Logit Wald")
return(CPLT)
}

##### 5. WALD_t - Coverage Probability
gcovpCTW<-function(n,alp,c,a,b,t1,t2)
{
####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCTW=0
qCTW=0
seCTW=0
LCTW=0
UCTW=0
DOF=0
cv=0
s=5000								#Simulation run to generate hypothetical p
cpCTW=matrix(0,k,s)
ctCTW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCTW=0
ctr=0

#MODIFIED_t-WALD METHOD
for(i in 1:k)
{
if(x[i]==0||x[i]==n)
{
pCTW[i]=(x[i]+2)/(n+4)
qCTW[i]=1-pCTW[i]
}else
{
pCTW[i]=x[i]/n
qCTW[i]=1-pCTW[i]
}
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOF[i]=2*((f1(pCTW[i],n))^2)/f2(pCTW[i],n)
cv[i]=qt(1-(alp/2), df=DOF[i])
seCTW[i]=cv[i]*sqrt(f1(pCTW[i],n))
LCTW[i]=pCTW[i]-(seCTW[i]+c)
UCTW[i]=pCTW[i]+(seCTW[i]+c)
if(LCTW[i]<0) LCTW[i]=0
if(UCTW[i]>1) UCTW[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCTW[i] && hp[j] < UCTW[i])
{
cpCTW[i,j]=dbinom(i-1, n,hp[j])
ctCTW[i,j]=1
}
}
cppCTW[j]=sum(cpCTW[,j])						#Coverage Probability
if(t1<cppCTW[j]&&cppCTW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
CPCTW=data.frame(hp,cp=cppCTW,method="Continuity corrected Wald-T")
# windows()
# plot(hp,cppCTW,type="l",xlab="p",ylab="Coverage Probability",main="Logit-Wald_CC")
# abline(h=1-(alp), lty=2)
return(CPCTW)
}

#############################################################################################################
#' Graphs  of Coverage Probability for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The graphs of the  Coverage Probability of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCAll(n,alp,c,a,b,t1,t2)
#' @export
##### 1.All methods - Coverage Probability 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
PlotcovpCAll<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")


  #### Calling functions and creating df
  df1    = gcovpCWD(n,alp,c,a,b,t1,t2)
  df2    = gcovpCAS(n,alp,c,a,b,t1,t2)
  df3    = gcovpCSC(n,alp,c,a,b,t1,t2)
  df4    = gcovpCLT(n,alp,c,a,b,t1,t2)
  df5    = gcovpCTW(n,alp,c,a,b,t1,t2)

  g.df= rbind(df1,df2,df3,df4,df5)

  ggplot2::ggplot(g.df, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected methods") +
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
#' Graphs  of Coverage Probability for continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Wald method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCWD(n,alp,c,a,b,t1,t2)
#' @export
PlotcovpCWD<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  df1    = gcovpCWD(n,alp,c,a,b,t1,t2)
  WaldcovpA.df    = covpCWD(n,alp,c,a,b,t1,t2)
  df1$mcp=WaldcovpA.df$mcpCW
  df1$micp=WaldcovpA.df$micpCW

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Wald method") +
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
#' Graphs  of Coverage Probability for continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected ArcSine method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCAS(n,alp,c,a,b,t1,t2)
#' @export
PlotcovpCAS<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  df1    = gcovpCAS(n,alp,c,a,b,t1,t2)
  ArcSinecovpA.df = covpCAS(n,alp,c,a,b,t1,t2)

  df1$mcp=ArcSinecovpA.df$mcpCA
  df1$micp=ArcSinecovpA.df$micpCA

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected ArcSine method") +
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
#' Graphs  of Coverage Probability for continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Score method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCSC(n,alp,c,a,b,t1,t2)
#' @export
PlotcovpCSC<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  df1    = gcovpCSC(n,alp,c,a,b,t1,t2)
  ScorecovpA.df   = covpCSC(n,alp,c,a,b,t1,t2)

  df1$mcp=ScorecovpA.df$mcpCS
  df1$micp=ScorecovpA.df$micpCS

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Score method") +
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
#' Graphs  of Coverage Probability for continuity corrected Logistic Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Logistic Wald method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCLT(n,alp,c,a,b,t1,t2)
#' @export
PlotcovpCLT<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  df1    = gcovpCLT(n,alp,c,a,b,t1,t2)
  WaldLcovpA.df   = covpCLT(n,alp,c,a,b,t1,t2)

  df1$mcp=WaldLcovpA.df$mcpCLT
  df1$micp=WaldLcovpA.df$micpCLT

  ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Logistic Wald method") +
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
#' Graphs  of Coverage Probability for continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The plot of the  Coverage Probability of continuity corrected Wald-T method
#' for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' PlotcovpCTW(n,alp,c,a,b,t1,t2)
#' @export
PlotcovpCTW<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  df1    = gcovpCTW(n,alp,c,a,b,t1,t2)

  AdWaldcovpA.df  = covpCTW(n,alp,c,a,b,t1,t2)

  df1$mcp=AdWaldcovpA.df$mcpCTW
  df1$micp=AdWaldcovpA.df$micpCTW

ggplot2::ggplot(df1, ggplot2::aes(x=hp, y=cp))+
    ggplot2::labs(title = "Coverage Probability of Continuity corrected Wald-T method") +
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
