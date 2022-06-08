## ---- echo=FALSE, results='asis'----------------------------------------------
# 00.Naming convention
Naming.convention = matrix( c("ci, ciA, ciC","Confidence Interval, adjusted CI and continuity corrected CI",
"covp, covpA, covpC","Coverage Probability, adjusted CP & continuity corrected CP",
"expl, explA, explC ","Expected Length, adjusted Expected Length & continuity corrected EL",
"length, lengthA, lengthC","Sum of Length, adjusted Sum of Length & continuity corrected sumLen",
"pCOpBI, pCOpBIA, pCOpBIC","p-Confidence & p-Bias, adjusted p-Conf & p-Bias ",
"","and continuity corrected p-Confidence & p-Bias",
"err, errA, errC","Error, adjusted error and continuity corrected error",
"AS","ArcSine",
"LR","Likelihood Ratio",
"LT","Logit Wald",
"SC","Score",
"TW","Wald-T",
"WD","Wald",
"All","6 base methods - Wald, Wald-T, Logit Wald, ArcSine, LR, Score",
"AAll","6 adj methods - Wald, Wald-T, Logit Wald, ArcSine, LR, Score",
"CAll","5 cont. corr. methods - Wald, Wald-T, Logit Wald, ArcSine,  Score",
"BA","Bayesian",
"EX" ,"Exact - setting e=0.5 gives mid-p and e=1 gives Clopper-Pearson"),byrow=TRUE, nrow=18,  ncol=2)

colnames(Naming.convention)=c("Abbrivation", "Expansion")

knitr::kable(Naming.convention,caption = "Naming convention used in functions")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 000.Function identification
Identify.Function = matrix( c("Plot","ci","A","AS","x","ci + A + AS + x =","ciAASx",
"","covp","C","SC","","Plot + ci + A + AS + x =","PlotciAASx",
"","expl","","BA","","Plot + covp + C + SC =","PlotcovpCSC",
"","length","","EX","","expl + A + TW =","explATW",
"","pCOpBI","","TW","","expl + A + TW + x =","explATWx",
"","err","","LT","","length + WD =","lenghtWD",
"","","","WD","","length + A + WD =","lengthAWD",
"","","","LR","","length + C + WD =","lengthCWD"),byrow=TRUE, nrow=8,  ncol=7)

colnames(Identify.Function)=c("Plot", "Concept","Modifications", "Name","Single x","Sample combination","Sample function")

knitr::kable(Identify.Function,caption = "Guide to identify core functions - Plot, Modifications and x are optional")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 1.ci
Numeric.ci = matrix( c("ciAS","ciASx","ciAAS","ciAASx","ciCAS","ciCASx",
"ciLR","ciLRx","ciALR","ciALRx","","",
"ciLT","ciLTx","ciALT","ciALTx","ciCLT","ciCLTx",
"ciSC","ciSCx","ciASC","ciASCx","ciCSC","ciCSCx",
"ciTW","ciTWx","ciATW","ciATWx","ciCTW","ciCTWx",
"ciWD","ciWDx","ciAWD","ciAWDx","ciCWD","ciCWDx",
"ciAll","ciAllx","ciAAll","ciAAllx","ciCAll","ciCAllx",
"ciBA","ciBAx","","","","",
"ciEX","ciEXx","","","",""),byrow=TRUE, nrow=9,  ncol=6)


colnames(Numeric.ci)=c("Basic","Basic-x","Adj","Adj-x","CC","CC-x")
rownames(Numeric.ci)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Numeric.ci,caption = "Confidence Interval")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 6. Plot ci
Plot.ci = matrix( c("PlotciAS","","PlotciAAS","","PlotciCAS","",
"PlotciLR","","PlotciALR","","","",
"PlotciLT","","PlotciALT","","PlotciCLT","",
"PlotciSC","","PlotciASC","","PlotciCSC","",
"PlotciTW","","PlotciATW","","PlotciCTW","",
"PlotciWD","","PlotciAWD","","PlotciCWD","",
"PlotciAllg","PlotciAllxg","PlotciAAllg","PlotciAAllxg","PlotciCAllg","PlotciCAllxg",
"PlotciAll","PlotciAllx","PlotciAAll","PlotciAAllx","PlotciCAll","PlotciCAllx",
"PlotciBA","","","","","",
"PlotciEX","PlotciEXx","","","",""),byrow=TRUE, nrow=10,  ncol=6)

colnames(Plot.ci)=c("Basic","Basic-x","Adj","Adj-x","CC","CC-x")
rownames(Plot.ci)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","Allg","All","Bayes","Exact")

knitr::kable(Plot.ci,caption = "Plotting functions of CI")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 2. covp
Numeric.covp = matrix( c("covpAS","covpAAS","covpCAS",
"covpLR","covpALR","",
"covpLT","covpALT","covpCLT",
"covpSC","covpASC","covpCSC",
"covpTW","covpATW","covpCTW",
"covpWD","covpAWD","covpCWD",
"covpAll","covpAAll","covpCAll",
"covpBA","","",
"covpEX","",""),byrow=TRUE, nrow=9,  ncol=3)


colnames(Numeric.covp)=c("Basic","Adjusted","Continuity corrected")
rownames(Numeric.covp)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Numeric.covp,caption = "Coverage Probability")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 7. Plot covp
Plot.covp = matrix( c("PlotcovpAS","PlotcovpAAS","PlotcovpCAS",
"PlotcovpLR","PlotcovpALR","",
"PlotcovpLT","PlotcovpALT","PlotcovpCLT",
"PlotcovpSC","PlotcovpASC","PlotcovpCSC",
"PlotcovpTW","PlotcovpATW","PlotcovpCTW",
"PlotcovpWD","PlotcovpAWD","PlotcovpCWD",
"PlotcovpAll","PlotcovpAAll","PlotcovpCAll",
"PlotcovpBA","","",
"PlotcovpEX","",""),byrow=TRUE, nrow=9,  ncol=3)

colnames(Plot.covp)=c("Basic","Adjusted","Continuity corrected")
rownames(Plot.covp)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Plot.covp,caption = "Plotting functions of Coverage Probability")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 3. Length
Numeric.length = matrix( c("lengthAS","lengthAAS","lengthCAS",
"lengthLR","lengthALR","",
"lengthLT","lengthALT","lengthCLT",
"lengthSC","lengthASC","lengthCSC",
"lengthTW","lengthATW","lengthCTW",
"lengthWD","lengthAWD","lengthCWD",
"lengthAll","lengthAAll","lengthCAll",
"lengthBA","","",
"lengthEX","",""),byrow=TRUE, nrow=9,  ncol=3)

colnames(Numeric.length)=c("SumLen","Adj-SumLen","CC-SumLen")
rownames(Numeric.length)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Numeric.length,caption = "Sum of length")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 8. Plot length
Plot.covp = matrix( c("PlotlengthAS","PlotexplAS","PlotlengthAAS","PlotexplAAS","PlotlengthCAS","PlotexplCAS",
"PlotlengthLR","PlotexplLR","PlotlengthALR","PlotexplALR","","",
"PlotlengthLT","PlotexplLT","PlotlengthALT","PlotexplALT","PlotlengthCLT","PlotexplCLT",
"PlotlengthSC","PlotexplSC","PlotlengthASC","PlotexplASC","PlotlengthCSC","PlotexplCSC",
"PlotlengthTW","PlotexplTW","PlotlengthATW","PlotexplATW","PlotlengthCTW","PlotexplCTW",
"PlotlengthWD","PlotexplWD","PlotlengthAWD","PlotexplAWD","PlotlengthCWD","PlotexplCWD",
"PlotlengthAll","PlotexplAll","PlotlengthAAll","PlotexplAAll","PlotlengthCAll","PlotexplCAll",
"PlotlengthBA","PlotexplBA","","","","",
"PlotlengthEX","PlotexplEX","","","",""),byrow=TRUE, nrow=9,  ncol=6)

#colnames(Plot.covp)=c("Basic","Adjusted","Continuity corrected")
colnames(Plot.covp)=c("SumLen","EL","Adj-SumLen","Adj-EL","CC-SumLen","CC-EL")
rownames(Plot.covp)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Plot.covp,caption = "Plotting functions of sum length and expected length (EL)")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 4. p-Confidence & p-Bias	
Numeric.pcpb = matrix( c("pCOpBIAS","pCOpBIAAS","pCOpBICAS",
"pCOpBILR","pCOpBIALR","",
"pCOpBILT","pCOpBIALT","pCOpBICLT",
"pCOpBISC","pCOpBIASC","pCOpBICSC",
"pCOpBITW","pCOpBIATW","pCOpBICTW",
"pCOpBIWD","pCOpBIAWD","pCOpBICWD",
"pCOpBIAll","pCOpBIAAll","pCOpBICAll",
"pCOpBIBA","","",
"pCOpBIEX","",""),byrow=TRUE, nrow=9,  ncol=3)

colnames(Numeric.pcpb)=c("Basic","Adjusted","Continuity corrected")
rownames(Numeric.pcpb)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Numeric.pcpb,caption = "p-Confidence & p-Bias")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 9. Plot p-Confidence & p-Bias
Plot.pcpb = matrix( c("PlotpCOpBIAS","PlotpCOpBIAAS","PlotpCOpBICAS",
"PlotpCOpBILR","PlotpCOpBIALR","",
"PlotpCOpBILT","PlotpCOpBIALT","PlotpCOpBICLT",
"PlotpCOpBISC","PlotpCOpBIASC","PlotpCOpBICSC",
"PlotpCOpBITW","PlotpCOpBIATW","PlotpCOpBICTW",
"PlotpCOpBIWD","PlotpCOpBIAWD","PlotpCOpBICWD",
"PlotpCOpBIAll","PlotpCOpBIAAll","PlotpCOpBICAll",
"PlotpCOpBIBA","","",
"PlotpCOpBIEX","",""),byrow=TRUE, nrow=9,  ncol=3)


colnames(Plot.pcpb)=c("Basic","Adjusted","Continuity corrected")
rownames(Plot.pcpb)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Plot.pcpb,caption = "Plotting functions for p-Confidence & p-Bias")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 5. Error and long term power
Numeric.error = matrix( c("errAS","errAAS","errCAS",
"errLR","errALR","",
"errLT","errALT","errCLT",
"errSC","errASC","errCSC",
"errTW","errATW","errCTW",
"errWD","errAWD","errCWD",
"errAll","errAAll","errCAll",
"errBA","","",
"errEX","",""),byrow=TRUE, nrow=9,  ncol=3)

colnames(Numeric.error)=c("Basic","Adjusted","Continuity corrected")
rownames(Numeric.error)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Numeric.error,caption = "Error and long term power")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 10. Plot error and long term power
Plot.error = matrix( c("PloterrAS","PloterrAAS","PloterrCAS",
"PloterrLR","PloterrALR","",
"PloterrLT","PloterrALT","PloterrCLT",
"PloterrSC","PloterrASC","PloterrCSC",
"PloterrTW","PloterrATW","PloterrCTW",
"PloterrWD","PloterrAWD","PloterrCWD",
"PloterrAll","PloterrAAll","PloterrCAll",
"PloterrBA","","",
"PloterrEX","",""),byrow=TRUE, nrow=9,  ncol=3)


colnames(Plot.error)=c("Basic","Adjusted","Continuity corrected")
rownames(Plot.error)=c("ArcSine", "LR","Logit","Score","Wald-T","Wald","All","Bayes","Exact")

knitr::kable(Plot.error,caption = "Plotting functions for error and long term power")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 11. Others
Other.functions = matrix( c("hypotestBAF1","covpGEN","lengthGEN","pCOpBIGEN","empericalBA","errGEN",
"hypotestBAF1x","PlotcovpGEN","PlotlengthGEN","PlotpCOpBIGEN","empericalBAx","",
"hypotestBAF2x","covpSIM","lengthSIM","","probPOSx","",
"hypotestBAF2","PlotcovpSIM","PlotlengthSIM","","probPOS","",
"hypotestBAF3x","","PlotexplGEN","","probPREx","",
"hypotestBAF3","","PlotexplSIM","","probPRE","",
"hypotestBAF4x","","","","","",
"hypotestBAF4","","","","","",
"hypotestBAF5x","","","","","",
"hypotestBAF5","","","","","",
"hypotestBAF6x","","","","","",
"hypotestBAF6","","","","",""),byrow=TRUE, nrow=12,  ncol=6)

colnames(Other.functions)=c("Hypothesis","covp","length","pCOpBI","Others","Error")

knitr::kable(Other.functions,caption = "Additional functions")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 12. Papers
Key.Papers = matrix( 
c("1","20","0","Newcombe","Wald ,Score,(both with","ciWDx,ciSCx","Methods such as Bayesian,Arcsine,",
"2","29","1","","and without CC) Exact","ciCWDx,","Logit Wald methods; Numerical",
"3","148","15","","and LR for CI","ciCSCx,","& graphical comparisons of methods",
"4","263","81","","","ciEXx,ciLRx","Use of general CC and adj. factor",
"5","10","10","Joseph 2005","Wald and Exact CI","ciWDx,ciEXx","Bayes factor",
"6","98","100","","","","",
"7","17","16","Zhou 2008","Wald, Score,","ciWDx, ciSCx","Other methods such as Bayesian,",
"8","14","12",""," Agresti-Coull &","ciAWDx","Arcsine Logit transformed methods",
"","","","","modified logit for CI","","Use of general CC and adj. factor",
"9","167","0","Wei 2012","Score, Agresti-Coull","ciSCx,","Other classical methods; Numerical",
"","","","","Bayesian(Jeffreys prior)","ciAWDx,","& graphical comparisons of methods", 
"","","","","& other two methods","ciBAx","Use of general CC and adj. factor",
"10","109","16","Tuyl 2008","Bayesian method with","ciBAx","Other classical methods; Numerical", 
"","","","","five different beta priors","","& graphical comparisons of methods",
"","","","","","","Use of general CC and adj. factor",
"11","NA","10","Vos 2005","p-confidence, p-bias","pCOpBIBA",""),byrow=TRUE, nrow=16,  ncol=7)

colnames(Key.Papers)=c("#","x","n","Paper","Methods","Function","Additional options")
#rownames(Key.Papers)=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15")

knitr::kable(Key.Papers,caption = "Additional functions")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 13. Paper 1
# Table 1 in document
library(proportion)
Paper1.1=ciAllx(x=0,n=20,alp=0.05)

knitr::kable(Paper1.1,caption = "Asymptotic methods CI using ciAllx(x=0,n=20,alp=0.05)")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 14. Paper 1
Paper1.2=ciEXx(x=0,n=20,alp=0.05,e=c(0.1,0.5,0.95,1)) 

knitr::kable(Paper1.2,caption = "Exact method CI using ciBAx(x=0,n=20,alp=0.05,e=c(0.1,0.5,0.95,1))")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 14. Paper 1
Paper1.31=ciBAx(x=0,n=20,alp=0.05,a=2,b=2) 
Paper1.32=ciBAx(x=0,n=20,alp=0.05,a=1,b=1) 
Paper1.33=ciBAx(x=0,n=20,alp=0.05,a=0.5,b=0.5) 
Paper1.34=ciBAx(x=0,n=20,alp=0.05,a=0.02,b=2) 
Paper1.31$desc = "Assuming Symmetry"
Paper1.32$desc = "Flat"
Paper1.33$desc = "Jeffreys"
Paper1.34$desc = "Near boundary"
Paper1.3.raw=rbind(Paper1.31,Paper1.32,Paper1.33,Paper1.34)

P.reorder=data.frame(Desc=Paper1.3.raw$desc,x=Paper1.3.raw$x,LBAQx=Paper1.3.raw$LBAQx,
                     UBAQx=Paper1.3.raw$UBAQx,LBAHx=Paper1.3.raw$LBAHx,UBAHx=Paper1.3.raw$UBAHx)

knitr::kable(P.reorder,caption = "Bayesian CI using ciBAx() with x=0,n=20,alp=0.05, varying a(2,1,0.05,0.02 and b(2,1,0.05,2)")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 15. Paper 1
Paper1.4=ciAAllx(x=0,n=20,alp=0.05,h=2)

knitr::kable(Paper1.4,caption = "Adding Pseudo constant using ciAAllx(x=0,n=20,alp=0.05,h=2)")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 16. Paper 1
Paper1.5=ciCAllx(x=0,n=20,alp=0.05,c=1/40)

knitr::kable(Paper1.5,caption = "Adding Continuity Correction, c = 1/(2n) & using ciCAllx(x=0,n=20,alp=0.05,c=1/40)")

## ---- echo=TRUE, results='asis'-----------------------------------------------
# 17. Paper 1
PlotciAllx(x=0,n=20,alp=0.05)


## ---- echo=TRUE, results='asis'-----------------------------------------------
# 18. Plot of sum of length of exact method
PlotlengthEX(n=10,alp=0.05,e=c(0.1,0.5,0.95,1),a=1,b=1) 


## ---- echo=FALSE, results='asis'----------------------------------------------
# 19. Paper 1
Paper1.6=covpAll(n=250,alp=0.05,a=1,b=1,t1=0.93,t2=0.97) 
knitr::kable(Paper1.6,caption = "Coverage probability using covpAll()")
PlotcovpAll(n=250,alp=0.05,a=1,b=1,t1=0.93,t2=0.97)

## ---- echo=TRUE, results='asis'-----------------------------------------------
# 20. Paper 1
PlotcovpAll(n=10,alp=0.05,a=1,b=1,t1=0.93,t2=0.97)

## ---- echo=TRUE, results='asis'-----------------------------------------------
# 21. Paper 2 - display the function
PlotciAllg(n=10,alp=0.05)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # 20. Paper 2 - display the function
#  hypotestBAF4(n=10, th0=0.9, a0=1,b0=1,a1=0.5,b1=0.5)

## ---- echo=FALSE, results='asis'----------------------------------------------
# 21. Paper 2
Paper2.1=hypotestBAF4(n=10, th0=0.9, a0=1,b0=1,a1=0.5,b1=0.5)
knitr::kable(Paper2.1,caption = "Hypothesis test, H0: p <= 0.9 vs. H1: p > 0.9")

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # 20. Function to evaluate ci varying the adding constant h
#  ciAAllx(x=16, n=17,alp = 0.05,h=0)
#  ciAAllx(x=16, n=109,alp = 0.05,h=0)

## ---- echo=FALSE, results='asis'----------------------------------------------
# 21. Paper 3&4
Paper2.2=ciAAllx(x=16, n=17,alp = 0.05,h=0)
Paper2.3=ciAAllx(x=16, n=109,alp = 0.05,h=0)
Paper2.4=ciAAllx(x=16, n=17,alp = 0.05,h=1)
Paper2.5=ciAAllx(x=16, n=109,alp = 0.05,h=1)
Paper2.6=ciAAllx(x=16, n=17,alp = 0.05,h=2)
Paper2.7=ciAAllx(x=16, n=109,alp = 0.05,h=2)
knitr::kable(Paper2.2,caption = "CI with x=16, n=17 & h=0")
knitr::kable(Paper2.3,caption = "CI with x=16, n=109 & h=0")
knitr::kable(Paper2.4,caption = "CI with x=16, n=17 & h=1")
knitr::kable(Paper2.5,caption = "CI with x=16, n=109 & h=1")
knitr::kable(Paper2.6,caption = "CI with x=16, n=17 & h=2")
knitr::kable(Paper2.7,caption = "CI with x=16, n=109 & h=2")

## ---- echo=TRUE, results='asis'-----------------------------------------------
# 21. Paper 3&4 - Plot of all the adjusted CI with h=1
PlotciAAllxg(x=16,n=17,alp=0.05,h=1)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # 22. Paper 3&4 - display the function
#  ciEXx(x=98, n=100,alp = 0.05,e=c(0.1,.5,0.95,1))

## ---- echo=FALSE, results='asis'----------------------------------------------
# 23. Paper 3&4
Paper3.2=ciEXx(x=98, n=100,alp = 0.05,e=c(0.1,.5,0.95,1))
knitr::kable(Paper3.2,caption = "CI-Exact with x=98, n=100")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 24. Paper 5
Paper5.1=probPREx(x=0,n=167,xnew=0,m=10,a1=1,a2=1)
Paper5.2=probPREx(x=0,n=167,xnew=0,m=50,a1=1,a2=1)
Paper5.3=probPREx(x=0,n=167,xnew=0,m=100,a1=1,a2=1)
Paper5.4=probPREx(x=0,n=167,xnew=0,m=150,a1=1,a2=1)
Paper5.5=probPREx(x=0,n=167,xnew=5,m=10,a1=1,a2=1)
Paper5.6=probPREx(x=0,n=167,xnew=25,m=50,a1=1,a2=1)
Paper5.7=probPREx(x=0,n=167,xnew=50,m=100,a1=1,a2=1)
Paper5.8=probPREx(x=0,n=167,xnew=75,m=150,a1=1,a2=1)
ndf=rbind(Paper5.1,Paper5.2,Paper5.3,Paper5.4,Paper5.5,Paper5.6,Paper5.7,Paper5.8)
knitr::kable(ndf,caption = "Predicted probability with x=0, n=167 varying xnew and m")

## ---- echo=FALSE, results='asis'----------------------------------------------
# 25. Paper 5
# Guidance for priors
Uniform.Prior=data.frame(Description="Uniform prior",a=1,b=1)
Jeffreys.Prior=data.frame(Description="Jeffreys prior",a=0.5,b=0.5)
Tuyl.sp1=data.frame(Description="Tuyl p1",a=0.042,b=27.96)
Tuyl.sp2=data.frame(Description="Tuyl p2",a=1,b=666)
Tuyl.sp3=data.frame(Description="Tuyl p3",a=1,b=398)
guidance.df=rbind(Uniform.Prior,Jeffreys.Prior,Tuyl.sp1,Tuyl.sp2,Tuyl.sp3)
knitr::kable(guidance.df,caption = "Guidance for priors used below")

# Data for table in paper
Pa.combo1.001= probPOSx(x=0,n=167,a=1,b=1,th=0.001)
Pa.combo2.001= probPOSx(x=0,n=167,a=0.5,b=0.5,th=0.001)
Pa.combo3.001= probPOSx(x=0,n=167,a=0.042,b=27.96,th=0.001)
Pa.combo4.001= probPOSx(x=0,n=167,a=1,b=666,th=0.001)
Pa.combo5.001= probPOSx(x=0,n=167,a=1,b=398,th=0.001)
Pa.combo1.01= probPOSx(x=0,n=167,a=1,b=1,th=0.01)
Pa.combo2.01= probPOSx(x=0,n=167,a=0.5,b=0.5,th=0.01)
Pa.combo3.01= probPOSx(x=0,n=167,a=0.042,b=27.96,th=0.01)
Pa.combo4.01= probPOSx(x=0,n=167,a=1,b=666,th=0.01)
Pa.combo5.01= probPOSx(x=0,n=167,a=1,b=398,th=0.01)
Pa.combo1.05= probPOSx(x=0,n=167,a=1,b=1,th=0.05)
Pa.combo2.05= probPOSx(x=0,n=167,a=0.5,b=0.5,th=0.05)
Pa.combo3.05= probPOSx(x=0,n=167,a=0.042,b=27.96,th=0.05)
Pa.combo4.05= probPOSx(x=0,n=167,a=1,b=666,th=0.05)
Pa.combo5.05= probPOSx(x=0,n=167,a=1,b=398,th=0.05)
Pa.combo1.1= probPOSx(x=0,n=167,a=1,b=1,th=0.1)
Pa.combo2.1= probPOSx(x=0,n=167,a=0.5,b=0.5,th=0.1)
Pa.combo3.1= probPOSx(x=0,n=167,a=0.042,b=27.96,th=0.1)
Pa.combo4.1= probPOSx(x=0,n=167,a=1,b=666,th=0.1)
Pa.combo5.1= probPOSx(x=0,n=167,a=1,b=398,th=0.1)
df.001=rbind(Pa.combo1.001,Pa.combo2.001,Pa.combo3.001,Pa.combo4.001,Pa.combo5.001)
df.01=rbind(Pa.combo1.01,Pa.combo2.01,Pa.combo3.01,Pa.combo4.01,Pa.combo5.01)
df.05=rbind(Pa.combo1.05,Pa.combo2.05,Pa.combo3.05,Pa.combo4.05,Pa.combo5.05)
df.1=rbind(Pa.combo1.1,Pa.combo2.1,Pa.combo3.1,Pa.combo4.1,Pa.combo5.1)
df001=t(df.001)
df01=t(df.01)
df1=t(df.1)
df05=t(df.05)
tr.df.001=df001[2,]
tr.df.01=df01[2,]
tr.df.1=df1[2,]
tr.df.05=df05[2,]
te.df=rbind(tr.df.001,tr.df.01,tr.df.1,tr.df.05)
rownames(te.df) <- c("th=0.001","th=0.01","th=0.1","th=0.5")
colnames(te.df) <- c("Uniform Prior","Jeffreys prior","Tuyl p1", "Tuyl  p2", "Tuyl p3")
knitr::kable(te.df,caption = "Posterior probability with x=0, n=167 varying th")


## ---- echo=FALSE, results='asis'----------------------------------------------
# 26. Paper 6
Paper6= pCOpBIBA(n=10,alp=0.05,a1=1,a2=1) 
knitr::kable(Paper6,caption = "p-Confidence & p-Bias of Bayesian method for n=10, a1=a2=1")

