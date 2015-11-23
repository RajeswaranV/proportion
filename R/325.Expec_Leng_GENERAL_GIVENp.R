#' Plot the expected length given hypothetical "p"
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param hp - Hypothetical "p"
#' @details  The  plot of the expected length for \code{n} given lower limit \code{LL} and  upper limit \code{UL}
#' @family Expected length
#' @examples
#' n= 5;
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' hp=seq(0,1,by=0.01)
#' PlotexplGEN(n,LL,UL,hp)
#' @export
##### 1.Expected Length - Graph
PlotexplGEN<-function(n,LL,UL,hp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(hp)) stop("'hp' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")
  if (any(hp>1) || any(hp<0)) stop("'hp' has to be between 0 and 1")

####INPUT n
x=0:n
k=n+1
s=length(hp)
ewi=matrix(0,k,s)						#Expected length quantity in sum
ew=0									#Expected Length
LE=0

for(i in 1:k)
{
LE[i]=UL[i]-LL[i]
}

####Expected Length

for (j in 1:s)
{
for(i in 1:k)
{
ewi[i,j]=LE[i]*dbinom(i-1, n,hp[j])
}
ew[j]=sum(ewi[,j])						#Expected Length
}
EL=data.frame(hp,ew)
ggplot2::ggplot(EL, ggplot2::aes(x=hp, y=ew))+
  ggplot2::labs(title = "Expected length given hypothetical 'p'") +
  ggplot2::labs(y = "Expected length") +
  ggplot2::labs(x = "p") +
  ggplot2::geom_vline(ggplot2::aes(xintercept=0.5,color="Cutoff"), linetype=2)+
  ggplot2::geom_line(ggplot2::aes(color="Expected length"))+
  ggplot2::geom_point(ggplot2::aes(color="EL Values"))+
  ggplot2::scale_colour_manual(name='Heading',
                               values=c('Cutoff'="brown",
                                        'Expected length'='red',
                                        'EL Values'='red'),
                               guide='legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(2,1,1),
                                                                     shape=c(NA, 16,NA))))

}
