#'Boundary crossing probabilities computation using exact test.
#'
#'Calculate boundary crossing probabilities of single-arm group sequential design with binary endpoint using binomial distribution
#'
#'This function is similar to \code{\link{asymprob}} except that the former uses binomial distribution and the latter
#'uses the normal asymptotic distribution. With \code{K=0}
#'(as default), \code{d} must be an object of class exactdesign. Meanwhile, other
#'arguments except for \code{p_1} will be inherited from \code{d} and the input values will be
#'ignored. With \code{K!=0}, the probabilities are derived from the input arguments. In
#'this circumstance, all the arguments except for \code{d} are required.
#'
#'The computation is based on the single-arm group sequential exact test
#'described in \code{\link{exactdesign}}. Therefore, for the output matrix of
#'upper bound crossing probabilities, the values for the first K-1 analyses are
#'zero since there is only one upper bound for the last analysis.
#'
#' @param K The maximum number of analyses, including the interim and the final. Should be an integer within (1,20]. K will be
#' rounded to the nearest whole number if it is not an integer. The default is 0.
#' @param p_0 The response rate or the probability of success under null hypothesis. Should be a scalar within (0,1).
#' @param p_1 A scalar or vector representing response rate or probability of success under the alternative hypothesis. The
#' value(s) should be within (p_0,1). It is a mandatory input.
#' @param n.I A vector of length K which contains sample sizes required at each analysis. Should be a positive and increasing
#' sequence.
#' @param u_K The upper boundary for the last analysis.
#' @param lowerbounds Non-decreasing lower boundaries for each analysis, in which each element is no less than -1 (no lower bound). With length K,
#' the last lower bound must be identical to u_K. With length K-1, the last element must be no greater than u_K and u_K will
#' be automatically added into the sequence. Note the lower bound must be less than the corresponding sample size.
#' @param d An object of the class exactdesign.
#'
#' @return An object of the class exactprob. This class contains:
#' \itemize{
#'   \item{p_0: As input with \code{d=NULL} or as in \code{d}.}
#'   \item{p_1: As input.}
#'   \item{K: K used in computation.}
#'   \item{n.I: As input with \code{d=NULL} or as in \code{d}.}
#'   \item{u_K: As input with \code{d=NULL} or as in \code{d}.}
#'   \item{lowerbounds: lowerbounds used in computation.}
#'   \item{problow: Probabilities of crossing the lower bounds at each analysis.}
#'   \item{probhi: Probability of crossing the upper bounds at each analysis.}
#' }
#'
#'@note The calculation of boundary crossing probabilities here borrowed strength from the
#'source code of function \code{gsBinomialExact} in package gsDesign and we really appreciate
#'their work.
#'
#'@section Reference: \itemize{ \item{Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
#'  Applications to Clinical Trials. Chapman and Hall/CRC, Boca Raton, FL, 2000.}
#'  \item{Keaven M. Anderson, Dan (Jennifer) Sun, Zhongxin (John) Zhang. gsDesign: An R
#'  Package for Designing Group Sequential Clinical Trials. R package version 3.0-1. }}
#'
#'@seealso \code{\link{exactdesign}}, \code{\link{exactcp}}, \code{\link{asymprob}}.
#'
#' @export
#'
#' @examples
#' I=c(0.2,0.4,0.6,0.8,0.99)
#' beta=0.2
#'betaspend=c(0.1,0.2,0.3,0.3,0.2)
#'alpha=0.05
#'p_0=0.3
#'p_1=0.5
#'K=4.6
#'tol=1e-6
#'tt1=asymdesign(I,beta,betaspend,alpha,p_0,p_1,K,tol)
#'tt2=exactdesign(tt1)
#'tt3=exactprob(p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),d=tt2)
#'tt3=exactprob(K=5,p_0=0.4,p_1=c(0.5,0.6,0.7,0.8),n.I=c(15,20,25,30,35),u_K=15,
#'lowerbounds=c(3,5,10,12,15))
exactprob<-function(K=0,p_0,p_1,n.I,u_K,lowerbounds,d=NULL){
  ##check validity of inputs
  if(K==0){##so the user input wrong parameter values,check if d belongs to class "exactdesign"
    if(!methods::is(d,"exactdesign"))
	  stop('d is not an object of class exactdesign')
	##if d is exactdesign class, then adopt the parameters except for p_1
	p_0=d$p_0
	if((min(p_1)<=p_0)|(max(p_1)>=1))
      stop('Please input p_1 that lies between p_0 and 1 (not including p_0 and 1).')
	K=d$K
	n.I=d$n.I
	u_K=d$u_K
	lowerbounds=d$lowerbounds
  }else{##in this case, we need to check the validity of inputs and adopt the inputs by the user.
    if(min(lowerbounds)<(-1))
      stop('Lowerbounds must be no less than -1.')
    u_K=round(u_K)
    temp1=check.prob(K,p_0,p_1,n.I,u_K,round(lowerbounds))
	K=temp1$K
	lowerbounds=temp1$lowerbounds
	n.I=temp1$n.I
	temp2=n.I-lowerbounds
	if(min(temp2)<=0)
	  stop('lower bound must be less than the corresponding sample size.')
  }

  probhi=matrix(0,1+length(p_1),K)
  problow=probhi

  ##compute boundary crossing probabilities under H0 and put them in the first row.
  temp1=exactprob1(n.I,lowerbounds,p_0,K,K)
  problow[1,]=temp1$plo
  probhi[1,K]=temp1$phi

  ##compute boundary crossing probabilities under H1
  for(i in 1:length(p_1)){
    temp1=exactprob1(n.I,lowerbounds,p_1[i],K,K)
    problow[(i+1),]=temp1$plo
    probhi[(i+1),K]=temp1$phi
  }

  problow=cbind(c(p_0,p_1),problow,rowSums(problow))
  probhi=cbind(c(p_0,p_1),probhi)
  colnames(problow)=c('p',1:K,'Total')
  colnames(probhi)=c('p',1:K)
  x=list(p_0=p_0,p_1=p_1,K=K,n.I=n.I,u_K=u_K,lowerbounds=lowerbounds,problow=problow,probhi=probhi)
  class(x)='exactprob'
  return(x)
}

