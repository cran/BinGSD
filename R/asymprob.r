#'Boundary crossing probabilities computation using asymptotic test.
#'
#'Calculate boundary crossing probabilities of single-arm group sequential
#'design with binary endpoint based on asymptotic test.
#'
#'This function calculates probabilities of crossing the upper or the lower
#'boundaries under null hypothesis and a set of alternative hypothese. With \code{K=0}
#'(as default), d must be an object of class asymdesign. Meanwhile, other
#'arguments except for \code{p_1} will be inherited from \code{d} and the input values will be
#'ignored. With \code{K!=0}, the probabilities are derived from the input arguments. In
#'this circumstance, all arguments except for \code{d} are required.
#'
#'The computation is based on the single-arm group sequential asymptotic test
#'described in \code{\link{asymdesign}}. Therefore, for the output matrix of
#'upper bound crossing probabilities, the values for the first K-1 analyses are
#'zero since there is only one upper bound for the last analysis.
#'
#'@param K The maximum number of analyses, including the interim and the final.
#'  Should be an integer within (1,20]. K will be rounded to its nearest whole
#'  number if it is not an integer. The default is 0.
#'@param p_0 The response rate or the probability of success under null
#'  hypothesis. Should be a scalar within (0,1).
#'@param p_1 A scalar or vector representing response rate or probability of
#'  success under the alternative hypothesis. The value(s) should be within
#'  (p_0,1). It is a mandatory input.
#'@param n.I A vector of length K which contains sample sizes required at each
#'  analysis. Should be a positive and increasing sequence.
#'@param u_K The upper boundary for the last analysis.
#'@param lowerbounds Non-decreasing lower boundaries for each analysis. With
#'  length K, the last lower bound must be identical to u_K. With length K-1,
#'  the last element must be no greater than u_K and u_K will be automatically
#'  added into the sequence.
#'@param d An object of the class asymdesign.
#'
#'@return An object of the class asymprob. This class contains: \itemize{
#'  \item{p_0: As input with \code{d=NULL} or as in \code{d}.} \item{p_1: As input.} \item{K:
#'  K used in computation.} \item{n.I: As input with \code{d=NULL} or as in \code{d}.}
#'  \item{u_K: As input with \code{d=NULL} or as in \code{d}.} \item{lowerbounds: lowerbounds
#'  used in computation.} \item{problow: Probabilities of crossing the lower
#'  bounds at each analysis.} \item{probhi: Probability of
#'  crossing the upper bounds at each analysis.} }
#'
#'
#'@section Reference: \itemize{ \item{Alan Genz et al. (2018). mvtnorm:
#'  Multivariate Normal and t Distributions. R package version 1.0-11.}}
#'
#'@seealso \code{\link{asymdesign}}, \code{\link{asymcp}}, \code{\link{exactprob}}.
#'
#'@export
#'
#' @examples
#' I=c(0.2,0.4,0.6,0.8,0.99)
#'beta=0.2
#'betaspend=c(0.1,0.2,0.3,0.3,0.2)
#'alpha=0.05
#'p_0=0.3
#'p_1=0.5
#'K=4.6
#'tol=1e-6
#'tt1=asymdesign(I,beta,betaspend,alpha,p_0,p_1,K,tol)
#'asymprob(p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),d=tt1)
#'asymprob(K=5,p_0=0.4,p_1=c(0.5,0.6,0.7,0.8),n.I=c(15,20,25,30,35),u_K=1.65,
#'lowerbounds=c(-1.2,-0.5,0.2,0.8,1.65))
asymprob<-function(K=0,p_0,p_1,n.I,u_K,lowerbounds,d=NULL){
  ##check validity of inputs
  if(K==0){##so the user input wrong parameter values,check if d belongs to class "asymdesign"
    if(!methods::is(d,"asymdesign"))
	  stop('d is not an object of class asymdesign')
	##if d is asymdesign class, then adopt the parameters except for p_1
	p_0=d$p_0
	if((min(p_1)<=p_0)|(max(p_1)>=1))
      stop('Please input p_1 that lies between p_0 and 1 (not including p_0 and 1).')
	K=d$K
	n.I=d$n.I
	u_K=d$u_K
	lowerbounds=d$lowerbounds
  }else{##in this case, we need to check the validity of inputs and adopt the inputs by the user.
    temp1=check.prob(K,p_0,p_1,n.I,u_K,lowerbounds)
	K=temp1$K
	lowerbounds=temp1$lowerbounds
	n.I=temp1$n.I
  }

  probhi=matrix(0,1+length(p_1),K)
  problow=probhi

  ##compute boundary crossing probabilities under H0 and put them in the first row.
  problow[1,]=asymprob2(n.I,lowerbounds,K)
  probhi[1,K]=1-sum(problow[1,])

  ##compute boundary crossing probabilities under H1
  for(i in 1:length(p_1)){
    problow[(i+1),]=asymprob1(n.I,lowerbounds,p_0,p_1[i],K)
  }
  probhi[2:(1+length(p_1)),K]=1-rowSums(problow[2:(1+length(p_1)),])

  problow=cbind(c(p_0,p_1),problow,(1-probhi[,K]))
  probhi=cbind(c(p_0,p_1),probhi)
  colnames(problow)=c('p',1:K,'Total')
  colnames(probhi)=c('p',1:K)
  x=list(p_0=p_0,p_1=p_1,K=K,n.I=n.I,u_K=u_K,lowerbounds=lowerbounds,problow=problow,probhi=probhi)
  class(x)='asymprob'
  return(x)
}
