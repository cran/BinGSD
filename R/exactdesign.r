#' Compute sample size and boundaries using exact binomial distribution
#'
#' Compute sample size and boundaries of single-arm group sequential design with binary endpoint using exact binomial distribution
#'
#'Suppose \eqn{X_{1}, X_{2}, \ldots} are binary outcomes following Bernoulli
#'distribution \eqn{b(1,p)}, in which 1 stands for the case that the subject
#'responds to the treatment and 0 otherwise. Consider a group sequential test
#'with \eqn{K} planned analyses, where the null and alternative hypotheses are
#'\eqn{H_0: p=p_0} and \eqn{H_1: p=p_1} respectively. Note that generally
#'\eqn{p_1} is greater than \eqn{p_0}. For \eqn{k<K}, the trial stops if and
#'only if the test statistic \eqn{Z_k} crosses the futility boundary, that is,
#'\eqn{Z_k<=l_k}. The lower bound for the last analysis \eqn{l_K} is set to be
#'equal to the last and only upper bound \eqn{u_K} to make a decision. At the
#'last analysis, the null hypothesis will be rejected if \eqn{Z_K>=u_K}.
#'
#'The computation of lower bounds except for the last one is implemented with
#'\eqn{u_K} fixed, thus the derived lower bounds are non-binding. Furthermore,
#'the overall type I error will not be inflated if the trial continues after
#'crossing any of the interim lower bounds, which is convenient for the purpose
#'of monitoring. Let the sequence of sample sizes required at each analysis be
#'\eqn{n_{1}, n_{2}, \ldots, n_{K}}. For binomial endpoint, the Fisher
#'information equals \eqn{n_k/p/(1-p)} which is proportional to \eqn{n_k}.
#'Accordingly, the information fraction available at each analysis is equivalent
#'to \eqn{n_k/n_K}.
#'
#'With exact test, the test statistic at analysis \eqn{k} is \eqn{Z_k=\sum_{s=1}^{n_k}X_s}
#'which follows binomial distribution \eqn{b(n_k,p)}. Actually, \eqn{Z_k} is the total
#'number of responses up to the kth analysis.
#'
#'Under the null hypothesis, \eqn{Z_k} follows a binomial distribution \eqn{b(n_k,p_0)}.
#'While under the alternative hypothesis, \eqn{Z_k} follows \eqn{b(n_k,p_1)}.
#'It may involve massive computation to simutaneously find proper \eqn{n_K} and \eqn{u_K}.
#'In fact, the sample sizes obtained from asymptotic test ought to be close to those from exact test.
#'Thus, we adopt \eqn{n_K} from asymptotic test as the starting value. The starting value of \eqn{u_K} is
#'computed given the \eqn{n_K}. Iteratively update \eqn{u_K} and \eqn{n_K} until errors are limited to
#'certain amount.
#'
#'Like \code{\link{asymdesign}}, the lower boundaries for the first \eqn{K-1} analyses are
#'sequentially determined by a search method. However, if the actual overall type II error exceeds the desired level,
#'not only sample sizes but also all the boundaries are updated, since the binomial distribution under \eqn{H_0}
#'involves with sample size.
#'
#'Due to the discreteness of binomial distribution, in exact test, the type I and
#'type II error actually spent at each analysis may not approximate the designated
#'amount. With the only one upper bound, the whole type I error is spent at the final analysis.
#' From some simulation studies, though not presented here, we found that carrying over
#'unused type II error has minor influence on the resulting boundaries and sample sizes.
#'However, in an attempt to reduce the false positive rate, we decided to recycle the unspent
#' amount of desired type II error. Thus, the elements of betaspend in an exactdesign object may be greater than
#' the amount pre-specified by the user.
#'
#' @param d An object of the class asymdesign.
#'
#' @return An object of the class exactdesign. This class contains:
#' \itemize{
#'   \item{I: I used in computation, as in d.}
#'   \item{beta: The desired overall type II error level, as in d.}
#'   \item{betaspend: The desired type II error spent at each analysis used in computation, as in d.}
#'   \item{alpha: The desired overall type I error level, as in d.}
#'   \item{p_0: The response rate or the probability of success under null hypothesis, as in d.}
#'   \item{p_1: The response rate or the probability of success under alternative hypothesis, as in d.}
#'   \item{K: K used in computation, as in d.}
#'   \item{n.I: A vector of length K which contains sample sizes required at each analysis to achieve desired type I and type
#'   II error requirements. n.I equals sample size for the last analysis times the vector of information fractions.}
#'   \item{u_K: The upper boundary for the last analysis.}
#'   \item{lowerbounds: A vector of length K which contains lower boundaries for each annalysis. Note that the lower
#'   boundaries are non-binding.}
#'   \item{problow: Probabilities of crossing the lower bounds under \eqn{H_1} or the actual type II error at each analysis.}
#'   \item{probhi: Probability of crossing the last upper bound under \eqn{H_0} or the actual type I error.}
#'   \item{power: power of the group sequential test with the value euqals 1-sum(problow).}
#' }
#'
#'
#'@section Reference: \itemize{ \item{Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
#'  Applications to Clinical Trials. Chapman and Hall/CRC, Boca Raton, FL, 2000.} }
#'
#'@seealso \code{\link{exactprob}}, \code{\link{exactcp}},
#'  \code{\link{asymdesign}}.
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
exactdesign<-function(d){
  if(!methods::is(d,"asymdesign"))
	stop('d is not an object of class asymdesign')
  I=d$I
  beta=d$beta
  betaspend=d$betaspend  ##the desired type II error spent at each analysis
  alpha=d$alpha
  p_0=d$p_0
  p_1=d$p_1
  K=d$K
  n_K=d$n.I[K]    ##adopt the maximum sample size output by asymdesign as starting value.
  s=0
  lowerbounds=rep(0,K) ##to store lower bounds
  problow=rep(0,K) ##to store actual lower bound crossing probabilities under H1

  while(s<=20){
    ##get starting value of u_K and check if the maximum sample size satisfies the two types of error requirements
    u_K=stats::qbinom(1-alpha,n_K,p_0)+1 ##starting value of u_K under H0
    temp1=stats::pbinom(u_K-1,n_K,p_1)
    t=0
    while((temp1>beta)&(t<=20)){  ##Step 0 of Algorithm 2
      n_K=n_K+1
      u_K=stats::qbinom(1-alpha,n_K,p_0)+1
	  temp1=stats::pbinom(u_K-1,n_K,p_1)
      t=t+1
    }
    if(temp1>beta)
      stop('cannot converge with the current setting of inputs. You may consider a larger beta or alpha.')

    n.I=ceiling(n_K*I)  ##the vector of sample sizes for all stages.

    temp1=c(stats::qbinom(betaspend[1],n.I[1],p_1)-1,stats::qbinom(betaspend[1],n.I[1],p_1))
	temp2=stats::pbinom(temp1,n.I[1],p_1)
	lowerbounds[1]=max(temp1[temp2<=betaspend[1]])  ##the first lower bound.
	if(lowerbounds[1]>=u_K)
	  stop('the first lower bound is larger than or equal to the last upper bound.')
	problow[1]=max(temp2[temp2<=betaspend[1]])  ##the actual type II error for first analysis
	betaspend[2]=betaspend[2]+betaspend[1]-problow[1] ##to carry over the unspent type II error

	##for design with more than 2 stages.
    if(K>2){
      for(k in 2:(K-1)){
	    temp1=bound2(k,lowerbounds[1:(k-1)],u_K,n.I[1:k],p_1,betaspend[k],K) ##get the lower bound for kth analysis.
		lowerbounds[k]=temp1$l_k
		problow[k]=temp1$error
		betaspend[(k+1)]=betaspend[(k+1)]+betaspend[k]-problow[k]  ##carryover the unspent type II error.
      }
    }
    ##make the last lower bound=u_K and calculate actual power
    lowerbounds[K]=u_K
	problow[K]=(exactprob1(n.I,lowerbounds,p_1,K,K)$plo)[K]

    ##check if the design satisfies the power constraint. if not ,increase the sample size.
	if(beta>=sum(problow))
       break
    ##if the actual power less than desired, increase the maximum sample size by 1.
    n_K=n_K+1	##increase sample size.
    s=s+1
 }

 if(beta<sum(problow)) ##the last while iteration is stopped due to s>20 but not convergence
    stop('cannot converge with the current tol.')

 probhi=exactprob1(n.I,lowerbounds,p_0,K,K)$phi ##the value is the probability of crossing the upper bound under H0.

 x=list(I=I,beta=beta,betaspend=d$betaspend,alpha=alpha,p_0=p_0,p_1=p_1,K=K,n.I=n.I,u_K=u_K,lowerbounds=lowerbounds,
   problow=problow,probhi=probhi,power=1-sum(problow))
 class(x)="exactdesign"

 return(x)
}
