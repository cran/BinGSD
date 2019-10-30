#'Boundary and sample size computation using asymptotic test.
#'
#'Calculate boundaries and sample sizes of single-arm group sequential design
#'with binary endpoint based on asymptotic test.
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
#'For a \eqn{p_0} not close to 1 or 0, with a large sample size, the test
#'statistic at analysis \eqn{k} is
#'\eqn{Z_k=\hat{\theta}_k\sqrt{n_k/p/(1-p)}=(\sum_{s=1}^{n_k}X_s/n_k-p_0)\sqrt{n_k/p/(1-p)}},
#'which follows the normal distribution \eqn{N(\theta \sqrt{n_k/p/(1-p)},1)}
#'with \eqn{\theta=p-p_0}. In practice, \eqn{p} in \eqn{Z_k} can be substituted
#'with the sample response rate \eqn{\sum_{s=1}^{n_k}X_s/n_k}.
#'
#'Under the null hypothesis, \eqn{\theta=0} and \eqn{Z_k} follows a standard
#'normal distribution. During the calculation, the only upper bound \eqn{u_K} is
#'firstly derived under \eqn{H_0}, without given \eqn{n_K}. Thus, there is no
#'need to adjust \eqn{u_K} for different levels of \eqn{n_K}. Following East,
#'given \eqn{u_K}, compute the maximum sample size \eqn{n_K} under \eqn{H_1}.
#'The rest sample sizes can be obtained by multiplying information fractions
#'and \eqn{n_K}. The lower boundaries for the first \eqn{K-1} analyses are
#'sequentially determined by a search method. The whole searching procedure
#'stops if the overall type II error does not excess the desired level or the
#'times of iteration excess 30. Otherwise, increase the sample sizes until
#'the type II error meets user's requirement.
#'
#'The multiple integrals of multivariate normal density functions are conducted with
#'\code{\link[mvtnorm]{pmvnorm}} in R package mvtnorm. Through a few transformations of the integral variables,
#'\code{\link[mvtnorm]{pmvnorm}} turns the multiple integral to the product of several
#'univariate integrals, which greatly reduces the computational burden of sequentially searching for
#'appropriate boundaries.
#'
#'@param I The information fractions at each analysis. For binary endpoints, the
#'  information fraction for analysis k is equal to n_k/n_K, where n_k is the
#'  sample size available at analysis k and n_K is the sample size available at
#'  the last analysis or the maximum sample size. Should be a positive
#'  increasing vector of length K or K-1. If I has K elements among which the
#'  last one is not 1, then I will be standardized so that the last information
#'  fraction is 1. If I has K-1 elements, the last element in I must be less
#'  than 1.
#'@param beta The desired overall type II error level. Should be a scalar within
#'  the interval (0,0.5]. Default value is 0.3, that is, power=0.7.
#'@param betaspend The proportions of beta spent at each analysis. Should be a
#'  vector of length K with all elements belong to [0,1]. If the sum of all
#'  elements in betaspend is not equal to 1, betaspend will be standardized.
#'@param alpha The desired overall type I error level. Should be a scalar within
#'  the interval (0,0.3]. Default is 0.05.
#'@param p_0 The response rate or the probability of success under null
#'  hypothesis. Should be a scalar within (0,1).
#'@param p_1 The response rate or the probability of success under alternative
#'  hypothesis. Should be a scalar within (p_0,1).
#'@param K The maximum number of analyses, including the interim and the final.
#'  Should be an integer within (1,20]. K will be rounded to its nearest whole
#'  number if it is not an integer.
#'@param tol The tolerance level which is essentially the maximum acceptable difference between
#'  the desired type II error spending and the actual type II error spending, when
#'  computing the boundaries using asymptotic test. Should be a positive scalar no
#'  more than 0.01. The default value is 1e-6.
#'
#'@return An object of the class asymdesign. This class contains:
#'\itemize{
#'  \item{I: I used in computation.}
#'  \item{beta: As input.}
#'  \item{betaspend: The desired type II error spent at each analysis used in computation.}
#'  \item{alpha: As input.}
#'  \item{p_0: As input.}
#'  \item{p_1: As input.}
#'  \item{K: K used in computation.}
#'  \item{tol: As input.}
#'  \item{n.I: A vector of length
#'  K which contains sample sizes required at each analysis to achieve desired
#'  type I and type II error requirements. n.I equals sample size for the last
#'  analysis times the vector of information fractions.}
#'  \item{u_K: The upper boundary for the last analysis.}
#'  \item{lowerbounds: A vector of length K
#'  which contains lower boundaries for each analysis. Note that the lower
#'  boundaries are non-binding.}
#'  \item{problow: Probabilities of crossing the
#'  lower bounds under \eqn{H_1} or the actual type II error at each analysis.}
#'  \item{probhi: Probability of crossing the last upper bound under \eqn{H_0} or the
#'  actual type I error.} \item{power: power of the group sequential test with
#'  the value equals 1-sum(problow).} }
#'
#'
#'@section Reference: \itemize{ \item{Cytel Inc. East Version 6.4.1 Manual.
#'  2017.}
#'  \item{Alan Genz et al. (2018). mvtnorm: Multivariate Normal and t Distributions. R package version 1.0-11.}}
#'
#'@seealso \code{\link{asymprob}}, \code{\link{asymcp}},
#'  \code{\link{exactdesign}}.
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
#'
asymdesign<-function(I,beta=0.3,betaspend,alpha=0.05,p_0,p_1,K,tol=1e-6){#I is timing/information fractions
  ##-------check the validity of all inputs-----
  temp1=check.asymdesign(I,beta,betaspend,alpha,p_0,p_1,K,tol)
  I=temp1$I
  betaspend=temp1$betaspend
  K=temp1$K
  lowerbounds=rep(0,K)  ##the vector for storing the lower bounds for each analysis.
  betaspend=betaspend*beta ##the vector of type II error spent at each analysis

  ##-------obtain the upper bound u_K and the initial value sample sizes n.I and l_1----
  u_K=stats::qnorm(1-alpha)  #the upper bound for the last analysis.
  n_K=ceiling(p_1*(1-p_1)*(((u_K-stats::qnorm(beta))/(p_1-p_0))^2)) ##n_K for initialization.
  n.I=ceiling(n_K*I)  ##the vector of sample sizes for all stages.
  lowerbounds[1]=stats::qnorm(betaspend[1])+(p_1-p_0)*sqrt(n.I[1]/p_1/(1-p_1)) ##lower boundary for the first analysis
  problow=betaspend   ##the vector for storing actual type II error achieved at each analysis under H1.

  ##for design with more than 2 stages.
  if(K>2){
    for(k in 2:(K-1)){
      temp1=bound1(k,lowerbounds[1:(k-1)],u_K,n.I[1:k],p_1,p_0,betaspend[k],tol) ##get the lower bound for kth analysis.
      flag=temp1$flag
      if(flag){##if l_k=u_K, then make the rest of lower bounds be u_K.
        lowerbounds[k:K]=u_K
        break
      }
      lowerbounds[k]=temp1$l_k
      problow[k]=temp1$error
    }
  }

  ##make the last lower bound=u_K and calculate actual power
  lowerbounds[K]=u_K
  if(!flag){## so only the type II error at the last analysis needs to calculate
    mean1=(p_1-p_0)*sqrt(n.I/p_1/(1-p_1))  ##the mean vector of the multivariate normal distribution.
    lowerlimits=c(lowerbounds[1:(K-1)],-Inf) ##the lower limits for computing the boundary crossing probability.
    sigma=matrix(0,K,K) #the covariance matrix of multivariate normal distribution.
    for(i in 1:K){
      for(j in 1:K){
        sigma[i,j]=mean1[min(i,j)]/mean1[max(i,j)]
      }
    }
    problow[K]=mvtnorm::pmvnorm(lower=lowerlimits,upper=c(rep(Inf,(K-1)),u_K),mean=mean1,sigma=sigma)[1]
  }else{##call for gsprob function to compute the boundary crossing probabilities.
    problow=asymprob1(n.I,lowerbounds,p_0,p_1,K)
  }

  ##check if the design satisfies the power constraint. if not ,increase the sample size.
  t=0
  ##if the actual power less than desired, increase sample size.
  while((beta<sum(problow))&(t<=30)){
    n_K=n_K+1	##increase sample size
    n.I=ceiling(n_K*I)  ##the vector of sample sizes for all stages.
    problow=asymprob1(n.I,lowerbounds,p_0,p_1,K)
    t=t+1
  }

  if(beta<sum(problow)) ##the last while iteration is stopped due to t>30 but not convergence
    stop('cannot converge with the current tol.')

  probhi=asymprob2(n.I,lowerbounds,K) ##the value is the probabilities of crossing the lower bounds under H0.
  probhi=1-sum(probhi)  ##the type I error is equal to 1-sum(lower bound crossing probabities under H0)

  x=list(I=I,beta=beta,betaspend=betaspend,alpha=alpha,p_0=p_0,p_1=p_1,K=K,tol=tol,n.I=n.I,u_K=u_K,lowerbounds=lowerbounds,
         problow=problow,probhi=probhi,power=1-sum(problow))
  class(x)="asymdesign"

  return(x)
}

