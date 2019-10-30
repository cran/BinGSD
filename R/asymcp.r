#'Conditional power computation using asymptotic test.
#'
#' Compute conditional power of single-arm group sequential design with binary endpoint based on asymptotic test, given the interim
#' result.
#'
#'Conditional power quantifies the conditional probability of crossing the upper bound given the interim result \eqn{z_i},
#'\eqn{1\le i<K}. Having inherited sample sizes and boundaries from \code{\link{asymdesign}} or \code{\link{asymprob}},
#'given the interim statistic at \eqn{i}th analysis \eqn{z_i}, the conditional power is defined as
#'
#'\eqn{\alpha _{i,K}(p|z_i)=P_{p}(Z_K\ge u_K, Z_{K-1}>l_{K-1}, \ldots, Z_{i+1}>l_{i+1}|Z_i=z_i)}
#'
#'With asymptotic test, the test
#'statistic at analysis \eqn{k} is
#'\eqn{Z_k=\hat{\theta}_k\sqrt{n_k/p/(1-p)}=(\sum_{s=1}^{n_k}X_s/n_k-p_0)\sqrt{n_k/p/(1-p)}},
#'which follows the normal distribution \eqn{N(\theta \sqrt{n_k/p/(1-p)},1)}
#'with \eqn{\theta=p-p_0}. In practice, \eqn{p} in \eqn{Z_k} can be substituted
#'with the sample response rate \eqn{\sum_{s=1}^{n_k}X_s/n_k}.
#'
#'The increment statistic \eqn{Z_k\sqrt{n_k/p/(1-p)}-Z_{k-1}\sqrt{n_{k-1}/p/(1-p)}} also follows a normal distribution independently
#'of \eqn{Z_{1}, \ldots, Z_{k-1}}. Then the conditional power can be easily obtained using a procedure similar
#'to that for unconditional boundary crossing probabilities.
#'
#' @param d An object of the class asymdesign or asymprob.
#' @param p_1 A scalar or vector representing response rate or probability of success under the alternative hypothesis. The
#' value(s) should be within (p_0,1).
#' @param i Index of the analysis at which the interim statistic is given. Should be an integer ranges from 1 to K-1. i will be
#' rounded to its nearest whole value if it is not an integer.
#' @param z_i The interim statistic at analysis i.
#'
#' @return A list with the elements as follows:
#' \itemize{
#'   \item{K: As in d.}
#'   \item{n.I: As in d.}
#'   \item{u_K: As in d.}
#'   \item{lowerbounds: As in d.}
#'   \item{i: i used in computation.}
#'   \item{z_i: As input.}
#'   \item{cp: A matrix of conditional powers under different response rates.}
#'   \item{p_1: As input.}
#'   \item{p_0: As input.}
#' }
#'
#'
#'@section Reference: \itemize{
#'  \item{Alan Genz et al. (2018). mvtnorm: Multivariate Normal and t Distributions. R package version 1.0-11.}}
#'
#'@seealso \code{\link{asymprob}}, \code{\link{asymdesign}},
#'  \code{\link{exactcp}}.
#'
#' @export
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
#'tt2=asymprob(p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),d=tt1)
#' asymcp(tt1,p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),1,2)
#' asymcp(tt2,p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),3,2.2)
asymcp<-function(d,p_1,i,z_i){
 ##check validity of inputs
 if(!(methods::is(d,"asymdesign")|methods::is(d,"asymprob")))
    stop("asymcp must be called with class of either asymdesign or asymprob.")
 p_0=d$p_0
 i=round(i)
 if((min(p_1)<=p_0)|(max(p_1)>=1))
    stop('Please input p_1 that lies between p_0 and 1 (not including p_0 and 1).')
 if((i>=K)|(i<1))
    stop('i is less than 1 or more than K-1.')
 if(abs(z_i)>=100)
    stop('invalid interim statistic: too large or too small.')
 K=d$K-i  ##the no. of analysis after the ith analysis.
 n.I=(d$n.I)[i:(i+K)]  ##sample sizes for i,i+1,...,K analysis
 lowerbounds=(d$lowerbounds)[(i+1):(i+K)]  ##lower bounds for i+1,i+2,...K analysis.
 cp=matrix(0,1+length(p_1),1) ##to store conditional power under different p

 #conditional power under H0
 infor=n.I/p_0/(1-p_0)  ##sequence of Fisher information under H0, for i, i+1,...K analysis.
 lowerbounds1=lowerbounds*sqrt(infor[-1])-z_i*sqrt(infor[1]) ##new lower bounds under H0,l_m^*
 sigma=matrix(0,K,K) #the covariance matrix of multivariate normal distribution, K here is essentially K-i.
 for(j1 in 1:K){
    for(j2 in 1:K){
        sigma[j1,j2]=infor[(min(j1,j2)+1)]  ##I_k, k is the smaller between j1 and j2.
	}
  }
 sigma=sigma-infor[1]   ##the covariance matrix of the joint normal dist. of Z_{i,i+1},Z_{i,i+2},...,Z_{i,K}
 cp[1,]=mvtnorm::pmvnorm(lower=lowerbounds1,upper=rep(Inf,K),mean=rep(0,K),sigma=sigma)[1] ##cp under H0

 ##conditional power under H1.
 for(s in 1:length(p_1)){
    p1=p_1[s]
    infor=n.I/p1/(1-p1)  ##sequence of Fisher information under H1, for i, i+1,...K analysis.
    lowerbounds1=lowerbounds*sqrt(infor[-1])-z_i*sqrt(infor[1]) ##new lower bounds under H1,l_m^*
    sigma=matrix(0,K,K) #the covariance matrix of multivariate normal distribution, K here is essentially K-i.
    for(j1 in 1:K){
       for(j2 in 1:K){
          sigma[j1,j2]=infor[(min(j1,j2)+1)]  ##I_k, k is the smaller between j1 and j2.
	    }
    }
    sigma=sigma-infor[1]   ##the covariance matrix of the joint normal dist. of Z_{i,i+1},Z_{i,i+2},...,Z_{i,K}
	mean1=(p1-p_0)*(infor[-1]-infor[1])
    cp[(s+1),]=mvtnorm::pmvnorm(lower=lowerbounds1,upper=rep(Inf,K),mean=mean1,sigma=sigma)[1] ##cp under H0
 }

 cp=cbind(c(p_0,p_1),cp)
 colnames(cp)=c('p','cp')
 return(list(K=d$K,n.I=d$n.I,u_K=d$u_K,lowerbounds=d$lowerbounds,i=i,z_i=z_i,cp=cp,p_1=p_1,p_0=p_0))
}
