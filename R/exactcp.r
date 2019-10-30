#'Conditional power computation using exact test.
#'
#'Compute conditional power of single-arm group sequential design with binary
#'endpoint based on binomial distribution.
#'
#'Conditional power quantifies the conditional probability of crossing the upper bound given the interim result \eqn{z_i},
#'\eqn{1\le i<K}. Having inherited sample sizes and boundaries from \code{\link{exactdesign}} or \code{\link{exactprob}},
#'given the interim statistic at \eqn{i}th analysis \eqn{z_i}, the conditional power is defined as
#'
#'\eqn{\alpha _{i,K}(p|z_i)=P_{p}(Z_K\ge u_K, Z_{K-1}>l_{K-1}, \ldots, Z_{i+1}>l_{i+1}|Z_i=z_i)}
#'
#'With exact test, the test statistic at analysis \eqn{k} is \eqn{Z_k=\sum_{s=1}^{n_k}X_s}
#'which follows binomial distribution \eqn{b(n_k,p)}. Actually, \eqn{Z_k} is the total
#'number of responses up to the kth analysis.
#'
#'The increment statistic \eqn{Z_k-Z_{k-1}} also follows a binomial distribution \eqn{b(n_k-n_{k-1},p)} independently
#'of \eqn{Z_{1}, \ldots, Z_{k-1}}. Then the conditional power can be easily obtained using the same procedure
#'for deriving unconditional boundary crossing probabilities.
#'
#'Note that \eqn{Z_{1}, \ldots, Z_{K}} is a non-decreasing sequence, thus the conditional power is 1 when the interim statistic
#' \eqn{z_i>=u_K}.
#'
#'@param d An object of the class exactdesign or exactprob.
#'@param p_1 A scalar or vector representing response rate or probability of
#'  success under the alternative hypothesis. The value(s) should be within
#'  (p_0,1).
#'@param i Index of the analysis at which the interim statistic is given. Should
#'  be an integer ranges from 1 to K-1. i will be rounded to its nearest whole
#'  value if it is not an integer.
#'@param z_i The interim statistic at analysis i.
#'
#'@return A list with the elements as follows: \itemize{ \item{K: As in d.}
#'  \item{n.I: As in d.} \item{u_K: As in d.} \item{lowerbounds: As in d.}
#'  \item{i: i used in computation.} \item{z_i: As input.} \item{cp: A matrix of
#'  conditional powers under different response rates.} \item{p_1: As input.}
#'  \item{p_0: As input.} }
#'
#'
#'@section Reference: \itemize{ \item{Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
#'  Applications to Clinical Trials. Chapman and Hall/CRC, Boca Raton, FL, 2000.} }
#'
#'@seealso \code{\link{exactprob}}, \code{\link{asymcp}},
#'  \code{\link{exactdesign}}.
#'
#'@export
#'
#' @examples
#'I=c(0.2,0.4,0.6,0.8,0.99)
#'beta=0.2
#'betaspend=c(0.1,0.2,0.3,0.3,0.2)
#'alpha=0.05
#'p_0=0.3
#'p_1=0.5
#'K=4.6
#'tol=1e-6
#'tt1=asymdesign(I,beta,betaspend,alpha,p_0,p_1,K,tol)
#'tt2=exactdesign(tt1)
#'tt3=exactprob(p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),d=tt2)
#'exactcp(tt2,p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),1,2)
#'exactcp(tt3,p_1=c(0.4,0.5,0.6,0.7,0.8,0.9),3,19)
exactcp<-function(d,p_1,i,z_i){
 ##check validity of inputs
 if(!(methods::is(d,"exactdesign")|methods::is(d,"exactprob")))
    stop("exactcp must be called with class of either exactdesign or exactprob.")
 p_0=d$p_0
 i=round(i)
 z_i=round(z_i)
 if((min(p_1)<=p_0)|(max(p_1)>=1))
    stop('Please input p_1 that lies between p_0 and 1 (not including p_0 and 1).')
 if((i>=K)|(i<1))
    stop('i is less than 1 or more than K-1.')
 if((z_i>=1000)|(z_i<=0))
    stop('invalid interim statistic: too large or too small.')
 K=d$K-i  ##the no. of analysis after the ith analysis.
 n.I=(d$n.I)[i:(i+K)]  ##sample sizes for i,i+1,...,K analysis
 lowerbounds=(d$lowerbounds)[(i+1):(i+K)]  ##lower bounds for i+1,i+2,...K analysis.
 lowerbounds=lowerbounds-z_i  ##new lower bounds for Z_{i,j}
 n.I=n.I[-1]-n.I[1]
 cp=matrix(0,1+length(p_1),1) ##to store conditional power under different p

 ##if z_i is euqal to or greater than u_K, cp is 1.
 if(lowerbounds[K]<=0){
     cp[,1]=1
     return(list(K=d$K,n.I=d$n.I,u_K=d$u_K,lowerbounds=d$lowerbounds,i=i,z_i=z_i,cp=cp,p_1=p_1,p_0=p_0))
 }

 ##if z_i is greater than one or some of the lower bounds.
 index1=which(lowerbounds<(-1))
 if(length(index1)>0){
     lowerbounds[1:max(index1)]=-1 #make these lower bounds equalling -1 so that the trial will continue even with 0 newly occurred events.
 }

 #conditional power under H0
 cp[1,]=exactprob1(n.I,lowerbounds,p_0,K,K)$phi ##cp under H0

 ##conditional power under H1.
 for(s in 1:length(p_1)){
    cp[(s+1),]=exactprob1(n.I,lowerbounds,p_1[s],K,K)$phi ##cp under H0
 }

 cp=cbind(c(p_0,p_1),cp)
 colnames(cp)=c('p','cp')
 return(list(K=d$K,n.I=d$n.I,u_K=d$u_K,lowerbounds=d$lowerbounds,i=i,z_i=z_i,cp=cp,p_1=p_1,p_0=p_0))
}

