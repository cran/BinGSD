#compute the lower boundary crossing probabilities given the design, under H1.
#asymprob1(n.I,lowerbounds,p_0,p_1,K)

asymprob1<-function(n.I,lowerbounds,p_0,p_1,K){
  mean1=(p_1-p_0)*sqrt(n.I/p_1/(1-p_1))  ##the mean vector of the multivariate normal distribution.
  sigma=matrix(0,K,K) #the covariance matrix of multivariate normal distribution.
  for(i in 1:K){
    for(j in 1:K){
        sigma[i,j]=mean1[min(i,j)]/mean1[max(i,j)]
	}	
  }
  problow=rep(0,K)
  problow[1]=stats::pnorm(lowerbounds[1]-mean1[1]) ##Z_1 follows a normal distribution.
  ##note the last (K-k) lower and upper integration would not influence the result since the intergrands are Inf or -Inf.
  for(k in 2:K){
    upperlimits=c(rep(Inf,(k-1)),lowerbounds[k],rep(Inf,(K-k)))
	lowerlimits=c(lowerbounds[1:(k-1)],rep(-Inf,(K-k+1)))
	problow[k]=mvtnorm::pmvnorm(lower=lowerlimits,upper=upperlimits,mean=mean1,sigma=sigma)[1]
  }
  return(problow)
}