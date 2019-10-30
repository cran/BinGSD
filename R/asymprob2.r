#compute the lower boundary crossing probabilities given the design, under H0.
#asymprob2(n.I,lowerbounds,K)

asymprob2<-function(n.I,lowerbounds,K){
  sigma=matrix(0,K,K) #the covariance matrix of multivariate normal distribution.
  for(i in 1:K){
    for(j in 1:K){
        sigma[i,j]=sqrt(n.I[min(i,j)]/n.I[max(i,j)])
	}	
  }
  problow=rep(0,K)
  problow[1]=stats::pnorm(lowerbounds[1]) ##Z_1 follows a standard normal distribution.
  ##note the last (K-k) lower and upper integration would not influence the result.
  for(k in 2:K){
    upperlimits=c(rep(Inf,k-1),lowerbounds[k],rep(Inf,(K-k)))
	lowerlimits=c(lowerbounds[1:(k-1)],rep(-Inf,(K-k+1)))
	problow[k]=mvtnorm::pmvnorm(lower=lowerlimits,upper=upperlimits,mean=rep(0,K),sigma=sigma)[1]
  }
  return(problow)
}