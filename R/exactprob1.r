#compute the lower boundary crossing probabilities given the design, under H0 or H1. using binomial distribution.
##this function borrows the idea of gsBinomialExact of package gsDesign.
##exactprob1(n.I,c(lowerbounds,lowerbounds[(k-1)]+1),p_1,k,K)
##n.I is the sample size for the first k analysis. lowerbounds are for the first k analysis.
##this function can be checked by using the results of gsBinomialExact with b=rep(n.I[1:(k-1)]+1,l_k) and a=c(lowerbounds[1:(k-1)],l_k-1)
exactprob1<-function(n.I,lowerbounds,p,k,K){
 m=c(n.I[1],diff(n.I))  ##m is the increment of sample size at each analysis,m=c(n_1,n_2-n_1)
 plo=rep(0,k)  ###store the probabilities of crossing the lower boundaries defined in (12.6)
 phi=NULL   ##when k=K, phi would be the upper boundary crossing probability
 c.mat=matrix(0,ncol=k,nrow=n.I[k]+1) ### c.mat is the recursive function defined in (12.5)
 c.mat[,1]=stats::dbinom(0:n.I[k],m[1],p) ##the probability of 0:n.I[k] responses occured at the first analysis
 plo[1]=sum(c.mat[(0:n.I[k])<=lowerbounds[1],1])
 for(i in 2:k){
    no.stop=((lowerbounds[(i-1)]+1):n.I[(i-1)])  #the number of responses occured that falls into the countine interval before ith analysis
    no.stop.mat=matrix(no.stop,byrow=T,nrow=n.I[k]+1,ncol=length(no.stop)) #j in (12.5) with each col stands for a different j.
    succ.mat=matrix(0:n.I[k],byrow=F,ncol=length(no.stop),nrow=n.I[k]+1)#y in (12.5) with each row stands for a different y.
    bin.mat=matrix(stats::dbinom(succ.mat-no.stop.mat,m[i],p),byrow=F,ncol=length(no.stop),nrow=n.I[k]+1)#B_{m_i}(y-j,p) in (12.5)
    c.mat[,i]=bin.mat%*%c.mat[no.stop+1,(i-1)] ##c_i(y,p)=sum{c_{i-1}(j,p)*B_{m_i}(y-j,p)} for each y.
	plo[i]=sum(c.mat[(0:n.I[k])<=lowerbounds[i],i])  ##r_k^l in (12.6)
	if(i==K){
	  plo[i]=sum(c.mat[(0:n.I[k])<lowerbounds[i],i]) #for the last analysis, lower bound crossing prob is P(z<l_K)
	  phi=sum(c.mat[(0:n.I[k])>=lowerbounds[i],i])   #for the last analysis, upper bound crossing prob is P(z>=l_K)
	}
 }
 return(list(plo=plo,phi=phi))
}



