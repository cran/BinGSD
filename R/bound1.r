#compute the lower bound for the kth analysis given:
#(1)lowerbounds: lower bounds before kth analysis, l_1,l_2,...,l_{k-1}.
#(2)n.I: sample sizes before and at kth analysis, n_1,n_2,...,n_k.
#(3)betak: desired type II error spent at kth analysis.
#(4)k: the index of the analysis at which the lower bound is computed.
#(5)p_1,p_0,the last upper bound u_K.
#NOTE:the computation of lower bounds are taken under H1!!
#bound1(k,lowerbounds[1:(k-1)],u_K,n.I[1:k],p_1,p_0,betaspend[k],tol)

bound1<-function(k,lowerbounds,u_K,n.I,p_1,p_0,betak,tol){
  mean1=(p_1-p_0)*sqrt(n.I/p_1/(1-p_1)) ##the mean vector of the multivariate normal distribution for analysis 1,...,k under H1.
  lowerlimits=c(lowerbounds,-Inf) ##the lower limits for computing the boundary crossing probability.
  sigma=matrix(0,k,k) #the covariance matrix of multivariate normal distribution.
  for(i in 1:k){
    for(j in 1:k){
        sigma[i,j]=mean1[min(i,j)]/mean1[max(i,j)]
	}	
  }
  
  uppermin=c(rep(Inf,(k-1)),lowerbounds[(k-1)])  ##use the lower bound at (k-1)th analysis as an initial lower searching boundary.  
  uppermax=c(rep(Inf,(k-1)),u_K) ##use the last upper bound as the initial upper searching boundary.
  errormin=mvtnorm::pmvnorm(lower=lowerlimits,upper=uppermin,mean=mean1,sigma=sigma)[1]##the type II error given l_k=l_{k-1}.
  errormax=mvtnorm::pmvnorm(lower=lowerlimits,upper=uppermax,mean=mean1,sigma=sigma)[1]##the type II error given l_k=u_K.
  if(((errormin+tol)>=betak)&(betak>=errormin)){ ##if l_k=l_{k-1} achieves the type II error requirement.
    return(list(l_k=lowerbounds[(k-1)],error=errormin,flag=0)) ##flag==0 means everything is fine.
  }
  if(((errormax+tol)>=betak)&(betak>=errormax)){ ##if l_k=u_K achieves the type II error requirement.
    return(list(l_k=u_K,error=errormax,flag=1)) ##flag==1 means lower bounds after kth analysis should be set to be u_K.
  }
 
  ##update the lower and upper searching bounds, uses bisection twice to narrow the searching interval
  boundmax=u_K
  boundmin=lowerbounds[(k-1)]
  for(i in 1:2){
    boundnew=(boundmin+boundmax)/2 
    uppernew=c(rep(Inf,(k-1)),boundnew)   
    errornew=mvtnorm::pmvnorm(lower=lowerlimits,upper=uppernew,mean=mean1,sigma=sigma)[1]##the type II error given boundnew.
    if(((errornew+tol)>=betak)&(betak>=errornew)){ ##if l_k=boundnew achieves the type II error requirement.
      return(list(l_k=boundnew,error=errornew,flag=0))
    }
    if(errornew>betak){##in this case,boundnew should be the new upper searching bound.
      boundmax=boundnew
	  errormax=errornew
    }else{
      boundmin=boundnew
	  errormin=errornew
    }
  }
	
  ##boundmin and boundmax did not meet the error requirement, so try with weighted average of boundmin and boundmax. 
  #the weights are chosen so that if (errormax-betak)>(betak-errormin),then the new bound would be closer to boundmin.
  boundnew=(boundmin*(errormax-betak)+boundmax*(betak-errormin))/(errormax-errormin)  
  uppernew=c(rep(Inf,(k-1)),boundnew)   
  errornew=mvtnorm::pmvnorm(lower=lowerlimits,upper=uppernew,mean=mean1,sigma=sigma)[1]##the type II error given boundnew.

  ##if none of the initial value of boundnew, u_K and l_{k-1} satisfies the type II error constraint.  
  t=0
  while((((errornew+tol)<betak)|(betak<errornew))&(t<=30)){
    if(betak>errornew){ ##in this case, errornew is too small, l_k should lie between boundnew and boundmax
       boundmin=boundnew
       errormin=errornew	   
	}else{##in this case,errornew is too large,l_k should lie between boundmin and boundnew
       boundmax=boundnew
	   errormax=errornew
	}
	boundnew=(boundmin*(errormax-betak)+boundmax*(betak-errormin))/(errormax-errormin) 
	uppernew=c(rep(Inf,(k-1)),boundnew)   
    errornew=mvtnorm::pmvnorm(lower=lowerlimits,upper=uppernew,mean=mean1,sigma=sigma)[1]
	t=t+1   
  }
 
  if(((errornew+tol)<betak)|(betak<errornew)) ##the last while iteration is stopped due to t>20 but not convergence
    stop('cannot converge with the current tol.')  
 
  return(list(l_k=boundnew,error=errornew,flag=0))  
}