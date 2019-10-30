##calculate lower bound for kth analysis given sample sizes and l_1,...l_{k-1} using binomial distribution
#(1)lowerbounds: lower bounds before kth analysis, l_1,l_2,...,l_{k-1}.
#(2)n.I: sample sizes before and at kth analysis, n_1,n_2,...,n_k.
#(3)betak: desired type II error spent at kth analysis.
#(4)k: the index of the analysis at which the lower bound is computed.
#(5)p_1,p_0,the last upper bound u_K.
#NOTE:the computation of lower bounds are taken under H1!!
#bound2(k,lowerbounds[1:(k-1)],u_K,n.I[1:k],p_1,betaspend[k],K)

bound2<-function(k,lowerbounds,u_K,n.I,p_1,betak,K){
  errormin=(exactprob1(n.I,c(lowerbounds,lowerbounds[(k-1)]+1),p_1,k,K)$plo)[k]##the type II error given l_k=l_{k-1}+1.
  errormax=(exactprob1(n.I,c(lowerbounds,u_K),p_1,k,K)$plo)[k]##the type II error given l_k=u_K.
  if(errormin>betak)
    stop(paste0('cannot find a lower bound for ',k,'th analysis.'))
  if(errormax<=betak){ 
    if(k<(K-1))
	   stop(paste0('cannot find a lower bound for ',k,'th analysis.'))
    return(list(l_k=u_K,error=errormax)) 
  }
  ##until now, betamin<=betak,betamax>betak and l_{k-1}+1<u_K
  if((u_K-lowerbounds[(k-1)])==2){ ##then l_k=l_{k-1}+1=u_K-1 is the largest integer satisifes the error requirement.
    return(list(l_k=lowerbounds[(k-1)]+1,error=errormin))
  }
  
  ##until now, (l_{k-1}+1)<(u_K-1)
  ##update the lower and upper searching bounds until the lower searching bound+1=the upper. 
  #uses bisection only to narrow the searching interval since the solution is integer. it converges fast.
  boundmax=u_K
  boundmin=lowerbounds[(k-1)]+1
  while(boundmax>(boundmin+1)){ ##iterate until boundmax=boundmin+1. Since errormax>betak and errormin<=betak, the final errormin is the solution.
    boundnew=floor((boundmin+boundmax)/2)
    errornew=(exactprob1(n.I,c(lowerbounds,boundnew),p_1,k,K)$plo)[k]##the type II error given boundnew.
    if(errornew>betak){##in this case,boundnew should be the new upper searching bound.
      boundmax=boundnew
	  errormax=errornew
    }else{
      boundmin=boundnew
	  errormin=errornew
    }
  }
 
  ##since errormin<=betak,errormax>betak,boundmin=boundmax-1, boundmin is the largest integer that makes type II error<=betak 
  return(list(l_k=boundmin,error=errormin))  
}