##check validity of inputs

check.prob<-function(K,p_0,p_1,n.I,u_K,lowerbounds){
 ##--round up K and check if K is no less than 1 and no more than 20----
 K=round(K)
 if((K<=1)|(K>20))
    stop('Please input a K that lies between 1 and 20 (not including 1).')
	  
 ##--check validity of p_1 and p_0
 if(length(p_0)>1)
   stop('p_0 is not a scalar.')
 if((p_0<=0)|(p_0>=1))
   stop('Please input a p_0 that lies between 0 and 1 (not including 0 and 1).')
 if((min(p_1)<=p_0)|(max(p_1)>=1))
    stop('Please input p_1 that lies between p_0 and 1 (not including p_0 and 1).') 
	
 ##--check the sequence of sample sizes--
 n.I=round(n.I)
 if(min(n.I)<=0)
   stop('element(s) of n.I should be positive.')
 temp1=length(n.I)
 temp2=n.I-c(0, n.I[1:(temp1-1)])  #temp2 is c(n_1,n_2-n_1,...,n_K-n_{K-1})
 if(min(temp2)<=0)
   stop('n.I should be an increasing vector.')  
 if(temp1!=K)##if length of n.I is not equal to K.
   stop('length of n.I is not K.')
    
 ##check lowerbounds
 temp1=length(lowerbounds)
 temp2=lowerbounds[-1]-lowerbounds[1:(temp1-1)] #temp2 is c(l_2-l_1,l_3-l_2,...,l_K-l_{K-1}).
 if(min(temp2)<0)
   stop('lowerbounds should be a non-decreasing vector.')
 if(temp1!=K){##if length of lowerbounds is not equal to K,check its values.
   if(temp1!=(K-1))
      stop('length of lowerbounds is neither K nor K-1.')
   if(lowerbounds[temp1]>u_K)
      stop('lower bound for analysis k-1 is greater than u_K.')
   lowerbounds=c(lowerbounds,u_K)
 }
 if(u_K!=lowerbounds[K])
   stop('u_K must equal the last lower bound.')
 
 return(list(K=K,lowerbounds=lowerbounds,n.I=n.I))
}