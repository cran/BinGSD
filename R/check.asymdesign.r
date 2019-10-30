##check validity of parameters
##I is an increasing vector of length K. and the last element should be 1 or I will be standardized. or of length K-1 and the 
##last element be less than 1.
##beta should be a value in the interval (0,0.5]. alpha be a scalar belongs to (0,0.3]. 
##betaspend be a vector with length K, with all elements greater than 0 and less than 1. and sum to 1, or the true betaspend 
##would be betaspend/sum(betaspend). p_1,p_0 be two numbers and p_1 greater than p_0.

check.asymdesign<-function(I,beta,betaspend,alpha,p_0,p_1,K,tol){
  ##--round up K and check if K is no less than 1 and no more than 20----
  K=round(K)
  if((K<=1)|(K>20))
    stop('Please input a K that lies between 1 and 20 (not including 1).')
  
  ##--check the sequence of information fractions I--
  if(min(I)<=0)
    stop('Element(s) of I should be positive.')
  temp1=length(I)
  temp2=I-c(0, I[1:(temp1-1)])  #temp2 is c(I[1],I[2]-I[1],I[3]-I[2],...)
  if(min(temp2)<=0)
    stop('I should be an increasing vector.')
  if(temp1==K){##if length of I is equal to K,check its values.
	if(I[K]!=1){
	   warning('I will be standardized so that the last element is 1.')
	   I=I/I[K]
	}
  }else{
    if(temp1!=(K-1))
	   stop('length of I is neither K nor K-1.')
    if(I[(K-1)]>=1)
	    stop('Please input I with I[K-1] less than 1.')
	I=c(I,1)  ##so that I if of length K.
  }
  
  ##--check beta----
  if(length(beta)>1)
    stop('beta is not a scalar.')
  if((beta<=0)|(beta>0.5))
    stop('Please input a beta that lies between 0 and 0.5 (not including 0).')

  ##--check betaspend----
  if(length(betaspend)!=K)
    stop('betaspend is not of length K.')
  if((min(betaspend)<0)|(max(betaspend)>1))
    stop('Element(s) in betaspend is less than 0 or greater than 1.')
  if(sum(betaspend)!=1){
    warning('betaspend will be standardized so that the total is 1.')
	betaspend=betaspend/sum(betaspend)
  }
	
  ##--check alpha--
  if(length(alpha)>1)
    stop('alpha is not a scalar.')
  if((alpha<=0)|(alpha>0.3))
    stop('Please input an alpha that lies between 0 and 0.3 (not including 0).')  
	
  ##--check p_0
  if(length(p_0)>1)
    stop('p_0 is not a scalar.')
  if((p_0<=0)|(p_0>=1))
    stop('Please input a p_0 that lies between 0 and 1 (not including 0 and 1).')
	
  ##--check p_1--
  if(length(p_1)>1)
    stop('p_1 is not a scalar.')
  if((p_1<=p_0)|(p_1>=1))
    stop('Please input a p_1 that lies between p_0 and 1 (not including p_0 and 1).') 
	
  ##check tol--
  if(length(tol)>1)
    stop('tol is not a scalar.')
  if(tol>0.01)
    stop('tol is greater than 0.01.')
  if(tol<=0)
    stop('tol is 0 or negative.')
	
  return(list(I=I,betaspend=betaspend,K=K))
}
