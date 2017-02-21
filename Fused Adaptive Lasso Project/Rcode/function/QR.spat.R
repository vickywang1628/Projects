QR.temp = function(y, tau, beta, lambda, maxit=200, toler=1e-3)
{
  n=length(y)
  phi=
  
  L=matrix(rep(1,n*n),n,n)
  L[lower.tri(L)] <- 0
  
  if(missing(beta)){
    #beta=lm(y~L-1)$coef
    beta=rep(1,n)
  }
  
  W = rep(1, n)
  beta1 = abs(QRMM(L,y,tau,beta,W,lambda,toler,maxit))^{-1}
  
  beta1[beta1>sqrt(n)] <- sqrt(n)
  W = beta1;
  
  betah = QRMM(L,y,tau, beta,W,lambda,toler,maxit)
  #beta = matrix(betah[-1])
  #b=betah[1]
  return(beta=betah);
}
