QR.temp = function(y, tau, beta, lambda, maxit=200, toler=1e-3)
{
  n=length(y)
  L=matrix(rep(1,n*n),n,n)
  L[lower.tri(L)] <- 0
  
 
  if(missing(beta)){
	  #beta=lm(y~L-1)$coef
	  beta=rep(1,n)
  }

  W = rep(1, n)
  beta1 = QRMM(L,y,tau,beta,W,lambda,toler,maxit)$theta
  qtau1 = t(L)%*%beta1
  
  beta1 <- abs(beta1)^{-1}
  beta1[beta1>sqrt(n)] <- sqrt(n)
  W = beta1;

  rst = QRMM(L,y,tau, beta,W,lambda,toler,maxit)
  betah = rst$theta
  error = rst$error
  f = rst$f
  qtau2 = t(L)%*%betah
  #beta = matrix(betah[-1])
  #b=betah[1]
  return(list(qtau1=qtau1, qtau2=qtau2, error = error, f = f));
}
