library(quantreg)
#similar to the temporal AFL
#neighbor.R: specify neighbors in the penalty term
source('neighbor.R')
adFAL.func <- function(y, tau=0.95, lambdas, plot.SIC=FALSE)
{
    y=as.vector(y)
    N = length(y)
    n = sqrt(N)
    # adaptive weights, Fedlity and Penalty matrices
    tmp.mat = design.mat.AF.func(y)
    Fid = tmp.mat$Fid
    R = tmp.mat$R
    #Use BIC to search for the appropriate lambda within lambdas
    ahat = SIC = NULL
    for(lam in lambdas)
    {
        Coef = as.vector(est.flasso(y,lam, tau=tau, Fid, R))
        uhat = y - Fid%*%Coef
        tmp.ahat = Coef[1:N]
        delta= round(tmp.ahat,4)
	    	edf=length(unique(delta))      
        tmp.SIC = log(sum(uhat*(tau-1*(uhat<0)))) + log(N)*edf
        SIC = c(SIC, tmp.SIC)
        ahat = rbind(ahat, tmp.ahat)
    }
    idx = order(SIC)[1]
    ahat = ahat[idx,]
    lam = lambdas[idx]
    cat("lam is", lam, "\n")
    Coef=ahat
    uhat = y - Fid%*%Coef
    if(plot.SIC){ plot(SIC~lambdas); abline(v=lam, col=2)}
    a.est = ahat
    return(list(lam=lam, a.est=a.est))
}



# obtain the arguments needed for rq.fit.sfn
design.mat.AF.func = function(y)
{
    N = length(y)
	  Z.csr=1
    X.csr = Z.csr %x% as(N,"matrix.diag.csr")
	  # penalization
    n=sqrt(N)
    R=as.matrix.csr(NearestNeighbors(n,-1))
    return(list(Fid=X.csr, R=R))
}


est.flasso = function(y, lam, tau=0.95, Fidelity, R)
{ 
  N = length(y)
	lambda = lam*rep(1,N)
  nlam = length(lambda)
	Penalty <- R* lambda
	D <- rbind(Fidelity, Penalty)
  y2 = c(y, rep(0, nlam))
	rhs = (1-tau) * (t(Fidelity)%*%rep(1,N)) +  0.5*(t(Penalty)%*%rep(1,nlam))
	# Sparse Regression Quantile Fitting
  coeff <- rq.fit.sfn(D, y2, rhs=rhs,tau=tau,sfn.control(nnzlmax=5*N^2))$coef
  return(coeff)
}

