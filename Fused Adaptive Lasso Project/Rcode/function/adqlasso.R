library(quantreg)
options(warn=-1)

# Fused adaptive Lasso quantile estimation
# tau: quantile level of interest
#Tuning parameter is selected by BIC
#lambda is a vector of tuning parameter candidates
# plot.SIC=TRUE: then plot SIC versus lambdas to help decide on the suitable range for choosing lambda
adFAL.func <- function(y, tau=0.95, lambdas, plot.SIC=FALSE)
{
    y=as.vector(y)
    N = length(y)
    # adaptive weights, Fedlity and R matrices
    tmp.mat = design.mat.AF.func(y,tau)
    Fid = tmp.mat$Fid
    R = tmp.mat$R
    wt = tmp.mat$wt
    # Use BIC to search for the appropriate lambda within lambdas
    ahat = SIC = NULL
    for(lam in lambdas)
    {
        Coef = as.vector(est.flasso(y,lam, tau=tau, Fid, R, wt))
        uhat = y - Fid%*%Coef
        tmp.ahat = Coef[1:N]
        delta = diff(tmp.ahat)
        edf = sum(abs(delta)>0.0001) 
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
# construct the design matrix, including Fidelity part and the penalty part
# tau: quantile level of interest
design.mat.AF.func = function(y,tau)
{
    N = length(y)
    Z.csr=1
    X.csr = Z.csr %x% as(N,"matrix.diag.csr")
	  # penalization
    R <- (as(1,"matrix.diag.csr") %x% as.matrix.csr(diff(diag(N))))
	  # adaptive weights
    b0.sfn = rq.fit.sfn(X.csr, y, tau=tau)$coef
    
    
    
    
    wt = 1/ pmax(abs(diff(b0.sfn)), 0.00001)
    wt =  wt[-N]
    return(list(Fid=X.csr, R=R, wt=wt))
}

# obtain the fused adaptive Lasso estimators for a given lambda value
# The estimation is obtained by using sfn: sparse fn
# y: observed response vector
# Fidelity: output from design.mat.csr
# R: output from design.mat.csr
# lam: tuning parameter
# tau: quantile level of interest, default is 0.95
# wt: the adaptive weights
est.flasso = function(y, lam, tau=0.95, Fidelity, R, wt )
{
    N = length(y)
    ### penalization
    lambda = lam*wt
    nlam = length(lambda)
    Penalty <- R*lambda
    D <- rbind(Fidelity, Penalty) #structure of the design matrix X stored in csr format
    y2 = c(y, rep(0, nlam)) # response vector
    #the right-hand-side of the dual problem
    rhs = (1-tau) * (t(Fidelity)%*%rep(1,N)) +  0.5*(t(Penalty)%*%rep(1,nlam))
    # Sparse Regression Quantile Fitting
    coeff <- rq.fit.sfn(D, y2, rhs=rhs, tau=tau)$coef
    return(coeff) 
}
