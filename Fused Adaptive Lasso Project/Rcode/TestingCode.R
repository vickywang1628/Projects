set.seed(922)
n=1000
z=rnorm(n)
tt=1:n
q1=tt<=n/3
q2=tt<=2*n/3
q3=tt<=n
sig=0.5*(q1+q2+q3)
y=z*sig

tau=0.95;
lambda = 18;

n=length(y)
x=matrix(rep(1,n*n),n,n)
x[lower.tri(x)] <- 0
theta=lm(y~x-1)$coef
#theta=rep(1,n)
W = rep(1, n)
maxit=200
toler=1e-3

xt = t(x);
epsilon=0.01;
##########################################
iteration=1;
error=1;
f_value = c();
alpha = 1;
while (iteration <= maxit && error[length(error)]>toler)
{
  thetaold = theta;    
  # below are MM steps
  r = y - xt %*% theta;  
  u = 1/(abs(r)+epsilon);
  #v = -1 + 2*tau + as.vector(y)/as.vector(abs(r)+epsilon);
  v = -1 + 2*tau + as.vector(r)/as.vector(abs(r)+epsilon);
  
  A = lambda/abs(W*theta);
  XUX = x %*% (as.vector(u)*as.matrix(xt));
  WAW = diag(as.vector(W*A*W));
  
  delta = solve(XUX + 2*WAW, -as.matrix(x) %*% as.vector(v) + 2*WAW %*% theta);   
  theta = theta - alpha*delta;
  alpha = alpha*0.5;
  
  error = c(error,sum(abs(theta - thetaold)));
  iteration = iteration + 1;
  r = y - xt %*% theta;
  r1 <- r;
  r1[r1>0] <- 0;
  f = sum(tau*r - r1) + lambda*sum(abs(t(W) %*% theta));
  f_value <- c(f_value, f)
}
##########################################
beta1 <- abs(theta)^{-1}
beta1[beta1>sqrt(n)] <- sqrt(n)
W = beta1;
##########################################
qtau2 = xt %*% theta
plot(y,pch=19,cex=0.5,xlab='t',ylab='Y',ylim = c(-4,4))
lines(qtau2 ,col='red')



