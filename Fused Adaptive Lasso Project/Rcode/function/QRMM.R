QRMM = function(x, y, tau, theta, W, lambda, toler=1e-3, maxit=100) 
{
  iteration=1;
  # x should be n*n matrix
  # W should be n vector
  n= dim(x)[1];
  p= dim(x)[2];
  xt = t(x);
  f_value = c();

  error=1;
  epsilon=0.01;########?? how much is appropriate??
  
  #toler=1e-3;
  #beta= solve(x*U*x.t()+2*W*A*W, x*v);
  
  while (iteration <= maxit && error[length(error)]>toler)
  {
    thetaold = theta;    
    # below are MM steps
    r = y - xt %*% theta;  
    u = 1/(abs(r)+epsilon);
    v = -1 + 2*tau + y/(abs(r)+epsilon);
       
    A = lambda/abs(W*theta);
    XUX = x %*% (as.vector(u)*as.matrix(xt));
    WAW = matrix(0,n, n)
    diag(WAW) = (W*A*W);
    
    theta = solve(XUX + 2*WAW , x %*% v); 
    
#    theta = theta+delta;
    error = c(error,sum(abs(theta - thetaold)));
    iteration = iteration + 1;
    r = y - xt %*% theta;
    r1 <- r;
    r1[r1>0] <- 0;
    f = sum(tau*r - r1) + lambda*sum(abs(t(W) %*% theta));
    f_value <- c(f_value, f)
    print(iteration)
  }
  return (list(theta = theta, error = error, f = f_value));
  
}
