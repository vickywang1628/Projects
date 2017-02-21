#include "RcppArmadillo.h"
// [[Rcpp::export]]

arma::vec QRMMCPP(arma::mat xr,arma::vec yr,arma::vec thetar,arma::vec Wr,double lambdar,double ta,double to,int m) 
{
// x should be n*n matrix
double toler=(to);
double tau=(ta);
double lambda=(lambdar);
int maxit=(m);
arma:: mat x=(xr), xt, xu, WAW;
arma:: vec y=(yr), r, W=(Wr), A, u, v; // delta
arma:: vec thetaold, theta=(thetar);
arma:: uvec order, index;

int n=x.n_rows;
int p=x.n_cols; // p = n!!
double error=10000, epsilon=0.9999;
//toler=1e-3;
int iteration=1;
r.zeros(n);

//beta= solve(x*U*x.t()+2*W*A*W, x*v);
WAW.zeros(n,n);
xu.zeros(n,p);
xt = x.t();

while (iteration<=maxit&& error>toler)
{
  thetaold = theta; 
  
  // below are MM steps
  r = y - x*theta;  
  u = 1/(arma::abs(r)+epsilon);
  v = -1 + 2*tau + y/(arma::abs(r)+epsilon);
  
  for (int i=0;i<n;i++)
  	{	A(i)=lambda/arma::abs(W(i)*theta(i));}

  for (int i=0;i<p;i++)
		{	xu.col(i)=x.col(i)*u(i);}
    
  for (int i=0;i<n;i++)
  	{	WAW[i, i] = W(i)*A(i)*W(i);}

 	theta = arma::solve(xu*xt + 2*WAW, x*v); 
  
	//theta = theta-delta;
  error = arma::sum(arma::abs(theta - thetaold));
  iteration++;
}

return theta;

}
