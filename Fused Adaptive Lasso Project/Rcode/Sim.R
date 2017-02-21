#https://github.com/cran/cqrReg


rm(list=ls())
#setwd("E:/Advanced Computing/790/Project/Rcode/Rcode")
setwd("/Volumes/TOSHIBA/Advanced Computing/790/Project/Rcode/Rcode")
#temporal adaptive fused Lasso
source('adqlasso.R')

#simulate time series with changing standard deviation 
set.seed(922)
n=1000
z=rnorm(n)
tt=1:n
q1=tt<=n/3
q2=tt<=2*n/3
q3=tt<=n
sig=0.5*(q1+q2+q3)
y=z*sig

#range of the tuning parameter lambda
lambda=seq(10,20,0.1)
tau=0.95
y=as.vector(y)
N = length(y)
# adaptive weights, Fedlity and R matrices
Z.csr=1
X.csr = Z.csr %x% as(N,"matrix.diag.csr")
Fid=X.csr



