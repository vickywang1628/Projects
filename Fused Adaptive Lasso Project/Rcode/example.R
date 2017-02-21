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
#temporal fused adaptive Lasso estimates of the 95% quantiles
m0=adFAL.func(y, tau=0.95,lambda, plot.SIC=T)
#Figure 1 in the paper
plot(y,pch=19,cex=0.5,xlab='t',ylab='Y')
lines(m0$a.est,col='red',lwd=2,lty=1)





rm(list=ls())
#spatial adaptive fused Lasso
source('square.R')

#simulate a spatial process on a 50 by 50 grid
library(SpatialExtremes)
set.seed(1023)
n=50
x <- seq(0, n, length = n)
coord <- cbind(x,x)
#two processes from Schlather's model
y1 <- rmaxstab(1, coord, cov.mod = "bessel", sill = 1, nugget=0, range = 0.01, smooth = 0.5, grid = T)
y2 <- 5+rmaxstab(1, coord, cov.mod = "bessel", sill = 1, nugget=0, range = 0.01, smooth = 0.5, grid = T)

grid=expand.grid(x,x)
d1=sqrt(grid[,1]^2+grid[,2]^2)
d2=sqrt((grid[,1]-n)^2+(grid[,2]-n)^2)
dm1=matrix(d1,n,n)+diag(n)*0.0001
dm2=matrix(d2,n,n)+diag(n)*0.0001
w1=1/dm1
w2=1/dm2
tot=w1+w2
#weighted average based on distance
y=(y1*w1+y2*w2)/tot
#plot simulated data
filled.contour(x, x, y, color.palette = terrain.colors,zlim=c(0,9))
title('Observations')

#range of lambda
tuning=seq(2,4,0.1)
#spatial fused adaptive Lasso estimates of the 95% quantiles
m0=adFAL.func(y, tau=0.95,tuning, plot.SIC=T)
m0mat=matrix(m0$a.est,n,n)
#plot quantile estimates
filled.contour(x, x, m0mat, color.palette = terrain.colors,zlim=c(0,9))
title('Quantile')

