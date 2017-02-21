#return neighboring structure for the penalty term in the spatial FAL method
NearestNeighbors<-function(n,r.par){

  n.grid<-n^2
  r=rep(0,n.grid^2)
  # construct gamma matrix
  index=NULL
  for(i in 1:(n.grid-n)){
	index=c(index,(i+n-1)*n.grid+i,(i-1)*n.grid+(i+n))
  }

  for(i in 1:(n^2-1)){
	# check boundary	
	if(i%%n!=0)	{
	 index=c(index,i*n.grid+i,(i-1)*n.grid+(i+1))
	}
  }
  r[index]=r.par
  r.matrix=as.matrix.csr(r,n.grid,n.grid)
  factor=diff(slot(r.matrix,"ia"))
 slot(r.matrix,"ra")=slot(r.matrix,"ra")/rep(factor,factor)  
  r.matrix+diag(n.grid)
  }