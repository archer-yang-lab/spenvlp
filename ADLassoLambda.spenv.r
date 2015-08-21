ADLassoLambda.spenv <- function(X, Y, u, lambda,eps=1e-10, maxit=1e4) {
  init = initial_value(X,Y,u)
	GEidx = GE(init)
	newY = Y[, GEidx]	
	newinit = init[GEidx,]
	m1 <- LassoLambda.spenv(X, newY, u, eps=eps, maxit=maxit, lambda, weight = rep(1,r-u),initial_value=newinit)
	#calculating weight
	w <- m1$Gammahat
	w_norm <- 1/sqrt(rowSums(w*w))[(u+1):r]
	#adaptive lasso step
  m2 <- LassoLambda.spenv(X, newY, u, eps=eps, maxit=maxit, lambda, weight = w_norm,initial_value=newinit)
  out=m2
  out$Gammahat=m2$Gammahat[order(GEidx),]
  out$betahat=m2$betahat[order(GEidx),]
  out$weight=w_norm
  out
}