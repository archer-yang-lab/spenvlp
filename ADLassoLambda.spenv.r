ADLassoLambda.spenv <- function(X, Y, u, eps=1e-10, maxit=1e4, lambda) {
	m1 <- LassoLambda.spenv(X, Y, u, eps=eps, maxit=maxit, lambda, weight = rep(1,r-u))
	#calculating weight
	w <- m1$Gammahat
	w_norm <- 1/sqrt(rowSums(w*w))[(u+1):r]
	#adaptive lasso step
	m2 <- LassoLambda.spenv(X, Y, u, eps=eps, maxit=maxit, lambda, weight = w_norm)
	m2
}