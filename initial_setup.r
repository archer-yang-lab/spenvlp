initial_setup <- function(X, Y, u) {

	a <- dim(Y)
	n <- a[1]
	r <- a[2]
	p <- ncol(X)
	sigY <- cov(Y) 
	invsigY <- chol2inv(chol(sigY)) * n / (n-1)
	invsigX <- chol2inv(chol(cov(X)))
	sigYX <- cov(Y, X)
	sigRes <- (sigY - sigYX %*% tcrossprod(invsigX, sigYX)) * (n-1) / n
	betaOLS <- sigYX %*% invsigX

	return(list(sigRes = sigRes, invsigY = invsigY, betaOLS = betaOLS))
}
	


