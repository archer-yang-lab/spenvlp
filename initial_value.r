initial_value <- function(X, Y, u) {

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

	tmp.y <- eigen(sigY)
	bsxb <- tcrossprod(betaOLS, sigYX)
	tmp2.y <- crossprod(tmp.y$vectors, bsxb)
	tmp3.y <- sort(diag(tcrossprod(tmp2.y, tmp2.y)), decreasing = TRUE, index.return = TRUE)
	init <- as.matrix(tmp.y$vectors[, tmp3.y$ix[1:u]]) 
	
	if (n > r + 1) {
		eig1 <- eigen(t(init) %*% sigRes %*% init)
		eig2 <- eigen(t(init) %*% invsigY %*% init)
		obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))
	
		tmp2 <- diag(1/sqrt(tmp.y$values)) %*% tmp2.y
		tmp3 <- sort(diag(tcrossprod(tmp2, tmp2)), decreasing = TRUE, index.return = TRUE)
		init.y <- as.matrix(tmp.y$vectors[, tmp3$ix[1:u]])
		e1 <- eigen(t(init.y) %*% sigRes %*% init.y)
		e2 <- eigen(t(init.y) %*% invsigY %*% init.y)
		obj2 <- sum(log(e1$values)) + sum(log(e2$values))		
		if (obj2 < obj1) {
			init <- init.y
			obj1 <- obj2
		}
	
		if (n > r + p + 1) {
			tmp.res <- eigen(sigRes * (n-1) / n)
			tmp2.res <- crossprod(tmp.res$vectors, bsxb)
			tmp3.res <- sort(diag(tcrossprod(tmp2.res, tmp2.res)), decreasing = TRUE, index.return = TRUE)	
			init.res <- as.matrix(tmp.res$vectors[, tmp3.res$ix[1:u]])
			e1 <- eigen(t(init.res) %*% sigRes %*% init.res)
			e2 <- eigen(t(init.res) %*% invsigY %*% init.res)
			obj3 <- sum(log(e1$values)) + sum(log(e2$values))			
			if (obj3 < obj1) {
				init <- init.res
				obj1 <- obj3
			}
		
			tmp2.res <- diag(1/sqrt(tmp.res$values)) %*% tmp2.res
			tmp3.res <- sort(diag(tcrossprod(tmp2.res, tmp2.res)), decreasing = TRUE, index.return = TRUE)
			init.res <- as.matrix(tmp.res$vectors[, tmp3.res$ix[1:u]])
			e1 <- eigen(t(init.res) %*% sigRes %*% init.res)
			e2 <- eigen(t(init.res) %*% invsigY %*% init.res)				
			obj4 <- sum(log(e1$values)) + sum(log(e2$values))			
			if (obj4 < obj1) {
				init <- init.res
				obj1 <- obj4
			}
		}
	}
  init
}
	


