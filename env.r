env <- function(X, Y, u) {

	a <- dim(Y)
	n <- a[1]
	r <- a[2]
	p <- ncol(X)

	if (u == 0) {
	
		Gammahat <- c()
		Gamma0hat <- diag(r)
		betahat <- matrix(0, r, p)
		etahat <- betahat
		return(list(betahat = betahat, etahat = etahat, Gammahat = Gammahat, Gamma0hat = Gamma0hat))
		
	} else if (u == r) {
	
		Yc <- as.matrix(scale(Y, center = T, scale = F))
		Xc <- as.matrix(scale(X, center = T, scale = F))
		invsigX <- chol2inv(chol(cov(X)))
		sigYX <- cov(Yc, Xc)
		betaOLS <- sigYX %*% invsigX
		Gammahat <- diag(r)
		Gamma0hat <- c()
		return(list(betahat = betaOLS, etahat = betaOLS, Gammahat = Gammahat, Gamma0hat = Gamma0hat))
		
	} else if (u == (r-1)) {
	
		maxiter = 100
		ftol = 1e-5
	
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

		Ginit <- init %*% solve(init[1:u, ])

		U1c2 <- array(0, dim = c(r-1, r-1))
		V1c2 <- array(0, dim = c(r-1, r-1))

		U1c2 <- sigRes[-r, -r] - as.matrix(sigRes[-r, r]) %*% sigRes[r, -r] / sigRes[r, r]
		V1c2 <- invsigY[-r, -r] - as.matrix(invsigY[-r, r]) %*% invsigY[r, -r] / invsigY[r, r]		


		i <- 1
		while (i < maxiter) {
			t2 <- sigRes[-r, r] / sigRes[r, r]
			t3 <- invsigY[-r, r] / invsigY[r, r]
			invC1 <- chol2inv(chol(U1c2))
			invC2 <- chol2inv(chol(V1c2))
				
			fobj <- function(x) {
				tmp2 <- x + t2
				tmp3 <- x + t3
				T2 <- invC1 %*% tmp2	
				T3 <- invC2 %*% tmp3
				-2 * log(1 + sum(x^2)) + log(1 + sigRes[r, r] * crossprod(tmp2, T2)) + log(1 + invsigY[r, r] * crossprod(tmp3, T3))
			}
	
			gobj <- function(x) {
				tmp2 <- x + t2
				tmp3 <- x + t3
				T2 <- invC1 %*% tmp2	
				T3 <- invC2 %*% tmp3
				-4 * x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / sigRes[r, r] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigY[r, r] + crossprod(tmp3, T3))	
			}
	
			res <- optim(Ginit[r, ], fobj, gobj, method = "BFGS")

			if (abs(fobj(Ginit[r,]) - fobj(res$par)) < ftol * fobj(Ginit[r,])) {
				Ginit[r,] <- res$par
				break
			} else {
				Ginit[r,] <- res$par
				i <- i + 1
			}
			
		}
		a <- qr.Q(qr(Ginit), complete = TRUE)
		Gammahat <- a[, 1:u]
		Gamma0hat <- a[, r]
		etahat <- crossprod(Gammahat, betaOLS)
		betahat <- Gammahat %*% etahat
		return(list(betahat = betahat, etahat = etahat, Gammahat = Gammahat, Gamma0hat = Gamma0hat))	
	
	} else {

		maxiter = 100
		ftol = 1e-5
		

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
		

		Ginit <- init %*% solve(init[1:u, ])


		GUG <- crossprod(Ginit, (sigRes %*% Ginit))	
		GVG <- crossprod(Ginit, (invsigY %*% Ginit))		

		
		t4 <- crossprod(Ginit[(u+1):r,], Ginit[(u+1):r, ]) + diag(u)
		i <- 1
		while (i < maxiter) {
			for (j in (u+1):r) {
				g <- as.matrix(Ginit[j, ])
				t2 <- crossprod(Ginit[-j, ], as.matrix(sigRes[-j, j])) / sigRes[j, j]
				t3 <- crossprod(Ginit[-j, ], as.matrix(invsigY[-j, j])) / invsigY[j, j]
				
				GUGt2 <- g + t2
				GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
				
				GVGt2 <- g + t3
				GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
				
				t4 <- t4 - tcrossprod(as.matrix(Ginit[j, ]), as.matrix(Ginit[j, ]))
				invC1 <- chol2inv(chol(GUG))
				invC2 <- chol2inv(chol(GVG))
				invt4 <- chol2inv(chol(t4))
						
				fobj <- function(x) {
					tmp2 <- x + t2
					tmp3 <- x + t3
					T1 <- invt4 %*% x
					T2 <- invC1 %*% tmp2	
					T3 <- invC2 %*% tmp3
					-2 * log(1 + x %*% T1) + log(1 + sigRes[j, j] * crossprod(tmp2, T2)) + log(1 + invsigY[j, j] * crossprod(tmp3, T3))
				}
			
				gobj <- function(x) {
					tmp2 <- x + t2
					tmp3 <- x + t3
					T1 <- invt4 %*% x
					T2 <- invC1 %*% tmp2	
					T3 <- invC2 %*% tmp3
					-4 	* T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / sigRes[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigY[j, j] + crossprod(tmp3, T3))	
				}
			
				res <- optim(Ginit[j,], fobj, gobj, method = "BFGS")
				Ginit[j, ] <- res$par
				g <- as.matrix(Ginit[j, ])
				t4 <- t4 + tcrossprod(g, g)
				GUGt2 <- g + t2
				GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
				
				GVGt2 <- g + t3
				GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
				
				
			}
	
			if (abs(fobj(Ginit[j,]) - res$value) < ftol * fobj(Ginit[j,])) {
				Ginit[j,] <- res$par
				break
			} else {
				Ginit[j,] <- res$par
				i <- i + 1
			}
		}
		a <- qr.Q(qr(Ginit), complete = TRUE)
		Gammahat <- a[, 1:u]
		Gamma0hat <- a[, (u+1):r]
		etahat <- crossprod(Gammahat, betaOLS)
		betahat <- Gammahat %*% etahat
		return(list(betahat = betahat, etahat = etahat, Gammahat = Gammahat, Gamma0hat = Gamma0hat))
	}
}
	


