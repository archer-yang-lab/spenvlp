spenv <- function(X, Y, u, eps=1e-5, maxit=3e2, lambda, weight, initial_value=NULL) {

	a <- dim(Y)
	n <- a[1]
	r <- a[2]
	p <- ncol(X)
	betahat = NA 
	etahat = NA 
	Gammahat = NA
	Gamma0hat = NA
	sigRes = NA
	invsigY = NA
	
	if (u == 0) {

		Gammahat <- NA
		Gamma0hat <- diag(r)
		betahat <- matrix(0, r, p)
		etahat <- betahat
		sigRes <- cov(Y) * (n - 1) / n
		invsigY <- chol2inv(chol(sigRes))
		
	} else if (u == r) {
	
		Yc <- as.matrix(scale(Y, center = T, scale = F))
		Xc <- as.matrix(scale(X, center = T, scale = F))
		invsigX <- chol2inv(chol(cov(X)))
		sigYX <- cov(Yc, Xc)
		betahat <- betaOLS <- sigYX %*% invsigX
		Gammahat <- diag(r)
		Gamma0hat <- NA
		sigY <- cov(Y) * (n-1) / n
		invsigY <- chol2inv(chol(sigY))
		sigRes <- sigY - sigYX %*% tcrossprod(invsigX, sigYX) * (n - 1) / n
		
	} else if (u == (r-1)) {
	
		if(missing(initial_value)) init <- initial_value(X,Y,u)
		else init=initial_value
		Ginit <- init %*% solve(init[1:u, ])

		# obj1 <- initial_value$obj1
		initial_setup_tmp <- initial_setup(X, Y, u)
		sigRes <- initial_setup_tmp$sigRes
		invsigY <- initial_setup_tmp$invsigY
		betaOLS <- initial_setup_tmp$betaOLS
			
		U1c2 <- array(0, dim = c(r-1, r-1))
		V1c2 <- array(0, dim = c(r-1, r-1))

		U1c2 <- sigRes[-r, -r] - as.matrix(sigRes[-r, r]) %*% sigRes[r, -r] / sigRes[r, r]
		V1c2 <- invsigY[-r, -r] - as.matrix(invsigY[-r, r]) %*% invsigY[r, -r] / invsigY[r, r]		

		i <- 1
		while (i < maxit) {
			t2 <- sigRes[-r, r] / sigRes[r, r]
			t3 <- invsigY[-r, r] / invsigY[r, r]
			invC1 <- chol2inv(chol(U1c2))
			invC2 <- chol2inv(chol(V1c2))
	
			res <- spenvlp(b2=drop(t2), b3=drop(t3), 
							A1=diag(r-1), A2=sigRes[r, r]*invC1, A3=invsigY[r, r]*invC2, 
							lambda=lambda, eps=eps, maxit=maxit, 
							weight=weight[r-u], 
							a_vec_init=drop(Ginit[r,]))
		
			old_Ginit <- Ginit[r, ]
			Ginit[r, ] <- res$a_vec
		
			if(sum((Ginit[r,]-old_Ginit)^2) < eps) break
			i <- i + 1		
		}
		a <- qr.Q(qr(Ginit), complete = TRUE)
		Gammahat <- a[, 1:u]
		Gamma0hat <- a[, r]
		etahat <- crossprod(Gammahat, betaOLS)
		betahat <- Gammahat %*% etahat	
	} 
	else {
		if(missing(initial_value)) init <- initial_value(X,Y,u)
		else init=initial_value
		Ginit <- init %*% solve(init[1:u, ])
    
		# obj1 <- initial_value$obj1
		initial_setup_tmp <- initial_setup(X, Y, u)
		sigRes <- initial_setup_tmp$sigRes
		invsigY <- initial_setup_tmp$invsigY
		betaOLS <- initial_setup_tmp$betaOLS

		GUG <- crossprod(Ginit, (sigRes %*% Ginit))	
		GVG <- crossprod(Ginit, (invsigY %*% Ginit))		

		t4 <- crossprod(Ginit[(u+1):r,], Ginit[(u+1):r, ]) + diag(u)
		i <- 1
		while (i < maxit) {
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

				res <- spenvlp(b2=drop(t2), b3=drop(t3), A1=invt4, A2=sigRes[j, j]*invC1, A3=invsigY[j, j]*invC2, lambda=lambda, eps=eps, maxit=maxit, weight=weight[j-u], a_vec_init=drop(Ginit[j,]))
				
				old_Ginit <- Ginit[j, ]
				Ginit[j, ] <- res$a_vec
				g <- as.matrix(Ginit[j, ])
				# print(g)
				t4 <- t4 + tcrossprod(g, g)
				GUGt2 <- g + t2
				GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
				
				GVGt2 <- g + t3
				GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
			}
			if(sum((Ginit[j,]-old_Ginit)^2) < eps) break
			i <- i + 1
		}
		a <- qr.Q(qr(Ginit), complete = TRUE)
		Gammahat <- a[, 1:u]
		Gamma0hat <- a[, (u+1):r]
		etahat <- crossprod(Gammahat, betaOLS)
		betahat <- Gammahat %*% etahat
	}
	cat('The number of iterations:',i,'.\n',sep='')
	if(!is.na(sum(Gammahat))) nonzero_index <- (rowSums(abs(Gammahat))>0)
	else nonzero_index <- NA
	return(list(betahat = betahat, etahat = etahat, Gammahat = Gammahat, Gamma0hat = Gamma0hat, sigRes = sigRes, invsigY = invsigY, r = r, n = n, nonzero_index = nonzero_index))
}
	


