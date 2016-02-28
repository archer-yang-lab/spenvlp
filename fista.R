fista <- function(b2, b3, A1, A2, A3, lambda, 
   eps, maxit, weight, a_vec_init) {
 #################################################################################
   #data setup
   lmax_A1 <- max(eigen(A1)$values)
	lmax_A2 <- max(eigen(A2)$values)
	lmax_A3 <- max(eigen(A3)$values)
	gamma <- 4*lmax_A1 + 2*lmax_A2 + 2*lmax_A3
	#################################################################################
	nobs <- NROW(A1)
	dif <- rep(NA, nobs)
	a_vec = a_vec_init
	olda_vec = a_vec
	npass <- 0		
	while(1){
		sk = 1
		olda_vec <- a_vec
		tmp1=A1%*%olda_vec
		tmp2=A2%*%(olda_vec+b2)
		tmp3=A3%*%(olda_vec+b3)
		U <- drop(4*tmp1/drop(1+olda_vec%*%tmp1)-					
		2*tmp2/drop(1+(olda_vec+b2)%*%tmp2)-		
		2*tmp3/drop(1+(olda_vec+b3)%*%tmp3))
		U_working <- (U + gamma * olda_vec)
		U_norm <- drop(sqrt(crossprod(U_working,U_working)))
		t <- U_norm - weight * lambda	
		if(t > 0){
			a_vec <- U_working * t / (gamma * U_norm)
		}else{
			a_vec <- rep(0,nobs)
		}
		sk_old = sk
		sk = (1+sqrt(1+4*sk_old^2))/2
		a_vec = a_vec + (sk_old-1)*(a_vec - olda_vec)/sk
	   dif <- a_vec - olda_vec
	   if(gamma*sum(dif^2) < eps) break
	   npass = npass + 1
	   if(npass > maxit) break
	}
  ################################################################################
   # output
   outlist <- list(a_vec = a_vec)
   outlist
} 
