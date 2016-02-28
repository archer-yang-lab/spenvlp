fista_backtracking <- function(b2, b3, A1, A2, A3, lambda, 
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
	L  <- 0.01
	eta <- 2 
	while(1){
		sk = 1
		olda_vec <- a_vec
		tmp1=A1%*%olda_vec
		tmp2=A2%*%(olda_vec+b2)
		tmp3=A3%*%(olda_vec+b3)
		U <- drop(4*tmp1/drop(1+olda_vec%*%tmp1)-					
		2*tmp2/drop(1+(olda_vec+b2)%*%tmp2)-		
		2*tmp3/drop(1+(olda_vec+b3)%*%tmp3))
		Ltmp = L
		ik = 1
		while(1){
			U_working <- (U + Ltmp * olda_vec)
			U_norm <- drop(sqrt(crossprod(U_working,U_working)))
			t <- U_norm - weight * lambda	
			if(t > 0){
				a_vec <- U_working * t / (Ltmp * U_norm)
			}else{
				a_vec <- rep(0,nobs)
			}
			dif <- F_fun(a_vec, b2, b3, A1, A2, A3, lambda) - Q_fun(a_vec, olda_vec, b2, b3, A1, A2, A3, lambda, U, Ltmp)
			if(dif > 0){
				Ltmp = Ltmp * eta^ik 
				ik <- ik + 1
			}else{
				L = Ltmp
				break
			}
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
   outlist <- list(a_vec = a_vec, L = L)
   outlist
} 

F_fun <- function(a, b2, b3, A1, A2, A3, lambda){
	-2*log(1+a%*%A1%*%a)+log(1+(a+b2)%*%A2%*%(a+b2))+log(1+(a+b3)%*%A3%*%(a+b3))
	+ lambda*sqrt(crossprod(a,a))
}

Q_fun <- function(a, b, b2, b3, A1, A2, A3, lambda, U, L){
 F_fun(b, b2, b3, A1, A2, A3, lambda) - (a-b)%*%U + L*crossprod(a-b,a-b)/2 + lambda*sqrt(crossprod(a,a)) 
}