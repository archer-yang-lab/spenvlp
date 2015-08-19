spenvlp <- function(b2, b3, A1, A2, A3, ulam, 
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
		olda_vec <- a_vec
		U <- drop(4*A1%*%olda_vec/drop(1+olda_vec%*%A1%*%olda_vec)-					
		2*A2%*%(olda_vec+b2)/drop(1+(olda_vec+b2)%*%A2%*%(olda_vec+b2))-		
		2*A3%*%(olda_vec+b3)/drop(1+(olda_vec+b3)%*%A3%*%(olda_vec+b3)))
		U_working <- (U + gamma * olda_vec)
		U_norm <- drop(sqrt(crossprod(U_working,U_working)))
		t <- U_norm - weight * ulam	
		if(t > 0){
			a_vec <- U_working * t / (gamma * U_norm)
		}else{
			a_vec <- rep(0,nobs)
		}
	    dif <- a_vec - olda_vec
	    if(gamma*sum(dif^2) < eps) break
	    npass = npass + 1
	    if(npass > maxit) break
	}
  ################################################################################
    # output
    outlist <- list(a_vec = a_vec)
    class(outlist) <- c("spenvlp")
    outlist
} 
