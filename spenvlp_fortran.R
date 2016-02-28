spenvlp_fortran <- function(lambda, weight, b2, b3, A1, A2, A3, eps, maxit, a_vec_init) {
 #################################################################################
    #data setup
    lmax_A1 <- max(eigen(A1)$values)
	lmax_A2 <- max(eigen(A2)$values)
	lmax_A3 <- max(eigen(A3)$values)
	gamma <- 4*lmax_A1 + 2*lmax_A2 + 2*lmax_A3
	#################################################################################
	nobs <- NROW(A1)
	a_vec = a_vec_init	
	fit <- .Fortran("spenvlp_engine",
			npass = integer(1),
			as.integer(maxit),
			as.integer(nobs),
			as.double(lambda), 
			as.double(weight),
			as.double(eps), 
			as.double(gamma),
			a_vec = as.double(a_vec), 
			as.double(b2), 
			as.double(b3), 
			as.double(A1),
			as.double(A2),
			as.double(A3))
 ################################################################################
    # output
    outlist <- list(a_vec = fit$a_vec)
    outlist
} 
