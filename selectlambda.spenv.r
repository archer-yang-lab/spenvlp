selectlambda.spenv <- function(X, Y, u, eps=1e-10, maxit=1e4, lambda) {

	initial_value_tmp = initial_value(X,Y,u)
	lambda <- seq(1,0.01,length.out=100)
	model_vec <- rep(NA, length(lambda))
	for(l in 1:length(lambda)){
		res <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=lambda[l], 
						weight=rep(1,r-u),initial_value=initial_value_tmp)
		loglik_tmp <- loglik(res$Gammahat, res$sigRes, res$invsigY, res$r, res$n)
		model_vec[l] <- -2 * loglik_tmp + log(res$n) * (sum(res$nonzero_index)-u) * u	
		initial_value_tmp <- res$Gammahat
	}
	model.min <- min(model_vec)
	idmin <- model_vec <= model.min
	lambda.min <- max(lambda[idmin])
	
}