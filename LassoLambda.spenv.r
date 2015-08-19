LassoLambda.spenv <- function(X, Y, u, eps=1e-10, maxit=1e4, lambda, weight = rep(1,r-u)) {
	initial_value_tmp = initial_value(X,Y,u)
	model_vec <- rep(NA, length(lambda))
	for(l in 1:length(lambda)){
		res <- spenv(X, Y, u, 
					eps=eps, 
					maxit=maxit, 
					ulam=lambda[l], 
					weight=weight,
					initial_value=initial_value_tmp)
		loglik_tmp <- loglik(res$Gammahat, res$sigRes, res$invsigY, res$r, res$n)
		model_vec[l] <- -2 * loglik_tmp + log(res$n) * (sum(res$nonzero_index)-u) * u	
		initial_value_tmp <- res$Gammahat
	}
	model.min <- min(model_vec)
	idmin <- model_vec <= model.min
	lambda.min <- max(lambda[idmin])
	lambda.min.id <- which(lambda==lambda.min)
	res <- spenv(X, Y, u, 
					eps=eps, maxit=maxit, ulam=lambda.min, 
					weight=weight,initial_value=initial_value_tmp)
	loglik.min <- loglik(res$Gammahat, res$sigRes, res$invsigY, res$r, res$n)
	Gammahat.min <- res$Gammahat
	nonzero_index.min <- res$nonzero_index
	out <- list(lambda.min=lambda.min, 
					lambda.min.id=lambda.min.id, 
					model_vec=model_vec,
					loglik.min=loglik.min,
					Gammahat.min=Gammahat.min,
					nonzero_index.min=nonzero_index.min)
	out
}