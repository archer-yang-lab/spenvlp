bic.spenv <- function(X, Y, u, eps=1e-10, maxit=1e4, ulam, weight) {

for(l in 1:length(lambda)){
	res <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=lambda[l], weight=rep(1,r-u))
	loglik <- loglik(res$Gammahat, res$sigRes, res$invsigY, res$r, res$n)
	-2 * loglik + log()
	
}

}