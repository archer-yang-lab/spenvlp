loglik <- function(Gammahat, sigRes, invsigY, r, n) {
	e1 <- eigen(t(Gammahat) %*% sigRes %*% Gammahat)
	e2 <- eigen(t(Gammahat) %*% invsigY %*% Gammahat)				
	-n / 2 * (sum(log(e1$values)) + sum(log(e2$values))) - n * r / 2 *log(2 * pi) - n * r / 2
}
	


