source("spenv.r", chdir = TRUE)
source("spenvlp.R", chdir = TRUE)
source("initial_value.r", chdir = TRUE)
source("loglik.R", chdir = TRUE)
source("initial_setup.r", chdir = TRUE)


n <- 2000
p <- 10
r <- 10
u <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- matrix(rnorm(n * r), n, r)

# res <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=1, weight=rep(1,r-u))
# Ginit <- initial_value(X,Y,u)

# lambda <- seq(0.1,0.01,length.out=100)
# model_vec <- rep(NA, length(lambda))
# for(l in 1:length(lambda)){
# 	res <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=lambda[l], weight=rep(1,r-u))
# 	loglik_tmp <- loglik(res$Gammahat, res$sigRes, res$invsigY, res$r, res$n)
# 	model_vec[l] <- -2 * loglik_tmp + log(res$n) * (sum(res$nonzero_index)-u) * u	
# }



res1 <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=10, 
				weight=rep(1,r-u),initial_value=initial_value_tmp)
loglik_tmp <- loglik(res1$Gammahat, res1$sigRes, res1$invsigY, res1$r, res1$n)
model_vec1 <- -2 * loglik_tmp + log(res1$n) * (sum(res1$nonzero_index)-u) * u	


res2 <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=1, 
				weight=rep(1,r-u),initial_value=initial_value_tmp)
loglik_tmp <- loglik(res2$Gammahat, res2$sigRes, res2$invsigY, res2$r, res2$n)
model_vec2 <- -2 * loglik_tmp + log(res2$n) * (sum(res2$nonzero_index)-u) * u	




