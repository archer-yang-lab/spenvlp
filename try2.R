source("spenv.r", chdir = TRUE)
source("spenvlp.R", chdir = TRUE)


n <- 2000
p <- 10
r <- 500
u <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- matrix(rnorm(n * r), n, r)


a = init()
anorm = norm(a)
for(i in nlam){
	a[,i] <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=0.2, weight=anorm)
}
best_a <- select_optimal(a)

best_a_norm = norm(best_a)
for(i in nlam){
a <- spenv(X, Y, u, eps=1e-10, maxit=1e4, ulam=0.2, weight=best_a_norm)
}