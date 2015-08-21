source("spenv.r", chdir = TRUE)
source("spenvlp.R", chdir = TRUE)
source("initial_value.r", chdir = TRUE)
source("loglik.R", chdir = TRUE)
source("initial_setup.r", chdir = TRUE)
source("grams.R", chdir = TRUE)
source("nulbasis.R", chdir = TRUE)
source("rref.R", chdir = TRUE)
source("LassoLambda.spenv.r", chdir = TRUE)
source("ADLassoLambda.spenv.r", chdir = TRUE)


library(mnormt)
set.seed(2015) 
p = 3
r = 10
q = 5
u = 2
n = 200

G_work=grams(matrix(rnorm(q*u),q,u));
G=rbind(G_work,matrix(0,r-q,u)) #r by u
G0 = grams(nulbasis(t(G)))

eta= 5*matrix(rnorm(u*p),u,p)
beta = G%*%eta; # r by p

sigma1 = 1
sigma2 = 5
sigmai = 2
Omega =  sigma1^2*diag(u);
Omega0 = diag( c(rep(sigma2^2, q - u),rep(sigmai^2, r-q)))
Sigma1 = G%*%Omega%*%t(G)
Sigma2 = G0%*%Omega0%*%t(G0)
Sigma = Sigma1 + Sigma2

B = matrix(rnorm(p*p),p,p)
SigmaX = t(B) %*% B

errMean = rep(0,r)
muX= rep(0, p)
alpha=matrix(rnorm(r),r,1)

# -----------------

X=rmnorm(n,muX,SigmaX)
Y=matrix(1,n,1)%*%t(alpha)+X%*%t(beta)+rmnorm(n,errMean,Sigma);

lambda_tmp <- seq(1,0.001,length.out=100)
lambda <- exp(seq(log(max(lambda_tmp)),log(min(lambda_tmp)),len=length(lambda_tmp)))

#adaptive lasso step
m1 <- ADLassoLambda.spenv(X, Y, u, eps=1e-10, maxit=1e4, lambda)

