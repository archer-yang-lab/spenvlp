load("spam.rda")
source("aobjects.R", chdir = TRUE)
source("kernelmatrix.R", chdir = TRUE)
source("kernels.R", chdir = TRUE)
source("spenvlp.R", chdir = TRUE)
source("spenv.r", chdir = TRUE)

rbf <- rbfdot(sigma = 0.05) 
## create two vectors 
x <- rnorm(10) 
y <- rnorm(10)
## calculate dot product 

## use the spam data 
dt <- as.matrix(spam[c(10:20,3000:3010),-58])
## initialize kernel function 
rbf <- rbfdot(sigma = 0.05) 
rbf

## calculate kernel matrix 
A1 = kernelMatrix(rbf, dt[,1:10])
A2 = kernelMatrix(rbf, dt[,11:20])
A3 = kernelMatrix(rbf, dt[,21:30])
b2 <- rnorm(22)
b3 <- rnorm(22)

ulam <- 0.5
weight <- 1
eps <- 1e-18
maxit <- 1e5
a_vec_init=rnorm(22)

olda_vec <- spenvlp(b2=b2, b3=b3, A1=A1, A2=A2, A3=A3, ulam=ulam,eps=eps, maxit=maxit, weight=weight, a_vec_init=a_vec_init)

olda_vec <- olda_vec$a_vec

lmax_A1 <- max(eigen(A1)$values)
lmax_A2 <- max(eigen(A2)$values)
lmax_A3 <- max(eigen(A3)$values)
gamma <- 4*lmax_A1 + 2*lmax_A2 + 2*lmax_A3

U <- drop(4*A1%*%olda_vec/drop(1+olda_vec%*%A1%*%olda_vec)-					
2*A2%*%(olda_vec+b2)/drop(1+(olda_vec+b2)%*%A2%*%(olda_vec+b2))-		
2*A3%*%(olda_vec+b3)/drop(1+(olda_vec+b3)%*%A3%*%(olda_vec+b3)))
U_working <- (U + gamma * olda_vec)
U_norm <- drop(sqrt(crossprod(U_working,U_working)))
t <- U_norm - weight * ulam

anorm <- sqrt(crossprod(olda_vec,olda_vec))
if(anorm!=0){
	tmp  <- -U + ulam * weight * olda_vec/anorm
	print(tmp)	
}else{
	print(t)
}




