load("spam.rda")
source("aobjects.R", chdir = TRUE)
source("kernelmatrix.R", chdir = TRUE)
source("kernels.R", chdir = TRUE)
source("spenvlp.R", chdir = TRUE)

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

ulam <- 1
weight <- 1
eps <- 1e-6
maxit <- 2000
a_vec_init=rnorm(22)

spenvlp(b2=b2, b3=b3, A1=A1, A2=A2, A3=A3, ulam=ulam,eps=eps, maxit=maxit, weight=weight, a_vec_init=a_vec_init)