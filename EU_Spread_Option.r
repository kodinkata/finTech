rm(list=ls())
# simulate a multiple dimension Geometric Brownian motion process
#Pricing an European spread Option on two stocks with the following characteristics:
#S1=S2 = 100, K = 1, r = 0.06, T = 1, sigma1 = 0.2, sigma2 = 0.3, and rho = 0.5.
# 1 - Properties for the 1st stock
# 2 - Properties for the 2nd stock
# S - current stock price
# K - Option Strike Price
# r - risk free rate
# T - time to expiration (1 = 1 year)
# sigma - volatility of the stock
# rho - correlation between the 2 stocks




# Define a function to simulate paths:
myGBM <- function(M,sigma,rho, u,x0,N=1,t=1){
  # Generate 3-dim X array for results:
  d <- length(u)
  B <- t(chol(Sigma))
  
  X <- array(NA,dim=c(M,N+1,d))
  # Set initial values = X0:
  X[,1,] <- x0
  dt <- t/N
  muMat <- matrix(rep(u,M),M,d,byrow=T)
  # Run the simulation:
  for (i in 1:N){
    Zmat <- matrix(rnorm(d*M),d,M)
    X[,i+1,] <- X[,i,]*exp(muMat*dt + sqrt(dt)*t(B%*%Zmat))
  }
  return(X)
}


x0 <- matrix(100,2,1)
r <- .06
u <- rep(r,2)
t <-1
M <- 10000
K <- 1
sigma <- c(0.2,0.3)
rho <- matrix(c(1,.5,.5,1),2,2,byrow = T)
Sigma <- diag(sigma) %*% rho %*% diag(sigma)
ans <- myGBM(M,Sigma,rho,u,x0)

soln <- pmax(ans[,2,1]-ans[,2,2]-K,0)
ans <- exp(-t*r)* mean(soln)
message(sprintf("The price of the option is %.4f",ans))