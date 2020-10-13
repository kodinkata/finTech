rm(list=ls())
# install.packages("quantmod")
#Perform VaR analysis (95%) over 1 year horizon on a naive portfolio (equal allocation of capital into each stock) of 5 stocks 
#Initial investment: $1,000,000
#Using historical returns for the last 3 years
#Using Geometric Brownian Motion to simulate stock prices paths

library("quantmod")
library(dplyr)
# Define few companies:
tics <- c("BX","BABA","NVDA","UPS","VZ")
P.list <- lapply(tics, function(tic)
  get(getSymbols(tic, from = "1980-01-01")))
sapply(P.list,nrow)
n <- length(tics)

for (i in 1:n){
  print(tics[i])
  print(tail(P.list[[i]],4))
  print("-----------------------------------------")
}
# Get the adjusted prices into a single object
P.adj <- lapply(P.list, function(p) p[,6])
# Merge the elements of the list
P <- Reduce(merge,P.adj)
names(P) <- tics
head(P)
date(P[1,])
#Last 3 years of data
P <- subset(P,date(P)>="2017-10-27" ) 

# Returns:
Rets <- P/(stats::lag(P)) - 1
# Define parameters for the optimization:
# 1. The returns: retain the mean returns for all non-missing observations
R <- apply(Rets,2,function(x) mean(x,na.rm = TRUE))
# 2. The covariance matrix: run pairwise correlations
Sigma <- var(Rets, use="complete.obs")
n <- length(R)
# Annualize data:
R <- (1 + R)^252 - 1
Sigma <- Sigma * 252
head(R,10)

myGBMD <- function(M,N,d,t,mu,X0,Sigma){
  # Generate 3-dim X array for results:
  X <- array(NA,dim=c(M,N+1,d))
  # Set initial values = X0:
  X[,1,] <- X0
  B <- t(chol(Sigma))
  dt <- t/N
  muMat <- matrix(rep(mu-0.5*diag(sigma)*diag(sigma),M),ncol=d,byrow=T)
  # Run the simulation:
  for (i in 1:N){
    Zmat <- matrix(rnorm(d*M),ncol=M)
    X[,i+1,] <- X[,i,]*exp(muMat*dt + sqrt(dt)*t(B%*%Zmat))
  }
  return(X)
}
rho <- cor(Rets,use="complete.obs")
print(rho)
sigma <- var(Rets,use="complete.obs")
print(sigma)
Sigma.test <- diag(diag(sqrt(sigma)))%*%rho%*%diag(diag(sqrt(sigma)))*252
print(Sigma.test)

# Test to check:
print(Sigma - Sigma.test)
head(R,10)

print(Sigma)

M <- 10000
N <- 252
d <- 5
t <- 1

X <- myGBMD(M,N,d,t,mu=R,X0=100,Sigma)
# x <- X[,,1]

dX <- log(X[,c(2:(N+1)),]) - log(X[,c(1:(N)),])

dX.long <- matrix(NA,ncol=5,nrow=M*N)
for (i in 1:d){
  dX.long[,i] <- matrix(dX[,,i],ncol=1)
}
print(cor(dX.long))
print(rho)
print(cov(dX.long)*N)
print(Sigma)

var <- 1000000* exp(quantile(apply(dX.long,1,sum),0.05))
message(sprintf("Value At Risk (VaR) over 1 year is %.4f", var))
