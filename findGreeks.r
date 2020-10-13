#Find metrics (Greeks: delta and vega) for various options
#Assumption: stock price evolves under a geometric Brownian Motion process 
#Given characteristics: sigma = 30%, S0=$50, K = $52, and three months to maturity.  Assume a risk free rate of 1%.
#sigma - stock's volatility, S0- current stock price, K - strike price

# Option A: European call option 
# Option B: Asian option with payoff equal to max{Average(S0,ST) - K, 0}
# Option C: Digital option whose payoff is equal to $1 if ST > K and $0 otherwise
# Option D: An up and out barrier option with a barrier of $60.  The payoff is max[S(T) - K,0] if max{S(t)} < barrier and $0 otherwise.


#Step 1: Using	Finite Difference Method and compare to Black Scholes Model
#Step 2: Using Pathwise Estimator to find greeks
#Step 3: Using Likelihood Ratio Method to find greeks


#~~~~~~~~~~~~Step 1~~~~~~~~~~~~
rm(list=ls())
M <- 10000 # number of paths - will keep constant
# Input params:
n <- 100
h <- seq(0.001,0.1,length=n)
sigma <- 0.3
t <- 3/12
r <- 0.01
S0 <- 50
K <- 52
num <- 100

# 
# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- rep(pnorm(d1),n)
vega.true <- rep(S0*dnorm(d1)*t,n)

# Plain vanilla call
f <- function(S,K){
  f <- pmax(S-K,0)
}
biasD <- rep(NA,length=n)
biasV <- rep(NA,length=n)
varD <- rep(NA,length=n)
varV <- rep(NA,length=n)
ss <- sample(100000,num*num) #sample using a different seed each time

# Forward Difference Estimators:
for (i in 1:n){
  fd  <- c(NA)
  for (j in 1:num){
    set.seed(ss[j])
    ST <- S0*exp((r - 0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(M))
    set.seed(ss[j])
    STh <- (S0+h[i])*exp((r - 0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(M))
    f0 <- f(ST,K)
    f0h <- f(STh,K)
    fd[j]<-exp(-r*t)*mean((f0h - f0)/h[i])
    
  }
  varD[i]<- var(fd)
  biasD[i] <- mean(fd)-delta.true[1]
  }

plot(y=biasD,x=h,main="Bias of Delta Forward Difference Variance") 
plot(y=varD,x=h,main="Variance of Delta Forward Difference Variance") 


biasD <- rep(NA,length=n)
biasV <- rep(NA,length=n)
varD <- rep(NA,length=n)
varV <- rep(NA,length=n)

# Central Difference :
for (i in 1:n){
  cd <- c(NA)
  # ss <- sample(ss)
  for (j in 1:num){
    set.seed(ss[j])
    ST <- (S0-h[i])*exp((r - 0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(M))
    set.seed(ss[j])
    STh <- (S0+h[i])*exp((r - 0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(M))
    f0 <- f(ST,K)
    f0h <- f(STh,K)
    cd[j] <- exp(-r*t)*mean((f0h - f0) / (2*h[i]))
  }
  varD[i] <- var(cd)
  biasD[i] <- mean(cd)-delta.true[1]
}
plot(y=biasD,x=h,main="Bias of Delta Forward Difference Variance") 
plot(y=varD,x=h,main="Variance of Delta Forward Difference Variance") 

#Vega
vega.true <- rep(S0*dnorm(d1)*sqrt(t),n)

# Plain vanilla call payoff function
f <- function(S,K){
  f <- pmax(S-K,0)
}

biasV  <- rep(NA,length=n)
varV  <- rep(NA,length=n)
ss <- sample(1000000,num*num)
# Forward Difference Estimators :

for (i in 1:n){
  fd <- c(NA)
  # ss <- sample(ss)
  for (j in 1:num){
    set.seed(ss[j])
    ST <- S0*exp((r - 0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(M))
    set.seed(ss[j])
    STh <- (S0)*exp((r - 0.5*(sigma+h[i])^2)*t+(sigma+h[i])*sqrt(t)*rnorm(M))
    f0 <- f(ST,K)
    f0h <- f(STh,K)
    fd[j] <- exp(-r*t)*mean((f0h - f0) / h[i])
  }
  varV[i] <- var(fd)
  biasV[i] <- mean(fd)-vega.true[1]
}
plot(y=biasV,x=h,main="Bias of Vega (Forward Diff)")

plot(y=varV,x=h,main="Variance of Vega (Forward Diff)")
biasV  <- rep(NA,length=n)
varV  <- rep(NA,length=n)

# Central Difference Estimator:
for (i in 1:n){
  cd <- c(NA)
  # ss <- sample(ss)
  for (j in 1:num){
    set.seed(ss[j])
    ST <- (S0)*exp((r - 0.5*(sigma-h[i])^2)*t+(sigma-h[i])*sqrt(t)*rnorm(M))
    set.seed(ss[j])
    STh <- (S0)*exp((r - 0.5*(sigma+h[i])^2)*t+(sigma+h[i])*sqrt(t)*rnorm(M))
    f0 <- f(ST,K)
    f0h <- f(STh,K)
    cd[j] <- exp(-r*t)*mean((f0h - f0) / (2*h[i]))
  }
  varV[i] <- var(cd)
  biasV[i] <- mean(cd)-vega.true[1]
}


plot(y=biasV,x=h,main="Bias of Vega (Central Diff)")
plot(y=varV,x=h,main="Variance of Vega (Central Diff)")

#Step 2
#Delta of Option A

rm(list=ls())
M <- 10000 # number of paths - will keep constant
# Input params:
S0 <- 50
sigma <- 0.3
t <- 3/12
r <- 0.01
K <- 52
# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- pnorm(d1)
vega.true<-S0*dnorm(d1)*sqrt(t)
f <- function(S,K){
  f <- rep(0,length=length(S))
  f[S>K] = 1
  return(f)
}
ST <- S0 * exp((r - 0.5*sigma^2)*t + sigma*sqrt(t)*rnorm(M))
delta.pw <- mean(exp(-r*t)*f(ST,K)*ST/S0)
message(sprintf("True delta = %5.4f",delta.true))
message(sprintf("PW delta = %5.4f",delta.pw))

#Step 2 

#Vega Option A
M <- 10000 # number of paths - will keep constant
# Input params:
S0 <- 50
sigma <- 0.3
t <- 3/12
r <- 0.01
K <- 52
# true vega:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
vega.true <- S0*dnorm(d1)*sqrt(t)

f <- function(S,K){
  f <- rep(0,length=length(S))
  f[S>K] = 1
  return(f)
}

ST <- S0 * exp((r - 0.5*sigma^2)*t + sigma*sqrt(t)*rnorm(M))
vega.pw <- mean(exp(-r*t)*f(ST,K)*ST*((log(ST/S0)-(r+0.5*sigma^2)*t)/(sigma)))
message(sprintf("True vega = %5.4f",vega.true))
message(sprintf("PW vega = %5.4f",vega.pw))

#Step 2 

rm(list=ls())
set.seed(100)

# Define GBM Function:
myGBM <- function(M,N,t,X0,mu,sigma){
  musim <- mu - 0.5*sigma*sigma
  dt <- t/N
  X <- matrix(NA,ncol=(N+1),nrow=M)
  X[,1] <- X0
  for (i in 1:N){
    Z <- rnorm(M)
    X[,i+1] <- X[,i]*exp(musim*dt + sigma*sqrt(dt)*Z)
  }
  return(X)
}
M <- 10000
N <- 13
t <- 3/12
S0 <- 50
r <- 0.01
sigma <- 0.3
K <- 52
X <- myGBM(M,N,t,S0,r,sigma)
# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- pnorm(d1)
f <- function(S,K){
  f <- rep(0,length=length(S))
  f[S>K] = 1
  return(f)
}
Smean <- rowMeans(X[,c(2:N+1)],dims=1)
# ST <- S0 * exp((r - 0.5*sigma^2)*t + sigma*sqrt(t)*rnorm(M))
# 
# 
# ST <- exp(-r*t)*pmax(Smean-K,0)
# 
# ST <- S0 * exp((r - 0.5*sigma^2)*t + sigma*sqrt(t)*rnorm(M))

delta.pw <- mean(exp(-r*t)*f(Smean,K)*Smean/S0)
message(sprintf("True delta = %5.4f",delta.true))
message(sprintf("PW delta = %5.4f",delta.pw))

#Step 2 
M <- 10000
N <- 52
t <- 3/12
S0 <- 50
r <- 0.01
sigma <- 0.3
K <- 52
# true vega:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
vega.true <- S0*dnorm(d1)*sqrt(t)
f <- function(S,K){
  f <- rep(0,length=length(S))
  f[S>K] = 1
  return(f)
}
vegaGBM <- function(M,N,t,X0,mu,sigma){
  musim <- mu - 0.5*sigma*sigma
  dt <- t/N
  X <- matrix(NA,ncol=(N+2),nrow=M)
  X[,1] <- X0
  change<- matrix(NA, ncol=(N+2),nrow=M)
  change[,1]<-0
  
  for (i in 1:N){
    Z <- rnorm(M)
    X[,i+1] <- X[,i]*exp((mu-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)
    change[,i+1] <- X[,i+1]*(log((X[,i+1])/X0)-(mu+0.5*sigma^2)*(dt*i))/sigma
  }
  
  X[,(ncol(X))]<- rowMeans(change[,c(1:N)],dim=1)
  return(X)
}
X <- vegaGBM(M,N,t,S0,r,sigma)
Smean <- rowMeans(X[,c(2:N+1)],dims=1)


vega.pw <- mean(exp(-r*t)*f(Smean,K)*X[,(ncol(X))])
message(sprintf("True vega = %5.4f",vega.true))
message(sprintf("PW vega = %5.4f",vega.pw))



##Step 3 
# Delta of A

rm(list=ls())
M <- 10000 # number of paths - will keep constant
N <- 52
t <- 3/12
S0 <- 50
r <- 0.01
sigma <- 0.3
K <- 52

# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- pnorm(d1)
f <- function(S,K){
  f <- pmax(S-K,0)
  return(f)
}
Z <- rnorm(M)
ST <- S0 * exp((r - 0.5*sigma^2)*t + sigma*sqrt(t)*Z)
score <- Z/(S0*sigma*sqrt(t))
delta.lr <- mean(exp(-r*t)*f(ST,K)*score)
message(sprintf("True delta = %5.4f",delta.true))
message(sprintf("LR delta Euro Call = %5.4f",delta.lr))


##Step 3
# b) Vega of A

rm(list=ls())
M <- 10000 # number of paths - will keep constant
# Input params:
S0 <- 50
sigma <- 0.3
t <- 3/12
r <- 0.01
K <- 52
# true vega:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
vega.true <- S0*dnorm(d1)*sqrt(t)
f <- function(S,K){
  f <- pmax(S-K,0)
  return(f)
}
Z <- rnorm(M)
ST <- S0 * exp((r - 0.5*sigma^2)*t + sigma*sqrt(t)*Z)
score <- (Z*Z - 1)/sigma - Z*sqrt(t)
vega.lr <- mean(exp(-r*t)*f(ST,K)*score)
message(sprintf("True vega = %5.4f",vega.true))
message(sprintf("LR vega Euro call = %5.4f",vega.lr))


##Step 3 
# c) Delta of B
rm(list=ls())
M <- 10000 # number of paths - will keep constant
N <- 52
t <- 3/12
S0 <- 50
r <- 0.01
sigma <- 0.3
K <- 52
myGBM <- function(M,N,t,X0,mu,sigma){
  musim <- mu - 0.5*sigma*sigma
  dt <- t/N
  X <- matrix(NA,ncol=(N+1),nrow=M)
  X[,1] <- X0
  for (i in 1:N){
    Z <- rnorm(M)
    X[,i+1] <- X[,i]*exp(musim*dt + sigma*sqrt(dt)*Z)
    if (i==1){
      D <- Z
    }
  }
  ls <- list(X,D)
  return(ls)
}
X <- myGBM(M,N,t,S0,r,sigma)
Z <- X[[2]]
X <- X[[1]]
Smean <- rowMeans(X[,c(2:N+1)],dims=1)
dt <- t/N
# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- pnorm(d1)
f <- function(S,K){
  f <- pmax(S-K,0)
  return(f)
}
# Z <- rnorm(M)
score <- Z/(S0*sigma*sqrt(dt))

delta.lr <- mean(exp(-r*t)*f(Smean,K)*score)

message(sprintf("True delta = %5.4f",delta.true))
message(sprintf("LR delta Asian = %5.4f",delta.lr))



##Step 3 
# Vega of B
rm(list=ls())
M <- 10000 # number of paths - will keep constant
N <- 52
t <- 3/12
S0 <- 50
r <- 0.01
sigma <- 0.3
K <- 52

# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- pnorm(d1)
vega.true <- S0*dnorm(d1)*sqrt(t)

f <- function(S,K){
  f <- pmax(S-K,0)
  return(f)
}
myGBM <- function(M,N,t,X0,mu,sigma){
  musim <- mu - 0.5*sigma*sigma
  dt <- t/N
  X <- matrix(NA,ncol=(N+1),nrow=M)
  Z <- matrix(NA,nrow=M,ncol=N)
  X[,1] <- X0
  for (i in 1:N){
    Z[,i] <- rnorm(M)
    X[,i+1] <- X[,i]*exp(musim*dt + sigma*sqrt(dt)*Z[,i])
  }
  
  ls <- list(X,Z)
  return(ls)
}
dt <- t/N

X <- myGBM(M,N,t,S0,r,sigma)
Z <- X[[2]]
X <- X[[1]]
Smean <- rowMeans(X[,c(2:N+1)],dims=1)
score <- apply((Z*Z-1)/sigma - Z*(sqrt(dt)),1,sum)
vega.lr <- mean(exp(-r*t)*f(Smean,K)*score)
message(sprintf("True vega = %5.4f",vega.true))
message(sprintf("LR Vega Asian = %5.4f",vega.lr))

##Step 3 
# Delta of C

rm(list=ls())
M <- 10000 # number of paths - will keep constant
N <- 52
t <- 3/12
S0 <- 50
r <- 0.01
sigma <- 0.3
K <- 52

# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- pnorm(d1)
f <- function(S,K){
  f <- rep(0,length=length(S))
  f[S>K] = 1
  return(f)
}
Z <- rnorm(M)
ST <- S0 * exp((r - 0.5*sigma^2)*t + sigma*sqrt(t)*Z)
score <- Z/(S0*sigma*sqrt(t))
delta.lr <- mean(exp(-r*t)*f(ST,K)*score)

delta.lr <- mean(exp(-r*t)*f(ST,K)*score)

message(sprintf("True delta = %5.4f",delta.true))
message(sprintf("LR delta 3e) = %5.4f",delta.lr))

##Step 3 
#
rm(list=ls())
M <- 10000 # number of paths - will keep constant
N <- 52
t <- 3/12
S0 <- 50
r <- 0.01
sigma <- 0.3
K <- 52
b <- 60
# true delta:
d1 <- (log(S0/K) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
delta.true <- pnorm(d1)
vega.true <- S0*dnorm(d1)*sqrt(t)

# Barrier Call
f <- function(ST,Smax,K,b){
  f <- pmax(ST-K,0)
  f[Smax>=b]=0
  return(f)
}
dt <- t/N

myGBM <- function(M,N,t,X0,mu,sigma){
  musim <- mu - 0.5*sigma*sigma
  dt <- t/N
  X <- matrix(NA,ncol=(N+1),nrow=M)
  X[,1] <- X0
  for (i in 1:N){
    Z <- rnorm(M)
    X[,i+1] <- X[,i]*exp(musim*dt + sigma*sqrt(dt)*Z)
    if (i==1){
      Z1 <- Z
    }
  }
  
  ls <- list(X,Z1)
  return(ls)
}
X <- myGBM(M,N,t,S0,r,sigma)

Z1 <-X[[2]]
X<- X[[1]]
Smax <- apply(X,1,max)
ST <- X[,N+1]
score <- Z1/(S0*sigma*sqrt(dt))
delta.lr <- mean(exp(-r*t)*f(ST,Smax,K,b)*score)

message(sprintf("True delta = %5.4f",delta.true))
message(sprintf("LR Barrier delta 3f) = %5.4f",delta.lr))