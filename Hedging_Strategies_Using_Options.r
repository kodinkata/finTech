#Stop Loss Strategy:
# Case Study: bank sold a Euro Call option on 100,000 shares of a non-dividend paying stock for $300,000.  
            #Assume S0 = 49, K = 50, r = 5%, sigma = 20%, and T = 20 weeks (0.3846 years), 
            #and the expected return on the stock is 13%.
# S0 - current stock price
# K - strike price
# r - risk free rate
# sigma - stock's volatility
# T - time to expiration

# Step 1:
  # Results: Simulation of a stop loss strategy using 10,000 simulated paths using differente time steps (5, 4, 2, 1, 0.5, and 0.25 weeks)

# Step 2:
  # Step 1 using delta hedging strategy and comparing the results

# Step 3:
  # Step 2 using Asian call option where the payoff is the mean weekly closing price of the stock over the life of the option

#calc the cost of the strategy (similar to the BLS price)
#Compare to the actual payouts and calc the volatility (st dev)

#~~~Step 1~~~~~~

rm(list=ls())
# Define a function to simulate paths:
myGBM <- function(M,N,S0,r,t,sigma){
  # Return a vector if only terminal prices are needed:
  dt <- t/N
  paths <- matrix(NA,ncol=(N+1),nrow=M)
  paths[,1] <- S0
  
  for (i in 1:N){
    paths[,i+1] <- paths[,i]* 
      exp((r-0.5*sigma*sigma)*dt+
            sigma*sqrt(dt)*rnorm(M))
    
  }
  return(paths)
}
sigma <- 0.2
M <- 10000
r <- 0.05
K <- 50
S0 <- 49
t <- 20/52
mu <- 0.13
N <- c(4,5,10,20,40,80,1000)

stopLoss <- function(M,N,S0,K,r,sigma,t,mu){
  d1 <- (log(S0/K) + (r + sigma*sigma/2)*t)/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  BLS <- S0*pnorm(d1) - K*exp(-r*t)*pnorm(d2)
  X <- myGBM(M,N,S0,r,t, sigma)
  # Generate a matrix for in our out of the money:
  I.O <- matrix(0,ncol=N+1,nrow=M)
  I.O[X > K] = 1
  
  # Matrix of changes:
  d.I.O <- matrix(0,ncol=N+1,nrow=M)
  if (S0 > K){
    d.I.O[,1] <- 1
  }
  d.I.O[,c(2:(N+1))] <- (I.O[,c(2:(N+1))] - I.O[,c(1:N)])
  
  # Compute the cost of the strategy:
  # 1. intermediate buy/sells:
  C <- -d.I.O*X
  
  # 2. clean up at Maturity:
  C[which(X[,N+1] > K),N+1] <- C[which(X[,N+1] > K),N+1] + K
  
  # 3. sum the costs:
  disc <- matrix(exp(-r*seq(0,t,length=N+1)),ncol=1)
  PV <- C%*%disc
  
  # compute performace
  H.perf <- sqrt(var(PV))/BLS
  outlist <- list("H.perf"=H.perf,"PV"=PV,"BLS"=BLS)
  return(outlist)
  
}
H <- c(NA)
PV <-c(NA)
for (j in 1:length(N)){
  tmp <- stopLoss(M,N[j],S0,K,r,sigma,t,mu)
  H[j] <- tmp$H.perf
  PV[j] <- mean(tmp$PV)
}
print(H)
print(PV)
print(tmp$BLS)


#~~~Step 2~~~
delta.hedge <- function(M,N,S0,K,r,sigma,t,mu){
  d1 <- (log(S0/K) + (r + sigma*sigma/2)*t)/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  BLS <- S0*pnorm(d1) - K*exp(-r*t)*pnorm(d2)
  delta <- function(S,K,r,t,sigma){
    d <- pnorm((log(S/K) + 
                  (r + sigma*sigma/2)*t)/(sigma*sqrt(t)))
    return(d)
  }
  
  X <- deltas <- matrix(NA, ncol=N+1,nrow=M)
  X[,1]<-S0
  dt <- t/N
  for (i in 1:N){
    X[,i+1] <- X[,i]*exp((mu-0.5*sigma*sigma)*dt + sigma*sqrt(dt)*rnorm(M))
    ttm <- t-dt*(i-1)
    deltas[,i] <- delta(X[,i],K,r,ttm,sigma)
  }
  
  #Fill in terminal deltas (1/0):
  deltas[,N+1] <- delta(X[,N+1],K,r,0,sigma)
  
  #GENERATE a matrix of positions:
  CF <- matrix(NA,ncol=N+1, nrow = M)
  CF[,1] <- -deltas[,1] * S0
  for (i in 1:(N-1)){
    CF[,i+1] <- -1*(deltas[,i+1] - deltas[,i])*X[,i+1]
  }
  IN <- which(X[,N+1] >K)
  CF[IN,N+1] <- K -(1-deltas[IN,N])*X[IN,N+1]
  CF[-IN, N+1] <- deltas[-IN,N]*X[-IN,N+1]
  
  #3.Sum the costs:
  disc <- matrix(exp(-r*seq(0,t,length=N+1)),ncol=1)
  PV <- CF%*% disc
  
  #compute performance
  H.perf <- sqrt(var(PV))/BLS
  outlist <- list("H.perf"=H.perf,"PV"=PV,"BLS"=BLS)
  return(outlist)
}

N <- c(4,5,10,20,40,80,1000)
H <- c(NA)
PV <- c(NA)
for (j in 1:length(N)){
  tmp <- delta.hedge(M,N[j],S0,K,r,sigma,t,mu)
  H[j] <- tmp$H.perf
  PV[j] <- mean(tmp$PV)
}

print(H)
print(PV)
print(tmp$BLS)

#~~~Step 3~~~
rm(list=ls())
# ================
# Delta Hedge
# ================

delta.hedge.PW <- function(M,N,S0,K,r,sigma,t,mu){
  print(N)
  
  # Simulate the paths and deltas...also want to save Z's:
  X <- deltas <- Zmat <- matrix(NA,ncol=N+1,nrow=M)
  X[,1] <- S0
  Zmat[,1] <- 0
  dt <- t/N
  for (i in 1:N){
    Z <- rnorm(M)
    X[,i+1] <- X[,i]*exp((r - 0.5*sigma*sigma)*dt + 
                           sigma*sqrt(dt)*Z)
    
    ttm <- t - dt*(i-1)
    Zmat[,i+1] <- Z
  }
  
  # Asian option price:
  Sbar.A <- apply(X,1,mean) - K
  tmp <- which(Sbar.A > 0)
  Sbar.A[-tmp] <- 0
  P.Asian <- mean(Sbar.A)
  
  
  # ======================
  # Pathwise Estimator
  # ======================
  sumS <- t(apply(X,1,cumsum))
  S <- X
  # 3. Payout function:
  f <- function(S,K){
    f <- pmax(S - K,0)
    return(f)
  }
  
  # 4. Compute discount factors (for the deltas):
  disc <- exp(-r*dt*rev(c(0:N)))
  
  # 5. Finally compute the deltas over time:
  deltas.A <- matrix(NA,ncol=N+1,nrow=M)
  for (i in 1:N){
    if (i < N){
      Zmat.sum <- matrix(NA,nrow=M,ncol=N+1)
      Zmat.sum[,c((i+1):(N+1))] <- t(apply(Zmat[,c((i+1):(N+1))],1,cumsum))
    }
    else{
      Zmat.sum <- matrix(NA,nrow=M,ncol=N+1)
      Zmat.sum[,c((i+1):(N+1))] <- Zmat[,N+1]
    }
    RHS <- matrix(NA,nrow=M,ncol=N+1)
    for (k in i:N){
      RHS[,k+1] <- (k-i+1)*(r - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*Zmat.sum[,k+1]
    }
    for (j in 1:M){
      S2 <- matrix(NA,nrow=M,ncol=N+1)
      S2 <- exp(RHS + log(S[j,i])) # Matrix of stock prices going forwward
      S2[,i] <- S[j,i] # Initial price
      Sbar <- (apply(S2[,c(i:(N+1))],1,sum) + sumS[j,i] - S[j,i])/(N+1)
      dSdJ <- (apply(S2[,c(i:(N+1))],1,sum))/(N+1)
      I <- rep(0,length=M)
      tmp <- which(Sbar > K)
      I[tmp] <- 1
      deltas.A[j,i] <- disc[i]*mean(I*dSdJ/S2[j,i])
    }
  }
  
  
  # Generate a matrix of positions:
  CF <- matrix(NA,ncol=N+1,nrow=M)
  CF[,1] <- -deltas.A[,1]*S0
  for (i in 1:(N-1)){
    CF[,i+1] <- -1*(deltas.A[,i+1] - deltas.A[,i])*X[,i+1]
  }
  
  Xbar <- apply(X,1,mean)
  IN <- which(Xbar > K)
  
  CF[IN,N+1] <- K - Xbar[IN] + deltas.A[IN,N]*X[IN,N+1]
  CF[-IN,N+1] <- deltas.A[-IN,N]*X[-IN,N+1]
  
  # 3. sum the costs:
  disc <- matrix(exp(-r*seq(0,t,length=N+1)),ncol=1)
  PV <- CF%*%disc
  
  # compute performace
  H.perf <- sqrt(var(PV))/P.Asian
  outlist <- list("H.perf"=H.perf,"PV"=PV,"P.Asian"=P.Asian,
                  "deltas.A"=deltas.A,"CF"=CF)
  return(outlist)
}


# ===============
# Test the function:
S0 <- 49
K <- 50
r <- 0.05
sigma <- 0.20
t <- 20/52
mu <- 0.13
M <- 500
N <- 12
h <- 0.1
set.seed(2500)
PW.1 <- delta.hedge.PW(M,N,S0,K,r,sigma,t,mu)
PW.means <- apply(PW.1$deltas.A,2,mean)
PW.sd <- sqrt(apply(PW.1$deltas.A,2,var))

print(PW.means)
print(PW.sd)
