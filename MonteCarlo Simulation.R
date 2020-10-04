rm(list=ls())

# Step 1: computing N Geometric Brownian Motion paths

# Step 2 & 3: using Monte Carlo simulation to price a European option (Step 2: Call, Step 3: Put): current price of $50, a strike price of $53, volatility of 30% and a risk free rate of 5%  

# Step 4: Using Monte Carlo simulation to price an Asian Option: current price of $40, strike price K of $40, time to maturity of 12 months, volatility of 40% and assume the risk free rate of 10%

# Step 5: Using Monte Carlo simulation to price Barrier Options:
  # a) Option 1: Down-and-out put barrier option on a stock with S0 = 50, K = 52, r = 0.1, T = 5/12, sigma = 0.4, and a barrier of 40
  # b) Option 2: Up-and-out call barrier option on a stock with S0 = 50, K = 52, r = 0.1, T = 5/12, sigma = 0.4, and a barrier of 60
  # c) Option 3: Down-and-in-put barrier option with S0 = $50, K = $52, sigma = 0.4, T = 5/12, r = 0.1, and a barrier of $46
# ---------------------------------
# Step 1:
library("ggplot2")
library("reshape")
myGBM <- function(M,N,S0,r,t,sigma){
  dt <- t/N
  # Return a vector if only terminal prices are needed:
  if (N == 1){
    paths <- S0*exp((r-0.5*sigma*sigma)*dt+
                      sigma*sqrt(dt)*rnorm(M))
  }
  else{
    paths <- matrix(NA,ncol=(N+1),nrow=M) 
    for (i in 1:N){
      paths[,1] <- S0
      paths[,i+1] <- paths[,i]* 
        exp((r-0.5*sigma*sigma)*dt+
              sigma*sqrt(dt)*rnorm(M))
    }
  }
  return(paths)
}

# Test the function:
S0 <- 50
r <- 0.05
t <- 1/2
N <- 126
M <- 10000
sigma <- 0.3
paths <- myGBM(M,N,S0,r,t,sigma)

time <- as.matrix(seq(0,t,length=N+1),ncol=1)
S <- t(paths)
df <- as.data.frame(cbind(time,S))

# Plot some paths:
mP <- 20 #limit the number of paths being plotted
Anames <- rep("",mP)
for (i in 1:mP){
  Anames[i] <- paste("A",i,sep="")
}
names(df) <- c("time", Anames)

# This creates a new data frame with columns x, variable and value
# x is the id, variable holds each of our timeseries designation
df.melted <- melt(df[,c(1:mP)], id = "time")

p <- ggplot(data = df.melted, aes(x = time, y = value, color = variable)) +
  geom_line()  + guides(fill=FALSE, color=FALSE)

p <- p + labs(x = "Time", y ="S(t)",title="GBM")
print(p)

# ---------------------------------
# Step 2 and 3:
# Compute call and puts together:
mcOpts <- function(S0,K,r,t,sigma,M,N=1,alpha=0.95){
  dt <- t
  z <- qnorm((1+alpha)/2)
  paths <- myGBM(M,N=1,S0,r,t,sigma)
  hc <- function(S){
    hc <- pmax(S-K,0)
  }
  hp <- function(S){
    hp <- pmax(K-S,0)
  }
  c0 <- exp(-r*t)*hc(paths)
  p0 <- exp(-r*t)*hp(paths)
  c <- mean(c0)
  p <- mean(p0)
  sdc <- sd(c0)
  sdp <- sd(p0)
  cCI <- c(c - sdc/sqrt(M)*z,c + sdc/sqrt(M)*z)
  pCI <- c(p - sdp/sqrt(M)*z,p + sdp/sqrt(M)*z)
  results <- list("c"=c,"p"=p,"cCI"=cCI,"pCI"=pCI)
  return(results)
}

S0 <- 50
K <- 53
sigma <- 0.30
t <- 0.5
r <- 0.05
N <- 1
M <- 10000

soln45 <- mcOpts(S0,K,r,t,sigma,M,N=1,alpha=0.95)

# Check Black-Scholes:
myBLS <- function(S0,K,r,t,sigma,print=1){
  opts <- c(2)
  d1 <- (log(S0/K) + (r + sigma^2/2)*t) / (sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  opts[1] <- S0 * pnorm(d1) - K*exp(-r*t)*pnorm(d2)
  opts[2] <- -S0*pnorm(-d1) + K*exp(-r*t)*pnorm(-d2)
  if (print == 1){
    message(sprintf("The price of the call option is $%5.4f",opts[1]))
    message(sprintf("The price of the put option is $%5.4f",opts[2]))
  }
  return(opts)
}
opts<-myBLS(S0,K,r,t,sigma,print=1)

# Run the expriment to see how SE's converge:
dat <- as.data.frame(matrix(NA,nrow=100,ncol=7))
M <- seq(10,10000,length=1000)
for (i in 1:length(M)){
  tmp <- mcOpts(S0,K,r,t,sigma,M[i],N=1,alpha=0.95)
  dat[i,] <- c(M[i],tmp$c,tmp$cCI,tmp$p,tmp$pCI)
}

mynames <- c("M","c","cLO","cHI","p","pLO","pHI")
names(dat) <- mynames

p <- ggplot(dat, aes(x=M,y = value, color = variable)) + 
  geom_line(aes(y=c, col="c"),color="black") + 
  geom_line(aes(y=cLO, col="cLO"),color="red") + 
  geom_line(aes(y=cHI, col="cHI"),color="red") + 
  
  geom_line(aes(y=p, col="p"),color="black") + 
  geom_line(aes(y=pLO, col="pLO"),color="blue") + 
  geom_line(aes(y=pHI, col="pHI"),color="blue") + 
  
  ggtitle("call and Put Option Estimates with 95% CI's")
print(p)

# Note that the CI's are changeing with the estimates.
# It is better to plot the "true CI's:
# We'll average across a huge sample to get the "true" CI's:
alpha <- 0.95
dat2 <- as.data.frame(matrix(NA,nrow=100,ncol=2))
M <- 1e6
for (i in 1:100){
  tmp <- mcOpts(S0,K,r,t,sigma,M,N=1,alpha=0.95)
  dat2[i,1] <- (tmp$cCI[2] - tmp$c)*sqrt(M)/qnorm((1+alpha)/2)
  dat2[i,2] <- (tmp$pCI[2] - tmp$p)*sqrt(M)/qnorm((1+alpha)/2)
}
sdc <- mean(dat2[,1])
sdp <- mean(dat2[,2])

dat$cLOt <- opts[1] - qnorm((1+alpha)/2)*sdc/sqrt(dat$M)
dat$cHIt <- opts[1] + qnorm((1+alpha)/2)*sdc/sqrt(dat$M)
dat$pLOt <- opts[2] - qnorm((1+alpha)/2)*sdp/sqrt(dat$M)
dat$pHIt <- opts[2] + qnorm((1+alpha)/2)*sdp/sqrt(dat$M)

mynames <- c("M","c","cLO","cHI","p","pLO","pHI",
             "cLOt","cHIt","pLOt","pHIt")
names(dat) <- mynames

p2 <- ggplot(dat, aes(x=M,y = value, color = variable)) + 
  geom_line(aes(y=c, col="c"),color="black") + 
  geom_line(aes(y=cLOt, col="cLOt"),color="red") + 
  geom_line(aes(y=cHIt, col="cHIt"),color="red") + 
  
  geom_line(aes(y=p, col="p"),color="black") + 
  geom_line(aes(y=pLOt, col="pLOt"),color="blue") + 
  geom_line(aes(y=pHIt, col="pHIt"),color="blue") + 
  
  ggtitle("call and Put Option Estimates with 95% CI's")
print(p2)


# ----------------------
# Step 4:
# Asian Option:
S0 <- 40
K <- 40
sigma <- 0.40
r <- 0.1
t <- 1
M <- 10000
N <- 12
paths <- myGBM(M,N,S0,r,t,sigma)

Smean <- rowMeans(paths[,c(2:(N+1))], dims = 1)
a0 <- exp(-r*t)*pmax(Smean - K,0)
a <- mean(a0)
aLO <- a - qnorm((1+alpha)/2)*sd(a0)/sqrt(M)
aHI <- a + qnorm((1+alpha)/2)*sd(a0)/sqrt(M)

message(sprintf("Asian Option: $%5.4f, ($%5.4f, $%5.4f)",a,aLO,aHI))

# Compare to BLS:
myBLS(S0,K,r,t,sigma,print=1)


# ------------------------
# Step 5:
S0 <- 50
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4
Sb <- 40
M <- 10000
N <- ceiling(252*5/12)
paths <- myGBM(M,N,S0,r,t,sigma)

# Part a:
h <- function(paths,K){
  Smin <- apply(paths,1,min)
  Smin[Smin <= Sb] <- 0
  Smin[Smin > Sb] <- 1
  St <- paths[,(N+1)]
  h <- pmax(K - St,0)*Smin
}

f0 <- exp(-r*t)*h(paths,K)
dop <- mean(f0)
dop_se <- sd(f0)/sqrt(M)
message(sprintf("Down-and-Out Put: %5.4f (%5.4f)",dop,dop_se))


# Part b:
Sb <- 60
h <- function(paths,K){
  Smax <- apply(paths,1,max)
  Smax[Smax < Sb] <- 1
  Smax[Smax >= Sb] <- 0
  St <- paths[,(N+1)]*Smax
  h <- pmax(St - K,0)
}

f0 <- exp(-r*t)*h(paths,K)
uoc <- mean(f0)
uoc_se <- sd(f0)/sqrt(M)
message(sprintf("Up-and-Out Call: %5.4f (%5.4f)",uoc,uoc_se))

# Part c:
Sb <- 40
h <- function(paths,K){
  Smin <- apply(paths,1,min)
  Smin[Smin <= Sb] <- 1
  Smin[Smin > Sb] <- 0
  St <- paths[,(N+1)]
  h <- pmax(K - St,0)*Smin
}

f0 <- exp(-r*t)*h(paths,K)
dip <- mean(f0)
dip_se <- sd(f0)/sqrt(M)
message(sprintf("Down-and-in Put: %5.4f (%5.4f)",dip,dip_se))

# Plain put:
h <- function(paths,K){
  St <- paths[,(N+1)]
  h <- pmax(K - St,0)
}

f0 <- exp(-r*t)*h(paths,K)
put <- mean(f0)
put_se <- sd(f0)/sqrt(M)
message(sprintf("Plain Put: %5.4f (%5.4f)",put,put_se))

