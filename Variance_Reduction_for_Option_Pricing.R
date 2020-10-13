#Price different options using variance reduction techniques

#Assumption: stock price evolves under a geometric Brownian Motion process 
#Given characteristics: sigma = 30%, S0=$50, K = $52, and three months to maturity.  Assume a risk free rate of 1%.
#sigma - stock's volatility, S0- current stock price, K - strike price

# Option A: European call option 
# Option B: Asian option with payoff equal to max{Average(S0,ST) - K, 0}
# Option C: Digital option whose payoff is equal to $1 if ST > K and $0 otherwise
# Option D: An up and out barrier option with a barrier of $60.  The payoff is max[S(T) - K,0] if max{S(t)} < barrier and $0 otherwise.


#Step 1: Using Control Variate Method & estimate the computational savings by implementing this method
  # a) Option A using the underlying stock as the control
  # b) Option B using the underlying stock as the control
  # c) Option B using a plain vanilla European option as the control (i.e., Black-Scholes)
  # d) Option B using geometric Asian option as the control
  # e) Option C using the underlying asset as the control
  # f) Option D using the underlying stock as the control
  # g) Option D using a plain vanilla European option as the control (i.e., Black-Scholes)
  
#Step 2: Using Antithetic Sampling Estimator & estimate the computational savings by implementing this method
  # a) Option A 
  # b) Option B
  # c) Option C
  # d) Option D

#Step 3: Stratified Random Sampling - stratified random sampling breaks down the population 
        # into subgroups and organizes them by similar traits, characteristics, and behavior


#~~~~~~~~~~Step 1 ~~~~~~~~~~~
#a)Vanilla Euro Call with underlying stock as the control
rm(list=ls())
M <- 10000
S0 <- 50
sigma <- 0.3
r <- 0.01
t <- 3/12
K <- 52

# Implment Control variate method using underlying asset to
# price a plain vanilla call:
h <- function(S,K){
  h <- pmax(S-K,0)
}
set.seed(100)
# 1. use a "pilot" sample to estimate b*:
m <- floor(0.05*M) #use 5% sample for b*
# Generate a vector X of the stock price at time T (i.e.,
# the control variable)
X <- S0*exp((r - 0.5*sigma*sigma)*t +
              sigma*sqrt(t)*rnorm(m))
# Gnenerate the option payouts to estimate b*:
Y <- exp(-r*t)*h(X,K)
# Plot Y vs. X:
plot(x=X,y=Y,main="Regression for K = 105")
abline(lm(Y ~ X))
# Estimate b*:
b.star <- -cov(X,Y)/var(X)
# Now compute the price:
X.new <- S0*exp((r - 0.5*sigma*sigma)*t +
                  sigma*sqrt(t)*rnorm(M-m))
Y.new <- exp(-r*t)*h(X.new,K)
CV <- Y.new + b.star*(X.new - S0*exp(r*t))
CV.c <- mean(CV) #option price
CV.CI <- c(CV.c - qnorm(0.975)*sqrt(var(CV)/(M-m)),
           CV.c + qnorm(0.975)*sqrt(var(CV)/(M-m)))
# Compare to MC estimator:
set.seed(100)
X.MC <- S0*exp((r - 0.5*sigma*sigma)*t +
                 sigma*sqrt(t)*rnorm(M-m))
Y.MC <- exp(-r*t)*h(X.MC,K)
CV.c.MC <- mean(Y.MC) #option price
CV.CI.MC <- c(CV.c.MC - qnorm(0.975)*sqrt(var(Y.MC)/(M)),
              CV.c.MC + qnorm(0.975)*sqrt(var(Y.MC)/(M)))
message(sprintf("Basic MC estimator: %5.4f (%5.4f, %5.4f)",
                CV.c.MC,CV.CI.MC[1],CV.CI.MC[2]))
message(sprintf("Control Variate estimator: %5.4f (%5.4f, %5.4f)",
                CV.c,CV.CI[1],CV.CI[2]))
message(sprintf("Variance reduction ratio: %5.4f",
                var(CV)/var(Y.MC)))
message(sprintf("Correlation rho: %5.4f",
                cor(X,Y)))

#b) Asian option with underlying stock


# Implment Control variate method using underlying asset to
# price a plain vanilla call:
h <- function(S,K){
  h <- pmax(S-K,0)
}
# Three potential Control variables:
set.seed(200)
m <- floor(0.05*M)
N <-13
# Generate paths first to get b*:
X <- matrix(NA,ncol=N+1,nrow=m)
X[,1] <- S0
dt <- t/N
for (i in 1:N){
  X[,i+1] <- X[,i]*exp((r - 0.5*sigma*sigma)*dt +
                         sigma*sqrt(dt)*rnorm(m))
}
# Save terminal prices (to make the notation cleaner):
XT <- X[,N+1]
# 1. Control Variate 1: underlying price:
CV.1 <- XT
E.CV.1 <- S0*exp(r*t) # expected value

#c) Asian option with Euro Call control
# 2. Control Variate 2: plain vanilla Black Scholes:
CV.2 <- exp(-r*t)*h(XT,K)
d1 <- (log(S0/K) + (r + sigma*sigma/2)*t)/(sigma*sqrt(t))
d2 <- d1 - sigma*sqrt(t)
E.CV.2 <- S0*pnorm(d1) - K*exp(-r*t)*pnorm(d2) # expected value

#d) Asian option with geometric Asian option control
# 3. Control Variate 3: Geometric Asian Option
Xbar.Geom <- apply(X,1,prod)^(1/(N+1))
CV.3 <- exp(-r*t)*h(Xbar.Geom,K)
tvec <- seq(0,t,length=N+1)
TT <- mean(tvec)
ivec <- c(1:(N+1))
sigma.bar.2 <- sigma^2/((N+1)^2*TT) * sum((2*ivec-1)*rev(tvec))
delta <- 0.5*(sigma^2 - sigma.bar.2)
d <- (log(S0/K) + (r - delta + 0.5*sigma.bar.2)*t)/
  (sqrt(sigma.bar.2)*sqrt(t))
E.CV.3 <- exp(-delta*TT)*S0*pnorm(d) -
  exp(-r*TT)*K*pnorm(d - sqrt(sigma.bar.2)*sqrt(TT))

# Generate the Arithmatic mean paths to get b*'s:
Xbar.Arit <- apply(X,1,mean)
Y <- exp(-r*t)*h(Xbar.Arit,K)
cor(Y,CV.3)


# Estimate b*'s for each control variable:
b.star <- c(NA,NA,NA)
b.star[1] <- -cov(CV.1,Y)/var(CV.1)
b.star[2] <- -cov(CV.2,Y)/var(CV.2)
b.star[3] <- -cov(CV.3,Y)/var(CV.3)
# Estimate the option values:
X.new <- matrix(NA,ncol=N+1,nrow=M-m)
X.new[,1] <- S0
dt <- t/N
for (i in 1:N){
  X.new[,i+1] <- X.new[,i]*exp((r - 0.5*sigma*sigma)*dt +
                                 sigma*sqrt(dt)*rnorm(m))
}
Xbar.Arit.new <- apply(X.new,1,mean)
Y.new <- exp(-r*t)*h(Xbar.Arit.new,K)

# COMPARE TO RAW MC ESTIMATOR:
Y.mean <- mean(Y.new)
Y.var <- var(Y.new)
# 1. underlying stock:
XT.new <- X.new[,N+1]
CV.1.new <- XT.new
theta.CV.1 <- Y.new + b.star[1]*(CV.1.new - E.CV.1)
CV.1.mean <- mean(theta.CV.1)
CV.1.var <- var(theta.CV.1)
CV.1.var.ratio <- CV.1.var/Y.var
message(sprintf("Control Variate 1:\n\tMean = %5.4f\n\tVariance = %5.4f\n\tVar. Ratio = %5.4f",
                CV.1.mean,CV.1.var,CV.1.var.ratio))
# 2. Black Scholes:
CV.2.new <- exp(-r*t)*h(XT.new,K)
theta.CV.2 <- Y.new + b.star[2]*(CV.2.new - E.CV.2)
CV.2.mean <- mean(theta.CV.2)
CV.2.var <- var(theta.CV.2)
CV.2.var.ratio <- CV.2.var/Y.var
message(sprintf("Control Variate 2:\n\tMean = %5.4f\n\tVariance = %5.4f\n\tVar. Ratio = %5.4f",
                CV.2.mean,CV.2.var,CV.2.var.ratio))
# 3. Geom Asian:
Xbar.Geom.new <- apply(X.new,1,prod)^(1/(N+1))
CV.3.new <- exp(-r*t)*h(Xbar.Geom.new,K)
theta.CV.3 <- Y.new + b.star[3]*(CV.3.new - E.CV.3)
CV.3.mean <- mean(theta.CV.3)
CV.3.var <- var(theta.CV.3)
CV.3.var.ratio <- CV.3.var/Y.var
message(sprintf("Control Variate 3:\n\tMean = %5.4f\n\tVariance = %5.4f\n\tVar. Ratio = %5.4f",
                CV.3.mean,CV.3.var,CV.3.var.ratio))

#e) Digital option with underlying stock
set.seed(100)
hDig <- function(S,K){
  h <- pmax(S-K,0)/(S-K)
}
# Estimate the option values:
X.MC <- S0*exp((r - 0.5*sigma*sigma)*t +
                 sigma*sqrt(t)*rnorm(M-m))
# Gnenerate the option payouts to estimate b*:
Y.MC <- exp(-r*t)*hDig(X.MC,K)

CV.c <- mean(Y.MC) #option price
CV.CI <- c(CV.c - qnorm(0.975)*sqrt(var(Y.MC)/(M)),
           CV.c + qnorm(0.975)*sqrt(var(Y.MC)/(M)))

message(sprintf("Digital Basic MC estimator: %5.4f (%5.4f, %5.4f)",
                CV.c,CV.CI[1],CV.CI[2]))

set.seed(100)
m <- floor(0.05*M)

# Generate a vector X of the stock price at time T (i.e.,
# the control variable)
X <- S0*exp((r - 0.5*sigma*sigma)*t +
              sigma*sqrt(t)*rnorm(m))
# Gnenerate the option payouts to estimate b*:
Y <- exp(-r*t)*hDig(X,K)

# Estimate b*:
b.star <- -cov(X,Y)/var(X)
# Now compute the price:
X.new <- S0*exp((r - 0.5*sigma*sigma)*t +
                  sigma*sqrt(t)*rnorm(M-m))
Y.new <- exp(-r*t)*hDig(X.new,K)
CV <- Y.new + b.star*(X.new - S0*exp(r*t))
CV.c <- mean(CV) #option price
CV.CI <- c(CV.c - qnorm(0.975)*sqrt(var(CV)/(M-m)),
           CV.c + qnorm(0.975)*sqrt(var(CV)/(M-m)))

message(sprintf("Digital option Control Variate estimator: %5.4f (%5.4f, %5.4f)",
                CV.c,CV.CI[1],CV.CI[2]))

message(sprintf("Digital Option Variance reduction ratio: %5.4f",
                var(CV)/var(Y.MC)))




#f)Barrier Option with underlying stock
bar <- 60
set.seed(100)
m <- floor(0.05*M)
h <- function(S,K){
  h <- pmax(S - K,0)
}

hBar <- function(S,ST,bar,K){
  mat <- matrix(0,nrow=length(S))
  mat[ST <bar]<-1
  h <- (pmax(S-K,0))*mat
}


# Generate paths first to get b*:
dt <- t/N

X <- matrix(NA,ncol=N+1,nrow=m)
X[,1] <- S0
for (i in 1:N){
  X[,i+1] <- X[,i]*exp((r - 0.5*sigma*sigma)*dt + 
                         sigma*sqrt(dt)*rnorm(m))
}
# Save terminal prices (to make the notation cleaner):
XT <- X[,N+1]
Smax <- apply(X,1,max)

# 1. Control Variate 1: underlying price:
CV.1 <- XT
E.CV.1 <- S0*exp(r*t) # expected value

# 2. Control Variate 2: plain vanilla Black Scholes:
CV.2 <- exp(-r*t)*h(XT,K)
d1 <- (log(S0/K) + (r + sigma*sigma/2)*t)/(sigma*sqrt(t))
d2 <- d1 - sigma*sqrt(t)
E.CV.2 <- S0*pnorm(d1) - K*exp(-r*t)*pnorm(d2) # expected value

# Estimate barrier option payoffs:
Y <- exp(-r*t)*hBar(XT,Smax,bar,K)

# Estimate b*'s for each control variable:
b.star <- c(NA,NA)
b.star[1] <- -cov(CV.1,Y)/var(CV.1)
b.star[2] <- -cov(CV.2,Y)/var(CV.2)


# Estimate the option values:
X.new <- matrix(NA,ncol=N+1,nrow=M-m)
X.new[,1] <- S0
dt <- t/N
for (i in 1:N){
  X.new[,i+1] <- X.new[,i]*exp((r - 0.5*sigma*sigma)*dt + 
                                 sigma*sqrt(dt)*rnorm(m))
}

XT <- X.new[,N+1]
Smax <- apply(X.new,1,max)
Y.new <- exp(-r*t)*hBar(XT,Smax,bar,K)

# COMPARE TO RAW MC ESTIMATOR:
Y.mean <- mean(Y.new)
Y.var <- var(Y.new)

#f) Underlying stock:
XT.new <- X.new[,N+1]
CV.1.new <- XT.new
theta.CV.1 <- Y.new + b.star[1]*(CV.1.new - E.CV.1)
CV.1.mean <- mean(theta.CV.1)
CV.1.var <- var(theta.CV.1)
CV.1.var.ratio <- CV.1.var/Y.var
message(sprintf("Control Variate 1:\n\tMean = %5.4f\n\tVariance = %5.4f\n\tVar. Ratio = %5.4f",
                CV.1.mean,CV.1.var,CV.1.var.ratio))

#g) EuroCall:
CV.2.new <- exp(-r*t)*h(XT.new,K)
theta.CV.2 <- Y.new + b.star[2]*(CV.2.new - E.CV.2)
CV.2.mean <- mean(theta.CV.2)
CV.2.var <- var(theta.CV.2)
CV.2.var.ratio <- CV.2.var/Y.var
message(sprintf("Control Variate 2:\n\tMean = %5.4f\n\tVariance = %5.4f\n\tVar. Ratio = %5.4f",
                CV.2.mean,CV.2.var,CV.2.var.ratio))

#~~~~~~~~~~ Step 2~~~~~~~~~~~~~~~~

#Euro Call
rm(list=ls())
# Antithetic Variates:
N <- 1
M <- 10000
S0 <- 50
sigma <- 0.3
r <- 0.01
t <- 3/12
K <- 52
# Implment Control variate method using underlying asset to
# price a plain vanilla call:
h <- function(S,K){
  h <- pmax(S-K,0)
}
# -----------------------------
# 1. Test the antithetic method for a plain vanilla call:
set.seed(100)
# Generate M Z values:
Z <- rnorm(M)
Z.A <- -Z
# Regular MC estimator:
start.time <- Sys.time()

Y <- exp(-r*t)*h(S0*exp((r - 0.5*sigma*sigma)*t + sigma*sqrt(t)*Z),K)
Y.var <- var(Y)
Y.mean <- mean(Y)
Y.se <- sqrt(Y.var/M)
stop.time <- Sys.time()
time <- stop.time-start.time
# Antithetic MC estimator:
start.time2 <- Sys.time()

Y.1 <- exp(-r*t)*h(S0*exp((r - 0.5*sigma*sigma)*t + sigma*sqrt(t)*Z),K)
Y.2 <- exp(-r*t)*h(S0*exp((r - 0.5*sigma*sigma)*t + sigma*sqrt(t)*Z.A),K)
Y.A <- (Y.1 + Y.2) / 2
Y.A.var <- var(Y.A)
Y.A.mean <- mean(Y.A)
Y.A.se <- sqrt(Y.A.var/M)
stop.time2 <- Sys.time()
time2 <- stop.time2-start.time2
message(sprintf("Speed MC\n\tRegular = %5.4f\n\tAntithetic = %5.4f",
                time,time2))
message(sprintf("Regular MC\n\tPrice = %5.4f\n\tVariance = %5.4f\n\tS.E. = %5.4f",
                Y.mean,Y.var,Y.se))
message(sprintf("Antithetic MC\n\tPrice = %5.4f\n\tVariance = %5.4f\n\tS.E. = %5.4f",
                Y.A.mean,Y.A.var,Y.A.se))



#Asian Option
rm(list=ls())
# Antithetic Variates:
N <- 13
M <- 10000
S0 <- 50
sigma <- 0.3
r <- 0.01

t <- 3/12
K <- 52
# Implment Control variate method using underlying asset to
# price a plain vanilla call:
h <- function(S,K){
  h <- pmax(S-K,0)
}
# -----------------------------
# 1. Test the antithetic method for a plain vanilla call:
set.seed(100)
# Generate M Z values:
Z <- matrix(rnorm(M*N),ncol=N,nrow=M)
Z.A <- -Z
X <- matrix(NA,ncol=N+1,nrow=M)
X[,1] <- S0
dt <- t/N
for (i in 1:N){
  X[,i+1] <- X[,i]*exp((r - 0.5*sigma*sigma)*dt +
                         sigma*sqrt(dt)*Z[,i])
}
# Regular MC estimator:
start.time <- Sys.time()

Y <- exp(-r*t)*h(apply(X,1,mean),K)
Y.var <- var(Y)
Y.mean <- mean(Y)
Y.se <- sqrt(Y.var/M)
stop.time <- Sys.time()
time <- stop.time-start.time
# Antithetic MC estimator:
start.time2 <- Sys.time()

X.neg <- matrix(NA,ncol=N+1,nrow=M)
X.neg[,1] <- S0
dt <- t/N
for (i in 1:N){
  X.neg[,i+1] <- X.neg[,i]*exp((r - 0.5*sigma*sigma)*dt +
                                 sigma*sqrt(dt)*Z.A[,i])
}
Y.1 <- exp(-r*t)*h(apply(X,1,mean),K)
Y.2 <- exp(-r*t)*h(apply(X.neg,1,mean),K)
Y.A <- (Y.1 + Y.2) / 2
Y.A.var <- var(Y.A)
Y.A.mean <- mean(Y.A)
Y.A.se <- sqrt(Y.A.var/M)
stop.time2 <- Sys.time()
time2 <- stop.time2-start.time2
message(sprintf("Speed\n\tRegular = %5.4f\n\tAntithetic = %5.4f",
                time,time2))

message(sprintf("Regular MC\n\tPrice = %5.4f\n\tVariance = %5.4f\n\tS.E. = %5.4f",
                Y.mean,Y.var,Y.se))

message(sprintf("Antithetic MC\n\tPrice = %5.4f\n\tVariance = %5.4f\n\tS.E. = %5.4f",
                Y.A.mean,Y.A.var,Y.A.se))


#Digital
set.seed(100)
# Generate M Z values:
Z <- matrix(rnorm(M),ncol=N,nrow=M)
Z.A <- -Z
X <- S0*exp((r - 0.5*sigma*sigma)*t + 
              sigma*sqrt(t)*Z)

h <- function(S,K){
  h <- pmax(S-K,0)/(S-K)
}

# Regular MC estimator:
start.time <- Sys.time()
Y <- exp(-r*t)*h(X,K)
Y.var <- var(Y)
Y.mean <- mean(Y)
Y.se <- sqrt(Y.var/M)
stop.time <- Sys.time()
time <- stop.time-start.time
# Antithetic MC estimator:
start.time2 <- Sys.time()
X.neg <- S0*exp((r - 0.5*sigma*sigma)*t + 
                  sigma*sqrt(t)*Z.A)

Y.1 <- exp(-r*t)*h(X,K)
Y.2 <- exp(-r*t)*h(X.neg,K)
Y.A <- (Y.1 + Y.2) / 2
Y.A.var <- var(Y.A)
Y.A.mean <- mean(Y.A)
Y.A.se <- sqrt(Y.A.var/M)
stop.time2 <- Sys.time()
time2 <- stop.time2-start.time2
message(sprintf("Speed\n\tRegular = %5.4f\n\tAntithetic = %5.4f",
                time,time2))

message(sprintf("Digital Regular MC\n\tPrice    = %5.4f\n\tVariance = %5.4f\n\tS.E.     = %5.4f",
                Y.mean,Y.var,Y.se))

message(sprintf("Digital Antithetic MC\n\tPrice    = %5.4f\n\tVariance = %5.4f\n\tS.E.     = %5.4f",
                Y.A.mean,Y.A.var,Y.A.se))

#d) Barrier
rm(list=ls())
N <- 13
M <- 10000
S0 <- 50
sigma <- 0.3
r <- 0.01

t <- 3/12
K <- 52
bar <- 60
hBar <- function(S,ST,bar,K){
  mat <- matrix(0,nrow=length(S))
  mat[ST <bar]<-1
  h <- (pmax(S-K,0))*mat
}

Z <- matrix(rnorm(M*N),ncol=N,nrow=M)
Z.A <- -Z
X <- matrix(NA,ncol=N+1,nrow=M)
X[,1] <- S0
dt <- t/N
for (i in 1:N){
  X[,i+1] <- X[,i]*exp((r - 0.5*sigma*sigma)*dt +
                         sigma*sqrt(dt)*Z[,i])
}
# Regular MC estimator:
start.time <- Sys.time()

Y <- exp(-r*t)*hBar(X[,N+1],apply(X,1,max),bar,K)
Y.var <- var(Y)
Y.mean <- mean(Y)
Y.se <- sqrt(Y.var/M)
stop.time <- Sys.time()
time <- stop.time-start.time
# Antithetic MC estimator:
start.time2 <- Sys.time()

X.neg <- matrix(NA,ncol=N+1,nrow=M)
X.neg[,1] <- S0
dt <- t/N
for (i in 1:N){
  X.neg[,i+1] <- X.neg[,i]*exp((r - 0.5*sigma*sigma)*dt +
                                 sigma*sqrt(dt)*Z.A[,i])
}
Y.1 <- exp(-r*t)*hBar(X[,N+1],apply(X,1,max),bar,K)
Y.2 <- exp(-r*t)*hBar(X.neg[,N+1],apply(X.neg,1,max),bar,K)
Y.A <- (Y.1 + Y.2) / 2
Y.A.var <- var(Y.A)
Y.A.mean <- mean(Y.A)
Y.A.se <- sqrt(Y.A.var/M)
stop.time2 <- Sys.time()
time2 <- stop.time2-start.time2
message(sprintf("Speed\n\tRegular = %5.4f\n\tAntithetic = %5.4f",
                time,time2))

message(sprintf("Barrier,Regular MC\n\tPrice = %5.4f\n\tVariance = %5.4f\n\tS.E. = %5.4f",
                Y.mean,Y.var,Y.se))

message(sprintf("Barrier,Antithetic MC\n\tPrice = %5.4f\n\tVariance = %5.4f\n\tS.E. = %5.4f",
                Y.A.mean,Y.A.var,Y.A.se))
#~~~~~~~~Step 3~~~~~~~~
#a)Standard Uniform
rm(list=ls())
set.seed(100)
# set number of bins:
k <- 100
# number of draws per bin:
n.A <- 5
# Total number of draws:
M <- k*n.A
# Regular simulation:
Z <- runif(M)
# Stratified Sampling:
# 1. Generate bin locations:
a <- seq(0,1,length=k+1)
# 2. Sample for each bin:
Z.S <- c(NA)
for (i in 1:k){
  u <- runif(n.A)
  v <- a[i] + u*(a[i+1] - a[i])
  Z.S <- c(Z.S,v)
}
# Plot the results:
par(mfrow=c(2,3))
hist(Z,25,main="Standard Uniform Regular")
hist(Z.S,25,main="Standard Uniform Strat")

#b) Standard Normal
# set number of bins:
k <- 100
# number of draws per bin:
n.A <- 5
# Total number of draws:
M <- k*n.A
# Regular simulation:
Z <- rnorm(M)
# Stratified Sampling:
# 1. Generate bin locations:
a <- seq(0,1,length=k+1)
# 2. Sample for each bin:
Z.S <- c(NA)
for (i in 1:k){
  u <- runif(n.A)
  v <- a[i] + u*(a[i+1] - a[i])
  Z.S <- c(Z.S,qnorm(v))
}
# Plot the results:
hist(Z,25,main="Standard Normal Regular")
hist(Z.S,25,main="Standard Normal Strat")
#c) Exponential

# set number of bins:
k <- 100
# number of draws per bin:
n.A <- 5
# Total number of draws:
M <- k*n.A
# Regular simulation:
Z <- rexp(M)
# Stratified Sampling:
# 1. Generate bin locations:
a <- seq(0,1,length=k+1)
# 2. Sample for each bin:
Z.S <- c(NA)
for (i in 1:k){
  u <- runif(n.A)
  v <- a[i] + u*(a[i+1] - a[i])
  Z.S <- c(Z.S,qexp(v))
}
# Plot the results:
hist(Z,25,main="Exponential Regular")
hist(Z.S,25,main="Exponential Strat")