#simulate Merton's jump diffusion model for stocks in the real world
#lambda <- 7: arrival of jumps follows a Poisson distribution with lambda = 7 (7 jumps per year)
#Intensity of the jumps: Yj ~ LN(a,b2) with a = 0.05 and b = 0.1.


# Simulate JD model - Fixed Dates:
set.seed(1234)
# Write a function to compare JD and GBM models:
myJD1 <- function(S0,M,N,mu,sigma,a,b,t,lambda){
  # Initialize:
  X.GBM <- X.JD <-matrix(NA,ncol=N+1,nrow=M)
  X.GBM[,1] <- log(S0)
  X.JD[,1] <- log(S0)
  dt <- t/N
  sqdt <- sqrt(dt)
  for (i in 1:N){
    # Draw a common Z for GBM and JD:
    Z <- matrix(rnorm(M),ncol=1)
    # Poiss for JD:
    NN <- matrix(rpois(M,lambda*dt),ncol=1)
    # Assume lognormal jumps:
    Z2 <- matrix(rnorm(M),ncol=1)
    MM <- a*NN + b*sqrt(NN)*Z2
    # Update:
    X.JD[,i+1] <- X.JD[,i] + (mu - 0.5*sigma^2)*dt +
      sigma*sqdt*Z + MM
    X.GBM[,i+1] <- X.GBM[,i] + (mu - 0.5*sigma^2)*dt +
      sigma*sqdt*Z
  }
  S.GBM <- exp(X.GBM)
  S.JD <- exp(X.JD)
  out <- list("GBM"=S.GBM,"JD"=S.JD)
  return(out)
}
S0 <- 50
N <- 252
r <- mu <- 0.1
sigma <- 0.4
a <- .05
b <- 0.1
lambda <- 7
M <- 10000
t <- 1
JD.test <- myJD1(S0,M,N,mu,sigma,a,b,t,lambda)
time <- as.matrix(seq(0,t,length=N+1),ncol=1)
S1 <- t(JD.test$GBM)
S2 <- t(JD.test$JD)
df1 <- as.data.frame(cbind(time,S1))
df2 <- as.data.frame(cbind(time,S2))
# Plot some paths:
mP <- 5 #limit the number of paths being plotted
Anames <- rep("",mP)
for (i in 1:mP){
  Anames[i] <- paste("A",i,sep="")
}
names(df1) <- c("time", Anames)
names(df2) <- c("time", Anames)
# This creates a new data frame with columns x, variable and value
# x is the id, variable holds each of our timeseries designation
library("reshape")
library("ggplot2")
df1.melted <- melt(df1[,c(1:mP)], id = "time")
df2.melted <- melt(df2[,c(1:mP)], id = "time")
p_Common <- ggplot() +
  geom_line(data = df1.melted, aes(x = time, y = value, factor=variable),color="red") +
  geom_line(data = df2.melted, aes(x = time, y = value, factor=variable),color="black") +
  xlab('time') +
  ylab('S')
print(p_Common)
