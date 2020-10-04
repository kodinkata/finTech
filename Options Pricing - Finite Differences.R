rm(list=ls())
# install.packages("plot3D")
library(plot3D)

#Option 1: European call with S0 (Current Stock Price) = 50, K (Strike Price) = 52, r (risk free rate) = 0.1, T (time to expiration)= 5/12, and sigma (stock's volatility) = 0.4.  
#Option 2: European put with S0 = 50, K = 52, r = 0.1, T = 5/12, and sigma = 0.4.  

#Step 1: Using Explicit (Forward) finite difference to price both options

#Step 2: Using implicit (Backward) finite difference to price both options

#Step 3: Using Crank Nicolson (Central) finite difference to price both options

#Step 4: Using Explicit (Forward) finite difference to price American Put Option

# American put with S0 = 50, K = 52, r = 0.1, T = 5/12, and sigma = 0.4.  

#~~~~~~~~~~Step 1~~~~~~~~~~
S0 <- 50
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4
M <- 101
N <- ceiling(t*M^2* sigma^2)

C <- matrix(NA,ncol=N,nrow=M)
Smin <- 0
Smax <-100
Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

discount <- exp(-r*rev(tvec))

h <- function(S,call=T){
  if (call==1){
  h <- pmax(S-K,0)
  }
  else {
    h <- pmax(K-S,0)
  }
}

C[1,] <- pmax(Smin - discount*K,0)
C[M,] <- pmax(Smax - discount*K,0)
C[,N] <- h(Svec, 1)

ivec <- seq(1,M-2)
ivec2 <- ivec*ivec

a <- dt/2*(sigma^2*ivec2 + r*ivec)
b <- 1-dt*(sigma^2*ivec2 +r)
c <- dt/2*(sigma^2*ivec2 - r*ivec)

#Solve in matrix 
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa +Ab + Ac
  return(A)
}

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)
# Solve each column
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j+1]*c[1]
  Kvec[M-2,1] <- C[M,j+1] * a[M-2]
  C[c(2:(M-1)),j] <- A%*%C[c(2:(M-1)),j+1] + Kvec
}

#Use linear interpolation 
price <- approx(x=Svec, y=C[,1],S0,method= "linear",
                rule=1, f=0, ties = mean)
message(sprintf("Call option price: $%5.4f",price$y))
# Plot the surface
plotting <- function(C){
  tmat <- matrix(rep(matrix(tvec,nrow=1),M),nrow=M,byrow=TRUE)
  Smat <- matrix(rep(matrix(Svec,ncol=1),N),ncol=N,byrow=FALSE)
  library(plot3D)
  surf3D(x=tmat,y=Smat,z=C,
         bty="b2",
         main="Call Option Surface",
         xlab="Time",
         ylab="Stock Price",
         zlab="Option Price")
  
}

plotting(C)

# Put Option
C[1,] <- pmax(discount*K-Smin ,0)
C[M,] <- pmax(discount*K-Smax ,0)
C[,N] <- h(Svec, 0)

ivec <- seq(2,M-1)-1
ivec2 <- ivec*ivec
a <- dt/2*(sigma^2*ivec2 + r*ivec)
b <- 1-dt*(sigma^2*ivec2 +r)
c <- dt/2*(sigma^2*ivec2 - r*ivec)

#Solve in matrix 
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa +Ab + Ac
  return(A)
}

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)
# Solve each column
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j+1]*c[1]
  Kvec[M-2,1] <- C[M,j+1] * a[M-2]
  C[c(2:(M-1)),j] <- A%*%C[c(2:(M-1)),j+1] + Kvec
}

#Use linear interpolation 
price <- approx(x=Svec, y=C[,1],S0,method= "linear",
                rule=1, f=0, ties = mean)
message(sprintf("Put option price: $%5.4f",price$y))
plotting(C)



#~~~~~~~~~~Step 2~~~~~~~~~~
# Put Option
C[1,] <- pmax(discount*K-Smin ,0)
C[M,] <- pmax(discount*K-Smax ,0)
C[,N] <- h(Svec, 0)

ivec <- seq(2,M-1)-1
ivec2 <- ivec*ivec
a <- -dt/2*(sigma^2*ivec2 + r*ivec)
b <- 1+dt*(sigma^2*ivec2 +r)
c <- -dt/2*(sigma^2*ivec2 - r*ivec)

#Solve in matrix 
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa +Ab + Ac
  return(A)
}

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j]*c[1]
  Kvec[M-2,1] <- C[M,j]*a[M-2]
  bvec <- C[c(2:(M-1)),j+1] - Kvec
  C[c(2:(M-1)),j] <- solve(A,bvec)
  }
price <- approx(x=Svec, y=C[,1],S0,method= "linear",
                rule=1, f=0, ties = mean)
message(sprintf("Put option price: $%5.4f",price$y))
plotting(C)


# Call Option
C[1,] <- pmax(Smin - discount*K,0)
C[M,] <- pmax(Smax - discount*K,0)
C[,N] <- h(Svec, 1)

ivec <- seq(2,M-1)-1
ivec2 <- ivec*ivec
a <- -dt/2*(sigma^2*ivec2 + r*ivec)
b <- 1+dt*(sigma^2*ivec2 +r)
c <- -dt/2*(sigma^2*ivec2 - r*ivec)

#Solve in matrix 
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa +Ab + Ac
  return(A)
}

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j]*c[1]
  Kvec[M-2,1] <- C[M,j]*a[M-2]
  bvec <- C[c(2:(M-1)),j+1] - Kvec
  C[c(2:(M-1)),j] <- solve(A,bvec)
}
price <- approx(x=Svec, y=C[,1],S0,method= "linear",
                rule=1, f=0, ties = mean)
message(sprintf("Call option price: $%5.4f",price$y))
plotting(C)


#~~~~~~~~~~Step 3~~~~~~~~~~
#Call Option
C[1,] <- pmax(Smin - discount*K,0)
C[M,] <- pmax(Smax - discount*K,0)
C[,N] <- h(Svec, 1)

ivec <- seq(2,M-1)-1
ivec2 <- ivec*ivec

a <- .25*dt*(sigma^2*ivec2 + r*ivec)
b <-  .5*dt*(sigma^2*ivec2 +r)
c <- .25*dt*(sigma^2*ivec2 - r*ivec)

#Solve in matrix 
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa +Ab + Ac
  return(A)
}

upper <- -diag(a[c(1:(M-3))])
lower <- -diag(c[c(2:(M-2))])
main <- 1+b
AL <- triD(main,upper,lower)

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- 1-b
AR <- triD(main,upper,lower)

# Solve each column
for (j in (N-1):1){
  Kvec1 <- matrix(0,ncol=1,nrow=M-2)
  Kvec1[1,1] <- C[1,j]*c[1]
  Kvec1[M-2,1] <- C[M,j] * a[M-2]
  
  Kvec2 <- matrix(0,ncol=1,nrow=M-2)
  Kvec2[1,1] <- C[1,j+1]*c[1]
  Kvec2[M-2,1] <- C[M,j+1]*a[M-2]
  
  bvec <- AR%*%C[c(2:(M-1)),j+1] + Kvec1+Kvec2
  C[c(2:(M-1)),j] <- solve(AL,bvec)
}

#Use linear interpolation 
price <- approx(x=Svec, y=C[,1],S0,method= "linear",
                rule=1, f=0, ties = mean)
message(sprintf("Call option price: $%5.4f",price$y))
plotting(C)

#Put Option
C[1,] <- pmax(discount*K - Smin ,0)
C[M,] <- pmax(discount*K-Smax ,0)
C[,N] <- h(Svec, 0)

ivec <- seq(2,M-1)-1
ivec2 <- ivec*ivec

a <- .25*dt*(sigma^2*ivec2 + r*ivec)
b <-  .5*dt*(sigma^2*ivec2 +r)
c <- .25*dt*(sigma^2*ivec2 - r*ivec)

#Solve in matrix 
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa +Ab + Ac
  return(A)
}

upper <- -diag(a[c(1:(M-3))])
lower <- -diag(c[c(2:(M-2))])
main <- 1+b
AL <- triD(main,upper,lower)

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- 1-b
AR <- triD(main,upper,lower)

# Solve each column
for (j in (N-1):1){
  Kvec1 <- matrix(0,ncol=1,nrow=M-2)
  Kvec1[1,1] <- C[1,j]*c[1]
  Kvec1[M-2,1] <- C[M,j] * a[M-2]
  
  Kvec2 <- matrix(0,ncol=1,nrow=M-2)
  Kvec2[1,1] <- C[1,j+1]*c[1]
  Kvec2[M-2,1] <- C[M,j+1]*a[M-2]
  
  bvec <- AR%*%C[c(2:(M-1)),j+1] + Kvec1+Kvec2
  C[c(2:(M-1)),j] <- solve(AL,bvec)
}

#Use linear interpolation 
price <- approx(x=Svec, y=C[,1],S0,method= "linear",
                rule=1, f=0, ties = mean)
message(sprintf("Call option price: $%5.4f",price$y))
plotting(C)



#~~~~Step 4~~~~~
rm(list=ls())
S0 <- 50
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4
M <- 101
N <- ceiling(t*M^2* sigma^2)

C <- matrix(NA,ncol=N,nrow=M)
Smin <- 0
Smax <-100
Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

discount <- exp(-r*rev(tvec))

h <- function(S,call=T){
  if (call==1){
    h <- pmax(S-K,0)
  }
  else {
    h <- pmax(K-S,0)
  }
}
C[1,] <- pmax(discount*K - Smin,0)
C[M,] <- pmax(discount*K - Smax ,0)
C[,N] <- h(Svec, 0)

ivec <- seq(2,M-1)-1
ivec2 <- ivec*ivec
a <- dt/2*(sigma^2*ivec2 + r*ivec)
b <- 1-dt*(sigma^2*ivec2 +r)
c <- dt/2*(sigma^2*ivec2 - r*ivec)

#Solve in matrix
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa +Ab + Ac
  return(A)
}
# A <- matrix(0, nrow=M-2,ncol=M-2)
# for(i in 1:(M-2)){
#   A[i,i] <- b[i]
# }
# for (i in 1:(M-3)){
#   A[i+1,i] <- c[i+1]
#   A[i,i+1] <- a[i]
# }

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)
# Solve each columnUS option

for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j+1]*c[1]
  Kvec[M-2,1] <- C[M,j+1] * a[M-2]
  intrin <- h(Svec[2:(M-1)],0)
  C[c(2:(M-1)),j] <- pmax(A%*%C[c(2:(M-1)),j+1] + Kvec,intrin)
}

#Use linear interpolation 
price <- approx(x=Svec, y=C[,1],S0,method= "linear",
                rule=1, f=0, ties = mean)
message(sprintf("Put option price: $%5.4f",price$y))
plotting(C)


