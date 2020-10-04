rm(list=ls())


# Step 1: Using Finite Differences to price Barrier Options:
  #a.	Option 1: Down-and-out put barrier option on a stock with S0 = 50, K = 52, r = 0.1, T = 5/12, sigma = 0.4, and a barrier of 40
  #b. b.	Option 2: Up-and-out call barrier option on a stock with S0 = 50, K = 52, r = 0.1, T = 5/12, sigma = 0.4, and a barrier of 60

# Step 2: Using Theta method to price  European options and compare it to the results from finite differences

# Step 3: Compare errors for the 3 methods
#~~~~~Step 1~~~~~~~~
#Explicit
S0 <- 50
K <- 52
BarA <- 40
r <- 0.1
t <- 5/12
sigma <- 0.4
Smin <- BarA
Smax <- 100
M <- 101
N <- ceiling(t*sigma^2*Smax^2/((Smax-BarA)/M)^2)

C <- matrix(NA,ncol=N,nrow=M)

Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

# Generate dS and dt:
dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

# Generate doiscount factors:
discount <- exp(-r*rev(tvec))

# Payoff function:
h <- function(S){
  h <- pmax(K - S,0)
}




# Set the boundary conditions:
C[,N] <- h(Svec)
C[1,] <- 0
C[M,] <- pmax(discount*K - Smax,0)


# 1. Function to get a tri-dagonal matrix of a', b's, and c's
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa + Ab + Ac
  return(A)
}

# Geneate ai, bi, ci:
ivec <- seq(1,M-2) + BarA/dS
ivec2 <- ivec*ivec
a <- dt/2*(sigma^2*ivec2 + r*ivec)
b <- 1 - dt*(sigma^2*ivec2 + r)
c <- dt/2*(sigma^2*ivec2 - r*ivec)

# Solve in matrix form:

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)

# 2. SOlve each column:
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j+1]*c[1]
  Kvec[M-2,1] <- C[M,j+1]*a[M-2]
  intrin <- h(Svec[2:(M-1)])*0
  C[c(2:(M-1)),j] <- pmax((A%*%C[c(2:(M-1)),j+1] + Kvec),intrin)
}

# Use linear interpolation to find the price on the grid:
price <- approx(x=Svec, y=C[,1], S0, method = "linear", 
                rule = 1, f = 0, ties = mean)
message(sprintf("a) Explicit Method: Down-and-Out Put Barrier: $%5.4f ",price$y))

#~~~~~~~~
rm(list=ls())

#Implicit
S0 <- 50
K <- 52
BarA <- 40
r <- 0.1
t <- 5/12
sigma <- 0.4
Smin <- BarA
Smax <- 100
M <- 101
N <- ceiling(t*sigma^2*Smax^2/((Smax-BarA)/M)^2)

C <- matrix(NA,ncol=N,nrow=M)

# Generate prices and times for the grid:
Smin <- BarA
Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

# Generate dS and dt:
dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

# Generate doiscount factors:
discount <- exp(-r*rev(tvec))

# Payoff function:
h <- function(S){
  h <- pmax(K - S,0)
}

# Set the boundary conditions:
C[,N] <- h(Svec)
C[1,] <- 0
C[M,] <- pmax(discount*K - Smax,0)


# 1. Function to get a tri-dagonal matrix of a', b's, and c's
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa + Ab + Ac
  return(A)
}


# Geneate ai, bi, ci:
ivec <- seq(1,M-2) + BarA/dS
ivec2 <- ivec*ivec
a <- -0.5*dt*(sigma*sigma*ivec2 + r*ivec)
b <- 1 + dt*(sigma*sigma*ivec2 + r)
c <- -0.5*dt*(sigma*sigma*ivec2 - r*ivec)

# Solve in matrix form:

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)

# 2. SOlve each column:
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j]*c[1]
  Kvec[M-2,1] <- C[M,j]*a[M-2]
  bvec <- C[c(2:(M-1)),j+1] - Kvec
  C[c(2:(M-1)),j] <- solve(A,bvec)
}

# Use linear interpolation to find the price on the grid:
price <- approx(x=Svec, y=C[,1], S0, method = "linear", 
                rule = 1, f = 0, ties = mean)
message(sprintf("a) Implicit Method: Down-and-Out Put Barrier: $%5.4f ",price$y))

#~~~~~~~~~~~~~~~~
#Crank-Nicolson
rm(list=ls())
S0 <- 50
K <- 52
BarA <- 40
r <- 0.1
t <- 5/12
sigma <- 0.4
Smin <- BarA
Smax <- 100
M <- 101
N <- ceiling(t*sigma^2*Smax^2/((Smax-BarA)/M)^2)

C <- matrix(NA,ncol=N,nrow=M)

# Generate prices and times for the grid:
Smin <- BarA
Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

# Generate dS and dt:
dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

# Generate doiscount factors:
discount <- exp(-r*rev(tvec))

# Payoff function:
h <- function(S){
  h <- pmax(K - S,0)
}

# Set the boundary conditions:
C[,N] <- h(Svec)
C[1,] <- 0
C[M,] <- pmax(discount*K - Smax,0)


# 1. Function to get a tri-dagonal matrix of a', b's, and c's
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa + Ab + Ac
  return(A)
}



# Geneate ai, bi, ci:
ivec <- seq(1,M-2) + BarA/dS
ivec2 <- ivec*ivec
a <- 0.25*dt*(sigma*sigma*ivec2 + r*ivec)
b <- 0.50*dt*(sigma*sigma*ivec2 + r)
c <- 0.25*dt*(sigma*sigma*ivec2 - r*ivec)

# Solve in matrix form:

upper <- -diag(a[c(1:(M-3))])
lower <- -diag(c[c(2:(M-2))])
main <- 1+b
AL <- triD(main,upper,lower)

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- 1-b
AR <- triD(main,upper,lower)

# 2. SOlve each column:
for (j in (N-1):1){
  Kvec1 <- matrix(0,ncol=1,nrow=M-2)
  Kvec1[1,1] <- C[1,j]*c[1]
  Kvec1[M-2,1] <- C[M,j]*a[M-2]
  
  Kvec2 <- matrix(0,ncol=1,nrow=M-2)
  Kvec2[1,1] <- C[1,j+1]*c[1]
  Kvec2[M-2,1] <- C[M,j+1]*a[M-2]
  
  bvec <- AR%*%C[c(2:(M-1)),j+1] + Kvec1 + Kvec2
  C[c(2:(M-1)),j] <- solve(AL,bvec)
}

# Use linear interpolation to find the price on the grid:
price <- approx(x=Svec, y=C[,1], S0, method = "linear", 
                rule = 1, f = 0, ties = mean)
message(sprintf("a) Crank Nicolson: Down-and-Out Put Barrier: $%5.4f ",price$y))

#~~~~~~~~~~~~
#Up-and-out Call
#Explicit
rm(list=ls())
S0 <- 50
BarB <- 60
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4

# Grid size:
M <- 101
N <- ceiling(t*M^2*sigma^2)


C <- matrix(NA,ncol=N,nrow=M)

# Generate prices and times for the grid:
Smin <- 0
Smax <- BarB
Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

# Generate dS and dt:
dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

# Generate doiscount factors:
discount <- exp(-r*rev(tvec))

# Payoff function:
h <- function(S){
  h <- pmax(S-K,0)
}

# Set the boundary conditions: (different from a usual call)
C[1,] <- pmax(Smin - discount*K,0)
C[M,] <- 0 # option value = $0 at the barrier
C[,N] <- h(Svec)

# 1. Function to get a tri-dagonal matrix of a', b's, and c's
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa + Ab + Ac
  return(A)
}
ivec <- seq(2,M-1) - 1
ivec2 <- ivec*ivec
a <- dt/2*(sigma*sigma*ivec2 + r*ivec)
b <- 1 - dt*(sigma*sigma*ivec2 + r)
c <- dt/2*(sigma*sigma*ivec2 - r*ivec)

# Solve in matrix form:

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)

# 2. SOlve each column:
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j+1]*c[1]
  Kvec[M-2,1] <- C[M,j+1]*a[M-2]
  intrin <- h(Svec[2:(M-1)])*0
  C[c(2:(M-1)),j] <- pmax((A%*%C[c(2:(M-1)),j+1] + Kvec),intrin)
}

price <- approx(x=Svec, y=C[,1], S0, method = "linear", 
                rule = 1, f = 0, ties = mean)
message(sprintf("Up-and-out Call Barrier: $%5.4f",price$y))

#Implicit
rm(list=ls())
S0 <- 50
BarB <- 60
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4

# Grid size:
M <- 101
N <- ceiling(t*M^2*sigma^2)


C <- matrix(NA,ncol=N,nrow=M)

# Generate prices and times for the grid:
Smin <- 0
Smax <- BarB
Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

# Generate dS and dt:
dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

# Generate doiscount factors:
discount <- exp(-r*rev(tvec))

# Payoff function:
h <- function(S){
  h <- pmax(S-K,0)
}

# Set the boundary conditions: (different from a usual call)
C[1,] <- pmax(Smin - discount*K,0)
C[M,] <- 0 # option value = $0 at the barrier
C[,N] <- h(Svec)

# 1. Function to get a tri-dagonal matrix of a', b's, and c's
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa + Ab + Ac
  return(A)
}
ivec <- seq(2,M-1) - 1
ivec2 <- ivec*ivec
a <- -0.5*dt*(sigma*sigma*ivec2 + r*ivec)
b <- 1 + dt*(sigma*sigma*ivec2 + r)
c <- -0.5*dt*(sigma*sigma*ivec2 - r*ivec)

# Solve in matrix form:

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- b
A <- triD(main,upper,lower)

# 2. SOlve each column:
for (j in (N-1):1){
  Kvec <- matrix(0,ncol=1,nrow=M-2)
  Kvec[1,1] <- C[1,j]*c[1]
  Kvec[M-2,1] <- C[M,j]*a[M-2]
  bvec <- C[c(2:(M-1)),j+1] - Kvec
  C[c(2:(M-1)),j] <- solve(A,bvec)
}

price <- approx(x=Svec, y=C[,1], S0, method = "linear", 
                rule = 1, f = 0, ties = mean)
message(sprintf("Up-and-out Call Barrier: $%5.4f",price$y))


#Crank-Nicolson
rm(list=ls())
S0 <- 50
BarB <- 60
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4

# Grid size:
M <- 101
N <- ceiling(t*M^2*sigma^2)


C <- matrix(NA,ncol=N,nrow=M)

# Generate prices and times for the grid:
Smin <- 0
Smax <- BarB
Svec <- seq(Smin,Smax,length=M)
tvec <- seq(0,t,length=N)

# Generate dS and dt:
dt <- t/(N-1)
dS <- Svec[2] - Svec[1]

# Generate doiscount factors:
discount <- exp(-r*rev(tvec))

# Payoff function:
h <- function(S){
  h <- pmax(S-K,0)
}

# Set the boundary conditions: (different from a usual call)
C[1,] <- pmax(Smin - discount*K,0)
C[M,] <- 0 # option value = $0 at the barrier
C[,N] <- h(Svec)

# 1. Function to get a tri-dagonal matrix of a', b's, and c's
triD <- function(main,upper,lower){
  m <- length(main)
  zeroCol <- matrix(0,ncol=1,nrow=m-1)
  zeroRow <- matrix(0,ncol=m,nrow=1)
  Aa <- rbind(cbind(zeroCol,upper),zeroRow)
  Ab <- diag(main)
  Ac <- rbind(zeroRow,cbind(lower,zeroCol))
  A <- Aa + Ab + Ac
  return(A)
}
ivec <- seq(2,M-1) -1
ivec2 <- ivec*ivec
a <- 0.25*dt*(sigma^2*ivec2 + r*ivec)
b <- 0.5 *dt*(sigma^2*ivec2 + r)
c <- 0.25*dt*(sigma^2*ivec2 - r*ivec)

# Solve in matrix form:
upper <- -diag(a[c(1:(M-3))])
lower <- -diag(c[c(2:(M-2))])
main <- 1+b
AL <- triD(main,upper,lower)

upper <- diag(a[c(1:(M-3))])
lower <- diag(c[c(2:(M-2))])
main <- 1-b
AR <- triD(main,upper,lower)

# 2. SOlve each column:
for (j in (N-1):1){
  Kvec1 <- matrix(0,ncol=1,nrow=M-2)
  Kvec1[1,1] <- C[1,j]*c[1]
  Kvec1[M-2,1] <- C[M,j]*a[M-2]
  
  Kvec2 <- matrix(0,ncol=1,nrow=M-2)
  Kvec2[1,1] <- C[1,j]*c[1]
  Kvec2[M-2,1] <- C[M,j]*a[M-2]
  
  bvec <- AR%*%C[c(2:(M-1)),j+1] + Kvec1 + Kvec2
  C[c(2:(M-1)),j] <- solve(AL,bvec)
}

price <- approx(x=Svec, y=C[,1], S0, method = "linear", 
                rule = 1, f = 0, ties = mean)
message(sprintf("Up-and-out Call Barrier: $%5.4f",price$y))

#~~~~~~~~~Step 2~~~~~~~~
rm(list=ls())
#Explicit Call
S0 <- 50
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4
Smin <- 0
Smax <- 100
M <- 101
N <- ceiling(t*M^2*sigma^2)
C<-matrix(ncol=N,nrow=M)
Svec<-seq(Smin,Smax,length=M)
tvec<-seq(0,t,length=N)
dt <- t/(N-1)
ds<-Svec[2]-Svec[1]

disc<-exp(-r*rev(tvec))

thetafun <- function(theta,call){
  if(call==1){
  h<-function(S){
    h<-pmax(S-K,0)
  }
  
  C[1,] <- pmax(Smin - disc*K,0) 
  C[M,] <- pmax(Smax - disc*K,0)
  C[,N] <- h(Svec)   
  }
  else{
    h<-function(S){
      h<-pmax(K-S,0)
    }
    
    C[1,] <- pmax(disc*K - Smin,0) 
    C[M,] <- pmax(disc*K - Smax,0)
    C[,N] <- h(Svec) 
  }
  ivec<-seq(2,M-1)-1
  ivec2<-ivec*ivec
  
  a <- dt/2*(sigma^2*ivec2 + r*ivec) 
  b <- dt*(sigma^2*ivec2 + r)
  c <- dt/2*(sigma^2*ivec2 - r*ivec)
  
  triD <- function(main,upper,lower){
    m <- length(main)
    zeroCol <- matrix(0,ncol=1,nrow=m-1) 
    zeroRow <- matrix(0,ncol=m,nrow=1)
    Aa <- rbind(cbind(zeroCol,upper),zeroRow) 
    Ab <- diag(main)
    Ac <- rbind(zeroRow,cbind(lower,zeroCol)) 
    A <- Aa + Ab + Ac
    return(A)
  }
  upper <- diag(a[c(1:(M-3))])*theta
  lower <- diag(c[c(2:(M-2))])*theta
  main <- -1-b*theta
  AL <- triD(main,upper,lower)
  
  upper <- diag(a[c(1:(M-3))])*(theta-1)
  lower <- diag(c[c(2:(M-2))])*(theta-1)
  main <- -1-b*(theta-1)
  AR <- triD(main,upper,lower)
  
  for (j in (N-1):1){
    Kvec <- matrix(0,ncol=1,nrow=M-2)
    Kvec[1,1] <- C[1,j+1]*c[1]*theta
    Kvec[M-2,1] <- C[M,j+1]*a[M-2]*theta
    
    Kvec2 <- matrix(0,ncol=1,nrow=M-2)
    
    Kvec2[1,1] <- C[1,j+1]*c[1]*(theta-1)
    Kvec2[M-2,1] <- C[M,j+1]*a[M-2]*(theta-1)
    bvec <- AR%*%C[c(2:(M-1)),j+1] - Kvec + Kvec2
    C[c(2:(M-1)),j] <- solve(AL,bvec)
  }
  price <- approx(x=Svec, y=C[,1], S0, method = "linear", rule = 1, f = 0, ties = mean)
}
thetaCall0<-thetafun(theta=0,call=1)
message(sprintf("Call option price: $%5.4f",thetaCall0$y))
thetaCall0.5<-thetafun(theta=0.5,call=1)
message(sprintf("Call option price: $%5.4f",thetaCall0.5$y))
thetaCall1<-thetafun(theta=1,call=1)
message(sprintf("Call option price: $%5.4f",thetaCall1$y))

thetaPut0<-thetafun(theta=0,call=0)
message(sprintf("Put option price: $%5.4f",thetaPut0$y))
thetaPut0.5<-thetafun(theta=0.5,call=0)
message(sprintf("Put option price: $%5.4f",thetaPut0.5$y))
thetaPut1<-thetafun(theta=1,call=0)
message(sprintf("Put option price: $%5.4f",thetaPut1$y))


#~~~~~~~~Step 3~~~~~~~~~~~
rm(list=ls())
#Explicit Call
S0 <- 50
K <- 52
r <- 0.1
t <- 5/12
sigma <- 0.4
Smin <- 0
Smax <- 100
M <- 101
N <- ceiling(t*M^2*sigma^2)
explicit <- function(S0,K,r,t,sigma,Smin,Smax,M,N,call){
  Svec <- seq(Smin,Smax,length=M)
  tvec <- seq(0,t,length=N)
  
  # Generate dS and dt:
  dt <- t/(N-1)
  dS <- Svec[2] - Svec[1]
  
  # Generate doiscount factors:
  discount <- exp(-r*rev(tvec))
  
  # Generate and empty matrix:
  P <- matrix(NA,ncol=N,nrow=M)
  if(call==1){
    hp <- function(S){
      hp <- pmax(S-K,0)
    }
    # Set the boundary conditions:
    P[1,] <- pmax(Smin - discount*K,0)
    P[M,] <- pmax(Smax - discount*K,0)
    P[,N] <- hp(Svec)
    }else{
    hp <- function(S){
      hp <- pmax(K-S,0)
    }
    # Set the boundary conditions:
    P[1,] <- pmax(discount*K - Smin,0)
    P[M,] <- pmax(discount*K - Smax,0)
    P[,N] <- hp(Svec)
    }
  
  # Geneate ai, bi, ci:
  ivec <- seq(1,M-2) 
  ivec2 <- ivec*ivec
  
  N <- max(N,ceiling(t*M^2*sigma^2)) # ensure convergence
  
  a <- dt/2*(sigma^2*ivec2 + r*ivec)
  b <- 1 - dt*(sigma^2*ivec2 + r)
  c <- dt/2*(sigma^2*ivec2 - r*ivec)
  
  # Solve in matrix form:
  # 1. Get a tri-dagonal matrix of a', b's, and c's
  triD <- function(main,upper,lower){
    m <- length(main)
    zeroCol <- matrix(0,ncol=1,nrow=m-1)
    zeroRow <- matrix(0,ncol=m,nrow=1)
    Aa <- rbind(cbind(zeroCol,upper),zeroRow)
    Ab <- diag(main)
    Ac <- rbind(zeroRow,cbind(lower,zeroCol))
    A <- Aa + Ab + Ac
    return(A)
  }
  upper <- diag(a[c(1:(M-3))])
  lower <- diag(c[c(2:(M-2))])
  main <- b
  A <- triD(main,upper,lower)
  
  # 2. SOlve each column:
  for (j in (N-1):1){
    Kvec <- matrix(0,ncol=1,nrow=M-2)
    Kvec[1,1] <- P[1,j+1]*c[1]
    Kvec[M-2,1] <- P[M,j+1]*a[M-2]
    P[c(2:(M-1)),j] <- A%*%P[c(2:(M-1)),j+1] + Kvec
  }
  # Use linear interpolation to find the price on the grid:
  price <- approx(x=Svec, y=P[,1], S0, method = "linear", 
                  rule = 1, f = 0, ties = mean)
  # message(sprintf("Put option price: $%5.4f",price$y))
  
  
  tmat <- matrix(rep(matrix(tvec,nrow=1),M),nrow=M,byrow=TRUE)
  Smat <- matrix(rep(matrix(Svec,ncol=1),N),ncol=N,byrow=FALSE)
  
  results <- list("tmat"=tmat,"Smat"=Smat,"f"=P,"f0"=price$y)
  return (results)
}

implicit <- function(S0,K,r,t,sigma,Smin,Smax,M,N,call){
  Svec <- seq(Smin,Smax,length=M)
  tvec <- seq(0,t,length=N)
  
  # Generate dS and dt:
  dt <- t/(N-1)
  dS <- Svec[2] - Svec[1]
  
  # Generate doiscount factors:
  discount <- exp(-r*rev(tvec))
  
  # Generate and empty matrix:
  P <- matrix(NA,ncol=N,nrow=M)
  if(call==1){
    hp <- function(S){
      hp <- pmax(S-K,0)
    }
    # Set the boundary conditions:
    P[1,] <- pmax(Smin - discount*K,0)
    P[M,] <- pmax(Smax - discount*K,0)
    P[,N] <- hp(Svec)
  }else{
    hp <- function(S){
      hp <- pmax(K-S,0)
    }
    # Set the boundary conditions:
    P[1,] <- pmax(discount*K - Smin,0)
    P[M,] <- pmax(discount*K - Smax,0)
    P[,N] <- hp(Svec)
  }
  
  # Geneate ai, bi, ci:
  ivec <- seq(1,M-2) 
  ivec2 <- ivec*ivec
  
  N <- ceiling(t*M^2*sigma^2)
  
  a <- -0.5*dt*(sigma*sigma*ivec2 + r*ivec)
  b <- 1 + dt*(sigma*sigma*ivec2 + r)
  c <- -0.5*dt*(sigma*sigma*ivec2 - r*ivec)
  
  # Solve in matrix form:
  # 1. Get a tri-dagonal matrix of a', b's, and c's
  triD <- function(main,upper,lower){
    m <- length(main)
    zeroCol <- matrix(0,ncol=1,nrow=m-1)
    zeroRow <- matrix(0,ncol=m,nrow=1)
    Aa <- rbind(cbind(zeroCol,upper),zeroRow)
    Ab <- diag(main)
    Ac <- rbind(zeroRow,cbind(lower,zeroCol))
    A <- Aa + Ab + Ac
    return(A)
  }
  upper <- diag(a[c(1:(M-3))])
  lower <- diag(c[c(2:(M-2))])
  main <- b
  A <- triD(main,upper,lower)
  
  # 2. SOlve each column:
  for (j in (N-1):1){
    Kvec <- matrix(0,ncol=1,nrow=M-2)
    Kvec[1,1] <- P[1,j]*c[1]
    Kvec[M-2,1] <- P[M,j]*a[M-2]
    bvec <- P[c(2:(M-1)),j+1] - Kvec
    P[c(2:(M-1)),j] <- solve(A,bvec)
  }
  # Use linear interpolation to find the price on the grid:
  price <- approx(x=Svec, y=P[,1], S0, method = "linear", 
                  rule = 1, f = 0, ties = mean)
  # message(sprintf("Put option price: $%5.4f",price$y))
  
  
  tmat <- matrix(rep(matrix(tvec,nrow=1),M),nrow=M,byrow=TRUE)
  Smat <- matrix(rep(matrix(Svec,ncol=1),N),ncol=N,byrow=FALSE)
  
  results <- list("tmat"=tmat,"Smat"=Smat,"f"=P,"f0"=price$y)
  return (results)
}
crank_nic <- function(S0,K,r,t,sigma,Smin,Smax,M,N,call){
  Svec <- seq(Smin,Smax,length=M)
  tvec <- seq(0,t,length=N)
  
  # Generate dS and dt:
  dt <- t/(N-1)
  dS <- Svec[2] - Svec[1]
  
  # Generate doiscount factors:
  discount <- exp(-r*rev(tvec))
  
  # Generate and empty matrix:
  P <- matrix(NA,ncol=N,nrow=M)
  if(call==1){
    hp <- function(S){
      hp <- pmax(S-K,0)
    }
    # Set the boundary conditions:
    P[1,] <- pmax(Smin - discount*K,0)
    P[M,] <- pmax(Smax - discount*K,0)
    P[,N] <- hp(Svec)
  }else{
    hp <- function(S){
      hp <- pmax(K-S,0)
    }
    # Set the boundary conditions:
    P[1,] <- pmax(discount*K - Smin,0)
    P[M,] <- pmax(discount*K - Smax,0)
    P[,N] <- hp(Svec)
  }
  
  # Geneate ai, bi, ci:
  ivec <- seq(1,M-2) 
  ivec2 <- ivec*ivec
  
  N <- ceiling(t*M^2*sigma^2)
  
  a <- 0.25*dt*(sigma*sigma*ivec2 + r*ivec)
  b <- 0.50*dt*(sigma*sigma*ivec2 + r)
  c <- 0.25*dt*(sigma*sigma*ivec2 - r*ivec)
  
  # Solve in matrix form:
  # 1. Get a tri-dagonal matrix of a', b's, and c's
  triD <- function(main,upper,lower){
    m <- length(main)
    zeroCol <- matrix(0,ncol=1,nrow=m-1)
    zeroRow <- matrix(0,ncol=m,nrow=1)
    Aa <- rbind(cbind(zeroCol,upper),zeroRow)
    Ab <- diag(main)
    Ac <- rbind(zeroRow,cbind(lower,zeroCol))
    A <- Aa + Ab + Ac
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
  
  # 2. SOlve each column:
  for (j in (N-1):1){
    Kvec1 <- matrix(0,ncol=1,nrow=M-2)
    Kvec1[1,1] <- P[1,j]*c[1]
    Kvec1[M-2,1] <- P[M,j]*a[M-2]
    
    Kvec2 <- matrix(0,ncol=1,nrow=M-2)
    Kvec2[1,1] <- P[1,j+1]*c[1]
    Kvec2[M-2,1] <- P[M,j+1]*a[M-2]
    
    bvec <- AR%*%P[c(2:(M-1)),j+1] + Kvec1 + Kvec2
    P[c(2:(M-1)),j] <- solve(AL,bvec)
  }
  # Use linear interpolation to find the price on the grid:
  price <- approx(x=Svec, y=P[,1], S0, method = "linear", 
                  rule = 1, f = 0, ties = mean)
  # message(sprintf("Put option price: $%5.4f",price$y))
  
  
  tmat <- matrix(rep(matrix(tvec,nrow=1),M),nrow=M,byrow=TRUE)
  Smat <- matrix(rep(matrix(Svec,ncol=1),N),ncol=N,byrow=FALSE)
  
  results <- list("tmat"=tmat,"Smat"=Smat,"f"=P,"f0"=price$y)
  return (results)
}
soln <- list(NA)
soln[[1]] <- explicit(S0,K,r,t,sigma,Smin,Smax,M,N,call=1)
soln[[2]] <- explicit(S0,K,r,t,sigma,Smin,Smax,M,N,call=0)
soln[[3]] <- implicit(S0,K,r,t,sigma,Smin,Smax,M,N,call=1)
soln[[4]] <- implicit(S0,K,r,t,sigma,Smin,Smax,M,N,call=0)
soln[[5]] <- crank_nic(S0,K,r,t,sigma,Smin,Smax,M,N,call=1)
soln[[6]] <- crank_nic(S0,K,r,t,sigma,Smin,Smax,M,N,call=0)

# Check the results versus Black-Scholes to see where the errors are:
# (Not required for HW 10)

# First compute a matrix of BLS prices:
checkErrors <- function(C,P,Smat,tmat,plot=0){
  myBLScall <- function(Smat,K,r,tmat,sigma){
    mydim <- dim(Smat)
    M <- mydim[1]
    N <- mydim[2]
    t <- tmat[M,N]
    tmat <- t - tmat
    
    c <- matrix(NA, nrow=M, ncol=N)
    
    d1 <- (log(Smat/K) + (r + sigma^2/2)*tmat) / (sigma*sqrt(tmat))
    d2 <- d1 - sigma*sqrt(tmat)
    c <- Smat * pnorm(d1) - K*exp(-r*tmat)*pnorm(d2)
    p <- -Smat * pnorm(-d1) + K*exp(-r*tmat)*pnorm(-d2)
    opts <- list(c=c,p=p)
    return(opts)
  }
  BLS <- myBLScall(Smat,K,r,tmat,sigma)
  
  # Calls:
  check_call <- C - BLS$c
  
  # Puts:
  check_put <- P - BLS$p
  
  results <- list("check_call"=check_call,"check_put"=check_put)
  return(results)
  if (plot == 1){
    pc <- plot_ly(x = tmat, y = Smat, z = check_call) %>%
      add_surface(colorscale="Portland") %>%
      layout(scene = list(xaxis = list(title = 'Time'),
                          yaxis = list(title = 'Stock Price'),
                          zaxis = list(title = 'BLS Call')))
    print(pc)
    
    pp <- plot_ly(x = tmat, y = Smat, z = check_put) %>%
      add_surface(colorscale="Portland") %>%
      layout(scene = list(xaxis = list(title = 'Time'),
                          yaxis = list(title = 'Stock Price'),
                          zaxis = list(title = 'BLS Put')))
    print(pp)
  }
}

Smat <- soln[[1]]$Smat
tmat <- soln[[1]]$tmat
# Explicit Method:
ex <- checkErrors(soln[[1]]$f,soln[[2]]$f,Smat,tmat)

# Implicit Method:
im <- checkErrors(soln[[3]]$f,soln[[4]]$f,Smat,tmat)

# Crank-Nicoloson Method:
cn <- checkErrors(soln[[5]]$f,soln[[6]]$f,Smat,tmat)

library(dplyr)
library(plotly)

# Call Options :
pc <- plot_ly() %>%
  add_surface(x = tmat, y = Smat, z = ex$check_call) %>%
  add_surface(x = tmat, y = Smat, z = .1+im$check_call) %>%
  add_surface(x = tmat, y = Smat, z = .2+cn$check_call) %>%
  layout(scene = list(xaxis = list(title = 'Time'),
                      yaxis = list(title = 'Stock Price'),
                      zaxis = list(title = 'BLS Call')))
print(pc)

# Put Options :
pc <- plot_ly() %>%
  add_surface(x = tmat, y = Smat, z = ex$check_put) %>%
  add_surface(x = tmat, y = Smat, z = .1+im$check_put) %>%
  add_surface(x = tmat, y = Smat, z = .2+cn$check_Put) %>%
  layout(scene = list(xaxis = list(title = 'Time'),
                      yaxis = list(title = 'Stock Price'),
                      zaxis = list(title = 'BLS Put')))
print(pc)


# Compare errors from Explicit and CN Methods:
pc <- plot_ly() %>%
  add_surface(x = tmat, y = Smat, z = ex$check_call - cn$check_call) %>%
  layout(scene = list(xaxis = list(title = 'Time'),
                      yaxis = list(title = 'Stock Price'),
                      zaxis = list(title = 'BLS Call)')))
print(pc)

# Compare errors from Explicit and CN Methods:
pc <- plot_ly() %>%
  add_surface(x = tmat, y = Smat, z = ex$check_put - cn$check_put) %>%
  layout(scene = list(xaxis = list(title = 'Time'),
                      yaxis = list(title = 'Stock Price'),
                      zaxis = list(title = 'BLS Put)')))
print(pc)