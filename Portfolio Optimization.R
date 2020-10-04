rm(list=ls())
#Mean-variance portfolio optimization for a set of 10 stocks
#Plotting the efficient frontier
#Step4- Portfolio optimization with a risk free asset

# An easy way to load multiple packages:
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ggplot2", "quantmod", "reshape", "alabama")
ipak(packages)

# ----------------------
# Step 1. Get the stock data
# Define some companies and also download the S&P 500 index 
tics <- c("TNDM","IBM","STZ","WFC","INTC",
          "CSCO","BABA","JPM","AMZN","MS",
          "SPY")

# Get market data over the past 10 years:
P.list <- lapply(tics, function(tic) 
  get(getSymbols(tic, from = "2008-10-01")))

sapply(P.list,nrow)

tail(P.list[[11]],10)

# Get the adjusted prices into a single object
P.adj <- lapply(P.list, function(p) p[,6])

# Merge the elements of the list
P <- Reduce(merge,P.adj)
names(P) <- tics

head(P,10)

tail(P,10)

# Compute returns:
Rets <- log(P/lag(P))

tail(P,5)
tail(Rets,5)

R <- apply(Rets,2,function(x) mean(x,na.rm = TRUE))

# ----------------------
# Step 2. Compute the expected returns:
# For this example, I will run the CAPM to get the mean return vector 


# First get some risk free rate data:
indices <- c("DGS3MO")
ind.list <- lapply(indices, function(tic) 
  get(getSymbols(tic, from = "2008-10-01", src = "FRED")) )
IND <- Reduce(merge,ind.list)
rf <- IND
# Make the interest rates daily:
IND <- log(1 + IND)/252
tail(IND,5)

# Merge the stock and index returns with the Treasury rates
dat <- merge(Rets,IND)

# Keep only the days with non missing data for all stocks and rates:
dat <- dat[complete.cases(dat),]
rf <- IND[length(IND)]

# Assume a 7% market risk premium:
RP <- 0.07
  
# Subtract off the risk free rate from all returns:
n <- length(tics) - 1
dat_rf <- dat
for (i in 1:(n+1)){
  dat_rf[,i] <- dat[,i] - dat[,(n+2)]
}

# Estimate the betas and mean return vector:
#  r_i - r_f = beta_i(r_m - r_f) + e_i
beta <- c(NA)
R <- c(NA)
for (i in 1:(n)){
  beta[i] <- lm(dat_rf[,i] ~ dat_rf[,(n+1)] - 1)$coefficients[1]
  R[i] <- rf + beta[i]*RP
}

# Esitmate the covariance matrix:
Sigma <- var(dat_rf[,1:10], use="pairwise") * 252

# ----------------------------------
# ----------------------
# Step 3. Plot the eficient frontier:
# ****NOTE****: generating the random portfolios to plot
# First define two functions to return the portfolio 
#   mean and variance given a vector of weights:
R_pf <- function(w,R=R){
  r <- t(w) %*% R
  return(r)
}

S_pf <- function(w,Sigma=Sigma){
  s <- t(w) %*% Sigma %*% w
  return(s)
}

# Define a function to generate random portfolios to plot:
randPfs <- function(R,Sigma){
  w <- rep(NA,length=n)
  w[1] <- runif(1)
  for (i in 2:n){
    w[i] <- runif(1,min=0,max=1 - sum(w,na.rm=TRUE))
  }
  w <- w / sum(w)
  w <- sample(w)
  R0 <- R_pf(w,R)
  S0 <- S_pf(w,Sigma)
  results <- list("R"=R0,"S"=S0,"w"=w)
  return(results)
}

# Call the above function and geenrate the random pf's
N <- n*500
pfs <- as.data.frame(matrix(NA,ncol=2+n,nrow=N))
for (i in 1:N){
  tmp <- randPfs(R,Sigma)
  pfs[i,] <- matrix(c(tmp$R,sqrt(tmp$S),t(tmp$w)))
}


# Plot the results of the random portoflios
p <- ggplot(data=pfs,aes(x=V2,y=V1)) + 
  geom_point() +
  ylab(label="R_pf") +
  xlab(label="sigma_pf") +
  ggtitle("Random Portfolios")  + 
  expand_limits(x=0, y=0)

# Plot the 100% portfolios:
for (i in 1:n){
  wtmp <- rep(0,length=n)
  wtmp[i]=1
  p <- p + geom_point(x=c(sqrt(S_pf(wtmp,Sigma))),y=c(R_pf(wtmp,R)),
                      color="green")
}
print(p)


# ------Minimum Variance Portfolio--------
# Inequality
hin <- function(x){
  # h <- x
  
  # # Or, alternatively:
  # b <- rep(0,length=length(x))
  # for (i in 1:j){
  #   h[j] <- x[j] - b[j]
  # }
  
  A <- diag(rep(1,length=length(x)))
  b <- rep(0,length=length(x))
  h <- as.vector(A%*%x - b)
  
  return(h)
}

# Equality
heq <- function(x){
  h <- sum(x) - 1
  return(h)
}

# ------Objective Function: Minimize Variance-----
eval_f <- function(x){
  s <- t(x) %*% Sigma %*% x
  return(s)
}

# Set the starting valae to a feasible portfolio:
x0 <- rep(1/n,length=n)

# Solve the optimization problem:
#   NOTE: the last two "trace" arguments in control.outer 
#   and countrol.optim turn off the output so it doesn't print
#   to the screen.  It is useful for the loops below.
tmp <- constrOptim.nl(par=x0, fn=eval_f,  
                      heq=heq, hin=hin,
                      "control.outer"=list("trace"=FALSE),
                      "control.optim"=list("reltol"=1e-12,
                                           "trace"=1)) 

# Print the return and risk of the estimated minimum variance portfolio:
Rmin <- tmp$par %*% R
print(Rmin)
Smin <- sqrt(S_pf(tmp$par,Sigma))
print(Smin)

message(sprintf("Minimum variance portfolio:\n\tR = %5.4f\n\tSigma=%5.4f",Rmin,Smin))

# ***************************************************
# Solve for the efficient frontier
# Constraints wi >= 0
# Set this in the loop because it changes every time

# sum(w) = 1
heq <- function(x){
  h <- sum(x) - 1
  return(h)
}

# Set the number of portoflios along [Rmin,Rmax]
npf <- 100
soln <- as.data.frame(matrix(NA,ncol=2+n,nrow=npf))
Rconst <- seq(Rmin,max(R),length=npf)

# Objective function (minimize variance)
eval_f <- function(x){
  s <- (t(x) %*% Sigma %*% x)
  return(s)
}

# Gradient of the objective function (if you do not
# define it, a numerical approximation will be used)
eval_g <- function(x){
  g <- 0.5*t(x)%*%Sigma
  return(g)
}

# Loop through to solve:
outMat <- as.data.frame(matrix(NA,ncol=9,nrow=npf))
outWts <- as.data.frame(matrix(NA,ncol=(n+1),nrow=npf))

# Set a starting value to a feasible portfolio 
# (Note that if the feasible region changes throughout the
# loop, then you will need to change x0 in the loop)
x0 <- rep(0,length=n)
x0[which.max(R)] <- 1

for (i in 1:npf){
  hin <- function(x){
    eps <- 1e-6
    h <- rep(NA,1)
    # Define no shorting constraints (the eps is to
    # mitigate the impact of rounding erros)
    for (j in 1:length(x)){
      h[j] <- x[j] + eps
    }
    # Define the return constraint (Note - this is
    # the only element that changes in the loop)
    h[length(x)+1] <- t(x)%*%R - Rconst[i] + eps
    return(h)
  }
  # Solve the problem (the tryCatch() assures the
  # code will run if a single iteration returns an
  # error)
  tryCatch({
    tmp <- constrOptim.nl(par=x0, fn=eval_f, 
                          gr=NULL,
                          heq=heq, hin=hin,
                          "control.outer"=list("trace"=FALSE),
                          "control.optim"=list("reltol"=1e-12,
                                               "abstol"=1e-12,
                                               "trace"=0))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # Save the resulting weight vector
  wtmp <- tmp$par
  
  # Save the results
  soln[i,2] <- sqrt(S_pf(wtmp,Sigma))
  soln[i,1] <- R_pf(wtmp,R)
  
  # ****THE FOLLOWING BLOCK IS SAVES AND PRINTS RESULTS
  soln[i,c(3:(n+2))] <- wtmp
  message(sprintf("%d. %6.5f\t%6.5f\t%6.5f\t%d\t%6.5f\t%8.7f\t%d\t%d",
                  i,Rconst[i],soln[i,1],soln[i,2],
                  tmp$outer.iterations,sum(wtmp),
                  Rconst[i]-soln[i,1],
                  tmp$convergence,
                  tmp$counts[1]))
  outMat[i,] <- c(i,Rconst[i],soln[i,1],soln[i,2],
                  tmp$outer.iterations,sum(wtmp),
                  Rconst[i]-soln[i,1],
                  tmp$convergence,
                  tmp$counts[1])
  outWts[i,] <- c(i,tmp$par)
  # **********END PRINTING BLOCK**********
}

# Add results to the plot:
p <- p + geom_point(data=soln, aes(x=soln[,2],y=soln[,1]),color="blue")
print(p)






# ---------------------------
# ***************************
# ---------------------------
# Step 4: Add in a risk free asset:

# 1. Return Vector
R <- c(R,rf)
tics <- c(tics,"RF")
print(R)

# 2. Covariance Matrix
Sigma <- rbind(Sigma,rep(0,length=n))
Sigma <- cbind(Sigma,rep(0,length=(n+1)))
print(Sigma)

# 3. Constraints
# sum(w) = 1
heq <- function(x){
  h <- sum(x) - 1
  return(h)
}

n<-length(R)
# Loop through to solve:
outMat <- as.data.frame(matrix(NA,ncol=9,nrow=npf))
outWts <- as.data.frame(matrix(NA,ncol=(n+1),nrow=npf))

# Set the number of portoflios along [Rmin,Rmax]
npf <- 100
soln <- as.data.frame(matrix(NA,ncol=2+n,nrow=npf))
Rmin <- 0
Rconst <- seq(Rmin,max(R),length=npf)

# Set a starting value to a feasible portfolio 
# (Note that if the feasible region changes throughout the
# loop, then you will need to change x0 in the loop)
x0 <- rep(0,length=n)
x0[which.max(R)] <- 1

for (i in 1:npf){
  hin <- function(x){
    eps <- 1e-6
    h <- rep(NA,1)
    # Define no shorting constraints (the eps is to
    # mitigate the impact of rounding erros)
    for (j in 1:length(x)){
      h[j] <- x[j] + eps
    }
    # Define the return constraint (Note - this is
    # the only element that changes in the loop)
    h[length(x)+1] <- t(x)%*%R - Rconst[i] + eps
    return(h)
  }
  # Solve the problem (the tryCatch() assures the
  # code will run if a single iteration returns an
  # error)
  tryCatch({
    tmp <- constrOptim.nl(par=x0, fn=eval_f, 
                          gr=NULL,
                          heq=heq, hin=hin,
                          "control.outer"=list("trace"=FALSE),
                          "control.optim"=list("reltol"=1e-12,
                                               "abstol"=1e-12,
                                               "trace"=0))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # Save the resulting weight vector
  wtmp <- tmp$par
  
  # Save the results
  soln[i,2] <- sqrt(S_pf(wtmp,Sigma))
  soln[i,1] <- R_pf(wtmp,R)
  
}

# Add results to the plot:
p <- p + geom_point(data=soln, aes(x=soln[,2],y=soln[,1]),color="yellow")
print(p)


# ---------------------------
# ***************************
# ---------------------------
# Step 5: Maximize Expected Utility: (with risk-free asset)
lambda <- 1.5
# Inequality
hin <- function(x){
  # h <- x
  
  # # Or, alternatively:
  # b <- rep(0,length=length(x))
  # for (i in 1:j){
  #   h[j] <- x[j] - b[j]
  # }
  
  A <- diag(rep(1,length=length(x)))
  b <- rep(0,length=length(x))
  h <- as.vector(A%*%x - b)
  
  return(h)
}

# Equality
heq <- function(x){
  h <- sum(x) - 1
  return(h)
}

# ------Objective Function: Minimize Variance-----
eval_f <- function(x){
  s <- -t(x)%*%R + lambda / 2 * t(x) %*% Sigma %*% x
  return(s)
}

# Set the starting valae to a feasible portfolio:
x0 <- rep(1/n,length=n)

# Solve the optimization problem:
#   NOTE: the last two "trace" arguments in control.outer 
#   and countrol.optim turn off the output so it doesn't print
#   to the screen.  It is useful for the loops below.
tmp <- constrOptim.nl(par=x0, fn=eval_f,  
                      heq=heq, hin=hin,
                      "control.outer"=list("trace"=FALSE),
                      "control.optim"=list("reltol"=1e-12,
                                           "trace"=1)) 

# Print the return and risk of the estimated minimum variance portfolio:
Rmin <- tmp$par %*% R
print(Rmin)
Smin <- sqrt(S_pf(tmp$par,Sigma))
print(Smin)

message(sprintf("Maximum utility Portfolio:\n\tR = %5.4f\n\tSigma=%5.4f",Rmin,Smin))
p1 <- p 
p2 <- p1 + geom_point() + geom_point(aes(x=Smin,y=Rmin),colour="red")
print(p2)

# Plot contour lines over the graph:
get_plot_limits <- function(plot) {
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}
lmts <- get_plot_limits(p2)
x.range <- seq(lmts$xmin,lmts$xmax,length=100)
y.range <- seq(lmts$ymin,lmts$ymax,length=100)
df.grid <- expand.grid(x=x.range,y=y.range)
df.grid$f <- df.grid$y - lambda / 2 * (df.grid$x)^2


p3 <- p2
p4 <- p3  +
  stat_contour(data = df.grid, aes(x = x, y = y, z = f, colour = ..level..),bins=50)
print(p4)
ggsave("lambda_eq_1_5.pdf", plot = last_plot(), 
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)



# ---------------
# ----REPEAT FOR LARGER LAMBDA

plot_utility <- function(lambda){
  # Inequality
  hin <- function(x){
    # h <- x
    
    # # Or, alternatively:
    # b <- rep(0,length=length(x))
    # for (i in 1:j){
    #   h[j] <- x[j] - b[j]
    # }
    
    A <- diag(rep(1,length=length(x)))
    b <- rep(0,length=length(x))
    h <- as.vector(A%*%x - b)
    
    return(h)
  }
  
  # Equality
  heq <- function(x){
    h <- sum(x) - 1
    return(h)
  }
  
  # ------Objective Function: Minimize Variance-----
  eval_f <- function(x){
    s <- -t(x)%*%R + lambda / 2 * t(x) %*% Sigma %*% x
    return(s)
  }
  
  # Set the starting valae to a feasible portfolio:
  x0 <- rep(1/n,length=n)
  
  # Solve the optimization problem:
  #   NOTE: the last two "trace" arguments in control.outer 
  #   and countrol.optim turn off the output so it doesn't print
  #   to the screen.  It is useful for the loops below.
  tmp <- constrOptim.nl(par=x0, fn=eval_f,  
                        heq=heq, hin=hin,
                        "control.outer"=list("trace"=FALSE),
                        "control.optim"=list("reltol"=1e-12,
                                             "trace"=1)) 
  
  # Print the return and risk of the estimated minimum variance portfolio:
  Rmin <- tmp$par %*% R
  print(Rmin)
  Smin <- sqrt(S_pf(tmp$par,Sigma))
  print(Smin)
  
  message(sprintf("Maximum utility Portfolio:\n\tR = %5.4f\n\tSigma=%5.4f",Rmin,Smin))
  p1 <- p 
  p2 <- p1 + geom_point() + geom_point(aes(x=Smin,y=Rmin),colour="red")
  print(p2)
  
  # Plot contour lines over the graph:
  lmts <- get_plot_limits(p2)
  x.range <- seq(lmts$xmin,lmts$xmax,length=100)
  y.range <- seq(lmts$ymin,lmts$ymax,length=100)
  df.grid <- expand.grid(x=x.range,y=y.range)
  df.grid$f <- df.grid$y - lambda / 2 * (df.grid$x)^2
  
  
  p3 <- p2
  p4 <- p3  +
    stat_contour(data = df.grid, aes(x = x, y = y, z = f, colour = ..level..),bins=50)
  print(p4)
  fname <- paste0("lambda_eq_",lambda,".pdf")
  ggsave(fname, plot = last_plot(), 
         scale = 1, width = 5, height = 4, units = "in",
         dpi = 300, limitsize = TRUE)
}
plot_utility(8)
