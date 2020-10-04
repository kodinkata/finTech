rm(list=ls())

#Step 1: one year put bermuda option on a stock with a current price of $40.00, a strike price of $43.50, volatility of 30%, and risk-free rate of 5%
#Can be excercised anytime after 6 months (European option until the first 6 months and then becomes American option)

#Step 2: Can be excercised at the end of each quarter (3,6,9 months) or at expiration

#Step 3: one-year American put option on a stock with a current price of $40.00, strike price of $43.50, volatility of 30%, and assume a risk-free rate of 5%.  
  # a) Option Pricing using Cox-Ross-Rubinstein parameters
  # b) Option Pricing if the stock pays dividens (5%) before expiration
  # c) Option Pricing if the stock pays dividens ($2) before expiration

#Step 4: one-year American put option on a stock with a current price of $40.00, strike price of $43.50, and volatility of 30%.
  # - Changing risk free rates over time

#Step 5: binomial lattices to find the price of an American spread option on two stocks with current prices S1=S2 = 100, K = 1, r = 0.06, T = 1, sigma1 = 0.2, sigma2 = 0.3, cor (correlation between the assets) = 0.5, q1 (Divident yield) = 0.03, and q2 (Divident yield) = 0.04

#~~~~Step 1~~~~~~~

#Parameters:
S0<-40
K<-43.50
sigma<-0.3
r<-0.05
t<-1
N<-12
dt<-t/N

lattice<-matrix(NA,ncol=N+1,nrow=N+1)
#Cox Ross Rubinstein parameterization
u<-exp(sigma*sqrt(dt))
d<-1/u
p<-(exp(r*dt)-d)/(u-d)

#fill up the prices

for(j in 1:(N+1))
{
  for(i in 1:j)
  {
    lattice[i,j]<-S0*u^(j-i)*d^(i-1)
  }
}
lattice
profit<-matrix(NA,nrow=N+1,ncol=N+1)

#payoff for puts
h<-function(S)
{
    h<-pmax(K-S,0)
  
  return(h)
}  
#
profit[,N+1]<-h(lattice[,N+1])

disc <- exp(-r*dt)

for (j in N:1){
  tj <- dt *j
  for (i in 1:j){
    if(tj>0.5){
      profit[i,j]<-max(disc*(p*(profit[i,j+1]) + (1-p)*profit[i+1,j+1]),h(lattice[i,j]))
      
    }
    else{
      profit[i,j]<-disc*(p*(profit[i,j+1]) + (1-p)*profit[i+1,j+1])
    }
  }}
message(sprintf("The value of option is $%5.2f",profit[1,1]))

#~~~~~~~Step 2~~~~~~~~

for (j in N:1){
  tj <- dt *j
  for (i in 1:j){
    if(tj==3/12 || tj==6/12 || tj==9/12){
      profit[i,j]<-max(disc*(p*(profit[i,j+1]) + (1-p)*profit[i+1,j+1]),h(lattice[i,j]))
      
    }
    else{
      profit[i,j]<-disc*(p*(profit[i,j+1]) + (1-p)*profit[i+1,j+1])
    }
  }}
message(sprintf("The value of option is $%5.2f",profit[1,1]))

#~~~~~~~~Step 3~~~~~~~~
# a)
for (j in N:1){
  for (i in 1:j){
      profit[i,j]<-max(disc*(p*profit[i,j+1] + (1-p)*profit[i+1,j+1]),h(lattice[i,j]))
      
   }
  }
message(sprintf("The value of option is $%5.2f",profit[1,1]))

# b)
tdiv <- 6/12
div <- 0.05
for (j in 1:(N+1)){
  tj <- dt *(j)
  for (i in 1:j){
    if(tj<=tdiv){
          lattice[i,j]<-S0*u^(j-i)*d^(i-1)
    }
    else{
      lattice[i,j]<-S0*u^(j-i)*d^(i-1)*(1-div)
      
    }
  }
}
profit<-matrix(NA,ncol=N+1,nrow=N+1)

profit[,N+1]<-h(lattice[,N+1])

for (j in N:1){
  for (i in 1:j){
    profit[i,j] <- max(disc*(p*profit[i,j+1]+(1-p)*profit[i+1,j+1]),h(lattice[i,j]))
  }
}
message(sprintf("The value of option is $%5.2f",profit[1,1]))

# c)
D <- 2

S <- S0 - D*exp(-r*tdiv)

for (j in 1:(N+1)){
  tj <- dt *j
  PVD <- D*exp(-r*(tdiv-tj))
  for (i in 1:j){
    if(tj<=tdiv){
      lattice[i,j]<-S*u^(j-i)*d^(i-1)+PVD
    }
    else{
      lattice[i,j]<-S*u^(j-i)*d^(i-1)
      
    }
  }
}

profit[,N+1]<-h(lattice[,N+1])

for (j in N:1){
  for (i in 1:j){
    profit[i,j] <- max(disc*(p*profit[i,j+1]+(1-p)*profit[i+1,j+1]),h(lattice[i,j]))
  }
}
message(sprintf("The value of option is $%5.2f",profit[1,1]))

#~~~~~~~Step 4~~~~~~~~~
r<-0.05

lattice<-matrix(NA,ncol=N+1,nrow=N+1)

for(j in 1:(N+1))
{
  for(i in 1:j)
  {
    lattice[i,j]<-S0*u^(j-i)*d^(i-1)
  }
}
profit<-matrix(NA,nrow=N+1,ncol=N+1)

profit[,N+1]<-h(lattice[,N+1])

disc <- exp(-r*dt)

for (j in N:1){
  tj <- dt *(j-1)
  
    if(tj<=3/12){
      r <- .02
    }
    else if (tj<=6/12 && tj>3/12) {
      r <- .0225
    }
    else if (tj<=9/12 && tj>6/12) {
      r <- .025
      }
    else if(tj>9/12 && tj <1){
      r <- .0275
    }
  disc <- exp(-r*dt)
  
    p<-(exp(r*dt)-d)/(u-d)
    for(i in 1:j){
    profit[i,j]<-max(disc*(p*(profit[i,j+1]) + (1-p)*profit[i+1,j+1])
                     ,h(lattice[i,j]))
  }
}


message(sprintf("The value of option is $%5.4f",profit[1,1]))

#~~~~~~~Step 5~~~~~~~~~~~
N <-3
S1 <- S2 <- 100
K <- 1
r <- 0.06
t <-1
sigma1 <- 0.2
sigma2 <- 0.3
rho <- 0.5
q1 <- .03
q2 <- .04
dt <- t/N
nu1 <- r-q1-0.5*sigma1^2
nu2 <- r-q2-0.5*sigma2^2

u1 <- exp(sigma1*sqrt(dt))
u2 <- exp(sigma2*sqrt(dt))
d1 <- 1/u1
d2 <- 1/u2
disc <-exp(-r*dt)
puu <- disc*0.25 * (1+sqrt(dt)* (nu1/sigma1 + nu2/sigma2) + rho)
pud <- disc*0.25 * (1+sqrt(dt)* (nu1/sigma1 - nu2/sigma2) - rho)
pdu <- disc*0.25 * (1+sqrt(dt)* (-nu1/sigma1 + nu2/sigma2) - rho)
pdd <- disc*0.25 * (1+sqrt(dt)* (-nu1/sigma1 - nu2/sigma2) + rho)

lattice <- matrix(NA,ncol=2*N+1,nrow=2*N+1)

# Generate two lattices to visualize the prices of each stock:
LS1 <- matrix(NA,ncol=2*N+1,nrow=2*N+1)
LS2 <- matrix(NA,ncol=2*N+1,nrow=2*N+1)
S1vec <- S1*u1^(seq(N,0))*d1^(seq(0,N))
S2vec <- S2*u2^(seq(N,0))*d2^(seq(0,N))
for (i in 1:(N+1)){
  for (j in 1:(N+1)){
    LS1[2*i-1,2*j-1] <- S1vec[i]
  }
}
for (i in 1:(N+1)){
  for (j in 1:(N+1)){
    LS2[2*i-1,2*j-1] <- S2vec[j]
  }
}
print(LS1,digits=4)
h<-function(S1,S2)
{
  h<-pmax(S1 - S2 -K,0)
  return(h)
  
}  
for (i in 1:(N+1)){
  for (j in 1:(N+1)){
    lattice[2*i-1,2*j-1] <- h(S1vec[i],S2vec[j])
  }
}

lattice2 <-  matrix(NA,ncol=N+1,nrow = N+1)

lattice2 <- lattice

#Solve recursively:
for(k in 1:N)
{
  for(i in seq(1+k,2*N+1-k,by=2))
  {
    for(j in seq(1+k,2*N+1-k,by=2))
    {
      fuu<-lattice[i-1,j-1]
      fud<-lattice[i-1,j+1]
      fdu<-lattice[i+1,j-1]
      fdd<-lattice[i+1,j+1]
      LS1[i,j] <- LS1[i-1,j-1]/u1
      LS2[i,j] <- LS2[i-1,j-1]/u2
      #American Option
      lattice[i,j]<- max(puu*fuu + pdu*fdu +pud*fud +pdd*fdd,h(LS1[i,j],LS2[i,j]))
    }
  }  
}  

message(sprintf("The value of option is $%5.4f",lattice[N+1,N+1]))


