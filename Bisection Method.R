#Bisection Method = a numerical method in Mathematics to find a root of a given function

myBis <- function(f,a,b,tol=1e-8,maxN=200){
  c <- (a + b)/2
  fa <- f(a)
  fb <- f(b)
  fc <- f(c)

  # Check to amke sure x* is on [a,b]
  if (fa*fb > 0) {
    stop("Solution not on [a,b]: widen search!")
  }
  
  if ((abs(a-b) < tol) | (abs(fc) < tol)){
    return(c)
    break
  }
  
  # Execute the main loop:
  for (i in 1:maxN){
    # print(c)
    # print(fc)
    if (fa*fc < 0){
      b <- c
      fb <- fc
    }
    else {
      a <- c
      fa <- fc
    }
    if ((abs(a-b) < tol) | (abs(fc) < tol)){
      message(sprintf("Solution achieved in %d iterations.",i))
      message(sprintf("\tx* = %5.4f",c))
      message(sprintf("\tf(x*) = %5.4f",fc))
      return(c)
      break
    }
    c <- (a + b)/2
    fc <- f(c)
  }
  message("Warning - reached maximum number of iterations!")
  return(c)
}