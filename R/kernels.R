
# Bartlett
bartlett <- function(x){
  if(abs(x)<1){
    k_x <- 1-abs(x)
  } else{
    k_x <- 0
  }
  return(k_x)
}

# Quadratic Spectral
qs <- function(x){
  p1 <- sin(6*pi*x/5)/(6*pi*x/5)
  p2 <- cos(6*pi*x/5)
  p3 <- 25/(12*pi^2*x^2)
  k_x <- p3*(p1-p2)
  if(x == 0){
    k_x <- 1
  }
  return(k_x)
}

# Parzen
parzen <- function(x){
  k_x <- 0
  if(abs(x)<=1){
    if(abs(x)<= 0.5){
      k_x <- 1 - 6*x^2 + 6*abs(x)^3
    } else {
      k_x <- 2*(1-abs(x))^3
    }
  }
  return(k_x)
}


# Tukey-Hanning
th <- function(x){
  k_x <- 0
  if(abs(x)<=1){
    k_x <- 0.5*(1 - cos(pi*x))
  }
  return(k_x)
}
