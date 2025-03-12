#' Lugsail Transformation (main)
#' 
#' @param x numeric, autocovariance values
#' @param lugsail_parameters list(numeric r, numeric c), lugsail parameters which
#' can also be obtained from [get_lugsail_parameters].
#' @param the_kernel function, the mother kernel to be used
lugsail <- function(x, lugsail_parameters, the_kernel= bartlett){
  r <- lugsail_parameters$r
  c <- lugsail_parameters$c

  # Actual lugsail
  y1 <- the_kernel(x)/(1-c)
  y2 <- 0

  # If QS, then always include the extra bit.
  if(deparse(substitute(the_kernel)) == "qs"){
    y2 <- the_kernel(x*r)*c/(1-c)
  }

  # If not QS, then you need to check
  if(abs(x) < 1/r){
    y2 <- the_kernel(x*r)*c/(1-c)
  }
  y <- y1- y2

  return(y)
}

#' Lugsail Support Function to get lugsail parameters
#' (default) b = Andrews (1991) Rule: 0.75*big_T^(-2*q/(2*q+1))
#' q for Bartlett, Parzen, and QS are 1,2,2, respectively.
#' 
#' @param big_T numeric, length of time series
#' @param q numeric, characteristic component (Parzen 1957). q for Bartlett, 
#' Parzen, and QS are 1,2,2, respectively.
#' @param method character, one of "Over", "Adaptive", or default "Zero" (not unique)
#' @param b numeric, bandwidth parameter, default is b = 0.75*big_T^(-2*q/(2*q+1))
get_lugsail_parameters <- function(big_T, q, method = "Zero",
                                   b = 0.75*big_T^(-2*q/(2*q+1))){

  if(method == "Over"){
    r <- 3
    c <- 2/(1+r^q)

  } else if(method == "Adaptive"){
    r <- 2
    M  <- big_T * b
    c_num <- (log(big_T) - log(M) + 1)
    c_den <- r^q*(log(big_T) - log(M)) + 1
    c <- c_num/c_den

  } else {
    # Zero or Manual lugsail
    r <- 2
    c <- r^(-q)

  }
  parameters <- list(r = r, c = round(c, 2))
  return(parameters)
}

