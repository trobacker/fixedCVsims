# mother kernel g_q
# zero lugsail g_q = 0
# over lugsail is the negative of mother
g_q <- list("bartlett" = 1, "parzen" = 6, "th" = pi^2/4, "qs" = 1.421223)



# Over bandwidth rule
b_rule <- function(rho, big_T, alpha, d, w_q, g_q, q=1, tau =  alpha*.15){
  try_b <- (1:(big_T/2))/big_T
  cv <- qchisq((1-alpha), d)

  # Type 1 error
  distortion <- dchisq(cv, d)*cv*2*rho^2*rho^(try_b*big_T)/(1 + rho)+ (try_b*big_T)^(-q)*dchisq(cv, d)*cv*g_q*w_q
  type_1 <- alpha + distortion
  opt_b <- try_b[which(abs(distortion) <= tau)]
  opt_b <- min(opt_b)
  if(opt_b == Inf){
    warning("No bandwidth meets criteria, it is recommended to increase tolerance level.")
    opt_b <- 0
  }
  return(opt_b)
}



# Main ---------------------------------------------------------

get_b <- function(the_data, alpha = 0.05, the_kernel ="bartlett", lugsail="Mother", tau = alpha*.15){
  the_data <- as.matrix(the_data)
  # dimensions
  big_T <- nrow(the_data)
  d <- ncol(the_data)

  # calculate the rho
  # average the AR(1) coefficient for all dimensions
  all_rhos <- rep(0, d)
  for(i in 1:d){
    all_rhos[i] <- stats::acf(the_data[,i], plot = F)$acf[2]
  }
  rho <- mean(all_rhos)

  # kernel statistic information
  q <- 1
  if(the_kernel != "bartlett"){q <- 2}
  if(q == 1){
    w_q <- 2*rho/(1-rho^2)
  } else {w_q <- 2*rho/(1-rho)^2 }
  g_q <- g_q[[the_kernel]]

  # g_q based on lugsail type
  if(lugsail == "Zero"){
    g_q <- 0
  } else if (lugsail == "Over"){
    g_q <- - g_q
  }

  if(lugsail == "Mother"){
    if(tau < 0.1*alpha){
      warning("Recommended tolerance level for mother setttings is 0.1*alpha or higher.")
    }
  }

  if(lugsail == "Over"){
    if(rho <0.9 | big_T <200){
      warning("Over lugsail is not recommended for small data sets, or data sets with low correlation.")
    }
  }

  b_opt <- b_rule(rho, big_T, alpha = 0.05, d =d, w_q, g_q, tau = tau)
  return(b_opt)
}



