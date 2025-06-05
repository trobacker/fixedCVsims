# mother kernel g_q
# zero lugsail g_q = 0
# over lugsail is the negative of mother
g_q <- list("bartlett" = 1, "parzen" = 6, "th" = pi^2/4, "qs" = 1.421223)



# Over bandwidth rule
b_rule <- function(rho, big_T, alpha, d, w_q, g_q, q=1, tau =  alpha*.15, auto_adjust = T){

  try_b <- (0:(big_T)/2)/big_T
  cv <- qchisq((1-alpha), d)

  # Type 1 error
  distortion <- dchisq(cv, d)*cv*2*rho^2*rho^(try_b*big_T)/(1 + rho)+
    (try_b*big_T)^(-q)*dchisq(cv, d)*cv*g_q*w_q
  type_1 <- alpha + distortion
  opt_b <- try_b[which(abs(distortion) <= tau)]

  if(length(opt_b) == 0 & auto_adjust){
    # If we don't find a suitable b and auto adjust is allowed,
    # then increase tolerance (tau) until we meet the optimal b.
    tau_new <- tau
    for(i in 1:1000){ # control how many times tau can increase
      tau_new  <- tau_new*1.25
      opt_b <- try_b[which(abs(distortion) <= tau_new)]
      if(length(opt_b) >0) {
        opt_b <- min(opt_b)
        break
      }
    }
    warning(paste("No bandwidth met the criteria using the tolerance level ",
                  round(tau, 4), ". Tolerance was increased to ", round(tau_new, 4),
                  " using the auto adjust feature.", sep = ""))
  } else if(length(opt_b) == 0 & auto_adjust == F){
  warning(paste("No bandwidth met the criteria using the tolerance level. It is recommended to increase the tolerance level or use the auto
                adjust feature to automatically compensate.", sep = ""))
    b_opt = 0
  } else{
    opt_b <- min(opt_b)
  }
  plot(try_b, abs(distortion), xlab = "Bandwidth (b)",
       ylab = "Abs. Distortion")
  abline(v= opt_b, col = "red", lwd = 2)
  legend("topright", col = "red", lwd = 2, legend = "Optimal", cex  =.5)
  return(opt_b)
}



# Main ---------------------------------------------------------

# tau = alpha * .15
get_b <- function(the_data, alpha = 0.05, the_kernel ="Bartlett", lugsail="Mother", tau = alpha*.15,
                  auto_adjust = T){
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
 # tau <- -1/ (big_T * log(rho)) # delete later

  # kernel statistic information
  q <- 1
  if(the_kernel != "bartlett"){q <- 2}
  if(q == 1){
    w_q <- 2*rho/(1-rho^2)
  } else {w_q <- 2*rho/(1-rho)^2 }


  # g_q based on lugsail type
  g_q <- g_q[[the_kernel]]
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

  #tau <- -.25/(big_T*log(rho)) # delete this line
  b_opt <- b_rule(rho, big_T, alpha = 0.05, d =d, w_q, g_q, tau = tau, auto_adjust = auto_adjust)

  return(b_opt)
}



