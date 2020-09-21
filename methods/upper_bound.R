# function of a linear uniform upper bound for sub-Gaussian
# k:     time point to obtain the upper bound
# alpha: Type 1 error level
# m:     parameter of the upper bound that change the tight region
ub_func_linear_normal = function(k, alpha, m = 10^4/4) {
  c = sqrt(-2*m*(log(alpha))) 
  ub_k = c + c/2/m*(k - m)
  return(ub_k)
}


# function of a curved uniform upper bound for sub-Gaussian using polynomial stitching
# k:     time point to obtain the upper bound
# alpha: Type 1 error level
ub_func_curve_poly = function(k, alpha) {
  ub_k = 1.7*sqrt(k*(log(log(2*k)) + 0.72*log(5.19/alpha)))
  return(ub_k)
}


# function of a curved uniform upper bound for sub-Gaussian using discrete mixture
# k:     time point to obtain the upper bound
# alpha: Type 1 error level
# f_mix: an intermediate function to compute the bound
ub_func_curve_mix = function(L, alpha, f_mix = f_mix) {
  eta = 1.1; s = 1.4
  lam_max <- sqrt(2*log(1/alpha))
  k_max <- ceiling(log(lam_max*(sqrt(5*L/log(1/alpha))), base = eta))
  
  lambda_seq <- lam_max / eta^(1:k_max + 1/2)
  w_seq <- lam_max*(eta - 1)*f_mix(lambda_seq*sqrt(eta), lam_max, s) / eta^(1:k_max + 1)
  
  p <- numeric(length = L)
  ub <- numeric(length = L)
  s = 0; expand = 1.001; shrink = 0.9
  for (t in 1:L) {
    s = max(sqrt(t), s)
    p[t] <- sum(w_seq*exp(lambda_seq*s - lambda_seq^2/2*t))
    while (p[t] > 1/alpha) {
      s = s*shrink
      p[t] <- sum(w_seq*exp(lambda_seq*s - lambda_seq^2/2*t))
    }
    while (p[t] < 1/alpha) {
      s = s*expand
      p[t] <- sum(w_seq*exp(lambda_seq*s - lambda_seq^2/2*t))
    }
    ub[t] <- s
  }
  return(ub)
}

f_mix <- function(x, lam_max, s){
  (s - 1) * (x <= lam_max) * (x >= 0) / x / log(exp(1)*lam_max/x)^s
}
#f_mix <- Vectorize(f_mix)


# a curved uniform upper bound for sub-Gaussian using inverted stitching
# L:     maximum index to compute the bound
# alpha: Type 1 error level
ub_func_curve_invert <- function(k, alpha){
  a = 2.42; b = 4.7
  ub = a*sqrt(k*log(log(exp(1)*k)) + b)
  return(ub)
}


# function of the bound in LORD
# k:     time point to obtain the upper bound
# alpha: Type 1 error level
ub_func_lord = function(k, alpha) {
  ub_k = alpha/3.524/k/(log(max(k,2)))^2
  return(ub_k)
} 
  

# function of a linear uniform upper bound for chi-square distribution
# k:     time point to obtain the upper bound
# alpha: Type 1 error level
# m:     parameter of the upper bound that change the tight region
ub_func_linear_chi <- function(k, alpha, m){
  #a = (sqrt(-2/m*(log(alpha))) + 1)*2*m
  #x = uniroot(function(x) {exp(-x/2 - m + m*log(x/2/m)) - alpha}, c(a, 2*a))$root
  x = uniroot(function(x) {exp(-x/2 - m + m*log(x/2/m)) - alpha}, c(0, 10^4))$root
  ub_k = x + (log(x/m/2) - 2)*(k - m) 
    #x + (log(x/m/2)/(0.5 - m/x) - 2)*(k - m) 
  return(ub_k)
}


# function of a linear uniform upper bound for sub-Exponential distribution
# k:     time point to obtain the upper bound
# alpha: Type 1 error level
# m:     parameter of the upper bound that change the tight region
ub_func_linear_exp <- function(k, alpha, m){
  c = sqrt(2); expand = 1.01
  x = sqrt(-2*m*(log(alpha)))
  p = exp(-x/c + m/c^2*log(1 + c*x/m))
  while (p > alpha) {    #expansion on x to make exact prob bound less than alpha
    x = expand*x
    p = exp(-x/c + m/c^2*log(1 + c*x/m))
  }
  # a = sqrt(-2*m*(log(alpha)))  
  # x = uniroot(function(x) {exp(-x/c + m/c^2*log(1 + c*x/m)) - alpha}, c(a, 2*a))$root
  b <- (1 + c*x/m)*log(1 + c*x/m) / c^2 / x*m - 1/c
  ub_scaled = x + b*(k - m) # a vector of length L
  ub_k = ub_scaled*2*sqrt(2) 
  return(ub_k)
}


# function of a curve uniform upper bound for sub-Gamma distribution
# k:     time point to obtain the upper bound
# alpha: Type 1 error level
ub_func_curve_gamma <- function(k, alpha){
  eta = 2.04; s = 1.4; c = sqrt(2)
  
  k1 = (eta^(1/4) + eta^(-1/4))/sqrt(2)
  k2 = sqrt(eta) + 1
  l = s*log(log(eta*k)) + log(zeta(s)/(alpha*(log(eta))^s))
  ub_scaled = k1*sqrt(k*l) + c*k2*l # a vector of length L
  ub_k = ub_scaled *2*sqrt(2) 
  return(ub_k)
}


## linear bound for variables bounded below using exponential 
# c: lower bound
ub_exp_para <- function(m, alpha, c = 2, expand = 1.01){
  x = sqrt(-2*m*(log(alpha)))   #initial x by taylor expansion to second degree
  p = exp(-x/c + m/c^2*log(1 + c*x/m))
  while (p > alpha) {    #expansion on x to make exact prob bound less than alpha
    x = expand*x
    p = exp(-x/c + m/c^2*log(1 + c*x/m))
  }
  b <- (1 + c*x/m)*log(1 + c*x/m) / c^2 / x*m - 1/c
  return(c(x, b))
}

ub_emp <- function(h_var, m, ub_para){
  ub <- ub_para[1] + ub_para[2]*(h_var - m)
  return(ub)
}



