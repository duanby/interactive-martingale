# indicator of an ordered sequence of N(mu,1) > 0,
#             and another indicator for |N(mu,1)| < |N(0,1)| 
# N1:        number of nonnulls
# mu:        alternative mean value
unit_order <- function(N1, mu){
  unorder <- rnorm(N1, mean = mu)
  ordered_ind <- order(abs(unorder), decreasing = TRUE)
  ordered <- unorder[ordered_ind]
  p <- ordered > 0
  q <- abs(rnorm(N1)) > abs(ordered)
  return(c(p, q))
}

# parameter given level alpha and total number of hypotheses n
CaN <- function(alpha, n){
  c <- 1.7*sqrt(log(log(2*n)) + 0.72*log(5.19/alpha))
  return(c)
}

# grid search for minimum alternative mean with targeted Type I and Type II error control
# N1:    number of nonnulls
# N0:    number of nulls
# alpha: Type I error control 
# beta:  Type II error control 
# R:     number of repetition for the Monte Carlo method
find_mu <- function(N1, N0, alpha, beta, R){
  n <- N1 + N0; C_alpha <- CaN(alpha, n); C_beta <- CaN(beta, n)
  mu = 0; rej <- FALSE
  while (!rej) {
    mu = mu + 0.1
    temp_fun <- function(x){ res <- unit_order(N1, mu = mu); return(res)}
    temp <- matrix(unlist(mclapply(1:R, temp_fun, mc.cores =  detectCores()/2)), ncol = R, byrow = FALSE)
    t <- rowMeans(temp); p <- t[1:N1]; q <- t[-(1:N1)]
    rej <- any(cumsum(2*p - 1) 
               - (C_alpha + C_beta)*sqrt(1:N1 + qbinom(beta/N1, N0, q, lower.tail = FALSE))
               >= 0)
  }
  return(mu)
}










