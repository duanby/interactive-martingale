library(parallel)
order_norm <- function(N1, mu, R){
  unorder <- matrix(rnorm(N1*R, mean = mu), nrow = R)
  sign <- unorder > 0
  ordered_abs <- t(apply(abs(unorder), 1, sort))[,N1:1] 
  ordered <- ordered_abs*(2*sign - 1)
  return(ordered)
}

unit_order <- function(N1, mu){
  unorder <- rnorm(N1, mean = mu)
  ordered_ind <- order(abs(unorder), decreasing = TRUE)
  ordered <- unorder[ordered_ind]
  p <- ordered > 0
  q <- abs(rnorm(N1)) > abs(ordered)
  return(c(p, q))
}


CaN <- function(alpha, n){
  c <- 1.7*sqrt(log(log(2*n)) + 0.72*log(5.19/alpha))
  return(c)
}
  

find_mu <- function(N1, N0, alpha, beta, R){
  n <- N1 + N0; C_alpha <- CaN(alpha, n); C_beta <- CaN(beta, n)
  mu = 0; rej <- FALSE
  while (!rej) {
    mu = mu + 0.1
    # ordered_mat <- order_norm(N1, mu, R)
    # 
    # p <- colMeans(ordered_mat > 0) #cannot move mu
    # q <- colMeans(abs(matrix(rnorm(N1*R), nrow = R)) > abs(ordered_mat))
    temp_fun <- function(x){ res <- unit_order(N1, mu = mu); return(res)}
    temp <- matrix(unlist(mclapply(1:R, temp_fun, mc.cores =  detectCores()/2)), ncol = R, byrow = FALSE)
    t <- rowMeans(temp); p <- t[1:N1]; q <- t[-(1:N1)]
    rej <- any(cumsum(2*p - 1) 
               - (C_alpha + C_beta)*sqrt(1:N1 + qbinom(beta/N1, N0, q, lower.tail = FALSE))
               >= 0)
  }
  return(mu)
}

R = as.numeric(commandArgs(TRUE))
alpha = 0.05; beta = 0.05/2
N1_seq <- round(10^(seq(2, 3, length.out = 10))); N0_seq <- round(10^(seq(2, 5, length.out = 10)))
mu_mat <- matrix(nrow = length(N1_seq), ncol = length(N0_seq))
i = 1; j = 1
for (N1 in N1_seq) {
  for (N0 in N0_seq) {
    mu_mat[i,j] <- find_mu(N1, N0, alpha, beta, R = R)
    print(c(N0, N1, mu_mat[i,j]))
    j = j + 1
  }
  i = i + 1; j = 1
}
save(mu_mat, file="rej-prob/mu_heatmap2.Rdata")

if(0){
  colnames(mu_mat) <- c(2, NA, NA, 3, NA, NA, 4, NA, NA, 5)
  rownames(mu_mat) <- round(seq(2, 3, length.out = 10),1)
  heatmap.2(mu_mat, Rowv=FALSE, Colv=FALSE, dendrogram="none",
            main="", xlab=expression(log(N[0])), ylab=expression(log(N[1])), density.info="none",
            #labCol = round(seq(2, 5, length.out = 10),1),
            #labRow = round(log(seq(100, 10^3, 100))/log(10), 1),
            srtCol=0,
            trace="none", key.xlab = expression(mu), key.title = "", cexRow = 1.2, cexCol = 1.2,
            lmat = rbind(c(3,4), c(2,1)))
  dev.off()
}







