# masked p-value generator
# P:         vector of p-values of length n
# mask_fun:  name of the masking function
# mask_para: a vector of parameters in the masking function
# g_num_map: mapping g computed numerically for a continous masking (object of g_map_gen)
mask_gen = function(P, mask_fun, mask_para, g_num_map = NULL, centered = FALSE){
  if (mask_fun == "tent"){
    h <- 1*(P < mask_para) + (-1)*(P >= mask_para)
    g <- pmin(P, mask_para/(1 - mask_para)*(1 - P))
  } else if (mask_fun == "railway"){
    h <- 1*(P < mask_para) + (-1)*(P >= mask_para)
    g <- P*(P < mask_para) +
      mask_para/(1 - mask_para)*(P - mask_para)*(P >= mask_para) 
  } else if (mask_fun == "continuous"){
    if(mask_para > 0){
      h <- log(mask_para*P^(mask_para - 1))
      if (centered) {h = h - (1 - mask_para + log(mask_para))}
      g <- g_num_map(P, mask_para)
    } else {
      ori_h <- (1 - P)/(log(P))^2/P + 1/(log(P))
      h <- log(ori_h); if (centered) {h = h - (-0.29)}
      g <- g_num_map(P, mask_para)
    }
  } else if (mask_fun == "Gaussians"){
    Z = qnorm(1 - P); Z_noise = rnorm(length(P))
    h = (Z + mask_para*Z_noise)/sqrt(1 + mask_para^2)
    g = (Z - 1/mask_para*Z_noise)/sqrt(1 + 1/mask_para^2); g = -g #smaller g indicate non-nulls
  }
  return(list(h = h, g = g))
}


# generate the mapping of g(P) for a continuous masking
g_map_gen = function(){
  for (eps in seq(0, 0.8, 0.2)) {
    if (eps > 0){
      p_star <- eps^(-1/(eps - 1))
      H <- function(p) {p^eps - p}
    } else {
      p_star <- uniroot(function(x){(1 - x)/(log(x))^2/x + 1/(log(x)) - 1}, c(0.1, 0.9))$root
      H <- function(p) {(p - 1)/log(p) - p}
    }
    g_num_val = sapply(seq(0, 0.999, length.out = 1000),
                       function(p) {uniroot(function(x){H(x) - H(p)}, c(0, p_star))$root})
    g_num_val = c(g_num_val, 0)
    save(g_num_val, file=paste("intermediate_results/g_num_val_",eps*10,".Rdata", sep = ""))
  }
}

# g(P) for a given p value and parameter eps
g_num_map = function(x, eps) {
  load(paste("intermediate_results/g_num_val_",eps*10,".Rdata", sep = ""))
  return(g_num_val[round(x*1000) + 1])
} 
g_num_map = Vectorize(g_num_map)

