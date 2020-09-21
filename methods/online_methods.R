# online martingale Stouffer test (MST-online)
# P:          vector of p-values
# alpha:      target Type 1 error level
# ub_func:    function for the uniform upper bound
# n_included: number of included hypotheses
# sum_stat:   sum statistics of included hypotheses
mst_online = function(P, alpha, ub_func, n_included, sum_stat){
  n_passed = length(P)
  cumsum_stat = sum_stat + cumsum(qnorm(1 - P))
  exceed_time = which(cumsum_stat > 
                        ub_func((n_included + 1):(n_included + n_passed), alpha = alpha))
  if(length(exceed_time)){
    rej = TRUE
    n_included = n_included + exceed_time[1]
    n_passed = exceed_time[1]
  } else{
    rej = FALSE
    n_included = n_included + n_passed
    sum_stat = cumsum_stat[n_passed]
  }
  return(list(rej = rej, n_included = n_included, n_passed = n_passed, sum_stat = sum_stat))
}


# online adaptive martingale test
# P:          vector of p-values 
# alpha:      target Type 1 error level
# n_included: number of included hypotheses
# h_sum:      sum of included h(P)
# mask_fun:   name of the masking function (tent/railway/continous)
# mask_para:  parameter in the masking function 
# ub_func:    function for the uniform upper bound
# threshold:  threshold for filtering the hypotheses
adaptive_mt_online = function(P, alpha, ub_func, n_included, h_sum,
                              mask_fun, mask_para, threshold){
  h_g = mask_gen(P, mask_fun, mask_para); h = h_g$h; g = h_g$g
  ind_select = which(g < threshold); n_select = length(ind_select)
  if (n_select == 0) {
    return(list(rej = FALSE, n_included = 0, n_passed = length(P), h_sum = h_sum))
  } else {
    #ind_select = ind_select[order(g[ind_select])] this is mini-batch setting
    cumsum_h = h_sum + cumsum(h[ind_select])
    exceed_time = which(cumsum_h > 
                          ub_func((n_included + 1):(n_included + n_select), alpha = alpha))
    if(length(exceed_time)){
      rej = TRUE
      n_included = n_included + exceed_time[1]
      n_passed = ind_select[exceed_time[1]]
    } else{
      rej = FALSE
      n_included = n_included + n_select
      n_passed = length(P)
      h_sum = cumsum_h[n_select]
    }
    return(list(rej = rej, n_included = n_included, n_passed = n_passed, h_sum = h_sum))
  }
}

#********** to be changed
# online interactive martingale test
# P:          vector of p-values 
# x:          side information 
# alpha:      target Type 1 error level
# n_included: number of included hypotheses
# h_sum:      sum of included h(P)
# mask_fun:   name of the masking function (tent/railway/continous)
# mask_para:  parameter in the masking function 
# ub_func:    function for the uniform upper bound
# threshold:  threshold for filtering the hypotheses
interactive_mt_online = function(P, x, alpha, ub_func, n_included, h_sum, structure,
                              mask_fun, mask_para, threshold){
  n_passed = length(P)
  h_g = mask_gen(P, mask_fun, mask_para); h = h_g$h; g = h_g$g
  m_set = rep(FALSE, n_passed)
  S = em_mixture(masked_P = g, x = x, m_set = m_set,
                 mask_fun = mask_fun, mask_para = mask_para, structure = structure)
  ind_select = which(S > threshold); n_select = length(ind_select)
  ind_select = ind_select[order(g[ind_select])] 
  cumsum_h = h_sum + cumsum(h[ind_select])
  exceed_time = which(cumsum_h > 
                        ub_func((n_included + 1):(n_included + n_select), alpha = alpha))
  if(length(exceed_time)){
    rej = TRUE
    n_included = n_included + exceed_time[1]
    n_passed = ind_select[exceed_time[1]]
  } else{
    rej = FALSE
    n_included = n_included + n_select
    h_sum = cumsum_h[n_select]
  }
  return(list(rej = rej, n_included = n_included, n_passed = n_passed, h_sum = h_sum))
}



# LORD (MST-online)
# P:          vector of p-values
# alpha:      target Type 1 error level
# ub_func:    function for the uniform upper bound
# n_included: number of included hypotheses
# sum_stat:   sum statistics of included hypotheses
lord = function(P, alpha, ub_func, n_included){
  n_passed = length(P)
  exceed_time = which(P < ub_func((n_included + 1):(n_included + n_passed), alpha = alpha))
  if(length(exceed_time)){
    rej = TRUE
    n_included = n_included + exceed_time[1]
    n_passed = exceed_time[1]
  } else{
    rej = FALSE
    n_included = n_included + n_passed
  }
  return(list(rej = rej, n_included = n_included, n_passed = n_passed))
}




