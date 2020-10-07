# experiment for offline settings: a sequence/grid/tree of hypotheses
expr_offline = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  if (structure == "grid") {
    n = D^2
    if ("weighted-HC" %in% methods) {
      guess_nonnull = nonnull_cluster(D = D, r = r, C = c(D/2, D/2))
      w0 = as.vector(1*guess_nonnull + 0.5*(!guess_nonnull)); w0 = w0/mean(w0)
      HC_dist = null_HC(w0)
      save(HC_dist, file="intermediate_results/HC_dist.Rdata")
    }
  } else if (structure == "tree-wide") {
    n = sum(20*n_child^(0:(levels-1))) + 1
  } else if (structure %in% c("tree-decrease", "tree-increase")) {
    n = sum(n_child^(0:(levels)))
  } 
  if(ub_type == "constant") {
    ub = ub_func_linear_normal(1:n, alpha = alpha, m = n/4)
    ub_eps = rep(log(1/alpha), n)
  } else if (ub_type == "linear") {
    ub = ub_func_linear_normal(1:n, alpha = alpha, m = m) 
    ub_eps = ub_func_linear_normal(1:n, alpha = alpha, m = m)
  }
  
  wrapper_func = function(i) {
    print(i)
    if (structure == "sequence_beta") {
      P = sim_seq(n = n, mu_1 = mu_1, mu_0 = mu_0,
                  pos_type = "beta", n_1 = n_1, beta_para = beta_para)
      dat = data.frame(P = P, x = 1:n); structure_para = NULL
    } else if (structure == "sequence_gather") {
      P = sim_seq(n = n, mu_1 = mu_1, mu_0 = mu_0,
                  pos_type = "gather", n_1 = n_1, sparsity = sparsity)
      dat = data.frame(P = P, x = 1:n); structure_para = NULL
    } else if (structure == "sequence_even") {
      P = sim_seq(n = n, mu_1 = mu_1, mu_0 = mu_0,
                  pos_type = "even", n_1 = n_1)
      dat = data.frame(P = P, x = 1:n); structure_para = NULL
    } else if (structure == "grid") {
      nonnull_ind = nonnull_cluster(D = D, r = r, C = C, mu_type = mu_type)
      dat = sim_cluster(nonnull_ind, mu_1 = mu_1, mu_0 = mu_0, rho = rho, mu_type = mu_type)
      structure_para = list(d = d, delta = delta)
    } else if (structure == "tree-wide") {
      tree_frame = tree_structure(levels = levels, n_child = n_child)
      dat = sim_tree(tree_frame, n_1 = n_nonnull, mu_1 = mu_1, mu_0 = mu_0,
                     decrease = decrease, tree_type = "wide")
      structure_para = tree_frame
    } else if (structure == "tree-decrease") {
      tree_frame = tree_structure(levels = levels, n_child = n_child, tree_type = "small")
      dat = sim_tree(tree_frame, n_1 = n_nonnull, mu_1 = mu_1, mu_0 = mu_0,
                     decrease = TRUE, tree_type = "small")
      structure_para = tree_frame
    } else if (structure == "tree-increase") {
      tree_frame = tree_structure(levels = levels, n_child = n_child, tree_type = "small")
      dat = sim_tree(tree_frame, n_1 = n_nonnull, mu_1 = mu_1, mu_0 = mu_0,
                     decrease = FALSE, tree_type = "small")
      structure_para = tree_frame
    }
    
    rejections = vector(length = length(methods)); names(rejections) = methods
    if ("Stouffer-batch" %in% methods) {
      rejections["Stouffer-batch"] =
        pnorm(sum(qnorm(1 - dat$P)), sd = sqrt(n)) > (1 - alpha)
    }
    if ("MST" %in% methods) {
      rejections["MST"] = 
        mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub,
            structure = structure, structure_para = structure_para)
    }
    if ("AMT-batch" %in% methods) {
      rejections["AMT-batch"] = 
        adaptive_mt(P = dat$P, alpha = alpha, ub = ub,
                    mask_fun = "tent", mask_para = 0.5)
    }
    if ("AMT-railway" %in% methods) {
      rejections["AMT-railway"] = 
        adaptive_mt(P = dat$P, alpha = alpha, ub = ub,
                    mask_fun = "railway", mask_para = 0.5)
    }
    if ("AMT-online" %in% methods) {
      rejections["AMT-online"] = 
        adaptive_mt_online(P = dat$P, alpha = alpha,
                           ub_func = function(k, alpha) ub_func_linear_normal(k, alpha, m = n/20), #n/25
                           n_included = 0, h_sum = 0, threshold = threshold,
                           mask_fun = "tent", mask_para = 0.5)$rej
    }
    if ("IMT-tent" %in% methods) {
      rejections["IMT-tent"] = 
        interactive_mt(P = dat$P, x = dat$x, alpha = alpha, ub = ub,
                       structure = structure, structure_para = structure_para,
                       mask_fun = "tent", mask_para = 0.5, S_model = smoothed)
    }
    if ("IMT-railway" %in% methods) {
      rejections["IMT-railway"] = 
        interactive_mt(P = dat$P, x = dat$x, alpha = alpha, ub = ub,
                       structure = structure, structure_para = structure_para,
                       mask_fun = "railway", mask_para = 0.5, S_model = smoothed)
    }
    if ("weighted-HC" %in% methods) {
      load(file = "intermediate_results/HC_dist.Rdata")
      pwval = cal_cdf(dat$P, w=w0)
      HC_val = hc_cal(pwval,t0ratio = 0.4)
      rejections["weighted-HC"] = HC_val > quantile(HC_dist, 0.95)
    }
    if ("AW-Fisher" %in% methods) {
      rejections["AW-Fisher"] =
        AWFisher_pvalue(matrix(dat$P, ncol = length(dat$P)))$pvalues < alpha
    }
    if ("IMT-Gaussians" %in% methods) {
      rejections["IMT-Gaussians"] = 
        interactive_mt(P = dat$P, x = dat$x, alpha = alpha, ub = ub,
                       structure = structure, structure_para = structure_para,
                       mask_fun = "Gaussians", mask_para = 0.1, S_model = FALSE)
    }
    for (eps in seq(0, 0.8, 0.2)) {
      if (paste("IMT-",eps*10, sep = "") %in% methods) {
        rejections[paste("IMT-",eps*10, sep = "")] =
          interactive_mt(P = dat$P, x = dat$x, alpha = alpha, ub = ub_eps,
                         structure = structure, structure_para = structure_para, S_model = FALSE,
                         mask_fun = "continuous", mask_para = eps, centered = centered)
      }
    }
    return(rejections)
  }
  #rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  rejections = lapply(1:R, wrapper_func)
  return(rejections)
}


expr_sparsity = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  ub = ub_func_linear_normal(1:n, alpha = alpha, m = n/4)
  if (ub_type == "linear" & "MST (m = n/4)" %in% methods) {
    ub_1 = ub_func_linear_normal(1:n, alpha = alpha, m = n/4)
    ub_2 = ub_func_linear_normal(1:n, alpha = alpha, m = n/2)
    ub_3 = ub_func_linear_normal(1:n, alpha = alpha, m = 3*n/4)
    ub_oracle = list()
    for(i in 1:length(sparsity)) {
      ub_oracle[[i]] = ub_func_linear_normal(1:n, alpha = alpha, m = sparsity[i]*n)
    }
    ub_eps = ub_func_linear_normal(1:n, alpha = alpha, m = n/4)
    #ub_exp_para(m = n/4, alpha = alpha, c = 2) 
  } else if (ub_type == "linear" & "MFT (m = n/4)" %in% methods) {
    ub_1 = ub_func_linear_exp(1:n, alpha = alpha, m = n/4)
    ub_2 = ub_func_linear_exp(1:n, alpha = alpha, m = n/2)
    ub_3 = ub_func_linear_exp(1:n, alpha = alpha, m = 3*n/4)
    ub_oracle = list()
    for(i in 1:length(sparsity)) {
      ub_oracle[[i]] = ub_func_linear_exp(1:n, alpha = alpha, m = sparsity[i]*n)
    }
  } else if (ub_type == "curve" & "MST (linear)" %in% methods) {
    ub_1 = ub_func_linear_normal(1:n, alpha = alpha, m = n/4)
    ub_2 = ub_func_curve_poly(1:n, alpha = alpha)
    ub_3 = ub_func_curve_mix(L = n, alpha = alpha, f_mix = f_mix)
    ub_4 = ub_func_curve_invert(1:n, alpha = alpha)
  } else if (ub_type == "curve" & "MFT (linear)" %in% methods) {
    ub_1 = ub_func_linear_exp(1:n, alpha = alpha, m = n/4)
    ub_2 = ub_func_curve_gamma(1:n, alpha = alpha)
  }
  
  wrapper_func = function(i) {
    print(i)
    P = sim_seq(n = n, mu_1 = mu_1, mu_0 = mu_0,
                pos_type = "gather_vary", n_1 = n_1, sparsity = sparsity)
    rejections = list()
    for(i in 1:length(sparsity)) {
      dat = data.frame(P = P[i,], x = 1:n); structure_para = NULL
      
      rejections[[as.character(sparsity[i])]] = vector(length = length(methods))
      names(rejections[[as.character(sparsity[i])]]) = methods

      if ("MST (m = n/4)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (m = n/4)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_1,
              structure = "sequence_gather", structure_para = NA)
      }
      if ("MST (m = n/2)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (m = n/2)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_2,
              structure = "sequence_gather", structure_para = NA)
      }
      if ("MST (m = 3n/4)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (m = 3n/4)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_3,
              structure = "sequence_gather", structure_para = NA)
      }
      if ("MST (m = l, oracle)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (m = l, oracle)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_oracle[[i]],
              structure = "sequence_gather", structure_para = NA)
      }
      
      
      if ("MST (linear)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (linear)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_1,
              structure = "sequence_gather", structure_para = NA)
      }
      if ("MST (poly)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (poly)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_2,
              structure = "sequence_gather", structure_para = NA)
      }
      if ("MST (discrete)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (discrete)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_3,
              structure = "sequence_gather", structure_para = NA)
      }
      if ("MST (invert)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST (invert)"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub_4,
              structure = "sequence_gather", structure_para = NA)
      }
      
      
      if ("MFT (m = n/4)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MFT (m = n/4)"] = 
          mft(P = dat$P, alpha = alpha, ub = ub_1)
      }
      if ("MFT (m = n/2)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MFT (m = n/2)"] = 
          mft(P = dat$P, alpha = alpha, ub = ub_2)
      }
      if ("MFT (m = 3n/4)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MFT (m = 3n/4)"] = 
          mft(P = dat$P, alpha = alpha, ub = ub_3)
      }
      if ("MFT (m = l, oracle)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MFT (m = l, oracle)"] = 
          mft(P = dat$P, alpha = alpha, ub = ub_oracle[[i]])
      }
      
      
      if ("MFT (linear)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MFT (linear)"] = 
          mft(P = dat$P, alpha = alpha, ub = ub_1)
      }
      if ("MFT (curve)" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MFT (curve)"] = 
          mft(P = dat$P, alpha = alpha, ub = ub_2)
      }
      
      
      if ("Stouffer-batch" %in% methods) {
        rejections[[as.character(sparsity[i])]]["Stouffer-batch"] =
          pnorm(sum(qnorm(1 - dat$P)), sd = sqrt(n)) > (1 - alpha)
      }
      if ("MST" %in% methods) {
        rejections[[as.character(sparsity[i])]]["MST"] = 
          mst(P = dat$P, x = dat$x, alpha = alpha, ub = ub,
              structure = "sequence_gather", structure_para = NA)
      }
      if ("AMT-batch" %in% methods) {
        rejections[[as.character(sparsity[i])]]["AMT-batch"] = 
          adaptive_mt(P = dat$P, alpha = alpha, ub = ub,
                      mask_fun = "tent", mask_para = 0.5)
      }
      if ("AMT-online" %in% methods) {
        rejections[[as.character(sparsity[i])]]["AMT-online"] = 
          adaptive_mt_online(P = dat$P, alpha = alpha,
                             ub_func = function(k, alpha) ub_func_linear_normal(k, alpha, m = n/20), #n/25
                             n_included = 0, h_sum = 0, threshold = threshold,
                             mask_fun = "tent", mask_para = 0.5)$rej
      }
    }
    return(rejections)
  }
  #rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  rejections = lapply(1:R, wrapper_func)
  return(rejections)
}




# experiment for online settings
expr_online = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  ub_online = function(k, alpha) {
    ub_func_linear_normal(k = k, alpha = alpha, m = m) #100
  }
  wrapper_func = function(i) {
    print(i)
    rejections = rep(FALSE, length = length(methods)); names(rejections) = methods
    n_included = rep(0, length = length(methods)); names(n_included) = methods
    n_total = rep(0, length = length(methods)); names(n_total) = methods
    sum_stat = rep(0, length = length(methods)); names(sum_stat) = methods
    if (structure == "tree-online") {prob_leaves = prob_ini}
    
    while (any(!rejections) & max(n_total) < exp(10)) {
      if (structure == "sequence") {
        P = sim_seq(n = n_unit, mu_1 = mu_1, mu_0 = mu_0,
                    pos_type = "even", prob = prob)
        dat = data.frame(P = P, x = 1:n_unit) #x to be changed
      } else if (structure == "seq-block") {
        P = sim_seq(n = n_unit, mu_1 = mu_1, mu_0 = mu_0,
                    pos_type = "block", prob = prob)
        dat = data.frame(P = P, x = 1:n_unit) #x to be changed
      } else if (structure == "tree-online") {
        if (n_total[1] == 0) {
          n = length(prob_leaves)
          nonnull_ind <- sapply(prob_leaves, function(x) rbinom(1,1,x)) 
          dat_set <- rnorm(n) + nonnull_ind*mu_1 + (1 - nonnull_ind)*mu_0
          P <- 1 - pnorm(dat_set)
          dat = data.frame(P = P, x = prob_leaves)
        } else {
          dat = sim_online_tree(parent_probs = dat$x, mu_1 = mu_1, mu_0 = mu_0,
                                n_child = n_child, child_fac = child_fac)
        }
      }
      
      if ("lord" %in% methods) {
        temp = lord(P = dat$P, alpha = alpha, ub_func = ub_func_lord,
                    n_included = n_included["lord"])
        n_included["lord"] = temp$n_included
        n_total["lord"] = n_total["lord"] + temp$n_passed
        rejections["lord"] = temp$rej
      }
      if ("MST" %in% methods) {
        temp = mst_online(P = dat$P, alpha = alpha, ub_func = ub_online,
                          n_included = n_included["MST"], sum_stat = sum_stat["MST"])
        n_included["MST"] = temp$n_included
        n_total["MST"] = n_total["MST"] + temp$n_passed
        sum_stat["MST"] = temp$sum_stat
        rejections["MST"] = temp$rej
      }
      if ("AMT-online" %in% methods) {
        temp = adaptive_mt_online(P = dat$P, alpha = alpha, ub_func = ub_online,
                          n_included = n_included["AMT-online"], h_sum = sum_stat["AMT-online"],
                          mask_fun = mask_fun, mask_para = mask_para, threshold = threshold_amt)
        n_included["AMT-online"] = temp$n_included
        n_total["AMT-online"] = n_total["AMT-online"] + temp$n_passed
        sum_stat["AMT-online"] = temp$h_sum
        rejections["AMT-online"] = temp$rej
      }
      if ("IMT-online" %in% methods) {
        if (structure == "seq-block") {
          threshold_seq = rep(threshold_imt, length(dat$P))
          for (i in 11:length(dat$P)) {
            if (all(dat$P[(i-10):(i-1)] < 0.1)) {threshold_seq[i] = threshold_imt*2}
            if (all(dat$P[(i-10):(i-1)] > 0.1)) {threshold_seq[i] = threshold_imt/4}
          }
          temp = adaptive_mt_online(P = dat$P, alpha = alpha, ub_func = ub_online,
                                    n_included = n_included["IMT-online"], h_sum = sum_stat["IMT-online"],
                                    mask_fun = mask_fun, mask_para = mask_para, threshold = threshold_seq)
        } else {
          temp = interactive_mt_online(
            P = dat$P, x = dat$x, alpha = alpha, ub_func = ub_online, structure = structure,
            n_included = n_included["IMT-online"], h_sum = sum_stat["IMT-online"],
            mask_fun = mask_fun, mask_para = mask_para, threshold = threshold_imt)
        }
        n_included["IMT-online"] = temp$n_included
        n_total["IMT-online"] = n_total["IMT-online"] + temp$n_passed
        sum_stat["IMT-online"] = temp$h_sum
        rejections["IMT-online"] = temp$rej
      }
    }
    return(n_total)
  }
  #rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  rejections = lapply(1:R, wrapper_func)
  return(rejections)
}


