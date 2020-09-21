# martingale Stouffer test (MST)
# P:         vector of p-values 
# x:         side information 
# alpha:     target Type 1 error level
# ub:        uniform upper bound for rejection
# structure: structure of the hypotheses (grid/tree)
# d, delta:  parameters of expanding M_k for the grid structure
mst = function(P, x, alpha, ub, structure, structure_para = NULL){
  n = length(P); h = qnorm(1 - P)
  if (structure %in% c("sequence_beta", "sequence_gather", "sequence_even",
                       "tree-wide", "tree-decrease", "tree-online")) { #in a tree structure P are traversed by level
    rej = any(cumsum(h) > ub)
  } else if (structure == "tree-increase") {
    rej = any(cumsum(h[n:1]) > ub)
  } else if (structure == "grid") {
    m_set = rep(FALSE, n)  # start with an empty set M_0
    iter = 0; rej = FALSE 
    while (!rej & any(!m_set)) {
      m_set = inclu_search_cluster(S = rep(1,n), x = x, m_set = m_set, random = TRUE,
                                   d = structure_para$d, delta = structure_para$delta)
      rej = sum(h[m_set]) > ub[sum(m_set)]
      iter = iter + 1
    }
  }
  return(rej)
}



# martingale Fisher test (MFT)
# P:         vector of p-values 
# alpha:     target Type 1 error level
# ub:        uniform upper bound for rejection
mft = function(P, alpha, ub){
  h = -2*log(P) - 2
  rej = any(cumsum(h) > ub)
  return(rej)
}



# adaptive martingale test (AMT)
# P:         vector of p-values
# alpha:     target Type 1 error level
# mask_fun:  name of the masking function (tent/railway/continous)
# mask_para: parameter in the masking function 
# ub:        uniform upper bound for rejection
adaptive_mt = function(P, alpha, mask_fun, mask_para, ub){
  # generate masked p-values and missing bits
  h_g = mask_gen(P, mask_fun, mask_para); h = h_g$h; g = h_g$g
  rej = any(cumsum(h[order(g)]) > ub)
  return(rej)
}


# interactive martingale test (IMT)
# P:         vector of p-values of length n
# x:         side information 
# alpha:     target Type 1 error level
# mask_fun:  name of the masking function
# mask_para: parameter in the masking function 
# ub:        uniform upper bound for rejection
# structure: structure of the hypotheses (grid/tree)
# S_model:   indicator of whether modelling the nonnull likelihood
# d, delta:  parameters of expanding M_k for the grid structure
# tree_obj:  the tree object when hypotheses form a tree
interactive_mt = function(P, x, alpha, mask_fun, mask_para, ub, structure, structure_para,
               S_model = TRUE, d = 5, delta = 0.05, tree_obj = NULL, centered = FALSE){
  if(mask_fun == "continuous" & length(ub) == 2) {
    ub_para = ub
  }
  # generate masked p-values and missing bits
  h_g = mask_gen(P, mask_fun, mask_para, centered); h = h_g$h; g = h_g$g
  
  n = length(P)
  if (structure %in% c("grid", "sequence_even", "tree-increase")) {
    m_set = rep(FALSE, n)  # start with an empty set M_0
  } else if (structure %in% c("tree-decrease", "tree-wide")) {
    m_set = rep(FALSE, n); m_set[1] = TRUE # start from the root
  }
  
  iter = 0; rej = FALSE 
  while (!rej & any(!m_set)) {
    masked_P = P*m_set + g*(1 - m_set)
    if (S_model) { 
      if (iter %% 100 == 0) { # update nonnull likelihood score S every 100 iterations
        S = em_mixture(masked_P, x, m_set, mask_fun, mask_para, structure)
      }
    } else {
      S = -masked_P
    }
    
    if(structure == "grid") {
      m_set = inclu_search_cluster(S = S, x = x, m_set = m_set, random = FALSE,
                                   d = structure_para$d, delta = structure_para$delta)
    } else if (structure %in% c("tree-decrease", "tree-increase", "tree-wide")) {
      m_set = inclu_search_tree(S = S, structure = structure, 
                                tree_obj = structure_para, m_set = m_set)
    } 
    if (mask_fun == "continuous" & length(ub) < 3) {
      ub = ub_emp(sum((h[m_set])^2), m = n/4, ub_para = ub_para)
      rej = sum(h[m_set]) > ub
    } else {
      rej = sum(h[m_set]) > ub[sum(m_set)]
      #print(c(sum(m_set), sum(h[m_set])))
    }
    iter = iter + 1
    #print(sum(m_set))
    if (0) {
      # if(iter == 0){
      #   record = data.frame(x = sum(m_set), y = sum(h[m_set]))
      # } else {
      #   record = rbind(record, data.frame(x = sum(m_set), y = sum(h[m_set])))
      # }
      # df = data.frame(x = 1:20, y = ub[1:20])
      # p = ggplot(data=df, aes(x=x, y=y, group=1)) + xlim(1,20) + ylim(0, 15) + 
      #   geom_line() + theme(legend.position = "none") + xlab("Number of included hypothesis") +
      #   ylab("") + 
      #   geom_point(data = record,
      #              mapping = aes(x = x, y = y, shape = "21", color = "red", size = 5))
      # ggsave(filename = paste("figure/line_", iter, ".png", sep = ""),
      #                plot = p, device = "png", width = 4, height = 4)
      # plot_pval = 1 - as.numeric(nonnull_ind)
      # longData = melt(matrix(plot_pval, nrow = D))
      # p = ggplot(longData, aes(x = Var2, y = Var1)) +
      #   geom_raster(aes(fill=value)) +
      #   scale_fill_gradient(low="violetred1", high= "slategray1", limits=c(0, 1)) +
      #   labs(title = bquote(t == .(iter)*","*~~~~~
      #                         "#cand" == .(sum(m_set))*","*~~~~~
      #                         sum(h) == .(sum(h[m_set])))) +
      #   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
      #                      axis.text.y=element_text(size=9),
      #                      axis.title.x=element_blank(), axis.title.y=element_blank(),
      #                      plot.title=element_text(size=11),
      #                      legend.position = "none",
      #                      panel.grid.major = element_blank(),
      #                      panel.grid.minor = element_blank()) +
      #   geom_vline(xintercept=seq(-0.5, 31.5, by=1), color = "grey") +
      #   geom_hline(yintercept=seq(-0.5, 31.5, by=1), color = "grey") +
      #   scale_x_continuous(breaks=seq(1, 10, 1)) +
      #   scale_y_continuous(breaks=seq(1, 10, 1)) +
      #   coord_cartesian(xlim = c(1, 10), ylim = c(1, 10))
      # plot(p)
      # ggsave(filename = paste("figure/true.png", sep = ""),
      #        plot = p, device = "png", width = 4, height = 4)
      # 
      # plot_pval = g*(!m_set) + P*m_set
      # longData = melt(matrix(plot_pval, nrow = D))
      # p = ggplot(longData, aes(x = Var2, y = Var1)) +
      #   geom_raster(aes(fill=value)) +
      #   scale_fill_gradient(low="violetred1", high= "slategray1", limits=c(0, 1)) +
      #   labs(title = bquote(t == .(iter)*","*~~~~~
      #                         "#cand" == .(sum(m_set))*","*~~~~~
      #                         sum(h) == .(sum(h[m_set])))) +
      #   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
      #                      axis.text.y=element_text(size=9),
      #                      axis.title.x=element_blank(), axis.title.y=element_blank(),
      #                      plot.title=element_text(size=11),
      #                      legend.position = "none",
      #                      panel.grid.major = element_blank(),
      #                      panel.grid.minor = element_blank()) +
      #   geom_vline(xintercept=seq(-0.5, 31.5, by=1), color = "grey") +
      #   geom_hline(yintercept=seq(-0.5, 31.5, by=1), color = "grey") +
      #   scale_x_continuous(breaks=seq(1, 10, 1)) +
      #   scale_y_continuous(breaks=seq(1, 10, 1)) +
      #   coord_cartesian(xlim = c(1, 10), ylim = c(1, 10))
      # plot(p)
      # ggsave(filename = paste("figure/gvals_", iter, ".png", sep = ""),
      #        plot = p, device = "png", width = 4, height = 4)
    }
  }
  return(rej)
}


# EM algorithm under mixture model to get the non-null likelihood score
# masked_P:  vector of masked p-value information
# x:         side information 
# rej_ind:   vector of indicators for whether a hypothesis is in the rejection set
# mask_fun:  name of the masking function
# mask_para: a vector of parameters in the masking function
# iter:      number of EM iterations, default to five
em_mixture = function(masked_P, x, m_set, mask_fun, mask_para, structure, decrease = FALSE,
                      iter = 5, df = 3){
  masked_P[masked_P < 10^(-10)] = 10^(-10)  ## avoid Inf in calculation
  masked_Z = qnorm(1 - masked_P) ## translate p-values into Z-scores
  if (mask_fun == "tent") {
    inv_Z = qnorm((1 - mask_para)/mask_para*(1 - pnorm(masked_Z)))
    drv = dnorm(masked_Z)/dnorm(inv_Z)*(1 - mask_para)/mask_para
  } else if (mask_fun == "railway") {
    inv_Z = qnorm((1 - mask_para)/mask_para*(pnorm(masked_Z) - 1 + mask_para)) 
    drv = dnorm(masked_Z)/dnorm(inv_Z)*(1 - mask_para)/mask_para
  } else if (mask_fun == "continuous") {
    
  }
  inv_Z[m_set] = 0
  
  if (structure == "tree-online") {
    pi_set = x
  } else {
    pi_set = rep(0.1, length(masked_P))
  }
  mu = 1 ## initial values for parameters in the EM algorithm
  para = vector(length = 4)
  for (i in 1:iter) {
    mu = mask_para %>% ifelse(mask_fun %in% c("tent", "railway"), ., .[1]) %>%
      max(mu, qnorm(1 - .) + 0.5)
    # likelihood of each case with respect to the latent labels w and q
    # a = pi_set*dnorm(masked_Z - mu); b = (1 - pi_set)*dnorm(masked_Z)
    # c = pi_set*dnorm(inv_Z - mu); d = (1 - pi_set)*dnorm(inv_Z)
    a = pi_set*dnorm(masked_Z - mu); b = (1 - pi_set)*dnorm(masked_Z)
    if (mask_fun %in% c("gap", "gap-railway")) {
      drv = dnorm(masked_Z)/dnorm(inv_Z)*(1 - mask_para[2])/mask_para[1]
    } else if (mask_fun %in% c("tent", "railway")) {
      drv = dnorm(masked_Z)/dnorm(inv_Z)*(1 - mask_para)/mask_para
    }
    #c = pi_set*dnorm(inv_Z - mu)*drv; d = (1 - pi_set)*dnorm(inv_Z)*drv
    
    a = ifelse(!m_set, a, a/(a+b)); a[is.na(a)] = 0.5
    b = ifelse(!m_set, b, 1 - a)
    c = ifelse(!m_set, pi_set*dnorm(inv_Z - mu)*drv, 0)
    d = ifelse(!m_set, (1 - pi_set)*dnorm(inv_Z)*drv, 0)
    sum_abcd = a + b + c + d
    q = (a+c)/sum_abcd
    a = a/sum_abcd; b = b/sum_abcd; c = c/sum_abcd; d = d/sum_abcd
    if (structure == "grid") {
      phi_x = bs(x[,1], df = df); phi_y = bs(x[,2], df); 
      phi = phi_x[,rep(1:df, each = df)] * phi_y[,rep(1:df, times = df)]
      pi_set = glm(q~phi, family = quasibinomial())$fitted.values
    } else if (structure == "tree-decrease") {
        pi_set = activeSet(isomat = x[,c(2,1)], mySolver = "LS", y = q,
                           weights = rep(1, length(masked_P)))$x
    } else if (structure == "tree-increase") {
      pi_set = activeSet(isomat = x, mySolver = "LS", y = q,
                         weights = rep(1, length(masked_P)))$x
    } else if (structure == "sequence") {
      phi = bs(x, df = df)
      pi_set = glm(q~phi, family = quasibinomial())$fitted.values
    } else if (structure == "tree-online") {
      pi_set = glm(q~x, family = quasibinomial())$fitted.values
    }
    mu = sum(a*masked_Z + c*inv_Z)/sum(a + c)
  }
  return(q)
}




# update M_k under clustered non-null structure 
# S:      non-null likelihood score
# x:      side information
# m_set:  vector of indicators for whether a hypothesis is in the previous M_k
# d:      number of cones
# delta:  proportion of hypothesis to be included in one cone
# random: indicator of whether include hypothesis from a random cone 
inclu_search_cluster = function(S, x, m_set, d, delta, random){
  if (sum(m_set) == 0 & random) {
    C = round(colMedians(x))
  } else if (sum(m_set) == 0 & !random) {
    S_mat = matrix(ncol = max(x[,1]), nrow = max(x[,2])); S_mat[x] = S
    aver_S <- apply(S_mat, 1, function(x, n = 2){stats::filter(x, rep(1/n, n), sides = 2)})
    aver_S <- apply(aver_S, 2, function(x, n = 2){stats::filter(x, rep(1/n, n), sides = 2)})
    C <- which(aver_S == max(aver_S,na.rm = T), arr.ind = TRUE); if(nrow(C) > 1){C = C[1,]}
  } else {
    if (sum(m_set) == 1) {
      C = x[m_set,]
    } else {
      C = round(colMedians(x[m_set,]))
    }
  } 
  # angle and distance from the center
  theta = atan((x[,1] - C[1]) / (x[,2] - C[2])) * (x[,1] - C[1] >= 0 & x[,2] - C[2] >= 0) +
    (atan((x[,1] - C[1]) / (x[,2] - C[2])) + pi) * (x[,2] - C[2] < 0) +
    (atan((x[,1] - C[1]) / (x[,2] - C[2])) + 2*pi) * (x[,1] - C[1] < 0 & x[,2] - C[2] >= 0)
  square_dist = (x[,1] - C[1])^2 + (x[,2] - C[2])^2
  
  S_cone = function(j) { #compute the averaged score for the candidate hypothesis in cone j
    cone_ind = which(theta >= 2*pi/d*(j - 1) & theta < 2*pi/d*j & !m_set)
    if (length(cone_ind) > 0){
      cand_ind =  order(square_dist[cone_ind]) %>%
        {.[1:round(delta*length(cone_ind))]} %>%
        cone_ind[.]
      s = mean(S[cand_ind])
    } else {
      s = -Inf; cand_ind = NA
    }
    return(list(s = s, cand_ind = cand_ind))
  } 
  
  eval_cones = lapply(1:d, S_cone)
  if (random) {
    random_ind = which(sapply(eval_cones, function(x) x$s) > -Inf)
    inclu_ind = eval_cones[[sample(random_ind,1)]]$cand_ind
  } else {
    inclu_ind = which.max(sapply(eval_cones, function(x) x$s)) %>% 
      eval_cones[[.]] %>% .$cand_ind
  }
  m_set[inclu_ind] = TRUE
  return(m_set)
}


# update M_k under a tree structure 
# S:        non-null likelihood score
# tree_obj: the tree structure where hypotheses lives
# m_set:    vector of indicators for whether a hypothesis is in the previous M_k
inclu_search_tree = function(S, structure, tree_obj, m_set) {
  if (structure %in% c("tree-wide", "tree-decrease")) {
    tree_obj$Set(cand = !m_set, traversal = "level")
    tree_obj$Set(evidence = S, traversal = "level")
    tree_obj$Set(active_cand = rep(FALSE, length(S)), traversal = "level")
    
    tree_obj$Do(function(node) node$active_cand =
                  ifelse(node$isRoot, FALSE, !node$parent$cand),
                filterFun = function(x) x$cand)
    evidence_cand <- tree_obj$Get('evidence', filterFun = function(x) {x$active_cand}, traversal = "level")
    if (length(evidence_cand) > 100) {
      threshold = quantile(evidence_cand, probs = c(0.9))
    } else {threshold = max(evidence_cand)}
    tree_obj$Do(function(node) {node$cand <- FALSE},
                filterFun = function(x) {x$active_cand & x$evidence >= threshold})
  } else if (structure == "tree-increase") {
    tree_obj$Set(cand = !m_set, traversal = "level")
    tree_obj$Set(evidence = S, traversal = "level")
    tree_obj$Set(active_cand = rep(FALSE, length(S)), traversal = "level")
    
    tree_obj$Do(function(node) node$active_cand =
                  ifelse(node$isLeaf, TRUE, all(!node$Get("cand")[-1])),
                filterFun = function(x) x$cand)
    evidence_cand <- tree_obj$Get('evidence', filterFun = function(x) {x$active_cand}, traversal = "level")
    if (length(evidence_cand) > 100) {
      threshold = quantile(evidence_cand, probs = c(0.9))
    } else {threshold = max(evidence_cand)}
    tree_obj$Do(function(node) {node$cand <- FALSE},
                filterFun = function(x) {x$active_cand & x$evidence >= threshold})
  }
  
  return(!tree_obj$Get('cand', traversal = "level"))
}


null_HC = function(w0, sample_size = 1000){
  n = length(w0)
  H_sample = foreach(i = 1:sample_size, .combine = c,
          .packages = c("foreach", "wHC")) %dopar% {
            print(i)
            P = runif(n)
            pwval=cal_cdf(P, w=w0)
            return(hc_cal(pwval,t0ratio =0.4))
          }
  return(H_sample)
}


