source("setup.R")
source("input.R", local = TRUE)
args <- commandArgs(trailingOnly = TRUE)

mode = args[1]
for(single_mode in strsplit(mode, split = ";")[[1]]){
  para_vary = strsplit(single_mode, split = ":")[[1]]
  if (para_vary[1] %in% c("structure", "ub_type", "mu_type")){
    assign(para_vary[1], para_vary[2])
  } else if (para_vary[1] %in% c("smoothed")){
    assign(para_vary[1], as.logical(para_vary[2]))
  } else {
    assign(para_vary[1], as.numeric(para_vary[2]))
  }
}
mu_1 = as.integer(args[2])*c
index = as.integer(args[3])
path = args[4]


if (structure == "grid") {
  n = D^2
  if ("weighted-HC" %in% methods) {
    guess_nonnull = nonnull_cluster(D = D, r = r, C = c(D/2, D/2))
    w0 = as.vector(1*guess_nonnull + 0.5*(!guess_nonnull)); w0 = w0/mean(w0)
    HC_dist = null_HC(w0)
    save(HC_dist, file="result/HC_dist.Rdata")
  }
} else if (structure == "tree-wide") {
  n = sum(20*n_child^(0:(levels-1))) + 1
} else if (structure %in% c("tree-decrease", "tree-increase")) {
  n = sum(n_child^(0:(levels)))
}
if(ub_type == "constant") {
  ub = rep(log(1/alpha), n)
  ub_eps = rep(log(1/alpha), n)
} else if (ub_type == "linear") {
  ub = ub_func_linear_normal(1:n, alpha = alpha, m = n/4) #n/25, 10 for tree
  ub_eps = ub_func_linear_normal(1:n, alpha = alpha, m = n/4)
    #ub_exp_para(m = n/4, alpha = alpha, c = 2) 
}


  if (structure == "sequence_beta") {
    P = sim_seq(n = n, mu_1 = mu_1, mu_0 = mu_0,
                pos_type = "beta", n_1 = n_1, beta_para = beta_para)
    dat = data.frame(P = P, x = 1:n); structure_para = NULL
  } else if (structure == "sequence_gather") {
    P = sim_seq(n = n, mu_1 = mu_1, mu_0 = mu_0,
                pos_type = "gather", n_1 = n_1, sparsity = sparsity)
    dat = data.frame(P = P, x = 1:n); structure_para = NULL
  } else if (structure == "grid") {
    nonnull_ind = nonnull_cluster(D = D, r = r, C = C)
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
  if ("AMT-online" %in% methods) {
    rejections["AMT-online"] = 
      adaptive_mt_online(P = dat$P, alpha = alpha,
                         ub_func = function(k, alpha) ub_func_linear_normal(k, alpha, m = n/25),
                         n_included = 0, h_sum = 0, threshold = threshold,
                         mask_fun = "tent", mask_para = 0.5)$rej
  }
  if ("IMT-tent" %in% methods) {
    rejections["IMT-tent"] = 
      interactive_mt(P = dat$P, x = dat$x, alpha = alpha, ub = ub,
                     structure = structure, structure_para = structure_para,
                     mask_fun = "tent", mask_para = 0.5, S_model = smoothed)
  }
  if ("weighted-HC" %in% methods) {
    load(file = "result/HC_dist.Rdata")
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
                       structure = structure, structure_para = structure_para,
                       mask_fun = "continuous", mask_para = eps, S_model = FALSE)
    }
  }


write.table(data.frame(rejections), file = sprintf("%s/%d_%d.txt", path, mu_1/c, index), 
            row.names = FALSE, col.names = FALSE)

