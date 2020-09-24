source("setup.R")


############### Figure 1 ###############
# (a) compare methods in the batch setting with a sequence of hypotheses 
sparsity_seq = c(0.01, seq(0.1, 1, length.out = 10))
para_vary = list(list(name = "methods",
                      value = c("Stouffer-batch", "MST", "AMT-batch", "AMT-online")),
                 list(name = "sparsity", value = sparsity_seq),
                 list(name = "n_1", value = 50),
                 list(name = "mu_1", value = 3),
                 list(name = "threshold", value = 0.05),
                 list(name = "R", value = 500))
result = expr_sparsity(para_vary)
save(result, file=paste(dirname(getwd()), "/results/seq_batch.Rdata", sep = ""))


# (b) compare methods in the online setting with a sequence of hypotheses 
mu_seq = seq(0, 4, length.out = 9)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("lord", "MST", "AMT-online")),
                   list(name = "structure", value = "sequence"),
                   list(name = "prob", value = 0.05),
                   list(name = "mu_1", value = mu),
                   list(name = "threshold_amt", value = "0.05"),
                   list(name = "m", value = 100),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_online(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/seq_online.Rdata", sep = ""))


############### Figure 2 ###############
# minimum alternative mean with targeted type-I and type-II error control
N1_seq <- round(10^(seq(2, 3, length.out = 10)))
N0_seq <- round(10^(seq(2, 5, length.out = 10)))
mu_mat <- matrix(nrow = length(N1_seq), ncol = length(N0_seq))
i = 1; j = 1
for (N1 in N1_seq) {
  for (N0 in N0_seq) {
    mu_mat[i,j] <- find_mu(N1, N0, alpha = 0.05, beta = 0.025, R = 1000)
    print(c(i,j))
    j = j + 1
  }
  i = i + 1; j = 1
}
save(mu_mat, file=paste(dirname(getwd()), "/results/mu_heatmap.Rdata", sep = ""))



############### Figure 4 ###############
# (a) compare methods in a grid of hypotheses with nonnulls in the center 
mu_seq = seq(0, 1.8, length.out = 7)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("Stouffer-batch", "MST", "IMT-tent")),
                   list(name = "structure", value = "grid"),
                   list(name = "mu_1", value = mu),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/grid_center.Rdata", sep = ""))
     

# (b) compare methods in a grid of hypotheses with nonnulls in the corner 
mu_seq = seq(0, 1.8, length.out = 7)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("Stouffer-batch", "MST", "IMT-tent")),
                   list(name = "structure", value = "grid"),
                   list(name = "mu_1", value = mu),
                   list(name = "C", value = c(20, 30)),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/grid_corner.Rdata", sep = ""))




############### Figure 5 ###############
# compare methods in a fixed tree of hypotheses with 801 hypotheses 
mu_seq = seq(0, 4, length.out = 9)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("Stouffer-batch", "MST", "IMT-tent")),
                   list(name = "structure", value = "tree-wide"),
                   list(name = "mu_1", value = mu),
                   list(name = "m", value = 10),
                   list(name = "smoothed", value = FALSE),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/tree_batch.Rdata", sep = ""))



############### Figure 6 ###############
# (a) compare methods in a fixed tree decreasing probability 
mu_seq = seq(0, 3, length.out = 7)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("MST", "IMT-tent")),
                   list(name = "structure", value = "tree-decrease"),
                   list(name = "levels", value = 4),
                   list(name = "n_child", value = 3),
                   list(name = "mu_1", value = mu),
                   list(name = "m", value = 10),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/tree_decrease.Rdata", sep = ""))


# (b) compare methods in a fixed tree increasing probability 
print("tree-increase")
mu_seq = seq(0, 3, length.out = 7)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("MST", "IMT-tent")),
                   list(name = "structure", value = "tree-increase"),
                   list(name = "levels", value = 4),
                   list(name = "n_child", value = 3),
                   list(name = "mu_1", value = mu),
                   list(name = "m", value = 10),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/tree_increase.Rdata", sep = ""))



############### Figure 7 ###############
# online interactive
mu_seq = seq(0, 4, length.out = 9)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("lord", "MST", "AMT-online", "IMT-online")),
                   list(name = "structure", value = "seq-block"),
                   list(name = "prob", value = 0.05),
                   list(name = "mu_1", value = mu),
                   list(name = "threshold_amt", value = 0.05),
                   list(name = "threshold_imt", value = 0.05),
                   list(name = "m", value = 100),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_online(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/seq_online_block.Rdata", sep = ""))


############### Figure 8 ###############
# compare methods in a growing tree of hypotheses
mu_seq = seq(0, 4, length.out = 9)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("lord", "MST", "AMT-online", "IMT-online")),
                   list(name = "structure", value = "tree-online"),
                   list(name = "mu_1", value = mu),
                   list(name = "threshold_amt", value = 0.05),
                   list(name = "threshold_imt", value = 0.6),
                   list(name = "m", value = 100),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_online(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/tree_online.Rdata", sep = ""))



############### Figure 9 ###############
# (a) compare methods in the batch setting with a sequence of hypotheses 
n_1 = 100; mu_1 = 1.5
mu_0_seq = seq(0, -4, length.out = 9)
result = list()
for (mu_0 in mu_0_seq) {
  para_vary = list(list(name = "methods",
                        value = c("Stouffer-batch", "MST", "AMT-batch", "AMT-railway")),
                   list(name = "structure", value = "sequence_even"),
                   list(name = "n", value = 1000),
                   list(name = "n_1", value = n_1),
                   list(name = "mu_1", value = mu_1),
                   list(name = "mu_0", value = mu_0),
                   list(name = "R", value = 500))
  result[[as.character(mu_0)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/conservative_", n_1, "_", mu_1, ".Rdata", sep = ""))



############### Figure 11 ###############
mu_seq = seq(0, 1.8, length.out = 7)
eps = seq(0, 0.8, 0.2)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("IMT-tent", paste("IMT-",eps*10, sep = ""))),
                   list(name = "structure", value = "grid"),
                   list(name = "mu_1", value = mu),
                   list(name = "smoothed", value = FALSE),
                   list(name = "ub_type", value = "constant"),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()),
                        "/results/continuous_masking_ubConstant.Rdata", sep = ""))


############### Figure 12 ###############
sparsity_seq = c(0.01, seq(0.1, 1, length.out = 10))
para_vary = list(list(name = "methods",
                      value = c("MST (m = n/4)", "MST (m = n/2)", "MST (m = 3n/4)", "MST (m = l, oracle)")),
                 list(name = "sparsity", value = sparsity_seq),
                 list(name = "R", value = 500))
result = expr_sparsity(para_vary)
save(result, file=paste(dirname(getwd()), "/results/mst_linear_bound.Rdata", sep = ""))


############### Figure 13 ###############
sparsity_seq = c(0.01, seq(0.1, 1, length.out = 10))
para_vary = list(list(name = "methods",
                      value = c("MST (linear)", "MST (poly)", "MST (discrete)", "MST (invert)")),
                 list(name = "ub_type", value = "curve"),
                 list(name = "sparsity", value = sparsity_seq),
                 list(name = "R", value = 500))
result = expr_sparsity(para_vary)
save(result, file=paste(dirname(getwd()), "/results/mst_curve_bound.Rdata", sep = ""))

############### Figure 14 ###############
sparsity_seq = c(0.01, seq(0.1, 1, length.out = 10))
para_vary = list(list(name = "methods",
                      value = c("MFT (m = n/4)", "MFT (m = n/2)", "MFT (m = 3n/4)", "MFT (m = l, oracle)")),
                 list(name = "sparsity", value = sparsity_seq),
                 list(name = "R", value = 500))
result = expr_sparsity(para_vary)
save(result, file=paste(dirname(getwd()), "/results/mft_linear_bound.Rdata", sep = ""))


############### Figure 15 ###############
sparsity_seq = c(0.01, seq(0.1, 1, length.out = 10))
para_vary = list(list(name = "methods",
                      value = c("MFT (linear)", "MFT (curve)")),
                 list(name = "ub_type", value = "curve"),
                 list(name = "sparsity", value = sparsity_seq),
                 list(name = "R", value = 500))
result = expr_sparsity(para_vary)
save(result, file=paste(dirname(getwd()), "/results/mft_curve_bound.Rdata", sep = ""))


############### Figure 16 ###############
# compare with adaptive weighted Fisher and weighted HC
mu_seq = seq(0, 1.8, length.out = 1)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("Stouffer-batch", "MST", "IMT-tent", "IMT-Gaussians",
                                  "AW-Fisher", "weighted-HC")),
                   list(name = "structure", value = "grid"),
                   list(name = "mu_1", value = mu),
                   list(name = "D", value = 10),
                   list(name = "r", value = 5), 
                   list(name = "C", value = c(5, 5)),
                   list(name = "m", value = 25),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()), "/results/grid_center_extra.Rdata", sep = ""))



###############Rebuttal Figure 1(b)
# heterogeneous non-nulls
mu_seq = seq(0, 1.8, length.out = 7)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("Stouffer-batch", "MST", "IMT-tent")),
                   list(name = "structure", value = "grid"),
                   list(name = "mu_type", value = "fading"),
                   list(name = "r", value = 70),
                   list(name = "mu_1", value = mu),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file="results/grid_center_heteroMu.Rdata")  


############### Rebuttal Figure 3 ###############
#(b)
mu_seq = seq(0, 1.8, length.out = 7)
eps = seq(0, 0.8, 0.2)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("IMT-tent", paste("IMT-",eps*10, sep = ""))),
                   list(name = "structure", value = "grid"),
                   list(name = "mu_1", value = mu),
                   list(name = "smoothed", value = FALSE),
                   list(name = "ub_type", value = "linear"),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()),
                        "/results/continuous_masking_ubLinear_zeroMean.Rdata", sep = ""))

#(c)
mu_seq = seq(0, 1.8, length.out = 7)
eps = seq(0, 0.8, 0.2)
result = list()
for (mu in mu_seq) {
  para_vary = list(list(name = "methods",
                        value = c("IMT-tent", paste("IMT-",eps*10, sep = ""))),
                   list(name = "structure", value = "grid"),
                   list(name = "mu_1", value = mu),
                   list(name = "smoothed", value = FALSE),
                   list(name = "ub_type", value = "linear"),
                   list(name = "centered", value = TRUE),
                   list(name = "R", value = 500))
  result[[as.character(mu)]] = expr_offline(para_vary)
}
save(result, file=paste(dirname(getwd()),
                        "/results/continuous_masking_ubLinear.Rdata", sep = ""))




