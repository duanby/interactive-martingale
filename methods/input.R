## problem setup
alpha = 0.05
mask_fun = "tent"
mask_para = 0.5
centered = FALSE

## a sequence of p-values
n = 10^4
n_1 = 100
mu_1 = 1.5
mu_0 = 0
sparsity = 0.2
beta_para = 20
prob = 0.01
n_unit = 10^4

## a grid of p-values
D = 100
r = 50
mu_1 = 1
mu_0 = 0
C = c(D/2, D/2)
d = 5; delta = 0.1
rho = 0
mu_type = "constant"

## a tree of p-values
levels = 4
n_child = 3
n_nonnull = 7
prob_ini = c(rep(0.9,5), rep(0.1,30), rep(0.9,5))
child_fac = c(0, 0.8, 1)
decrease = TRUE

## methods para
threshold_amt = 0.05
smoothed = TRUE
ub_type = "linear"
m = n/4
methods = c("Stouffer-batch", "MST", "IMT-tent")
