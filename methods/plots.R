source("setup.R")
rowSes = function(m) {apply(m, 1, sd)/sqrt(ncol(m))}
allcolor_name = c("deepskyblue2", "purple", "green3", "yellow4", "red1", "magenta", "tan1")


figure_generator = function(mode, legend_name, color_name, expr_type = "offline", 
                           alpha = 0.05, alpha_bar = TRUE, error_bar = FALSE,
                           exclude_methods = c(), legend_pos = c(0.8,0.7), legend_size = 8,
                           y_scale = NA, x_factor = NA, x_name = NA,
                           save = FALSE){
  load(file = paste(dirname(getwd()),"/results/", mode,".Rdata", sep = ""))
  load(file = paste(dirname(getwd()),"/results/", mode,".Rdata", sep = ""))
  if(expr_type == "offline") {
    if(is.data.frame(result[[1]][[1]])){
      power = sapply(result, function(x) {
        rowMeans(matrix(unlist(x), nrow = nrow(x[[1]])),na.rm = TRUE)})
    } else {
      power = sapply(result, function(x) {
        rowMeans(matrix(unlist(x), nrow = length(x[[1]])),na.rm = TRUE)})
    }
    sd_power = sapply(result, function(x) {rowSes(matrix(unlist(x), nrow = length(x[[1]])))})
    mu_seq = rep(as.numeric(colnames(power)), each = nrow(power)); neg = FALSE
    if(any(mu_seq < 0)) {
      mu_seq = -mu_seq; neg = TRUE
    }
    if(!is.na(x_factor)) {mu_seq = x_factor*mu_seq}
    if(is.na(x_name)) {x_name = "alternative mean"}
    y_name = "power"
  } else if(expr_type == "online") {
    power = sapply(result, function(x) {
      rowMeans(matrix(unlist(x)/10^3, nrow = length(x[[1]])),na.rm = TRUE)})
    sd_power = sapply(result, function(x) {rowSes(matrix(unlist(x)/10^3, nrow = length(x[[1]])))})
    mu_seq = rep(as.numeric(colnames(power)), each = nrow(power))
    x_name = "alternative mean"; y_name = expression(paste("detection time/", 10^{3}))
  } else if (expr_type == "sparsity") {
    sparsity_seq = c(0.01, seq(0.1, 1, length.out = 10))
    power = matrix(nrow = length(result[[1]][[1]]), ncol = length(result[[1]])); sd_power = power
    for(i in 1:length(result[[1]])){
      power[,i] = rowMeans(sapply(result, function(x) {x[[i]]}))
      sd_power[,i] = rowSes(sapply(result, function(x) {x[[i]]}))
    }
    mu_seq = rep(sparsity_seq, each = nrow(power))
    x_name = "sparsity"; y_name = "power"
  }
  
  df_power = data.frame(mu_seq = mu_seq,
                        power = as.vector(power),
                        sd = as.vector(sd_power),
                        grp = rep(legend_name,
                                  ncol(power)))
  p = ggplot(data = subset(df_power, !(grp %in% exclude_methods)),
             aes(x = mu_seq, y = power, group = grp, fill = grp)) +
    geom_line(aes(linetype = grp, color = grp), size = 0.8) +
    geom_point(aes(shape = grp, color = grp), size = 2.5) +
    scale_color_manual(values = c(color_name)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = legend_pos, legend.text = element_text(size = legend_size)) +
    xlab(x_name) + ylab(y_name) 
  if(expr_type %in% c("offline", "sparsity") & is.na(y_scale)) {
    p = p + scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))
  } 
  if(alpha_bar) {
    p = p + geom_hline(yintercept = alpha)
  }
  if(error_bar) {
    p = p + geom_errorbar(aes(ymin=power-sd, ymax=power+sd), width=.1, 
                          position=position_dodge(0.05)) 
  }
  if(!is.na(y_scale)){
    p = p + ylab("Type 1 error") +
      scale_y_continuous(breaks = seq(0, y_scale, length.out = 6), limits = c(0,y_scale))
  }
  if(neg){
    p = p + scale_x_continuous(labels = 0:-4)
  }
  plot(p)
  if(save){
    ggsave(filename = paste(dirname(getwd()),"/figures/", mode, ".pdf", sep = ""),
           plot = p, device = "pdf", width = 4, height = 3.6)
  }
}

bound_figure_generator = function(mode, bound, legend_name, color_name,
                                  legend_pos = c(0.8,0.2), save = FALSE, add_point = TRUE,
                                  x_name = "k", y_name = "bound", width = NA, height = 3.6){
  n = nrow(bound)
  legend_name = factor(legend_name, levels = legend_name)
  if (add_point) {k = rep(1:n, ncol(bound))} else {k = rep((1:n)/(n+1), ncol(bound))}
  ub_df = data.frame(k = k, bound = as.vector(bound),
                     type = rep(legend_name, each = n))
  p = ggplot(ub_df, aes(x=k, y=bound, group=type, fill = type)) +
    geom_line(aes(linetype = type, color = type), size = 0.8) +
    scale_color_manual(values = color_name) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = legend_pos, legend.text = element_text(size = 8)) + 
    xlab(x_name) + ylab(y_name) 
  if (add_point) {
    p = p + geom_point(data = dplyr::filter(ub_df, k == 2000),
                       aes(shape = type, color = type), size = 2.5) 
  }
  if (!is.na(width)){
    p = p + theme(legend.key.width=unit(width,"cm"))
  }
  plot(p)
  if(save){
    ggsave(filename = paste(dirname(getwd()),"/figures/", mode, ".pdf", sep = ""),
           plot = p, device = "pdf", width = 4, height = height)
  }
}

######Figure 1 
#(a)
methods = c("Stouffer (batch)", "MST", "AMT (batch)", "AMT (online)")
figure_generator(mode = "seq_batch", expr_type = "sparsity",
                 legend_name = factor(methods, levels = methods[c(3,4,2,1)]),
                 color_name = allcolor_name[2:5], alpha_bar = FALSE)
#(b)
methods = c("Bonferroni", "MST", "AMT (online)")
figure_generator(mode = "seq_online", expr_type = "online",
                 legend_name = factor(methods, levels = methods[c(3,2,1)]),
                 color_name = allcolor_name[3:5], alpha_bar = FALSE)

#######Figure 2
load(file = paste(dirname(getwd()),"/results/mu_heatmap.Rdata", sep = ""))
longData = melt(mu_mat)
N1_seq <- round(10^(seq(2, 3, length.out = 10)))
N0_seq <- round(10^(seq(2, 5, length.out = 10)))
p = ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="yellow", high= "red", limits=c(0, 4), name = expression(mu)) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     #axis.title.x=element_blank(), axis.title.y=element_blank(),
                     plot.title=element_text(size=11),
                     legend.position = "right",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  geom_vline(xintercept=seq(-0.5, 31.5, by=1), color = "grey") +
  geom_hline(yintercept=seq(-0.5, 31.5, by=1), color = "grey") +
  xlab(expression(log[10](N[0]))) + ylab(expression(log[10](N[1]))) +
  scale_x_continuous(labels = round(log10(N0_seq[seq(2, 10, 2)]),1), breaks = seq(2, 10, 2)) +
  scale_y_continuous(labels = round(log10(N1_seq[seq(2, 10, 2)]),1), breaks = seq(2, 10, 2)) +
  coord_cartesian(xlim = c(1, 10), ylim = c(1, 10)) 
plot(p)
ggsave(filename = paste(dirname(getwd()),"/figures/mu_heatmap.pdf", sep = ""),
       plot = p, device = "pdf", width = 5, height = 4)



#######Figure 3
nonnull_ind = nonnull_cluster(D = 100, r = 30, C = c(30,30))
dat = sim_cluster(nonnull_ind, mu_1 = 2, mu_0 = 0, rho = 0); n = length(dat$P)
structure_para = list(d = 10, delta = 0.1)
ub = ub_func_linear_normal(1:n, alpha = 0.05, m = n/4)
example_expr = interactive_mt(P = dat$P, x = dat$x, alpha = 0.05, ub = ub,
                 structure = "grid", structure_para = structure_para,
                 mask_fun = "tent", mask_para = 0.5, S_model = TRUE, figure = TRUE)


######Figure 4
methods = c("Stouffer (batch)", "MST", "IMT")
#(a)
figure_generator(mode = "grid_center",
                 legend_name = factor(methods, levels = methods[c(3,2,1)]),
                 legend_pos = c(0.8, 0.2),
                 color_name = allcolor_name[c(1,4,5)])
#(b)
figure_generator(mode = "grid_corner",
                 legend_name = factor(methods, levels = methods[c(3,2,1)]),
                 legend_pos = c(0.8, 0.2),
                 color_name = allcolor_name[c(1,4,5)])

######Figure 5
methods = c("Stouffer (batch)", "MST", "IMT")
figure_generator(mode = "tree_batch", 
                 legend_name = factor(methods, levels = methods[c(3,2,1)]),
                 legend_pos = c(0.2, 0.8), x_factor = 0.5,
                 color_name = allcolor_name[c(1,4,5)])

######Figure 6
methods = c("Stouffer (batch)", "MST", "IMT")
#(a)
figure_generator(mode = "tree_decrease", 
                 legend_name = factor(methods, levels = methods[c(3,2,1)]),
                 exclude_methods = methods[1], legend_pos = c(0.2, 0.8), x_factor = 0.5,
                 color_name = allcolor_name[c(1,4)])
#(b)
figure_generator(mode = "tree_increase", 
                 legend_name = factor(methods, levels = methods[c(3,2,1)]),
                 exclude_methods = methods[1], legend_pos = c(0.2, 0.8), x_factor = 0.5,
                 color_name = allcolor_name[c(1,4)])

######Figure 7
methods = c("Bonferroni", "MST", "AMT (online)", "IMT(online)")
figure_generator(mode = "seq_online_block", expr_type = "online",
                 legend_name = factor(methods, levels = methods[c(4,3,2,1)]),
                 legend_pos = c(0.8, 0.8),
                 color_name = allcolor_name[c(1,3,4,5)], error_bar = TRUE, alpha_bar = FALSE)

######Figure 8
methods = c("Bonferroni", "MST", "AMT (online)", "IMT(online)")
figure_generator(mode = "tree_online", expr_type = "online",
                 legend_name = factor(methods, levels = methods[c(4,3,2,1)]),
                 legend_pos = c(0.8, 0.8),
                 color_name = allcolor_name[c(1,3,4,5)], alpha_bar = FALSE)

######Figure 9
methods = c("Stouffer (batch)", "MST", "IMT (tent)", "IMT (railway)")
figure_generator(mode = "conservative_100_1.5", x_name = "null mean",
                 legend_name = factor(methods, levels = methods[c(4,3,2,1)]),
                 legend_pos = c(0.75, 0.6), legend_size = 12, 
                 color_name = allcolor_name[c(2,1,4,5)])

######Figure 10
grid_P = seq(0.001, 0.999, length.out = 999); eps = c(seq(0.2, 0.8, 0.2), 0)
h_g = lapply(eps, function(x){
  mask_gen(grid_P, mask_fun = "continuous", mask_para = x)})
h_g = c(list(mask_gen(grid_P, mask_fun = "tent", mask_para = 0.5)), h_g)
#(a)
h_mat = sapply(h_g, function(x){x$h})
legend_name = c("h", paste("log(f_",eps[-5],")", sep = ""), "log(f_m)")
bound_figure_generator(mode = "h_varyc", bound = h_mat, legend_name = legend_name,
                       color_name = allcolor_name[1:6], add_point = FALSE, width = 1.5, height = 3,
                       legend_pos = c(0.8,0.7), x_name = "p-value", y_name = "missing bits")
#(b)
g_mat = sapply(h_g, function(x){x$g})
legend_name = c("g_ori", paste("g_",eps[-5], sep = ""), "g_m")
bound_figure_generator(mode = "g_varyc", bound = g_mat, legend_name = legend_name,
                       color_name = allcolor_name[1:6], add_point = FALSE, width = 1.5, height = 3,
                       legend_pos = c(0.8,0.7), x_name = "p-value", y_name = "masked p-values")


######Figure 11
eps = seq(0, 0.8, 0.2)
methods = c("h", "f_m", paste("f_",eps[-1], sep = ""))
figure_generator(mode = "continuous_masking_ubConstant", 
                 legend_name = factor(methods, levels = methods[c(1, 3:6, 2)]),
                 legend_pos = c(0.2, 0.7), x_factor = 0.3, 
                 color_name = allcolor_name[1:6])

######Figure 12
#(a)
methods = c("m = n/4", "m = n/2", "m = 3n/4", "m = l (oracle)")
figure_generator(mode = "mst_linear_bound", expr_type = "sparsity",
                 legend_name = factor(methods, levels = methods[1:4]),
                 legend_pos = c(0.8, 0.6), alpha_bar = FALSE,
                 color_name = allcolor_name[c(1,5,6,7)])

#(b)
n = 10^4
ub_mat = sapply(1:3, function(x){ub_func_linear_normal(1:n, alpha = 0.05, m = x*n/4)})
legend_name = c("m = n/4", "m = n/2", "m = 3n/4")
bound_figure_generator(mode = "mst_linear_bound_plot", bound = ub_mat, legend_name = legend_name,
                       color_name = allcolor_name[c(1,5,6)])


##############Figure 13
#(a)
methods = c("linear", "curve (poly)", "curve (discrete)", "curve (inverted)")
figure_generator(mode = "mst_curve_bound", expr_type = "sparsity",
                 legend_name = factor(methods, levels = methods[1:4]),
                 legend_pos = c(0.8, 0.6), alpha_bar = FALSE,
                 color_name = allcolor_name[c(1,5,3,7)])

#(b)
n = 10^4
ub_mat = cbind(ub_func_linear_normal(1:n, alpha = alpha, m = n/4), 
               ub_func_curve_poly(1:n, alpha = alpha),
               ub_func_curve_mix(L = n, alpha = alpha, f_mix = f_mix),
               ub_func_curve_invert(1:n, alpha = alpha))
legend_name = c("linear", "curve (poly)", "curve (discrete)", "curve (inverted)")
bound_figure_generator(mode = "mst_curve_bound_plot", bound = ub_mat, legend_name = legend_name,
                       color_name = allcolor_name[c(1,5,3,7)])

###############Figure 14
#(a)
methods = c("m = n/4", "m = n/2", "m = 3n/4", "m = l (oracle)")
figure_generator(mode = "mft_linear_bound", expr_type = "sparsity",
                 legend_name = factor(methods, levels = methods[1:4]),
                 legend_pos = c(0.8, 0.6), alpha_bar = FALSE,
                 color_name = allcolor_name[c(1,5,6,7)])

#(b)
n = 10^4
ub_mat = sapply(1:3, function(x){ub_func_linear_exp(1:n, alpha = alpha, m = x*n/4)})
legend_name = c("m = n/4", "m = n/2", "m = 3n/4")
bound_figure_generator(mode = "mft_linear_bound_plot", bound = ub_mat, legend_name = legend_name,
                       color_name = allcolor_name[c(1,5,6)])


######Figure 15
methods = c("linear", "curve")
#(a)
figure_generator(mode = "mft_curve_bound", expr_type = "sparsity",
                 legend_name = factor(methods, levels = methods[1:2]),
                 legend_pos = c(0.8, 0.6), alpha_bar = FALSE,
                 color_name = allcolor_name[c(1,5)])

#(b)
n = 10^4
ub_mat = cbind(ub_func_linear_exp(1:n, alpha = alpha, m = n/4),
               ub_func_curve_gamma(1:n, alpha = alpha))
legend_name = c("linear", "curve")
bound_figure_generator(mode = "mft_curve_bound_plot", bound = ub_mat, legend_name = legend_name,
                       color_name = allcolor_name[c(1,5)])



###### Figure 16
methods = c("Stouffer-batch", "MST", "IMT", "IMT-Gaussians",
            "AW-Fisher", "weighted-HC")
figure_generator(mode = "grid_center_extra", 
                 legend_name = factor(methods, levels = methods[1:6]),
                 exclude_methods = methods[c(1,2,4)], legend_pos = c(0.2, 0.8), 
                 color_name = allcolor_name[c(1,6,7)])






#########rebuttal Figure 1 (b)
methods = c("Stouffer (batch)", "MST", "IMT")
figure_generator(mode = "grid_center_heteroMu", 
                 legend_name = factor(methods, levels = methods[3:1]),
                 legend_pos = c(0.8, 0.2),
                 color_name = allcolor_name[c(1,4,5)], x_name = "average alternative mean")

######### rebuttal Figure 3
eps = seq(0, 0.8, 0.2)
methods = c("h", "f_m", paste("f_",eps[-1], sep = ""))
#(b)
figure_generator(mode = "continuous_masking_ubLinear", 
                 legend_name = factor(methods, levels = methods[c(1, 3:6, 2)]),
                 legend_pos = c(0.2, 0.7),
                 color_name = allcolor_name[1:6])
#(c)
figure_generator(mode = "continuous_masking_ubLinear_zeroMean", 
                 legend_name = factor(methods, levels = methods[c(1, 3:6, 2)]),
                 legend_pos = c(0.2, 0.7), x_factor = 0.3,
                 color_name = allcolor_name[1:6])





