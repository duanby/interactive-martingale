# simulate a sequence of p-values
# n:         number of hypotheses
# mu_1:      non-null mean value
# mu_0:      null mean value
# pos_type:  type of position of non-nulls
# n_1:       number of non-nulls if pos_type = gather
# sparsity:  parameter encoding the sparisty of the non-nulls if pos_type = gather
# beta_para: parameter encoding the position of non-nulls if pos_type = beta
# prob:      probability of being non-null if pos_type = even
sim_seq = function(n, mu_1, mu_0, pos_type, n_1 = NA,
                   sparsity = NA, beta_para = NA, prob = NA){
  if (pos_type == "gather") {
    mu_set = rep(mu_0, n); mu_set[sample(sparsity*n,n_1)] = mu_1
  } else if (pos_type == "even") {
    if (!is.na(prob)) {
      nonnull_ind <- rbinom(n, 1, prob) 
    } else {
      nonnull_ind <- sample(n, n_1)
    }
    mu_set = rep(mu_0, n); mu_set[nonnull_ind == 1] = mu_1
  } else if (pos_type == "beta") {
    mu_set = rep(mu_0, n); mu_set[unique(round(rbeta(n_1, 2, beta_para)*n))] = mu_1
  } else if (pos_type == "block") {
    if(is.na(n_1)) {n_1 = n*prob}
    start_pos = sample(1:(n-n_1 + 1), 1)
    mu_set = rep(mu_0, n); mu_set[start_pos:(start_pos + n_1 - 1)] = mu_1
  } else if (pos_type == "gather_vary") {
    dat_null = rnorm(n - n_1); dat_nonnull = rnorm(n_1) + mu_1
    dat_set = matrix(nrow = length(sparsity), ncol = n)
    for(i in 1:length(sparsity)){
      nonnull_ind = sample(sparsity[i]*n,n_1)
      dat_set[i,nonnull_ind] = dat_nonnull; dat_set[i,-nonnull_ind] = dat_null
    }
  }
  
  if (pos_type != "gather_vary") {
    dat_set = rnorm(n) + mu_set
  } 
  
  P = 1 - pnorm(dat_set)
  return(P)
}




# generate a grid of non-null indicators with non-nulls clustered.
# D: size of grid D*D
# d: size of non-null cluster (= radius^2)
# C: coordinates of the center of the cluster
nonnull_cluster = function(D, r, C, mu_type = "constant"){
  nonnull_ind = matrix(FALSE, ncol = D, nrow = D)
  bd_x <- floor(sqrt(r))
  if (mu_type == "fading") {
    for (i in (C[1] - bd_x):(C[1] + bd_x)) {
      bd_y <- floor(sqrt(r - (C[1] - i)^2))
      for (j in (C[2] - bd_y):(C[2] + bd_y)) {
        nonnull_ind[i,j] = 1 - ((i - C[1])^2 + (j - C[2])^2)/r^2
      }
    }
  } else {
    for (i in (C[1] - bd_x):(C[1] + bd_x)) {
      bd_y <- floor(sqrt(r - (C[1] - i)^2))
      for (j in (C[2] - bd_y):(C[2] + bd_y)) {
        nonnull_ind[i,j] = TRUE
      }
    }
  }
  return(nonnull_ind)
}


# simulate a vector p-values and side information for cluster structure.
# nonnull_ind: a grid of non-null indicators
# mu_1:        non-null mean value
# mu_0:        null mean value
# rho:         correlations among data
sim_cluster = function(nonnull_ind, mu_1, mu_0, rho, mu_type = "constant"){
  n_col = ncol(nonnull_ind); n_row = nrow(nonnull_ind)
  x <- cbind(rep(1:n_col, each = n_row), rep(1:n_row, times = n_col))
  
  if (mu_type == "constant"){
    mu <- as.vector(nonnull_ind)*mu_1 + as.vector(!nonnull_ind)*mu_0
  } else if (mu_type == "varying") {
    mu_1_seq = mu_1*(abs(sin(1:n_col*n_row)) + 1/2)
    mu <- as.vector(nonnull_ind)*mu_1_seq + as.vector(!nonnull_ind)*mu_0
  } else if (mu_type == "fading") {
    mu <- as.vector(nonnull_ind)*mu_1 + as.vector(nonnull_ind == 0)*mu_0
  }
  
  if(rho == 0){
    dat <- rnorm(n_col*n_row) + mu
  } else {
    n = n_col*n_row
    Sigma = matrix(rho, nrow = n, ncol = n); diag(Sigma) = 1
    dat = mvrnorm(n = 1, mu = mu, Sigma = Sigma)
  }
  P <- 1 - pnorm(dat)
  return(list(Z = dat, P = as.vector(P), x = x))
}




# generate a tree structure.
# levels:  number of generations in the tree
# n_child: number of children for each parent node (from the third generation)
tree_structure = function(levels, n_child, parent = NULL, tree_type = "wide"){
  if (is.null(parent)) {
    parent <- Node$new("1")
    if(tree_type == "wide") {
      for (i in 1:20) child <- parent$AddChild(i)  ##20 children on the first level
      if (levels > 1) for (i in 1:20) tree_structure(n_child = n_child, levels = levels - 1, 
                                                     parent = parent$children[[i]])
    }
  }
  for (i in 1:n_child) child <- parent$AddChild(i)
  if (levels > 1) for (i in 1:n_child) tree_structure(n_child = n_child, levels = levels - 1, 
                                                      parent = parent$children[[i]])
  return(parent)
}


# generate the side information for a tree structure: edge matrix
# levels:      number of generations in the tree
# n_child_seq: a vector of number of children for the parent on each level
edge_gen = function(levels, n_child_seq) {
  edge_mat = matrix(ncol = 2, nrow = 0)
  children_count = cumprod(n_child_seq)
  node_count = cumsum(c(1, children_count))
  for (level in 1:levels) {
    temp_mat = matrix(ncol = 2, nrow = children_count[level])
    if (level == 1){
      temp_mat[,1] = rep(1:node_count[level],
                         each = n_child_seq[level])
    } else {
      temp_mat[,1] = rep((node_count[level - 1] + 1):node_count[level],
                         each = n_child_seq[level])
    }
    temp_mat[,2] = (node_count[level] + 1):(node_count[level + 1])
    edge_mat = rbind(edge_mat, temp_mat)
  }
  return(edge_mat)
}

# simulate a tree of p-values with a subtree of non-nulls.
# tree-frame: a tree object generated by tree_structure
# n_1:        number of non-nulls
# mu_1:       non-null mean value
# mu_0:       null mean value
sim_tree = function(tree_frame, n_1, mu_1, mu_0, decrease, tree_type){
  n_node <- tree_frame$totalCount
  #insert non-nulls by a depth-first ordering
  tree_frame$Set(nonnull_ind = c(rep(1, n_1), rep(0, n_node - n_1)),
                 traversal = ifelse(decrease, "pre-order", "post-order")) 
  Sort(tree_frame, attribute = "nonnull_ind", decreasing = FALSE) 
  #arrange non-nulls at later-visited node
  tree_frame$Do(function(node) {
    mean_val = node$nonnull_ind*mu_1 + (1 - node$nonnull_ind)*mu_0
    node$pvalue <- 1 - pnorm(rnorm(1) + mean_val)
  })
  
  P = tree_frame$Get('pvalue', traversal = "level")
  levels = tree_frame$height - 1
  if (tree_type == "wide") {
    edge_mat = edge_gen(levels = levels, n_child_seq = c(20, rep(3, levels - 1))) #assume each node have three children!!!!!
  } else if (tree_type == "small") {
    edge_mat = edge_gen(levels = levels, n_child_seq = rep(length(tree_frame$children), levels))
  }
  if (decrease) {edge_mat = edge_mat[,c(2,1)]}
  return(list(tree_obj = tree_frame, P = P, x = edge_mat))
}

# simulate p-values for the new generation of a growing tree
# parent_probs: a vector of probability of being non-null of the nodes on the previous level
# mu_1:         nonnull mean value
# mu_0:         null mean value
# n_child:      number of children from each node
# child_fac:    factor to decrease the parent_prob for each chil node
sim_online_tree = function(parent_probs, mu_1, mu_0, n_child, child_fac){
  n = length(parent_probs)
  child_probs = rep(parent_probs, each = n_child)*rep(child_fac, n)
  nonnull_ind <- sapply(child_probs, function(x) rbinom(1,1,x)) 
  mu_set = nonnull_ind*mu_1 + (1 - nonnull_ind)*mu_0
  
  dat_set <- rnorm(3*n) + mu_set
  P <- 1 - pnorm(dat_set)
  return(list(P = P, x = child_probs))
}



