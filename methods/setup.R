source("mask.R")
source("interactive.R")
source("online_methods.R")
source("upper_bound.R")
source("sim_dat.R")
source("find_mu.R")
source("single_experiment.R")

suppressPackageStartupMessages({
  require(stats); require(graphics)
  library(doParallel)
  library(parallel)
  
  library(robustbase)
  library(splines)
  library(zoo)
  library(pracma)
  
  library(igraph)
  library(isotone)
  library(data.tree)
  
  library(rootSolve)
  
  library(devtools)
  library(AWFisher)
  library(wHC)
  
  library(ggplot2)
  library(reshape2)
  library(dplyr)
})