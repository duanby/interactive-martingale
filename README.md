# Interactive martingale test
This code is an accompaniment to the paper titled "Interactive Martingale Tests for the Global Null", and produce all plots in the paper.

## Overview
The interactive martingale test controls the global Type 1 error in multiple testing. It achieves high power against structured alterantives by allowing human interaction to incorporates prior knowledge, covariate information, non-null structures, in a highly flexbile manner. We also propose martingale variants of classical tests such as the Fisher test and the Stouffer test, which can have high power if the hypotheses have a good prior ordering (mostly non-nulls at front). 

To reproduce the results:
```
$Rscript reproduce.R
```
To generate the figures:
```
$Rscript plots.R
```
Reproducing the results/ folder is time-consuming, but plotting the results to generate the figures/ folder is quick.

## Functions and files
methods/interactive.R: interactively ordered martingale test; martingale Stouffer test; and other comparing methods.

methods/online_methods.R: online version of the interactive test and other comparing methods.

methods/upper_bound.R: martingale boundaries to decide rejection.

methods/mask.R: several forms of masking p-values (i.e. decompose a p-value into two parts and hide one part from the analyst).

methods/sim_dat.R: simulate p-values of various structures.

results/.: rejection decision of each method in each experiment.

figures/.: figures in the paper.

## Dependencies
The code was tested using R (version 3.6.3) and the following packages:
stats, graphics, robustbase, splines, zoo, pracma, igraph, isotone, data.tree, rootSolve, devtools, AWFisher, wHC, ggplot2, reshape2, dpylr
