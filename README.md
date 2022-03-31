# SC Forest for R

R package for "Survey Scale Forests: Estimating Valid Latent Variable Scores from Conditionally Causal Models"

We develop a survey scale forest (SC Forest) algorithm for the estimation of latent variable scores from conditionally causal models with one or more latent variables. SC Forest establishes conditional causality in confirmatory factor analysis (CFA) models with ordinal and/or numerical response variables. Through parametric model restrictions paired with a non-parametric tree-based machine learning approach, SC Forest estimates latent variables scores that fulfill the main criteria for construct validity.

Necessary Arguments:<br>
`data`: Dataset to be analyzed. Note that the observed responses of your model cannot contain missing data. <br>
`model`: Model description in `lavaan` terminology.<br>
`input`: Character vector of partitioning variables. Note that every partitioning variable must either be defineds as `vector` (if it should be treted as categorical) or `numeric`. <br>

Optional Arguments:
`ordered`: Character vector of ordinal observed variables of the model. Default = `NULL.`<br>
`ntrees`: Number of trees to be computed. Default = `100`<br>
`split`: Number of partitioning variales to be selected for random split selection at every split point within a tree (`mtry` argument in `ctree` function). Default = `2`.<br>
`minsize`: Minimum size of teminal nodes of trees. Needs to be big enough to estimate the model parmeters. Default = `300`. <br>
`cutoff_rmsea`: Cutoff value for model fit evaluation. Default = `.05`.<br>
`cutoff_loading`: Cutoff value for multiplicative parameters in model. Note that model with parameter estimates not significantly different form 0 are excluded per default. Default = `0`.<br>
`direct`: No double sampling or bagging used. Variance within forest only comes from random split selection. Default = `FALSE` <br>
`dbsamp`: Double sampling. Default = `TRUE`.<br>
`bagging`: Numeric vector defining the proportion of original data used for tree growing (first number in vector) and re-fitting terminal nodes (second number in vector). Default = `NULL`.<br>
`std.lv`: Standardization of latent variable variances in model. Default = `FALSE`.<br>
`ctree_control`: Control settings for `ctree` function. Default = `ctree_control(minbucket=minsize, mtry= split)`.
