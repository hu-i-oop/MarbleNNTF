# MarbleNNTF
A package to implement non-negative tensor factorization for EHR phenotyping in R, according to the Marble algorithm by Joyce Ho et al., (2014).
(https://doi.org/10.1145/2623330.2623658)


# Directions for End-Users

### Installing the Package

• Ensure that `git` is installed in your device

• Run the command `git clone https://github.com/hu-i-oop/MarbleNNTF/` in a terminal.

• In R or Rstudio, set the working directory to the cloned folder using the `setwd(...)` function.

• Type `source("MarbleAlgorithm.R")` to load the function.

Then, you may call the `Marble(...)` function in R to decompose X into a bias tensor and interaction tensor.

To do so, you must have: (1) An array/tensor object in R, named `X`, (2) the number of ranks `R` of the outputted interaction tensor, and (3) Values set for the `alpha` and `gamma` hyperparameters.

### How to Use the Package

For example, for a Rank-3 Tensor X: `Marble(X = X, R = 5, alpha = 0.01, gamma = c(0.1, 0.1, 0.1))`

This function outputs a list of arrray objects: the bias tensor "C" and rank-(R = 5) interaction tensors "V", respectively.
