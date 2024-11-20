library(Matrix)
library(rTensor)

########################
### Marble Algorithm ###
########################

Marble <- function(X, R, Alpha, Gamma, MaxIterations = 100, KKT_Tolerance = 1e-4,
                   Convergence = 1e-4, Reproducible = TRUE, Verbose = TRUE) {
  if (Reproducible == TRUE){
    set.seed(915)
  }
  Dimensions <- dim(X)
  N <- length(Dimensions)
  Factor <- Initialize(Dimensions, R)
  u <- Factor$u
  A <- Factor$A
  Lambda <- runif(R)
  Lambda <- Lambda/sum(runif(R))
  Xi <- 0
  LatestObjective <- Inf
  for (Iteration in 1:MaxIterations) {
    for (n in 1:N) {
      Psi <- as.matrix(Equation9(u[-n]))
      B <- Equation6(A[[n]], Lambda)
      Pi <- Equation7(A[-n])
      XMatrix <- Matricize(X, n)
      for (l in 1:10) {
        Phi <- Equation12(XMatrix, Alpha, u[[n]], Psi, B, Pi)
        if (all(min(B, 1 - Phi) <= KKT_Tolerance)) break
        B <- B * Phi
      }
      B <- pmax(B, Xi * Gamma[n])
      Lambda <- colSums(B)
      A[[n]] <- sweep(B, 2, Lambda, `/`)
      for (l in 1:10) {
        Z <- Equation13(XMatrix, Alpha, u[[n]], Psi, B, Pi)
        if (all(min(u[[n]], 1 - Z) <= KKT_Tolerance)) break
        u[[n]] <- u[[n]] * Z
      }
      u[[n]] <- u[[n]] / sum(u[[n]])
    }
    C <- array(Alpha * (u[[1]] %*% Psi), dim = Dimensions)
    V <- array(B %*% Pi, dim = Dimensions)        
    M <- C + V
    head(C)
    if (any(M <= 0)) M[M <= 0] <- 1e-10
    Objective <- sum(X * log(M) - M)
    if (Verbose == TRUE){
      cat("Iteration:", Iteration, " --  Objective:", abs(Objective), "\n")
    }
    if (abs(LatestObjective - Objective) < Convergence) break
    LatestObjective <- Objective
    Xi <- Equation11(Xi, Objective, LatestObjective)
  }
  list(V = V, C = C)
}

########################
### Helper Functions ###
########################

Initialize <- function(Dimensions, R) {
  u <- lapply(Dimensions, function(d) {
    v <- runif(d)
    v / sum(v)
  })
  A <- lapply(Dimensions, function(d) {
    PopulateA <- matrix(runif(d * R), nrow = d, ncol = R)
    sweep(PopulateA, 2, colSums(PopulateA), `/`)
  })
  list(u = u, A = A)
}

Matricize <- function(X, mode) {
  Dimensions <- dim(X)          
  N <- length(Dimensions)        
  Permutation <- c(mode, setdiff(1:N, mode))  
  Reshaping <- aperm(X, Permutation)          
  matrix(Reshaping, nrow = Dimensions[mode], ncol = prod(Dimensions[-mode]))
}

Equation6 <- function(A, Lambda) {
  LambdaDiagonals <- diag(Lambda)
  B <- A %*% LambdaDiagonals
  return(B)
}

Equation7 <- function(A) {
  Result <- A[[1]]
  for (i in 2:length(A)) {
    Result <- khatri_rao(Result, A[[i]])
  }
  return(t(Result))
}

Equation9 <- function(u) {
  u <- lapply(u, function(v) {
    if (is.null(dim(v))) matrix(v, ncol = 1) else v
  })
  Result <- u[[1]]
  for (i in 2:length(u)) {
    Result <- khatri_rao(Result, u[[i]])
  }
  return(t(Result))
}

Equation11 <- function(Xi, CurrentF, PreviousF) {
  Kappa <- 1 - abs(PreviousF - CurrentF) / abs(PreviousF)
  SteppedXi <- max(Xi, 0.5 * Xi + 0.5 * Kappa)
  return(SteppedXi)
}

Equation12 <- function(MatrixX, Alpha, u, Psi, B, Pi) {
  Denominator <- Alpha * (u %*% Psi) + (B %*% Pi)
  HalfResult <- MatrixX / (Denominator + 1e-10)
  Phi <- HalfResult %*% t(Pi)
  return(Phi)
}

Equation13 <- function(MatrixX, Alpha, u, Psi, B, Pi) {
  Denominator <- Alpha * (u %*% Psi) + (B %*% Pi)
  HalfResult <- MatrixX / (Denominator + 1e-10)
  Z <- HalfResult %*% t(Psi)
  return(Z)
}
