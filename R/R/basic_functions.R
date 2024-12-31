###### Packages
# basic matrix function such as eye, ones, zeros, dot
library(pracma)
library(latex2exp)
library(glmnet)
library(dplyr)
### Kendall' tau
library(Kendall)
library(igraph)
library(parallel)
library(nloptr)
library(MASS)
library(stats) 

###### Solving Optimal Proportion
###Fisher-Info Utilities
fisher_logit_rank_single <- function(theta, phi, a) {
  # theta : numeric vector of parameter values (theta_1, theta_2, ..., theta_p)
  # phi   : scalar, the current value of phi
  # a     : integer vector c(i,j)
  #
  # Returns the 2x2 Fisher info matrix for a single observation X^(i,j)
  # under the logit "Dav" model.
  i <- a[1]
  j <- a[2]
  dtheta <- theta[i] - theta[j]
  v <- exp(phi)
  
  csh <- cosh(dtheta / 2)
  ssh <- sinh(dtheta / 2)
  denom <- (v + 2 * csh)^2
  
  I11 <- (v * csh + 2) / (2 * denom)     # I_{Delta theta, Delta theta}
  I12 <- -(v * ssh) / denom             # I_{Delta theta, phi}
  I22 <- (2 * v * csh) / denom          # I_{phi, phi}
  
  I_mat <- matrix(c(I11, I12,
                    I12, I22),
                  nrow = 2, byrow = TRUE)
  I_mat
}

construct_U <- function(p) {
  # p: dimension (integer > 1).
  # Returns p x (p-1) matrix U, whose columns are orthonormal to (1/sqrt(p))*1_p
  # and each has length=1
  U <- matrix(0, nrow = p, ncol = p - 1)
  for (k in 1:(p - 1)) {
    U[1:k, k] <- 1
    U[k + 1, k] <- -k
    norm_k <- sqrt(k * (k + 1))
    U[, k] <- U[, k] / norm_k
  }
  U
}

unit_vector <- function(p, i) {
  # returns e_i in R^p
  v <- numeric(p)
  v[i] <- 1
  v
}

fisher_info_block <- function(theta_vals, phi_val, pair_a, k_B, b) {
  # theta_vals : length p
  # phi_val    : scalar
  # pair_a     : c(i, j)
  # k_B        : total # of phi
  # b          : which phi => 1..k_B
  p <- length(theta_vals)
  i <- pair_a[1]
  j <- pair_a[2]
  
  # 1) single 2x2 from fisher_logit_rank_single
  I_single <- fisher_logit_rank_single(theta_vals, phi_val, pair_a)
  
  # 2) build e_ij
  e_i <- unit_vector(p, i)
  e_j <- unit_vector(p, j)
  eij <- e_i - e_j
  
  # 3) embed 2x2 into (p+k_B)x(p+k_B)
  M <- matrix(0, nrow = p + k_B, ncol = p + k_B)
  M[1:p, 1:p] <- I_single[1,1] * (eij %*% t(eij))
  M[1:p, p + b] <- I_single[1,2] * eij
  M[p + b, 1:p] <- I_single[2,1] * t(eij)
  M[p + b, p + b] <- I_single[2,2]
  
  M
}
### Moore-Penrose under rank=$d-1$ assumption, with small eps regularization}
moore_penrose_custom <- function(I_pi, p, k_B, eps=1e-8) {
  # I_pi : (d x d) matrix
  # rank(I_pi)=d-1
  # p,k_B => d=p+k_B
  # Steps:
  #   1) Q => d x (d-1), top-block=U_p, bottom-block= I_{k_B}
  #   2) I_pi_low_rank= t(Q)%*%I_pi%*%Q => (d-1)x(d-1)
  #   3) invert with regularization
  #   4) MPinv= Q * inv_low_rank * t(Q)
  d <- p + k_B
  U_p <- construct_U(p)
  Q <- matrix(0, nrow=d, ncol=d-1)
  
  # top block => (p x (p-1))= U_p
  Q[1:p, 1:(p-1)] <- U_p
  # bottom => identity across last k_B columns
  for (i in seq_len(k_B)) {
    Q[p + i, (p-1) + i] <- 1
  }
  
  I_pi_low_rank <- t(Q) %*% I_pi %*% Q
  # add eps diag for stability
  inv_low_rank <- solve(I_pi_low_rank + eps*diag(d-1))
  MPinv <- Q %*% inv_low_rank %*% t(Q)
  MPinv
}

### Objective function and gradient
F_pi_par <- function(pi_mat, par, connected_pair_matrix, k_B, H) {
  # pi_mat => (k x k_B)
  # par => c(theta, phi)
  # returns trace( H (I^pi)^+ )
  p <- length(par) - k_B
  theta_vals <- par[1:p]
  phi_vals   <- par[(p+1):(p+k_B)]
  
  k <- nrow(connected_pair_matrix)
  d <- p + k_B
  I_pi <- matrix(0, d, d)
  
  for (a in seq_len(k)) {
    pair_a <- connected_pair_matrix[a, ]
    for (b in seq_len(k_B)) {
      w_ab <- pi_mat[a,b]
      if (abs(w_ab) > 1e-15) {
        I_ab <- fisher_info_block(theta_vals, phi_vals[b], pair_a, k_B, b)
        I_pi <- I_pi + w_ab * I_ab
      }
    }
  }
  # moore penrose
  MPinv <- moore_penrose_custom(I_pi, p, k_B, eps=epss*sum(diag(I_pi))/p)
  val <- sum(diag(H %*% MPinv))
  val
}

gradF_pi_par <- function(pi_mat, par, connected_pair_matrix, k_B, H) {
  # partial => - trace( H (I^pi)^+ I_{a,b} (I^pi)^+ )
  p <- length(par) - k_B
  theta_vals <- par[1:p]
  phi_vals   <- par[(p+1):(p+k_B)]
  k <- nrow(connected_pair_matrix)
  d <- p + k_B
  
  # build I_pi
  I_pi <- matrix(0, d, d)
  for (a in seq_len(k)) {
    pair_a <- connected_pair_matrix[a,]
    for (b in seq_len(k_B)) {
      w_ab <- pi_mat[a,b]
      #if (abs(w_ab)>1e-15) {
      I_ab <- fisher_info_block(theta_vals, phi_vals[b], pair_a, k_B, b)
      I_pi <- I_pi + w_ab*I_ab
      #}
    }
  }
  MPinv <- moore_penrose_custom(I_pi, p, k_B, eps=epss*sum(diag(I_pi))/p)
  
  grad_mat <- matrix(0, k, k_B)
  for (a in seq_len(k)) {
    pair_a <- connected_pair_matrix[a,]
    for (b in seq_len(k_B)) {
      I_ab  <- fisher_info_block(theta_vals, phi_vals[b], pair_a, k_B, b)
      M_ab  <- MPinv %*% I_ab %*% MPinv
      val_ab<- sum(diag(H %*% M_ab))
      grad_mat[a,b] <- - val_ab
    }
  }
  grad_mat
}
### Flatten $\pi => x$, set eq constraints, box constraints, use nloptr
# equality: sum_{a=1..k} pi(a,b)= v_b  for each b
eq_constraints <- function(x, v_vec, k, k_B) {
  out <- numeric(k_B)
  for (b in seq_len(k_B)) {
    idx_b <- ((b-1)*k + 1):(b*k)
    sum_b <- sum(x[idx_b])
    out[b] <- sum_b - v_vec[b]
  }
  out
}

eq_jacobian <- function(x, v_vec, k, k_B) {
  # partial of out[b] wrt each x => 1 if x in that col
  nvar <- k*k_B
  jac  <- matrix(0, nrow=k_B, ncol=nvar)
  for (b in seq_len(k_B)) {
    idx_b <- ((b-1)*k +1):(b*k)
    jac[b, idx_b] <- 1
  }
  jac
}

# objective => f_objective, gradient => f_gradient
f_objective <- function(x, par, connected_pair_matrix, k_B, H, k) {
  # reshape
  pi_mat <- matrix(x, nrow=k, ncol=k_B)
  val <- F_pi_par(pi_mat, par, connected_pair_matrix, k_B, H)
  val
}

f_gradient <- function(x, par, connected_pair_matrix, k_B, H, k) {
  pi_mat <- matrix(x, nrow=k, ncol=k_B)
  g_mat  <- gradF_pi_par(pi_mat, par, connected_pair_matrix, k_B, H)
  grad_vec <- as.vector(g_mat)
  grad_vec
}

solve_pi_nloptr <- function(par, connected_pair_matrix, k_B, H,
                            v_vec, pi_init,
                            maxeval=10000,
                            method="NLOPT_LD_SLSQP") {
  k <- nrow(connected_pair_matrix)
  x0 <- as.vector(pi_init)
  
  # define lower, upper => x>=0
  lb <- rep(0, k*k_B)
  ub <- rep(Inf, k*k_B)
  
  res <- nloptr::nloptr(
    x0 = x0,
    eval_f      = function(x) f_objective(x, par, connected_pair_matrix, k_B, H, k),
    eval_grad_f = function(x) f_gradient(x, par, connected_pair_matrix, k_B, H, k),
    eval_g_eq    = function(x) eq_constraints(x, v_vec, k, k_B),
    eval_jac_g_eq= function(x) eq_jacobian(x, v_vec, k, k_B),
    lb = lb,
    ub = ub,
    opts = list(
      algorithm   = method,        # e.g. "NLOPT_LD_SLSQP"
      xtol_rel    = 1e-8,
      maxeval     = maxeval,
      print_level = 1
    )
  )
  
  # result
  x_sol <- res$solution
  pi_sol<- matrix(x_sol, nrow=k, ncol=k_B)
  list(
    status  = res$status,
    message = res$message,
    Fval    = res$objective,
    pi      = pi_sol
  )
}
### Graph generation
init_edge_generate_RA <- function(p, method = c("full", "regular", "min"), reg = 2) {
  #
  # p      : number of items (vertices)
  # method : one of "full", "regular", "min"
  # reg    : integer for 'regular' method => create random reg-regular graph
  #
  # returns: (k x 2) integer matrix connected_pair_matrix => each row is (i,j)
  #           with 1 <= i<j <= p (undirected). 
  method <- match.arg(method)
  if (p < 2) {
    stop("p must be >= 2 to form edges.")
  }
  if (method == "full") {
    ## 1) FULLY CONNECTED: all edges {i<j}, i=1..p-1, j=i+1..p
    edges_list <- list()
    idx <- 1
    for (i in seq_len(p-1)) {
      for (j in seq(i+1, p)) {
        edges_list[[idx]] <- c(i, j)
        idx <- idx + 1
      }
    }
    connected_pair_matrix <- do.call(rbind, edges_list)
  } else if (method == "regular") {
    ## 2) s-regular: use igraph::sample_k_regular
    ##    sample_k_regular(p, reg, directed=FALSE, multiple=FALSE).
    ##    This yields a random k-regular graph on p vertices.
    if (reg >= p) {
      stop("For a k-regular graph, reg must be < p in most typical scenarios.")
    }
    g_reg <- sample_k_regular(p, k=reg, directed=FALSE, multiple=FALSE)
    edges_ij <- as_edgelist(g_reg)
    # edges_ij is a (E x 2) matrix with vertex ids in {1..p}
    # typically iGraph might not ensure i<j, so we fix that:
    edges_sorted <- t(apply(edges_ij, 1, function(v) sort(v)))
    # remove duplicates if any (shouldn't happen for an undirected simple graph)
    edges_unique <- unique(edges_sorted)
    connected_pair_matrix <- edges_unique
    
  } else if (method == "min") {
    ## 3) random MST: 
    ##    we create a complete graph with random weights,
    ##    then pick the MST edges => subgraph with p-1 edges
    # step a) fully connected
    g_full <- make_full_graph(p, directed=FALSE)  # all edges
    # step b) assign random edge weights
    w_rand <- runif(ecount(g_full))
    g_full <- set_edge_attr(g_full, name="weight", value=w_rand)
    # step c) find MST => e.g. using 'mst' in igraph
    g_mst  <- mst(g_full, weights = E(g_full)$weight)
    edges_ij <- as_edgelist(g_mst)
    edges_sorted <- t(apply(edges_ij, 1, function(v) sort(v)))
    connected_pair_matrix <- edges_sorted
  }
  
  # ensure each row is i<j
  # no duplication
  # return 
  mode(connected_pair_matrix) <- "integer"
  connected_pair_matrix
}
### Generate a random parameter and inital proportion
init_par_test <- function(p, k_B, seed = NULL) {
  # p   : number of theta parameters
  # k_B : number of phi parameters
  # seed: optional integer to set.seed(...) if reproducibility is desired
  #
  # returns a numeric vector of length (p + k_B).
  #   par_test[1..p]   ~ runif( p, -3, 3)
  #   par_test[p+1..p+k_B] ~ runif( k_B, -3, 1)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  theta_vals <- runif(p, min=-3, max=3)
  phi_vals   <- runif(k_B, min=-3, max=1)
  par_test <- c(theta_vals, phi_vals)
  par_test
}
random_col <- function(k, total) {
  #
  # k     : integer > 0
  # total : numeric, the desired column sum
  #
  # returns a length-k vector x >= 0 with sum(x)= total
  #
  x <- runif(k)          # k uniform(0,1)
  x <- x / sum(x)        # normalize to sum=1
  x * total              # scale to sum= 'total'
}
init_pi_mat <- function(k, v_vec) {
  #
  # k      : row dimension (the 'k' in pi(a,b))
  # v_vec  : numeric vector of length k_B => each v_b >= 0, sum(v_vec)=1
  # returns a (k x k_B) matrix whose b-th column sums to v_vec[b]
  #
  k_B <- length(v_vec)
  pi_init <- matrix(0, nrow = k, ncol = k_B)
  for (b in seq_len(k_B)) {
    pi_init[, b] <- random_col(k, v_vec[b])
  }
  pi_init
}
### Generate $H$ matrix
build_H <- function(p, k_B, eta_theta = 1, eta = 0.5) {
  # build diagonal matrix of size d=p+k_B
  # first p diag entries =1, last k_B diag entries= eta
  d <- p + k_B
  diag_entries <- c(rep(eta_theta, p), rep(eta, k_B))
  diag(diag_entries, nrow=d, ncol=d)
}

build_soft_ranking_H <- function(par_true, p, eta = 0.5, delta = 1) {
  # Input:
  #   theta_star : numeric vector (the "true" parameter theta^*)
  #   delta      : nonnegative scalar
  # Output:
  #   H_mat : a length(theta_star) x length(theta_star) matrix
  # Formula:
  #   H_mat = sum_{1 <= i < j <= p} (|theta_star[i] - theta_star[j]| + delta)^2* (e_i - e_j)(e_i - e_j)^T
  # where p = length(theta_star).
  theta_star=par_true[1:p]
  d=length(par_true)
  k_B=d-p
  stopifnot(delta >= 0)
  # Initialize p x p matrix
  H_mat <- matrix(0, nrow = p, ncol = p)
  H_mat_full <- matrix(0, nrow = d, ncol = d)
  for(i in (p+1):d){
    H_mat_full[i, i]=eta
  }
  for (i in seq_len(p - 1)) {
    for (j in seq(i + 1, p)) {
      diff_ij <- abs(theta_star[i] - theta_star[j]) + delta
      val <- diff_ij^2
      # e_i - e_j => a length-p vector with +1 at i, -1 at j
      eij <- numeric(p)
      eij[i] <- 1
      eij[j] <- -1
      # outer product => eij %*% t(eij)
      H_mat <- H_mat + val * (eij %*% t(eij))
    }
  }
  H_mat_full[1:p,1:p]=H_mat/p
  return(H_mat_full)
}

###### SLSQP for MLE solver
### Data generator
expit <- function(z) {
  1 / (1 + exp(-z))
}

data_random_generator_roundrobin_BTD <- function(
    p,
    connected_pair_matrix,
    k_B,
    theta_min = -2,
    theta_max =  2,
    phi_min   = -3,
    phi_max   =  3,
    par_true  = NULL  # new optional argument
) {
  #
  # p  : dimension of theta
  # k_B: number of categories for tie parameter
  #
  # If par_true is not NULL and length(par_true)== p+k_B, we use that for theta,phi
  # Else we sample randomly in user-specified ranges.
  #
  # output => list( data_random=..., theta=..., phi=... )
  # data_random has columns X1..Xp, Y, cat_1..cat_k_B
  k <- nrow(connected_pair_matrix)
  nRows <- max(k, k_B)
  
  # 1) Decide on theta, phi
  use_random <- TRUE
  if (!is.null(par_true)) {
    # check length
    if (length(par_true) == p + k_B) {
      # use par_true
      use_random <- FALSE
      theta <- par_true[1:p]
      phi   <- par_true[(p+1) : (p + k_B)]
    }
  }
  
  if (use_random) {
    # sample randomly
    theta <- runif(p, min=theta_min, max=theta_max)
    phi   <- runif(k_B, min=phi_min, max=phi_max)
  }
  
  # 2) Build design + sample outcomes row by row
  design_matrix <- matrix(0, nrow = nRows, ncol = p)
  cat_matrix    <- matrix(0, nrow = nRows, ncol = k_B)
  Y_vec         <- integer(nRows)
  
  for (row_i in seq_len(nRows)) {
    i_pair <- ((row_i - 1) %% k) + 1
    i <- connected_pair_matrix[i_pair, 1]
    j <- connected_pair_matrix[i_pair, 2]
    
    b_star <- ((row_i - 1) %% k_B) + 1
    cat_matrix[row_i, b_star] <- 1
    
    # design => +1 at i, -1 at j
    design_matrix[row_i, i] <- +1
    design_matrix[row_i, j] <- -1
    
    # sample outcome
    Delta_ij <- theta[i] - theta[j]
    f_ij     <- expit(Delta_ij)
    oneMinus_f_ij <- expit(-Delta_ij)  # stable => 1 - f_ij
    v_star   <- exp(phi[b_star])
    sqrt_part<- sqrt(f_ij * oneMinus_f_ij)
    denom    <- 1 + v_star * sqrt_part
    
    p_plus  <-  f_ij         / denom
    p_minus <-  oneMinus_f_ij/ denom
    p_zero  <- (v_star * sqrt_part) / denom
    
    randu <- runif(1)
    if (randu < p_plus) {
      Y_vec[row_i] <- +1L
    } else if (randu < p_plus + p_minus) {
      Y_vec[row_i] <- -1L
    } else {
      Y_vec[row_i] <- 0L
    }
  }
  
  # 3) Build final df
  df_out <- as.data.frame(design_matrix)
  colnames(df_out) <- paste0("X", seq_len(p))
  df_out$Y <- Y_vec
  
  cat_colnames <- paste0("cat_", seq_len(k_B))
  colnames(cat_matrix) <- cat_colnames
  df_out <- cbind(df_out, cat_matrix)
  
  list(
    data_random = df_out,
    theta       = theta,
    phi         = phi
  )
}
### Negative log-likelihood and gradient for BTD model
neg_loglik_btd_ball_unified_noloop <- function(par, df_data, p, k_B) {
  #
  # par => c(theta, phi), length= p+k_B
  # df_data => columns X1..Xp, Y in {+1,-1,0}, cat_1..cat_k_B
  #
  theta <- par[1:p]
  phi   <- par[(p+1):(p+k_B)]
  
  Xmat   <- as.matrix(df_data[, paste0("X", seq_len(p)), drop=FALSE])
  Yvec   <- df_data$Y
  catMat <- as.matrix(df_data[, paste0("cat_", seq_len(k_B)), drop=FALSE])
  n      <- nrow(Xmat)
  
  bIndex <- catMat %*% seq_len(k_B)
  DeltaVec <- as.numeric(Xmat %*% theta)
  
  fVec       <- expit(DeltaVec)
  oneMinus_f <- expit(-DeltaVec)
  phiVec     <- phi[bIndex]
  vVec       <- exp(phiVec)
  
  sqrtPartVec <- sqrt(fVec* oneMinus_f)
  DenVec      <- 1 + vVec* sqrtPartVec
  
  p1 <- fVec       / DenVec
  p2 <- oneMinus_f / DenVec
  p3 <- (vVec* sqrtPartVec)/ DenVec
  
  Iplus  <- (Yvec==+1)
  Iminus <- (Yvec==-1)
  Izero  <- (Yvec==0)
  
  pChosen<- p1*Iplus + p2*Iminus + p3*Izero
  -sum(log(pChosen ))
}

grad_neg_loglik_btd_ball_unified_noloop <- function(par, df_data, p, k_B) {
  theta <- par[1:p]
  phi   <- par[(p+1):(p+k_B)]
  
  Xmat  <- as.matrix(df_data[, paste0("X", seq_len(p)), drop=FALSE])
  Yvec  <- df_data$Y
  catMat<- as.matrix(df_data[, paste0("cat_", seq_len(k_B)), drop=FALSE])
  n     <- nrow(Xmat)
  
  bIndex <- catMat %*% seq_len(k_B)
  DeltaVec <- as.numeric(Xmat %*% theta)
  
  eNegDelta <- exp(-DeltaVec)
  phiVec    <- phi[bIndex]
  vVec      <- exp(phiVec)
  
  DenVec    <- 1 + eNegDelta + vVec* exp(-DeltaVec/2)
  
  # partial wrt Delta => dNLL_dDelta
  # cX(Y)= (1 -Y)/2 => 0 if Y=+1, 1 if Y=-1, 0.5 if Y=0
  cXVec <- (1 - Yvec)/2
  dDen_dx <- - eNegDelta - 0.5* vVec* exp(-DeltaVec/2)
  dNLL_dDeltaVec <- cXVec + dDen_dx/ DenVec
  
  # 1) gradient wrt theta => X^T %*% dNLL_dDeltaVec
  grad_theta <- as.numeric(t(Xmat) %*% dNLL_dDeltaVec)
  
  # 2) gradient wrt phi
  # partial wrt phi => (dDen_dphi/DenVec) - i0, where i0= 1-Y^2 => 1 if Y=0, else0
  i0 <- 1 - (Yvec^2)
  dDen_dphiVec <- exp(-DeltaVec/2)* vVec
  dNLL_dphiVec <- (dDen_dphiVec / DenVec) - i0
  
  grad_phi <- as.numeric(rowsum(dNLL_dphiVec, group= bIndex))
  
  c(grad_theta, grad_phi)
}

### MLE solver
fit_btd_mle_ball_noloop_zerosum <- function(
    df_data, p, k_B, R = 2,
    phi_lower = -3, phi_upper = 3,
    init_par   = NULL,
    method     = "NLOPT_LD_SLSQP",
    maxeval    = 1000,
    verbose    = FALSE
) {
  #
  # df_data => columns X1..Xp, Y in {+1,-1,0}, cat_1..cat_k_B
  # p,k_B => dimension => length(par)= p+k_B
  # R => sum(theta^2)<= R^2
  # sum(theta)=0
  # phi in [phi_lower,phi_upper]
  # init_par => if null => random
  # method => nloptr method => "NLOPT_LD_SLSQP"
  # maxeval => iteration limit
  # verbose => if TRUE => print each iteration's objective
  #
  # returns => nloptr result with an added $obj_vals => the entire obj. history
  
  stopifnot(requireNamespace("nloptr", quietly=TRUE))
  
  nParam <- p + k_B
  
  # 1) Build init if needed
  if (is.null(init_par)) {
    # random then project inside ball => sum(theta^2)<= R^2
    init_theta <- runif(p, min=-1, max=1)
    nt <- sum(init_theta^2)
    if (nt > R^2) {
      init_theta <- init_theta*( R / sqrt(nt))
    }
    # also shift => sum(theta)=0
    shift_val<- mean(init_theta)
    init_theta<- init_theta - shift_val
    
    # random phi in [phi_lower, phi_upper]
    init_phi<- runif(k_B, min= phi_lower, max= phi_upper)
    init_par<- c(init_theta, init_phi)
  } else {
    if (length(init_par)!= nParam) {
      stop("init_par length mismatch => must be p + k_B.")
    }
  }
  
  # environment to track iteration # and objective
  histEnv <- new.env()
  histEnv$iterCount <- 0
  histEnv$objVals   <- numeric(0)
  
  # define objective + gradient
  obj_fn <- function(x) {
    val <- neg_loglik_btd_ball_unified_noloop(x, df_data, p, k_B)
    histEnv$iterCount <- histEnv$iterCount + 1
    histEnv$objVals[ histEnv$iterCount ] <- val
    if (verbose) {
      cat(sprintf("Iter=%d, Obj=%.6f\n", histEnv$iterCount, val))
    }
    val
  }
  grad_fn<- function(x) {
    grad_neg_loglik_btd_ball_unified_noloop(x, df_data, p, k_B)
  }
  
  # eq => sum(theta)=0
  eq_fn <- function(x) {
    the<- x[1:p]
    c(sum(the))
  }
  eq_grad_fn <- function(x) {
    c(rep(1,p), rep(0,k_B))
  }
  
  # ineq => sum(theta^2)- R^2 <=0
  ineq_fn <- function(x) {
    the<- x[1:p]
    c(sum(the^2)- R^2)
  }
  ineq_grad_fn<- function(x) {
    the<- x[1:p]
    c(2*the, rep(0,k_B))
  }
  
  # box => phi in [phi_lower, phi_upper], no constraints on theta
  lb <- c(rep(-Inf,p), rep(phi_lower,k_B))
  ub <- c(rep( Inf,p), rep(phi_upper,k_B))
  
  # call nloptr => SLSQP
  res <- nloptr::nloptr(
    x0= init_par,
    eval_f      = obj_fn,
    eval_grad_f = grad_fn,
    eval_g_eq    = eq_fn,
    eval_jac_g_eq= eq_grad_fn,
    eval_g_ineq  = ineq_fn,
    eval_jac_g_ineq= ineq_grad_fn,
    lb= lb,
    ub= ub,
    opts= list(
      algorithm   = method,
      xtol_rel    = 1e-8,
      maxeval     = maxeval,
      print_level = 0
    )
  )
  # store the entire objective history
  res$obj_vals <- histEnv$objVals
  res
}
###### Adaptive sampling
### Adaptive Sampling Utilities
one_step_sample_b <- function(v, vec_b) {
  v <- v/sum(v)
  k_B <- length(v)
  counts<- numeric(k_B)
  for (catv in vec_b) {
    counts[catv]<- counts[catv]+1
  }
  n_prev<- length(vec_b)
  if(n_prev>0) {
    f<- counts/n_prev
  } else {
    f<- rep(0,k_B)
  }
  gap<- v - f
  b_star<- which.max(gap)
  b_star
}

compute_pi_n_general <- function(df_data, connected_pair_matrix, p, k_B) {
  k <- nrow(connected_pair_matrix)
  pi_mat<- matrix(0, k, k_B)
  
  # helper => find row in cpm that matches (i,j) or (j,i)
  find_row_in_cpm <- function(i,j) {
    hits<- which(
      (connected_pair_matrix[,1]== i & connected_pair_matrix[,2]== j) |
        (connected_pair_matrix[,1]== j & connected_pair_matrix[,2]== i)
    )
    if(length(hits)<1) {
      return(NA_integer_)
    } else {
      return(hits[1])
    }
  }
  
  n_rows<- nrow(df_data)
  for(r in seq_len(n_rows)) {
    rowX<- as.numeric(df_data[r, paste0("X", seq_len(p)), drop=FALSE])
    i_candidates<- which(rowX== +1)
    j_candidates<- which(rowX== -1)
    if(length(i_candidates)!=1 || length(j_candidates)!=1) {
      next
    }
    i<- i_candidates[1]
    j<- j_candidates[1]
    
    rowCats<- as.numeric(df_data[r, paste0("cat_", seq_len(k_B)), drop=FALSE])
    b_idx<- which(rowCats== 1)
    if(length(b_idx)!=1) {
      next
    }
    a_idx<- find_row_in_cpm(i,j)
    if(is.na(a_idx)) {
      next
    }
    pi_mat[a_idx,b_idx]<- pi_mat[a_idx,b_idx]+1
  }
  total_count<- sum(pi_mat)
  if(total_count>0) {
    pi_mat<- pi_mat/ total_count
  }
  pi_mat
}

select_action_a <- function(
    par,   # c(theta, phi)
    H,     # (p+k_B)x(p+k_B)
    pi_n,  # (k x k_B)
    b,     # chosen category in {1..k_B}
    connected_pair_matrix, # (k x 2)
    p, k_B,
    eps=1e-8
) {
  d<- p+k_B
  I_pi<- matrix(0,d,d)
  
  # build I^pi
  for(a_row in seq_len(nrow(connected_pair_matrix))) {
    pair_a <- connected_pair_matrix[a_row,]
    for(b_col in seq_len(k_B)) {
      w_ab<- pi_n[a_row,b_col]
      if(abs(w_ab)>1e-15) {
        # build fisher_info_block => I_{a_row,b_col}
        I_ab<- fisher_info_block(
          theta_vals= par[1:p],
          phi_val   = par[(p+1):(p+k_B)][b_col],
          pair_a    = pair_a,
          k_B       = k_B,
          b         = b_col
        )
        I_pi<- I_pi + w_ab* I_ab
      }
    }
  }
  
  # invert => moore_penrose_custom(I_pi, p, k_B, eps)
  MPinv<- moore_penrose_custom(I_pi, p, k_B, eps= 0)
  
  # partial wrt each a => a= a_row
  partial_vals<- numeric(nrow(connected_pair_matrix))
  for(a_row in seq_len(nrow(connected_pair_matrix))) {
    pair_ab<- connected_pair_matrix[a_row,]
    # fisher block => I_{a_row,b}
    I_ab<- fisher_info_block(
      theta_vals= par[1:p],
      phi_val   = par[(p+1):(p+k_B)][b],
      pair_a    = pair_ab,
      k_B       = k_B,
      b         = b
    )
    M_ab<- MPinv %*% I_ab %*% MPinv
    val_ab<- sum(diag(H %*% M_ab))
    partial_vals[a_row]<-  val_ab
  }
  a_star<- which.max(partial_vals)
  a_star
}
build_vec_b_from_df <- function(df_data, k_B) {
  # df_data: must have columns cat_1..cat_k_B
  # k_B: total categories
  #
  # returns: an integer vector vec_b_init with length = nrow(df_data),
  #          where each entry is the category b in {1..k_B} for that row.
  n_rows <- nrow(df_data)
  vec_b <- integer(n_rows)
  
  for (i in seq_len(n_rows)) {
    # read the cat_ columns
    row_cat <- as.numeric(df_data[i, paste0("cat_", seq_len(k_B)), drop=FALSE])
    b_index <- which(row_cat == 1)
    if (length(b_index) == 1) {
      vec_b[i] <- b_index
    } else {
      # if 0 or multiple 1's => skip or handle error
      # for safety, we can do:
      vec_b[i] <- NA_integer_
    }
  }
  vec_b
}
###### Design and estimation
single_round_adaptive <- function(
    df_data,
    connected_pair_matrix,
    p, k_B,
    R=3,
    phi_lower=-3, phi_upper=3,
    init_par=NULL,
    method="NLOPT_LD_SLSQP",
    maxeval=2000,
    verbose=FALSE,
    v,          # vector for one_step_sample_b
    vec_b,      # history of categories
    par_true,   # real param => c(theta_true,phi_true)
    H,          # for select_action_a
    eps=1e-8,
    pi_n=NULL   # if not NULL => we skip computing from df_data
) {
  
  # Step 2) Fit MLE:
  res_fit <- fit_btd_mle_ball_noloop_zerosum(
    df_data= df_data, p= p, k_B= k_B,
    R= R, phi_lower= phi_lower, phi_upper= phi_upper,
    init_par= init_par,
    method= method, maxeval= maxeval,
    verbose= verbose
  )
  par_mle <- res_fit$solution   # c(theta*, phi*)
  
  # Step 3(a) => pick b_star
  b_star <- one_step_sample_b(v, vec_b)
  
  # Step 3(b) => pi_n
  if (is.null(pi_n)) {
    # compute from df_data => more general approach
    pi_n <- compute_pi_n_general(df_data, connected_pair_matrix, p, k_B)
  }
  
  # select a => using select_action_a(par=..., H=..., pi_n=..., b=..., ...)
  a_star <- select_action_a(
    par= par_mle,
    H= H,
    pi_n= pi_n,
    b= b_star,
    connected_pair_matrix= connected_pair_matrix,
    p= p, k_B= k_B,
    eps= eps
  )
  
  # Step 4) sample new outcome => from par_true => 
  i_sel <- connected_pair_matrix[a_star, 1]
  j_sel <- connected_pair_matrix[a_star, 2]
  theta_true <- par_true[1:p]
  phi_true   <- par_true[(p+1):(p+k_B)]
  v_sel <- exp(phi_true[b_star])
  
  Delta_sel<- theta_true[i_sel] - theta_true[j_sel]
  f_sel    <- expit(Delta_sel)
  oneMinus_f_sel<- expit(-Delta_sel)
  sqrt_part<- sqrt(f_sel* oneMinus_f_sel)
  denom_sel<- 1+ v_sel* sqrt_part
  
  p_plus  <- f_sel         / denom_sel
  p_minus <- oneMinus_f_sel/ denom_sel
  p_zero  <- (v_sel* sqrt_part)/ denom_sel
  
  randu <- runif(1)
  Y_new <- 0L
  if (randu< p_plus) {
    Y_new<- +1L
  } else if (randu< p_plus + p_minus) {
    Y_new<- -1L
  } else {
    Y_new<- 0L
  }
  
  # Step 5) build new row => X1..Xp, Y, cat_1..cat_k_B
  new_design<- rep(0, p)
  new_design[i_sel]<- +1
  new_design[j_sel]<- -1
  new_cat<- rep(0, k_B)
  new_cat[b_star]<- 1
  
  new_row_data<- c(new_design, Y_new, new_cat)
  coln_data   <- c(paste0("X", seq_len(p)), "Y", paste0("cat_", seq_len(k_B)))
  new_dfrow   <- as.data.frame(t(new_row_data))
  colnames(new_dfrow)<- coln_data
  
  df_data_new<- rbind(df_data, new_dfrow)
  
  # update vec_b => add b_star
  vec_b_new<- c(vec_b, b_star)
  
  # also update pi_n => we do the incremental approach:
  # if old pi_n was based on nrow(df_data) samples => total = n
  # new row => n+1 => we do
  # pi_n_new= n*pi_n + delta_{(a_star,b_star)} => then / (n+1)
  n_old <- nrow(df_data)
  pi_n_new <- pi_n * n_old  # scale back to counts
  pi_n_new[a_star, b_star] <- pi_n_new[a_star, b_star] + 1
  n_new <- n_old + 1
  pi_n_new<- pi_n_new / n_new
  
  # return
  list(
    df_data   = df_data_new,
    pi_n      = pi_n_new,
    vec_b     = vec_b_new,
    par_mle   = par_mle,
    res_mle   = res_fit,
    a_star    = a_star,
    b_star    = b_star,
    Y_new     = Y_new
  )
}


###############################################################################
### multiple_rounds_adaptive_dynamicH => re-computes H each iteration
### and stores par_mle in par_mle_tensor
###############################################################################
multiple_rounds_adaptive_dynamicH <- function(
    T,
    df_data_init,
    pi_n_init = NULL,
    connected_pair_matrix,
    p, k_B,
    R = 3,
    phi_lower = -3, phi_upper = 3,
    init_par = NULL,
    method = "NLOPT_LD_SLSQP",
    maxeval = 2000,
    verbose = FALSE,
    v,                 # target proportion => one_step_sample_b
    vec_b_init,        # initial category history
    par_true,          # real param => c(theta_true, phi_true)
    H_dynamic_fun,     # function(par_mle) => build dynamic H
    eps = 1e-8
) {
  # 1) Initialize
  df_data <- df_data_init
  pi_n    <- pi_n_init
  vec_b   <- vec_b_init
  
  if (is.null(pi_n)) {
    # compute from df_data
    pi_n <- compute_pi_n_general(df_data, connected_pair_matrix, p, k_B)
  }
  
  # total new samples to add:
  n_init <- nrow(df_data)
  max_steps <- T - n_init
  if (max_steps < 0) {
    stop("Initial data size > T.")
  }
  
  # We'll store pi_n in a 3D array => (k, k_B, max_steps+1).
  # We'll also store MLE in a 2D array => (p + k_B, max_steps+1).
  # Because we do an extra final MLE after the loop, let's allocate (max_steps+1).
  k <- nrow(connected_pair_matrix)
  pi_tensor <- array(0, dim=c(k, k_B, max_steps + 1))
  pi_tensor[,,1] <- pi_n
  
  d_param <- p + k_B
  par_mle_tensor <- matrix(0, nrow=d_param, ncol=max_steps + 1)
  
  # Possibly do initial MLE if you want the MLE for the *initial data*
  # Let's do it so that par_mle_tensor[,1] is MLE from df_data_init
  if (n_init > 0) {
    # Fit MLE on initial data
    res_mle_init <- fit_btd_mle_ball_noloop_zerosum(
      df_data  = df_data,
      p        = p,
      k_B      = k_B,
      R        = R,
      phi_lower= phi_lower,
      phi_upper= phi_upper,
      init_par = init_par,
      method   = method,
      maxeval  = maxeval,
      verbose  = verbose
    )
    par_mle_init <- res_mle_init$solution
    par_mle_tensor[, 1] <- par_mle_init
  } else {
    # If there is no initial data, or we want to skip that, we can do random or store init_par
    if (!is.null(init_par)) {
      par_mle_tensor[,1] <- init_par
    }
  }
  
  # We'll keep track of the current MLE to build dynamic H each iteration
  current_par_mle <- par_mle_tensor[,1]
  
  # 2) main loop => each iteration adds exactly 1 new sample
  for (step_i in seq_len(max_steps)) {
    # (a) Build dynamic H from the *current* MLE
    H_dynamic <- H_dynamic_fun(current_par_mle)
    
    # (b) pick b_star => largest gap from v
    b_star <- one_step_sample_b(v, vec_b)
    
    # (c) select action a_star => row
    a_star <- select_action_a(
      par = current_par_mle,
      H   = H_dynamic,
      pi_n= pi_n,
      b   = b_star,
      connected_pair_matrix= connected_pair_matrix,
      p= p, k_B= k_B,
      eps= eps
    )
    
    # (d) sample new outcome from par_true => (theta_true, phi_true)
    i_sel <- connected_pair_matrix[a_star,1]
    j_sel <- connected_pair_matrix[a_star,2]
    theta_true <- par_true[1:p]
    phi_true   <- par_true[(p+1):(p+k_B)]
    v_sel <- exp(phi_true[b_star])
    
    Delta_sel<- theta_true[i_sel] - theta_true[j_sel]
    f_sel    <- expit(Delta_sel)
    oneMinus_f_sel<- expit(-Delta_sel)
    sqrt_part<- sqrt(f_sel * oneMinus_f_sel)
    denom_sel<- 1 + v_sel * sqrt_part
    
    p_plus  <- f_sel / denom_sel
    p_minus <- oneMinus_f_sel/ denom_sel
    p_zero  <- (v_sel * sqrt_part)/ denom_sel
    
    randu <- runif(1)
    Y_new <- 0L
    if (randu < p_plus) {
      Y_new<- +1L
    } else if (randu < p_plus + p_minus) {
      Y_new<- -1L
    } else {
      Y_new<- 0L
    }
    
    # (e) append new row => df_data
    new_design<- rep(0, p)
    new_design[i_sel]<- +1
    new_design[j_sel]<- -1
    new_cat<- rep(0, k_B)
    new_cat[b_star]<- 1
    
    new_row_data<- c(new_design, Y_new, new_cat)
    coln_data   <- c(paste0("X", seq_len(p)), "Y", paste0("cat_",seq_len(k_B)))
    new_dfrow   <- as.data.frame(t(new_row_data))
    colnames(new_dfrow)<- coln_data
    
    df_data <- rbind(df_data, new_dfrow)
    vec_b   <- c(vec_b, b_star)
    
    # (f) update pi_n => incremental
    n_old <- n_init + step_i - 1  # or simply nrow(df_data)-1
    pi_n <- pi_n * n_old
    pi_n[a_star, b_star] <- pi_n[a_star, b_star] + 1
    pi_n <- pi_n / (n_old + 1)
    
    pi_tensor[,, step_i + 1] <- pi_n
    
    # (g) Fit MLE again => so next iteration uses new MLE as warm start
    res_mle <- fit_btd_mle_ball_noloop_zerosum(
      df_data  = df_data,
      p        = p,
      k_B      = k_B,
      R        = R,
      phi_lower= phi_lower,
      phi_upper= phi_upper,
      init_par = current_par_mle,
      method   = method,
      maxeval  = maxeval,
      verbose  = verbose
    )
    current_par_mle <- res_mle$solution
    
    # store in par_mle_tensor => step_i+1
    par_mle_tensor[, step_i + 1] <- current_par_mle
  }
  
  # 3) after we have T rows => final MLE if we want it:
  # We already did MLE after the last sample in the loop => current_par_mle
  # is consistent with the final dataset. So no extra step is needed
  # unless you want an additional solve. We'll just say this is final MLE.
  final_par_mle <- current_par_mle
  
  # finalize => how many columns we used
  used_slices <- max_steps + 1  # we have up to that many
  pi_tensor_used <- pi_tensor[,, seq_len(used_slices), drop=FALSE]
  par_mle_tensor_used <- par_mle_tensor[, seq_len(used_slices), drop=FALSE]
  
  list(
    df_data         = df_data,
    pi_n            = pi_n,
    vec_b           = vec_b,
    final_par_mle   = final_par_mle,
    pi_tensor       = pi_tensor_used,
    par_mle_tensor  = par_mle_tensor_used
  )
}



ranking_eval_kendall <- function(par_mle_tensor, p, theta_true) {
  # par_mle_tensor : (p + k_B) x N matrix 
  #                  each column => MLE = (theta_1..theta_p, phi_1..phi_kB)
  # p : integer, dimension for theta
  # theta_true : numeric vector of length p => the true theta
  #
  # returns: numeric vector of length N => Kendall's tau for each column
  
  # 1) load or ensure the "Kendall" function is available
  #    typically from library(Kendall)
  if (!requireNamespace("Kendall", quietly = TRUE)) {
    stop("Package 'Kendall' must be installed to compute Kendall's tau.")
  }
  
  # 2) figure out how many columns => number of iterations
  N <- ncol(par_mle_tensor)
  
  # 3) initialize output
  tau_vec <- numeric(N)
  
  # 4) for each column => extract MLE’s theta part => compare to theta_true
  for (i in seq_len(N)) {
    # MLE’s theta => first p entries of column i
    mle_theta <- par_mle_tensor[1:p, i]
    # compute Kendall’s tau
    # 'Kendall' returns a list with $tau => the correlation
    tau_val <- Kendall::Kendall(mle_theta, theta_true)$tau
    tau_vec[i] <- tau_val
  }
  
  # 5) return the entire vector of Kendall’s tau
  tau_vec
}



###############################################################################
### data_random_generator_MST_BTD
###   - Takes p, connected_pair_matrix, k_B, plus optional par_true
###   - Builds MST from connected_pair_matrix => edges_mst
###   - Let nEdges_mst = p - 1
###   - Then we produce nRows = max(nEdges_mst, k_B)
###   - Round-robin: for row_i in 1..nRows:
###        i_edge = ((row_i-1) mod nEdges_mst) + 1   => MST edge
###        i_cat  = ((row_i-1) mod k_B) + 1          => tie category
###     Then sample outcome with either par_true or random (theta, phi).
###############################################################################
data_random_generator_MST_BTD <- function(
    p,
    connected_pair_matrix,
    k_B,
    theta_min = -2,
    theta_max =  2,
    phi_min   = -3,
    phi_max   =  3,
    par_true  = NULL  # optional (theta, phi)
) {
  #
  # p  : dimension of theta
  # connected_pair_matrix: (k x 2) integer => each row => (i,j)
  # k_B: number of categories for tie parameter
  #
  # If par_true is not NULL and length(par_true)== p+k_B, we use that for theta,phi.
  # Else we sample randomly in user-specified ranges for theta, phi.
  #
  # Then build an MST from connected_pair_matrix:
  #   1) we treat connected_pair_matrix as edges => build a graph in igraph
  #   2) if the graph is not yet fully connected, or if you want random weights,
  #      we can do so => see code below
  #   3) edges_mst => p-1 edges
  #
  # nEdges_mst = p-1
  # nRows = max(nEdges_mst, k_B)
  # for row_i=1..nRows:
  #   i_edge = ((row_i-1) %% nEdges_mst)+1
  #   i_cat  = ((row_i-1) %% k_B)+1
  #   build design => +1 at i, -1 at j
  #   sample outcome from the 3-outcome BTD model
  #
  # output => list( data_random=..., theta=..., phi=... )
  # data_random has columns X1..Xp, Y in {+1,-1,0}, cat_1..cat_k_B
  #
  # The MST approach means we only use p-1 edges from connected_pair_matrix
  # then do a round-robin among those MST edges.
  
  # ensure igraph
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' must be installed to build MST from connected_pair_matrix.")
  }
  
  # 1) Decide on theta, phi
  use_random <- TRUE
  if (!is.null(par_true)) {
    if (length(par_true) == p + k_B) {
      # use par_true
      use_random <- FALSE
      theta <- par_true[1:p]
      phi   <- par_true[(p+1):(p+k_B)]
    }
  }
  if (use_random) {
    # sample randomly
    theta <- runif(p, min=theta_min, max=theta_max)
    phi   <- runif(k_B, min=phi_min, max=phi_max)
  }
  
  # 2) Build MST from connected_pair_matrix
  #    step (a) => create an igraph with p vertices => labeled 1..p
  #    step (b) => edges come from connected_pair_matrix => each row => (i,j)
  #    step (c) => assign random weights or set all = runif?
  #    step (d) => compute MST => yields (p-1) edges => edges_mst
  
  # Make a graph from connected_pair_matrix => edges in {1..p}
  # We can store random weights for MST or just do uniform
  g0 <- igraph::make_empty_graph(n=p, directed=FALSE)
  # edges => a vector c(i1, j1, i2, j2, ...)
  edges_vec <- t(connected_pair_matrix)
  edges_vec <- as.vector(edges_vec)
  g_full <- igraph::add_edges(g0, edges= edges_vec)
  
  # assign random weights => runif(ecount)
  ecount_full <- igraph::ecount(g_full)
  w_rand <- runif(ecount_full)
  g_full <- igraph::set_edge_attr(g_full, name="weight", value=w_rand)
  
  # build MST => p-1 edges
  g_mst  <- igraph::mst(g_full, weights= igraph::E(g_full)$weight)
  edges_ij_mst <- igraph::as_edgelist(g_mst)
  # edges_ij_mst is a matrix => (#(p-1) x 2)
  # each row => c(i, j) in {1..p}
  # ensure i<j for consistency
  edges_sorted <- t(apply(edges_ij_mst, 1, function(v) sort(v)))
  mode(edges_sorted) <- "integer"
  edges_mst <- edges_sorted
  
  nEdges_mst <- nrow(edges_mst) # should be p-1
  
  # 3) nRows => max(nEdges_mst, k_B)
  nRows <- max(nEdges_mst, k_B)
  
  # 4) Construct design_matrix, cat_matrix, Y
  design_matrix <- matrix(0, nrow=nRows, ncol=p)
  cat_matrix    <- matrix(0, nrow=nRows, ncol=k_B)
  Y_vec         <- integer(nRows)
  
  # define expit => 1/(1+exp(-z))
  expit <- function(z) 1 / (1 + exp(-z))
  
  for (row_i in seq_len(nRows)) {
    # a) pick which edge => i_edge
    i_edge <- ((row_i - 1) %% nEdges_mst) + 1
    i <- edges_mst[i_edge,1]
    j <- edges_mst[i_edge,2]
    
    # b) pick category => i_cat
    i_cat <- ((row_i - 1) %% k_B) + 1
    cat_matrix[row_i, i_cat] <- 1
    
    # design => +1 at i, -1 at j
    design_matrix[row_i, i] <- +1
    design_matrix[row_i, j] <- -1
    
    # c) sample outcome => from tie-augmented BTD
    #    Delta_ij = theta[i] - theta[j]
    #    f_ij = expit(Delta_ij)
    Delta_ij <- theta[i] - theta[j]
    f_ij     <- expit(Delta_ij)
    oneMinus_f_ij <- expit(-Delta_ij)
    # tie param => phi[i_cat]
    v_star   <- exp(phi[i_cat])
    sqrt_part<- sqrt(f_ij * oneMinus_f_ij)
    denom    <- 1 + v_star* sqrt_part
    
    p_plus  <- f_ij / denom
    p_minus <- oneMinus_f_ij / denom
    p_zero  <- (v_star* sqrt_part)/ denom
    
    randu <- runif(1)
    if (randu < p_plus) {
      Y_vec[row_i] <- +1L
    } else if (randu < (p_plus + p_minus)) {
      Y_vec[row_i] <- -1L
    } else {
      Y_vec[row_i] <- 0L
    }
  }
  
  # 5) build final data frame
  df_out <- as.data.frame(design_matrix)
  colnames(df_out) <- paste0("X", seq_len(p))
  df_out$Y <- Y_vec
  colnames(cat_matrix) <- paste0("cat_", seq_len(k_B))
  df_out <- cbind(df_out, cat_matrix)
  
  list(
    data_random = df_out,
    theta       = theta,
    phi         = phi,
    edges_mst   = edges_mst    # optionally return the MST edges used
  )
}


### 1) Preliminary: Helper to compute 3-outcome distribution & Shannon Entropy
three_outcome_entropy <- function(
    a,             # row index in connected_pair_matrix
    b_star,        # chosen category
    par_mle,       # c(theta, phi)
    connected_pair_matrix,
    p, k_B
) {
  # 1) parse par_mle => c(theta, phi)
  theta_hat <- par_mle[1:p]
  phi_hat   <- par_mle[(p+1):(p+k_B)]
  
  # 2) extract the edge => (i,j)
  i <- connected_pair_matrix[a,1]
  j <- connected_pair_matrix[a,2]
  
  # 3) distribution:
  #    Delta_ij = theta_hat[i] - theta_hat[j]
  #    f_ij = expit(Delta_ij)
  #    v_b  = exp(phi_hat[b_star])
  expit <- function(z) 1/(1+exp(-z))
  
  Delta_ij <- theta_hat[i] - theta_hat[j]
  f_ij     <- expit(Delta_ij)
  oneMinus_f_ij <- expit(-Delta_ij)
  v_b      <- exp(phi_hat[b_star])
  sqrt_part<- sqrt(f_ij * oneMinus_f_ij)
  denom    <- 1 + v_b * sqrt_part
  
  p_plus  <- f_ij          / denom
  p_minus <- oneMinus_f_ij / denom
  p_zero  <- (v_b*sqrt_part)/ denom
  
  # 4) Shannon entropy => - sum(p(y)*log p(y)) for y in {+1,-1,0}, ignoring zeros
  probs <- c(p_plus, p_minus, p_zero)
  probs_pos <- probs[ probs>1e-15 ]
  # 
  ent <- - sum( probs_pos * log(probs_pos) )   # natural log
  
  ent
}

### 2) Function: multiple_rounds_adaptive_unif
multiple_rounds_adaptive_unif <- function(
    n,  # final total sample size
    df_data_init,
    pi_n_init = NULL,
    connected_pair_matrix,
    p, k_B,
    R = 3,
    phi_lower = -3, phi_upper = 3,
    init_par = NULL,
    method = "NLOPT_LD_SLSQP",
    maxeval = 2000,
    verbose = FALSE,
    v,
    vec_b_init,
    par_true,
    eps = 1e-8
) {
  # 1) Initialize
  df_data <- df_data_init
  pi_n    <- pi_n_init
  vec_b   <- vec_b_init
  if (is.null(pi_n)) {
    pi_n <- compute_pi_n_general(df_data, connected_pair_matrix, p, k_B)
  }
  
  k <- nrow(connected_pair_matrix)
  n_init <- nrow(df_data)
  max_steps <- n - n_init
  if (max_steps < 0) {
    stop("Initial data size > n.")
  }
  
  # We'll store pi_n in a 3D array => (k, k_B, max_steps+1)
  pi_tensor <- array(0, dim=c(k, k_B, max_steps+1))
  pi_tensor[,,1] <- pi_n
  
  # We'll also store MLE in a 2D array => (p + k_B, max_steps+1)
  d_param <- p + k_B
  par_mle_tensor <- matrix(0, nrow=d_param, ncol=max_steps+1)
  
  # Possibly do initial MLE
  if (n_init>0) {
    res_mle_init <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=init_par,
      method=method, maxeval=maxeval, verbose=verbose
    )
    par_mle_tensor[,1] <- res_mle_init$solution
  } else {
    if (!is.null(init_par)) {
      par_mle_tensor[,1] <- init_par
    }
  }
  
  current_par_mle <- par_mle_tensor[,1]
  
  # 2) main loop
  for (step_i in seq_len(max_steps)) {
    # (a) pick b_star => one_step_sample_b
    b_star <- one_step_sample_b(v, vec_b)
    
    # (b) pick a_star => uniform among 1..k
    a_star <- sample.int(k, size=1)
    
    # (c) sample outcome from par_true => same approach
    i_sel <- connected_pair_matrix[a_star,1]
    j_sel <- connected_pair_matrix[a_star,2]
    
    theta_true <- par_true[1:p]
    phi_true   <- par_true[(p+1):(p+k_B)]

    Delta_ij <- theta_true[i_sel] - theta_true[j_sel]
    f_ij     <- expit(Delta_ij)
    oneMinus_f_ij<- expit(-Delta_ij)
    v_sel    <- exp(phi_true[b_star])
    sqrt_part<- sqrt(f_ij*oneMinus_f_ij)
    denom_sel<- 1+ v_sel* sqrt_part
    
    p_plus  <- f_ij/ denom_sel
    p_minus <- oneMinus_f_ij/ denom_sel
    p_zero  <- (v_sel* sqrt_part)/ denom_sel
    
    randu <- runif(1)
    Y_new <- 0L
    if(randu < p_plus) {
      Y_new <- +1L
    } else if(randu < p_plus + p_minus) {
      Y_new <- -1L
    } else {
      Y_new <- 0L
    }
    
    # (d) append new row => df_data
    new_design<- rep(0, p)
    new_design[i_sel]<- +1
    new_design[j_sel]<- -1
    new_cat<- rep(0, k_B)
    new_cat[b_star]<-1
    
    row_data<- c(new_design, Y_new, new_cat)
    coln_data<- c(paste0("X", seq_len(p)), "Y", paste0("cat_", seq_len(k_B)))
    new_dfrow<- as.data.frame(t(row_data))
    colnames(new_dfrow)<- coln_data
    
    df_data <- rbind(df_data, new_dfrow)
    vec_b   <- c(vec_b, b_star)
    
    # (e) update pi_n => incremental
    n_old <- n_init + step_i - 1
    pi_n <- pi_n * n_old
    pi_n[a_star,b_star] <- pi_n[a_star,b_star] + 1
    pi_n <- pi_n/(n_old+1)
    pi_tensor[,, step_i+1] <- pi_n
    
    # (f) Fit MLE => warm start
    res_mle <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=current_par_mle,
      method=method, maxeval=maxeval, verbose=verbose
    )
    current_par_mle <- res_mle$solution
    par_mle_tensor[, step_i+1] <- current_par_mle
  }
  
  final_par_mle <- current_par_mle
  used_slices <- max_steps+1
  pi_tensor_used <- pi_tensor[,, seq_len(used_slices), drop=FALSE]
  par_mle_tensor_used <- par_mle_tensor[, seq_len(used_slices), drop=FALSE]
  
  list(
    df_data        = df_data,
    pi_n           = pi_n,
    vec_b          = vec_b,
    final_par_mle  = final_par_mle,
    pi_tensor      = pi_tensor_used,
    par_mle_tensor = par_mle_tensor_used
  )
}

multiple_rounds_adaptive_uncertainty <- function(
    n,
    df_data_init,
    pi_n_init = NULL,
    connected_pair_matrix,
    p, k_B,
    R = 3,
    phi_lower = -3, phi_upper = 3,
    init_par = NULL,
    method = "NLOPT_LD_SLSQP",
    maxeval = 2000,
    verbose = FALSE,
    v,
    vec_b_init,
    par_true,
    eps = 1e-8
) {
  # 1) Initialize
  df_data <- df_data_init
  pi_n    <- pi_n_init
  vec_b   <- vec_b_init
  if (is.null(pi_n)) {
    pi_n <- compute_pi_n_general(df_data, connected_pair_matrix, p, k_B)
  }
  
  k <- nrow(connected_pair_matrix)
  n_init <- nrow(df_data)
  max_steps <- n - n_init
  if (max_steps < 0) {
    stop("Initial data size > n.")
  }
  
  pi_tensor <- array(0, dim=c(k, k_B, max_steps+1))
  pi_tensor[,,1] <- pi_n
  
  d_param <- p + k_B
  par_mle_tensor <- matrix(0, nrow=d_param, ncol=max_steps+1)
  
  if (n_init>0) {
    res_mle_init <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=init_par,
      method=method, maxeval=maxeval, verbose=verbose
    )
    par_mle_tensor[,1] <- res_mle_init$solution
  } else {
    if (!is.null(init_par)) {
      par_mle_tensor[,1] <- init_par
    }
  }
  
  current_par_mle <- par_mle_tensor[,1]
  
  # 2) main loop
  for (step_i in seq_len(max_steps)) {
    # (a) pick b_star => largest gap from v
    b_star <- one_step_sample_b(v, vec_b)
    
    # (b) among a=1..k => find the row that max. Shannon entropy
    #    using the "current" MLE => c(theta, phi)
    #    we define a helper to compute => see "three_outcome_entropy" above
    best_a <- NA
    best_ent <- -Inf
    for (a in seq_len(k)) {
      ent_a <- three_outcome_entropy(
        a      = a,
        b_star = b_star,
        par_mle= current_par_mle,
        connected_pair_matrix= connected_pair_matrix,
        p= p, k_B= k_B
      )
      if (ent_a > best_ent) {
        best_ent <- ent_a
        best_a   <- a
      }
    }
    a_star <- best_a
    
    # (c) sample outcome from par_true => same approach
    i_sel <- connected_pair_matrix[a_star,1]
    j_sel <- connected_pair_matrix[a_star,2]
    
    theta_true <- par_true[1:p]
    phi_true   <- par_true[(p+1):(p+k_B)]
    
    Delta_ij <- theta_true[i_sel] - theta_true[j_sel]
    f_ij     <- expit(Delta_ij)
    oneMinus_f_ij<- expit(-Delta_ij)
    v_sel    <- exp(phi_true[b_star])
    sqrt_part<- sqrt(f_ij*oneMinus_f_ij)
    denom_sel<- 1+ v_sel* sqrt_part
    
    p_plus  <- f_ij/ denom_sel
    p_minus <- oneMinus_f_ij/ denom_sel
    p_zero  <- (v_sel* sqrt_part)/ denom_sel
    
    randu <- runif(1)
    Y_new <- 0L
    if(randu < p_plus) {
      Y_new <- +1L
    } else if(randu < p_plus + p_minus) {
      Y_new <- -1L
    } else {
      Y_new <- 0L
    }
    
    # (d) append new row => df_data
    new_design<- rep(0, p)
    new_design[i_sel]<- +1
    new_design[j_sel]<- -1
    new_cat<- rep(0, k_B)
    new_cat[b_star]<-1
    
    row_data<- c(new_design, Y_new, new_cat)
    coln_data<- c(paste0("X", seq_len(p)), "Y", paste0("cat_", seq_len(k_B)))
    new_dfrow<- as.data.frame(t(row_data))
    colnames(new_dfrow)<- coln_data
    
    df_data <- rbind(df_data, new_dfrow)
    vec_b   <- c(vec_b, b_star)
    
    # (e) update pi_n => incremental
    n_old <- n_init + step_i - 1
    pi_n <- pi_n * n_old
    pi_n[a_star,b_star] <- pi_n[a_star,b_star] + 1
    pi_n <- pi_n/(n_old+1)
    pi_tensor[,, step_i+1] <- pi_n
    
    # (f) Fit MLE => warm start
    res_mle <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=current_par_mle,
      method=method, maxeval=maxeval, verbose=verbose
    )
    current_par_mle <- res_mle$solution
    par_mle_tensor[, step_i+1] <- current_par_mle
  }
  
  final_par_mle <- current_par_mle
  used_slices <- max_steps+1
  pi_tensor_used <- pi_tensor[,, seq_len(used_slices), drop=FALSE]
  par_mle_tensor_used <- par_mle_tensor[, seq_len(used_slices), drop=FALSE]
  
  list(
    df_data        = df_data,
    pi_n           = pi_n,
    vec_b          = vec_b,
    final_par_mle  = final_par_mle,
    pi_tensor      = pi_tensor_used,
    par_mle_tensor = par_mle_tensor_used
  )
}


###########################################################################
# Implements the heterogeneous soft ranking loss:
hetero_soft_ranking_loss <- function(
    theta,          # numeric vector of length p
    theta_star,     # numeric vector of length p (true parameter)
    gamma,          # numeric or list for the tie parameters (estimated)
    gamma_star,     # numeric or list for the tie parameters (true)
    C         = NULL,   # function(x) => link 'C'. If NULL => C(x) = exp(x).
    Cprime0   = NULL,   # numeric => value of C'(0). If NULL & C=exp => =1.
    u         = NULL,   # function(x) => 'u(x)'. If NULL => u(x)= x+delta
    delta     = 1,      # numeric >=0 => shift for u(0)= delta
    eta       = 1       # numeric >0 => weighting for tie-parameter penalty
) {
  ###########################################################################
  # Implements the heterogeneous soft ranking loss:
  #
  #   L_{h-soft}^{C,u}(\theta,theta^*, gamma, gamma^*)
  #    = (1/(2p)) * sum_{1<=i<j<=p} [ h_ij^{C,u}+ h_ji^{C,u} ]
  #      + eta * sum_b || gamma_b - gamma_b^* ||^2
  #
  # with:
  #    h_ij^{C,u} = hC( - ((Delta_ij - Delta_ij^*) * sign(Delta_ij^*) * u(|Delta_ij^*|)) )
  #    hC(x) = C(x) - C'(0)* x - C(0)
  # and sign(0) => +1.
  ###########################################################################
  
  # 1) Dimensions & basic checks
  p <- length(theta)
  stopifnot(length(theta_star) == p)
  
  # 2) Default C => exp(x)
  #    => C'(0)=1, C(0)=1
  if (is.null(C)) {
    C <- function(x) exp(x)
    if (is.null(Cprime0)) {
      Cprime0 <- 1.0
    }
    C0_value <- 1.0
  } else {
    # user-supplied C => must also have C'(0)
    if (is.null(Cprime0)) {
      stop("If you provide a custom C, you must also provide 'Cprime0' = C'(0).")
    }
    C0_value <- C(0)
  }
  
  # 3) Bregman divergence => hC(x) = C(x) - C'(0)* x - C(0)
  hC <- function(x) {
    C(x) - Cprime0*x - C0_value
  }
  
  # 4) if u is NULL => use u(x)= x + delta
  if (is.null(u)) {
    u <- function(x) x + delta
  }
  
  # 5) special sign => sign(0)= +1
  safe_sign <- function(x) {
    if (x > 0) {
      1
    } else if (x < 0) {
      -1
    } else {
      1  # sign(0)= +1
    }
  }
  
  # 6) compute the ranking sum => sum_{i<j} [ h_ij + h_ji ] / 2
  ranking_sum <- 0.0
  for (i in seq_len(p-1)) {
    for (j in (i+1):p) {
      d_ij     <- theta[i] - theta[j]
      dstar_ij <- theta_star[i] - theta_star[j]
      s_ij_star <- safe_sign(dstar_ij)
      mag_ij_star<- abs(dstar_ij)
      
      arg_ij <- -((d_ij - dstar_ij)* s_ij_star * u(mag_ij_star))
      hij    <- hC(arg_ij)
      
      # symmetrical => Delta_ji= -d_ij, Delta_ji^*= -dstar_ij
      arg_ji <- -((-(d_ij) - (-(dstar_ij))) * safe_sign(-dstar_ij) * u(abs(-dstar_ij)))
      hji    <- hC(arg_ji)
      
      ranking_sum <- ranking_sum + 0.5*(hij + hji)
    }
  }
  # then multiply by 1/p
  ranking_term <- ranking_sum / p
  
  # 7) penalty => eta * sum( (gamma - gamma_star)^2 )
  # flatten if list
  if (is.list(gamma)) {
    gamma_vec <- unlist(gamma)
  } else {
    gamma_vec <- gamma
  }
  if (is.list(gamma_star)) {
    gamma_star_vec <- unlist(gamma_star)
  } else {
    gamma_star_vec <- gamma_star
  }
  stopifnot(length(gamma_vec) == length(gamma_star_vec))
  penalty_sum <- sum((gamma_vec - gamma_star_vec)^2)
  
  # 8) final
  loss_value <- ranking_term + eta * penalty_sum
  loss_value
}


# Suppose p=3, k_B=2 => theta,gamma each dimension
theta       <- c(1.0, -0.3, 0.8)
theta_star  <- c(0.9, -0.2, 0.5)
gamma       <- c(0.1, 0.2)   # tie param (estimated)
gamma_star  <- c(0.0, 0.0)   # tie param (true)

# Evaluate the default with sign(0)= +1, delta=1, eta=0.5
val_loss <- hetero_soft_ranking_loss(
  theta      = theta,
  theta_star = theta_star,
  gamma      = gamma,
  gamma_star = gamma_star,
  delta      = 1,
  eta        = 0.5
)


### (a) Uncertainty sampling function
multiple_rounds_uncertainty <- function(
    T,
    df_data_init,
    pi_n_init = NULL,
    connected_pair_matrix,
    p, k_B,
    R = R_theta,
    phi_lower = -3, phi_upper = 3,
    init_par = NULL,
    method = "NLOPT_LD_SLSQP",
    maxeval = 100,
    verbose = FALSE,
    v,
    vec_b_init,
    par_true,
    eps = 1e-8
) {
  # Similar to multiple_rounds_adaptive_dynamicH, but
  # we pick b_star from one_step_sample_b,
  # then pick a_star by maximizing the Shannon entropy
  # (like you had in your previous 'uncertainty' approach).
  # We'll store pi_n, par_mle_tensor, etc.
  
  # 1) Initialization
  df_data <- df_data_init
  pi_n    <- pi_n_init
  vec_b   <- vec_b_init
  if (is.null(pi_n)) {
    pi_n <- compute_pi_n_general(df_data, connected_pair_matrix, p, k_B)
  }
  
  k <- nrow(connected_pair_matrix)
  n_init <- nrow(df_data)
  max_steps <- T - n_init
  if (max_steps < 0) {
    stop("Initial data size > T.")
  }
  
  # We'll store pi_n in a 3D array => (k, k_B, max_steps+1)
  pi_tensor <- array(0, dim=c(k, k_B, max_steps+1))
  pi_tensor[,,1] <- pi_n
  
  # We'll also store MLE in a 2D array => (p + k_B, max_steps+1)
  d_param <- p + k_B
  par_mle_tensor <- matrix(0, nrow=d_param, ncol=max_steps+1)
  
  # Possibly do initial MLE
  if (n_init>0) {
    res_mle_init <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=init_par,
      method=method, maxeval=maxeval, verbose=verbose
    )
    par_mle_tensor[,1] <- res_mle_init$solution
  } else {
    if (!is.null(init_par)) {
      par_mle_tensor[,1] <- init_par
    }
  }
  
  current_par_mle <- par_mle_tensor[,1]
  
  # Helper for Shannon entropy
  three_outcome_entropy <- function(a, b_star, par_mle) {
    # a => row in connected_pair_matrix => (i,j)
    i_sel <- connected_pair_matrix[a,1]
    j_sel <- connected_pair_matrix[a,2]
    theta_hat <- par_mle[1:p]
    phi_hat   <- par_mle[(p+1):(p+k_B)]
    Delta_ij  <- theta_hat[i_sel] - theta_hat[j_sel]
    f_ij      <- expit(Delta_ij)
    oneMinus_f_ij <- expit(-Delta_ij)
    v_sel     <- exp(phi_hat[b_star])
    sqrt_part <- sqrt(f_ij * oneMinus_f_ij)
    denom     <- 1 + v_sel* sqrt_part
    
    p_plus  <- f_ij/denom
    p_minus <- oneMinus_f_ij/denom
    p_zero  <- (v_sel* sqrt_part)/denom
    probs   <- c(p_plus, p_minus, p_zero)
    
    # Shannon entropy => -sum_{y} p(y) log p(y)
    # ignoring p(y)<1e-15 if needed
    probs_pos <- probs[probs>1e-15]
    ent <- -sum(probs_pos * log(probs_pos))
    ent
  }
  
  # 2) main loop
  for (step_i in seq_len(max_steps)) {
    # pick b_star => largest gap from v
    b_star <- one_step_sample_b(v, vec_b)
    
    # among a=1..k => pick that row that yields maximum entropy
    best_a <- 1
    best_entropy <- -Inf
    for (a_candidate in seq_len(k)) {
      ent_val <- three_outcome_entropy(a_candidate, b_star, current_par_mle)
      if (ent_val > best_entropy) {
        best_entropy <- ent_val
        best_a <- a_candidate
      }
    }
    a_star <- best_a
    
    # sample outcome from par_true
    i_sel <- connected_pair_matrix[a_star,1]
    j_sel <- connected_pair_matrix[a_star,2]
    
    # same logic => compute p_plus, p_minus, p_zero
    Delta_ij<- par_true[1:p][i_sel] - par_true[1:p][j_sel]
    f_ij    <- expit(Delta_ij)
    oneMinus_f_ij<- expit(-Delta_ij)
    v_sel   <- exp(par_true[(p+1):(p+k_B)][b_star])
    sqrt_part<- sqrt(f_ij* oneMinus_f_ij)
    denom_sel<- 1+ v_sel* sqrt_part
    
    p_plus  <- f_ij/ denom_sel
    p_minus <- oneMinus_f_ij/ denom_sel
    p_zero  <- (v_sel* sqrt_part)/ denom_sel
    
    randu <- runif(1)
    Y_new <- 0L
    if(randu < p_plus) {
      Y_new<- +1L
    } else if(randu < p_plus + p_minus) {
      Y_new<- -1L
    } else {
      Y_new<- 0L
    }
    
    # append row => df_data
    new_design<- rep(0, p)
    new_design[i_sel]<- +1
    new_design[j_sel]<- -1
    new_cat<- rep(0, k_B)
    new_cat[b_star]<-1
    
    row_data<- c(new_design, Y_new, new_cat)
    coln_data<- c(paste0("X", seq_len(p)), "Y", paste0("cat_", seq_len(k_B)))
    new_dfrow<- as.data.frame(t(row_data))
    colnames(new_dfrow)<- coln_data
    
    df_data <- rbind(df_data, new_dfrow)
    vec_b   <- c(vec_b, b_star)
    
    # update pi_n => incremental
    n_old <- n_init + step_i - 1
    pi_n <- pi_n * n_old
    pi_n[a_star,b_star] <- pi_n[a_star,b_star] + 1
    pi_n <- pi_n/(n_old+1)
    pi_tensor[,, step_i+1] <- pi_n
    
    # Fit MLE => warm start
    res_mle <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=current_par_mle,
      method=method, maxeval=maxeval, verbose=verbose
    )
    current_par_mle <- res_mle$solution
    par_mle_tensor[, step_i+1] <- current_par_mle
  }
  
  final_par_mle <- current_par_mle
  used_slices <- max_steps+1
  pi_tensor_used <- pi_tensor[,, seq_len(used_slices), drop=FALSE]
  par_mle_tensor_used <- par_mle_tensor[, seq_len(used_slices), drop=FALSE]
  
  list(
    df_data        = df_data,
    pi_n           = pi_n,
    vec_b          = vec_b,
    final_par_mle  = final_par_mle,
    pi_tensor      = pi_tensor_used,
    par_mle_tensor = par_mle_tensor_used
  )
}


### (b) Uniform sampling
multiple_rounds_uniform <- function(
    T,
    df_data_init,
    pi_n_init = NULL,
    connected_pair_matrix,
    p, k_B,
    R = 3,
    phi_lower = -3, phi_upper = 3,
    init_par = NULL,
    method = "NLOPT_LD_SLSQP",
    maxeval = 2000,
    verbose = FALSE,
    v,
    vec_b_init,
    par_true,
    eps = 1e-8
) {
  # Each iteration picks b_star from one_step_sample_b,
  # picks a_star uniformly from 1..k,
  # then sample => new row => re-fit MLE.
  # store par_mle_tensor, etc.
  
  df_data <- df_data_init
  pi_n    <- pi_n_init
  vec_b   <- vec_b_init
  if (is.null(pi_n)) {
    pi_n <- compute_pi_n_general(df_data, connected_pair_matrix, p, k_B)
  }
  
  k <- nrow(connected_pair_matrix)
  n_init <- nrow(df_data)
  max_steps <- T - n_init
  if (max_steps < 0) {
    stop("Initial data size > T.")
  }
  
  pi_tensor <- array(0, dim=c(k, k_B, max_steps+1))
  pi_tensor[,,1] <- pi_n
  
  d_param <- p + k_B
  par_mle_tensor <- matrix(0, nrow=d_param, ncol=max_steps+1)
  
  if (n_init>0) {
    res_mle_init <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=init_par,
      method=method, maxeval=maxeval, verbose=verbose
    )
    par_mle_tensor[,1] <- res_mle_init$solution
  } else {
    if (!is.null(init_par)) {
      par_mle_tensor[,1] <- init_par
    }
  }
  
  current_par_mle <- par_mle_tensor[,1]
  
  for (step_i in seq_len(max_steps)) {
    # pick b_star => largest gap
    b_star <- one_step_sample_b(v, vec_b)
    
    # pick a_star => uniform
    a_star <- sample.int(k, size=1)
    
    # sample outcome from par_true
    i_sel <- connected_pair_matrix[a_star,1]
    j_sel <- connected_pair_matrix[a_star,2]
    Delta_ij <- par_true[1:p][i_sel] - par_true[1:p][j_sel]
    f_ij     <- expit(Delta_ij)
    oneMinus_f_ij<- expit(-Delta_ij)
    v_sel    <- exp(par_true[(p+1):(p+k_B)][b_star])
    sqrt_part<- sqrt(f_ij*oneMinus_f_ij)
    denom_sel<- 1+ v_sel* sqrt_part
    
    p_plus  <- f_ij/ denom_sel
    p_minus <- oneMinus_f_ij/ denom_sel
    p_zero  <- (v_sel* sqrt_part)/ denom_sel
    
    randu <- runif(1)
    Y_new <- 0L
    if(randu < p_plus) {
      Y_new<- +1L
    } else if(randu < p_plus + p_minus) {
      Y_new<- -1L
    } else {
      Y_new<- 0L
    }
    
    new_design<- rep(0, p)
    new_design[i_sel]<- +1
    new_design[j_sel]<- -1
    new_cat<- rep(0, k_B)
    new_cat[b_star]<-1
    
    row_data<- c(new_design, Y_new, new_cat)
    coln_data<- c(paste0("X", seq_len(p)), "Y", paste0("cat_", seq_len(k_B)))
    new_dfrow<- as.data.frame(t(row_data))
    colnames(new_dfrow)<- coln_data
    
    df_data <- rbind(df_data, new_dfrow)
    vec_b   <- c(vec_b, b_star)
    
    n_old <- n_init + step_i -1
    pi_n <- pi_n * n_old
    pi_n[a_star,b_star] <- pi_n[a_star,b_star] + 1
    pi_n <- pi_n/(n_old+1)
    pi_tensor[,, step_i+1] <- pi_n
    
    # fit MLE => warm start
    res_mle <- fit_btd_mle_ball_noloop_zerosum(
      df_data=df_data, p=p, k_B=k_B,
      R=R, phi_lower=phi_lower, phi_upper=phi_upper,
      init_par=current_par_mle,
      method=method, maxeval=maxeval, verbose=verbose
    )
    current_par_mle <- res_mle$solution
    par_mle_tensor[, step_i+1] <- current_par_mle
  }
  
  final_par_mle <- current_par_mle
  used_slices <- max_steps+1
  pi_tensor_used <- pi_tensor[,, seq_len(used_slices), drop=FALSE]
  par_mle_tensor_used <- par_mle_tensor[, seq_len(used_slices), drop=FALSE]
  
  list(
    df_data        = df_data,
    pi_n           = pi_n,
    vec_b          = vec_b,
    final_par_mle  = final_par_mle,
    pi_tensor      = pi_tensor_used,
    par_mle_tensor = par_mle_tensor_used
  )
}
