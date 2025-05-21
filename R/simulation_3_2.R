setwd("/panfs/jay/groups/27/shenx/zhao1118/soft_ranking_loss")

############################################################
## Section 1 – Functions                                  ##
############################################################

# sample_X_ab() : draw one observation X^{a,b}.
sample_X_ab <- function(theta, sigma2_vec, a, b, d_mat) {
  if (!(a %in% seq_len(ncol(d_mat)))) stop("`a` must be in 1:k.")
  if (!(b %in% seq_along(sigma2_vec))) stop("`b` must be in 1:k'.")
  mu <- sum(d_mat[, a] * theta)
  sd <- sqrt(sigma2_vec[b])
  rnorm(1, mean = mu, sd = sd)
}

# sample_X_vec() : vectorised version, returns X_1,…,X_n.
sample_X_vec <- function(theta, sigma2_vec, a_vec, b_vec, d_mat) {
  if (length(a_vec) != length(b_vec))
    stop("`a_vec` and `b_vec` must have the same length.")
  if (!all(a_vec %in% seq_len(ncol(d_mat))))
    stop("elements of `a_vec` must be in 1:k.")
  if (!all(b_vec %in% seq_along(sigma2_vec)))
    stop("elements of `b_vec` must be in 1:k'.")
  mu     <- as.vector(crossprod(theta, d_mat[, a_vec, drop = FALSE]))
  sd_vec <- sqrt(sigma2_vec[b_vec])
  rnorm(length(mu), mean = mu, sd = sd_vec)
}

# fit_mle_bca() : constrained MLE via block-coordinate ascent.
fit_mle_bca <- function(d_mat, a_vec, b_vec, X_vec,
                        max_iter = 50, radius = 5,
                        lower_bd = 0.1, upper_bd = 10,
                        tol = 1e-10) {
  
  if (length(a_vec) != length(b_vec) || length(b_vec) != length(X_vec))
    stop("`a_vec`, `b_vec`, and `X_vec` must have the same length.")
  if (!all(a_vec %in% seq_len(ncol(d_mat))))
    stop("elements of `a_vec` must be in 1:k.")
  
  k_prime  <- max(b_vec)
  theta_vec <- rep(0, nrow(d_mat))
  res       <- X_vec
  var_vec   <- tapply(res^2, b_vec, mean)
  var_vec   <- pmax(pmin(var_vec, upper_bd), lower_bd)
  if (length(var_vec) < k_prime)
    var_vec <- rep(var_vec, length.out = k_prime)
  
  D_samp <- d_mat[, a_vec, drop = FALSE]
  
  for (iter in seq_len(max_iter)) {
    w           <- 1 / var_vec[b_vec]
    G           <- D_samp %*% diag(w) %*% t(D_samp)
    g           <- D_samp %*% (w * X_vec)
    theta_new   <- solve(G, g)
    
    nrm <- sqrt(sum(theta_new^2))
    if (nrm > radius) theta_new <- theta_new * (radius / nrm)
    
    res     <- X_vec - as.vector(crossprod(theta_new, D_samp))
    var_new <- tapply(res^2, b_vec, mean)
    var_new <- pmax(pmin(var_new, upper_bd), lower_bd)
    if (length(var_new) < k_prime)
      var_new <- rep(var_new, length.out = k_prime)
    
    if (max(abs(theta_new - theta_vec), abs(var_new - var_vec)) < tol) {
      theta_vec <- theta_new
      var_vec   <- var_new
      break
    }
    theta_vec <- theta_new
    var_vec   <- var_new
  }
  
  list(theta_hat = as.numeric(theta_vec),
       var_hat   = as.numeric(var_vec),
       iterations = iter,
       converged  = (iter < max_iter))
}


fisher_single <- function(a, b, d_mat, var_vec) {
  p <- nrow(d_mat)
  k_prime <- length(var_vec)
  out <- matrix(0, p + k_prime, p + k_prime)
  d  <- d_mat[, a, drop = FALSE]
  out[1:p, 1:p] <- (1 / var_vec[b]) * (d %*% t(d))
  out[p + b, p + b] <- 1 / (2 * var_vec[b]^2)
  out
}

sequential_GI1 <- function(theta, sigma2_vec, b_vec, d_mat, a_vec_init,
                           radius = 5, lower_bd = 0.1, upper_bd = 10) {
  k        <- ncol(d_mat)
  k_prime  <- length(sigma2_vec)
  n        <- length(b_vec)
  m0       <- length(a_vec_init)
  if (m0 < 2)
    stop("`a_vec_init` must have length ≥ 2.")
  if (length(unique(a_vec_init)) < 2)
    stop("`a_vec_init` must contain at least three distinct design indices.")
  if (length(unique(b_vec[seq_len(m0)])) == 1)
    stop("The first `m0` entries of `b_vec` must include at least two different tier values.")
  if (!all(a_vec_init %in% seq_len(k)))
    stop("`a_vec_init` contains indices outside 1:k.")
  if (!all(b_vec %in% seq_len(k_prime)))
    stop("`b_vec` contains indices outside 1:k'.")
  if (!all(seq_len(k_prime) %in% b_vec[seq_len(m0)]))
    stop("The first `m0` entries of `b_vec` must include every tier 1,…,k'.")
  a_vec <- integer(n)
  X_vec <- numeric(n)
  a_vec[1:m0] <- a_vec_init
  X_vec[1:m0] <- sample_X_vec(theta, sigma2_vec, a_vec_init, b_vec[1:m0], d_mat)
  for (t in m0:(n - 1)) {
    mle <- fit_mle_bca(d_mat, a_vec[1:t], b_vec[1:t], X_vec[1:t],
                       radius = radius, lower_bd = lower_bd,
                       upper_bd = upper_bd)
    var_hat <- mle$var_hat
    p <- nrow(d_mat)
    info_tot <- matrix(0, p + k_prime, p + k_prime)
    for (j in 1:t){info_tot <- info_tot+fisher_single(a_vec[j], b_vec[j], d_mat, var_hat)}
    Sigma_hat <- solve(info_tot)
    
    scores <- numeric(k)
    for (a in 1:k) {
      Icand <- fisher_single(a, b_vec[t + 1], d_mat, var_hat)
      scores[a] <- sum(diag(Sigma_hat %*% Icand %*% Sigma_hat))
    }
    a_next <- which.max(scores)
    a_vec[t + 1] <- a_next
    X_vec[t + 1] <- sample_X_ab(theta, sigma2_vec,
                                a_next, b_vec[t + 1], d_mat)
  }
  
  list(a_vec = a_vec, X_vec = X_vec)
}

############################################################
## Section 2 – Problem setting                            ##
############################################################
p       <- 2
k       <- 3
k_prime <- 2

theta       <- c(1, -1)
sigma2_vec  <- c(1, 2)
set.seed(2025)

d_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)

############################################################
## Section 3 – Simulation & MLE                           ##
############################################################
set.seed(123)
n  <- 100
a_init <- c(1, 2, 3)          # five initial designs (≥3 distinct)
b_vec <- rep(1:2, length.out = n) 

res <- sequential_GI1(theta, sigma2_vec, b_vec, d_mat, a_init)
head(res$a_vec)
head(res$X_vec)
final_mle <- fit_mle_bca(d_mat, res$a_vec, b_vec, res$X_vec)


############################################################
## Section 4 – Optimal design π* (minimise tr{I(π)⁻¹}) ##
############################################################

library(nloptr)

# solve_pi_opt: solve min_{π≥0} tr{[I(π)+ridge·I]⁻¹}
# subject to ∑ₐ π_{a,b}=v_b  and  ∑_{a,b}π_{a,b}=1
solve_pi_opt <- function(d_mat,
                         sigma2_vec,
                         v_vec,
                         maxeval = 5000,
                         ridge   = 1e-8,
                         lb_eps  = 1e-6) {
  
  p       <- nrow(d_mat)
  k       <- ncol(d_mat)
  k_prime <- length(sigma2_vec)
  
  stopifnot(length(v_vec) == k_prime,
            abs(sum(v_vec) - 1) < 1e-10)
  
  # Fisher information block for design a, tier b
  J_block <- function(a, b) {
    J <- matrix(0, p + k_prime, p + k_prime)
    d <- d_mat[, a, drop = FALSE]
    J[1:p, 1:p]     <- (1 / sigma2_vec[b]) * (d %*% t(d))
    J[p + b, p + b] <- 1 / (2 * sigma2_vec[b]^2)
    J
  }
  
  # Assemble full Fisher information I(π)
  I_from_pi <- function(pi_vec) {
    I <- matrix(0, p + k_prime, p + k_prime)
    for (i in seq_along(pi_vec)) {
      w <- pi_vec[i]
      if (w <= 0) next
      a <- ((i - 1) %% k) + 1
      b <- ((i - 1) %/% k) + 1
      I <- I + w * J_block(a, b)
    }
    I
  }
  
  # Objective: trace of the inverse (with ridge stabilization)
  f_objective <- function(x) {
    Iinv <- solve(I_from_pi(x) + diag(ridge, p + k_prime))
    sum(diag(Iinv))
  }
  
  # Gradient: ∂/∂π_i tr(I⁻¹) = -tr(I⁻¹ J_i I⁻¹)
  f_gradient <- function(x) {
    Iinv <- solve(I_from_pi(x) + diag(ridge, p + k_prime))
    grad <- numeric(length(x))
    for (i in seq_along(x)) {
      a <- ((i - 1) %% k) + 1
      b <- ((i - 1) %/% k) + 1
      J <- J_block(a, b)
      grad[i] <- -sum(diag(Iinv %*% J %*% Iinv))
    }
    grad
  }
  
  # Equality constraints: column sums = v_vec; total sum = 1
  eq_constraints <- function(x) {
    M <- matrix(x, nrow = k, ncol = k_prime)
    c(colSums(M) - v_vec,
      sum(M) - 1)
  }
  
  # Jacobian of the equality constraints
  eq_jacobian <- function(x) {
    Jac <- matrix(0, k_prime + 1, k * k_prime)
    for (b in seq_len(k_prime)) {
      Jac[b, ((b - 1) * k + 1):(b * k)] <- 1
    }
    Jac[k_prime + 1, ] <- 1
    Jac
  }
  
  # Initial uniform guess: π_{a,b} = v_b / k
  x0 <- as.vector(outer(rep(1 / k, k), v_vec))
  
  # Solve via Augmented Lagrangian with SLSQP inner solver
  res <- nloptr(
    x0            = x0,
    eval_f        = f_objective,
    eval_grad_f   = f_gradient,
    eval_g_eq     = eq_constraints,
    eval_jac_g_eq = eq_jacobian,
    lb            = rep(lb_eps, k * k_prime),
    ub            = rep(1,   k * k_prime),
    opts          = list(
      algorithm   = "NLOPT_LD_AUGLAG",
      local_opts  = list(
        algorithm = "NLOPT_LD_SLSQP",
        xtol_rel  = 1e-8
      ),
      maxeval     = maxeval,
      xtol_rel    = 1e-8,
      print_level = 0
    )
  )
  
  list(
    status  = res$status,
    message = res$message,
    Fval    = res$objective,
    pi_hat  = matrix(res$solution, nrow = k, ncol = k_prime)
  )
}

## Example run for v = (0.5, 0.5)ᵀ
v_vec   <- c(0.5, 0.5) 
opt_pi <- solve_pi_opt(d_mat, sigma2_vec, v_vec)
opt_pi$status    # optimizer status (≥ 1 indicates success)
opt_pi$pi_hat    # estimated optimal weights π̂_{a,b}
opt_pi$Fval      # minimal tr{I(π)⁻¹}


############################################################
## Section 5 – Adaptive sampling with random stopping    ##
############################################################

run_adaptive_stop <- function(theta,
                              sigma2_vec,
                              d_mat,
                              a_init,
                              c_val,
                              v_vec,
                              radius   = 5,
                              lower_bd = 0.1,
                              upper_bd = 10,
                              tol      = 1e-10) {
  
  p        <- length(theta)
  k_prime  <- length(sigma2_vec)
  k        <- ncol(d_mat)
  m0       <- length(a_init)
  
  a_seq <- a_init
  b_seq <- rep(1:k_prime, length.out = m0)
  X_seq <- sample_X_vec(theta, sigma2_vec, a_seq, b_seq, d_mat)
  
  repeat {
    m       <- length(a_seq)
    mle     <- fit_mle_bca(d_mat, a_seq, b_seq, X_seq,
                           radius   = radius,
                           lower_bd = lower_bd,
                           upper_bd = upper_bd,
                           tol      = tol)
    var_hat <- mle$var_hat
    
    info_tot <- matrix(0, p + k_prime, p + k_prime)
    for (j in seq_len(m)) {
      info_tot <- info_tot +
        fisher_single(a_seq[j], b_seq[j], d_mat, var_hat)
    }
    
    Iinv <- solve(info_tot) 
    
    crit <- sum(diag(Iinv[1:p, 1:p])) 
    if (crit <= c_val) {
      psi_hat <- c(mle$theta_hat, mle$var_hat)
      diff    <- psi_hat - c(theta, sigma2_vec)
      D_val   <- as.numeric(diff %*% info_tot %*% diff)
      break
    }
    
    b_next <- ((m) %% k_prime) + 1
    scores <- numeric(k)
    for (a in seq_len(k)) {
      Icand     <- fisher_single(a, b_next, d_mat, var_hat)
      scores[a] <- sum(diag(Iinv %*% Icand %*% Iinv))
    }
    a_next <- which.max(scores)
    
    a_seq <- c(a_seq, a_next)
    b_seq <- c(b_seq, b_next)
    X_seq <- c(X_seq,
               sample_X_ab(theta, sigma2_vec, a_next, b_next, d_mat))
  }
  
  list(
    tau       = m,
    D         = D_val,
    theta_hat = mle$theta_hat,
    var_hat   = mle$var_hat,
    a_seq     = a_seq,
    b_seq     = b_seq,
    X_seq     = X_seq
  )
}

############################################################
## Example & empirical vs. optimal comparison            ##
############################################################

# 1. set parameters and run adaptive sampling
c_val  <- 0.01
v_vec  <- c(0.5, 0.5)
a_init <- rep(c(1, 2, 3),1)

res <- run_adaptive_stop(theta, sigma2_vec, d_mat, a_init, c_val, v_vec)

# 2. display stopping results
cat("Stopping time τ =", res$tau, "\n")
cat("MLE θ at stop =", res$theta_hat, "\n")
cat("MLE γ at stop =", res$var_hat, "\n")
cat("Mahalanobis D =", res$D, "\n\n")

# 3. empirical allocation π_emp
k       <- ncol(d_mat)
k_prime <- length(sigma2_vec)
counts  <- table(
  factor(res$a_seq, 1:k),
  factor(res$b_seq, 1:k_prime)
)
pi_emp <- counts / res$tau

# 4. optimal allocation π_opt
opt_pi <- solve_pi_opt(d_mat, sigma2_vec, v_vec)
pi_opt <- opt_pi$pi_hat
F_opt  <- opt_pi$Fval

# 5. compare allocations
cat("π_emp - π_opt:\n")
print(pi_emp - pi_opt)

# 6. evaluate objective at π_emp
I_from_pi_mat <- function(pi_mat) {
  p <- nrow(d_mat)
  I <- matrix(0, p + k_prime, p + k_prime)
  for (a in 1:k) for (b in 1:k_prime) {
    w <- pi_mat[a, b]
    if (w > 0) {
      d <- d_mat[, a, drop = FALSE]
      I[1:p,1:p]     <- I[1:p,1:p] + w/sigma2_vec[b] * (d %*% t(d))
      I[p+b,p+b]     <- I[p+b,p+b] + w/(2*sigma2_vec[b]^2)
    }
  }
  I
}
I_emp <- I_from_pi_mat(pi_emp)
F_emp <- sum(diag(solve(I_emp)))

# 7. compare objective values
cat("\nF(pi_emp) =", F_emp, "\n")
cat("F(pi_opt)  =", F_opt, "\n")
cat("Difference =", F_emp - F_opt, "\n")



############################################################
## Section 6 – Monte Carlo simulation of (τ, D)          ##
############################################################
############################################################
## Section 6 – Monte-Carlo study for several c values     ##
############################################################
c_grid <- c(0.1,0.05, 0.025)   # thresholds to investigate
N      <- 10000               # replications per c
a_init <- rep(c(1, 2, 3), 2)    # initial designs
df     <- length(theta) + length(sigma2_vec)
p      <- length(theta)
k_prime<- length(sigma2_vec)
n_bins <- 20                    # bins for both histograms

## ---------- deterministic τ_det is c-dependent ----------------------------
get_tau_det <- function(c_val) {
  opt_pi <- solve_pi_opt(d_mat, sigma2_vec, v_vec)
  pi_opt <- opt_pi$pi_hat
  I_det  <- matrix(0, df, df)
  for (a in seq_len(ncol(d_mat)))
    for (b in seq_len(k_prime))
      I_det <- I_det + pi_opt[a, b] *
    fisher_single(a, b, d_mat, sigma2_vec)
  tr_Iinv_theta <- sum(diag(solve(I_det)[1:p, 1:p]))
  ceiling(tr_Iinv_theta / c_val)
}

## ---------- loop over c values -------------------------------------------
set.seed(2025)
for (c_val in c_grid) {
  
  ## 1. Monte-Carlo sampling of (τc , D)
  tau_vec <- integer(N)
  D_vec   <- numeric(N)
  for (i in seq_len(N)) {
    res        <- run_adaptive_stop(theta, sigma2_vec, d_mat,
                                    a_init, c_val, v_vec)
    tau_vec[i] <- res$tau
    D_vec[i]   <- res$D
  }
  
  ## 2. Save results (optional)
  write.csv(
    data.frame(run = seq_len(N), tau = tau_vec, D = D_vec),
    file = sprintf("tau_D_c%g.csv", c_val),
    row.names = FALSE
  )
  
  ## 3. deterministic stopping time
  tau_det <- get_tau_det(c_val)
  
