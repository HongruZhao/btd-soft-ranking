# -1) source the external R file
source("basic_functions.R")

###############################################################################
### 0) Preliminary Setup
###############################################################################
# We'll define expit if not already:
if (!exists("expit")) {
  expit <- function(z) 1/(1 + exp(-z))
}

###############################################################################
### 1) Parameter Settings
###############################################################################
set.seed(2025)

p        <- 10        # dimension of theta
reg      <- 4         # "4-regular" or "full"
T_final  <- 200       # final sample size per replicate
M        <- 1000        # replicate count
v_vec    <- c(0.3, 0.4, 0.3)
k_B      <- length(v_vec)
epss     <- 1e-8

# For building H in "Case1" (soft ranking):
eta_val   <- 0.1
delta_val <- 0.1

cat("p=", p, "k_B=", k_B, "\neta_val=", eta_val, " delta_val=", delta_val,"\n")

###############################################################################
### 2) Build connected_pair_matrix => e.g. "full" or "regular"
###############################################################################
connected_pair_matrix <- init_edge_generate_RA(
  p      = p,
  method = "full",   # or "regular", "min"
  reg    = reg
)

cat("connected_pair_matrix dim=", dim(connected_pair_matrix), "\n")

###############################################################################
### 3) Generate MST-based initial data => data_random_generator_MST_BTD(...)
###############################################################################
initial_out <- data_random_generator_MST_BTD(
  p                    = p,
  connected_pair_matrix= connected_pair_matrix,
  k_B                  = k_B,
  theta_min            = -2,
  theta_max            =  2,
  phi_min              = -3,
  phi_max              =  0,
  par_true             = NULL   # random
)

df_data_init <- initial_out$data_random
cat("Initial data size =", nrow(df_data_init), "\n")

theta_true <- initial_out$theta - mean(initial_out$theta)
phi_true   <- initial_out$phi
par_true_full <- c(theta_true, phi_true)

cat("True param =>\n")
print(par_true_full)
sum(theta_true^2)

###############################################################################
### 4) Define dynamic H for Case1 and Case2
###############################################################################
H_dynamic_fun_case1 <- function(par_mle) {
  build_soft_ranking_H(
    par_true = par_mle,
    p        = p,
    eta      = eta_val,
    delta    = delta_val
  )
}
H_dynamic_fun_case2 <- function(par_mle) {
  build_H(
    p       = p,
    k_B     = k_B,
    eta_theta = 1,
    eta       = eta_val
  )
}

###############################################################################
### 5) Solve pi_soft_limit for "Case1"
###############################################################################
H_case1_true <- build_soft_ranking_H(
  par_true = par_true_full,
  p        = p,
  eta      = eta_val,
  delta    = delta_val
)
pi_init_case1 <- init_pi_mat(nrow(connected_pair_matrix), v_vec)
res_case1_soft <- solve_pi_nloptr(
  par   = par_true_full,
  connected_pair_matrix = connected_pair_matrix,
  k_B   = k_B,
  H     = H_case1_true,
  v_vec = v_vec,
  pi_init= pi_init_case1,
  maxeval= 2000
)
pi_soft_limit_case1 <- res_case1_soft$pi
cat("Case1 => pi_soft_limit solved.\n")

###############################################################################
### 6) We'll define "run_one_replicate()" => run 4 selection approaches
###############################################################################
R_theta=4
run_one_replicate <- function(R_theta) {
  # 1) Case1
  res_case1 <- multiple_rounds_adaptive_dynamicH(
    T              = T_final,
    df_data_init   = df_data_init,
    pi_n_init      = NULL,
    connected_pair_matrix= connected_pair_matrix,
    p= p,  k_B= k_B,
    R= R_theta,                    # sum(theta^2)<=4 constraint
    phi_lower= -3, phi_upper= 3,
    init_par= NULL,
    method= "NLOPT_LD_SLSQP",
    maxeval= 200,
    verbose= FALSE,
    v= v_vec,
    vec_b_init= integer(0),
    par_true= par_true_full,
    H_dynamic_fun= H_dynamic_fun_case1
  )
  
  # 2) Case2
  res_case2 <- multiple_rounds_adaptive_dynamicH(
    T              = T_final,
    df_data_init   = df_data_init,
    pi_n_init      = NULL,
    connected_pair_matrix= connected_pair_matrix,
    p= p,  k_B= k_B,
    R= 4,
    phi_lower= -3, phi_upper= 3,
    init_par= NULL,
    method= "NLOPT_LD_SLSQP",
    maxeval= 200,
    verbose= FALSE,
    v= v_vec,
    vec_b_init= integer(0),
    par_true= par_true_full,
    H_dynamic_fun= H_dynamic_fun_case2
  )
  
  # 3) Uniform
  res_unif <- multiple_rounds_uniform(
    T              = T_final,
    df_data_init   = df_data_init,
    pi_n_init      = NULL,
    connected_pair_matrix= connected_pair_matrix,
    p= p,  k_B= k_B,
    R= 4,
    phi_lower= -3, phi_upper= 3,
    init_par= NULL,
    method= "NLOPT_LD_SLSQP",
    maxeval= 200,
    verbose= FALSE,
    v= v_vec,
    vec_b_init= integer(0),
    par_true= par_true_full
  )
  
  # 4) Uncertainty
  res_uncert <- multiple_rounds_uncertainty(
    T              = T_final,
    df_data_init   = df_data_init,
    pi_n_init      = NULL,
    connected_pair_matrix= connected_pair_matrix,
    p= p,  k_B= k_B,
    R= 4,
    phi_lower= -3, phi_upper= 3,
    init_par= NULL,
    method= "NLOPT_LD_SLSQP",
    maxeval= 200,
    verbose= FALSE,
    v= v_vec,
    vec_b_init= integer(0),
    par_true= par_true_full
  )
  
  list(case1   = res_case1,
       case2   = res_case2,
       unif    = res_unif,
       uncert  = res_uncert)
}

###############################################################################
### 7) We define metrics => (A) F_n difference, (B) Kendall's tau,
###    (C) hetero_soft_ranking_loss, (D) hetero_soft_ranking_loss * n
###############################################################################
compute_diff_vec_sameH <- function(
    pi_tensor, pi_soft_limit, par_true, H_same,
    connected_pair_matrix, k_B
) {
  n_iter <- dim(pi_tensor)[3]
  diff_vec <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    val_i <- F_pi_par(
      pi_mat = pi_tensor[,, i],
      par    = par_true,
      connected_pair_matrix= connected_pair_matrix,
      k_B= k_B,
      H= H_same
    )
    val_soft <- F_pi_par(
      pi_mat = pi_soft_limit,
      par    = par_true,
      connected_pair_matrix= connected_pair_matrix,
      k_B= k_B,
      H= H_same
    )
    diff_vec[i] <- val_i - val_soft
  }
  diff_vec
}

compute_hsr_loss_vec <- function(
    par_mle_tensor,
    p, k_B,
    theta_true, phi_true,
    delta=0.1, eta=0.1
) {
  n_iter <- ncol(par_mle_tensor)
  out_vec <- numeric(n_iter)
  
  for (i in seq_len(n_iter)) {
    par_est   <- par_mle_tensor[, i]
    est_theta <- par_est[1:p]
    est_phi   <- par_est[(p+1):(p+k_B)]
    val_loss  <- hetero_soft_ranking_loss(
      theta      = est_theta,
      theta_star = theta_true,
      gamma      = est_phi,
      gamma_star = phi_true,
      delta      = delta,
      eta        = eta
    )
    out_vec[i] <- val_loss
  }
  out_vec
}

###############################################################################
### 8) Main Monte Carlo loop => M=10 replicates
###############################################################################
n_init <- nrow(df_data_init)
max_steps <- T_final - n_init
n_iter <- max_steps + 1

diff_mat_case1  <- matrix(NA, M, n_iter)
diff_mat_case2  <- matrix(NA, M, n_iter)
diff_mat_unif   <- matrix(NA, M, n_iter)
diff_mat_uncert <- matrix(NA, M, n_iter)

tau_mat_case1   <- matrix(NA, M, n_iter)
tau_mat_case2   <- matrix(NA, M, n_iter)
tau_mat_unif    <- matrix(NA, M, n_iter)
tau_mat_uncert  <- matrix(NA, M, n_iter)

hsr_mat_case1   <- matrix(NA, M, n_iter)
hsr_mat_case2   <- matrix(NA, M, n_iter)
hsr_mat_unif    <- matrix(NA, M, n_iter)
hsr_mat_uncert  <- matrix(NA, M, n_iter)

hsrN_mat_case1  <- matrix(NA, M, n_iter)
hsrN_mat_case2  <- matrix(NA, M, n_iter)
hsrN_mat_unif   <- matrix(NA, M, n_iter)
hsrN_mat_uncert <- matrix(NA, M, n_iter)

cat("\n=== Starting M=", M, "replicates...\n")
for (m_i in seq_len(M)) {
  cat("Replicate =", m_i,"\n")
  out_m <- run_one_replicate(R_theta)
  
  # =========== CASE1 ============
  diff_case1 <- compute_diff_vec_sameH(
    pi_tensor      = out_m$case1$pi_tensor,
    pi_soft_limit  = pi_soft_limit_case1,
    par_true       = par_true_full,
    H_same         = H_case1_true,
    connected_pair_matrix= connected_pair_matrix,
    k_B= k_B
  )
  diff_mat_case1[m_i, ] <- diff_case1
  
  tau_case1 <- ranking_eval_kendall(
    par_mle_tensor = out_m$case1$par_mle_tensor,
    p= p,
    theta_true= par_true_full[1:p]
  )
  tau_mat_case1[m_i, ] <- tau_case1
  
  hsr_case1 <- compute_hsr_loss_vec(
    par_mle_tensor= out_m$case1$par_mle_tensor,
    p= p, k_B= k_B,
    theta_true= par_true_full[1:p],
    phi_true  = par_true_full[(p+1):(p+k_B)],
    delta= delta_val,
    eta= eta_val
  )
  hsr_mat_case1[m_i, ] <- hsr_case1
  
  sample_size_vec <- (n_init + seq_len(n_iter) - 1)
  hsrN_case1 <- hsr_case1 * sample_size_vec
  hsrN_mat_case1[m_i, ] <- hsrN_case1
  
  # =========== CASE2 ============
  diff_case2 <- compute_diff_vec_sameH(
    out_m$case2$pi_tensor, pi_soft_limit_case1,
    par_true_full, H_case1_true, connected_pair_matrix, k_B
  )
  diff_mat_case2[m_i, ] <- diff_case2
  
  tau_case2 <- ranking_eval_kendall(
    out_m$case2$par_mle_tensor, p, par_true_full[1:p]
  )
  tau_mat_case2[m_i, ] <- tau_case2
  
  hsr_case2 <- compute_hsr_loss_vec(
    out_m$case2$par_mle_tensor,
    p, k_B,
    par_true_full[1:p], par_true_full[(p+1):(p+k_B)],
    delta= delta_val, eta= eta_val
  )
  hsr_mat_case2[m_i, ] <- hsr_case2
  hsrN_case2 <- hsr_case2 * sample_size_vec
  hsrN_mat_case2[m_i, ] <- hsrN_case2
  
  # =========== UNIFORM ============
  diff_unif <- compute_diff_vec_sameH(
    out_m$unif$pi_tensor, pi_soft_limit_case1,
    par_true_full, H_case1_true, connected_pair_matrix, k_B
  )
  diff_mat_unif[m_i, ] <- diff_unif
  
  tau_unif <- ranking_eval_kendall(
    out_m$unif$par_mle_tensor, p, par_true_full[1:p]
  )
  tau_mat_unif[m_i, ] <- tau_unif
  
  hsr_unif <- compute_hsr_loss_vec(
    out_m$unif$par_mle_tensor,
    p, k_B,
    par_true_full[1:p], par_true_full[(p+1):(p+k_B)],
    delta= delta_val, eta= eta_val
  )
  hsr_mat_unif[m_i, ] <- hsr_unif
  hsrN_unif <- hsr_unif * sample_size_vec
  hsrN_mat_unif[m_i, ] <- hsrN_unif
  
  # =========== UNCERTAINTY ========
  diff_uncert <- compute_diff_vec_sameH(
    out_m$uncert$pi_tensor, pi_soft_limit_case1,
    par_true_full, H_case1_true, connected_pair_matrix, k_B
  )
  diff_mat_uncert[m_i, ] <- diff_uncert
  
  tau_uncert <- ranking_eval_kendall(
    out_m$uncert$par_mle_tensor, p, par_true_full[1:p]
  )
  tau_mat_uncert[m_i, ] <- tau_uncert
  
  hsr_uncert <- compute_hsr_loss_vec(
    out_m$uncert$par_mle_tensor,
    p, k_B,
    par_true_full[1:p], par_true_full[(p+1):(p+k_B)],
    delta= delta_val, eta= eta_val
  )
  hsr_mat_uncert[m_i, ] <- hsr_uncert
  hsrN_uncert <- hsr_uncert * sample_size_vec
  hsrN_mat_uncert[m_i, ] <- hsrN_uncert
}

###############################################################################
### 9) Summaries => 2.5%, 50%, 97.5% across replicates for each metric
###############################################################################
Q_lower <- 0.05
Q_upper <- 0.95


get_quantiles <- function(mat) {
  # mat: M x n_iter numeric matrix
  # We want:
  #   row1 => 5% quantile
  #   row2 => mean
  #   row3 => 95% quantile
  
  q5  <- apply(mat, 2, quantile, probs = Q_lower, na.rm = TRUE)
  q95 <- apply(mat, 2, quantile, probs = Q_upper, na.rm = TRUE)
  mu  <- colMeans(mat, na.rm = TRUE)
  
  # Return a 3 x n_iter matrix:
  #   row1 = 5% quantile
  #   row2 = mean
  #   row3 = 95% quantile
  rbind(q5, mu, q95)
}




q_diff_case1  <- get_quantiles(diff_mat_case1)
q_diff_case2  <- get_quantiles(diff_mat_case2)
q_diff_unif   <- get_quantiles(diff_mat_unif)
q_diff_uncert <- get_quantiles(diff_mat_uncert)

q_tau_case1   <- get_quantiles(tau_mat_case1)
q_tau_case2   <- get_quantiles(tau_mat_case2)
q_tau_unif    <- get_quantiles(tau_mat_unif)
q_tau_uncert  <- get_quantiles(tau_mat_uncert)

q_hsr_case1   <- get_quantiles(hsr_mat_case1)
q_hsr_case2   <- get_quantiles(hsr_mat_case2)
q_hsr_unif    <- get_quantiles(hsr_mat_unif)
q_hsr_uncert  <- get_quantiles(hsr_mat_uncert)

q_hsrN_case1  <- get_quantiles(hsrN_mat_case1)
q_hsrN_case2  <- get_quantiles(hsrN_mat_case2)
q_hsrN_unif   <- get_quantiles(hsrN_mat_unif)
q_hsrN_uncert <- get_quantiles(hsrN_mat_uncert)

###############################################################################
### 10) a helper to plot 4 curves => median lines => 2.5..97.5 band
###     plus new argument legend_pos for specifying legend position
###############################################################################
### new
plot_four_approaches <- function(
    xvals,
    q_case1, q_case2, q_unif, q_uncert,  # each is a 3 x n_iter matrix => row1=5%, row2=mean, row3=95%
    ylim_vec       = NULL,               # user-specified y-limits or NULL
    main_title     = "",
    ylab_str       = "",
    legend_labels  = c("HSR loss optimal",
                       expression("l2 loss optimal"),
                       "Uniform sampling",
                       "Uncertainty sampling"),
    line_colors    = c("blue","brown","grey","red"),
    hline0         = TRUE,
    legend_pos     = "topright"
) {
  # xvals: numeric vector, length = n_iter
  # q_caseX: a 3 x n_iter matrix:
  #    row1 => 5%  => lty=3, lwd=1
  #    row2 => mean => lty=1, lwd=2
  #    row3 => 95% => lty=3, lwd=1
  # ylim_vec: if NULL => auto from data; else c(ymin, ymax)
  # legend_labels: char or expression vector of length 4
  # line_colors: color for each approach
  # hline0: draw horizontal line at y=0 if TRUE
  # legend_pos: legend location, e.g. "topright", "bottomleft", etc.
  
  n_iter  <- length(xvals)
  all_vals <- c(q_case1, q_case2, q_unif, q_uncert)
  y_min  <- min(all_vals, na.rm=TRUE)
  y_max  <- max(all_vals, na.rm=TRUE)
  
  # Determine y-limits
  if (is.null(ylim_vec)) {
    ylim_use <- c(y_min, y_max)
  } else {
    ylim_use <- ylim_vec
  }
  
  # Base plot
  plot(
    x = xvals,
    y = rep(NA, n_iter),  # placeholders
    type = "n",
    main = main_title,
    xlab = "Sample Size n",
    ylab = ylab_str,
    ylim = ylim_use
  )
  
  # Optionally a horizontal line at 0
  if (hline0) {
    abline(h = 0, col="black", lwd=1)
  }
  
  # Helper to add lines for each approach
  add_lines <- function(q_mat, col_main) {
    # row1 => lower 5% => lty=3, lwd=1
    lines(xvals, q_mat[1,], col=col_main, lty=3, lwd=1)
    # row2 => mean => lty=1, lwd=2
    lines(xvals, q_mat[2,], col=col_main, lty=1, lwd=2)
    # row3 => upper 95% => lty=3, lwd=1
    lines(xvals, q_mat[3,], col=col_main, lty=3, lwd=1)
  }
  
  # Draw lines for each approach
  add_lines(q_case1, line_colors[1])
  add_lines(q_case2, line_colors[2])
  add_lines(q_unif,  line_colors[3])
  add_lines(q_uncert,line_colors[4])
  
  # Legend => no box => bty="n"
  legend(
    legend_pos,
    legend = legend_labels,
    col    = line_colors,
    lty    = 1,      # show the main line type (mean) in the legend
    lwd    = 2,
    bty    = "n"     # remove border box
  )
}

### old
plot_four_approaches <- function(
    xvals,
    ylim_vec,
    q_case1, q_case2, q_unif, q_uncert,
    main_title      = "",
    ylab_str        = "",
    legend_labels   = c(
      "HSR loss optimal",
      expression("l2 loss optimal"),  # \ell_2
      "Uniform sampling",
      "Uncertainty sampling"
    ),
    line_colors     = c("blue","brown","grey","red"),
    hline0          = TRUE,
    legend_pos      = "topright"  # new argument
) {
  n_iter <- length(xvals)
  all_vals <- c(q_case1, q_case2, q_unif, q_uncert)
  y_min <- min(all_vals, na.rm=TRUE)
  y_max <- max(all_vals, na.rm=TRUE)
  
  plot(
    xvals,
    rep(NA, n_iter),
    #    type="n",
    main= main_title,
    xlab= "Sample Size n",
    ylab= ylab_str,
    ylim= ylim_vec
  )
  
  draw_band <- function(x, low_vec, high_vec, col_band) {
    polygon(c(x, rev(x)), c(low_vec, rev(high_vec)), border=NA, col=col_band)
  }
  
  # approach1 => case1
  draw_band(xvals, q_case1[1,], q_case1[3,], adjustcolor(line_colors[1], alpha.f=0.2))
  lines(xvals, q_case1[2,], col=line_colors[1], lty=1, lwd=1.5)
  
  # approach2 => case2
  draw_band(xvals, q_case2[1,], q_case2[3,], adjustcolor(line_colors[2], alpha.f=0.2))
  lines(xvals, q_case2[2,], col=line_colors[2], lty=2, lwd=2)
  
  # approach3 => unif
  draw_band(xvals, q_unif[1,], q_unif[3,], adjustcolor(line_colors[3], alpha.f=0.2))
  lines(xvals, q_unif[2,], col=line_colors[3], lty=3, lwd=2)
  
  # approach4 => uncert
  draw_band(xvals, q_uncert[1,], q_uncert[3,], adjustcolor(line_colors[4], alpha.f=0.2))
  lines(xvals, q_uncert[2,], col=line_colors[4], lty=4, lwd=2)
  
  if (hline0) {
    abline(h=0, col="black", lwd=1)
  }
  
  legend(
    legend_pos,
    legend = legend_labels,
    col    = line_colors,
    lty    = 1:4,
    lwd    = 2,
    bty    = "n"     # remove border box
  )
}

###############################################################################
### 11) produce 4 separate plots:
###   (1) F_n difference
###   (2) Kendall's tau
###   (3) hetero_soft_ranking_loss
###   (4) hetero_soft_ranking_loss * n
###############################################################################
xvals <- seq_len(n_iter)+T_final-n_iter

#xvals[-(1:12)]

# (1) F_n difference
png("plot_Fn_diff.png", width=1200/3*2, height=900/3*2, res=120)
plot_four_approaches(
  xvals,
  ylim_vec=c(0,1000),
  q_case1 = q_diff_case1,
  q_case2 = q_diff_case2,
  q_unif  = q_diff_unif,
  q_uncert= q_diff_uncert,
  main_title = NULL,
  ylab_str   = expression(F[n]),
  legend_labels= c(
    "HSR loss optimal",
    expression("l2 loss optimal"),
    "Uniform sampling",
    "Uncertainty sampling"
  ),
  line_colors = c("blue","brown","darkgreen","red"),
  hline0      = TRUE,
  legend_pos  = "topright"  # specify legend position
)
dev.off()



cat("\nAll done! Created 4 PNG plots:\n",
    "1) plot_Fn_diff.png\n")

