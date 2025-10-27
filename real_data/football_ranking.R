# Helper: locate project root (assumes script is in "real data")
find_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- NULL
  idx <- grep(file_arg, args)
  if (length(idx) > 0) {
    script_path <- normalizePath(sub(file_arg, "", args[idx]))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    script_path <- normalizePath(sys.frames()[[1]]$ofile)
  }
  if (is.null(script_path)) {
    return(normalizePath(".."))
  }
  normalizePath(file.path(dirname(script_path), ".."))
}

project_root <- find_project_root()

results_path <- file.path(project_root, "real data", "results.csv")
raw_results <- read.csv(results_path, stringsAsFactors = FALSE)
raw_results$date <- as.Date(raw_results$date)

# ---- Minimal BTD utilities (no external tidyverse dependencies) ----
if (!requireNamespace("nloptr", quietly = TRUE)) {
  stop("Package 'nloptr' is required. Please install it before running this script.")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required. Please install it before running this script.")
}

# Hyper-parameters (user adjustable)
MC_REPS <- 100L
MAX_TEAMS_TO_KEEP <- 64L
MIN_MATCHES_THRESHOLD <- 30L
TOTAL_ADDITIONAL_SAMPLES <- 200L
WITH_REPLACEMENT <- TRUE

mc_override_env <- Sys.getenv("MC_REPS_OVERRIDE")
if (nzchar(mc_override_env)) {
  MC_REPS <- as.integer(mc_override_env)
  if (is.na(MC_REPS) || MC_REPS < 1L) {
    stop("Environment variable MC_REPS_OVERRIDE must be a positive integer.")
  }
}

with_replacement_env <- Sys.getenv("WITH_REPLACEMENT_OVERRIDE")
if (nzchar(with_replacement_env)) {
  WITH_REPLACEMENT <- !(tolower(with_replacement_env) %in% c("0", "false", "no"))
}

total_samples_override <- Sys.getenv("TOTAL_ADDITIONAL_SAMPLES_OVERRIDE")
if (nzchar(total_samples_override)) {
  TOTAL_ADDITIONAL_SAMPLES <- as.integer(total_samples_override)
  if (is.na(TOTAL_ADDITIONAL_SAMPLES) || TOTAL_ADDITIONAL_SAMPLES < 0L) {
    stop("Environment variable TOTAL_ADDITIONAL_SAMPLES_OVERRIDE must be a non-negative integer.")
  }
}

MASTER_SEED <- 3077L
N_PLOT_INITIAL <- 20L
set.seed(MASTER_SEED)

script_start_time <- proc.time()

expit <- function(z) {
  1 / (1 + exp(-z))
}

construct_U <- function(p) {
  U <- matrix(0, nrow = p, ncol = p - 1)
  if (p <= 1) return(U)
  for (k in 1:(p - 1)) {
    U[1:k, k] <- 1
    U[k + 1, k] <- -k
    norm_k <- sqrt(k * (k + 1))
    U[, k] <- U[, k] / norm_k
  }
  U
}

unit_vector <- function(p, i) {
  v <- numeric(p)
  v[i] <- 1
  v
}

fisher_logit_rank_single <- function(theta, phi, pair_idx) {
  i <- pair_idx[1]
  j <- pair_idx[2]
  dtheta <- theta[i] - theta[j]
  v <- exp(phi)
  csh <- cosh(dtheta / 2)
  ssh <- sinh(dtheta / 2)
  denom <- (v + 2 * csh)^2
  I11 <- (v * csh + 2) / (2 * denom)
  I12 <- -(v * ssh) / denom
  I22 <- (2 * v * csh) / denom
  matrix(c(I11, I12, I12, I22), nrow = 2, byrow = TRUE)
}

fisher_info_block <- function(theta_vals, phi_val, pair_a, k_B, b_idx) {
  p <- length(theta_vals)
  I2 <- fisher_logit_rank_single(theta_vals, phi_val, pair_a)
  e_i <- unit_vector(p, pair_a[1])
  e_j <- unit_vector(p, pair_a[2])
  eij <- e_i - e_j
  d <- p + k_B
  M <- matrix(0, nrow = d, ncol = d)
  M[1:p, 1:p] <- I2[1, 1] * (eij %*% t(eij))
  M[1:p, p + b_idx] <- I2[1, 2] * eij
  M[p + b_idx, 1:p] <- I2[2, 1] * t(eij)
  M[p + b_idx, p + b_idx] <- I2[2, 2]
  M
}

moore_penrose_custom <- function(I_pi, p, k_B, eps = 1e-8) {
  d <- p + k_B
  if (p <= 1) return(matrix(0, d, d))
  U_p <- construct_U(p)
  Q <- matrix(0, nrow = d, ncol = d - 1)
  Q[1:p, 1:(p - 1)] <- U_p
  if (k_B > 0) {
    for (i in seq_len(k_B)) {
      Q[p + i, (p - 1) + i] <- 1
    }
  }
  I_low <- t(Q) %*% I_pi %*% Q
  inv_low <- solve(I_low + eps * diag(d - 1))
  Q %*% inv_low %*% t(Q)
}

build_H <- function(p, k_B, eta_theta = 1, eta = 0.5) {
  d <- p + k_B
  diag(c(rep(eta_theta, p), rep(eta, k_B)), nrow = d, ncol = d)
}

build_soft_ranking_H <- function(par_est, p, eta = 0.5, delta = 1) {
  theta_star <- par_est[1:p]
  d <- length(par_est)
  k_B <- d - p
  H_mat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      diff_ij <- abs(theta_star[i] - theta_star[j]) + delta
      val <- diff_ij^2
      eij <- numeric(p)
      eij[i] <- 1
      eij[j] <- -1
      H_mat <- H_mat + val * (eij %*% t(eij))
    }
  }
  H_full <- diag(d)
  H_full[1:p, 1:p] <- H_mat / max(p, 1)
  if (k_B > 0) {
    for (idx in 1:k_B) {
      H_full[p + idx, p + idx] <- eta
    }
  }
  H_full
}

select_action_restricted <- function(
    par,
    H,
    pi_n,
    b_idx,
    allowed_pairs,
    connected_pair_matrix,
    p,
    k_B,
    eps = 1e-8) {
  if (length(allowed_pairs) == 0) {
    return(NA_integer_)
  }
  d <- p + k_B
  k <- nrow(connected_pair_matrix)
  I_pi <- matrix(0, nrow = d, ncol = d)
  for (a_row in seq_len(k)) {
    for (b_col in seq_len(k_B)) {
      weight <- pi_n[a_row, b_col]
      if (weight <= 0) next
      I_pi <- I_pi + weight * fisher_info_block(
        theta_vals = par[1:p],
        phi_val = par[(p + 1):(p + k_B)][b_col],
        pair_a = connected_pair_matrix[a_row, ],
        k_B = k_B,
        b_idx = b_col
      )
    }
  }
  MPinv <- moore_penrose_custom(I_pi, p, k_B, eps)
  best_val <- -Inf
  best_pair <- allowed_pairs[1]
  for (a_row in allowed_pairs) {
    I_ab <- fisher_info_block(
      theta_vals = par[1:p],
      phi_val = par[(p + 1):(p + k_B)][b_idx],
      pair_a = connected_pair_matrix[a_row, ],
      k_B = k_B,
      b_idx = b_idx
    )
    M_ab <- MPinv %*% I_ab %*% MPinv
    val_ab <- sum(diag(H %*% M_ab))
    if (val_ab > best_val) {
      best_val <- val_ab
      best_pair <- a_row
    }
  }
  best_pair
}

three_outcome_entropy <- function(
    pair_idx,
    b_idx,
    par_mle,
    connected_pair_matrix,
    p,
    k_B) {
  teams_pair <- connected_pair_matrix[pair_idx, ]
  theta_hat <- par_mle[1:p]
  phi_hat <- par_mle[(p + 1):(p + k_B)]
  delta_ij <- theta_hat[teams_pair[1]] - theta_hat[teams_pair[2]]
  f_plus <- expit(delta_ij)
  f_minus <- expit(-delta_ij)
  v_val <- exp(phi_hat[b_idx])
  sqrt_part <- sqrt(f_plus * f_minus)
  denom <- 1 + v_val * sqrt_part
  probs <- c(f_plus / denom, (v_val * sqrt_part) / denom, f_minus / denom)
  probs_pos <- probs[probs > 1e-15]
  -sum(probs_pos * log(probs_pos))
}

neg_loglik_btd <- function(par, df_data, p, k_B) {
  theta <- par[1:p]
  phi <- par[(p + 1):(p + k_B)]
  Xmat <- as.matrix(df_data[, seq_len(p), drop = FALSE])
  Yvec <- df_data$Y
  catMat <- as.matrix(df_data[, (p + 2):(p + 1 + k_B), drop = FALSE])
  bIndex <- catMat %*% seq_len(k_B)
  DeltaVec <- as.numeric(Xmat %*% theta)
  fVec <- expit(DeltaVec)
  oneMinus_f <- expit(-DeltaVec)
  phiVec <- phi[bIndex]
  vVec <- exp(phiVec)
  sqrtPart <- sqrt(fVec * oneMinus_f)
  denom <- 1 + vVec * sqrtPart
  p1 <- fVec / denom
  p2 <- oneMinus_f / denom
  p3 <- (vVec * sqrtPart) / denom
  Iplus <- (Yvec == 1L)
  Iminus <- (Yvec == -1L)
  Izero <- (Yvec == 0L)
  probs <- p1 * Iplus + p2 * Iminus + p3 * Izero
  -sum(log(probs))
}

grad_neg_loglik_btd <- function(par, df_data, p, k_B) {
  theta <- par[1:p]
  phi <- par[(p + 1):(p + k_B)]
  Xmat <- as.matrix(df_data[, seq_len(p), drop = FALSE])
  Yvec <- df_data$Y
  catMat <- as.matrix(df_data[, (p + 2):(p + 1 + k_B), drop = FALSE])
  bIndex <- catMat %*% seq_len(k_B)
  DeltaVec <- as.numeric(Xmat %*% theta)
  eNegDelta <- exp(-DeltaVec)
  phiVec <- phi[bIndex]
  vVec <- exp(phiVec)
  DenVec <- 1 + eNegDelta + vVec * exp(-DeltaVec / 2)
  cXVec <- (1 - Yvec) / 2
  dDen_dx <- -eNegDelta - 0.5 * vVec * exp(-DeltaVec / 2)
  dNLL_dDelta <- cXVec + dDen_dx / DenVec
  grad_theta <- as.numeric(t(Xmat) %*% dNLL_dDelta)
  i0 <- 1 - (Yvec^2)
  dDen_dphi <- exp(-DeltaVec / 2) * vVec
  dNLL_dphi <- (dDen_dphi / DenVec) - i0
  grad_phi <- as.numeric(rowsum(dNLL_dphi, group = bIndex))
  c(grad_theta, grad_phi)
}

hetero_soft_ranking_loss <- function(
    theta,
    theta_star,
    gamma,
    gamma_star,
    delta = 1,
    eta = 1) {
  p <- length(theta)
  loss_sum <- 0
  if (p >= 2) {
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        d_ij <- theta[i] - theta[j]
        d_star <- theta_star[i] - theta_star[j]
        sign_star <- if (d_star > 0) 1 else if (d_star < 0) -1 else 1
        u_val <- abs(d_star) + delta
        hC <- function(x) exp(x) - x - 1
        loss_sum <- loss_sum +
          0.5 * (hC(-((d_ij - d_star) * sign_star * u_val)) +
                   hC(-(((-d_ij) - (-d_star)) * (-sign_star) * u_val)))
      }
    }
  }
  ranking_term <- loss_sum / max(p, 1)
  gamma_vec <- as.numeric(gamma)
  gamma_star_vec <- as.numeric(gamma_star)
  ranking_term + eta * sum((gamma_vec - gamma_star_vec)^2)
}

fit_btd_mle <- function(
    df_data,
    p,
    k_B,
    R = 3,
    phi_lower = -3,
    phi_upper = 3,
    init_par = NULL,
    method = "NLOPT_LD_SLSQP",
    maxeval = 1000,
    verbose = FALSE) {
  nParam <- p + k_B
  if (is.null(init_par)) {
    init_theta <- runif(p, min = -1, max = 1)
    nt <- sum(init_theta^2)
    if (nt > R^2) {
      init_theta <- init_theta * (R / sqrt(nt))
    }
    init_theta <- init_theta - mean(init_theta)
    init_phi <- runif(k_B, min = phi_lower, max = phi_upper)
    init_par <- c(init_theta, init_phi)
  } else if (length(init_par) != nParam) {
    stop("init_par must have length p + k_B")
  }

  hist_list <- list(vals = numeric(0))
  obj_fn <- function(x) {
    val <- neg_loglik_btd(x, df_data, p, k_B)
    hist_list$vals <- c(hist_list$vals, val)
    if (verbose) {
      cat("Objective:", val, "\n")
    }
    val
  }
  grad_fn <- function(x) {
    grad_neg_loglik_btd(x, df_data, p, k_B)
  }
  eq_fn <- function(x) {
    sum(x[1:p])
  }
  eq_grad <- function(x) {
    c(rep(1, p), rep(0, k_B))
  }
  ineq_fn <- function(x) {
    sum(x[1:p]^2) - R^2
  }
  ineq_grad <- function(x) {
    c(2 * x[1:p], rep(0, k_B))
  }
  lb <- c(rep(-Inf, p), rep(phi_lower, k_B))
  ub <- c(rep(Inf, p), rep(phi_upper, k_B))
  res <- nloptr::nloptr(
    x0 = init_par,
    eval_f = obj_fn,
    eval_grad_f = grad_fn,
    eval_g_eq = eq_fn,
    eval_jac_g_eq = eq_grad,
    eval_g_ineq = ineq_fn,
    eval_jac_g_ineq = ineq_grad,
    lb = lb,
    ub = ub,
    opts = list(
      algorithm = method,
      xtol_rel = 1e-8,
      maxeval = maxeval,
      print_level = if (verbose) 1 else 0
    )
  )
  res$obj_vals <- hist_list$vals
  res
}

# Summarise tournaments
years_all <- as.integer(format(raw_results$date, "%Y"))
split_years <- split(years_all, raw_results$tournament)
matches_count <- sort(table(raw_results$tournament), decreasing = TRUE)
first_year <- vapply(split_years, function(v) min(v, na.rm = TRUE), integer(1))
last_year <- vapply(split_years, function(v) max(v, na.rm = TRUE), integer(1))
tournament_summary <- data.frame(
  tournament = names(matches_count),
  matches = as.integer(matches_count[names(matches_count)]),
  first_year = first_year[names(matches_count)],
  last_year = last_year[names(matches_count)],
  row.names = NULL,
  stringsAsFactors = FALSE
)
cat("\nDistinct tournaments with match counts (top 25):\n")
print(head(tournament_summary, 25))

# Keep last 50 years
max_year <- max(years_all, na.rm = TRUE)
cutoff_year <- max_year - 49
recent_idx <- which(years_all >= cutoff_year)
results_recent <- raw_results[recent_idx, , drop = FALSE]

# Target tournaments (indices 1-10, 15,16,17,19,21 from summary)
target_tournaments <- c(
  "Friendly",
  "FIFA World Cup qualification",
  "UEFA Euro qualification",
  "African Cup of Nations qualification",
  "FIFA World Cup",
  "Copa América",
  "AFC Asian Cup qualification",
  "African Cup of Nations",
  "UEFA Nations League",
  "CECAFA Cup",
  "AFC Asian Cup",
  "Gold Cup",
  "Gulf Cup",
  "UEFA Euro",
  "COSAFA Cup"
)

results_major <- results_recent[results_recent$tournament %in% target_tournaments, , drop = FALSE]
results_major$home_score <- as.integer(results_major$home_score)
results_major$away_score <- as.integer(results_major$away_score)

if (nrow(results_major) == 0L) {
  stop("No matches remain after filtering to target tournaments.")
}

cat("\nMatches per selected tournament (last 100 years):\n")
print(sort(table(results_major$tournament), decreasing = TRUE))

# Iteratively remove teams with no wins or too few matches
prune_sparse_teams <- function(
    df,
    min_matches = 20L,
    min_wins = 1L,
    min_losses = 1L,
    max_teams = NULL
) {
  current_df <- df
  repeat {
    matches_counts <- table(c(current_df$home_team, current_df$away_team))
    winners <- c(
      current_df$home_team[current_df$home_score > current_df$away_score],
      current_df$away_team[current_df$away_score > current_df$home_score]
    )
    wins_counts <- table(winners)
    losers <- c(
      current_df$home_team[current_df$home_score < current_df$away_score],
      current_df$away_team[current_df$away_score < current_df$home_score]
    )
    losses_counts <- table(losers)
    teams_all <- sort(unique(c(names(matches_counts), names(wins_counts), names(losses_counts))))
    matches_vec <- as.numeric(matches_counts[teams_all])
    matches_vec[is.na(matches_vec)] <- 0
    wins_vec <- as.numeric(wins_counts[teams_all])
    wins_vec[is.na(wins_vec)] <- 0
    losses_vec <- as.numeric(losses_counts[teams_all])
    losses_vec[is.na(losses_vec)] <- 0
    summary_df <- data.frame(
      team = teams_all,
      matches = as.integer(matches_vec),
      wins = as.integer(wins_vec),
      losses = as.integer(losses_vec),
      stringsAsFactors = FALSE
    )
    keep_mask <- summary_df$matches >= min_matches & summary_df$wins >= min_wins
    keep_mask <- keep_mask & (summary_df$losses >= min_losses)
    keep_candidates <- summary_df[keep_mask, , drop = FALSE]
    if (!is.null(max_teams) && nrow(keep_candidates) > max_teams) {
      keep_candidates <- keep_candidates[order(-keep_candidates$matches, -keep_candidates$wins), , drop = FALSE]
      keep_candidates <- keep_candidates[seq_len(max_teams), , drop = FALSE]
    }
    keep_teams <- keep_candidates$team
    filtered_df <- current_df[
      current_df$home_team %in% keep_teams &
        current_df$away_team %in% keep_teams,
      , drop = FALSE
    ]
    if (nrow(filtered_df) == nrow(current_df)) {
      summary_keep <- summary_df[summary_df$team %in% keep_teams, , drop = FALSE]
      return(list(data = filtered_df, summary = summary_keep))
    }
    if (length(keep_teams) < 2L) {
      stop("Filtering removed too many teams. Consider lowering thresholds.")
    }
    current_df <- filtered_df
  }
}

min_matches_threshold <- MIN_MATCHES_THRESHOLD
max_teams_to_keep <- MAX_TEAMS_TO_KEEP
pruned <- prune_sparse_teams(
  results_major,
  min_matches = min_matches_threshold,
  min_wins = 1L,
  min_losses = 1L,
  max_teams = max_teams_to_keep
)
results_pruned <- pruned$data
team_summary <- pruned$summary[order(-pruned$summary$matches), , drop = FALSE]

cat("\nTeam summary after pruning (top 20 by matches):\n")
print(head(team_summary, 20))

if (nrow(results_pruned) == 0L) {
  stop("No matches remain after pruning.")
}

graph_edges_char <- cbind(as.character(results_pruned$home_team), as.character(results_pruned$away_team))
g_results <- igraph::graph_from_edgelist(graph_edges_char, directed = FALSE)
comp_results <- igraph::components(g_results)
if (comp_results$no > 1) {
  giant_idx <- which.max(comp_results$csize)
  keep_teams_component <- names(comp_results$membership)[comp_results$membership == giant_idx]
  results_pruned <- results_pruned[
    results_pruned$home_team %in% keep_teams_component &
      results_pruned$away_team %in% keep_teams_component,
    , drop = FALSE
  ]
  team_summary <- team_summary[team_summary$team %in% keep_teams_component, , drop = FALSE]
  cat("Filtered to largest connected component with", length(keep_teams_component), "teams.\n")
}

teams <- sort(unique(c(results_pruned$home_team, results_pruned$away_team)))
p <- length(teams)
cat("\nNumber of teams retained:", p, "\n")

categories <- sort(unique(results_pruned$tournament))
cat("\nTournament categories retained:\n")
print(categories)

team_index <- setNames(seq_along(teams), teams)
home_idx <- unname(team_index[results_pruned$home_team])
away_idx <- unname(team_index[results_pruned$away_team])
n_matches <- nrow(results_pruned)

to_design_matrix <- function(df_subset) {
  n_obs <- nrow(df_subset)
  X <- matrix(0, nrow = n_obs, ncol = p)
  if (n_obs > 0) {
    home_id <- unname(team_index[df_subset$home_team])
    away_id <- unname(team_index[df_subset$away_team])
    X[cbind(seq_len(n_obs), home_id)] <- 1
    X[cbind(seq_len(n_obs), away_id)] <- -1
  }
  colnames(X) <- paste0("X", seq_len(p))
  Y <- integer(n_obs)
  if (n_obs > 0) {
    Y[df_subset$home_score > df_subset$away_score] <- 1L
    Y[df_subset$home_score < df_subset$away_score] <- -1L
  }
  cat_factor <- factor(df_subset$tournament, levels = categories)
  cat_ids <- as.integer(cat_factor)
  cat_matrix <- matrix(0, nrow = n_obs, ncol = length(categories))
  if (n_obs > 0) {
    cat_matrix[cbind(seq_len(n_obs), cat_ids)] <- 1
  }
  colnames(cat_matrix) <- paste0("cat_", seq_len(length(categories)))
  data.frame(X, Y = Y, cat_matrix, check.names = FALSE)
}

compute_pair_probabilities <- function(df) {
  pair_counts <- list()
  tournament_totals <- table(df$tournament)
  for (idx in seq_len(nrow(df))) {
    tour <- df$tournament[idx]
    home <- df$home_team[idx]
    away <- df$away_team[idx]
    if (identical(home, away)) {
      next
    }
    score_h <- df$home_score[idx]
    score_a <- df$away_score[idx]
    pair_names <- sort(c(home, away))
    pair_key <- paste(pair_names, collapse = "||")
    bucket <- pair_counts[[tour]]
    if (is.null(bucket)) {
      bucket <- list()
    }
    if (is.null(bucket[[pair_key]])) {
      bucket[[pair_key]] <- c(team1_wins = 0L, team2_wins = 0L, ties = 0L, total = 0L)
    }
    counts_vec <- bucket[[pair_key]]
    counts_vec["total"] <- counts_vec["total"] + 1L
    if (score_h > score_a) {
      # home win
      if (home == pair_names[1]) {
        counts_vec["team1_wins"] <- counts_vec["team1_wins"] + 1L
      } else {
        counts_vec["team2_wins"] <- counts_vec["team2_wins"] + 1L
      }
    } else if (score_h < score_a) {
      if (away == pair_names[1]) {
        counts_vec["team1_wins"] <- counts_vec["team1_wins"] + 1L
      } else {
        counts_vec["team2_wins"] <- counts_vec["team2_wins"] + 1L
      }
    } else {
      counts_vec["ties"] <- counts_vec["ties"] + 1L
    }
    bucket[[pair_key]] <- counts_vec
    pair_counts[[tour]] <- bucket
  }
  pair_summary_list <- list()
  global_pairs <- list()
  for (tour in names(pair_counts)) {
    bucket <- pair_counts[[tour]]
    pair_keys <- names(bucket)
    tour_df <- data.frame(
      tournament = character(0),
      team1 = character(0),
      team2 = character(0),
      total = integer(0),
      win_team1 = integer(0),
      win_team2 = integer(0),
      ties = integer(0),
      p_team1 = numeric(0),
      p_team2 = numeric(0),
      p_tie = numeric(0),
      stringsAsFactors = FALSE
    )
    for (key in pair_keys) {
      counts_vec <- bucket[[key]]
      teams_split <- strsplit(key, "\\|\\|")[[1]]
      total_games <- counts_vec["total"]
      if (total_games <= 0) next
      win1 <- counts_vec["team1_wins"]
      win2 <- counts_vec["team2_wins"]
      tie <- counts_vec["ties"]
      probs <- c(win1, win2, tie) / total_games
      tour_df <- rbind(
        tour_df,
        data.frame(
          tournament = tour,
          team1 = teams_split[1],
          team2 = teams_split[2],
          total = total_games,
          win_team1 = win1,
          win_team2 = win2,
          ties = tie,
          p_team1 = probs[1],
          p_team2 = probs[2],
          p_tie = probs[3],
          stringsAsFactors = FALSE
        )
      )
      global_pairs[[key]] <- TRUE
    }
    pair_summary_list[[tour]] <- tour_df
  }
  list(
    tournament_pairs = pair_summary_list,
    all_pair_keys = sort(names(global_pairs)),
    tournament_totals = tournament_totals
  )
}

pair_stats <- compute_pair_probabilities(results_pruned)
pair_summary_list <- pair_stats$tournament_pairs
all_pair_keys <- pair_stats$all_pair_keys
tournament_totals <- pair_stats$tournament_totals
pair_summary_all <- do.call(
  rbind,
  lapply(names(pair_summary_list), function(tnm) {
    df <- pair_summary_list[[tnm]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
  })
)
if (is.null(pair_summary_all)) {
  pair_summary_all <- data.frame(
    tournament = character(0),
    team1 = character(0),
    team2 = character(0),
    total = integer(0),
    win_team1 = integer(0),
    win_team2 = integer(0),
    ties = integer(0),
    p_team1 = numeric(0),
    p_team2 = numeric(0),
    p_tie = numeric(0),
    stringsAsFactors = FALSE
  )
}

cat("\nExample pair summary (first tournament):\n")
if (length(pair_summary_list) > 0) {
  show_tour <- names(pair_summary_list)[1]
  print(head(pair_summary_list[[show_tour]], 10))
}

n_pairs <- length(all_pair_keys)
connected_pair_matrix <- matrix(0L, nrow = n_pairs, ncol = 2)
pair_index_map <- setNames(seq_len(n_pairs), all_pair_keys)
for (key in all_pair_keys) {
  idx <- pair_index_map[[key]]
  teams_split <- strsplit(key, "\\|\\|")[[1]]
  connected_pair_matrix[idx, 1] <- team_index[[teams_split[1]]]
  connected_pair_matrix[idx, 2] <- team_index[[teams_split[2]]]
}
colnames(connected_pair_matrix) <- c("team1", "team2")

k_B <- length(categories)
prob_win1 <- matrix(0, nrow = n_pairs, ncol = k_B)
prob_win2 <- matrix(0, nrow = n_pairs, ncol = k_B)
prob_tie <- matrix(0, nrow = n_pairs, ncol = k_B)
pairs_count_matrix <- matrix(0, nrow = n_pairs, ncol = k_B)
pairs_by_tournament <- vector("list", k_B)
names(pairs_by_tournament) <- categories

for (b_idx in seq_along(categories)) {
  tour <- categories[b_idx]
  tour_df <- pair_summary_list[[tour]]
  if (is.null(tour_df) || nrow(tour_df) == 0) {
    pairs_by_tournament[[b_idx]] <- integer(0)
    next
  }
  pair_ids <- as.integer(pair_index_map[paste(tour_df$team1, tour_df$team2, sep = "||")])
  pairs_by_tournament[[b_idx]] <- pair_ids
  prob_win1[pair_ids, b_idx] <- tour_df$p_team1
  prob_win2[pair_ids, b_idx] <- tour_df$p_team2
  prob_tie[pair_ids, b_idx] <- tour_df$p_tie
  pairs_count_matrix[pair_ids, b_idx] <- tour_df$total
  pair_summary_list[[tour]]$pair_index <- pair_ids
}

tournament_count_vec <- as.numeric(tournament_totals[categories])
tournament_count_vec[is.na(tournament_count_vec)] <- 0
if (sum(tournament_count_vec) == 0) {
  stop("No tournament counts available after pruning.")
}
tournament_probabilities <- tournament_count_vec / sum(tournament_count_vec)
pair_total_counts <- rowSums(pairs_count_matrix)
pair_tournaments <- vector("list", n_pairs)
for (b_idx in seq_along(categories)) {
  pair_ids <- pairs_by_tournament[[b_idx]]
  if (!length(pair_ids)) next
  for (pid in pair_ids) {
    pair_tournaments[[pid]] <- c(pair_tournaments[[pid]], b_idx)
  }
}

pair_lookup_matrix <- matrix(0L, nrow = p, ncol = p)
for (pid in seq_len(n_pairs)) {
  i <- connected_pair_matrix[pid, 1]
  j <- connected_pair_matrix[pid, 2]
  pair_lookup_matrix[i, j] <- pid
  pair_lookup_matrix[j, i] <- pid
}

choose_tournament_for_pair <- function(pair_idx) {
  tours <- pair_tournaments[[pair_idx]]
  if (is.null(tours) || length(tours) == 0) {
    return(NA_integer_)
  }
  counts <- pairs_count_matrix[pair_idx, tours]
  counts[is.na(counts)] <- 0
  tours[which.max(counts)]
}

sample_outcome_empirical <- function(pair_idx, b_idx) {
  probs <- c(
    prob_win1[pair_idx, b_idx],
    prob_tie[pair_idx, b_idx],
    prob_win2[pair_idx, b_idx]
  )
  if (sum(probs) <= 0) {
    probs <- c(0.45, 0.1, 0.45)
  }
  probs <- probs / sum(probs)
  sample(c(1L, 0L, -1L), size = 1, prob = probs)
}

make_design_row <- function(pair_idx, b_idx, outcome) {
  row_vec <- numeric(p + 1 + k_B)
  teams_pair <- connected_pair_matrix[pair_idx, ]
  row_vec[teams_pair[1]] <- 1
  row_vec[teams_pair[2]] <- -1
  row_vec[p + 1] <- outcome
  row_vec[p + 1 + b_idx] <- 1
  row_vec
}

create_initial_samples <- function() {
  edge_matrix <- connected_pair_matrix
  edge_matrix_char <- matrix(as.character(edge_matrix), ncol = 2)
  g <- igraph::graph_from_edgelist(edge_matrix_char, directed = FALSE)
  igraph::E(g)$pair_idx <- seq_len(n_pairs)
  igraph::E(g)$weight <- 1 / (pair_total_counts + 1e-6)
  if (!igraph::is_connected(g)) {
    warning("Pair graph not fully connected; using spanning forest.")
    comps <- igraph::components(g)
    mst_pairs <- integer(0)
    for (comp_id in seq_len(comps$no)) {
      verts <- which(comps$membership == comp_id)
      sub_g <- igraph::induced_subgraph(g, verts)
      if (igraph::ecount(sub_g) == 0) next
      mst_sub <- igraph::mst(sub_g, weights = igraph::E(sub_g)$weight)
      if (igraph::ecount(mst_sub) > 0) {
        mst_pairs <- c(mst_pairs, igraph::E(mst_sub)$pair_idx)
      }
    }
    mst_pairs <- unique(mst_pairs)
  } else {
    mst_sub <- igraph::mst(g, weights = igraph::E(g)$weight)
    mst_pairs <- igraph::E(mst_sub)$pair_idx
  }
  records <- data.frame(pair_idx = integer(0), b_idx = integer(0), result = integer(0))
  row_list <- list()
  tournaments_used <- integer(0)
  for (pid in mst_pairs) {
    b_idx <- choose_tournament_for_pair(pid)
    if (is.na(b_idx)) next
    outcome <- sample_outcome_empirical(pid, b_idx)
    row_list[[length(row_list) + 1]] <- make_design_row(pid, b_idx, outcome)
    records <- rbind(records, data.frame(pair_idx = pid, b_idx = b_idx, result = outcome))
    tournaments_used <- c(tournaments_used, b_idx)
  }
  missing_tournaments <- setdiff(seq_along(categories), unique(tournaments_used))
  for (b_idx in missing_tournaments) {
    available_pairs <- pairs_by_tournament[[b_idx]]
    if (!length(available_pairs)) next
    counts <- pairs_count_matrix[available_pairs, b_idx]
    counts[is.na(counts)] <- 0
    pid <- available_pairs[which.max(counts)]
    outcome <- sample_outcome_empirical(pid, b_idx)
    row_list[[length(row_list) + 1]] <- make_design_row(pid, b_idx, outcome)
    records <- rbind(records, data.frame(pair_idx = pid, b_idx = b_idx, result = outcome))
  }
  if (length(row_list) == 0) {
    stop("Failed to create initial samples; no valid pairs.")
  }
  design_matrix <- do.call(rbind, row_list)
  colnames(design_matrix) <- c(paste0("X", seq_len(p)), "Y", paste0("cat_", seq_len(k_B)))
  list(
    df = as.data.frame(design_matrix, check.names = FALSE),
    records = records,
    mst_pairs = mst_pairs
  )
}

simulate_strategy <- function(
    strategy,
    total_samples,
    initial_df,
    initial_records,
    par_init,
    tournament_probabilities,
    par_global,
    theta_global,
    phi_global,
    soft_eta = 0.1,
    soft_delta = 0.1,
    R_theta = 5,
    phi_bounds = c(-3, 3),
    verbose = FALSE,
    with_replacement = TRUE) {
  n_init <- nrow(initial_df)
  if (total_samples < n_init) {
    stop("total_samples must be >= size of initial sample.")
  }
  df_current <- initial_df
  par_current <- par_init
  pi_counts <- matrix(0, nrow = n_pairs, ncol = k_B)
  if (nrow(initial_records) > 0) {
    for (row_idx in seq_len(nrow(initial_records))) {
      rec <- initial_records[row_idx, ]
      pi_counts[rec$pair_idx, rec$b_idx] <- pi_counts[rec$pair_idx, rec$b_idx] + 1
    }
  }
  available_pairs_list <- NULL
  if (!with_replacement) {
    available_pairs_list <- lapply(pairs_by_tournament, function(x) x)
  }
  get_allowed_pairs <- function(b_idx) {
    if (with_replacement) {
      pairs_by_tournament[[b_idx]]
    } else {
      available_pairs_list[[b_idx]]
    }
  }
  metrics_length <- total_samples - n_init + 1
  mse_vec <- numeric(metrics_length)
  hsr_vec <- numeric(metrics_length)
  tau_vec <- numeric(metrics_length)
  compute_metrics <- function(par_est) {
    theta_est <- par_est[1:p]
    phi_est <- par_est[(p + 1):(p + k_B)]
    mse_val <- mean((theta_est - theta_global)^2)
    hsr_val <- hetero_soft_ranking_loss(theta_est, theta_global, phi_est, phi_global, delta = soft_delta, eta = soft_eta)
    tau_val <- suppressWarnings(cor(theta_est, theta_global, method = "kendall"))
    c(mse_val, hsr_val, tau_val)
  }
  initial_metrics <- compute_metrics(par_current)
  mse_vec[1] <- initial_metrics[1]
  hsr_vec[1] <- initial_metrics[2]
  tau_vec[1] <- initial_metrics[3]
  current_total <- n_init
  for (step_idx in 2:metrics_length) {
    # select tournament
    b_star <- NA_integer_
    attempt <- 0
    while (is.na(b_star)) {
      attempt <- attempt + 1
      if (attempt > 200) {
        stop("Unable to sample a tournament with available pairs.")
      }
      b_candidate <- sample.int(k_B, size = 1, prob = tournament_probabilities)
      allowed_pairs_candidate <- get_allowed_pairs(b_candidate)
      if (!is.null(allowed_pairs_candidate) && length(allowed_pairs_candidate) > 0) {
        b_star <- b_candidate
      }
    }
    allowed_pairs <- get_allowed_pairs(b_star)
    pi_n <- if (current_total > 0) pi_counts / current_total else matrix(0, n_pairs, k_B)
    if (strategy == "soft") {
      H_dynamic <- build_soft_ranking_H(par_current, p, eta = soft_eta, delta = soft_delta)
      a_star <- select_action_restricted(
        par_current, H_dynamic, pi_n, b_star, allowed_pairs,
        connected_pair_matrix, p, k_B
      )
    } else if (strategy == "l2") {
      H_static <- build_H(p, k_B, eta_theta = 1, eta = soft_eta)
      a_star <- select_action_restricted(
        par_current, H_static, pi_n, b_star, allowed_pairs,
        connected_pair_matrix, p, k_B
      )
    } else if (strategy == "uniform") {
      a_star <- sample(allowed_pairs, size = 1)
    } else if (strategy == "uncertainty") {
      best_val <- -Inf
      a_star <- allowed_pairs[1]
      for (candidate in allowed_pairs) {
        ent_val <- three_outcome_entropy(candidate, b_star, par_current, connected_pair_matrix, p, k_B)
        if (ent_val > best_val) {
          best_val <- ent_val
          a_star <- candidate
        }
      }
    } else {
      stop("Unknown strategy")
    }
    if (is.na(a_star) || length(a_star) == 0) {
      next
    }
    if (!with_replacement) {
      available_pairs_list[[b_star]] <- setdiff(available_pairs_list[[b_star]], a_star)
    }
    outcome <- sample_outcome_empirical(a_star, b_star)
    new_row <- make_design_row(a_star, b_star, outcome)
    new_df <- as.data.frame(matrix(new_row, nrow = 1), stringsAsFactors = FALSE)
    colnames(new_df) <- c(paste0("X", seq_len(p)), "Y", paste0("cat_", seq_len(k_B)))
    df_current <- rbind(df_current, new_df)
    rownames(df_current) <- NULL
    pi_counts[a_star, b_star] <- pi_counts[a_star, b_star] + 1
    current_total <- current_total + 1
    fit_res <- fit_btd_mle(
      df_data = df_current,
      p = p,
      k_B = k_B,
      R = R_theta,
      phi_lower = phi_bounds[1],
      phi_upper = phi_bounds[2],
      init_par = par_current,
      maxeval = 500,
      verbose = FALSE
    )
    par_current <- fit_res$solution
    metrics_vals <- compute_metrics(par_current)
    mse_vec[step_idx] <- metrics_vals[1]
    hsr_vec[step_idx] <- metrics_vals[2]
    tau_vec[step_idx] <- metrics_vals[3]
  }
  list(
    mse = mse_vec,
    hsr = hsr_vec,
    tau = tau_vec,
    final_par = par_current
  )
}

k_B <- length(categories)
df_btd <- to_design_matrix(results_pruned)
Y_vec <- df_btd$Y

# Dataset with tournament/date/result
match_outcomes <- data.frame(
  date = results_pruned$date,
  tournament = results_pruned$tournament,
  result = Y_vec,
  stringsAsFactors = FALSE
)

cat("\nSample of retained match outcomes (tournament/date/result):\n")
print(head(match_outcomes, 20))

set.seed(MASTER_SEED + 1L)
initial_sample <- create_initial_samples()
df_initial <- initial_sample$df
initial_records <- initial_sample$records
mst_pair_indices <- initial_sample$mst_pairs

mst_edge_matrix <- connected_pair_matrix[mst_pair_indices, , drop = FALSE]
mst_graph <- igraph::graph_from_edgelist(matrix(as.character(mst_edge_matrix), ncol = 2), directed = FALSE)
mst_connected <- if (nrow(mst_edge_matrix) > 0) igraph::is_connected(mst_graph) else FALSE
cat("MST edges count:", nrow(mst_edge_matrix), "Connected:", mst_connected, "\n")

initial_edge_matrix <- connected_pair_matrix[initial_records$pair_idx, , drop = FALSE]
initial_graph <- igraph::graph_from_edgelist(matrix(as.character(initial_edge_matrix), ncol = 2), directed = FALSE)
initial_connected <- if (nrow(initial_edge_matrix) > 0) igraph::is_connected(initial_graph) else FALSE
cat("Initial design edges count:", nrow(initial_edge_matrix), "Connected:", initial_connected, "\n")

R_theta_constraint <- 5
initial_fit <- fit_btd_mle(
  df_data = df_initial,
  p = p,
  k_B = k_B,
  R = R_theta_constraint,
  phi_lower = -3,
  phi_upper = 3,
  init_par = NULL,
  maxeval = 1000,
  verbose = FALSE
)
par_initial <- initial_fit$solution

years_pruned <- as.integer(format(results_pruned$date, "%Y"))
global_cutoff_year <- max(years_pruned) - 49
results_global <- results_pruned[years_pruned >= global_cutoff_year, , drop = FALSE]
if (nrow(results_global) < p) {
  results_global <- results_pruned
}
df_global <- to_design_matrix(results_global)
global_fit <- fit_btd_mle(
  df_data = df_global,
  p = p,
  k_B = k_B,
  R = R_theta_constraint,
  phi_lower = -3,
  phi_upper = 3,
  init_par = par_initial,
  maxeval = 1500,
  verbose = FALSE
)
par_global <- global_fit$solution
theta_global <- par_global[1:p]
phi_global <- par_global[(p + 1):(p + k_B)]

total_samples_target <- nrow(df_initial) + TOTAL_ADDITIONAL_SAMPLES
strategies <- c("soft", "l2", "uniform", "uncertainty")
mc_reps <- MC_REPS
metrics_len <- total_samples_target - nrow(df_initial) + 1
strategy_results <- lapply(strategies, function(x) {
  list(
    mse = matrix(0, nrow = mc_reps, ncol = metrics_len),
    hsr = matrix(0, nrow = mc_reps, ncol = metrics_len),
    tau = matrix(0, nrow = mc_reps, ncol = metrics_len)
  )
})
names(strategy_results) <- strategies

for (rep_idx in seq_len(mc_reps)) {
  set.seed(MASTER_SEED + 1000L + rep_idx)
  for (s_idx in seq_along(strategies)) {
    strategy_name <- strategies[s_idx]
    sim_res <- simulate_strategy(
      strategy = strategy_name,
      total_samples = total_samples_target,
      initial_df = df_initial,
      initial_records = initial_records,
      par_init = par_initial,
      tournament_probabilities = tournament_probabilities,
      par_global = par_global,
      theta_global = theta_global,
      phi_global = phi_global,
      soft_eta = 0.1,
      soft_delta = 0.1,
      R_theta = R_theta_constraint,
      phi_bounds = c(-3, 3),
      verbose = FALSE,
      with_replacement = WITH_REPLACEMENT
    )
    strategy_results[[strategy_name]]$mse[rep_idx, ] <- sim_res$mse
    strategy_results[[strategy_name]]$hsr[rep_idx, ] <- sim_res$hsr
    strategy_results[[strategy_name]]$tau[rep_idx, ] <- sim_res$tau
  }
}

sample_sizes <- seq(nrow(df_initial), total_samples_target)
strategy_means <- lapply(strategy_results, function(res) {
  list(
    mse = colMeans(res$mse),
    hsr = colMeans(res$hsr),
    tau = colMeans(res$tau)
  )
})

plot_colors <- c(soft = "steelblue", l2 = "sienna", uniform = "darkgreen", uncertainty = "purple")

plot_metric_curve <- function(metric_name, y_label, file_name) {
  png(file.path(project_root, "real data", file_name), width = 1200, height = 800, res = 150)
  on.exit(dev.off(), add = TRUE)
  y_values <- sapply(strategies, function(name) strategy_means[[name]][[metric_name]])
  start_idx <- which.min(abs(sample_sizes - (nrow(df_initial) + N_PLOT_INITIAL - 1L)))
  if (start_idx < 1) start_idx <- 1
  y_slice <- y_values[start_idx:nrow(y_values), , drop = FALSE]
  x_slice <- sample_sizes[start_idx:length(sample_sizes)]
  y_min <- min(y_slice, na.rm = TRUE)
  y_max <- max(y_slice, na.rm = TRUE)
  plot(x_slice, y_slice[, 1], type = "l", lwd = 2, col = plot_colors[strategies[1]],
       xlab = "Sample size", ylab = y_label, ylim = c(y_min, y_max), main = paste("Average", y_label))
  if (length(strategies) > 1) {
    for (idx in 2:length(strategies)) {
      lines(x_slice, y_slice[, idx], lwd = 2, col = plot_colors[strategies[idx]])
    }
  }
  legend("topright", legend = strategies, col = plot_colors[strategies], lwd = 2, bty = "n")
}

plot_metric_curve("mse", "Mean Squared Error", "mc_average_mse.png")
plot_metric_curve("hsr", "Soft Ranking Loss", "mc_average_soft_ranking_loss.png")
plot_metric_curve("tau", "Kendall's Tau", "mc_average_kendall_tau.png")

metric_summary_df <- do.call(rbind, lapply(strategies, function(name) {
  data.frame(
    strategy = name,
    sample_size = sample_sizes,
    mse = strategy_means[[name]]$mse,
    soft_ranking_loss = strategy_means[[name]]$hsr,
    kendall_tau = strategy_means[[name]]$tau,
    stringsAsFactors = FALSE
  )
}))
write.csv(metric_summary_df, file.path(project_root, "real data", "mc_average_metrics.csv"), row.names = FALSE)

# Per-tournament team summaries
summarize_teams <- function(df) {
  teams_all <- sort(unique(c(df$home_team, df$away_team)))
  matches <- integer(length(teams_all))
  wins <- integer(length(teams_all))
  losses <- integer(length(teams_all))
  draws <- integer(length(teams_all))
  names(matches) <- names(wins) <- names(losses) <- names(draws) <- teams_all
  for (i in seq_len(nrow(df))) {
    h <- df$home_team[i]
    a <- df$away_team[i]
    hs <- df$home_score[i]
    as <- df$away_score[i]
    matches[h] <- matches[h] + 1L
    matches[a] <- matches[a] + 1L
    if (hs > as) {
      wins[h] <- wins[h] + 1L
      losses[a] <- losses[a] + 1L
    } else if (hs < as) {
      wins[a] <- wins[a] + 1L
      losses[h] <- losses[h] + 1L
    } else {
      draws[h] <- draws[h] + 1L
      draws[a] <- draws[a] + 1L
    }
  }
  data.frame(
    team = teams_all,
    matches = matches,
    wins = wins,
    losses = losses,
    draws = draws,
    stringsAsFactors = FALSE
  )
}

tournament_team_summaries <- lapply(split(results_pruned, results_pruned$tournament), summarize_teams)

cat("\nPer-tournament team summaries (example: first tournament listed):\n")
first_tournament <- categories[1]
print(head(tournament_team_summaries[[first_tournament]], 10))

tournament_team_summary_df <- do.call(rbind, lapply(names(tournament_team_summaries), function(tnm) {
  df <- tournament_team_summaries[[tnm]]
  df$tournament <- tnm
  df
}))
tournament_team_summary_df <- tournament_team_summary_df[, c("tournament", "team", "matches", "wins", "losses", "draws")]

# Overall summary across tournaments
overall_summary <- summarize_teams(results_pruned)
overall_summary <- overall_summary[order(-overall_summary$matches), , drop = FALSE]

cat("\nOverall team summary across selected tournaments (top 20):\n")
print(head(overall_summary, 20))

write.csv(match_outcomes, file.path(project_root, "real data", "selected_match_outcomes.csv"), row.names = FALSE)
write.csv(tournament_team_summary_df, file.path(project_root, "real data", "tournament_team_summaries.csv"), row.names = FALSE)
write.csv(overall_summary, file.path(project_root, "real data", "overall_team_summary.csv"), row.names = FALSE)
write.csv(pair_summary_all, file.path(project_root, "real data", "pair_empirical_summary.csv"), row.names = FALSE)
write.csv(initial_records, file.path(project_root, "real data", "initial_sample_records.csv"), row.names = FALSE)
write.csv(df_initial, file.path(project_root, "real data", "initial_sample_design.csv"), row.names = FALSE)

cat("\nFitting BTD model with heterogeneous tie parameters...\n")
R_theta <- 8
fit_res <- fit_btd_mle(
  df_data = df_btd,
  p = p,
  k_B = k_B,
  R = R_theta,
  phi_lower = -3,
  phi_upper = 3,
  init_par = par_initial,
  method = "NLOPT_LD_SLSQP",
  maxeval = 2000,
  verbose = TRUE
)

theta_hat <- fit_res$solution[1:p]
phi_hat <- fit_res$solution[(p + 1):(p + k_B)]

latent_scores <- data.frame(
  team = teams,
  theta = theta_hat,
  stringsAsFactors = FALSE
)
latent_scores <- latent_scores[order(-latent_scores$theta), , drop = FALSE]

tie_parameters <- data.frame(
  tournament = categories,
  phi = phi_hat,
  tie_weight = exp(phi_hat),
  stringsAsFactors = FALSE
)

cat("\nTop 20 latent team scores:\n")
print(head(latent_scores, 20))

cat("\nTie parameters by tournament:\n")
print(tie_parameters)

write.csv(latent_scores, file.path(project_root, "real data", "latent_team_scores.csv"), row.names = FALSE)
write.csv(tie_parameters, file.path(project_root, "real data", "tie_parameters.csv"), row.names = FALSE)

cat("\nResults written to:\n")
cat("  - real data/selected_match_outcomes.csv\n")
cat("  - real data/tournament_team_summaries.csv\n")
cat("  - real data/overall_team_summary.csv\n")
cat("  - real data/pair_empirical_summary.csv\n")
cat("  - real data/initial_sample_records.csv\n")
cat("  - real data/initial_sample_design.csv\n")
cat("  - real data/latent_team_scores.csv\n")
cat("  - real data/tie_parameters.csv\n")
cat("  - real data/mc_average_mse.png\n")
cat("  - real data/mc_average_soft_ranking_loss.png\n")
cat("  - real data/mc_average_kendall_tau.png\n")
cat("  - real data/mc_average_metrics.csv\n")
cat("Adaptive sampling was performed ", if (WITH_REPLACEMENT) "with" else "without", " replacement.\n", sep = "")

elapsed_time <- proc.time() - script_start_time
cat("Total runtime for MC_REPS =", mc_reps, "was", elapsed_time["elapsed"], "seconds.\n")
if (mc_reps > 0) {
  per_rep <- elapsed_time["elapsed"] / mc_reps
  cat("Approximate time per replicate:", per_rep, "seconds.\n")
  cat("Estimated runtime for 1000 replicates:", per_rep * 1000, "seconds.\n")
}
