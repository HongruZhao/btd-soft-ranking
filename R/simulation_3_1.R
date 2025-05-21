set.seed(123)       # for reproducibility

# Function to simulate one stopped‐sample mean
simulate_barX_tau <- function(theta, c, max_n = 1e8) {
  n     <- 0
  sum_x <- 0
  repeat {
    n     <- n + 1
    x     <- rnorm(1, mean = theta, sd = 1)
    sum_x <- sum_x + x
    # stopping rule: n * \bar{X}_n >= c
    if (n * (sum_x / n - theta) >= c || n >= max_n) break
  }
  
  barX_tau <- sum_x / n          # sample mean at stopping time
  list(
    mean             = barX_tau,               # \bar{X}_τ
    tau              = n,                      # stopping time τ
    scaled_deviation = sqrt(n) * (barX_tau - theta)  # √τ ( \bar{X}_τ − θ )
  )
}


# Parameters
theta  <- 1    # true θ; change as desired
c_val  <- 3   # threshold
nsim   <- 10000  # number of replications

# Run the simulation

sim_mat   <- replicate(nsim, simulate_barX_tau(theta, c_val))  # 2 × nsim matrix
tau_vals  <- as.numeric(sim_mat["tau",          ])     # fast column access
barX_vals <- as.numeric(sim_mat["scaled_deviation",   ])
# Quick diagnostics
mean(barX_vals)    # empirical mean of barX_{τ_c}
sd(barX_vals)      # empirical SD of barX_{τ_c}

## save current settings so you can restore them later
op <- par(no.readonly = TRUE)



par(mar = c(5, 5, 1, 1))           # optional: wider left margin
hist(
  barX_vals,
  breaks = 25,
  freq   = FALSE,
  main   = "",
  xlab   = expression(sqrt(tau[c]) * (bar(X)[tau[c]] - theta)),
  col    = "skyblue",
  border = "white"
)




par(mar = c(5, 5, 1, 1)) 
## --- QQ-plot against N(0,1) -----------------------------------------------
qqnorm(barX_vals,
       main = "",
       pch  = 19,           # solid points
       col  = "steelblue")
qqline(barX_vals, col = "red", lwd = 2)



par(mar = c(5, 5, 1, 1))            # ample left margin
hist(
  log(tau_vals),
  freq   = FALSE,                   # density scale
  breaks = 30,
  main   = "",
  xlab   = expression(log(tau[c])),
  col    = "skyblue",
  border = "white"
)











results_df <- data.frame(
  tau  = tau_vals,
  barX = barX_vals
)

# Write to CSV (no row numbers)
write.csv(results_df, file = "fail.csv", row.names = FALSE)








# -------------------------------------------------------------------
# Two-panel PDF:  (L)  sqrt(tau_c)*(X̄_τc−θ)   |   (R)  log τ_c
# -------------------------------------------------------------------
pdf(
  file      = "hist_failure_combined.pdf",  # one output file
  width     = 12,  # twice the single-plot width (6 in × 2 panels)
  height    = 4.5,
  pointsize = 12
)

par(mfrow = c(1, 2))            # 1 row, 2 columns
par(mar  = c(5, 5, 1, 1))       # uniform margins for both panels

## ---- (1) Histogram of scaled deviation ---------------------------------
hist(
  barX_vals,
  breaks = 25,
  freq   = FALSE,
  main   = "",
  xlab   = expression(sqrt(tau[c]) * (bar(X)[tau[c]] - theta)),
  col    = "skyblue",
  border = "white"
)

## ---- (2) Histogram of log(tau_c) ---------------------------------------
hist(
  log(tau_vals),
  breaks = 30,
  freq   = FALSE,
  main   = "",
  xlab   = expression(log(tau[c])),
  col    = "skyblue",
  border = "white"
)

dev.off()                       # close the PDF device

