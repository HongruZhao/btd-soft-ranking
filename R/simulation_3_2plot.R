setwd("/panfs/jay/groups/27/shenx/zhao1118/soft_ranking_loss")

############################################################
## Load results and build six histograms (tau and capped D)
############################################################

dat01  <- read.csv("tau_D_c0.1.csv")     # c = 0.10
dat05  <- read.csv("tau_D_c0.05.csv")    # c = 0.05
dat025 <- read.csv("tau_D_c0.025.csv")   # c = 0.025

tau01  <- dat01$tau
tau05  <- dat05$tau
tau025 <- dat025$tau

## 99th chi^2 quantile for D
df      <- length(theta) + length(sigma2_vec)   # here df = 5
D_cap99 <- qchisq(0.99, df = df)

D01  <- pmin(dat01$D,  D_cap99)
D05  <- pmin(dat05$D,  D_cap99)
D025 <- pmin(dat025$D, D_cap99)

## Common breaks
breaks_tau <- seq(min(c(tau01, tau05, tau025)),
                  max(c(tau01, tau05, tau025)),
                  length.out = 20)

breaks_D <- seq(0, D_cap99, length.out = 15)

## Deterministic tau_det for each c
tau_det01  <- get_tau_det(0.10)
tau_det05  <- get_tau_det(0.05)
tau_det025 <- get_tau_det(0.025)

## Helper: histogram plus chi^2 overlay on a unified y scale
plot_hist_with_chisq <- function(z, breaks, df,
                                 xlab, col = "skyblue", border = "white") {
  
  ## 1. bin counts (no 'freq' argument while plot = FALSE)
  h0       <- hist(z, breaks = breaks, plot = FALSE)
  counts   <- h0$counts
  binwidth <- diff(h0$breaks)
  dens     <- counts / (sum(counts) * binwidth)
  
  ## 2. chi^2 density on a dense grid
  xg       <- seq(min(breaks), max(breaks), length.out = 2000)
  y_chi    <- dchisq(xg, df = df)
  y_max    <- max(dens, y_chi)
  
  ## 3. final plot
  hist(z, breaks = breaks, freq = FALSE, ylim = c(0, y_max),
       main = "", xlab = xlab, col = col, border = border)
  curve(dchisq(x, df = df), from = min(breaks), to = max(breaks),
        add = TRUE, col = "black", lwd = 2)
}

## ------------------------------------------------------------
## Write the 3 x 2 histogram panel to a PDF file
## ------------------------------------------------------------
pdf(
  file      = "tau_D_histograms.pdf",
  width     = 6,                  # inches
  height    = 6 * 1.176,          # 7.06 inches, keeps the original ratio
  pointsize = 12
)

par(mfrow = c(3, 2), mar = c(4, 4, 1, 1))

## tau (c = 0.10)
hist(tau01, breaks = breaks_tau, freq = FALSE,
     main = "", xlab = expression(tau[c] ~ "(c = 0.10)"),
     col  = "skyblue", border = "white")
abline(v = tau_det01, col = "black", lwd = 2)
legend("topright", legend = expression(tau[det]),
       bty = "n", col = "black", lwd = 2)

## D (c = 0.10)
plot_hist_with_chisq(D01, breaks_D, df,
                     xlab = expression(D ~ "(c = 0.10)"))

## tau (c = 0.05)
hist(tau05, breaks = breaks_tau, freq = FALSE,
     main = "", xlab = expression(tau[c] ~ "(c = 0.05)"),
     col  = "skyblue", border = "white")
abline(v = tau_det05, col = "black", lwd = 2)
legend("topright", legend = expression(tau[det]),
       bty = "n", col = "black", lwd = 2)

## D (c = 0.05)
plot_hist_with_chisq(D05, breaks_D, df,
                     xlab = expression(D ~ "(c = 0.05)"))

## tau (c = 0.025)
hist(tau025, breaks = breaks_tau, freq = FALSE,
     main = "", xlab = expression(tau[c] ~ "(c = 0.025)"),
     col  = "skyblue", border = "white")
abline(v = tau_det025, col = "black", lwd = 2)
legend("topright", legend = expression(tau[det]),
       bty = "n", col = "black", lwd = 2)

## D (c = 0.025)
plot_hist_with_chisq(D025, breaks_D, df,
                     xlab = expression(D ~ "(c = 0.025)"))

par(mfrow = c(1, 1))   # restore defaults
dev.off()               # save the PDF

