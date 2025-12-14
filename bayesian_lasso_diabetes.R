rm(list = ls())
set.seed(100)

# ============================================================
# Bayesian Lasso (Park & Casella, 2008) — Diabetes (lars)
# Reproduz:
#   (i) LS, Lasso-CV e Ridge-CV em dados padronizados
#   (ii) Gibbs sampler do Bayesian Lasso com mistura Normal–Exponencial
#   (iii) Matching de lambda via ||beta||_1 (BL mediana ~ Lasso-CV)
#   (iv) Diagnósticos MCMC e figuras de comparação
# ============================================================

# ============================================================
# 0) Dependencies
# ============================================================
pkgs <- c("lars", "MASS", "statmod", "coda")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(lars)
library(MASS)
library(statmod)
library(coda)

# ============================================================
# 1) Data: Diabetes (Efron et al.) + standardization
# ============================================================
data(diabetes, package = "lars")
X_raw <- diabetes$x
y_raw <- diabetes$y

n <- nrow(X_raw)
p <- ncol(X_raw)
cat("n =", n, " | p =", p, "\n")

# Variable names used in the lars::diabetes dataset
var_names <- c("age","sex","bmi","map","tc","ldl","hdl","tch","ltg","glu")

# Standardize X (mean 0, sd 1) and center y
X <- scale(X_raw, center = TRUE, scale = TRUE)
y <- y_raw - mean(y_raw)
colnames(X) <- var_names

# Sanity checks (should be ~0 up to floating-point error)
cat("Check X: max|mean| =", max(abs(colMeans(X))), "\n")
cat("Check X: max|sd-1| =", max(abs(apply(X, 2, sd) - 1)), "\n")
cat("Check y: mean =", mean(y), "\n")

# Precompute cross-products
XtX <- crossprod(X)     # X'X
Xty <- crossprod(X, y)  # X'y

# ============================================================
# 2) Bayesian Lasso (Park & Casella): Gibbs sampler
# ============================================================
# Model:
#   y | beta, sig2 ~ N(X beta, sig2 I)
#
# Priors:
#   beta_j | sig2, tau2_j ~ N(0, sig2 * tau2_j)
#   tau2_j ~ Exp(rate = lambda^2 / 2)
#   sig2 ~ IG(a0, b0)  (implemented as 1/Gamma)
#
# Full conditionals:
#   beta | sig2, tau2, y ~ N_p(A^{-1} X'y, sig2 A^{-1}), A = X'X + D^{-1}
#   sig2 | beta, tau2, y ~ IG(a0+(n+p)/2, b0+0.5*(||y-Xb||^2 + sum(b_j^2/tau2_j)))
#   w_j = 1/tau2_j | beta_j, sig2 ~ InvGaussian(mean=mu_j, shape=lambda^2),
#     mu_j = sqrt((lambda^2 * sig2) / beta_j^2)
#
# Reference: Park & Casella (2008), "The Bayesian Lasso"
# ============================================================

a0 <- 1e-2
b0 <- 1e-2

# Numerical guard for beta_j ~ 0 in mu_j
eps_beta <- 1e-12

sample_beta <- function(XtX, Xty, sig2, tau2) {
  A <- XtX + diag(1 / tau2, nrow = length(tau2))
  R <- chol(A)  # A = R'R
  
  mean_beta <- backsolve(R, forwardsolve(t(R), Xty))
  z <- rnorm(length(tau2))
  v <- backsolve(R, z)  # v ~ N(0, A^{-1})
  
  as.vector(mean_beta + sqrt(sig2) * v)
}

sample_sig2 <- function(y, X, beta, tau2, a0, b0) {
  n <- length(y)
  p <- length(beta)
  res <- as.vector(y - X %*% beta)
  
  quad <- sum(res^2) + sum(beta^2 / tau2)
  shape <- a0 + 0.5 * (n + p)
  rate  <- b0 + 0.5 * quad
  
  1 / rgamma(1, shape = shape, rate = rate)
}

sample_tau2 <- function(beta, sig2, lambda2, eps_beta = 1e-12) {
  p <- length(beta)
  tau2_new <- numeric(p)
  
  for (j in 1:p) {
    bj2 <- max(beta[j]^2, eps_beta)
    mu_j <- sqrt((lambda2 * sig2) / bj2)
    w <- statmod::rinvgauss(1, mean = mu_j, shape = lambda2)  # w = 1/tau2
    tau2_new[j] <- 1 / w
  }
  tau2_new
}

run_bl_mcmc <- function(lambda_val, X, y, XtX, Xty,
                        a0, b0,
                        M = 12000, burn = 2000, thin = 10,
                        init = NULL,
                        eps_beta = 1e-12,
                        verbose = TRUE) {
  
  p <- ncol(X)
  lambda2 <- lambda_val^2
  
  if (is.null(init)) {
    beta <- rep(0, p)
    sig2 <- var(y)
    tau2 <- rep(1, p)
  } else {
    beta <- init$beta
    sig2 <- init$sig2
    tau2 <- init$tau2
  }
  
  beta_chain <- matrix(NA_real_, nrow = M, ncol = p)
  sig2_chain <- numeric(M)
  
  beta_chain[1,] <- beta
  sig2_chain[1]  <- sig2
  
  for (k in 2:M) {
    beta <- sample_beta(XtX, Xty, sig2, tau2)
    sig2 <- sample_sig2(y, X, beta, tau2, a0, b0)
    tau2 <- sample_tau2(beta, sig2, lambda2, eps_beta)
    
    beta_chain[k,] <- beta
    sig2_chain[k]  <- sig2
    
    if (verbose && (k %% 2000 == 0)) cat("BL MCMC:", k, "/", M, "\n")
  }
  
  idx <- seq(burn + 1, M, by = thin)
  beta_post <- beta_chain[idx, , drop = FALSE]
  sig2_post <- sig2_chain[idx]
  
  list(
    lambda = lambda_val,
    beta_chain = beta_chain,
    sig2_chain = sig2_chain,
    beta_post = beta_post,
    sig2_post = sig2_post,
    last_state = list(beta = beta, sig2 = sig2, tau2 = tau2)
  )
}

summarize_bl <- function(beta_post, var_names) {
  beta_median <- apply(beta_post, 2, median)
  beta_ci <- t(apply(beta_post, 2, quantile, probs = c(0.025, 0.975)))
  data.frame(
    var = var_names,
    BL_median = beta_median,
    BL_low = beta_ci[,1],
    BL_high = beta_ci[,2],
    row.names = NULL
  )
}

# ============================================================
# 3) Baselines: LS, Lasso-CV, Ridge-CV
# ============================================================

# Least Squares on standardized X / centered y
beta_ls <- as.vector(solve(XtX, Xty))

# Lasso-CV with explicit folds (reproducible)
set.seed(100)
K <- 10
fold_id <- sample(rep(1:K, length.out = n))

lasso_cv_fraction <- function(X, y, fold_id, s_grid,
                              normalize = FALSE, intercept = FALSE) {
  
  K <- max(fold_id)
  mse_mat <- matrix(NA_real_, nrow = length(s_grid), ncol = K)
  
  for (kk in 1:K) {
    test <- which(fold_id == kk)
    train <- setdiff(seq_len(nrow(X)), test)
    
    Xtr <- X[train, , drop = FALSE]
    ytr <- y[train]
    Xte <- X[test, , drop = FALSE]
    yte <- y[test]
    
    fit <- lars(x = Xtr, y = ytr, type = "lasso",
                normalize = normalize, intercept = intercept)
    
    coef_mat <- predict(fit, s = s_grid, type = "coefficients", mode = "fraction")$coefficients
    pred_mat <- coef_mat %*% t(Xte)
    
    mse_mat[, kk] <- rowMeans((pred_mat - matrix(yte, nrow = length(s_grid),
                                                 ncol = length(yte), byrow = TRUE))^2)
  }
  
  mse_mean <- rowMeans(mse_mat)
  s_star <- s_grid[which.min(mse_mean)]
  
  list(s_star = s_star, mse_mean = mse_mean, mse_by_fold = mse_mat)
}

s_grid_cv <- seq(0, 1, length.out = 200)
cv_lasso_out <- lasso_cv_fraction(X, y, fold_id, s_grid_cv, normalize = FALSE, intercept = FALSE)
s_cv <- cv_lasso_out$s_star
cat("s_cv (Lasso CV) =", s_cv, "\n")

lasso_path <- lars(x = X, y = y, type = "lasso", normalize = FALSE, intercept = FALSE)
beta_lasso_cv <- as.vector(predict(lasso_path, s = s_cv, type = "coefficients", mode = "fraction")$coefficients)

L1_lasso_cv <- sum(abs(beta_lasso_cv))
cat("||beta||_1 (Lasso-CV) =", L1_lasso_cv, "\n")

# Ridge-CV (10-fold) on a grid of k
k_grid <- exp(seq(log(1e-3), log(1e3), length.out = 200))
I_p <- diag(p)
cv_mse_ridge <- numeric(length(k_grid))

for (ii in seq_along(k_grid)) {
  k <- k_grid[ii]
  mse_fold <- numeric(K)
  
  for (kk in 1:K) {
    test <- which(fold_id == kk)
    train <- setdiff(1:n, test)
    
    Xtr <- X[train, , drop = FALSE]
    ytr <- y[train]
    Xte <- X[test, , drop = FALSE]
    yte <- y[test]
    
    XtX_tr <- crossprod(Xtr)
    Xty_tr <- crossprod(Xtr, ytr)
    
    beta_k <- solve(XtX_tr + k * I_p, Xty_tr)
    pred <- as.vector(Xte %*% beta_k)
    
    mse_fold[kk] <- mean((yte - pred)^2)
  }
  cv_mse_ridge[ii] <- mean(mse_fold)
}

k_star <- k_grid[which.min(cv_mse_ridge)]
cat("k* (Ridge-CV) =", k_star, "\n")

ridge_coef_mat <- matrix(NA_real_, nrow = length(k_grid), ncol = p)
for (ii in seq_along(k_grid)) {
  k <- k_grid[ii]
  ridge_coef_mat[ii, ] <- solve(XtX + k * I_p, Xty)
}
beta_ridge_star <- as.vector(solve(XtX + k_star * I_p, Xty))

# ============================================================
# 4) Lambda grid (Bayesian Lasso) + L1 matching
# ============================================================

M_grid    <- 6000
burn_grid <- 1000
thin_grid <- 5

lambda_grid <- exp(seq(log(0.05), log(300), length.out = 45))

# Warm-start: last state at lambda_i initializes lambda_{i+1}
res_list <- vector("list", length(lambda_grid))
init_state <- NULL

t0 <- Sys.time()
for (i in seq_along(lambda_grid)) {
  cat("BL grid:", i, "/", length(lambda_grid), " | lambda =", lambda_grid[i], "\n")
  
  out_i <- run_bl_mcmc(lambda_grid[i], X, y, XtX, Xty, a0, b0,
                       M = M_grid, burn = burn_grid, thin = thin_grid,
                       init = init_state, eps_beta = eps_beta, verbose = FALSE)
  
  beta_med_i <- apply(out_i$beta_post, 2, median)
  
  res_list[[i]] <- list(
    lambda = lambda_grid[i],
    beta_median = beta_med_i,
    L1_median = sum(abs(beta_med_i))
  )
  
  init_state <- out_i$last_state
}
cat("BL grid total time:", Sys.time() - t0, "\n")

lambda_vec <- sapply(res_list, `[[`, "lambda")
L1_vec     <- sapply(res_list, `[[`, "L1_median")
beta_med_mat <- do.call(rbind, lapply(res_list, `[[`, "beta_median"))
colnames(beta_med_mat) <- var_names

# lambda* chosen by matching ||median(beta)||_1 to ||beta_lasso_cv||_1
i_star <- which.min(abs(L1_vec - L1_lasso_cv))
lambda_star <- lambda_vec[i_star]

cat("lambda* (L1 matching) =", lambda_star, "\n")
cat("||median(beta_BL)||_1 =", L1_vec[i_star], " | ||beta_LassoCV||_1 =", L1_lasso_cv, "\n")

# ============================================================
# 5) Final long chain at lambda* (posterior summaries + diagnostics)
# ============================================================

M_final    <- 12000
burn_final <- 2000
thin_final <- 10

fit_bl_star <- run_bl_mcmc(lambda_star, X, y, XtX, Xty, a0, b0,
                           M = M_final, burn = burn_final, thin = thin_final,
                           init = NULL, eps_beta = eps_beta, verbose = TRUE)

beta_chain <- fit_bl_star$beta_chain
sig2_chain <- fit_bl_star$sig2_chain
beta_post  <- fit_bl_star$beta_post
sig2_post  <- fit_bl_star$sig2_post

cat("Posterior draws (after burn/thin):", nrow(beta_post), "\n")

# ============================================================
# 6) MCMC diagnostics: trace, density, ACF, ESS
# ============================================================

beta_med_tmp <- apply(beta_post, 2, median)
ord_beta <- order(abs(beta_med_tmp), decreasing = TRUE)
j1 <- ord_beta[1]; j2 <- ord_beta[2]

par(mfrow = c(2,2))
plot(sig2_chain, type = "l", main = expression("Trace: " * sigma^2),
     xlab = "Iteration", ylab = expression(sigma^2))
plot(density(sig2_post, adjust = 2), main = expression("Posterior density: " * sigma^2),
     xlab = expression(sigma^2), lwd = 2)
plot(beta_chain[, j1], type = "l", main = paste0("Trace: ", var_names[j1]),
     xlab = "Iteration", ylab = var_names[j1])
plot(beta_chain[, j2], type = "l", main = paste0("Trace: ", var_names[j2]),
     xlab = "Iteration", ylab = var_names[j2])

par(mfrow = c(2,2))
acf(sig2_post, main = expression("ACF: " * sigma^2))
acf(beta_post[, j1], main = paste0("ACF: ", var_names[j1]))
acf(beta_post[, j2], main = paste0("ACF: ", var_names[j2]))

mcmc_sig2  <- mcmc(sig2_post)
mcmc_betas <- mcmc(beta_post)

ess_sig2 <- effectiveSize(mcmc_sig2)
ess_beta <- effectiveSize(mcmc_betas)

cat("ESS(sigma^2) =", ess_sig2, "\n")
cat("ESS(beta):\n")
print(setNames(as.numeric(ess_beta), var_names))

par(mfrow = c(1,1))

# ============================================================
# 7) Posterior table + coefficient CI plot
# ============================================================

summary_bl <- summarize_bl(beta_post, var_names)

compare_tab <- data.frame(
  var = var_names,
  BL_median = summary_bl$BL_median,
  BL_low    = summary_bl$BL_low,
  BL_high   = summary_bl$BL_high,
  LS        = beta_ls,
  Lasso_CV  = beta_lasso_cv,
  Ridge_CV  = beta_ridge_star,
  row.names = NULL
)

print(summary_bl)
print(compare_tab)

ord2 <- order(summary_bl$BL_median)
m <- length(ord2)
xlim_all <- range(c(summary_bl$BL_low, summary_bl$BL_high, beta_ls, beta_lasso_cv))

plot(NA, xlim = xlim_all, ylim = c(1, m),
     xlab = "Standardized coefficients",
     ylab = "Variable (ordered)",
     main = "Bayesian Lasso (lambda*): median and 95% CI (vs LS and Lasso-CV)",
     yaxt = "n")

axis(2, at = 1:m, labels = var_names[ord2], las = 1)
abline(v = 0, lty = 2)

for (i in 1:m) {
  j <- ord2[i]
  segments(summary_bl$BL_low[j], i, summary_bl$BL_high[j], i, lwd = 2)
  points(summary_bl$BL_median[j], i, pch = 4,  cex = 1.2)  # BL median
  points(beta_ls[j],              i, pch = 1,  cex = 1.0)  # LS
  points(beta_lasso_cv[j],        i, pch = 17, cex = 1.0)  # Lasso-CV
}

legend("bottomright",
       legend = c("BL median", "LS", "Lasso (CV)"),
       pch = c(4, 1, 17), bty = "n")

# ============================================================
# 8) Coefficient paths: Lasso vs Bayesian Lasso vs Ridge
# x-axis: ||beta||_1 / max(||beta||_1)
# ============================================================

# Lasso path
s_grid_plot <- seq(0, 1, length.out = 200)
coef_lasso_path <- predict(lasso_path, s = s_grid_plot, type = "coefficients",
                           mode = "fraction")$coefficients

L1_lasso_path <- apply(coef_lasso_path, 1, function(b) sum(abs(b)))
s_lasso_rel <- L1_lasso_path / max(L1_lasso_path)

beta_lasso_cv_check <- as.vector(
  predict(lasso_path, s = s_cv, type = "coefficients", mode = "fraction")$coefficients
)
s_lasso_cv_rel <- sum(abs(beta_lasso_cv_check)) / max(L1_lasso_path)

# Bayesian Lasso path (medians on lambda grid)
s_bl <- L1_vec / max(L1_vec)
ord_s <- order(s_bl)
s_bl_ord <- s_bl[ord_s]
beta_bl_ord <- beta_med_mat[ord_s, , drop = FALSE]
s_bl_star <- s_bl[i_star]

# Ridge path
L1_ridge <- apply(ridge_coef_mat, 1, function(b) sum(abs(b)))
s_ridge_rel <- L1_ridge / max(L1_ridge)
s_ridge_star_rel <- sum(abs(beta_ridge_star)) / max(L1_ridge)

par(mfrow = c(1,3), mar = c(4,4,3,1))

matplot(s_lasso_rel, coef_lasso_path, type = "l", lty = 1,
        xlab = expression("||beta||"[1] / max("||beta||"[1])),
        ylab = "Coefficients",
        main = "(a) Lasso")
abline(v = s_lasso_cv_rel, lty = 2)

matplot(s_bl_ord, beta_bl_ord, type = "l", lty = 1,
        xlab = expression("||beta||"[1] / max("||beta||"[1])),
        ylab = "Coefficients",
        main = "(b) Bayesian Lasso")
abline(v = s_bl_star, lty = 2)

matplot(s_ridge_rel, ridge_coef_mat, type = "l", lty = 1,
        xlab = expression("||beta||"[1] / max("||beta||"[1])),
        ylab = "Coefficients",
        main = "(c) Ridge")
abline(v = s_ridge_star_rel, lty = 2)

par(mfrow = c(1,1))

# ============================================================
# 9) Optional: predictive CV (includes BL and can be expensive)
# ============================================================
run_full_cv <- FALSE

if (run_full_cv) {
  set.seed(100)
  fold_id2 <- fold_id
  
  mse_ls <- numeric(K)
  mse_lasso <- numeric(K)
  mse_ridge <- numeric(K)
  mse_bl <- numeric(K)
  
  M_cv <- 4000
  burn_cv <- 1000
  thin_cv <- 10
  
  for (kk in 1:K) {
    test <- which(fold_id2 == kk)
    train <- setdiff(1:n, test)
    
    Xtr <- X[train, , drop = FALSE]
    ytr <- y[train]
    Xte <- X[test, , drop = FALSE]
    yte <- y[test]
    
    XtX_tr <- crossprod(Xtr)
    Xty_tr <- crossprod(Xtr, ytr)
    
    b_ls <- as.vector(solve(XtX_tr, Xty_tr))
    mse_ls[kk] <- mean((yte - as.vector(Xte %*% b_ls))^2)
    
    b_ridge <- as.vector(solve(XtX_tr + k_star * I_p, Xty_tr))
    mse_ridge[kk] <- mean((yte - as.vector(Xte %*% b_ridge))^2)
    
    lasso_tr <- lars(x = Xtr, y = ytr, type = "lasso", normalize = FALSE, intercept = FALSE)
    b_lasso <- as.vector(predict(lasso_tr, s = s_cv, type = "coefficients", mode = "fraction")$coefficients)
    mse_lasso[kk] <- mean((yte - as.vector(Xte %*% b_lasso))^2)
    
    fit_bl_cv <- run_bl_mcmc(lambda_star, Xtr, ytr, XtX_tr, Xty_tr, a0, b0,
                             M = M_cv, burn = burn_cv, thin = thin_cv,
                             init = NULL, eps_beta = eps_beta, verbose = FALSE)
    b_bl <- apply(fit_bl_cv$beta_post, 2, median)
    mse_bl[kk] <- mean((yte - as.vector(Xte %*% b_bl))^2)
    
    cat("CV fold", kk, "done\n")
  }
  
  cat("\nCV-MSE (mean ± sd)\n")
  cat("LS    :", mean(mse_ls),    "±", sd(mse_ls), "\n")
  cat("Lasso :", mean(mse_lasso), "±", sd(mse_lasso), "\n")
  cat("Ridge :", mean(mse_ridge), "±", sd(mse_ridge), "\n")
  cat("BL    :", mean(mse_bl),    "±", sd(mse_bl), "\n")
}

# ============================================================
# 10) Console summary (for quick reporting)
# ============================================================

cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")

cat(sprintf("n = %d, p = %d\n", n, p))
cat(sprintf("lambda* (L1 matching) = %.6f\n", lambda_star))
cat(sprintf("||median(beta_BL)||_1 = %.4f | ||beta_LassoCV||_1 = %.4f\n",
            L1_vec[i_star], L1_lasso_cv))
cat(sprintf("k* (Ridge-CV) = %.6f\n", k_star))

ess_beta_num <- as.numeric(ess_beta)
names(ess_beta_num) <- var_names
cat(sprintf("ESS(sigma^2) = %.2f\n", as.numeric(ess_sig2)))
cat(sprintf("ESS(beta): min = %.2f | median = %.2f | max = %.2f\n",
            min(ess_beta_num), median(ess_beta_num), max(ess_beta_num)))

active_idx <- which(summary_bl$BL_low > 0 | summary_bl$BL_high < 0)
if (length(active_idx) == 0) {
  cat("95% credible intervals not crossing 0: none\n")
} else {
  act <- summary_bl[active_idx, , drop = FALSE]
  act <- act[order(-abs(act$BL_median)), , drop = FALSE]
  cat("95% credible intervals not crossing 0:\n")
  for (i in 1:nrow(act)) {
    cat(sprintf(" - %s: median = %.4f | CI95%% = [%.4f, %.4f]\n",
                act$var[i], act$BL_median[i], act$BL_low[i], act$BL_high[i]))
  }
}

sign_cmp <- data.frame(
  var = compare_tab$var,
  sign_BL = sign(compare_tab$BL_median),
  sign_LS = sign(compare_tab$LS),
  sign_Lasso = sign(compare_tab$Lasso_CV)
)

sgn_txt <- function(s) ifelse(s > 0, "+", ifelse(s < 0, "-", "0"))
sign_cmp$sign_BL <- sgn_txt(sign_cmp$sign_BL)
sign_cmp$sign_LS <- sgn_txt(sign_cmp$sign_LS)
sign_cmp$sign_Lasso <- sgn_txt(sign_cmp$sign_Lasso)

cat("\nSign check (BL vs LS vs Lasso-CV):\n")
print(sign_cmp, row.names = FALSE)

cat("============================================================\n")
