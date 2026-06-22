# --- 1. LIBRARIES & DATA ---
library(fitdistrplus)
library(goftest)
library(actuar)

# Your specific dataset
data <- c(0.39, 0.81, 0.85, 0.98, 1.08, 1.12, 1.17, 1.18, 1.22, 1.25, 1.36, 1.41,
          1.47, 1.57, 1.57, 1.59, 1.59, 1.61, 1.61, 1.69, 1.69, 1.71, 1.73, 1.80,
          1.84, 1.84, 1.87, 1.89, 1.92, 2.00, 2.03, 2.03, 2.05, 2.12, 2.17, 2.17, 2.17,
          2.35, 2.38, 2.41, 2.43, 2.48, 2.48, 2.50, 2.53, 2.55, 2.55, 2.56, 2.59, 
          2.67, 2.73, 2.74, 2.76, 2.77, 2.79, 2.81, 2.81, 2.82, 2.83, 2.85, 2.87, 
          2.88, 2.93)
n <- length(data)

# --- 2. THETA-WMWPF FUNCTIONS (WITH STABILITY FIX) ---
p_wmwpf <- function(q, theta, lambda, k, alpha, beta) {
  p_val <- 1 - (exp(-(1 - theta) * (lambda * q)^k) * (1 - (q / beta)^alpha)^theta)
  p_val[q >= beta] <- 1
  p_val[q <= 0] <- 0
  # FIX: Ensure strictly in [0, 1]
  return(pmin(pmax(p_val, 0), 1))
}

ll_wmwpf <- function(p, x) {
  th <- p[1]; lam <- p[2]; k <- p[3]; al <- p[4]; be <- p[5]
  if(th <= 0 || th >= 1 || lam <= 0 || k <= 0 || al <= 0 || be <= max(x)) return(1e15)
  t1 <- -(1 - th) * (lam * x)^k
  t2 <- th * log(pmax(1 - (x / be)^al, 1e-15))
  t3 <- log(pmax((1 - th) * k * (lam^k) * (x^(k - 1)) + 
                   th * (al * x^(al - 1)) / (be^al - x^al), 1e-15))
  return(-sum(t1 + t2 + t3))
}

# --- 3. WEIBULL-LOMAX FUNCTIONS (WITH STABILITY FIX) ---
p_wlx <- function(q, theta, lambda, k, alpha, beta_l) {
  term1 <- exp((theta - 1) * (lambda * q)^k)
  term2 <- (1 + beta_l * q)^(-theta * alpha)
  p_val <- 1 - (term1 * term2)
  p_val[q <= 0] <- 0
  # FIX: Ensure strictly in [0, 1]
  return(pmin(pmax(p_val, 0), 1))
}

ll_wlx <- function(p, x) {
  th <- p[1]; lam <- p[2]; k <- p[3]; al <- p[4]; be_l <- p[5]
  if(any(p <= 0)) return(1e15)
  part1 <- exp((th - 1) * (lam * x)^k) * (1 + be_l * x)^(-th * al)
  part2 <- (1 - th) * lam * k * (lam * x)^(k - 1) + (th * al * be_l) / (1 + be_l * x)
  pdf_vals <- pmax(part1 * part2, 1e-15)
  return(-sum(log(pdf_vals)))
}

# --- 4. FIT MODELS ---
f_w  <- fitdist(data, "weibull")
f_g  <- fitdist(data, "gamma")
f_ln <- fitdist(data, "lnorm")
f_ll <- fitdist(data, "llogis")

# Improved starting values for your data
fit_th <- optim(par = c(0.4, 0.2, 1.5, 1.0, max(data) + 1), fn = ll_wmwpf, x = data,
                method = "L-BFGS-B", lower = c(0.01, 0.001, 0.1, 0.1, max(data) + 0.01))

fit_wlx <- optim(par = c(1.1, 0.2, 1.1, 1.0, 0.5), fn = ll_wlx, x = data,
                 method = "L-BFGS-B", lower = rep(0.001, 5))

# --- 5. DYNAMIC CALCULATION FUNCTION ---
calc_stats <- function(dist_name, estimates, logL, k_num, data) {
  aic  <- 2 * k_num - 2 * logL
  bic  <- k_num * log(n) - 2 * logL
  caic <- aic + (2 * k_num * (k_num + 1)) / (n - k_num - 1)
  hqic <- -2 * logL + 2 * k_num * log(log(n))
  
  arg_list <- as.list(estimates)
  ks  <- suppressWarnings(do.call(ks.test, c(list(data, dist_name), arg_list)))
  cvm <- do.call(cvm.test, c(list(data, dist_name), arg_list))$statistic
  ad  <- do.call(ad.test, c(list(data, dist_name), arg_list))$statistic
  
  return(c(LogL=logL, AIC=aic, BIC=bic, CAIC=caic, HQIC=hqic, W_star=cvm, A_star=ad, KS=ks$statistic, P=ks$p.value))
}

# --- 6. COMPILE TABLES ---
results_table <- rbind(
  Weibull      = calc_stats("pweibull", f_w$estimate, f_w$loglik, 2, data),
  Gamma        = calc_stats("pgamma", f_g$estimate, f_g$loglik, 2, data),
  Lognormal    = calc_stats("plnorm", f_ln$estimate, f_ln$loglik, 2, data),
  Log_Logistic = calc_stats("pllogis", f_ll$estimate, f_ll$loglik, 2, data),
  Weibull_Lomax = calc_stats("p_wlx", fit_wlx$par, -fit_wlx$value, 5, data),
  theta_WMWPf  = calc_stats("p_wmwpf", fit_th$par, -fit_th$value, 5, data)
)

# --- 7. PRINT OUTPUT ---
cat("\n--- MLE PARAMETERS ---\n")
cat("Weibull (shape, scale):", round(f_w$estimate, 4), "\n")
cat("Gamma (shape, rate):", round(f_g$estimate, 4), "\n")
cat("Lognormal (meanlog, sdlog):", round(f_ln$estimate, 4), "\n")
cat("Log-Logistic (shape, scale):", round(f_ll$estimate, 4), "\n")
cat("Weibull-Lomax (th, lam, k, al, be):", round(fit_wlx$par, 4), "\n")
cat("theta-WMWPf (th, lam, k, al, be):", round(fit_th$par, 4), "\n")
cat("\n--- GOODNESS OF FIT COMPARISON ---\n")
print(round(results_table, 4))