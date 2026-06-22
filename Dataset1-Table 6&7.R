# --- 1. LIBRARIES & DATA ---
library(fitdistrplus)
library(goftest)
library(actuar)

data <- c(0.1, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 6.0,
          7.0, 11.0, 12.0, 18.0, 18.0, 18.0, 18.0, 18.0, 21.0, 32.0,
          36.0, 40.0, 45.0, 45.0, 47.0, 50.0, 55.0, 60.0, 63.0, 63.0,
          67.0, 67.0, 67.0, 67.0, 72.0, 75.0, 79.0, 82.0, 82.0, 83.0,
          84.0, 84.0, 84.0, 85.0, 85.0, 85.0, 85.0, 85.0, 86.0, 86.0)
n <- length(data)

# --- 2. THETA-WMWPF FUNCTIONS ---
p_wmwpf <- function(q, theta, lambda, k, alpha, beta) {
  p_val <- 1 - (exp(-(1 - theta) * (lambda * q)^k) * (1 - (q / beta)^alpha)^theta)
  p_val[q >= beta] <- 1
  p_val[q <= 0] <- 0
  return(p_val)
}

ll_wmwpf <- function(p, x) {
  th <- p[1]; lam <- p[2]; k <- p[3]; al <- p[4]; be <- p[5]
  if(th <= 0 || th >= 1 || lam <= 0 || k <= 0 || al <= 0 || be <= max(x)) return(1e15)
  t1 <- -(1 - th) * (lam * x)^k
  t2 <- th * log(1 - (x / be)^al)
  t3 <- log((1 - th) * k * (lam^k) * (x^(k - 1)) + 
              th * (al * x^(al - 1)) / (be^al - x^al))
  return(-sum(t1 + t2 + t3))
}

# --- 3. WEIBULL-LOMAX FUNCTIONS ---
p_wlx <- function(q, theta, lambda, k, alpha, beta_l) {
  # CDF: 1 - exp((theta-1)*(lambda*q)^k) * (1 + beta_l*q)^(-theta*alpha)
  term1 <- exp((theta - 1) * (lambda * q)^k)
  term2 <- (1 + beta_l * q)^(-theta * alpha)
  p_val <- 1 - (term1 * term2)
  p_val[q <= 0] <- 0
  return(p_val)
}

ll_wlx <- function(p, x) {
  th <- p[1]; lam <- p[2]; k <- p[3]; al <- p[4]; be_l <- p[5]
  # CRITICAL FIX: Ensure theta is strictly between 0 and 1, and others are positive
  if(th <= 0 || th >= 1 || lam <= 0 || k <= 0 || al <= 0 || be_l <= 0) return(1e15)
  
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

# Fitting theta-WMWPf
fit_th <- optim(par = c(0.4, 0.02, 1.5, 1.0, max(data)), fn = ll_wmwpf, x = data,
                method = "L-BFGS-B", lower = c(0.01, 0.001, 0.1, 0.1, max(data) + 0.1))

# Fitting Weibull-Lomax
# Corrected Fitting for Weibull-Lomax with strict bounds on theta
fit_wlx <- optim(par = c(0.5, 0.02, 0.5, 1.2, 0.1), # Changed starting theta to 0.5
                 fn = ll_wlx, 
                 x = data,
                 method = "L-BFGS-B", 
                 lower = c(0.001, 0.001, 0.001, 0.001, 0.001), # Lower boundaries
                 upper = c(0.999,   Inf,   Inf,   Inf,   Inf)) # Upper boundary enforced for theta!

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
cat("\n--- MLE PARAMETERS (WLx vs WMWPf) ---\n")
# --- 6. PRINT OUTPUT ---
cat("\n--- MLE PARAMETERS ---\n")
cat("Weibull (shape, scale):", round(f_w$estimate, 4), "\n")
cat("Gamma (shape, rate):", round(f_g$estimate, 4), "\n")
cat("Lognormal (meanlog, sdlog):", round(f_ln$estimate, 4), "\n")
cat("Log-Logistic (shape, scale):", round(f_ll$estimate, 4), "\n")
cat("Weibull-Lomax (th, lam, k, al, be):", round(fit_wlx$par, 4), "\n")
cat("theta-WMWPf (th, lam, k, al, be):", round(fit_th$par, 4), "\n")

cat("\n--- GOODNESS OF FIT COMPARISON ---\n")
print(round(results_table, 4))