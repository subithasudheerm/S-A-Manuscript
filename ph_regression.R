# =================================================================
# 1. LOAD DATA & PRE-PROCESSING
# =================================================================
if (!require("survival")) install.packages("survival")
library(survival)

data(jasa)

# Prepare dataframe
my_data <- data.frame(
  time       = jasa$futime,
  age_s      = as.numeric(scale(jasa$age)), 
  surgery    = jasa$surgery,
  transplant = jasa$transplant
)

my_data <- na.omit(my_data)
my_data <- my_data[my_data$time > 0, ]

# =================================================================
# 2. DEFINE THE THETA-WMWPF FUNCTIONS
# =================================================================

S0_wmwpf <- function(t, p) {
  theta <- p[1]; lambda <- p[2]; k <- p[3]; alpha <- p[4]; beta <- p[5]
  if(theta <= 0 || theta >= 1 || any(t >= beta) || any(t <= 0)) return(1e-10)
  S_val <- theta * (1 - (t/beta)^k) + (1 - theta) * exp(-(lambda * t)^alpha)
  return(max(S_val, 1e-10))
}

h0_wmwpf <- function(t, p) {
  theta <- p[1]; lambda <- p[2]; k <- p[3]; alpha <- p[4]; beta <- p[5]
  num <- (theta * k * (t^(k-1)) / (beta^k)) + 
    ((1 - theta) * alpha * (lambda^alpha) * (t^(alpha-1)))
  den <- S0_wmwpf(t, p)
  return(max(num/den, 1e-10))
}

# =================================================================
# 3. FULL LOG-LIKELIHOOD ESTIMATION
# =================================================================

log_lik_reg <- function(params, data) {
  psi <- params[1:5]
  gamma <- params[6:8]
  
  z <- as.matrix(data[, c("age_s", "surgery", "transplant")])
  LP <- z %*% gamma
  
  res <- tryCatch({
    h0_vals <- sapply(data$time, h0_wmwpf, p = psi)
    H0_vals <- sapply(data$time, function(x) -log(S0_wmwpf(x, psi)))
    
    val <- sum(log(h0_vals) + LP - (H0_vals * exp(LP)))
    
    if(!is.finite(val)) return(1e10)
    return(-val) 
  }, error = function(e) 1e10)
  return(res)
}

# Optimization setup
init_p <- c(0.5, 0.001, 1.2, 0.8, max(my_data$time) + 100, 0.1, 0.1, 0.1)

fit_final <- optim(par = init_p, fn = log_lik_reg, data = my_data, 
                   method = "L-BFGS-B", 
                   lower = c(0.01, 1e-6, 0.1, 0.1, max(my_data$time) + 1, -5, -5, -5),
                   upper = c(0.99, 0.1, 15, 15, 10000, 5, 5, 5), 
                   hessian = TRUE)

# =================================================================
# 4. FINAL TABLE GENERATION (Table 5)
# =================================================================

se_full <- sqrt(abs(diag(solve(fit_final$hessian + diag(1e-5, 8)))))

dist_params <- c("theta", "lambda", "k", "alpha", "beta")
dist_est    <- fit_final$par[1:5]
dist_se     <- se_full[1:5]

B      <- fit_final$par[6:8]
SE_B   <- se_full[6:8]
Wald   <- (B / SE_B)^2
Sig    <- 1 - pchisq(Wald, df = 1)
Exp_B  <- exp(B)
L_CI   <- exp(B - 1.96 * SE_B)
U_CI   <- exp(B + 1.96 * SE_B)

# --- Output Formatting with Exp(B) column ---
cat("\nTable 5: MLE for theta-WMWPf Parameters and Covariates\n")
cat(paste(rep("-", 120), collapse=""), "\n")
cat(sprintf("%-15s | %-10s | %-10s | %-10s | %-8s | %-10s | %-16s\n", 
            "Parameter", "Estimate", "Std. Error", "Wald/Z", "Sig.", "Exp(B)", "95% CI (Exp B)"))
cat(paste(rep("-", 120), collapse=""), "\n")

# Distribution Parameters
for(i in 1:5) {
  cat(sprintf("%-15s | %-10.4f | %-10.4f | %-10s | %-8s | %-10s | %-16s\n", 
              dist_params[i], dist_est[i], dist_se[i], "-", "-", "-", "-"))
}

cat(paste(rep(".", 120), collapse=""), "\n")

# Regression Coefficients
labels <- c("Gamma1", "Gamma2", "Gamma3")
for(i in 1:3) {
  cat(sprintf("%-15s | %-10.4f | %-10.4f | %-10.3f | %-8.4f | %-10.4f | (%.3f, %.3f)\n", 
              labels[i], B[i], SE_B[i], Wald[i], Sig[i], Exp_B[i], L_CI[i], U_CI[i]))
}
cat(paste(rep("-", 120), collapse=""), "\n")