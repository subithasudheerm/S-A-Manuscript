# --- 1. CLEAN ENVIRONMENT & REPRODUCIBILITY ---
rm(list = ls())
set.seed(123)

# --- 2. DATA GENERATION FUNCTION (Direct CDF Inverse Transform) ---
r_wmwpf <- function(n, theta, lambda, k, alpha, beta) {
  
  # Step 3: Define g(x) using direct CDF formulation
  target_fn <- function(x, u, theta, lambda, k, alpha, beta) {
    # Numerical boundary guards: guarantees standard opposite signs at exact endpoints
    if (x <= 0) return(-u)
    if (x >= beta) return(1 - u)
    
    # Step 2: Cumulative Distribution Function F(x)
    cdf_val <- 1 - (exp(-(1 - theta) * (lambda * x)^k) * (1 - (x / beta)^alpha)^theta)
    
    # g(x) = F(x) - u = 0
    return(cdf_val - u)
  }
  
  samples <- numeric(n)
  for (i in 1:n) {
    # Step 1: Generate independent Ui ~ Uniform(0,1)
    u <- runif(1)
    
    # Step 4: Numerical root-finding over the exact physical support bounds [0, beta]
    res <- uniroot(target_fn, interval = c(0, beta), u = u, 
                   theta = theta, lambda = lambda, k = k, 
                   alpha = alpha, beta = beta)
    samples[i] <- res$root
  }
  return(samples)
}

# --- 3. NEGATIVE LOG-LIKELIHOOD FUNCTION (Hardened against boundary NaN walls) ---
neg_log_lik <- function(p, data) {
  theta <- p[1]; lambda <- p[2]; k <- p[3]; alpha <- p[4]; beta_p <- p[5]
  
  # Parameter constraints for mathematical validity
  if(theta <= 0 || theta >= 1 || lambda <= 0 || k <= 0 || alpha <= 0) return(1e15)
  if(any(data >= beta_p)) return(1e15)
  
  term1 <- -(1 - theta) * (lambda * data)^k
  term2 <- theta * log(pmax(1 - (data / beta_p)^alpha, 1e-15)) # Protected log(0)
  
  denom <- pmax(beta_p^alpha - data^alpha, 1e-15)            # Protected denominator
  term3 <- log(pmax((1 - theta) * k * (lambda^k) * (data^(k - 1)) + 
                      theta * (alpha * data^(alpha - 1)) / denom, 1e-15))
  
  return(-sum(term1 + term2 + term3))
}

# --- 4. SIMULATION SETUP ---
n_size <- 500            # Manually update to 50, 100, 500 for your other tables
replications <- 1000
true_params <- c(theta=0.60, lambda=0.15, k=2.00, alpha=1.20, beta=5.00)

results <- matrix(NA, nrow = replications, ncol = 5)
colnames(results) <- names(true_params)

# --- 5. EXECUTION LOOP ---
cat("Starting simulation for n =", n_size, "...\n")
for (i in 1:replications) {
  # Generate sample using the hardened root finder
  sample_data <- r_wmwpf(n_size, true_params[1], true_params[2], 
                         true_params[3], true_params[4], true_params[5])
  
  # Initial Guess adjustments
  start_vals <- true_params * 1.15
  start_vals[5] <- max(sample_data) + 0.1  
  
  fit <- optim(par = start_vals, fn = neg_log_lik, data = sample_data,
               method = "L-BFGS-B", 
               lower = c(0.01, 0.01, 0.01, 0.01, max(sample_data) + 1e-4), 
               upper = c(0.99, 10.0, 10.0, 10.0, 20.0))
  
  if (fit$convergence == 0) {
    results[i, ] <- fit$par
  }
}

# --- 6. CALCULATE PERFORMANCE METRICS ---
clean_results <- results[complete.cases(results), ]

# Average MLE: Mean of the estimates
mle_means <- colMeans(clean_results)

# Standard Error (SE): Standard deviation of the estimates across replications
se_vals <- apply(clean_results, 2, sd)

# Mean Square Error (MSE): Mean of (Estimate - True)^2
mse_vals <- colMeans(sweep(clean_results, 2, true_params)^2)

# --- 7. FINAL TABLE FOR PUBLICATION ---
summary_table <- data.frame(
  Parameter   = names(true_params),
  True_Value  = as.numeric(true_params),
  Average_MLE = round(as.numeric(mle_means), 4),
  SE          = round(as.numeric(se_vals), 4),
  MSE         = round(as.numeric(mse_vals), 4)
)

cat("\nSimulation Results for n =", n_size, "\n")
print(summary_table)