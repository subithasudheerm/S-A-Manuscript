# =====================================================================
# 1. LOAD DATA & PRE-PROCESSING (Stanford Heart Transplant Data)
# =====================================================================
# Correct: Standard single-nested quotes
if (!require("survival")) install.packages("survival")

library(survival)

data(jasa)

# Prepare clean data frame matching your regression variables
my_data <- data.frame(
  time       = jasa$futime,
  age_s      = as.numeric(scale(jasa$age)), 
  surgery    = jasa$surgery,
  transplant = jasa$transplant # This acts as our Treatment indicator (0 or 1)
)

# Clean missing records and drop non-positive times
my_data <- na.omit(my_data)
my_data <- my_data[my_data$time > 0, ]

# =====================================================================
# 2. FIT COVARIATE-ADJUSTED COX MODEL (To generate survival matrices)
# =====================================================================
# Fit a Cox model adjusting for baseline demographics
cox_fit <- coxph(Surv(time) ~ transplant + age_s + surgery, data = my_data)

# Create counterfactual datasets:
# What if EVERYONE received a transplant vs. What if NO ONE received a transplant
data_all_treated <- my_data
data_all_treated$transplant <- 1

data_all_control <- my_data
data_all_control$transplant <- 0

# Generate the adjusted multi-column survival curves
surv_treated <- survfit(cox_fit, newdata = data_all_treated)
surv_control <- survfit(cox_fit, newdata = data_all_control)

# =====================================================================
# 3. DEFINE TIME POINTS OF INTEREST
# =====================================================================
time_points_interest <- c(100, 500, 1000)

# Extract survival summaries at the specified intervals
summary_S1 <- summary(surv_treated, times = time_points_interest)
summary_S0 <- summary(surv_control, times = time_points_interest)

# =====================================================================
# --- YOUR SAFE DUAL-FUNCTION HELPER ---
# =====================================================================
extract_and_pad_surv <- function(summary_obj, requested_times) {
  # Guard against completely empty summary objects
  if (is.null(summary_obj$time) || length(summary_obj$time) == 0) {
    return(rep(NA, length(requested_times)))
  }
  
  # Handle 2D matrix (adjusted curves) vs 1D vector (Kaplan-Meier) dynamically
  if (!is.null(ncol(summary_obj$surv)) && ncol(summary_obj$surv) > 1) {
    raw_values <- rowMeans(summary_obj$surv)
  } else {
    raw_values <- as.vector(summary_obj$surv)
  }
  
  # Align values to requested time points to fix any row mismatch errors
  aligned_vector <- rep(NA, length(requested_times))
  matching_indices <- match(summary_obj$time, requested_times)
  
  valid_matches <- !is.na(matching_indices)
  if (any(valid_matches)) {
    aligned_vector[matching_indices[valid_matches]] <- raw_values[valid_matches]
  }
  
  return(aligned_vector)
}

# =====================================================================
# 4. CALCULATE CAUSAL EFFECT (ATE) SAFELY
# =====================================================================
S1_values     <- extract_and_pad_surv(summary_S1, time_points_interest)
S0_values     <- extract_and_pad_surv(summary_S0, time_points_interest)
causal_effect <- S1_values - S0_values

# =====================================================================
# 5. CREATE AND DISPLAY THE FORMATTED DATA FRAME
# =====================================================================
causal_table <- data.frame(
  `Days (t)`   = time_points_interest,
  `ˆ S(1)(t)`  = round(S1_values, 4),
  `ˆ S(0)(t)`  = round(S0_values, 4),
  `ATE`        = round(causal_effect, 4),
  check.names  = FALSE # Preserves your exact academic math column names
)

cat("\nTable: Calculated Causal Effect at Specific Time Intervals (Stanford Heart Data)\n")
print(causal_table, row.names = FALSE)