# S-A-Manuscript
## Software Environment & System Dependencies
* **R Version:** 4.5.1 
* **Required CRAN Packages:** * `survival` (>= 3.5-5)
* ## Causal Analysis (`causal_analysis.R`)

### Description
This script implements a G-computation (Causal Standardization) analysis pipeline using the baseline demographics of the Stanford Heart Transplant dataset (`jasa`). It evaluates the Average Treatment Effect (ATE) of receiving a heart transplant on patient survival probability across standardized, discrete time horizons ($t = 100, 500,$ and $1000$ days).

### Key Features & Workflow:
1. **Data Pre-processing:** Extracts the survival timeline, scales the age covariate, and removes missing observations to build a clean analytical dataframe.
2. **Counterfactual Matrix Generation:** Fits a covariate-adjusted Cox Proportional Hazards regression model (`coxph`). It then clones the baseline data to create two separate counterfactual scenarios:
   * **Treated Scenario ($S^{(1)}$):** A structural paradigm where *every* patient hypothetically receives a transplant.
   * **Control Scenario ($S^{(0)}$):** A structural paradigm where *no* patient receives a transplant.
3. **Robust Array Alignment (`extract_and_pad_surv`):** Features an error-insulated helper function designed to safely aggregate individual counterfactual curves via marginal means (`rowMeans`) and align them to exact time milestones, preventing vector length mismatches.
4. **ATE Calculation:** Compares the marginal survival probabilities to compute the absolute risk difference ($\text{ATE} = \hat{S}^{(1)}(t) - \hat{S}^{(0)}(t)$).

### Required R Packages
* **`survival` (v3.5-5 or higher):** Core package used for loading the `jasa` dataset, fitting the Cox model via `coxph()`, and extracting adjusted survival pathways using `survfit()`.
* **`stats` (Base R):** Internal package utilized for data cleaning via `na.omit()` and standardizing covariates with `scale()`.
* **`base` (Base R):** Utilized for operational array checking (`rowMeans()`, `match()`) and terminal printing.
## Regression Analysis and Maximum Likelihood Estimation (`ph_regression.R`)

### Description
This script implements a parametric proportional hazards regression framework using the baseline demographics of the Stanford Heart Transplant dataset (`jasa`). It maximizes the full log-likelihood function under a baseline $\theta$-WMWPf survival distribution to simultaneously estimate baseline distributional parameters ($\theta, \lambda, k, \alpha, \beta$) and covariate effects ($\gamma_1, \gamma_2, \gamma_3$).

### Key Features & Workflow:
1. **Data Integration & Transformation:** Loads the baseline `jasa` cohort data. It isolates the operational timelines (`futime`), transforms and standardizes the continuous age covariate (`scale`), and matches structural indicator dummy variables for history of surgery and transplant status.
2. **Parametric Formulation:** Code definitions are written directly for the structural baseline survival function ($S_0(t)$) and its derivative baseline hazard mechanism ($h_0(t)$) dictated by the parameters bounded inside the valid space $(0 < \theta < 1, \lambda > 0, k > 0, \alpha > 0, \beta > \max(t))$.
3. **Log-Likelihood Maximization:** Establishes a joint log-likelihood interface mapping the linear prediction engine ($z^T\gamma$) against cumulative hazard matrices ($H_0(t) = -\log S_0(t)$). Numerical optimization is handled using the bounded L-BFGS-B algorithm (`optim`), integrating safe boundary handlers (`tryCatch`) to insulate calculation steps against invalid parameters.
4. **Statistical Inference & Wald Tests:** Infers asymptotic standard errors using a regularized observed Hessian inversion matrix (`solve`). It computes standard Wald chi-square scores ($\text{Wald} = (\hat{\gamma}/\text{SE})^2$) and yields significant $p$-values, hazard ratios ($\text{Exp}(\beta)$), and symmetric 95% exponential confidence bounds.
5. **Table Generation (Table 10):** Automatically formats and prints an academic-standard summary text block matching the precise structural presentation of Table 10 in the manuscript.

### Required R Packages
* **`survival` (v3.5-5 or higher):** Core package used for downloading the historical `jasa` dataset.

* **`stats` (Base R):** Provides the mathematical optimization interface via `optim(method = "L-BFGS-B")` and standard chi-squared
* distribution computing via `pchisq()`.
* **`base` (Base R):** Utilized for structural linear algebra calculations (`%*%`, `diag`, `solve`) and element mappings (`sapply`).
* ## Sampling and Monte Carlo Simulation (`sampling_and_simulation.R`)

### Description
This script implements an automated Monte Carlo simulation framework to evaluate the finite-sample performance of Maximum Likelihood Estimators (MLE) for the $\theta$-WMWPf distribution. It establishes a random variate generator, evaluates parameter estimation across a series of replications, and outputs core performance benchmarks (Average MLE, Standard Error, and Mean Squared Error).

### Key Features & Workflow:
1. **Inverse Transform Sampling Framework (`r_wmwpf`):** Generates independent random realizations from the $\theta$-WMWPf distribution using the Probability Integral Transform. Since the cumulative distribution function (CDF) cannot be inverted analytically, it implements an independent random uniform draw $U_i \sim \text{Uniform}(0,1)$ and executes numerical root-finding using `uniroot` bounded across the exact support interval $[0, \beta]$.
2. **Hardened Estimation Objective Matrix (`neg_log_lik`):** Evaluates the negative joint log-likelihood function. It includes strict numerical boundary guards and protected logarithm thresholds (`pmax`) to insulate the execution loop against boundary $\text{NaN}$ traps and division-by-zero errors when sample values draw close to the parameter ceiling $\beta$.
3. **Execution Loop & Optimization Setup:** Simulates data generations over a specified number of independent replications (e.g., $1,000$ iterations) given a target baseline parameter vector $\Theta_0 = (\theta = 0.60, \lambda = 0.15, k = 2.00, \alpha = 1.20, \beta = 5.00)$. Parameter optimization is calculated via the bounded L-BFGS-B algorithm (`optim`).
4. **Computational Efficiency Protocol:** To optimize computational efficiency, the Monte Carlo simulations are executed sequentially for each sample size ($n = 20, 50, 100,$ and $500$) rather than simultaneously. This serial configuration prevents memory saturation and mitigates excessive computational overhead on standard hardware structures.
5. **Statistical Metrics & Output Compiler:** Filters out un-converged optimization runs using strict complete-case filtering. It processes the empirical sampling matrix to compile statistical metrics:
   * **Average MLE:** $\frac{1}{B}\sum_{b=1}^{B} \hat{\psi}_b$
   * **Standard Error (SE):** The empirical standard deviation of parameter estimations across valid iterations.
   * **Mean Squared Error (MSE):** $\frac{1}{B}\sum_{b=1}^{B} (\hat{\psi}_b - \psi_{\text{true}})^2$
   It automatically formats and outputs a publication-ready data frame containing these metrics.

### Required R Packages
* **`stats` (Base R):** Handles inverse uniform generation (`runif`), physical root evaluations (`uniroot`), bounded optimization processing (`optim(method = "L-BFGS-B")`), and structural statistical metrics (`sd`).
* **`base` (Base R):** Utilized for structural array loops, array dimension cleaning (`complete.cases`), matrix operations (`colMeans`, `sweep`), and formatted console outputs.
* ##  Empirical Real Data Application (`Dataset1- Table 6&7.R`)

### Description
This script conducts a comprehensive real-data application and comparative goodness-of-fit analysis on an empirical dataset. It evaluates the fitting flexibility of the proposed bounded $\theta$-WMWPf distribution against five competing models: the standard baseline classical models (Weibull, Gamma, Lognormal, Log-Logistic) and a advanced five-parameter generalization (Weibull-Lomax distribution).

### Key Features & Workflow:
1. **Data Initialization:** Loads and analyzes a continuous survival/reliability dataset ($n = 50$). 
2. **Probability Distribution Modeling:** Code definitions are written directly for the cumulative distribution functions (CDF) and log-likelihood engines of both the bounded $\theta$-WMWPf distribution and the complex Weibull-Lomax ($\text{WLx}$) model.
3. **Bounded Multi-Model Estimation:** Leverages the `fitdistrplus` architecture to fit baseline distributions via classical maximum likelihood estimation. For the highly non-linear, parameter-constrained $\theta$-WMWPf and Weibull-Lomax configurations, custom optimizations are computed using the `L-BFGS-B` algorithm (`optim`), integrating rigid boundary thresholds to restrict the mixing parameter space strictly within $0 < \theta < 1$.
4. **Goodness-of-Fit Metric Compilations (`calc_stats`):** A unified automated helper routine that parses the resulting maximized log-likelihood configurations to calculate seven critical academic selection benchmarks:
   * **Information Criteria:** Akaike Information Criterion ($\text{AIC}$), Bayesian Information Criterion ($\text{BIC}$), Consistent Akaike Information Criterion ($\text{CAIC}$), and Hannan-Quinn Information Criterion ($\text{HQIC}$).
   * **Empirical Distance Tests:** Cramer-von Mises statistic ($W^*$), Anderson-Darling statistic ($A^*$), and the Kolmogorov-Smirnov ($\text{KS}$) distance alongside its corresponding $p$-value.
5. **Table Generation (Tables 6 & 7):** Merges the statistical metrics into an empirical evaluation matrix. It automatically formats and outputs the exact values required to compile both **Table 6 (MLE Parameters)** and **Table 7 (Goodness of Fit Metrics)** in the final manuscript.

### Required R Packages
* **`fitdistrplus` (v1.1-11 or higher):** Utilized for numerical parameter estimations of standard continuous baseline alternatives (`fitdist`).
* **`actuar` (v3.3-3 or higher):** Extends baseline sampling by providing the Log-Logistic (`llogis`) operational primitives.
* **`goftest` (v1.2-3 or higher):** Provides the mathematical engines for advanced empirical distance testing via Cramer-von Mises (`cvm.test`) and Anderson-Darling (`ad.test`).
* **`stats` (Base R):** Handles core optimizations, the Kolmogorov-Smirnov infrastructure (`ks.test`), and distribution computations (`pweibull`, `pgamma`, `plnorm`).
* ## 6. Empirical Application on Carbon Fiber Strengths (`Dataset2- Table 8&9.R`)

### Description
This script handles the second empirical application in the manuscript, fitting the bounded $\theta$-WMWPf distribution and competing models to a reliability dataset representing the single-filament strengths of carbon fibers ($n = 63$). It compares the bounded distribution's fitting capabilities against classical models (Weibull, Gamma, Lognormal, Log-Logistic) and the generalized five-parameter Weibull-Lomax distribution.

### Key Features & Workflow:
1. **Data Initialization:** Loads the empirical carbon fiber tensile strength vector ($n = 63$).
2. **Probability Distribution Framework with Stability Guards:** Implements the custom cumulative distribution functions (CDF) and log-likelihood calculation engines. The functions feature numerical bounding constraints (`pmin` and `pmax` down to `1e-15`) to handle precision floors and ensure stability when working with small fractions close to zero or the parameter boundary $\beta$.
3. **Parametric Model Fitting:** Executes automated parametric fitting. Standard continuous distributions are calculated via the `fitdistrplus` package. Non-linear, complex multi-parameter spaces ($\theta$-WMWPf and Weibull-Lomax) are solved using the bounded `L-BFGS-B` optimization routine (`optim`) with updated safe initial guesses tailored specifically for this dataset.
4. **Comprehensive Goodness-of-Fit Comparison (`calc_stats`):** Automatically computes seven essential model-selection selection criteria based on the optimized, maximized log-likelihood value:
   * **Information Criteria:** Akaike Information Criterion ($\text{AIC}$), Bayesian Information Criterion ($\text{BIC}$), Consistent Akaike Information Criterion ($\text{CAIC}$), and Hannan-Quinn Information Criterion ($\text{HQIC}$).
   * **Distance and Hypothesis Tests:** Cramer-von Mises ($W^*$), Anderson-Darling ($A^*$), and Kolmogorov-Smirnov ($\text{KS}$) distance metrics along with the resulting asymptotic $p$-values.
5. **Table Generation (Tables 8 & 9):** Combines the parameters and metrics into a summary matrix. It outputs the numerical results required to construct **Table 8 (MLE Parameters for Dataset 2)** and **Table 9 (Goodness of Fit Metrics for Dataset 2)** in the final manuscript.

### Required R Packages
* **`fitdistrplus` (v1.1-11 or higher):** Used to perform automatic maximum likelihood estimations on standard continuous distributions (`fitdist`).
* **`actuar` (v3.3-3 or higher):** Extends baseline distribution families to include the Log-Logistic (`llogis`) baseline distribution.
* **`goftest` (v1.2-3 or higher):** Computes empirical distance tests via Cramer-von Mises (`cvm.test`) and Anderson-Darling (`ad.test`) functions.
* **`stats` (Base R):** Coordinates bounded multivariate optimization routines, basic data summary measures, and hypothesis checking (`ks.test`).
