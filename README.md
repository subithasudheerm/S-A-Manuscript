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
5. **Table Generation (Table 5):** Automatically formats and prints an academic-standard summary text block matching the precise structural presentation of Table 5 in the manuscript.

### Required R Packages
* **`survival` (v3.5-5 or higher):** Core package used for downloading the historical `jasa` dataset.
* **`stats` (Base R):** Provides the mathematical optimization interface via `optim(method = "L-BFGS-B")` and standard chi-squared distribution computing via `pchisq()`.
* **`base` (Base R):** Utilized for structural linear algebra calculations (`%*%`, `diag`, `solve`) and element mappings (`sapply`).
