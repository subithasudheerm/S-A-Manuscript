# S-A-Manuscript
## Software Environment & System Dependencies
* **R Version:** 4.5.1 
* **Required CRAN Packages:** * `survival` (>= 3.5-5)
* ## 4. Causal Analysis (`causal_analysis.R`)

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
