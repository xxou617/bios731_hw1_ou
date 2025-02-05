README
================
Xiaxian Ou

# Folder description

*simulations*

`simulations/run_simulations_parallel.R`

- set simulation design elements
- simulation function `simulate_scenario` for one scenario
- parallel calculation by `future, furrr`
- save the simulation results in folder `data`

*source*

`source/01_simulate_data.R`

- `get_simdata` function to get simulation data after setting
  requirements

`source/02_apply_methods.R`

`fit_model`: Use `lm()` when epsilon follows normal distribution and
`glm()` with `gaussian(log)` for lognormal distribution

functions to calculate confidence interval by different methods

- `standard_wald`: Wald confidence intervals
- `boot_quantile`: Nonparametric bootstrap percentile intervals
- `boot_t`: Nonparametric bootstrap t intervals; `boot_quantile` is also
  used inside this function to calculate $\hat{se}(\hat{\theta}^*_b)$
  for each bootstrap estimate

*data*

- simulation results (.Rdata) from each scenario
- the overall results (all_scenarios.Rdata)

*analysis*

- Rmarkdown file to analyze the simulation results
- pdf file rendered from Rmarkdown

*results*

- plots

# Workflow

## main simulation

In `simulations/run_simulations_parallel.R`: Source
`source/01_simulate_data.R` and `source/02_apply_methods.R` to load the
necessary functions into the environment.

### 1. Define Simulation Parameters

The simulation begins by specifying the design elements. These include:

- Sample Sizes: 10, 50, and 500.
- True Treatment Effects ($\beta_{treatment}$): 0, 0.5, and 2.
- Error Distributions: Normally distributed and log-normal.
- CI Methods: Wald, boot quantile, and boot t intervals.

These parameters are combined into a full factorial grid, resulting in
54 unique scenarios. Each scenario is assigned a unique identifier for
tracking.

### 2. Simulation Function

The function `simulate_scenario` is defined in a specific scenario to:

- Simulate Data: Generate datasets tailored to the scenario’s parameters
  using `get_simdata()`.
- Fit Model: Fit a regression model to estimate the treatment effect by
  `lm()` or `glm()` with `gaussian('log)`.
- Compute CIs: Use one of the three CI methods: `standard_wald()`,
  `boot_quantile()` or `boot_t()`. Meanwhile, use `toc` and `tic` to
  calculate the computation time.
- Summarise results: Get bias and coverage
- Save Result: The scenario-specific results are saved as .Rdata files
  for later analysis.

### 3. Parallel Processing

- The future and furrr packages are used to distribute the workload
  across available CPU cores.
- Each scenario is simulated independently in parallel, leveraging
  future_map.

### 4. Run Simulations

Each scenario is simulated with 475 iterations. Random seeds are
generated to ensure reproducibility. The bootstrapping parameters nboot
= 500 for nonparametric bootstrap percentile intervals; nboot_t = 500
for nonparametric bootstrap t intervals with nboot_se = 200 for nested
bootstrapping.

### 5. Combine and Save Results

The results from all scenarios are combined into a single data frame
(final_results) for analysis. Scenario-specific results are stored as
.Rdata files, while the complete dataset is saved as
all_scenarios.Rdata.

## Analysis and Visualization

- The saved results can be used to compare bias, coverage, and
  computation time across methods.
- Visualize performance usnig ggplot2 to identify trends and differences
  across scenarios.
