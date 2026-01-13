# Load necessary packages
library(readxl)    # For reading Excel files
library(dplyr)     # For data manipulation
library(brms)      # For Bayesian multilevel modeling
library(tidybayes) # For working with Bayesian posteriors
library(ggplot2)   # For visualization
library(emmeans)   # For estimated marginal means (integrates with brms)
library(ggridges)  # For plotting posterior distributions of group means (ridge plots)
library(tidyr)     # For tidyr::separate() and tidyr::pivot_longer
library(bayesplot) # For posterior predictive checks (PPCs)
library(tibble)    # For as_tibble
# install.packages("patchwork") # Uncomment to install if not already installed
library(patchwork) # For combining plots, if desired (optional)

# --- Install and set up cmdstanr if you plan to use backend = "cmdstanr" ---
# This is often faster but requires a separate installation and setup.
# If you prefer to use the default rstan backend, you can comment out
# the 'backend = "cmdstanr"' line in the brm() calls below.

# Step 1: Install cmdstanr package (only run once if not installed)
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# Step 2: Load cmdstanr and install CmdStan (the Stan compiler)
# This will download and compile CmdStan. It may take a few minutes and requires a C++ compiler.
# You only need to run install_cmdstan() once.
# library(cmdstanr) # Uncomment this line to load the package
# install_cmdstan() # Uncomment this line to install CmdStan

# If you still get "CmdStan path has not been set yet" after running install_cmdstan(),
# you might need to explicitly set the path. Replace "path/to/CmdStan" with the actual path
# where CmdStan was installed (e.g., cmdstan_path() will show it after successful install).
# set_cmdstan_path("path/to/CmdStan")


# Set a random seed for reproducibility for all analyses
set.seed(123)

# --- Define Priors (Common for all models) ---
message("Defining common weakly informative priors for all models.")
fixed_effects_prior <- prior(normal(0, 5), class = b)
intercept_prior <- prior(normal(0, 5), class = Intercept)
random_effects_sd_prior <- prior(cauchy(0, 2.5), class = sd)
# sigma_prior is removed for ordinal models as 'sigma' parameter is not directly applicable.
# Priors for 'cutpoints' are handled automatically by brms for the 'cumulative' family.
priors_with_fixed_effects <- c(fixed_effects_prior, intercept_prior, random_effects_sd_prior)
priors_model0 <- c(intercept_prior, random_effects_sd_prior)


# --- Function to Run I-MAIHDA Analysis for a Given Stratum Definition ---
# This function is now more generic, taking a vector of dimension names and a list of their levels.
run_imaihd_analysis <- function(df_original, dimensions_vec, levels_list, analysis_label) {
  # dimensions_vec: A character vector of 3 column names from df_original (e.g., c("Gender", "Position.Level", "Sexuality"))
  # levels_list: A list of 3 character vectors, where each vector contains the levels for the corresponding dimension in dimensions_vec
  
  message(paste0("\n--- Starting Analysis for: ", analysis_label, " ---"))
  
  # Create a copy of the original data and select only relevant columns
  df_data <- df_original %>%
    dplyr::select(all_of(c(dimensions_vec, "Conflict Risk"))) # Select the dimensions and target variable
  
  # Convert categorical variables to factors with specified levels
  for (i in 1:length(dimensions_vec)) {
    col_name <- dimensions_vec[i]
    levels <- levels_list[[i]]
    
    # Handle "Unknown" or "Prefer Not To Answer" by converting to NA and then dropping
    # This ensures these values are not treated as valid levels for analysis
    df_data[[col_name]] <- na_if(df_data[[col_name]], "Unknown")
    df_data[[col_name]] <- na_if(df_data[[col_name]], "Prefer Not To Answer")
    
    # Convert to factor with explicit levels.
    df_data[[col_name]] <- factor(df_data[[col_name]], levels = levels)
  }
  
  # --- CRITICAL CHANGE FOR ORDINAL REGRESSION: Ensure Conflict.Risk is an ORDERED FACTOR ---
  # Assuming Conflict Risk values are numeric codes for ordered categories (e.g., 1, 2, 3, 4, 5)
  # It's crucial that these are treated as ordered, not just nominal or continuous.
  # If Conflict Risk is already factors like "Low", "Medium", "High", then `ordered()` will respect that order.
  df_data$Conflict.Risk <- ordered(as.numeric(df_data$`Conflict Risk`))
  
  # Get the numeric version of Conflict.Risk for PPC stats that require it (e.g., mean, sd)
  # This column is used for calculations and plots that treat the ordinal categories numerically.
  df_data$Conflict.Risk_numeric <- as.numeric(df_data$`Conflict Risk`)
  
  # Create the unique Intersectional_ID for this specific combination
  # Use !!!syms() for unquoting column names in dplyr verbs for interaction
  df_data <- df_data %>%
    mutate(
      Intersectional_ID = interaction(!!!syms(dimensions_vec), sep = "_", drop = TRUE)
    )
  
  # --- CRITICAL FIX: Explicitly remove NAs for modeling and PPCs ---
  # Filter out rows where Conflict.Risk or any of the dimensions are NA.
  # This ensures that the data used for fitting and the 'y' for PPCs are consistent.
  df_cleaned <- df_data %>%
    drop_na(all_of(c("Conflict.Risk", dimensions_vec)))
  
  message(paste0("\nData Preparation for ", analysis_label, ":"))
  message("Number of observations per intersectional stratum (after NA handling):")
  print(df_cleaned %>% group_by(Intersectional_ID) %>% summarise(n = n(), .groups = 'drop'))
  message(paste0("Total observations used for modeling (after NA removal): ", nrow(df_cleaned)))
  
  # Define formulas dynamically based on the dimensions_vec
  formula_base <- "Conflict.Risk ~ 1"
  formula_random_effect <- " + (1 | Intersectional_ID)"
  
  # Formulas for individual fixed effects
  formula_dim1 <- paste0("Conflict.Risk ~ ", dimensions_vec[1], formula_random_effect)
  formula_dim2 <- paste0("Conflict.Risk ~ ", dimensions_vec[2], formula_random_effect)
  formula_dim3 <- paste0("Conflict.Risk ~ ", dimensions_vec[3], formula_random_effect)
  
  # Full additive formula
  formula_additive <- paste0("Conflict.Risk ~ ", paste(dimensions_vec, collapse = " + "), formula_random_effect)
  
  # Full interaction formula
  # This dynamically creates the full three-way interaction term: Dim1 * Dim2 * Dim3
  formula_interaction <- paste0("Conflict.Risk ~ ", paste(dimensions_vec, collapse = " * "), formula_random_effect)
  
  # --- Fit Models (using df_cleaned) ---
  # Model 0: Null Model
  message(paste0("\nFitting Model 0 for ", analysis_label, "..."))
  model0 <- brm(
    as.formula(paste0(formula_base, formula_random_effect)),
    data = df_cleaned, # Use df_cleaned
    family = cumulative("probit"), # --- CHANGED TO PROBIT FAMILY ---
    prior = priors_model0,
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 123,
    backend = "cmdstanr", control = list(adapt_delta = 0.99)
  )
  message("Model 0 fitted.")
  
  # Model 1a: Adjusting for Dim1
  message(paste0("\nFitting Model 1a (", dimensions_vec[1], ") for ", analysis_label, "..."))
  model1a_dim1 <- brm(
    as.formula(formula_dim1),
    data = df_cleaned, # Use df_cleaned
    family = cumulative("probit"), # --- CHANGED TO PROBIT FAMILY ---
    prior = priors_with_fixed_effects,
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 123, backend = "cmdstanr",
    control = list(adapt_delta = 0.99)
  )
  message("Model 1a fitted.")
  
  # Model 1b: Adjusting for Dim2
  message(paste0("\nFitting Model 1b (", dimensions_vec[2], ") for ", analysis_label, "..."))
  model1b_dim2 <- brm(
    as.formula(formula_dim2),
    data = df_cleaned, # Use df_cleaned
    family = cumulative("probit"), # --- CHANGED TO PROBIT FAMILY ---
    prior = priors_with_fixed_effects,
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 123, backend = "cmdstanr",
    control = list(adapt_delta = 0.99)
  )
  message("Model 1b fitted.")
  
  # Model 1c: Adjusting for Dim3
  message(paste0("\nFitting Model 1c (", dimensions_vec[3], ") for ", analysis_label, "..."))
  model1c_dim3 <- brm(
    as.formula(formula_dim3),
    data = df_cleaned, # Use df_cleaned
    family = cumulative("probit"), # --- CHANGED TO PROBIT FAMILY ---
    prior = priors_with_fixed_effects,
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 123, backend = "cmdstanr",
    control = list(adapt_delta = 0.99)
  )
  message("Model 1c fitted.")
  
  # Model 2: Full Additive Model
  message(paste0("\nFitting Model 2 (Full Additive) for ", analysis_label, "..."))
  model2_additive <- brm(
    as.formula(formula_additive),
    data = df_cleaned, # Use df_cleaned
    family = cumulative("probit"), # --- CHANGED TO PROBIT FAMILY ---
    prior = priors_with_fixed_effects,
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 123, backend = "cmdstanr",
    control = list(adapt_delta = 0.99)
  )
  message("Model 2 fitted.")
  
  # Model 3: Full Interaction Model (NEW)
  message(paste0("\nFitting Model 3 (Full Interaction) for ", analysis_label, "..."))
  model3_interaction <- brm(
    as.formula(formula_interaction),
    data = df_cleaned,
    family = cumulative("probit"),
    prior = priors_with_fixed_effects,
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 123, backend = "cmdstanr",
    control = list(adapt_delta = 0.99)
  )
  message("Model 3 fitted.")
  
  # --- Extract Key Metrics ---
  # ICC from Model 0 (Interpreted on the latent scale for ordinal models)
  posterior_samples_model0 <- as_draws_df(model0)
  
  sigma_u0_samples <- posterior_samples_model0$`sd_Intersectional_ID__Intercept`
  variance_u0_samples <- sigma_u0_samples^2
  
  # For cumulative('probit'), the residual variance on the latent scale is fixed at 1
  # This is a common convention to ensure identifiability when using the probit link.
  variance_e0_fixed <- 1 
  
  # Calculate ICC based on the variance of the random intercept and the fixed residual variance
  icc_samples <- variance_u0_samples / (variance_u0_samples + variance_e0_fixed)
  
  # PCV for Full Additive Model
  posterior_samples_model2 <- as_draws_df(model2_additive)
  sigma_u_additive_samples <- posterior_samples_model2$`sd_Intersectional_ID__Intercept`
  variance_u_additive_samples <- sigma_u_additive_samples^2
  pcv_full_additive_samples <- (variance_u0_samples - variance_u_additive_samples) / variance_u0_samples
  
  # Remaining Between-Stratum Variance (from Model 2)
  remaining_variance_samples <- variance_u_additive_samples
  
  # --- Calculate RMSE and Bayesian R-squared for Model 3 (Interaction Model) ---
  message(paste0("\nCalculating RMSE and Bayesian R-squared for Model 3 (", analysis_label, "))..."))
  
  # Get observed y (numeric coding) for calculations
  y_observed_for_metrics <- as.numeric(model3_interaction$data$Conflict.Risk)
  
  # Generate posterior predictions for RMSE calculation (mean predicted category)
  # type="response" gives the predicted category directly (e.g., 1, 2, 3...)
  yrep_for_metrics_categories <- posterior_predict(model3_interaction, newdata = model3_interaction$data, type = "response")
  
  # For RMSE, we'll use the mean predicted category (which is a numeric value)
  # colMeans because yrep_for_metrics_categories is draws x observations
  y_pred_mean_category <- colMeans(yrep_for_metrics_categories)  
  
  # Calculate RMSE based on the predicted numeric categories vs observed numeric categories
  rmse_value <- sqrt(mean((y_observed_for_metrics - y_pred_mean_category)^2))
  
  # Calculate Bayesian R-squared (Conditional R-squared is typically most relevant for multilevel models)
  # For ordinal models, bayes_R2 is interpreted on the latent scale.
  bayes_r2_results <- bayes_R2(model3_interaction)
  
  # --- MODIFIED R-SQUARED EXTRACTION LOGIC (more robust) ---
  if (is.numeric(bayes_r2_results) && length(bayes_r2_results) > 0) {
    r2_values <- bayes_r2_results
    r2_conditional_mean <- mean(r2_values, na.rm = TRUE)
    r2_conditional_lower <- quantile(r2_values, 0.10, na.rm = TRUE)
    r2_conditional_upper <- quantile(r2_values, 0.90, na.rm = TRUE)
  } else if (is.data.frame(bayes_r2_results) && "R2" %in% names(bayes_r2_results)) {
    r2_values <- bayes_r2_results$R2
    r2_conditional_mean <- mean(r2_values, na.rm = TRUE)
    r2_conditional_lower <- quantile(r2_values, 0.10, na.rm = TRUE)
    r2_conditional_upper <- quantile(r2_values, 0.90, na.rm = TRUE)
  } else {
    message("Warning: bayes_R2 did not return a valid numeric vector or data frame with 'R2' column. R-squared will be NA.")
    r2_conditional_mean <- NA
    r2_conditional_lower <- NA
    r2_conditional_upper <- NA
  }
  # --- END MODIFIED R-SQUARED EXTRACTION LOGIC ---
  
  message(paste0("  RMSE (on numeric categories): ", round(rmse_value, 3)))
  message(paste0("  Bayesian R-squared (Conditional, latent scale): ", round(r2_conditional_mean, 3), 
                 " [", round(r2_conditional_lower, 3), ", ", round(r2_conditional_upper, 3), "]"))
  
  # Store results
  results <- list(
    label = analysis_label,
    icc_mean = mean(icc_samples),
    icc_ci_lower = quantile(icc_samples, 0.10),
    icc_ci_upper = quantile(icc_samples, 0.90),
    pcv_full_mean = mean(pcv_full_additive_samples),
    pcv_full_ci_lower = quantile(pcv_full_additive_samples, 0.10),
    pcv_full_ci_upper = quantile(pcv_full_additive_samples, 0.90),
    remaining_var_mean = mean(remaining_variance_samples),
    remaining_var_ci_lower = quantile(remaining_variance_samples, 0.10),
    remaining_var_ci_upper = quantile(remaining_variance_samples, 0.90),
    
    # Add new metrics (from interaction model)
    rmse = rmse_value,
    r2_conditional_mean = r2_conditional_mean,
    r2_conditional_lower = r2_conditional_lower,
    r2_conditional_upper = r2_conditional_upper,
    
    # --- NEW: Add the data with predicted values to the results list ---
    # Store the predicted category for each individual (from interaction model)
    predicted_data = df_cleaned %>%
      mutate(Predicted_Conflict_Risk = factor(
        round(y_pred_mean_category), # Round to nearest integer category
        levels = levels(df_cleaned$Conflict.Risk), # Use original factor levels
        ordered = TRUE
      ))
  )
  
  # --- Posterior Predictive Checks (PPCs) for Model 3 (Interaction Model) ---
  message(paste0("\n--- Performing Posterior Predictive Checks (PPCs) for Model 3 (", analysis_label, ") ---"))
  
  # Debugging: Print number of rows in df_cleaned and model's internal data
  message(paste0("  Number of rows in df_cleaned (original cleaned data): ", nrow(df_cleaned)))
  message(paste0("  Number of rows in model3_interaction$data (data used by model): ", nrow(model3_interaction$data)))
  
  # Use the data actually used by the model for the observed 'y' values
  # For ppc_bars_grouped, y needs to be numeric representing the categories.
  y_observed_for_ppc_numeric <- as.numeric(model3_interaction$data$Conflict.Risk)
  
  # Generate posterior predictions (a few draws for plotting). type = "response" gives categories.
  yrep_for_bayesplot_categories <- posterior_predict(model3_interaction, newdata = model3_interaction$data, draws = 100, type = "response")
  
  # For PPC stat plots (mean, sd), convert yrep categories back to numeric
  yrep_for_bayesplot_numeric <- apply(yrep_for_bayesplot_categories, MARGIN = 2, FUN = as.numeric)
  
  # Debugging: Print dimensions of yrep_for_bayesplot_categories and length of y_observed_for_ppc_numeric
  message(paste0("  Dimensions of yrep_for_bayesplot_categories (rows x columns): ", nrow(yrep_for_bayesplot_categories), " x ", ncol(yrep_for_bayesplot_categories)))
  message(paste0("  Length of y_observed_for_ppc_numeric: ", length(y_observed_for_ppc_numeric)))
  
  if (ncol(yrep_for_bayesplot_categories) != length(y_observed_for_ppc_numeric)) {
    stop(paste0("CRITICAL ERROR: Mismatch in observation counts for PPC. ",
                "Columns in yrep (observations): ", ncol(yrep_for_bayesplot_categories),
                ", Length of observed y: ", length(y_observed_for_ppc_numeric), ". ",
                "This indicates a fundamental issue with posterior_predict or data handling. Please check your data and model fit."))
  }
  
  # Initialize a list to store plots for this analysis label
  current_plots <- list()
  
  # 1. Plot: Distribution of observed vs. predicted (Bars for Ordinal Categories)
  message("  Plotting observed vs. predicted category distribution...")
  plot_bars_grouped <- ppc_bars_grouped(
    y = y_observed_for_ppc_numeric, # --- FIXED: Pass numeric version of observed categories ---
    yrep = yrep_for_bayesplot_categories,
    group = model3_interaction$data$Intersectional_ID, # Groups for plotting
    prob = TRUE # Show proportions
  ) +
    labs(title = paste0("PPC: Observed vs. Predicted Category Proportions (", analysis_label, ")"),
         x = "Conflict Risk Category", y = "Proportion") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels if too crowded
  current_plots[["category_proportions"]] <- plot_bars_grouped
  
  # 2. Plot: Observed vs. Predicted Means (on numeric coding)
  message("  Plotting observed vs. predicted means (on numeric coding)...")
  plot_mean <- ppc_stat(y = y_observed_for_ppc_numeric, yrep = yrep_for_bayesplot_numeric, stat = "mean") +
    labs(title = paste0("PPC: Observed vs. Predicted Mean (", analysis_label, ")"),
         x = "Mean Conflict Risk (Numeric Code)", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  current_plots[["predicted_mean"]] <- plot_mean
  
  # 3. Plot: Observed vs. Predicted Standard Deviations (on numeric coding)
  message("  Plotting observed vs. predicted standard deviations (on numeric coding)...")
  plot_sd <- ppc_stat(y = y_observed_for_ppc_numeric, yrep = yrep_for_bayesplot_numeric, stat = "sd") +
    labs(title = paste0("PPC: Observed vs. Predicted Standard Deviation (", analysis_label, ")"),
         x = "Standard Deviation of Conflict Risk (Numeric Code)", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  current_plots[["predicted_sd"]] <- plot_sd
  
  # 4. Plot: Observed vs. Predicted Group Means (by Intersectional_ID, on numeric coding)
  message("  Plotting observed vs. predicted group means (by Intersectional_ID, on numeric coding)...")
  
  # Calculate observed group means from the data used by the model (model3_interaction$data)
  observed_group_means <- model3_interaction$data %>% # Use model's internal data for observed group means
    group_by(Intersectional_ID) %>%
    summarise(Observed_Mean = mean(as.numeric(Conflict.Risk)), .groups = 'drop') # Use numeric version for mean
  
  # Create a data frame that combines Intersectional_ID with predicted values
  # Each row of yrep_for_bayesplot_categories is a draw of predictions for all observations.
  predicted_data_long <- lapply(1:nrow(yrep_for_bayesplot_categories), function(i) { # Iterate over draws (rows of yrep)
    tibble(
      Intersectional_ID = model3_interaction$data$Intersectional_ID, # Use model's internal data's Intersectional_ID
      Predicted_Value = as.numeric(yrep_for_bayesplot_categories[i, ]), # Select the i-th draw (row), convert to numeric
      Draw_ID = paste0("Draw_", i)
    )
  }) %>%
    bind_rows()
  
  predicted_group_means_long <- predicted_data_long %>%
    group_by(Intersectional_ID, Draw_ID) %>%
    summarise(Predicted_Group_Mean_Per_Draw = mean(Predicted_Value), .groups = 'drop') %>%
    group_by(Intersectional_ID) %>%
    summarise(Mean_of_Predicted_Means = mean(Predicted_Group_Mean_Per_Draw), .groups = 'drop')
  
  # Merge observed and predicted means for plotting
  group_means_plot_data <- left_join(observed_group_means, predicted_group_means_long, by = "Intersectional_ID")
  
  plot_group_means <- ggplot(group_means_plot_data, aes(x = Observed_Mean, y = Mean_of_Predicted_Means)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") + # Ideal line
    labs(title = paste0("PPC: Observed vs. Predicted Group Means (", analysis_label, ")"),
         x = "Observed Group Mean Conflict Risk (Numeric Code)",
         y = "Mean of Predicted Group Mean Conflict Risk (Numeric Code)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_equal() # Ensures x and y axes have the same scale
  current_plots[["group_means"]] <- plot_group_means
  
  message(paste0("--- PPCs complete for Model 3 (", analysis_label, ") ---"))
  
  # Return both results and the list of plots
  return(list(results = results, plots = current_plots))
}

# --- Main Script Execution ---
message("\n--- Loading Original Data for All Analyses ---")
tryCatch({
  # Ensure your Excel file is named "test_data_2.xlsx" and the sheet is "Sheet3"
  df_original_data <- read_excel("test_data_2.xlsx", sheet = "Sheet3")
  message("Original data loaded successfully from 'test_data_2.xlsx' sheet 'Sheet3'.")
}, error = function(e) {
  stop("Error loading data. Make sure 'test_data_2.xlsx' is in your working directory and 'Sheet3' sheet exists, or provide the correct path. Details: ", e$message)
})

# Rename original columns to be R-friendly for consistent passing to function
df_original_data <- df_original_data %>%
  rename(
    Gender = Gender,
    Birth.Country = `Birth Region`, # Note the space here if it exists in your Excel
    Position.Level = Position,
    Sexuality = Sexuality,
    Disability = Disability,
    `Conflict Risk` = `Conflict Risk` # Keep original name for Conflict Risk as it's passed directly
  )


# --- Run Analysis for Each Stratum Definition ---
all_results <- list()
all_plots_list <- list() # New list to store all plots

# 1. Gender x Birth.Country x Position.Level
analysis_output <- run_imaihd_analysis(
  df_original = df_original_data,
  dimensions_vec = c("Gender", "Birth.Country", "Position.Level"),
  levels_list = list(
    c("Man", "Woman"),
    c("American", "Non-American"), # Assuming these are the levels in your data
    c("Junior/Mid Staff", "Management Staff")
  ),
  analysis_label = "Gender x Birth.Country x Position.Level"
)
all_results[["Gender_BirthCountry_PositionLevel"]] <- analysis_output$results
all_plots_list[["Gender_BirthCountry_PositionLevel"]] <- analysis_output$plots

# 2. Gender x Birth.Country x Sexuality
analysis_output <- run_imaihd_analysis(
  df_original = df_original_data,
  dimensions_vec = c("Gender", "Birth.Country", "Sexuality"),
  levels_list = list(
    c("Man", "Woman"),
    c("American", "Non-American"), # Assuming these are the levels in your data
    c("Heterosexual / Straight", "Minority/LGBTQ+") # Ensure these match your data
  ),
  analysis_label = "Gender x Birth.Country x Sexuality"
)
all_results[["Gender_BirthCountry_Sexuality"]] <- analysis_output$results
all_plots_list[["Gender_BirthCountry_Sexuality"]] <- analysis_output$plots

# 3. Gender x Birth.Country x Disability
analysis_output <- run_imaihd_analysis(
  df_original = df_original_data,
  dimensions_vec = c("Gender", "Birth.Country", "Disability"),
  levels_list = list(
    c("Man", "Woman"),
    c("American", "Non-American"), # Assuming these are the levels in your data
    c("Yes", "No") # Ensure these match your data
  ),
  analysis_label = "Gender x Birth.Country x Disability"
)
all_results[["Gender_BirthCountry_Disability"]] <- analysis_output$results
all_plots_list[["Gender_BirthCountry_Disability"]] <- analysis_output$plots

# 4. Gender x Position.Level x Sexuality
analysis_output <- run_imaihd_analysis(
  df_original = df_original_data,
  dimensions_vec = c("Gender", "Position.Level", "Sexuality"),
  levels_list = list(
    c("Man", "Woman"),
    c("Junior/Mid Staff", "Management Staff"),
    c("Heterosexual / Straight", "Minority/LGBTQ+") # Ensure these match your data
  ),
  analysis_label = "Gender x Position.Level x Sexuality"
)
all_results[["Gender_PositionLevel_Sexuality"]] <- analysis_output$results
all_plots_list[["Gender_PositionLevel_Sexuality"]] <- analysis_output$plots

# 5. Gender x Position.Level x Disability
analysis_output <- run_imaihd_analysis(
  df_original = df_original_data,
  dimensions_vec = c("Gender", "Position.Level", "Disability"),
  levels_list = list(
    c("Man", "Woman"),
    c("Junior/Mid Staff", "Management Staff"),
    c("Yes", "No") # Ensure these match your data
  ),
  analysis_label = "Gender x Position.Level x Disability"
)
all_results[["Gender_PositionLevel_Disability"]] <- analysis_output$results
all_plots_list[["Gender_PositionLevel_Disability"]] <- analysis_output$plots


# --- Save Derived Features to CSV for Python ---
message("\n--- Saving Derived I-MAIHDA Features to CSV for Python Integration ---")

for (name in names(all_results)) {
  # Save EMMs if they exist (assuming you might re-add them later)
  if (!is.null(all_results[[name]]$emms_data)) {
    emms_file_name <- paste0("emms_", name, ".csv")
    write.csv(all_results[[name]]$emms_data, file = emms_file_name, row.names = FALSE)
    message(paste0("Saved EMMs for ", name, " to ", emms_file_name))
  }
  # Save RIs if they exist (assuming you might re-add them later)
  if (!is.null(all_results[[name]]$ris_data)) {
    ris_file_name <- paste0("ris_", name, ".csv")
    write.csv(all_results[[name]]$ris_data, file = ris_file_name, row.names = FALSE)
    message(paste0("Saved RIs for ", name, " to ", ris_file_name))
  }
}

# --- NEW: Export individual predicted risk values to CSV ---
message("\n--- Saving Individual Predicted Conflict Risk to CSV ---")

for (name in names(all_results)) {
  if (!is.null(all_results[[name]]$predicted_data)) {
    predicted_data_file_name <- paste0("individual_predictions_", name, ".csv")
    write.csv(all_results[[name]]$predicted_data, file = predicted_data_file_name, row.names = FALSE)
    message(paste0("Saved individual predictions for ", name, " to ", predicted_data_file_name))
  }
}
message("All individual predicted conflict risk data exported.")

# --- Export all plots to a single PDF file ---
pdf("I-MAIHDA_PPC_Visualizations.pdf", width = 10, height = 8) # Adjust width/height as needed for readability
message("\n--- Exporting all visualizations to 'I-MAIHDA_PPC_Visualizations.pdf' ---")

for (analysis_label in names(all_plots_list)) {
  message(paste0("  Exporting plots for: ", analysis_label))
  plots_for_label <- all_plots_list[[analysis_label]]
  for (plot_name in names(plots_for_label)) {
    print(plots_for_label[[plot_name]]) # Print each ggplot object to the PDF
  }
}
dev.off() # Close the PDF device
message("All visualizations exported to 'I-MAIHDA_PPC_Visualizations.pdf'.")

# --- Comparative Analysis of Results (as before) ---
message("\n--- Comparative Analysis of All Intersectional Strata ---")

# Convert list of results to a data frame for easy comparison
# Filter out the 'emms_data' and 'ris_data' elements before binding rows
comparison_df <- bind_rows(lapply(all_results, function(x) {
  # Select only the relevant metrics for the comparison table
  x[c("label", "icc_mean", "icc_ci_lower", "icc_ci_upper",
      "pcv_full_mean", "pcv_full_ci_lower", "pcv_full_ci_upper",
      "remaining_var_mean", "remaining_var_ci_lower", "remaining_var_ci_upper",
      "rmse", "r2_conditional_mean", "r2_conditional_lower", "r2_conditional_upper")]
})) %>%
  mutate(
    ICC = paste0(round(icc_mean, 3), " [", round(icc_ci_lower, 3), ", ", round(icc_ci_upper, 3), "]"),
    PCV_Full_Additive = paste0(round(pcv_full_mean, 3), " [", round(pcv_full_ci_lower, 3), ", ", round(pcv_full_ci_upper, 3), "]"),
    Remaining_Variance = paste0(round(remaining_var_mean, 3), " [", round(remaining_var_ci_lower, 3), ", ", round(remaining_var_ci_upper, 3), "]"),
    RMSE = round(rmse, 3), # Format RMSE
    R2_Conditional = paste0(round(r2_conditional_mean, 3), " [", round(r2_conditional_lower, 3), ", ", round(r2_conditional_upper, 3), "]") # Format R-squared
  ) %>%
  dplyr::select(label, ICC, PCV_Full_Additive, Remaining_Variance, RMSE, R2_Conditional) # Select all columns for the final table

message("\nSummary Table of Key I-MAIHDA Metrics Across Different Strata:")
print(comparison_df)

# --- Interpretation of Non-Additive Intersectionality ---
message("\n--- Interpretation: Which Stratum Exhibits More Non-Additive Intersectionality? ---")
message("To determine which stratum exhibits more non-additive intersectionality, we look for combinations that show:")
message("1. A relatively high ICC from the Null Model (indicating substantial overall group-level variance).")
message("2. A low or negative PCV for the Full Additive Model (meaning additive effects explain little of that group-level variance).")
message("3. A substantial 'Remaining Between-Stratum Variance' (the 'pure' intersectional effect) after accounting for additive effects.")
message("\nBased on the results above, evaluate the 'PCV_Full_Additive' and 'Remaining_Variance' columns:")
message("- A PCV closer to 0 or a negative PCV, especially with a credible interval largely below 0, suggests that the additive effects do NOT explain the group-level variance. This points to strong non-additivity.")
message("- A larger 'Remaining_Variance' (and particularly if its lower CI bound is confidently above zero) indicates a more substantial 'pure' intersectional effect that is not reducible to additive components.")
message("\nPerform a final comparison of the numerical values in the table above to draw your conclusion.")
message("For example, if 'Gender x Birth.Country x Position.Level' has a PCV close to 0 or negative, and a relatively high remaining variance, it would suggest stronger non-additive intersectionality compared to a stratum where PCV is higher (closer to 1) and remaining variance is lower.")
