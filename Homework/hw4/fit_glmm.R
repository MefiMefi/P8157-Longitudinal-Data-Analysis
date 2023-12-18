
# Description -------------------------------------------------------------
# This code is a sequence of customized functions for fitting 
# generalized linear mixed effect model using maximum likelihood 
# with Gauss-Hermite approximation.

# Directions --------------------------------------------------------------
# Press Ctrl/Command + Shift + Enter to load all functions into the environment,
# and excute the custom model fitting 


# Author ------------------------------------------------------------------
# Ryan Wei from P8157 Fall 2023 at Mailman School of Public Health,
# Columbia University


load("../../datasets/ICHS/ICHS.RData")
# Using glmer to validate my method
true_res = lme4::glmer(infection ~ age + xerop +(1|id) ,data = ichs, family = binomial, nAGQ = 3)

sample.dmtx = model.matrix(infection ~ age + xerop ,data = ichs)

# Function to calculate the Hermite polynomial coefficients using the recurrence relation
hermite_poly_coefs <- function(n) {
  if (n == 0) {
    return(c(1))
  } else if (n == 1) {
    return(c(2, 0))
  } else {
    H_n_minus_two <- c(1)
    H_n_minus_one <- c(2, 0)
    for (i in 2:n) {
      H_n <- c(2 * H_n_minus_one, 0) - c(0, 0, 2 * (i - 1) * H_n_minus_two)
      H_n_minus_two <- H_n_minus_one
      H_n_minus_one <- H_n
    }
    return(H_n)
  }
}

# Function to calculate the roots of the Hermite polynomial
hermite_roots <- function(n) {
  coefs <- hermite_poly_coefs(n)
  # Reverse the coefficients for polyroot function
  coefs <- rev(coefs)
  # Find the roots of the polynomial
  roots <- polyroot(coefs)
  # Return only the real parts of the roots
  return(Re(roots))
}

# Function to calculate the weights for the Gauss-Hermite quadrature
gauss_hermite_weights <- function(m) {
  # Calculate the roots of the Hermite polynomial of degree m
  roots <- hermite_roots(m)
  weights <- numeric(m)
  for (j in 1:m) {
    H_m_derivative <- hermite_poly_coefs(m - 1)
    H_m_value <- sum(rev(H_m_derivative) * roots[j]^(0:(m-1)))
    weights[j] <- 2^(m-1) * factorial(m) * sqrt(pi) / (m^2 * H_m_value^2)
  }
  return(weights)
}

# Wrapper
gauss_hermite_point <- function(m){
  x = hermite_roots(m)
  w = gauss_hermite_weights(m)
  return(list(x = x, w = w))
} 

# using fastGHQuad package to check the validity of the weights
# fastGHQuad::gaussHermiteData(5)
# gauss_hermite_point(5)
# two outputs are the same

# Logit function
logit <- function(linear_predictor){
  p = 1/(1 + exp(linear_predictor))
  return(p)
}

# Conditional log_likelihood given random effect, fixed effect for a given subject
conditional_log_likelihood <- function(random_effect, fixed_effects, response, X_matrix) {
  linear_predictor <- X_matrix %*% fixed_effects + random_effect
  p <- logit(linear_predictor)
  log_lik <- dbinom(response, size = rep(1, length(response)), prob = p, log = T)
  return(sum(log_lik))
}

# Gauss-Hermite approximation for the contribution for likelihood of each subject, given estimate of variance component sigma_gamma
gh_approx <- function(fixed_effects, response, X_matrix, gh_points, sigma_gamma){
  integral <- 0
  gh_x = gh_points$x
  gh_w = gh_points$w
  for (i in seq_along(gh_x)) {
    # Scale the node by the random effect standard deviation
    node_scaled <- gh_x[i] * sqrt(2) * sigma_gamma
    # Calculate the log-likelihood contribution for the current node and weight
    ll_contrib <- conditional_log_likelihood(node_scaled, fixed_effects, response, X_matrix)
    # Update the integral approximation
    integral <- integral + exp(ll_contrib) * gh_w[i]
    
  }
  integral_scaled <- integral/ sqrt(pi)
  return(integral_scaled)
}

# log-likelihood for the whole data
# parameters to be optimized
## fixed_effects : beta
## sigma_gamma : variance component
calculate_log_likelihood <- function(data, idVar, responseVar, fixedVars, fixed_effects, sigma_gamma, gh_points) {
  unique_ids <- unique(data[[idVar]])
  total_log_likelihood <- 0
  
  # Ensure fixed_effects is a numeric vector
  fixed_effects <- as.numeric(fixed_effects)
  
  # Pre-calculate the model matrix for fixed effects
  formula <- as.formula(paste(responseVar, "~", paste(fixedVars, collapse = "+")))
  X_matrix <- model.matrix(formula, data)
  
  # Loop over each subject
  for (subj_id in unique_ids) {
    # Subset the data and model matrix for the subject
    subj_indices <- which(data[[idVar]] == subj_id)
    subject_data <- data[subj_indices, ]
    subj_X_matrix <- X_matrix[subj_indices, , drop = FALSE]
    #print(subj_X_matrix)
    #print(fixed_effects)
    
    response <- subject_data[[responseVar]]
    
    # Calculate subject-specific likelihood contribution using gh_approx
    subject_log_likelihood <- gh_approx(fixed_effects, response, subj_X_matrix, gh_points, sigma_gamma)
    
    # Sum the log-likelihood contributions
    total_log_likelihood <- total_log_likelihood + log(subject_log_likelihood)  # Log is taken here as gh_approx returns the scaled integral, not the log likelihood
  }
  
  return(total_log_likelihood)
}

# Optimizer for maximize likelihood
optimize_parameters <- function(data, idVar, responseVar, fixedVars, gh_points) {
  # Initial guesses for parameters (fixed for simplicity)
  initial_fixed_effects <- rep(0, length(fixedVars)+1)
  initial_sigma_gamma <- 1
  #print(initial_fixed_effects)
  # Objective function to minimize (negative log-likelihood)
  objective_func <- function(params) {
    #print(params)
    #print(fixedVars)
    fixed_effects <- params[1:(length(fixedVars)+1)]
    #print(fixed_effects)
    sigma_gamma <- params[(length(fixedVars) + 2)]
    -calculate_log_likelihood(data, idVar, responseVar, fixedVars, fixed_effects, sigma_gamma, gh_points)
  }
  
  # Combine initial guesses into a single vector
  initial_params <- c(initial_fixed_effects, initial_sigma_gamma)
  #print(length(initial_fixed_effects))
  
  # Use optim to maximize the log-likelihood (minimize its negative)
  optim_results <- optim(initial_params, objective_func, method = "BFGS")
  
  return(optim_results)
}

# Formula wrapper
parse_formula <- function(formula) {
  # Extract the response variable
  responseVar <- as.character(formula[[2]])
  
  # Extract the fixed effects variables
  fixed_effects <- attr(terms(formula), "term.labels")
  fixedVar <- fixed_effects[!grepl("\\|", fixed_effects)]
  
  # Identify the random effect and ID variable
  random_effects <- fixed_effects[grepl("\\|", fixed_effects)]
  if (length(random_effects) > 0) {
    randomVar <- stringr::str_trim(gsub("(.*)\\|.*", "\\1", random_effects))
    idVar <- stringr::str_trim(gsub(".*\\|(.*)", "\\1", random_effects))
  } else {
    randomVar <- NA
    idVar <- NA
  }
  
  return(list(responseVar = responseVar, fixedVar = fixedVar, randomVar = randomVar, idVar = idVar))
}


# Model fiting options
fit_glmm <- function(formula, data, nAGQ){
  modelVars = parse_formula(formula)
  response_var = modelVars$responseVar
  fixed_vars = c(modelVars$fixedVar)
  random_var = as.numeric(modelVars$randomVar) # only accepts random intercepts model
  id_var = modelVars$idVar
  
  subject_ids <- unique(data$id_var)
  random_effects <- rep(0, length(subject_ids))
  
  gh_points = gauss_hermite_point(nAGQ)
  
  # Minimizing the negative log-likelihood, need to put a minus sign in front of all estimates
  optim_res = optimize_parameters(data, id_var, response_var, fixed_vars, gh_points)
  fit_res = list(
    fixed_effects_est = -optim_res$par[1:(length(fixed_vars)+1)],
    variance_component_est = optim_res$par[(length(fixed_vars) + 2)],
    loglik = -optim_res$value,
    data = data,
    responseVar = response_var,
    fixedVars = fixed_vars,
    design_mtx = model.matrix(as.formula(paste(response_var, "~", paste(fixed_vars, collapse = "+"))), data),
    idVar = id_var,
    nAGQ = nAGQ,
    gh_points = gh_points,
    formula = formula
  )

  return(fit_res)
}

# Fisher info for estimating variance of fixed effect estimates
calculate_fisher_information <- function(fixed_effects, data, idVar, responseVar, fixedVars, gh_points, sigma_gamma) {
  # Number of fixed effects
  p <- length(fixed_effects)
  
  # Initialize the Fisher Information matrix
  fisher_info <- matrix(0, nrow = p, ncol = p)
  #print(fixed_effects); print(sigma_gamma)
  params = c(fixed_effects)
  # Calculate the Hessian matrix by numerical differentiation
  hessian_func <- function(params) {
    fixed_effects <- params[1:p]
    #print(params)
    calculate_log_likelihood(data, idVar, responseVar, fixedVars, fixed_effects, sigma_gamma, gh_points)
  }
  
  # Calculate the Fisher Information matrix
  # Here, we are using the observed Hessian (second derivative of the log-likelihood)
  fisher_info <- -numDeriv::hessian(hessian_func, params)
  
  return(fisher_info)
}

# This function would be used after obtaining the maximum likelihood estimates
generate_summary <- function(fit_res) {
  fixed_effects_est = fit_res$fixed_effects_est
  variance_component_est = fit_res$variance_component_est
  loglik = fit_res$loglik
  data = fit_res$data
  idVar = fit_res$idVar
  fixedVars = fit_res$fixedVars
  responseVar = fit_res$responseVar
  nAGQ = fit_res$nAGQ
  gh_points = fit_res$gh_points
  formula = fit_res$formula
  # Calculate the Fisher Information matrix
  fisher_info <- calculate_fisher_information(fixed_effects_est, data, idVar, responseVar, fixedVars, gh_points, variance_component_est)
  
  # The estimated covariance matrix of the fixed effects
  #print(fisher_info)
  cov_matrix <- solve(fisher_info)
  #print(cov_matrix)
  # Standard errors of the fixed effects
  std_errors <- sqrt(diag(cov_matrix))
  
  # Calculate AIC and BIC
  n <- nrow(data)
  p <- length(fixed_effects_est)
  AIC <- -2 * loglik + 2 * p
  BIC <- -2 * loglik + log(n) * p
  
  # Calculate deviance as -2 * log-likelihood
  deviance <- -2 * loglik
  
  # Calculate z-values
  z_values <- fixed_effects_est / std_errors
  
  # Calculate p-values
  p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
  
  # Calculate scaled residuals
  formula_fix <- as.formula(paste(responseVar, "~", paste(fixedVars, collapse = "+")))
  X_matrix <- model.matrix(formula_fix, data)
  linear_predictor <- X_matrix %*% fixed_effects_est
  fitted_probs <- plogis(linear_predictor)
  response <- data[[responseVar]]
  raw_residuals <- response - fitted_probs
  variance_function <- fitted_probs * (1 - fitted_probs)
  scaled_residuals <- raw_residuals / sqrt(variance_function)
    
  # Random effects summary
  random_effects_summary <- list(
    Groups = 'id',
    Name = 'Intercept',
    Variance = variance_component_est^2,
    Std.Dev. = variance_component_est
  )
  
  # Fixed effects summary
  fixed_effects_summary <- cbind(Estimate = fixed_effects_est, 
                                 `Std. Error` = std_errors, 
                                 `z value` = z_values, 
                                 `Pr(>|z|)` = p_values)
  
  # Significance codes
  sig_codes <- symnum(p_values, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),  symbols = c("***", "**", "*", ".", " "))
  
  # Correlation of Fixed Effects
  correlation_matrix <- cov_matrix
  rownames(correlation_matrix) <- names(fixed_effects_est)
  colnames(correlation_matrix) <- names(fixed_effects_est)
  
  # Prepare the output similar to summary.glmer
  summary_list <- list(
    call = formula,
    AIC = AIC,
    BIC = BIC,
    logLik = loglik,
    deviance = deviance,
    df.resid = n - p,
    scaled_residuals = scaled_residuals,
    random_effects = random_effects_summary,
    fixed_effects = fixed_effects_summary,
    correlation_fixed_effects = correlation_matrix,
    signif_codes = sig_codes
  )
  
  # Function to print the summary
  print_summary <- function(x) {
    cat("Generalized linear mixed model fit by custom maximum likelihood\n")
    cat(" Family: binomial  ( logit )\n")
    cat("Formula: ", deparse(formula), "\n")
    cat("   Data: ", deparse(substitute(data)), "\n\n")
    
    cat("     AIC      BIC   logLik deviance df.resid \n")
    cat(sprintf("  %6.1f  %6.1f %7.1f  %6.1f  %7d \n\n", AIC, BIC, loglik, deviance, n - p))
    
    cat("Random effects:\n")
    cat(" Groups Name        Variance Std.Dev.\n")
    cat(sprintf(" id     (Intercept) %.4f   %.4f  \n", variance_component_est^2, variance_component_est))
    cat("Number of obs:", n, ", groups:  id, ", length(unique(data$id)), "\n\n")
    
    cat("Fixed effects:\n")
    print(fixed_effects_summary, digits = 3)
    cat("\n")
    cat("Correlation Matrix of Fixed Effects:\n")
    print(correlation_matrix, digits = 3)
    cat("---\n")
    cat("Signif. codes: ", paste(names(attr(sig_codes, "legend")), collapse = " "), "\n")
  }
  
  class(summary_list) <- "custom_glmm_summary"
  attr(summary_list, "print_summary") <- print_summary
  
  return(summary_list)
}


# Print attribute
print.custom_glmm_summary <- function(x) {
  attr(x, "print_summary")(x)
}


# Application
# Inputs are similar to glmer() except for family argument
myfit = fit_glmm(infection ~ age + xerop +(1|id) ,data = ichs, nAGQ = 3)
# A summary function to generate results similar to summary.glmer
mysummary = generate_summary(myfit)
# Model summary
mysummary
