load("../../datasets/ICHS/ICHS.RData")
# function for separating formula
lme4::glmer(infection ~ age + xerop +(1|id) ,data = ichs, family = binomial, nAGQ = 3)

sample.dmtx = model.matrix(infection ~ age + xerop ,data = ichs)

hermite.poly <- function(n) {
  Hn <- function(x) {
    (-1)^n * exp(x^2) * stats::deriv(expr = exp(-x^2), name = "x", order = n)$value
  }
  return(Hn)
}

hermite.poly(2)

# Function to calculate the nodes and weights for Gauss-Hermite quadrature
gh_nodes_weights <- function(n) {
  # Hermite polynomial of degree n
  Hn <- function(x) {
    hermite.poly(n)(x)
  }
  
  # Derivative of Hermite polynomial of degree n
  Hn_prime <- function(x) {
    hermite.poly(n)(x)
  }
  
  # Find roots of Hermite polynomial, which are the nodes
  nodes <- as.numeric(polyroot(c(rep(0, n), 1)))
  
  # Calculate weights
  # w_j = 2^(n-1) * n! * sqrt(pi) / (n^2 * (H'_n(x_j))^2)
  weights <- 2^(n-1) * factorial(n) * sqrt(pi) / (n^2 * (Hn_prime(nodes))^2)
  
  return(list(nodes = nodes, weights = weights))
}

gh_nodes_weights(4)

gauss_hermite_quadrature <- function(h, n) {
  # Generate the Hermite polynomial of degree n
  Hn <- function(x) {
    (-1)^n * exp(x^2) * stats::deriv(expr = exp(-x^2), name = "x", order = n)$value
  }
  
  # Find the roots of the Hermite polynomial
  # Hermite polynomials are symmetric about the origin, so we only need to find positive roots
  hermite_roots <- Re(polyroot(c(rep(0, n), 1)))
  hermite_roots <- hermite_roots[hermite_roots >= 0]
  
  # Define Hermite polynomial of degree n-1
  Hn_minus_one <- function(x) {
    (-1)^(n - 1) * exp(x^2) * stats::deriv(expr = exp(-x^2), name = "x", order = n - 1)$value
  }
  
  # Calculate the weights for each root
  weights <- 2^(n - 1) * factorial(n) * sqrt(pi) / (n^2 * Hn_minus_one(hermite_roots)^2)
  
  # Approximate the integral using the roots and weights
  integral_approx <- sum(sapply(hermite_roots, function(x) h(x) * exp(x^2))) * weights
  
  return(integral_approx)
}

# Define your function h(x)
h <- function(x) { 1 }  # Example function, replace with your actual function h(x)

# Specify the number of roots (degree of Hermite polynomial)
n <- 10

# Perform Gauss-Hermite quadrature to approximate the integral of h(x)
integral_approximation <- gauss_hermite_quadrature(h, n)
print(integral_approximation)


# Define the negative log-likelihood function for a logistic regression with random intercept
neg_log_likelihood <- function(beta, y, X, groups) {
  # Calculate the linear predictors
  eta <- X %*% beta
  
  # Calculate the log-likelihood of the fixed effects
  loglik_fixed <- sum(y * eta - log(1 + exp(eta)))
  
  # Placeholder for the random intercept effects
  # In a real implementation, you would integrate over the random effects distribution
  # For simplicity, we will just assume a single random intercept with a mean of zero
  # and a standard deviation estimated from the data
  
  
  sigma <- sd(groups)
  random_intercepts <- rnorm(length(unique(groups)), 0, sigma)
  loglik_random <- sum(dnorm(random_intercepts, 0, sigma, log = TRUE))
  
  # Calculate the combined log-likelihood
  loglik <- loglik_fixed + loglik_random
  
  # Return the negative log-likelihood
  return(-loglik)
}

# Define the function to fit the GLMM
fit_glmm <- function(y, X, groups) {
  # Starting values for the fixed effects parameters
  start_val <- rep(0, ncol(X))
  
  # Optimize the negative log-likelihood function
  opt_res <- optim(par = start_val, fn = neg_log_likelihood, y = y, X = X, groups = groups)
  
  # Return the optimization results
  return(opt_res)
}

# Example usage:
# y is the binary outcome vector
# X is the model matrix for the fixed effects (including the intercept)
# groups is the vector representing the grouping factor for the random intercepts

# Example data
y <- c(0, 1, 1, 0, 1)  # Binary outcomes
X <- cbind(1, matrix(rnorm(10), ncol = 2))  # Model matrix with intercept and two covariates
groups <- gl(5, 1)  # Random grouping factor

# Fit the model
model <- fit_glmm(y, X, groups)

# Check the results
print(model$par)  # Estimated coefficients