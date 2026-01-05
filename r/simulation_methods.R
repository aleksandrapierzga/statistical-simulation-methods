# Task 1
# Acceptance–rejection method
generate_normal_rejection <- function(mu = 0, sigma = 1, n = 1) {
  x <- numeric(n)
  i <- 1
  while (i <= n) {
    u <- runif(1, min = -pi/2, max = pi/2)  # Generate U from a uniform distribution
    v <- runif(1)  # Generate V from a uniform distribution
    x_proposed <- tan(u)  # Transform to a proposed value from the Cauchy distribution
    density_proposed <- dnorm(x_proposed, mean = mu, sd = sigma)
    if (v <= density_proposed) {
      x[i] <- x_proposed
      i <- i + 1
    }
  }
  return(mu + sigma * x)
}

# Example usage
set.seed(123)  # Set seed for reproducibility

# Measure execution time
system.time({
  sample_rejection <- generate_normal_rejection(mu = 0, sigma = 1, n = n)
})

# Display the first few values
head(sample_rejection)

# Test sample normality
test_rejection <- shapiro.test(sample_rejection)
list(test_rejection_result = test_rejection)

# Inverse CDF method
generate_normal_inverse_cdf <- function(mu = 0, sigma = 1, n = 1) {
  # Generate n random numbers from a uniform distribution
  u <- runif(n)

  # Transform to the normal distribution
  x <- qnorm(u, mean = mu, sd = sigma)

  return(x)
}

# Example usage
set.seed(123)  # Set seed for reproducibility

# Measure execution time
system.time({
  sample_inverse_cdf <- generate_normal_inverse_cdf(mu = 0, sigma = 1, n = 1000)
})

# Display the first few values
head(sample_inverse_cdf)

# Test sample normality
test_inverse_cdf <- shapiro.test(sample_inverse_cdf)
list(test_inverse_cdf_result = test_inverse_cdf)

# Transformation method (built-in generator)
generate_normal_transform <- function(n, mu = 0, sigma = 1) {
  rnorm(n, mean = mu, sd = sigma)
}

# Example usage
set.seed(123)
n <- 1000

# Measure execution time
system.time({
  sample_transform <- generate_normal_transform(n, mu = 0, sigma = 1)
})

# Display the first few samples
head(sample_transform)

# Test sample normality
test_transform <- shapiro.test(sample_transform)
list(test_transform_result = test_transform)

# Box–Muller method
generate_normal_box_muller <- function(n, mu = 0, sigma = 1) {
  m <- ceiling(n / 2)
  u1 <- runif(m)  # Generate m random numbers U1
  u2 <- runif(m)  # Generate m random numbers U2
  z1 <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)  # First Z variable
  z2 <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)  # Second Z variable
  z <- c(z1, z2)[1:n]  # Trim to desired length n
  return(mu + sigma * z)  # Scale to desired mu and sigma
}

# Example usage
set.seed(123)
n <- 1000

# Measure execution time and generate sample
system.time({
  sample_box_muller <- generate_normal_box_muller(n, mu = 0, sigma = 1)
})

# Display the first few values
head(sample_box_muller)

# Test sample normality
test_box_muller <- shapiro.test(sample_box_muller)
list(test_box_muller_result = test_box_muller)

# Task 2

library(MASS)
library(Matrix)

# Generate a matrix from the Wishart distribution
generate_wishart <- function(d, n, sigma) {
  # Check whether sigma is symmetric and positive definite
  if (!(Matrix::isSymmetric(sigma) && all(eigen(sigma)$values > 0))) {
    stop("Sigma must be symmetric and positive definite.")
  }

  # Generate Z matrix of size d x n from the normal distribution
  Z <- matrix(rnorm(d * n), nrow = d, ncol = n)

  # Cholesky decomposition of sigma to obtain matrix A
  A <- chol(sigma)

  # Scale Z by A to obtain AZ
  AZ <- A %*% Z

  # Compute M as AZ multiplied by the transpose of AZ
  M <- AZ %*% t(AZ)

  return(M)
}

# Define parameters for the Wishart distribution
d <- 3
n <- 5
sigma <- diag(1, d)

# Generate a single Wishart sample
wishart_sample <- generate_wishart(d, n, sigma)

# Print the generated sample
print(wishart_sample)

# Verify correctness by generating many samples
samples <- lapply(1:10000, function(i) generate_wishart(d, n, sigma))

# Compute empirical mean and covariance
empirical_mean <- Reduce("+", samples) / length(samples)
empirical_cov <- matrix(0, nrow = d^2, ncol = d^2)
for (s in samples) {
  vec_s <- as.vector(s)
  empirical_cov <- empirical_cov + (vec_s %*% t(vec_s))
}
empirical_cov <- empirical_cov / length(samples)

theoretical_mean <- n * sigma
theoretical_cov <- diag(2 * n, d^2)

# Print results for comparison
print(paste("Empirical mean of generated samples:", empirical_mean))
print(paste("Theoretical mean of the Wishart distribution:", theoretical_mean))

# Task 3

# Function for Monte Carlo method to estimate pi
monte_carlo_pi <- function(method_function, n) {
  pi_estimate <- method_function(n)
  return(pi_estimate)
}

# Example Monte Carlo methods for pi
dartboard_method <- function(n) {
  x <- runif(n, -1, 1)
  y <- runif(n, -1, 1)
  inside_circle <- sum(x^2 + y^2 <= 1)
  return(4 * inside_circle / n)
}

integration_method <- function(n) {
  x <- runif(n)
  y <- runif(n)
  inside_circle <- sum(x^2 + y^2 <= 1)
  return(4 * inside_circle / n)
}

# Set the number of iterations
n <- 1000000

print(monte_carlo_pi(dartboard_method, n))
print(monte_carlo_pi(integration_method, n))

# Check execution time for Dartboard Method
system.time(pi_dartboard <- monte_carlo_pi(dartboard_method, n))

# Check execution time for Integration Method
system.time(pi_integration <- monte_carlo_pi(integration_method, n))

# True value of pi
true_pi <- pi

# Calculate estimates for both methods
pi_dartboard <- monte_carlo_pi(dartboard_method, n)
pi_integration <- monte_carlo_pi(integration_method, n)

# Calculate absolute differences from true pi
abs_diff_dartboard <- abs(pi_dartboard - true_pi)
abs_diff_integration <- abs(pi_integration - true_pi)

# Print absolute differences
cat(abs_diff_dartboard)
cat(abs_diff_integration)
