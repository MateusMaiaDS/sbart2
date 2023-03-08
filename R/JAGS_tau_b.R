# Building a JAGS model for splines with the penalized
rm(list=ls())
source("R/other_functions.R")
source("R/wrap_bart.R")
source("R/bayesian_simulation.R")
library(R2jags)
library(MASS) # Useful for mvrnorm function
library(splines) # Useful for creating the B-spline basis functions
# Maths -------------------------------------------------------------------

# Notation:
# y(t): Response variable at time t, defined on continuous time
# y: vector of all observations
# B: design matrix of spline basis functions
# beta; spline weights
# sigma: residual standard deviation parameter (sometimes known in the GP world as the nugget)
# sigma_b: spline random walk parameter

# Likelihood:
# y ~ N(B%*%beta, sigma)
# beta_j - beta_{j-1} ~ N(0, sigma_b)

# Priors
# sigma ~ cauchy(0, 10)
# sigma_b ~ cauchy(0, 10)

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
set.seed(42)
N <- 100 # Number of observations
x <- sort(runif(N, 0, 10)) # Create some covariate values
nIknots <- 100
knots <- quantile(x,seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]
B <- cbind(1,splines::ns(x = x,knots = knots,intercept = FALSE))
tau_b <- 1 # Parameters as above
tau <- 10
beta <- cumsum(c(1, rnorm(ncol(B) - 1, 0, sqrt(tau_b^-1))))
y <- rnorm(N, mean = B %*% beta, sd = sqrt(tau^-1))
plot(x, y)
lines(x, B %*% beta, col = "red") # True line

# Scaling y
min_y <- min(y)
max_y <- max(y)

y_scale <- normalize_bart(y = y,a = min_y,b = max_y)

# Setting other parameters
nsigma <- naive_sigma(x = x,y = y_scale)
df <- 3
# Calculating tau hyperparam
a_tau <- df/2
sigquant <- 0.9

tau_b_0 <- 4*(2^2)*1
# Calculating lambda
qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
lambda <- (nsigma*nsigma*qchi)/df
d_tau <- (lambda*df)/2


# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code <- "
model
{
  # Likelihood
  for (t in 1:N) {
    y[t] ~ dnorm(inprod(B[t,], beta), tau)
  }

  # RW prior on beta
  beta[1] ~ dnorm(0, tau_b_0)
  for (i in 2:N_knots) {
    beta[i] ~ dnorm(beta[i-1], tau_b)
  }

  # Priors on beta values
  tau ~ dgamma(a_tau, d_tau)
  tau_b ~ dgamma(0.5*nu,0.5*delta*nu)
  delta ~ dgamma(a_delta,d_delta)

}
"

# Set up the data
model_data <- list(N = N, y = y, B = B,
                   N_knots = ncol(B),
                   a_tau = a_tau,
                   d_tau = d_tau,
                   a_delta = 0.01,
                   d_delta = 0.01,
                   tau_b_0 = tau_b_0,
                   nu  = 20)

# Choose the parameters to watch
model_parameters <- c("beta", "tau", "tau_b","delta")

# Run the model - can be slow
model_run <- jags(
     data = model_data,
     parameters.to.save = model_parameters,
     model.file = textConnection(model_code)
)

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run)

# Get the posterior betas and 50% CI
beta_post <- model_run$BUGSoutput$sims.list$beta
beta_quantile <- apply(beta_post, 2, quantile, prob = c(0.25, 0.5, 0.75))

# Plot the output with uncertainty bands
plot(x, y)
lines(x, B %*% beta, col = "red") # True line
lines(x, B %*% beta_quantile[2, ], col = "blue") # Predicted line
lines(x, B %*% beta_quantile[1, ], col = "blue", lty = 2) # Predicted low
lines(x, B %*% beta_quantile[3, ], col = "blue", lty = 2) # Predicted high
legend("topleft", c(
     "True line",
     "Posterior lines (with 50% CI)",
     "Data"
),
lty = c(1, 1, -1),
pch = c(-1, -1, 1),
col = c("red", "blue", "black")
)

# Create some new predictions on a grid of new values
# Needs to be in the same range as the previous values (if not you need to go back to the creation of B above)
# x_new <- seq(min(x), max(x), length = 1000)
# nots <- quantile(x,seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]
# B_new <- cbind(1,splines::ns(x = x_new,knots = knots,intercept = FALSE))
# plot(x, y)
# lines(x_new, B_new %*% beta_quantile[2, ], col = "blue") # Beautifully smooth

# Real example ------------------------------------------------------------

# Data wrangling and jags code to run the model on a real data set in the data directory


# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks
