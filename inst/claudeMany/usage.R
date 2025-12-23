library(Rcpp)
library(RcppEigen)

# Compile the C++ code
sourceCpp("fastLmMany.cpp")
source("manyreg2.R")

# Generate example data
set.seed(123)
n <- 100
p <- 5
m <- 20  # 20 different response variables

X <- matrix(rnorm(n * p), n, p)
Y <- matrix(rnorm(n * m), n, m)

# Add some real effects
true_coefs <- matrix(rnorm(p * m), p, m)
Y <- X %*% true_coefs + Y

# Fit all models at once
result <- fastLmMany_R(X, Y, intercept = TRUE)

# View results
dim(result$coefficients)  # (p+1) x m
dim(result$se)            # (p+1) x m

# Compare with lm() for first response
lm_result <- lm(Y[,1] ~ X)
cbind(
    fastLm = result$coefficients[,1],
    lm = coef(lm_result),
    se_fastLm = result$se[,1],
    se_lm = summary(lm_result)$coefficients[,2]
)

