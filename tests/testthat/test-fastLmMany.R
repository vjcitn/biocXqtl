test_that("fastLmMany_R produces correct results", {
    set.seed(123)
    n <- 50
    p <- 3
    m <- 5
    
    X <- matrix(rnorm(n * p), n, p)
    Y <- matrix(rnorm(n * m), n, m)
    
    result <- fastLmMany_R(X, Y)
    
    # Compare with lm for first response
    lm_fit <- lm(Y[,1] ~ X)
    
    expect_equal(result$coefficients[,1], coef(lm_fit), tolerance = 1e-10)
    expect_equal(result$se[,1], summary(lm_fit)$coefficients[,2], 
                 tolerance = 1e-10)
})
