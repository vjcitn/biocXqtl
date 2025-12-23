#' Fit multiple linear models with common design matrix
#'
#' @param X Design matrix (n x p)
#' @param Y Response matrix (n x m), each column is a separate response
#' @param intercept Logical, add intercept column to X?
#' @return List with coefficients and standard errors matrices
fastLmMany_R <- function(X, Y, intercept = TRUE) {
    # Convert to matrices if needed
    if (!is.matrix(X)) X <- as.matrix(X)
    if (!is.matrix(Y)) Y <- as.matrix(Y)
    
    # Store original X column names before modification
    orig_X_names <- colnames(X)
    orig_p <- ncol(X)
    
    # Add intercept if requested
    if (intercept) {
        X <- cbind(1, X)
    }
    
    # Check dimensions
    if (nrow(X) != nrow(Y)) {
        stop("X and Y must have the same number of rows")
    }
    
    # Call C++ function
    result <- fastLmMany(X, Y)
    
    # Create coefficient names
    if (intercept) {
        if (is.null(orig_X_names)) {
            coef_names <- c("(Intercept)", paste0("X", seq_len(orig_p)))
        } else {
            coef_names <- c("(Intercept)", orig_X_names)
        }
    } else {
        if (is.null(orig_X_names)) {
            coef_names <- paste0("X", seq_len(orig_p))
        } else {
            coef_names <- orig_X_names
        }
    }
    
    # Create response names
    if (is.null(colnames(Y))) {
        resp_names <- paste0("Y", seq_len(ncol(Y)))
    } else {
        resp_names <- colnames(Y)
    }
    
    # Add row and column names
    rownames(result$coefficients) <- coef_names
    rownames(result$se) <- coef_names
    colnames(result$coefficients) <- resp_names
    colnames(result$se) <- resp_names
    
    return(result)
}

