//  from claude dec 23
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

// [[Rcpp::export]]
Rcpp::List fastLmMany(const Eigen::Map<Eigen::MatrixXd> X,
                      const Eigen::Map<Eigen::MatrixXd> Y) {
    
    const int n = X.rows();
    const int p = X.cols();
    const int m = Y.cols();  // number of response variables
    
    // QR decomposition of X (only done once!)
    Eigen::ColPivHouseholderQR<MatrixXd> qr(X);
    const int rank = qr.rank();
    
    // Initialize output matrices
    MatrixXd coefs(p, m);
    MatrixXd se(p, m);
    
    // Solve for each response variable
    for (int i = 0; i < m; i++) {
        VectorXd y = Y.col(i);
        
        // Solve for coefficients
        VectorXd beta = qr.solve(y);
        coefs.col(i) = beta;
        
        // Calculate residuals
        VectorXd resid = y - X * beta;
        
        // Residual standard error
        double rss = resid.squaredNorm();
        int df = n - rank;
        double sigma2 = rss / df;
        
        // Standard errors: sqrt(diag((X'X)^-1) * sigma2)
        // Using QR: (X'X)^-1 = (R'R)^-1 = R^-1 * R^-T
        MatrixXd R = qr.matrixQR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
        MatrixXd Rinv = R.inverse();
        VectorXd XtX_inv_diag = (Rinv * Rinv.transpose()).diagonal();
        
        // Handle rank deficiency
        VectorXd se_vec = VectorXd::Constant(p, NA_REAL);
        Eigen::PermutationMatrix<Eigen::Dynamic> P = qr.colsPermutation();
        
        for (int j = 0; j < rank; j++) {
            int idx = P.indices()[j];
            se_vec[idx] = std::sqrt(XtX_inv_diag[j] * sigma2);
        }
        
        se.col(i) = se_vec;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("coefficients") = coefs,
        Rcpp::Named("se") = se,
        Rcpp::Named("rank") = rank
    );
}
