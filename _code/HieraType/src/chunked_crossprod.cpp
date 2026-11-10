#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

//' Chunked weighted cross-product for semi-dense sparse matrices
//'
//' Processes rows of X in dense chunks so each t(X)*X call hits BLAS dsyrk/dgemm,
//' avoiding CHOLMOD's sparse path that is slow at >10% density.
//' Peak memory is O(chunk_size * p) instead of O(N * p).
//'
//' @param X sparse N x p dgCMatrix (Xsigned assembled on the R side)
//' @param y_r dense N x m response matrix (double)
//' @param w_r prior_level_weights vector of length N, or length-0 numeric(0) for unweighted
//' @param chunk_size number of rows per BLAS chunk (default 10000)
//' @return named list: a = p x p (X'WX), b = p x m (X'Wy)
// [[Rcpp::export]]
Rcpp::List chunked_weighted_crossprod(
    const Eigen::MappedSparseMatrix<double>& X,
    Rcpp::NumericMatrix                      y_r,
    Rcpp::NumericVector                      w_r,
    int                                      chunk_size = 10000
) {
  typedef Eigen::MappedSparseMatrix<double>::InnerIterator SpIter;

  int N = X.rows(), p = X.cols(), m = y_r.ncol();
  bool weighted = ((int)w_r.size() == N);

  Eigen::Map<Eigen::MatrixXd> y(y_r.begin(), N, m);
  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(p, p);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p, m);

  for(int start = 0; start < N; start += chunk_size) {
    int end = std::min(start + chunk_size, N);
    int cnt = end - start;

    // Materialise sparse chunk as dense — from here t(chunk)*chunk calls BLAS
    Eigen::MatrixXd chunk = Eigen::MatrixXd::Zero(cnt, p);
    for(int col = 0; col < p; ++col) {
      for(SpIter it(X, col); it; ++it) {
        int r = it.row();
        if(r >= end) break;        // dgCMatrix rows are sorted within each column
        if(r >= start) chunk(r - start, col) = it.value();
      }
    }

    Eigen::MatrixXd yc = y.middleRows(start, cnt);

    if(weighted) {
      Eigen::VectorXd wc   = Eigen::Map<Eigen::VectorXd>(w_r.begin() + start, cnt);
      Eigen::VectorXd sqwc = wc.array().sqrt();
      // Row-scale: Xw(i,j) = chunk(i,j) * sqrt(w[i])
      Eigen::MatrixXd Xw = sqwc.asDiagonal() * chunk;
      a.noalias() += Xw.transpose() * Xw;
      b.noalias() += chunk.transpose() * (wc.asDiagonal() * yc);
    } else {
      a.noalias() += chunk.transpose() * chunk;
      b.noalias() += chunk.transpose() * yc;
    }
  }

  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("b") = b);
}
