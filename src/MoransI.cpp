#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MoransI(SEXP mat1, SEXP r1) {

  NumericMatrix Mat1(mat1);
  int R1 = as<int>(r1);
  int Rows = Mat1.nrow();
  int Cols = Mat1.ncol();
  int n1 = Rows * Cols;
  int k;
  int lowerp;
  int lowerq;

  // containers for distance neighbour classifier
  double dist;
  double neigh;

  // containers for mean and sum of squares
  double mu;
  double ssq;

  // dummy variables for row and column indices
  int m;
  int n;

  // computing mean and sum of squared difference
  // from the mean and then mean centering the matrix
  mu = sum(Mat1)/n1;
  ssq = sum(pow(Mat1 - mu, 2));
  Mat1 = Mat1 - mu;

  double MoransIout = 0.0;
  double WtSum = 0.0;
  for(int j = 0; j < n1; ++j) {
    // creating the neighbourhood classifier
    // first figure out a smaller subset
    // of locations that are in the neighbourhood.
    n = floor(j / Rows);    // column index
    m = j - n * Rows;       // row index
    if (m - R1 < 0) {
      lowerp = 0;
    } else {
      lowerp = m - R1;
    }
    if (n - R1 < 0) {
      lowerq = 0;
    } else {
      lowerq = n - R1;
    }
    for(int p = lowerp; p < m + R1 + 1; ++p){
      for(int q = lowerq; q < n + R1 + 1; ++q){
        // computing the corresponding vector indices
        // again I preserve consistency with R's
        // column major approach.
        k = q*Rows + p;

        // calculating the distance matrix/vector
        if (k >= 0 && k < n1) {
          dist = sqrt(pow(m - p, 2) + pow(n - q, 2));
        }

        // filling in the neighbourhood matrix/vector
        if (k >= 0 && k < n1 && dist <= R1) {
          neigh = 1.0;
        }
        if (k >= 0 && k < n1 && dist == 0) {
          neigh = 0.0;
        }
        if (k >= 0 && k < n1 && dist > R1) {
          neigh = 0.0;
        }

        // adding to the cross-product
        if (k >= 0 && k < n1 && p < Rows && q < Cols) {
          MoransIout += Mat1(m, n) * Mat1(p, q) * neigh;
          WtSum += neigh;
        }
      }
    }
  }

  // re-scaling
  if (ssq > 0) {
    MoransIout = MoransIout*n1/WtSum/ssq;
  } else {
    MoransIout = -999.0;
  }

  return MoransIout;
}
