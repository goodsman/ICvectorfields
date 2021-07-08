#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MoransI(SEXP mat1, SEXP r1) {

  NumericMatrix Mat1(mat1);
  int R1 = as<int>(r1);
  int Rows = Mat1.nrow();
  int Cols = Mat1.ncol();
  int n1 = Rows * Cols;

  IntegerVector RowIndex(n1);
  IntegerVector ColIndex(n1);

  // containers for the distance matrix
  // and neighbour matrix in vector form
  NumericVector dvec(n1);
  NumericVector nvec(n1);

  // containers for products and for
  // the vectorized matrix
  NumericVector prod1(n1);
  NumericVector mvec(n1);

  // containers for differences in first
  // and second dimension.
  IntegerVector diff1(n1);
  IntegerVector diff2(n1);

  // containers for mean and sum of squares
  double mu;
  int ssq;

  // dummy variables for row and column indices
  int m;
  int n;

  // computing mean and sum of squared difference
  // from the mean and then mean centering the matrix
  mu = sum(Mat1)/n1;
  ssq = sum(pow(Mat1 - mu, 2));
  Mat1 = Mat1 - sum(Mat1)/n1;

  // computing row and column indices
  // and filling in mvec
  for(int i = 0; i < n1; ++i) {
    ColIndex[i] = ceil(i / Rows);
    RowIndex[i] = i - (ColIndex[i] - 1) * Rows;

    // vectorizing the matrix
    m = RowIndex[i] - 1;
    n = ColIndex[i] - 1;
    mvec[i] = Mat1[m, n];
  }

  double MoransIout = 0.0;
  double WtSum = 0.0;
  for(int j = 0; j < n1; ++j) {

    diff1 = RowIndex[j] - RowIndex;
    diff2 = ColIndex[j] - ColIndex;
    dvec = sqrt(pow(diff1, 2) + pow(diff2, 2));

    nvec[dvec <= R1] = 1.0;
    nvec[dvec == 0.0] = 0.0;
    nvec[dvec > R1] = 0.0;

    prod1 = mvec[j] * mvec;

    MoransIout += sum(prod1 * nvec);
    WtSum += sum(nvec);
  }

  // re-scaling
  if (ssq > 0) {
    MoransIout = MoransIout*n1/WtSum/ssq;
  } else {
    MoransIout = -999.0;
  }

  return MoransIout;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
MoransI(runif(9) nrow = 3), 1)
*/
