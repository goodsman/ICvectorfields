#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector MoransI(NumericMatrix mat1, IntegerVector rad1) {
  return mat1 * 2 * rad1;

  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
MoransI(matrix(rep(1, 9), nrow = 3), 3)
*/
