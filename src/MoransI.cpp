#include <Rcpp.h>
using namespace Rcpp;
//' Efficiently compute Moran's I statistic
//'
//' Compute Moran's I for a matrix. A fast implementation of Moran's I for
//' gridded data, with neighbours defined based on a radial distance. Note
//' that when using radius to define the neighbourhood, a radius of one
//' corresponds to the rook's neibhourhood. There is currently no
//' equivalent to queen's neighbourhood.
//'
//' @param mat1 a matrix of values; NA/Inf values must be coded as NA and are ignored
//' @param r1 an integer representing the distance (radius), within which nearby
//'     cells are considered neighbours in units of rows/columns
//'
//' @return a single numeric value for Moran's I
//'
//' @export
//'
//' @examples
//' (TestMat <- matrix(c(1, 0, 1, 0, 1,
//'                      0, 1, 0, 1, 0,
//'                      1, 0, 1, 0, 1,
//'                      0, 1, 0, 1, 0,
//'                      1, 0, 1, 0, 1),
//'                    nrow =5))
//' # the code below should return -1
//' MoransI(TestMat, r1 = 1)
//[[Rcpp::export]]
double MoransI(SEXP mat1, SEXP r1) {
  NumericMatrix Mat1(mat1);
  int R1 = as<int>(r1);
  int Rows = Mat1.rows();
  int Cols = Mat1.cols();
  int n1 = Rows * Cols;
  int n2 = 0;
  int k;
  int lowerp;
  int lowerq;

  // containers for distance neighbour classifier
  double dist;
  double neigh = 0.0;

  // container for mean
  double mu = 0.0;

  // dummy variables for row and column indices
  int m;
  int n;

  // computing mean while dropping NA/Inf
  // values coded as NA
  for(int i = 0; i < n1; ++i) {
    n = floor(i / Rows);    // column index
    m = i - n * Rows;       // row index
    if (NumericVector::is_na(Mat1(m, n)) == 0) {
      mu += Mat1(m, n);
      n2 += 1;
    }
  }

  // normalizing the sum to obtain the mean
  mu = mu/n2;

  double MoransIout = 0.0;
  double WtSum = 0.0;
  double ssq = 0.0;
  for(int j = 0; j < n1; ++j) {
    Rcpp::checkUserInterrupt();
    n = floor(j / Rows);    // column index
    m = j - n * Rows;       // row index

    // computing the sum of squared differences from the mean
    if (NumericVector::is_na(Mat1(m, n)) == 0) {
      ssq += pow(Mat1(m, n) - mu, 2);
    }

    // creating the neighbourhood classifier
    // first figure out a smaller subset
    // of locations that are in the neighbourhood.
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
        if (k >= 0 && k < n1 && p < Rows && q < Cols && NumericVector::is_na(Mat1(m, n)) == 0 && NumericVector::is_na(Mat1(p, q)) == 0) {
          MoransIout += (Mat1(m, n) - mu) * (Mat1(p, q) - mu) * neigh;
          WtSum += neigh;
        }
      }
    }
  }

  // re-scaling
  if (ssq > 0) {
    // it is important to use n2 rather than n1 below, which is
    // the number of observations with NA/Inf values flagged
    // as -999.0 and removed.
    MoransIout = MoransIout/ssq/WtSum*n2;
  } else {
    MoransIout = -999.0;
  }

  return MoransIout;
}
