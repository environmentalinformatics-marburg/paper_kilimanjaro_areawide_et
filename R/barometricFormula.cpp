#include <Rcpp.h>
using namespace Rcpp;

//// difference -----
// [[Rcpp::export]]
NumericVector difference(double x, NumericVector y) {
  
  int len = y.size();
  NumericVector out(len);
  for (int i = 0; i < len; i++) {
    out[i] = x - y[i];
  }
  
  return abs(out);
}

//// is.na() -----
// [[Rcpp::export]]
LogicalVector isNA(NumericVector x) {
  int n = x.size();
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
  }
  return out;
}

//// which.min() -----
// [[Rcpp::export]]
NumericVector whichMin(NumericVector x) {
  
  int len = x.size();
  bool isna, out = FALSE;

  int n = 0;
  NumericVector tmp;
  for (int i = 0; i < len; i++) {
    isna = isNA(x[i]);
    if (!isna) {
      tmp[n] = x[i];
      n += 1;
    }
  }
  
  // int i = 0;
  // while (!out | i == 0) {
  //   std::cout << i << "\n";
  //   isna = isNA(x[i]);
  //   if (!isna)
  //     out = (min(x) == x[i]);
  //   i++;
  // }
  
  return tmp;
}

//// barometric formula -----
// [[Rcpp::export]]
int barometricFormula(double z, NumericVector gp, NumericVector ta, 
                         IntegerVector p) {
  
  NumericVector differences = difference(z, gp);
  int id = whichMin(differences);

  // double h0, h1, t, p0;
  // 
  // bool isna = isNA(gp[id + 1]);
  // if ((differences[id - 1] < differences[id + 1]) | isna) {
  //   h0 = gp[id];
  //   h1 = z;
  //   t = ta[id];
  //   p0 = p[id];
  // } else {
  //   h0 = gp[id + 1];
  //   h1 = z;
  //   t = ta[id + 1];
  //   p0 = p[id + 1];
  // }
  // 
  // // apply barometric formula, directly taken from 
  // // https://en.wikipedia.org/wiki/Barometric_formula
  // double dh, g, M, R;
  // 
  // dh = h1 - h0;
  // g = 9.80665;
  // M = 0.0289644;
  // R = 8.31432;
  // 
  // double out;
  // out = p0 * exp((-1) * (g * M * dh) / (R * (t + 273.15)));
  
  return id;
}
