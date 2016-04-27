#include <RcppArmadillo.h>
using namespace arma;

// Weighted samplling with replacement. Input x is a vector of weight.
// [[Rcpp::export]]
arma::vec sample2 (arma::vec x, int N) {
  int n = x.n_elem;
  x = x/arma::sum(x);
  arma::vec p(n + 1);
  p(0) = 0;
  for (int i = 1; i < n + 1; i++) p(i) = p(i - 1) + x(i - 1);

  arma::vec runif = arma::randu<vec>(N);
  arma::vec out(N);
  for (int i = 0; i < N; i++) {
    int a = 0; int b = n;
    while (b - a > 1) {
      if  (runif(i) >= p((a + b)/2)) {
        a = (a + b)/2;
      } else {
        b = (a + b)/2;
      }
    }
    out(i) = a;
  }
  return out;
}

// Generate polynomial matrix: 1 x x^2 ......
// [[Rcpp::export]]
arma::mat poly2(arma::vec x, int r) {
  arma::mat out(x.n_elem, r + 1);
  out.col(0).fill(1.0);
  for (int i = 1; i < r + 1 ; i++) out.col(i) = arma::pow(x, i);
  return out;
}

// OLS and return the constant term: predicted value at cutoff point
// [[Rcpp::export]]
double lpreg1(arma::mat X, arma::mat W, arma::vec y) {
  arma::vec coef = arma::inv_sympd(X.t()*W*X)*X.t()*W*y;
  return coef(0);
}

// OLS and return fitted values + residuals
// [[Rcpp::export]]
arma::mat lpreg2(arma::mat X, arma::mat W, arma::vec y) {
  arma::mat out(y.n_elem, 2);
  arma::vec coef = arma::inv_sympd(X.t()*W*X)*X.t()*W*y;
  out.col(0) = X*coef;
  out.col(1) = y - out.col(0);
  return out;
}

// Sharp RD: estimate treatment effect tau
// [[Rcpp::export]]
double srd(arma::mat yxwl, arma::mat yxwr, int r) { // only wh

  arma::mat Xl = poly2(yxwl.col(1), r);
  arma::mat Wl = arma::diagmat(yxwl.col(2));

  arma::mat Xr = poly2(yxwr.col(1), r);
  arma::mat Wr = arma::diagmat(yxwr.col(2));

  double tau = lpreg1(Xr, Wr, yxwr.col(0)) - lpreg1(Xl, Wl, yxwl.col(0));
  return tau;
}


// Sharp RD: residual sampling
// [[Rcpp::export]]
arma::mat srdsampling(arma::vec yres, arma::vec w, int Nbc) {

  int n = yres.n_elem;
  int q = 0;
  arma::vec random = sample2(w, n*Nbc);

  arma::mat out(n, Nbc);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Nbc; j++) {
      out(i, j) = yres(random(q)); q++;
    }
  }

  return out;
}


// Sharp RD: estimate bias-corrected treatment effect tau
// [[Rcpp::export]]
arma::vec srdbc(arma::mat yxwl, arma::mat yxwr, int Nbc, int p, int q, int dist) { // wh and wb

  arma::uvec indexh(3), indexb(3), indexlh, indexrh;
  indexh << 0 << 1 << 2;
  indexb << 0 << 1 << 3;
  indexlh = arma::find(yxwl.col(2) > 0);
  indexrh = arma::find(yxwr.col(2) > 0);

  double ptau = srd(yxwl.submat(indexlh, indexh), yxwr.submat(indexrh, indexh), p);
  double qtau = srd(yxwl.cols(indexb), yxwr.cols(indexb), q);

  // left model of order q
  arma::mat Xl = poly2(yxwl.col(1), q);
  arma::mat Wl = arma::diagmat(yxwl.col(3));
  arma::mat yleft = lpreg2(Xl, Wl, yxwl.col(0));

  // right model of order q
  arma::mat Xr = poly2(yxwr.col(1), q);
  arma::mat Wr = arma::diagmat(yxwr.col(3));
  arma::mat yright = lpreg2(Xr, Wr, yxwr.col(0));

  // bootstrap new data
  arma::mat bootleft = srdsampling(yleft.col(1), yxwl.col(3), Nbc);
  arma::mat bootright = srdsampling(yright.col(1), yxwr.col(3), Nbc);

  // run regression with new data
  arma::vec Bptau(Nbc);
  for (int i = 0; i < Nbc; i++) {

    arma::mat newl = yxwl;
    newl.col(0) = yleft.col(0) + bootleft.col(i);

    arma::mat newr = yxwr;
    newr.col(0) = yright.col(0) + bootright.col(i);

    Bptau(i) = srd(newl.submat(indexlh, indexh), newr.submat(indexrh, indexh), p);
  }

  if (dist == 0) {
    arma::vec out(3);
    out(0) = ptau;
    out(1) = qtau;
    out(2) = ptau + qtau - mean(Bptau);
    return out;
  } else {
    arma::vec out(3 + Nbc);
    out(0) = ptau;
    out(1) = qtau;
    out(2) = ptau + qtau - mean(Bptau);
    out.tail(Nbc) = Bptau;
    return out;
  }
}


// Sharp RD: bootstrap bias-corrected treatment effect tau
// [[Rcpp::export]]
arma::vec srdbcboot(arma::mat yxwl, arma::mat yxwr, int Nbc, int Nci, int p, int q) {

  arma::uvec indexh(3), indexb(3);
  indexh << 0 << 1 << 2;
  indexb << 0 << 1 << 3;

  // left model of order q
  arma::mat Xl = poly2(yxwl.col(1), q);
  arma::mat Wl = arma::diagmat(yxwl.col(3));
  arma::mat yleft = lpreg2(Xl, Wl, yxwl.col(0));

  // right model of order q
  arma::mat Xr = poly2(yxwr.col(1), q);
  arma::mat Wr = arma::diagmat(yxwr.col(3));
  arma::mat yright = lpreg2(Xr, Wr, yxwr.col(0));

  // bootstrap new data
  arma::mat bootleft = srdsampling(yleft.col(1), yxwl.col(3), Nci);
  arma::mat bootright = srdsampling(yright.col(1), yxwr.col(3), Nci);

  // run regression with new data
  arma::vec Btaubc(Nci);
  for (int i = 0; i < Nci; i++) {

    arma::mat newl = yxwl;
    newl.col(0) = yleft.col(0) + bootleft.col(i);

    arma::mat newr = yxwr;
    newr.col(0) = yright.col(0) + bootright.col(i);

    Btaubc(i) = srdbc(newl, newr, Nbc, p, q, 0)(2);
  }

  return Btaubc;
}

