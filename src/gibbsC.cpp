#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @param n the parameter of the conditional distribution(binomal) of x
//' @param a the first parameter of the conditional distribution(beta) of y
//' @param b the second parameter of the conditional distribution(beta) of y
//' @return generated samples \code{mat}
//' @examples
//' \dontrun{
//'	a <- b <- 5;n <- 20;N <- 1000
//' gibR <- gibbsR(N, 10, n, a, b)
//' par(mfrow=c(2,1));
//' plot(gibR[,1],type='l')
//' plot(gibR[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin, int n, int a, int b) {
  NumericMatrix mat(N, 2);
  double x = 10, y = 0.5; 
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, x + a, n - x + b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
