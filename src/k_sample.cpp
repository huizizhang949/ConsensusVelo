#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
int draw_cat(NumericVector prob){

  // number of categories
  int k = prob.size();

  // the function returns a one-hot encoding
  // set first argument to 1 for a Categorical distribution
  IntegerVector outcome(k);
  rmultinom(1, prob.begin(), k, outcome.begin());

  // return the label
  // we need to target the label '1' (others are all 0)


  // if use match function, the first element should be a vector,
  // and it returns type IntegerVector. No need to increment the index by 1.
  //IntegerVector target(1,1);
  //IntegerVector ix = match(target, Z);

  int ix = which_max(outcome) + 1;

  return ix;

}


// [[Rcpp::export]]
IntegerVector k_sample_cpp(NumericMatrix Prob){

  int C = Prob.nrow();

  IntegerVector k(C);

  for (int c = 0; c < C; ++c){
    NumericVector prob_c = Prob.row(c);
    k[c] = draw_cat(prob_c);
  }

  return k;

}


