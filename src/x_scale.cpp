
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::uword i, j, n_vars;

arma::mat
  x_input,
  y_input;

arma::uvec
  cols,
  rows;

arma::uvec::iterator iit, jit;

arma::vec w_input;

// [[Rcpp::export]]
arma::mat x_scale_wtd(NumericMatrix&  x,
                      NumericVector&  w,
                      arma::ivec& cols_,
                      arma::ivec& rows_){


  x_input = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
  w_input = arma::vec(w.begin(), w.size(), false);
  cols = arma::conv_to<arma::uvec>::from(cols_);
  rows = arma::conv_to<arma::uvec>::from(rows_);

  n_vars = cols.size();

  arma::mat x_node = x_input(rows, cols);
  arma::vec w_node = w_input(rows);

  arma::mat out(n_vars, 2);
  arma::vec means = out.unsafe_col(0);   // Reference to column 1
  arma::vec scales = out.unsafe_col(1);  // Reference to column 2

  double w_node_sum = arma::sum(w_node);

  for(i = 0; i < n_vars; i++) {

    arma::vec x_i = x_node.unsafe_col(i);

    means.at(i) = arma::sum( w_node % x_i ) / w_node_sum;

    x_i -= means.at(i);

    scales.at(i) = arma::sum(w_node % arma::abs(x_i));

    if(scales(i) > 0)
      scales.at(i) = w_node_sum / scales.at(i);
    else
      scales.at(i) = 1.0; // rare case of constant covariate;

    x_i *= scales.at(i);

  }


  return(out);

}

// [[Rcpp::export]]
arma::mat x_scale_wtd_iter(NumericMatrix&  x,
                           NumericVector&  w,
                           arma::ivec& cols_,
                           arma::ivec& rows_){


  x_input = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
  w_input = arma::vec(w.begin(), w.size(), false);
  cols = arma::conv_to<arma::uvec>::from(cols_);
  rows = arma::conv_to<arma::uvec>::from(rows_);

  arma::mat out(n_vars, 2);
  arma::vec means = out.unsafe_col(0);   // Reference to column 1
  arma::vec scales = out.unsafe_col(1);  // Reference to column 2

  double w_node_sum = arma::sum(w_input(rows));

  i=0;

  for(iit = cols.begin(); iit != cols.end(); ++iit) {

    arma::vec x_i = x_input.unsafe_col(*iit);

    means.at(i) = arma::sum( w_input(rows) % x_i(rows) ) / w_node_sum;

    x_i(rows) -= means.at(i);

    scales.at(i) = arma::sum(w_input(rows) % arma::abs(x_i(rows)));

    if(scales(i) > 0)
      scales.at(i) = w_node_sum / scales.at(i);
    else
      scales.at(i) = 1.0; // rare case of constant covariate;

    x_i(rows) *= scales.at(i);

    i++;

  }


  return(out);

}

//
// // this is the version currently deployed.
// // It is very readable and very fast
//
// // [[Rcpp::export]]
// arma::mat x_scale_wtd_1(arma::mat& x_node,
//                         arma::vec& w_node){
//
//  // set aside memory for outputs
//  // first column holds the mean values
//  // second column holds the scale values
//
//  arma::uword n_vars = x_node.n_cols;
//
//  arma::mat out(n_vars, 2);
//  arma::vec means = out.unsafe_col(0);   // Reference to column 1
//  arma::vec scales = out.unsafe_col(1);  // Reference to column 2
//
//  double w_node_sum = arma::sum(w_node);
//
//  for(i = 0; i < n_vars; i++) {
//
//   arma::vec x_i = x_node.unsafe_col(i);
//
//   means.at(i) = arma::sum( w_node % x_i ) / w_node_sum;
//
//   x_i -= means.at(i);
//
//   scales.at(i) = arma::sum(w_node % arma::abs(x_i));
//
//   if(scales(i) > 0)
//    scales.at(i) = w_node_sum / scales.at(i);
//   else
//    scales.at(i) = 1.0; // rare case of constant covariate;
//
//   x_i *= scales.at(i);
//
//  }
//
//
//  return(out);
//
// }
//
// // this version is less readable and not any faster
//
// // [[Rcpp::export]]
// arma::mat x_scale_wtd_2(arma::mat& x_node,
//                         arma::vec& w_node){
//
//  // set aside memory for outputs
//  // first column holds the mean values
//  // second column holds the scale values
//
//  arma::uword n_vars = x_node.n_cols;
//
//  arma::mat out(n_vars, 2);
//  arma::vec means = out.unsafe_col(0);   // Reference to column 1
//  arma::vec scales = out.unsafe_col(1);  // Reference to column 2
//
//  double w_node_sum = arma::sum(w_node);
//
//  for(i = 0; i < n_vars; i++) {
//
//   out.at(i, 0) = arma::sum( w_node % x_node.col(i) ) / w_node_sum;
//
//   x_node.col(i) -= out.at(i, 0);
//
//   out.at(i,1) = arma::sum(w_node % arma::abs(x_node.col(i)));
//
//   if(scales(i) > 0)
//    out.at(i,1) = w_node_sum / out.at(i,1);
//   else
//    out.at(i,1) = 1.0; // rare case of constant covariate;
//
//   x_node.col(i) *= out.at(i,1);
//
//  }
//
//
//  return(out);
//
// }
