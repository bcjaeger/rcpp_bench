#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double log_rank_test_wtd(arma::mat& y_node,
                         arma::vec& w_node,
                         arma::vec& group){

 double temp1, temp2, deaths;
 double n_risk=0, observed=0, expected=0, V=0, g_risk=0;
 bool break_loop = false;

 arma::uword j = y_node.n_rows-1;

 while (!break_loop){

  temp1 = y_node.at(j, 0);

  deaths = 0;

  for (; y_node.at(j, 0) == temp1; j--) {

   n_risk += w_node[j];
   deaths += y_node.at(j, 1) * w_node[j];
   g_risk += group[j] * w_node[j];
   observed += y_node.at(j, 1) * group[j] * w_node[j];

   if(j == 0){
    break_loop = true;
    break;
   }

  }

  // should only do these calculations if deaths > 0,
  // but turns out its faster to multiply by 0 than
  // it is to check whether deaths is > 0

  temp2 = g_risk / n_risk;
  expected += deaths * temp2;

  // update variance if n_risk > 1 (if n_risk == 1, variance is 0)
  // definitely check if n_risk is > 1 b/c otherwise divide by 0
  if (n_risk>1){
   temp1 = deaths * temp2 * (n_risk-deaths) / (n_risk-1);
   V += temp1 * (1 - temp2);
  }

 }


 return(pow(expected-observed, 2) / V);


}
