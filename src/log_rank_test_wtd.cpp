#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double log_rank_test_wtd(arma::mat& y_node,
                         arma::vec& w_node,
                         arma::vec& group){

 double temp1, temp2, n_events;
 double n_at_risk=0, observed=0, expected=0, V=0, g_risk=0;
 bool break_loop = false;

 arma::uword j = y_node.n_rows-1;

 while (!break_loop){

  temp1 = y_node.at(j, 0);

  n_events = 0;

  for (; y_node.at(j, 0) == temp1; j--) {

   n_at_risk += w_node[j];
   n_events += y_node.at(j, 1) * w_node[j];
   g_risk += group[j] * w_node[j];
   observed += y_node.at(j, 1) * group[j] * w_node[j];

   if(j == 0){
    break_loop = true;
    break;
   }

  }

  // should only do these calculations if n_events > 0,
  // but turns out its faster to multiply by 0 than
  // it is to check whether n_events is > 0

  temp2 = g_risk / n_at_risk;
  expected += n_events * temp2;

  // update variance if n_at_risk > 1 (if n_at_risk == 1, variance is 0)
  // definitely check if n_at_risk is > 1 b/c otherwise divide by 0
  if (n_at_risk>1){
   temp1 = n_events * temp2 * (n_at_risk-n_events) / (n_at_risk-1);
   V += temp1 * (1 - temp2);
  }

 }


 return(pow(expected-observed, 2) / V);


}
