#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat lrt_multi_v1(arma::mat& y_node,
                       arma::vec& XB,
                       double min_obs,
                       double min_events){

 double temp1, temp2, deaths, n_risk, observed, expected, V, g_risk;
 double n_events=0, n_obs=0;
 bool break_loop;

 arma::uvec XB_sort = arma::sort_index(XB);

 arma::uword i, i_lower;

 for(i = 0; i < y_node.n_rows; i++){
  n_events += y_node.at(XB_sort[i], 1);
  n_obs++;
  if(n_events >= min_events && n_obs >= min_obs) break;
 }

 i_lower = i;

 Rcout << "i_lower: " << i_lower << std::endl;

 n_events=0;
 n_obs=0;
 arma::vec group(y_node.n_rows);
 group.fill(0);

 for(i = y_node.n_rows - 1; ; i--){
  n_events += y_node.at(XB_sort[i], 1);
  group[XB_sort[i]] = 1;
  n_obs++;
  if(n_events >= min_events && n_obs >= min_obs) break;
 }


 arma::uword total_steps = i - i_lower;
 arma::uword g_index = i;
 arma::uword t;
 Rcout << "g_index: " << g_index << std::endl;

 arma::mat out(total_steps, 2);

 for(t = 0; t < total_steps-1; t++){

  n_risk=0;
  g_risk=0;

  observed=0;
  expected=0;

  V=0;

  break_loop = false;

  i = y_node.n_rows-1;

  while (!break_loop){

   temp1 = y_node.at(i, 0);

   deaths = 0;

   for (; y_node.at(i, 0) == temp1; i--) {

    n_risk++;
    deaths += y_node.at(i, 1);
    g_risk += group[i];
    observed += y_node.at(i, 1) * group[i];

    if(i == 0){
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

  Rcout << "V: " << V << std::endl;

  out(t,0) = XB[XB_sort(g_index)];
  out(t,1) = pow(expected-observed, 2) / V;
  Rcout << "LRT chisq:  " << out(t,1) << std::endl;

  Rcout << "g_index: " << g_index << std::endl;

  group(XB_sort[g_index]) = 1;
  g_index--;

  Rcout << "donezo" << std::endl;

 }

 return(out);

}
