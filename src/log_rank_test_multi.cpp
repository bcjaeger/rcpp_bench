#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double lrt_multi_v1(arma::mat& y_node,
                    arma::vec& XB,
                    arma::vec& group,
                    arma::uword n_cps,
                    double min_obs,
                    double min_events){

  // about this function - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // this function returns a cutpoint obtaining the local maximum
  // of the log-rank test (lrt) statistic. The default value (+Inf)
  // is really for diagnostic purposes. Put another way, if the
  // return value is +Inf (an impossible value for a cutpoint),
  // that means that we didn't find any valid cut-points and
  // the node cannot be grown with the current XB.
  //
  // if there is a valid cut-point, then the main side effect
  // of this function is to modify the group vector, which
  // will be used to assign observations to the two new nodes.
  //
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  double temp1, temp2, observed, expected, V, g_risk;
  double deaths=0, n_risk=0, cp = R_PosInf;
  bool break_loop;

  // initialize at the lowest possible LRT stat value
  double lrt_current, lrt_max = 0;

  // sort XB- we need to iterate over the sorted indices
  arma::uvec XB_sort = arma::sort_index(XB);

  // speaking of which, here is that iterator.
  arma::uvec::iterator it, it_best;

  // unsafe columns point to specific cols in y_node.
  // this makes the code more readable and doesn't copy data
  arma::vec status = y_node.unsafe_col(1);
  arma::vec time = y_node.unsafe_col(0);

  // unsigned integers used by loops
  arma::uword
    // i is used within the lrt (inner) loop.
    i,
    // j is used in the outer loop.
    j,
    // k and p are used to determine the range of j.
    k=0, p=0;

  // first determine the lowest value of XB that will
  // be a valid cut-point to split a node. A valid cut-point
  // is one that, if used, will result in at least min_obs
  // and min_events in both the left and right node.
  for(it = XB_sort.begin(); it != XB_sort.end(); ++it){
    deaths += status[*it];
    n_risk++; // keep separate from k - n_risk will be weighted
    k++;
    if(deaths >= min_events && n_risk >= min_obs) break;
  }

  // got to reset these before finding the upper limit of j
  deaths=0; n_risk=0;

  // group should be initialized as all 0s
  group.fill(0);

  for(it = XB_sort.end()-1; it >= XB_sort.begin(); --it){
    deaths += status[*it];
    group[*it] = 1;
    n_risk++;
    p++;
    if(deaths >= min_events && n_risk >= min_obs) break;
  }

  // this is just done to avoid compilation warnings
  // (it_best should be initialized before it is used)
  it_best = it;

  // what happens if we don't have enough events or obs to split?
  // the first valid lower cut-point (at XB_sort[k]) is > the first
  // valid upper cutpoint (current value of n_risk). Put another way,
  // k (the number of steps taken from beginning of the XB vec)
  // will be > n_rows - p, where the difference on the RHS is
  // telling us where we are after taking p steps from the end
  // of the XB vec. Returning the infinite cp is a red flag.

  // re-form p to be our location in XB when we take p steps
  // away from the end of XB.
  p = y_node.n_rows - p;

  if ( k > p ) return(cp);

  // reset p to tell us how many total steps we need to take
  // in the j-loop that we are about to start.
  p = p - k;

  // reset k to tell us when we are at a landmark value of j
  k = ceil(p / n_cps);

  // Rcout << "upper index: " << y_node.n_rows - n_risk - 1 << std::endl;
  // Rcout << "lower index: " << steps_lower << std::endl;
  // Rcout << "total steps: " << total_steps << std::endl;
  // Rcout << "sum(g==1) = " << arma::sum(group) << std::endl;

  Rcout << "k: " << k << std::endl;
  Rcout << "p: " << p << std::endl;

  // begin outer loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(j = 0; j < p; j++){

    // begin inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - -
    // if j modulo q = 0, it means we are q steps into the next
    // batch of j values, which means its time to compute lrt.
    if(j % k == 0){

      n_risk=0;
      g_risk=0;

      observed=0;
      expected=0;

      V=0;

      break_loop = false;

      i = y_node.n_rows-1;

      while (!break_loop){

        temp1 = time[i];

        deaths = 0;

        for ( ; time[i] == temp1; i--) {

          n_risk++;
          deaths += status[i];
          g_risk += group[i];
          observed += status[i] * group[i];

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

      lrt_current = pow(expected-observed, 2) / V;

      Rcout << "lrt_current: " << lrt_current << std::endl;

      if(lrt_current > lrt_max){
        it_best = it;
        lrt_max = lrt_current;
        cp = XB[*it];
        Rcout << "cp" << cp << std::endl;
      }

    }
    // end inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

    // drop down one spot on XB
    --it;
    // update group to be consistent with where we are on XB
    group[*it] = 1;
    // end outer loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

  }

  // rewind it until it is back where it was when we got the
  // best lrt stat. While rewinding it, also reset the group
  // values so that group is as it was when we got the best
  // lrt stat.
  while(it <= it_best){
    group[*it] = 0;
    ++it;
  }

  return(cp);

}

// [[Rcpp::export]]
arma::mat lrt_multi_v2(arma::mat& y_node,
                       arma::vec& XB,
                       double min_obs,
                       double min_events){

  double temp1, deaths, observed, expected, V;
  double n_events=0, n_obs=0, n_risk=0, g_risk=0;
  bool break_loop=false;

  arma::uvec XB_sort = arma::sort_index(XB);

  arma::uword i, j, i_lower;

  for(i=0 ; i < y_node.n_rows; i++){
    n_events += y_node.at(XB_sort[i], 1);
    n_obs++;
    if(n_events >= min_events && n_obs >= min_obs) break;
  }

  i_lower = i;

  // Rcout << "i_lower: " << i_lower << std::endl;

  n_events=0; n_obs=0;
  arma::vec group(y_node.n_rows);
  group.fill(0);

  for(i = y_node.n_rows - 1; ; i--){
    n_events += y_node.at(XB_sort[i], 1);
    group[XB_sort[i]] = 1;
    n_obs++;
    if(n_events >= min_events && n_obs >= min_obs) break;
  }

  arma::uword g = i;
  arma::uword total_steps = g - i_lower;

  i = y_node.n_rows-1;
  j = 0;

  arma::mat out(total_steps, 2);
  arma::mat lrt(y_node.n_rows, 4);

  while (!break_loop){

    temp1 = y_node.at(i, 0);

    deaths = 0;

    for (; y_node.at(i, 0) == temp1; i--) {

      n_risk++;
      deaths   += y_node.at(i, 1);
      g_risk   += group[i];

      if(i == 0){
        break_loop = true;
        break;
      }

    }

    lrt(j, 0) = i+1;
    lrt(j, 1) = g_risk;
    lrt(j, 2) = n_risk;
    lrt(j, 3) = deaths;
    j++;

  }

  lrt.resize(j,4);

  //Rcout << "lrt first row: " << lrt.row(0) << std::endl;

  arma::vec lrt_index  = lrt.unsafe_col(0);
  arma::vec lrt_g_risk = lrt.unsafe_col(1);
  arma::vec lrt_n_risk = lrt.unsafe_col(2);
  arma::vec lrt_deaths = lrt.unsafe_col(3);
  arma::vec status     = y_node.unsafe_col(1);

  arma::vec v2 = lrt_g_risk / lrt_n_risk;
  arma::vec v1 = lrt_deaths % v2 % (lrt_n_risk-lrt_deaths) / (lrt_n_risk-1);
  arma::vec lrt_V = v1 % (1 - v2);

  expected = arma::sum(lrt_deaths % v2);
  observed = arma::sum(status % group);

  if(isnan(lrt_V[0])) lrt_V[0] = 0;

  V = arma::sum(lrt_V);

  j = 0;

  out(j,0) = XB[XB_sort(g)];
  out(j,1) = pow(expected-observed, 2) / V;

  j++;
  g--;

  arma::uword k;

  while(g > i_lower){

    for(k = lrt.n_rows - 1; ; k--) {
      if (lrt_index[k] == XB_sort[g]) {
        break;
      } else if (lrt_index[k] > XB_sort[g]) {
        k = k+1;
        break;
      }
    }

    observed += status[XB_sort[g]];

    for(; k < lrt.n_rows; k++){

      expected += lrt_deaths[k] / lrt_n_risk[k];

      v2[k] += 1 / lrt_n_risk[k];

      v1[k] += lrt_deaths[k] * (lrt_n_risk[k]-lrt_deaths[k]) /
        (lrt_n_risk[k] * (lrt_n_risk[k]-1));

      lrt_V[k] = v1[k] * (1-v2[k]);

    }

    V = arma::sum(lrt_V);

    out(j,0) = XB[XB_sort[g]];
    out(j,1) = pow(expected-observed, 2) / V;

    j++;
    g--;


  }




  return(out);

}
