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
  // this function returns a cutpoint obtaining a local maximum
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
  double n_events=0, n_risk=0;
  bool break_loop;

  // initialize at the lowest possible LRT stat value
  double stat_current, stat_best = 0;

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
    k=0, p=y_node.n_rows;

  // first determine the lowest value of XB that will
  // be a valid cut-point to split a node. A valid cut-point
  // is one that, if used, will result in at least min_obs
  // and min_events in both the left and right node.
  for(it = XB_sort.begin(); it != XB_sort.end(); ++it){
    n_events += status[*it];
    n_risk++; // keep separate from k - n_risk will be weighted
    k++;
    if(n_events >= min_events && n_risk >= min_obs) break;
  }

  // got to reset these before finding the upper limit of j
  n_events=0; n_risk=0;

  // group should be initialized as all 0s
  group.fill(0);

  for(it = XB_sort.end()-1; it >= XB_sort.begin(); --it){
    n_events += status[*it];
    group[*it] = 1;
    n_risk++;
    // p: location in XB when we take p steps away from the end of XB.
    p--;
    if(n_events >= min_events && n_risk >= min_obs) break;
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

  if ( k > p ) return(R_PosInf);

  // adjust p to indicate steps taken in the outer loop.
  p -= k;

  // adjust k to tell us when we are at a landmark value of j.
  k = round(p / (n_cps));

  // Rcout << "upper index: " << p << std::endl;
  // Rcout << "landmarks at each " << k << " steps" << std::endl;
  // arma::uvec tmp = arma::linspace<arma::uvec>(0, p, n_cps);
  // Rcout << tmp << std::endl;

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

        n_events = 0;

        for ( ; time[i] == temp1; i--) {

          n_risk++;
          n_events += status[i];
          g_risk += group[i];
          observed += status[i] * group[i];

          if(i == 0){
            break_loop = true;
            break;
          }

        }

        // should only do these calculations if n_events > 0,
        // but turns out its faster to multiply by 0 than
        // it is to check whether n_events is > 0

        temp2 = g_risk / n_risk;
        expected += n_events * temp2;

        // update variance if n_risk > 1 (if n_risk == 1, variance is 0)
        // definitely check if n_risk is > 1 b/c otherwise divide by 0
        if (n_risk>1){
          temp1 = n_events * temp2 * (n_risk-n_events) / (n_risk-1);
          V += temp1 * (1 - temp2);
        }

      }

      stat_current = pow(expected-observed, 2) / V;

      // Rcout << "stat_current: " << stat_current << std::endl;
      // Rcout << "J: " << j << std::endl;

      if(stat_current > stat_best){
        it_best = it;
        stat_best = stat_current;
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
  while(it < it_best){
    group[*it] = 0;
    ++it;
  }

  return(XB[*it_best]);

}

// [[Rcpp::export]]
double lrt_multi_v2(arma::mat& y_node,
                     arma::vec& XB,
                     arma::vec& group,
                     arma::uword n_cps,
                     double min_obs,
                     double min_events){

  // about this function - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // this function returns a cutpoint obtaining a local maximum
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
  double n_events=0, n_risk=0;
  bool break_loop;

  // initialize at the lowest possible LRT stat value
  double stat_current, stat_best = 0;

  // sort XB- we need to iterate over the sorted indices
  arma::uvec XB_sort = arma::sort_index(XB);

  // speaking of which, here is that iterator.
  arma::uvec::iterator iit, jit, iit_best;

  // unsafe columns point to specific cols in y_node.
  // this makes the code more readable and doesn't copy data
  arma::vec status = y_node.unsafe_col(1);
  arma::vec time = y_node.unsafe_col(0);

  // unsigned integers used by loops
  arma::uword
    // i is used within the lrt (inner) loop.
    i,
    // j is used to keep track of group assignments
    j=0,
    // p is used to determine the range of j.
    k=y_node.n_rows;

  // first determine the lowest value of XB that will
  // be a valid cut-point to split a node. A valid cut-point
  // is one that, if used, will result in at least min_obs
  // and min_events in both the left and right node.
  for(iit = XB_sort.begin(); iit != XB_sort.end(); ++iit){
    n_events += status[*iit];
    n_risk++;
    j++;
    if(n_events >= min_events && n_risk >= min_obs) break;
  }

  // got to reset these before finding the upper limit of j
  n_events=0; n_risk=0;

  // group should be initialized as all 0s
  group.fill(0);

  for(iit = XB_sort.end()-1; iit >= XB_sort.begin(); --iit){
    n_events += status[*iit];
    group[*iit] = 1;
    n_risk++;
    // k: location in XB when we take p steps away from the end of XB.
    k--;
    if(n_events >= min_events && n_risk >= min_obs) break;
  }

  // this is just done to avoid compilation warnings
  // (iit_best should be initialized before iit is used)
  iit_best = iit;

  // what happens if we don't have enough events or obs to split?
  // the first valid lower cut-point (at XB_sort[k]) is > the first
  // valid upper cutpoint (current value of n_risk). Put another way,
  // k (the number of steps taken from beginning of the XB vec)
  // will be > n_rows - p, where the difference on the RHS is
  // telling us where we are after taking p steps from the end
  // of the XB vec. Returning the infinite cp is a red flag.

  if ( j > k ) return(R_PosInf);

  // adjust p to indicate steps taken in the outer loop.
  k -= j;

  j = 0;

  // Rcout << "upper index: " << p << std::endl;
  // Rcout << "landmarks at each " << k << " steps" << std::endl;
  arma::uvec jvals = arma::linspace<arma::uvec>(0, k, n_cps);
  // Rcout << tmp << std::endl;

  // begin outer loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(jit = jvals.begin(); jit != jvals.end(); ++jit){

    // drop down one spot on XB
    // Rcout << "iit points to" << *iit << std::endl;
    // Rcout << "jit points to" << *jit << std::endl;

    for( ; j < *jit; j++){
      group[*iit] = 1;
      --iit;
    }

    n_risk=0;
    g_risk=0;

    observed=0;
    expected=0;

    V=0;

    break_loop = false;

    i = y_node.n_rows-1;

    // begin inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - -
    while (!break_loop){

      temp1 = time[i];

      n_events = 0;

      for ( ; time[i] == temp1; i--) {

        n_risk++;
        n_events += status[i];
        g_risk += group[i];
        observed += status[i] * group[i];

        if(i == 0){
          break_loop = true;
          break;
        }

      }

      // should only do these calculations if n_events > 0,
      // but turns out its faster to multiply by 0 than
      // it is to check whether n_events is > 0

      temp2 = g_risk / n_risk;
      expected += n_events * temp2;

      // update variance if n_risk > 1 (if n_risk == 1, variance is 0)
      // definitely check if n_risk is > 1 b/c otherwise divide by 0
      if (n_risk>1){
        temp1 = n_events * temp2 * (n_risk-n_events) / (n_risk-1);
        V += temp1 * (1 - temp2);
      }

    }
    // end inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

    stat_current = pow(expected-observed, 2) / V;

    // Rcout << "stat_current: " << stat_current << std::endl;
    // Rcout << "J: " << j << std::endl;

    if(stat_current > stat_best){
      iit_best = iit;
      stat_best = stat_current;
    }


    // end outer loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

  }

  // rewind iit until it is back where it was when we got the
  // best lrt stat. While rewinding iit, also reset the group
  // values so that group is as it was when we got the best
  // lrt stat.
  while(iit < iit_best){
    group[*iit] = 0;
    ++iit;
  }

  return(XB[*iit_best]);

}


// [[Rcpp::export]]
double lrt_multi_v3b(arma::mat& y_node,
                     arma::vec& XB,
                     arma::vec& group,
                     arma::uword n_cps,
                     double min_obs,
                     double min_events){

  // about this function - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // this function returns a cutpoint obtaining a local maximum
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
  double n_events=0, n_risk=0;
  bool break_loop;

  // initialize at the lowest possible LRT stat value
  double stat_current, stat_best = 0;

  // sort XB- we need to iterate over the sorted indices
  arma::uvec XB_sort = arma::sort_index(XB);

  // speaking of which, here is that iterator.
  arma::uvec::iterator iit, jit, iit_best;

  // unsafe columns point to specific cols in y_node.
  // this makes the code more readable and doesn't copy data
  arma::vec status = y_node.unsafe_col(1);
  arma::vec time = y_node.unsafe_col(0);

  // unsigned integers used by loops
  arma::uword
    // i is used within the lrt (inner) loop.
    i,
    // j is used to keep track of group assignments
    j=0,
    // p is used to determine the range of j.
    k=y_node.n_rows;

  // first determine the lowest value of XB that will
  // be a valid cut-point to split a node. A valid cut-point
  // is one that, if used, will result in at least min_obs
  // and min_events in both the left and right node.
  for(iit = XB_sort.begin(); iit != XB_sort.end(); ++iit){
    n_events += status[*iit];
    n_risk++;
    j++;
    if(n_events >= min_events && n_risk >= min_obs) break;
  }

  // got to reset these before finding the upper limit of j
  n_events=0; n_risk=0;

  // group should be initialized as all 0s
  group.fill(0);

  for(iit = XB_sort.end()-1; iit >= XB_sort.begin(); --iit){
    n_events += status[*iit];
    group[*iit] = 1;
    n_risk++;
    // k: location in XB when we take p steps away from the end of XB.
    k--;
    if(n_events >= min_events && n_risk >= min_obs) break;
  }

  // this is just done to avoid compilation warnings
  // (iit_best should be initialized before iit is used)
  iit_best = iit;

  // what happens if we don't have enough events or obs to split?
  // the first valid lower cut-point (at XB_sort[k]) is > the first
  // valid upper cutpoint (current value of n_risk). Put another way,
  // k (the number of steps taken from beginning of the XB vec)
  // will be > n_rows - p, where the difference on the RHS is
  // telling us where we are after taking p steps from the end
  // of the XB vec. Returning the infinite cp is a red flag.

  if ( j > k ) return(R_PosInf);

  // adjust p to indicate steps taken in the outer loop.
  k -= j;

  j = 0;

  // Rcout << "upper index: " << p << std::endl;
  // Rcout << "landmarks at each " << k << " steps" << std::endl;
  arma::uvec jvals = arma::linspace<arma::uvec>(0, k, n_cps);
  // Rcout << tmp << std::endl;

  // begin outer loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(jit = jvals.begin(); jit != jvals.end(); ++jit){

    // drop down one spot on XB
    // Rcout << "iit points to" << *iit << std::endl;
    // Rcout << "jit points to" << *jit << std::endl;

    for( ; j < *jit; j++){
      group[*iit] = 1;
      --iit;
    }

    n_risk=0;
    g_risk=0;

    observed=0;
    expected=0;

    V=0;

    break_loop = false;

    i = y_node.n_rows-1;

    // begin inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - -
    for (; ;){

      temp1 = time[i];

      n_events = 0;

      for ( ; time[i] == temp1; i--) {

        n_risk++;
        n_events += status[i];
        g_risk += group[i];
        observed += status[i] * group[i];

        if(i == 0){
          break_loop = true;
          break;
        }

      }

      // should only do these calculations if n_events > 0,
      // but turns out its faster to multiply by 0 than
      // it is to check whether n_events is > 0

      temp2 = g_risk / n_risk;
      expected += n_events * temp2;

      // update variance if n_risk > 1 (if n_risk == 1, variance is 0)
      // definitely check if n_risk is > 1 b/c otherwise divide by 0
      if (n_risk>1){
        temp1 = n_events * temp2 * (n_risk-n_events) / (n_risk-1);
        V += temp1 * (1 - temp2);
      }

      if(break_loop) break;

    }
    // end inner loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

    stat_current = pow(expected-observed, 2) / V;

    // Rcout << "stat_current: " << stat_current << std::endl;
    // Rcout << "J: " << j << std::endl;

    if(stat_current > stat_best){
      iit_best = iit;
      stat_best = stat_current;
    }


    // end outer loop  - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

  }

  // rewind iit until it is back where it was when we got the
  // best lrt stat. While rewinding iit, also reset the group
  // values so that group is as it was when we got the best
  // lrt stat.
  while(iit < iit_best){
    group[*iit] = 0;
    ++iit;
  }

  return(XB[*iit_best]);

}
