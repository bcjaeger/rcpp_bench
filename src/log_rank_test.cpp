
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
void survdiff2(int    nn,
               int    nngroup,
               double rho,
               NumericVector& time,
               IntegerVector& status,
               IntegerVector& group,
               NumericVector& obs,
               NumericVector& exp,
               NumericVector& var,
               NumericVector& risk,
               NumericVector& kaplan){

  register int i,j,k;
  int kk;
  int n, ngroup, ntot;
  int istart, koff;
  double km, nrisk, wt, tmp;
  double deaths;

  ntot = nn;
  n = ntot;
  ngroup = nngroup;

  istart=0;
  koff=0;

  for (i=0; i< ngroup*ngroup; i++)  var[i]=0;

  for (i=0; i< ngroup; i++) {
    obs[i]=0;
    exp[i]=0;
  }

  while (istart < ntot) {

    for (i=0; i<ngroup; i++) risk[i]=0;

    /*
     ** Compute the k-m, which is only needed if rho!=0
     **   We want it set up as a left-continuous function (unusual)
     */
    if (rho !=0){

      km = 1;

      for (i = istart; i<n; ) {

        kaplan[i] = km;
        nrisk = n - i;
        deaths = status[i];

        for (j=i+1; j<n && time[j]==time[i]; j++) {

          kaplan[j] = km;
          deaths += status[j];

        }

        km = km * (nrisk - deaths) / nrisk;

        i = j;

      }

    }

    /*
     ** Now for the actual test
     */

    for (i = n-1; i >= istart; i--){

      if (rho ==0){
        wt = 1;
      } else {
        wt = pow(kaplan[i], rho);
      }

      deaths = 0;

      for (j = i; j >= istart && time[j] == time[i]; j--) {
        k = group[j]-1;
        deaths += status[j];
        risk[k] += 1;
        obs[k + koff] += status[j] * wt;
      }

      i = j + 1;
      nrisk = n - i;

      if (deaths>0) {
        /* a death time */
        for (k=0; k<ngroup; k++)
          exp[k+koff] += wt * deaths * risk[k] / nrisk;

        // Rcout << "deaths: " << deaths << std::endl;
        // Rcout << "risk: " << risk[1] << std::endl;
        // Rcout << "nrisk: " << nrisk << std::endl;
        // Rcout << "exp: " << exp[1] << std::endl;
        // Rcout << "-----" << std::endl;

        if (nrisk>1){
          kk =0;
          wt = wt*wt;
          for (j=0; j<ngroup; j++) {
            tmp = wt* deaths* risk[j]* (nrisk-deaths)/(nrisk *(nrisk-1));
            var[kk+j] += tmp;
            for (k=0; k<ngroup; k++) {
              var[kk] -= tmp * risk[k] / nrisk;
              kk++ ;
            }
          }
        }

      }

    }

    istart = n;
    koff += ngroup;

  }
}

// [[Rcpp::export]]
double log_rank_test_v1(arma::mat& y,
                        arma::vec& group,
                        double rho){


  arma::uword i, j, n;
  double km, nrisk, wt, tmp, deaths;

  arma::vec time = y.unsafe_col(0);
  arma::vec status = y.unsafe_col(1);

  n = y.n_rows;

  double obs=0, exp=0, V=0, risk=0;

  arma::vec kaplan(n);

  /*
   ** Compute the k-m, which is only needed if rho!=0
   **   We want it set up as a left-continuous function (unusual)
   */
  if (rho !=0){

    km = 1;

    for (i = 0; i < n; ) {

      kaplan[i] = km;
      nrisk = n - i;
      deaths = status[i];

      for (j = i+1; j < n && time[j] == time[i]; j++) {

        kaplan[j] = km;
        deaths += status[j];

      }

      km = km * (nrisk - deaths) / nrisk;

      i = j;

    }

  }

  // Now for the actual test

  bool break_loop = false;

  i = n;

  nrisk = 0;

  while (!break_loop){

    i--;


    if(i == 0) break_loop = true;

    if (rho ==0){
      wt = 1;
    } else {
      wt = pow(kaplan[i], rho);
    }

    deaths = 0;

    j = i;

    for (; time[j] == time[i]; j--) {

      nrisk++;

      deaths += status[j];

      if(group[j] == 2){
        risk += 1;
        obs += status[j] * wt;
      }

      if(j == 0){
        break_loop = true;
        break;
      }

    }

    i = j+1;

    if (deaths>0) {

      exp += wt * deaths * risk / nrisk;

      // Rcout << "e1: " << exp;
      // Rcout << "; o1: " << deaths;
      // Rcout << "; Y1: " << nrisk;
      // Rcout << "; Y: " << risk;
      // Rcout << std::endl;

      // Rcout << "deaths: " << deaths << std::endl;
      // Rcout << "risk: " << risk << std::endl;
      // Rcout << "nrisk: " << nrisk << std::endl;
      // Rcout << "exp: " << exp << std::endl;
      // Rcout << "i: " << i << std::endl;
      // Rcout << "j: " << j << std::endl;
      // Rcout << "break_loop: " << break_loop << std::endl;
      // Rcout << "-----" << std::endl;

      // update variance if nrisk > 1 (if nrisk == 1, variance is 0)
      if (nrisk>1){

        tmp = wt*wt * deaths * risk * (nrisk-deaths) / (nrisk * (nrisk-1));
        V = V + tmp * (1 - risk/nrisk);

      }

    }

  }

  //Rcout << "obs: " << obs << std::endl;
  //Rcout << "V:   " << V << std::endl;

  return(pow(exp-obs, 2) / V);


}

// [[Rcpp::export]]
double log_rank_test_v2(arma::mat& y_node,
                        arma::vec& group,
                        double rho){

  arma::uword i, j;
  double n, km, nrisk, wt, tmp, deaths;
  double observed=0, expected=0, V=0, risk=0;
  bool break_loop;

  n = y_node.n_rows;

  break_loop = false;
  i = n;
  nrisk = 0;

  if(rho == 0){

    while (!break_loop){

      i--;

      if(i == 0) break_loop = true;

      deaths = 0;

      j = i;

      for (; y_node.at(j, 0) == y_node.at(i, 0); j--) {

        nrisk++;

        deaths += y_node.at(j, 1);

        if(group[j] == 2){
          risk += 1;
          observed += y_node.at(j, 1);
        }

        if(j == 0){
          break_loop = true;
          break;
        }

      }

      i = j+1;

      if (deaths>0) {

        expected += deaths * risk / nrisk;
        // Rcout << "deaths: " << deaths << std::endl;
        // Rcout << "risk: " << risk << std::endl;
        // Rcout << "nrisk: " << nrisk << std::endl;
        // Rcout << "expected: " << expected << std::endl;
        // Rcout << "i: " << i << std::endl;
        // Rcout << "j: " << j << std::endl;
        // Rcout << "break_loop: " << break_loop << std::endl;
        // Rcout << "-----" << std::endl;

        // update variance if nrisk > 1 (if nrisk == 1, variance is 0)
        if (nrisk>1){

          tmp = deaths * risk * (nrisk-deaths) / (nrisk * (nrisk-1));
          V = V + tmp * (1 - risk/nrisk);

        }

      }

    };
  } else {

    // Compute the k-m, which is only needed if rho!=0
    //   We want it set up as a left-continuous function (unusual)

    arma::vec kaplan(n);

    km = 1;

    for (i = 0; i < n; ) {

      kaplan[i] = km;
      nrisk = n - i;
      deaths = y_node(i, 1);

      for (j = i+1; j < n && y_node.at(j,0) == y_node.at(i,0); j++) {

        kaplan[j] = km;
        deaths += y_node.at(j, 1);

      }

      km = km * (nrisk - deaths) / nrisk;

      i = j;

    }

    nrisk = 0;

    while (!break_loop){

      i--;

      if(i == 0) break_loop = true;

      wt = kaplan[i];

      deaths = 0;

      j = i;

      for (; y_node.at(j, 0) == y_node.at(i, 0); j--) {

        nrisk++;

        deaths += y_node.at(j, 1);

        if(group[j] == 2){
          risk += 1;
          observed += y_node.at(j, 1) * wt;
        }

        if(j == 0){
          break_loop = true;
          break;
        }

      }

      i = j+1;

      if (deaths>0) {

        expected += wt * deaths * risk / nrisk;
        // Rcout << "deaths: " << deaths << std::endl;
        // Rcout << "risk: " << risk << std::endl;
        // Rcout << "nrisk: " << nrisk << std::endl;
        // Rcout << "expected: " << expected << std::endl;
        // Rcout << "i: " << i << std::endl;
        // Rcout << "j: " << j << std::endl;
        // Rcout << "break_loop: " << break_loop << std::endl;
        // Rcout << "-----" << std::endl;

        // update variance if nrisk > 1 (if nrisk == 1, variance is 0)
        if (nrisk>1){

          tmp = wt*wt * deaths*risk * (nrisk-deaths) / (nrisk*(nrisk-1));
          V = V + tmp * (1 - risk/nrisk);

        }

      }

    }
  }

  //Rcout << "observed: " << observed << std::endl;
  //Rcout << "V:   " << V << std::endl;

  return(pow(expected-observed, 2) / V);


}

// [[Rcpp::export]]
double log_rank_test_v3(arma::mat& y_node,
                        arma::vec& group){

  double temp1, temp2, deaths;
  double n_risk=0, observed=0, expected=0, V=0, g_risk=0;
  bool break_loop = false;

  arma::uword j = y_node.n_rows-1;

  while (!break_loop){

    temp1 = y_node.at(j, 0);

    deaths = 0;

    for (; y_node.at(j, 0) == temp1; j--) {

      n_risk++;
      deaths += y_node.at(j, 1);
      g_risk += group[j];
      observed += y_node.at(j, 1) * group[j];

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

// [[Rcpp::export]]
arma::mat log_rank_test_v4(arma::mat& y_node,
                           arma::vec& group){

  double temp1, deaths;
  double n_risk=0, g_risk=0;
  bool break_loop = false;

  arma::uword j = y_node.n_rows-1;
  arma::mat out(j, 4);
  arma::uword i = 0;


  while (!break_loop){

    temp1 = y_node.at(j, 0);

    deaths = 0;

    for (; y_node.at(j, 0) == temp1; j--) {

      n_risk++;
      deaths   += y_node.at(j, 1);
      g_risk   += group[j];

      if(j == 0){
        break_loop = true;
        break;
      }

    }

    // should only do these calculations if deaths > 0,
    // but turns out its faster to multiply by 0 than
    // it is to check whether deaths is > 0

    out(i, 0) = j+1;
    out(i, 1) = g_risk;
    out(i, 2) = n_risk;
    out(i, 3) = deaths;
    i++;

  }


  out.resize(i, 4);

  if(out.at(i-1,0)==1) out.at(i-1,0)=0;

  return(out);


}

// [[Rcpp::export]]
double log_rank_test_v0(arma::mat& y,
                        arma::vec& g){

  arma::uword i, lwr, upr, count;
  arma::uword n = y.n_rows;
  arma::uword Y = n;

  double Y1 = arma::sum(g);

  lwr = 0;
  upr = 0;

  // i dont need two stopping conditions for this loop
  // because there is guaranteed to be at least 1 event
  // if this function is called.

  i=0;

  while(y.at(i, 1)==0) {upr++;}

  // for(i=0; i<n; i++){
  //
  //   if(y.at(i, 1) == 0){
  //
  //     upr++;
  //
  //   } else {
  //
  //     break;
  //
  //   }
  //
  // }

  count = upr + 1; // starts at 1 to mimic size

  double d = 0;
  double d1 = 0;

  for(i = lwr; i <= upr; i++){
    d += y.at(i, 1);
    d1 += y.at(i, 1) * g[i];
  }

  double e1 = Y1 * d / Y;
  double o1 = d1;

  double V = (Y-Y1) * Y1 * d * (Y-d) / (pow(Y,2) * (Y-1));

  Y -= count;

  Y1 -= arma::sum(g.subvec(lwr, upr));

  //for(i = lwr; i <= upr; i++){Y1 -= g[i];}

  lwr=upr+1;

  for( ; ; ){

    upr = lwr;
    count = 1;

    while( (y.at(upr, 1) == 0) & (upr < n-1) ){
      upr++;
      count++;
    }

    if(upr == n-1){

      // Rcout << "upper == n-1" << std::endl;

      if( y.at(upr, 1) == 0 ){

        break;

      } else {

        d = 0; d1 = 0;

        for(i = lwr; i <= upr; i++){
          d += y.at(i, 1);
          d1 += y.at(i, 1) * g[i];
        }

        e1 += (Y1 * d/Y);
        o1 += d1;

        V += (Y-Y1) * Y1 * d * (Y-d) / (pow(Y, 2) * (Y-1));
        Y -= count;

        Y1 -= arma::sum(g.subvec(lwr, upr));
        //for(i = lwr; i <= upr; i++) Y1 -= g[i];

        break;

      }

    }

    d = 0; d1 = 0;

    for(i = lwr; i <= upr; i++){
      d += y.at(i, 1);
      d1 += y.at(i, 1) * g[i];
    }

    // Rcout << "d: " << d << std::endl;
    // Rcout << "d1: " << d1 << std::endl;

    e1 += (Y1*d / Y);
    o1 += d1;

    V += (Y-Y1) * Y1 * d * (Y-d) / (pow(Y,2) * (Y-1));
    Y -= count;

    Y1 -= arma::sum(g.subvec(lwr, upr));
    //for(i = lwr; i <= upr; i++) Y1 -= g[i];

    // Rcout << "indices: " << lwr << " - " << upr << std::endl;
    lwr=upr+1;

    // Rcout << "e1: " << e1 << "; ";
    // Rcout << "o1: " << o1 << "; ";
    // Rcout << "Y1: " << Y1 << "; ";
    // Rcout << "Y:  " << Y  << "; ";
    // Rcout << "V:  " << V  <<  std::endl;
    // Rcout << std::endl;

    if(Y==1){

      // Rcout << "lwr: " << lwr << std::endl;

      e1 += (Y1 * y.at(lwr, 1) / Y);
      o1 += y.at(lwr, 1) * g[lwr];

      break;

    }

  }


  return pow(o1-e1,2) / V;

}
