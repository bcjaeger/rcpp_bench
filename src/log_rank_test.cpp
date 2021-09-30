
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
void survdiff2(int    nn,
               int    nngroup,
               double rho,
               NumericVector time,
               IntegerVector status,
               IntegerVector group,
               NumericVector obs,
               NumericVector exp,
               NumericVector var,
               NumericVector risk,
               NumericVector kaplan){

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
double log_rank_test(arma::mat& y,
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

    while(time[j] == time[i]) {

      deaths += status[j];

      if(group[j] == 2){
        risk += 1;
        obs += status[j] * wt;
      }

      if(j==0){break;} else {j--;}

    }

    i = j+1;
    nrisk = n - i;

    if (deaths>0) {


      exp += wt * deaths * risk / nrisk;
      Rcout << "deaths: " << deaths << std::endl;
      Rcout << "risk: " << risk << std::endl;
      Rcout << "nrisk: " << nrisk << std::endl;
      Rcout << "exp: " << exp << std::endl;

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
