#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(binsize);
  DATA_INTEGER(nwc);
  DATA_VECTOR(freq);
  DATA_SCALAR(n);
  DATA_SCALAR(epsilon_a);
  DATA_SCALAR(epsilon_r);
  DATA_SCALAR(A);
  DATA_SCALAR(eta_m);
  DATA_SCALAR(meanloga);
  DATA_SCALAR(sdloga);
  DATA_INTEGER(isSurvey);
  DATA_INTEGER(usePois);
  DATA_SCALAR(totalYield);
  PARAMETER(loga);
  PARAMETER(x);
  PARAMETER(logFm);
  PARAMETER(logWinf);
  PARAMETER(logWfs);
  PARAMETER(logSigma);
  PARAMETER(logeta_S);
  Type sigma = exp(logSigma);
  Type u = 10.0;
  vector<Type> Nvec(nwc);
  vector<Type> Weight(nwc);
  Type Fm = exp(logFm);
  Type Winf = exp(logWinf);
  Type Wfs = exp(logWfs);
  Type a = exp(loga);
  Type eta_S = exp(logeta_S);
  Type cumsum, nc, w, psi_m, psi_F, psi_S, g, m, N, wr, alpha, ssb, rmax, R;
  wr = 0.001;
  cumsum=0.0;
  nc = 0.0;
  ssb = 0.0;
  Type Y = 0.0;
  for(int j=0; j<nwc; j++) {
    w = binsize * (j + 0.5);
    Weight(j) = w;
    psi_m = 1 / (1 + pow(w / (Winf * eta_m), -10));
    psi_F = 1 / (1 + pow(w / (Wfs), -u));
    psi_S = 1 / (1 + pow(w / (Winf * eta_S), -u));
    g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m));
    m = a * A * pow(w, n-1);
    cumsum += (m + Fm *  psi_F) / g * binsize; 
    N = exp(-cumsum) / g;
    alpha = epsilon_r * (1 - epsilon_a) * A * pow(Winf, n-1) / wr;
    ssb += psi_m  * N * w * binsize;
    Nvec(j) = N * (isSurvey == 0 ? psi_F : psi_S);
    Y +=  Fm * N * psi_F * w * binsize;
    nc += Nvec(j); // * binsize;
  }
  Type Rrel = 1 - (pow(Winf, 1-n) * wr) / (epsilon_r * (1 - epsilon_a) * A * ssb);
  N = N * Rrel;
  Y = Y * Rrel;
  rmax = totalYield / Y;
  ssb = ssb * Rrel * rmax;
  R = Rrel * rmax;
  Type ff = freq.sum() ; // * binsize;
  Type nll=0.0;
  for(int i=0; i<nwc; i++) {
    if(usePois) {
      if(freq(i) > 0) {
        //        Type mean = Nvec(i) / nc * ff;
        //        Type size = pow(mean, 2) / (pow(sigma, 2) - mean);
        //        Type prob = size / (size + mean);
        //        nll -= dnbinom(freq(i), size, prob, true);
        nll -= dpois(freq(i), Nvec(i) / nc * sigma, true);
        //nll += pow(sigma, 2);
      }
    } else {
      if(freq(i) > 0) {
        nll -= dnorm(log(freq(i) / ff), log(Nvec(i) / nc), sigma, true);
      }
    }
  }
  nll -= dnorm(loga, meanloga, sdloga, true);
  nll += pow(x, 2);
  vector<Type> residuals(nwc);
  residuals = Nvec / nc * ff - freq;
  
  
  ADREPORT(Fm);
  ADREPORT(Winf);
  ADREPORT(Wfs);
  ADREPORT(eta_S);
  ADREPORT(a);
  ADREPORT(sigma);
  ADREPORT(Y);
  ADREPORT(ssb);
  ADREPORT(R);
  
  REPORT(Nvec);
  REPORT(Weight);
  REPORT(freq);
  REPORT(residuals);
  REPORT(Rrel);
  REPORT(rmax);
  return nll;
}

