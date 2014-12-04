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
  PARAMETER(loga);
  PARAMETER(x);
  PARAMETER(logFm);
  PARAMETER(logWinf);
  PARAMETER(logWfs);
  PARAMETER(logSigma);
  Type sigma = exp(logSigma);
  Type u = 10.0;
  vector<Type> Nvec(nwc);
  Type Fm = exp(logFm);
  Type Winf = exp(logWinf);
  Type Wfs = exp(logWfs);
  Type a = exp(loga);
  ADREPORT(Fm);
  ADREPORT(Winf);
  ADREPORT(Wfs);
  ADREPORT(a);
  ADREPORT(sigma);
  Type cumsum, nc, w, psi_m, psi_F, psi_S, g, m, N, wr, alpha, ssb;
  wr = 0.001;
  cumsum=0.0;
  nc = 0.0;
  ssb = 0.0;
  for(int j=0; j<nwc; j++) {
    w = binsize * (j + 1);
    psi_m = 1 / (1 + pow(w / (Winf * eta_m), -10));
    psi_F = 1 / (1 + pow(w / (Wfs), -u));
    psi_S = 1 / (1 + pow(w / (Winf * 0.001), -u));
    g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m));
    m = a * A * pow(w, n-1);
    cumsum += (m + Fm *  psi_F) / g * binsize; 
    N = exp(-cumsum) / g / wr;
    alpha = epsilon_r * (1 - epsilon_a) * A * pow(Winf, n-1) / wr;
    ssb += psi_m  * N * w * binsize;
    Nvec(j) = N * Fm * (isSurvey == 0 ? psi_F : psi_S);
    nc += Nvec(j) * binsize;
  }
  Type Rrel = 1 - (pow(Winf, 1-n) * wr) / (epsilon_r * (1 - epsilon_a) * A * ssb);
  N = N * Rrel;
  Type ff = freq.sum() * binsize;
  Type nll=0.0;
  for(int i=0; i<nwc; i++) {
    //if(freq(i) > 0) 
    {
      nll -= dnorm(freq(i) / ff, 
                   Nvec(i) / nc, sigma, true);
    }
  }
  nll -= dnorm(loga, meanloga, sdloga, true);
  nll += pow(x, 2);
  vector<Type> residuals(nwc);
  residuals = Nvec / nc - freq / ff;
  REPORT(Nvec);
  REPORT(freq);
  REPORT(residuals);
  REPORT(Rrel);
  return nll;
}

