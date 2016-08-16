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
  DATA_INTEGER(estFmsy);
  PARAMETER(loga);
  PARAMETER(logFm);
  PARAMETER(logWinf);
  PARAMETER(logWfs);
  PARAMETER(logSigma);
  PARAMETER(logeta_S);
  PARAMETER(logu);
  PARAMETER(logFmsy);
  Type Fmsy = exp(logFmsy);
  Type u = exp(logu);
  Type sigma = exp(logSigma);
  Type Fm = exp(logFm);
  Type Winf = exp(logWinf);
  Type Wfs = exp(logWfs);
  Type a = exp(loga);
  Type eta_S = exp(logeta_S);
  vector<Type> Nvec(nwc);
  vector<Type> Weight(nwc);
  vector<Type> residuals(nwc);
  Type cumsum, nc, w, psi_m, psi_F, psi_S, g, m, N, wr, alpha, ssb, rmax, R, Rp, Yfmsy;
  wr = 0.001;
  cumsum=0.0;
  nc = 0.0;
  ssb = 0.0;
  Type Y = 0.0;
  for(int j = 0; j < nwc; j++) {
    w = binsize * (j + 0.5);
    Weight(j) = w;
    psi_m = 1 / (1 + pow(w / (Winf * eta_m), -10));
    psi_F = 1 / (1 + pow(w / (Wfs), -u));
    psi_S = 1 / (1 + pow(w / (Winf * eta_S), -u));
    g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m));
    m = a * A * pow(w, n - 1);
    cumsum += (m + Fm *  psi_F) / g * binsize; 
    N = exp(-cumsum) / g;
    ssb += psi_m  * N * w * binsize;
    Nvec(j) = N * (isSurvey == 0 ? psi_F : psi_S);
    Y +=  Fm * N * psi_F * w * binsize;
    nc += Nvec(j); // * binsize;
  }
  Type Rrel = 1 - (pow(Winf, 1-n) * wr) / (epsilon_r * (1 - epsilon_a) * A * ssb);
  Y = Y * Rrel;
  rmax = totalYield / Y;
  ssb = ssb * Rrel * rmax;
  R = Rrel * rmax;
  alpha = epsilon_r * (1 - epsilon_a) * A * pow(Winf, n-1) / wr;
  Rp = alpha * ssb;
  Type nll=0.0;
  for(int i=0; i<nwc; i++) {
    if(usePois) {
      if(freq(i) > 0) {
        nll -= dpois(freq(i), Nvec(i) / nc * sigma, true);
        residuals(i) = qnorm(ppois(freq(i), Nvec(i) / nc * sigma));
      }
    } else {
      if(freq(i) > 0) {
        nll -= dnorm(log(freq(i) / freq.sum()), log(Nvec(i) / nc), sigma, true);
        residuals(i) = freq(i) / freq.sum() - Nvec(i) / nc;
      }
    }
  }
  nll -= dnorm(loga, meanloga, sdloga, true);
  
  if(estFmsy) {
    Type w_r = 0.001;
    Type B = 0.0;
    Type psifNwdelta = 0.0;
    Type delta;
    w = w_r;
    g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * epsilon_a);
    m = a * A * pow(w, n - 1);
    N = 1 / g; 
    psi_m = 0;
    cumsum = 0.0;
    int M = 1000;
    for(int i = 1; i < M; i++) {
      w = exp(log(w_r) + i * (log(Winf) - log(w_r)) / (M - 1.0));
      delta = w - exp(log(w_r) + (i - 1) * (log(Winf) - log(w_r)) / (M - 1.0));
      psi_m = 1 / (1 + pow(w / (eta_m * Winf), -10));
      psi_F = 1 / (1 + pow(w / Wfs, -u));
      cumsum += m / g * delta;
      N = exp(-cumsum) / g;
      psifNwdelta += psi_F * N * w * delta;
      B += psi_m * N * w * delta;
      m = a * A * pow(w, n - 1) + Fmsy * psi_F;
      g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m));
    } 
    Rrel = 1 - (pow(Winf, 1 - n) * w_r / (epsilon_r * (1 - epsilon_a) * A * B));
    Yfmsy =  Fmsy * Rrel * psifNwdelta;
    
    nll -= Yfmsy;
  }
  ADREPORT(Fm);
  ADREPORT(Winf);
  ADREPORT(Wfs);
  ADREPORT(eta_S);
  ADREPORT(a);
  ADREPORT(sigma);
  ADREPORT(Y);
  ADREPORT(ssb);
  ADREPORT(R);
  ADREPORT(u);
  ADREPORT(Rrel);
  ADREPORT(Rp);
  ADREPORT(rmax);
  ADREPORT(Fmsy);
  ADREPORT(Fm/Fmsy);
  ADREPORT(Yfmsy);  
  
  REPORT(residuals);
  REPORT(Weight);
  REPORT(freq);
  REPORT(Nvec);

  return nll;
}

