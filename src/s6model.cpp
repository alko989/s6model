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
  PARAMETER(loga);
  PARAMETER(logFm);
  PARAMETER(logWinf);
  PARAMETER(logWfs);
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
  Type cumsum, totalFishing, delta, w, psi_m, psi_F, g, m, N;
  cumsum=0.0;
  totalFishing = 0.0;  
  delta = binsize;
  for(int j=0; j<nwc; j++) {
    w = delta * (j + 1);
    psi_m = 1 / (1 + pow(w / (Winf * eta_m), -10));
    psi_F = 1 / (1 + pow(w / (Wfs), -u));
    g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m));
    m = a * A * pow(w, n-1);
    cumsum += (m + Fm *  psi_F) / g * delta; 
    N = exp(-cumsum) / g;
    Nvec(j) = N * Fm * psi_F;
    totalFishing += N * psi_F * Fm * delta;
  }
  Type nll=0.0;
  for(int i=0; i<nwc; i++) {
    if(Nvec(i) > 0)
      nll += -log(Nvec(i)/totalFishing) * freq(i);
  }
  nll -= dnorm(loga, meanloga, sdloga, true);
  return nll;
}

