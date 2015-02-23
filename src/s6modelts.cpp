#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(binsize);
  DATA_INTEGER(nwc);
  DATA_MATRIX(freq);
  int nyrs = freq.transpose().col(1).size();
  DATA_SCALAR(n);
  DATA_SCALAR(epsilon_a);
  DATA_SCALAR(epsilon_r);
  DATA_SCALAR(A);
  DATA_SCALAR(eta_m);
  DATA_SCALAR(meanloga);
  DATA_SCALAR(sdloga);
  DATA_INTEGER(isSurvey);
  DATA_INTEGER(usePois);
  DATA_VECTOR(totalYield);
  PARAMETER(loga);
  PARAMETER(x);
  PARAMETER_VECTOR(logFm);
  PARAMETER(logWinf);
  PARAMETER_VECTOR(logWfs);
  PARAMETER_VECTOR(logSigma);
  PARAMETER_VECTOR(logeta_S);
//  PARAMETER(logsdFm);
//  PARAMETER(logsdWfs);
  Type u = 10.0;
  vector<Type> Fm = exp(logFm);
  Type Winf= exp(logWinf);
  vector<Type> Wfs = exp(logWfs);
  vector<Type> eta_S = exp(logeta_S);
  vector<Type> sigma = exp(logSigma);
  Type a = exp(loga);
//  Type sdFm = exp(logsdFm);
//  Type sdWfs = exp(logsdWfs);
  Type cumsum, nc, w, psi_m, psi_F, psi_S, g, m, N, wr, alpha, Y, Rrel;
  vector<Type> ssb(nyrs);
  vector<Type> rmax(nyrs);
  vector<Type> R(nyrs);
  wr = 0.001;
  Type nll = 0.0;
  
  for(int yr = 0; yr < nyrs; ++yr) {
     cumsum=0.0;
     nc = 0.0;
     ssb(yr) = 0.0;
     Y = 0.0;
     vector<Type> Nvec(nwc);
     Rrel = 0.0;
    for(int j=0; j<nwc; j++) {
      w = binsize * (j + 0.5);
      psi_m = 1 / (1 + pow(w / (Winf * eta_m), -10));
      psi_F = 1 / (1 + pow(w / (Wfs(yr)), -u));
      psi_S = 1 / (1 + pow(w / (Winf * eta_S(yr)), -u));
      g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m));
      m = a * A * pow(w, n-1);
      cumsum += (m + Fm(yr) *  psi_F) / g * binsize; 
      N = exp(-cumsum) / g;
      alpha = epsilon_r * (1 - epsilon_a) * A * pow(Winf, n-1) / wr;
      ssb(yr) += psi_m  * N * w * binsize;
      Nvec(j) = N * (isSurvey ? psi_S : psi_F);
      Y +=  Fm(yr) * N * psi_F * w * binsize;
      nc += Nvec(j);
    }
    
    Rrel = 1 - (pow(Winf, 1-n) * wr) / (epsilon_r * (1 - epsilon_a) * A * ssb(yr));
    N = N * Rrel;
    Y = Y * Rrel;
    rmax(yr) = totalYield(yr) / Y;
    ssb(yr) *=  Rrel * rmax(yr);
    R(yr) = Rrel * rmax(yr);
    for(int i=0; i<nwc; i++) {
      if(usePois) {
        if(freq(i, yr) > 0) 
        {
          nll -= dpois(freq(i, yr), Nvec(i) / nc * sigma(yr), true);  
        }
      } else {
        nll -= dnorm(freq(i, yr) / freq.sum(), Nvec(i) / nc, sigma(yr), true);
      }
    }
  }
  nll -= dnorm(loga, meanloga, sdloga, true);
  nll += pow(x, 2);
//  for(int i=0; i<Fm.size()-1; ++i){
//    nll -= dnorm(Fm(i+1), Fm(i), sdFm, true);
//    nll -= dnorm(Wfs(i+1), Wfs(i), sdWfs, true);
//  }

  ADREPORT(Fm);
  ADREPORT(Winf);
  ADREPORT(Wfs);
  ADREPORT(eta_S);
  ADREPORT(a);
  ADREPORT(sigma);
  ADREPORT(Y);
  ADREPORT(ssb);
  ADREPORT(R);
  // ADREPORT(sdFm);
  // ADREPORT(sdWfs);
  
  REPORT(Rrel);
  REPORT(rmax);
  return nll;
}

