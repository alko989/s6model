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
  PARAMETER_VECTOR(logFm);
  PARAMETER(logWinf);
  PARAMETER_VECTOR(logWfs);
  PARAMETER_VECTOR(logSigma);
  PARAMETER(logeta_S);
  PARAMETER(logu);
  Type u = exp(logu);
//  PARAMETER(logsdFm);
//  PARAMETER(logsdWfs);
  vector<Type> Fm = exp(logFm);
  Type Winf = exp(logWinf);
  vector<Type> Wfs = exp(logWfs);
  Type eta_S = exp(logeta_S);
  vector<Type> sigma = exp(logSigma);
  matrix<Type> residuals(nwc,nyrs);
  matrix<Type> Weight(nwc,nyrs);
  matrix<Type> Nvec(nwc,nyrs);
  Type a = exp(loga);
//  Type sdFm = exp(logsdFm);
//  Type sdWfs = exp(logsdWfs);
  Type cumsum, nc, w, psi_m, psi_F, psi_S, g, m, N, wr, alpha;
  vector<Type> ssb(nyrs);
  vector<Type> rmax(nyrs);
  vector<Type> R(nyrs);
  vector<Type> Rrel(nyrs);
  vector<Type> Y(nyrs);
  vector<Type> Rp(nyrs);
  wr = 0.001;
  Type nll = 0.0;
  
  for(int yr = 0; yr < nyrs; ++yr) {
     cumsum=0.0;
     nc = 0.0;
     ssb(yr) = 0.0;
     Y(yr) = 0.0;
     Rrel(yr) = 0.0;
    for(int j=0; j<nwc; j++) {
      w = Weight(j, yr) = binsize * (j + 0.5);
      psi_m = 1 / (1 + pow(w / (Winf * eta_m), -10));
      psi_F = 1 / (1 + pow(w / (Wfs(yr)), -u));
      psi_S = 1 / (1 + pow(w / (Winf * eta_S), -u));
      g = A * pow(w, n) * (1 - pow(w / Winf, 1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m));
      m = a * A * pow(w, n-1);
      cumsum += (m + Fm(yr) *  psi_F) / g * binsize; 
      N = exp(-cumsum) / g;
      ssb(yr) += psi_m  * N * w * binsize;
      Nvec(j, yr) = N * (isSurvey ? psi_S : psi_F);
      Y(yr) +=  Fm(yr) * N * psi_F * w * binsize;
      nc += Nvec(j, yr);
    }
    Rrel(yr) = 1 - (pow(Winf, 1-n) * wr) / (epsilon_r * (1 - epsilon_a) * A * ssb(yr));
    Y(yr) = Y(yr) * Rrel(yr);
    rmax(yr) = totalYield(yr) / Y(yr);
    ssb(yr) *=  Rrel(yr) * rmax(yr);
    R(yr) = Rrel(yr) * rmax(yr);
    alpha = epsilon_r * (1 - epsilon_a) * A * pow(Winf, n-1) / wr;
    Rp(yr)  = alpha * ssb(yr);
    for(int i=0; i<nwc; i++) {
      if(usePois) {
        if(freq(i, yr) > 0) 
        {
          nll -= dpois(freq(i, yr), Nvec(i) / nc * sigma(yr), true);
          residuals(i, yr) = ppois(freq(i, yr), Nvec(i) / nc * sigma(yr)); // freq(i, yr) - Nvec(i) / nc * sigma(yr);
        }
      } else {
        nll -= dnorm(freq(i, yr) / freq.sum(), Nvec(i) / nc, sigma(yr), true);
        residuals(i,yr) = freq(i, yr) / freq.sum() - Nvec(i) / nc;
      }
    }
  }
  nll -= dnorm(loga, meanloga, sdloga, true);
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
  ADREPORT(u);
  ADREPORT(Rrel);
  ADREPORT(Rp);
  ADREPORT(rmax);
  // ADREPORT(sdFm);
  // ADREPORT(sdWfs);
  
  REPORT(residuals);
  REPORT(Weight);
  REPORT(freq);
  REPORT(Nvec);
  return nll;
}

