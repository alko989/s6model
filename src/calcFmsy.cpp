#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{ 
  DATA_SCALAR(n);
  DATA_SCALAR(epsilon_a);
  DATA_SCALAR(epsilon_r);
  DATA_SCALAR(A);
  DATA_SCALAR(eta_m);
  DATA_SCALAR(a);
  DATA_SCALAR(Winf);
  DATA_SCALAR(Wfs);
  DATA_INTEGER(md); //From Ken's paper notation
  DATA_SCALAR(u);
  PARAMETER(logF);
  Type Fmsy = exp(logF);
  ADREPORT(Fmsy);
  vector<Type> g(md);
  vector<Type> m(md);
  vector<Type> psi_F(md);
  vector<Type> psi_m(md);
  vector<Type> N(md);
  vector<Type> ww(md);
  vector<Type> delta(md);
  Type cumsum;
  Type w_r;
  Type B;
  Type Rrel;
  Type Y;
  
  w_r=0.001;
  ww(0) = w_r;
  g(0)=A * pow(ww(0),n) * (1 - pow(ww(0)/Winf, 1-n) * (epsilon_a ));
  m(0) = a * A * pow(ww(0), n - 1);
  N(0)=1/g(0); 
  psi_m(0) = 0;

  for(int i=1; i<md; i++) {
    ww(i) = exp(log(w_r) + i * (log(Winf) - log(w_r)) / (md - 1.0));
    delta(i) = ww(i) - ww(i-1);
    psi_m(i) = 1 / (1 + pow(ww(i)/(eta_m * Winf),-10));
    psi_F(i) = 1 / (1 + pow(ww(i)/(Wfs),-u));
    m(i) = a * A * pow(ww(i), n - 1) + Fmsy * psi_F(i); // exp(-(ww(i)/Winf)) + 0.4 * (1 - ww(i) / Winf) + Fmsy * psi_F(i);
    g(i) = A * pow(ww(i),n) * (1 - pow(ww(i)/Winf, 1-n) * (epsilon_a + (1-epsilon_a)*psi_m(i)));
    cumsum += (m(i-1))/g(i-1) * delta(i);
    N(i)=exp(-cumsum)/g(i);
    
  } 
  delta(md-1) = 0;
  B = (psi_m * N * ww * delta).sum();
  Rrel = 1 - (pow(Winf,1-n) * w_r / (epsilon_r * (1 - epsilon_a) * A * B));
  Y =  Fmsy * Rrel * (psi_F * N * ww * delta).sum();
  return - Y;
}
