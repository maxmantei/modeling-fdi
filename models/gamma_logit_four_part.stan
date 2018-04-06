data{
  int N;
  int N_hostyear;
  int N_sourceyear;
  int N_dyad;
  vector[N] FDI_AM;
  int<lower=0,upper=N_hostyear> hostyear[N];
  int<lower=0,upper=N_sourceyear> sourceyear[N];
  int<lower=0,upper=N_dyad> dyad[N];
  int<lower=0,upper=1> priors;
}
transformed data{
  int FDI_AM_zero[N];
  for(n in 1:N)
    FDI_AM_zero[n] = (FDI_AM[n] < 0.01 && FDI_AM[n] > -0.01) ? 1 : 0;
}
parameters{
  real<lower=0,upper=1> theta_zero;
  real<lower=0,upper=1> theta_positive;
  
  real<lower=0> phi;
  
  real alpha_neg;
  
  real alpha;
  vector[N_hostyear-2] alpha_hostyear;
  vector[N_sourceyear-2] alpha_sourceyear;
  vector[N_dyad-2] alpha_dyad;
}
transformed parameters{
  vector<lower=0>[N] mu;

  vector[N_hostyear] fe_hostyear;
  vector[N_sourceyear] fe_sourceyear;
  vector[N_dyad] fe_dyad;

  fe_hostyear[1:2] = rep_vector(0, 2);
  fe_hostyear[3:N_hostyear] = alpha_hostyear;
  fe_sourceyear[1:2] = rep_vector(0, 2);
  fe_sourceyear[3:N_sourceyear] = alpha_sourceyear;
  fe_dyad[1:2] = rep_vector(0, 2);
  fe_dyad[3:N_dyad] = alpha_dyad;
  
  for(n in 1:N){
    if(FDI_AM_zero[n] == 1){
      mu[n] = 0;
  } else{
      if(FDI_AM[n] > 0){
        mu[n] = exp(alpha + fe_hostyear[hostyear[n]] + fe_sourceyear[sourceyear[n]] + fe_dyad[dyad[n]]);
      } else{
        mu[n] = exp(alpha_neg)/exp(alpha + fe_hostyear[hostyear[n]] + fe_sourceyear[sourceyear[n]] + fe_dyad[dyad[n]]);
      }
    }
  }

}
model{
  
  for(n in 1:N){
    if(FDI_AM_zero[n] == 1){
      target += bernoulli_lpmf(1 | theta_zero);
    } else{
      target += bernoulli_lpmf(0 | theta_zero);
      if(FDI_AM[n] > 0){
        target += bernoulli_lpmf(1 | theta_positive);
        target += gamma_lpdf( FDI_AM[n] | phi, phi/mu[n]);
      } else{
        target += bernoulli_lpmf(0 | theta_positive);
        target += gamma_lpdf(-FDI_AM[n] | phi, phi/mu[n]);
      }
    }
  }
  
  if(priors){
    phi ~ normal(0, 1);
  }

}
generated quantities{
}
