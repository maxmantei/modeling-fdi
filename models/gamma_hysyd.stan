data{
  int N;
  int N_hostyear;
  int N_sourceyear;
  int N_dyad;
  vector<lower=0>[N] FDI_AM;
  int hostyear[N];
  int sourceyear[N];
  int dyad[N];
  int<lower=0,upper=1> priors;
}
transformed data{
  //vector[N] log_FDI_AM = log(FDI_AM);
}
parameters{
  real<lower=0> phi;
  
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
  
  mu = exp(alpha + fe_hostyear[hostyear] + fe_sourceyear[sourceyear] + fe_dyad[dyad]);
}
model{
  
  if(priors){
    phi ~ cauchy(0, 1);
  
    alpha ~ normal(0, 10);
  
    alpha_hostyear ~ normal(0, 5);
    alpha_sourceyear ~ normal(0, 5);
    alpha_dyad ~ normal(0, 5);
  } 
  
  target += gamma_lpdf(FDI_AM | phi, rep_vector(phi, N)./mu);

}
generated quantities{
}
