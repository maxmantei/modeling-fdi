functions{
  real inv_gaussian_lpdf(vector x, vector mu, real lambda){
    real n_half = num_elements(x)*0.5;
    return -n_half*lambda*sum(((x-mu) .* (x-mu)) ./ (x .* (mu .* mu))) -n_half*log(2*pi()) -1.5*log(sum(x)) +n_half*log(lambda);
  }
}
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
  vector[N] FDI_AM_transformed = FDI_AM ./ 1000;
}
parameters{
  real log_phi;
  
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
    log_phi ~ normal(-20, 5);
    
    alpha ~ normal(0, 10);
    
    alpha_hostyear ~ normal(0, 5);
    alpha_sourceyear ~ normal(0, 5);
    alpha_dyad ~ normal(0, 5);
  } 
  
  target += inv_gaussian_lpdf(FDI_AM_transformed | mu, exp(log_phi));
  
}
generated quantities{
}
