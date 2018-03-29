functions{
  real inv_gaussian_lpdf(vector x, vector mu, real lambda){
    real n_half = num_elements(x)*0.5;
    return -n_half*lambda*sum(((x-mu) .* (x-mu)) ./ (x .* (mu .* mu))) -n_half*log(2*pi()) -1.5*log(sum(x)) +n_half*log(lambda);
  }
}
data{
  int N;
  int N_host;
  int N_source;
  int N_year;
  vector<lower=0>[N] FDI_AM;
  int host[N];
  int source[N];
  int year[N];
  int<lower=0,upper=1> priors;
}
transformed data{
  vector[N] FDI_AM_transformed = FDI_AM ./ 1000;
}
parameters{
  real<lower=0> phi;
  
  real alpha;
  vector[N_host-2] alpha_host;
  vector[N_source-1] alpha_source;
  vector[N_year-1] alpha_year;
}
transformed parameters{
  vector<lower=0>[N] mu;
  vector[N_host] fe_host;
  vector[N_source] fe_source;
  vector[N_year] fe_year;
  
  fe_host[1:2] = rep_vector(0, 2);
  fe_host[3:N_host] = alpha_host;
  fe_source[1] = 0;
  fe_source[2:N_source] = alpha_source;
  fe_year[1] = 0;
  fe_year[2:N_year] = alpha_year;
  
  mu = exp(alpha + fe_host[host] + fe_source[source] + fe_year[year]);
}
model{

  if(priors){
    phi ~ double_exponential(0, 1);
    
    alpha ~ normal(0, 10);
    
    alpha_host ~ normal(0, 5);
    alpha_source ~ normal(0, 5);
    alpha_year ~ normal(0, 5);
  } 
  
  target += inv_gaussian_lpdf(FDI_AM_transformed | mu, phi);
  
}
generated quantities{
}
