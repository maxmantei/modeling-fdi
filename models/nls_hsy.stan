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
  //vector[N] log_FDI_AM = log(FDI_AM);
}
parameters{
  real<lower=0> sigma;
  
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
    (sigma/1000) ~ cauchy(0, 2.5);
  
    alpha ~ normal(0, 10);
  
    alpha_host ~ normal(0, 5);
    alpha_source ~ normal(0, 5);
    alpha_year ~ normal(0, 5);
  } 
  
  target += normal_lpdf(FDI_AM | mu, sigma);

}
generated quantities{
}