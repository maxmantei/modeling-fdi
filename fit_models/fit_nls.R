library(tidyverse)
library(rstan)
source("functions/calc_se.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fdi_data <- readRDS("data/fdi_minimum_data.rds")

fdi_data_subset <- fdi_data %>% 
  filter(FDI_AM > 0)

################################################################################
#------------------- HOST, SOURCE, YEAR ----------------------------------------

fdi_data_subset_fe <- fdi_data_subset %>%
  group_by(host) %>% 
  mutate(n_host = n(), 
         fe_host = if_else(n_host < 2, "_restriction", host)) %>%
  group_by(source) %>% 
  mutate(n_source = n(), 
         fe_source = if_else(n_source < 2, "_restriction", source)) %>% 
  group_by(year) %>% 
  mutate(n_year = n()) %>% 
  ungroup() %>% 
  mutate(fe_host_int = as.integer(factor(fe_host)),
         fe_source_int = as.integer(factor(fe_source)),
         fe_year_int = as.integer(factor(year)))

standata <- list(
  N = nrow(fdi_data_subset_fe),
  N_host = max(fdi_data_subset_fe$fe_host_int),
  N_source = max(fdi_data_subset_fe$fe_source_int),
  N_year = max(fdi_data_subset_fe$fe_year_int),
  FDI_AM = fdi_data_subset_fe$FDI_AM,
  host = fdi_data_subset_fe$fe_host_int,
  source = fdi_data_subset_fe$fe_source_int,
  year = fdi_data_subset_fe$fe_year_int
)

standata$priors <- 0
empty_stanmodel <- stan("models/nls_hsy.stan", data = standata, iter = 0, chains = 0)

standata$priors <- 0
mle_fit <- optimizing(get_stanmodel(empty_stanmodel), data = standata, 
                      verbose = TRUE, iter = 1e4, 
                      as_vector = FALSE, hessian = TRUE)


#------------------- HOSTYEAR, SOURCEYEAR, DYAD --------------------------------

fdi_data_subset_fe <- fdi_data_subset %>% 
  mutate(hostyear = paste(host, year, sep = "_"),
         sourceyear = paste(source, year, sep = "_")) %>%
  group_by(hostyear) %>% 
  mutate(n_hostyear = n(), 
         fe_hostyear = if_else(n_hostyear < 2, "_restriction", hostyear)) %>%
  group_by(sourceyear) %>% 
  mutate(n_sourceyear = n(), 
         fe_sourceyear = if_else(n_sourceyear < 2, "_restriction", sourceyear)) %>% 
  group_by(dyad) %>% 
  mutate(n_dyad = n(), 
         fe_dyad = if_else(n_dyad < 2, "_restriction", dyad)) %>% 
  ungroup() %>% 
  mutate(fe_hostyear_int = as.integer(factor(fe_hostyear)),
         fe_sourceyear_int = as.integer(factor(fe_sourceyear)),
         fe_dyad_int = as.integer(factor(fe_dyad)))

standata <- list(
  N = nrow(fdi_data_subset_fe),
  N_hostyear = max(fdi_data_subset_fe$fe_hostyear_int),
  N_sourceyear = max(fdi_data_subset_fe$fe_sourceyear_int),
  N_dyad = max(fdi_data_subset_fe$fe_dyad_int),
  FDI_AM = fdi_data_subset_fe$FDI_AM,
  hostyear = fdi_data_subset_fe$fe_hostyear_int,
  sourceyear = fdi_data_subset_fe$fe_sourceyear_int,
  dyad = fdi_data_subset_fe$fe_dyad_int
)

#------------------------ FITTING MODEL ----------------------------------------

standata$priors <- 0
empty_stanmodel <- stan("models/nls_hysyd.stan", data = standata, iter = 0, chains = 0)

standata$priors <- 0
mle_fit <- optimizing(get_stanmodel(empty_stanmodel), data = standata, 
                      verbose = TRUE, iter = 1e5, 
                      as_vector = FALSE, hessian = TRUE)

results <- tibble(parameter = names(unlist(mle_fit$par)),
       estimate = unlist(mle_fit$par),
       se = calc_se(mle_fit, stan_fit = empty_stanmodel, n_sims = 10000))

results_extended <- results %>% 
  mutate(z_score = abs(estimate/se), 
         p_value = 2*pnorm(z_score, lower.tail = FALSE),
         ci_90_lo = estimate - se*qnorm(0.05, lower.tail = FALSE),
         ci_90_hi = estimate + se*qnorm(0.05, lower.tail = FALSE))

results_extended %>% filter(parameter)
