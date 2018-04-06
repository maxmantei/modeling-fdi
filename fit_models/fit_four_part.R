library(tidyverse)
library(rstan)
source("functions/calc_se.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fdi_data <- readRDS("data/fdi_minimum_data.rds") %>%
  mutate(FDI_AM = if_else(is.na(FDI_AM), 0, FDI_AM))

# Prepare data

num_covariates <- 1

fdi_data_fe <- fdi_data %>% 
  filter(!between(FDI_AM, -0.01, 0.01)) %>% 
  mutate(hostyear = paste(host, year, sep = "_"),
         sourceyear = paste(source, year, sep = "_")) %>%
  group_by(hostyear) %>% 
  mutate(n_hostyear = n(), 
         fe_hostyear = if_else(n_hostyear < num_covariates + 1, "_restriction", hostyear)) %>%
  group_by(sourceyear) %>% 
  mutate(n_sourceyear = n(), 
         fe_sourceyear = if_else(n_sourceyear < num_covariates + 1, "_restriction", sourceyear)) %>% 
  group_by(dyad) %>% 
  mutate(n_dyad = n(), 
         fe_dyad = if_else(n_dyad < num_covariates + 1, "_restriction", dyad)) %>% 
  ungroup() %>% 
  mutate(fe_hostyear_int = as.integer(factor(fe_hostyear)),
         fe_sourceyear_int = as.integer(factor(fe_sourceyear)),
         fe_dyad_int = as.integer(factor(fe_dyad))) %>%
  right_join(fdi_data)
fdi_data_fe[is.na(fdi_data_fe)] <- 0

standata <- list(
  N = nrow(fdi_data_fe),
  N_hostyear = max(fdi_data_fe$fe_hostyear_int),
  N_sourceyear = max(fdi_data_fe$fe_sourceyear_int),
  N_dyad = max(fdi_data_fe$fe_dyad_int),
  FDI_AM = fdi_data_fe$FDI_AM,
  hostyear = fdi_data_fe$fe_hostyear_int,
  sourceyear = fdi_data_fe$fe_sourceyear_int,
  dyad = fdi_data_fe$fe_dyad_int
)

# Fitting the model

standata$priors <- 0
empty_stanmodel <- stan(file = "C:/Users/win7/Desktop/gamma_logit_four_part.stan",
                        data = standata, chains = 0, iter = 0)

standata$priors <- 0
mle_fit <- optimizing(get_stanmodel(empty_stanmodel), data = standata, 
                      verbose = TRUE, iter = 1e4, 
                      as_vector = FALSE, hessian = TRUE)

results <- inference(mle_fit, stan_fit = empty_stanmodel, 
                     ci_level = 0.95, n_sims = 10000,
                     transformations = list(theta_zero = plogis,
                                            theta_positive = plogis,
                                            phi = exp))

resid <- fdi_data_fe %>%
  mutate(estimate = mle_fit$par$mu,
         estimate = if_else(FDI_AM < 0, -estimate, estimate),
         resid = FDI_AM - estimate,
         resid_pearson = resid/sqrt(((estimate^2)/mle_fit$par$phi)))

resid %>% ggplot(aes(x = FDI_AM, y = estimate)) + 
  geom_point(alpha = 0.8, color = "azure3") + 
  geom_abline(intercept = 0, slope = 1, color = "cornflowerblue", size = 1) +
  geom_smooth(color = "coral4") + 
  coord_cartesian(xlim = c(min(resid$FDI_AM), max(resid$FDI_AM)), ylim = c(min(resid$FDI_AM), max(resid$FDI_AM))) +
  xlab(expression(FDI)) +
  ylab(expression(hat(mu))) +
  theme_minimal() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
resid %>% ggplot(aes(x = asinh(FDI_AM), y = asinh(estimate))) + 
  geom_point(alpha = 0.8, color = "azure3") + 
  geom_abline(intercept = 0, slope = 1, color = "cornflowerblue", size = 1) +
  geom_smooth(color = "coral4") +
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) +
  xlab(expression(log(FDI))) +
  ylab(expression(log(hat(mu)))) +
  theme_minimal() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

# only valid for the part of non-zero flows
resid %>% ggplot(aes(x = asinh(estimate), y = resid_pearson)) +
  geom_point(alpha = 0.8, color = "azure3") +
  geom_smooth(color = "coral4") +
  coord_cartesian(ylim = c(-50, 50)) +
  xlab(expression(asinh(hat(mu)))) +
  ylab(expression(hat(epsilon)/hat(theta)*hat(mu)^{2})) +
  theme_minimal() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
resid %>% ggplot(aes(x = asinh(estimate), y = sqrt(abs(resid_pearson)))) +
  geom_point(alpha = 0.8, color = "azure3") +
  geom_smooth(color = "coral4") +
  xlab(expression(asinh(hat(mu)))) +
  ylab(expression(sqrt(abs(hat(epsilon)/hat(theta)*hat(mu)^{2})))) +
  theme_minimal() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))