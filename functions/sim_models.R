library(tidyverse)
library(broom)

fdi_data <- read_rds("data/fdi_minimum_data.rds")

iterations <- 99
out <- NULL

fdi_data_pos <- fdi_data %>% filter(FDI_AM > 0)

for(i in 1:iterations){
dt <- fdi_data_pos %>% mutate(x = rnorm(n()))

cat("\nfitting:")

cat(" log-normal")
out <- bind_rows(out,
  glm(FDI_AM ~ x + host + source + factor(year), data = dt, family = gaussian(link = "log")) %>% 
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-normal"))

cat(" | log-gamma")
out <- bind_rows(out,
  glm(FDI_AM ~ x + host + source + factor(year), data = dt, family = Gamma(link = "log")) %>%
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-gamma"))

cat(" | log-poisson")
out <- bind_rows(out,
  glm(FDI_AM ~ x + host + source + factor(year), data = dt, family = poisson(link = "log")) %>%
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-poisson"))

cat(" | log-quasipoisson")
out <- bind_rows(out,
  glm(FDI_AM ~ x + host + source + factor(year), data = dt, family = quasipoisson(link = "log")) %>%
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-quasipoisson"))

cat(" | log-negative-binomial")
out <- bind_rows(out,
  MASS::glm.nb(FDI_AM ~ x + host + source + factor(year), data = dt) %>% 
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-negative-binomial"))

cat(" | log-quasi-const")
out <- bind_rows(out,
                 glm(FDI_AM ~ x + host + source + factor(year), data = dt, family = quasi(link = "log", variance = "constant")) %>% 
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-quasi-const"))

cat(" | log-quasi-mu")
out <- bind_rows(out,
  glm(FDI_AM ~ x + host + source + factor(year), data = dt, family = quasi(link = "log", variance = "mu")) %>%
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-quasi-mu"))

cat(" | log-quasi-mu^2")
out <- bind_rows(out,
  glm(FDI_AM ~ x + host + source + factor(year), data = dt, family = quasi(link = "log", variance = "mu^2")) %>%
    tidy() %>% filter(term == "x") %>% mutate(iter = i, model = "log-quasi-mu^2"))

cat("\niteration:", i, "of", iterations,"complete.\n")
}; rm(i)

out <- out %>% mutate(date = date()) %>% bind_rows(read_rds("functions/temp/sim_out_pos.rds"))
write_rds(out, "functions/temp/sim_out_pos.rds")