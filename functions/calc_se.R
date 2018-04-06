inference <- function(mle_fit, stan_fit, ci_level = 0.9, n_sims = 10000,
                      transformations = list()){
  
require(bdsmatrix)
  
ci_value <- (1 - ci_level)/2
  
cat("Extract the vector of estimates and the Hessian for the unconstrained parameters\n")
upars <- unconstrain_pars(empty_stanmodel, relist(unlist(mle_fit$par), mle_fit$par))
Hessian <- mle_fit$hessian

cat("\nExtract the Cholesky decomposition of the (negative) Hessian and invert\n")
cat("This might take a while...\n")
cat("  - doing Cholesky decomposition...\n")
R <- gchol(-Hessian)

cat("  - inverting...\n")
V <- solve(R, full = FALSE)
rownames(V)

cat("Produce a matrix with some specified number of simulation draws from a multinormal\n")
len <- length(upars)
cat("This might take a while...\n")
cat("  - producing", n_sims, "by", len, "simulation matrix\n")

#unconstrained <- upars + t(chol(V)) %*% matrix(rnorm(sims * len), nrow = len, ncol = sims)
sims <- mvnfast::rmvn(n_sims, upars, V, isChol = TRUE) #this is not very fast...
colnames(sims) <- colnames(Hessian)

stats <- tibble(parameter = names(unlist(mle_fit$par)[1:len]),
                estimate = upars, 
                se = apply(sims, 2, sd)) %>% 
  mutate(ci_low = estimate - se*qnorm(ci_value, lower.tail = FALSE), 
         ci_high = estimate + se*qnorm(ci_value, lower.tail = FALSE))

for(par in names(transformations)){
  stats[stats$parameter == par, "se"] <- sd(transformations[[par]](sims[,par]))
  for(val in c("estimate", "ci_low", "ci_high")){
    stats[stats$parameter == par, val] <- transformations[[par]](as.numeric(stats[stats$parameter == par, val]))
  }
}

#cat("  - convert back to constrained parameter space...")
#constrained <- matrix(nrow = length(theta_hat), ncol = sims)
# for loop is a few seconds faster than apply
#for(column in 1:ncol(unconstrained)){
#  constrained[,column] <- unlist(constrain_pars(empty_stanmodel, unconstrained[,column]))
#  if(column == round(ncol(unconstrained)/2)) cat("50% done!\n")
#}
#
#cat("Produce estimated standard errors for the constrained parameters\n")
#stats <- stats[1:len,]
#stats[which(!(upars %in% theta_hat)),"se"] <- apply(constrained[which(!(upars %in% theta_hat)),], 1, sd)

# set problematic SEs to NA
stats$se[as.vector(diag(-Hessian) <= 0)] <- NA
stats <- stats %>% 
  mutate(parameter = names(unlist(mle_fit$par)[1:len]),
         ci_low = if_else(is.na(se), NA_real_, ci_low),
         ci_high = if_else(is.na(se), NA_real_, ci_high))

return(stats)

}