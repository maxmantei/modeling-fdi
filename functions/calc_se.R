calc_se <- function(mle_fit, stan_fit, n_sims = 10000){
  
cat("Extract the vector of estimates and the Hessian for the unconstrained parameters\n")
theta_hat <- unlist(mle_fit$par)
upars <- unconstrain_pars(empty_stanmodel, relist(theta_hat, mle_fit$par))
Hessian <- mle_fit$hessian
cat("There are", sum(diag(-Hessian) <= 0), "problematic parameters!")

cat("\nExtract the Cholesky decomposition of the (negative) Hessian and invert\n")
cat("This might take a while...\n")
cat("  - doing Cholesky decomposition...\n")
R <- as.matrix(bdsmatrix::gchol(-Hessian))

cat("  - inverting...\n")
V <- chol2inv(R)
rownames(V) <- colnames(V) <- colnames(Hessian)

cat("Produce a matrix with some specified number of simulation draws from a multinormal\n")
sims <- n_sims
len <- length(upars)
cat("This might take a while...\n")
cat("  - producing", len, "by", sims, "simulation matrix\n")

robust <- FALSE
if(robust){
  t_scaling <- 1/sqrt(rchisq(sims * len, 22688-5691)/(22688-5691))
} else{
  t_scaling <- 1
}

unconstrained <- upars + t(chol(V)) %*% matrix(rnorm(sims * len)*t_scaling, nrow = len, ncol = sims)

cat("  - convert back to constrained parameter space\n")
theta_sims <- t(apply(unconstrained, 2, FUN = 
                        function(upars) {
                          unlist(constrain_pars(empty_stanmodel, upars))
                        }
))

cat("Produce estimated standard errors for the constrained parameters\n")
se <- apply(theta_sims, 2, sd)
se_estimated_params <- se[names(unlist(mle_fit$par)) %in% str_remove(colnames(Hessian), "[.]")]
se_estimated_params[as.vector(diag(-Hessian) <= 0)] <- NA
se[names(unlist(mle_fit$par)) %in% str_remove(colnames(Hessian), "[.]")] <- se_estimated_params
se[which(se == 0)] <- NA
return(se)
}