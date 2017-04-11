## data to pass to JAGS
dat_jags <- c("dat_age","ln_dat_esc","dat_harv","scl_cvrs","n_cov",
              "n_yrs","A","age_min","age_max","age_skip","n_fore") 

## 2. model params/states for JAGS to return
par_jags <- c("alpha","E_Rkr_a","mu_Rkr_a","Rkr_b","ln_Rkr_a",
              "cc",
              "Sp","Rec","tot_ln_Rec","ln_RS","p_vec",
              "var_Qr","var_Rs","res_ln_Rec")

## 3. MCMC control params
# MCMC parameters
mcmc_chains <- 4
mcmc_length <- 5e5
mcmc_burn <- 2.5e5
mcmc_thin <- 500
# total number of MCMC samples
mcmc_samp <- (mcmc_length-mcmc_burn)*mcmc_chains/mcmc_thin

## function to create JAGS inits
init_vals <- function() {
  list(mu_Rkr_a=1, cc=rnorm(n_cov,0,0.1),
       Rkr_b=1/exp(mean(ln_dat_esc, na.rm=TRUE)),
       p_pi=1, p_mu=rep(1,A),
       p_vec=matrix(c(0.01,0.3,0.48,0.15,0.05,0.01),n_yrs-age_min+n_fore,A,byrow=TRUE),
       Rec_mu=log(1000),
       Rec_sig=0.1,
       tot_ln_Rec=rep(log(1000),n_yrs-age_min+n_fore),
       innov_1=0,
       phi=0.5)
}

## list of model info for JAGS
mod_jags <- list(data=dat_jags,
                 inits=init_vals,
                 parameters.to.save=par_jags,
                 model.file=fn_jags,
                 n.chains=as.integer(mcmc_chains),
                 n.iter=as.integer(mcmc_length),
                 n.burnin=as.integer(mcmc_burn),
                 n.thin=as.integer(mcmc_thin),
                 DIC=TRUE)

## fit the model in JAGS & store results
mod_fit <- do.call(jags.parallel, mod_jags)
