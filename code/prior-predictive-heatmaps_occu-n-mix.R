#################################
##  Prior-predictive heatmaps  ##
#################################

library(tidyverse)

# sampling model for binomial RV
gen_binom <- function(z, p, nvisit){
  out <- rbinom(n = nvisit, size = 1, prob = z*p)
}

################################
## Standard occupancy model #### 
################################
# function to generate prior predictive data under 
# Z_i ~ Bern(psi), i = 1,2,3,...,nsite
# Y_{i,1:j}|Z_i ~ Binom(n = nvisit, prob = p*z_i), j = 1,2,...,nvisit
# psi ~ beta(a_psi,b_psi)
# p ~ beta(a_p, b_p)

prior_pred_y <- function(a_psi = 1, b_psi = 1, 
                         a_p = 1, b_p = 1, 
                         nsite = 10, nvisit = 4){
  # priors
  psi <- rbeta(n = 1, shape1 = a_psi, shape2 = b_psi)
  p <- rbeta(n = 1, shape1 = a_p, shape2 = b_p)
  
  # generate data 
  z <- rbinom(n = nsite, size = 1, prob = psi)
  y <- matrix(nrow = nsite, ncol = nvisit)
  y <- t(apply(cbind(z), 1, gen_binom, p = p, nvisit = nvisit))
  out <- list(dg_psi = psi, dg_p = p, y = y, z = z)
  return(out)
}

# prior predictive datasets for user-provided priors on psi and p
prior_pred_occu <- function(nsim, 
                            a_psi = 1, b_psi = 1, 
                            a_p = 1, b_p = 1, 
                            nsite = 10, nvisit = 4){
  pp_occu <- rep(NA, nsim)
  pp_det <- rep(NA, nsim)
  plot_df <- tibble::tibble()
  for(i in 1:nsim){
    prior_gen <- prior_pred_y(a_psi = a_psi, 
                              b_psi = b_psi, 
                              a_p = a_p, b_p = b_p, 
                              nsite = nsite, 
                              nvisit = nvisit)
    
    # these are not correct ests, just thinking about viz 
    ppy <- data.frame(prior_gen$y)
    names(ppy) <- paste0("V", 1:nvisit)
    naive_det <- apply(ppy, 1, sum)/nvisit
    pp_det[i] <- mean(naive_det[which(naive_det > 0)])
    pp_occu[i] <- mean(apply(ppy, 1, sum) > 0)
    
    df_pp <- ppy |>
      dplyr::mutate(site = factor(paste0("S", 1:nsite), 
                                  levels = paste0("S", 1:nsite)), 
                    iter = i, 
                    z_state = prior_gen$z) |>
      tidyr::pivot_longer(cols = 1:nvisit, names_to = "visit", 
                          values_to = "y")
    plot_df <- rbind(plot_df, df_pp)
  }
  out <- list(plot_df = plot_df,
              pp_occu = pp_occu,
              pp_det = pp_det)
  return(out)
}

# function for plotting heatmaps
pp_plot <- function(df){
  df |> ggplot2::ggplot(ggplot2::aes(x = visit, y = site, fill = factor(y))) +
    ggplot2::geom_tile() + scale_fill_viridis_d() +
    labs(fill = "Det/non-det") + 
    facet_wrap(~iter) 
}

# Recreate figure 1
set.seed(7292025)
pp_heat_occu <- prior_pred_occu(20)

pp_heat_occu$plot_df |> 
  pp_plot() + 
  theme_bw(base_size = 20) 

# Recreate figure 2 
# function for summarizing the proportion of 
# sites without variation in the detection history for 
# a given set of prior predictive datasets
prior_pred_novar <- function(plot_df, ...){
  site_visit <- plot_df |> 
    group_by(site, iter) |> 
    summarise(num_y = length(unique(y)))
  df_novar <- site_visit |> 
    group_by(iter) |> 
    summarise(prop_novar = table(num_y)[1]/length(unique(site)))
  out <- df_novar |> tibble()
  return(out)
}


set.seed(72920252)
# generate 100 pp datasets with n = 10, J = 4
# p ~ beta(1,1) and psi ~ beta(1,1) (default settings for args in prior_pred_occu)
pp_occu_novar <- prior_pred_occu(100)

# run prior-predictive assesment with T = # of sites with dh of (0,0,0,0) or (1,1,1,1)
novar <- prior_pred_novar(pp_occu_novar$plot_df)

# create visualization (Fig 2)
novar |> ggplot(aes(x = prop_novar)) + 
  geom_histogram(bins = 10) + 
  theme_minimal(base_size = 20) +
  xlab("Proportion of sites with no variation (all 0 or all 1)") + 
  xlim(c(0,1)) 


######################### 
## N-mixture models. ####
#########################

# sampling model for typical N-mixture models
# the paper focuses on the two options for the distribution 
# of the state process model Poisson/Binomial and NegativeBinomial/Binomial 
# N-mixture models.
# this function is general enough to also consider Poisson/Poisson N-mixture models
# Poisson/binomial 
# N ~ Poisson(lambda)
# Y|N= n ~ Binomial(n, p)

# Negative Binomial state 
# N ~ NegativeBinomial(mu sigma)
# Y|N= n ~ Binomial(n, p)

# Poisson/Poisson 
# N ~ Piosson(lamba)
# Y|N = n ~ Piosson(phi*lambda)
gen_nmix <- function(nsite = 10, nvisit = 4, p = 0.5, phi = 0.5, sig_nb = 1, 
                     lambda = 12, lhood = "pois", obs = "binom"){
  if (lhood == "pois"){ 
    n <- rpois(n  = nsite, lambda = lambda)
    if (obs == "binom"){
      ct <- apply(cbind(n),1, function(x){rbinom(n = nvisit, size = x, prob = p)}) 
      out <- list(y = t(ct), n = n, p = p, lambda = lambda)
    } else{
      ct <- apply(cbind(n),1, function(x){rpois(n = nvisit, lambda = x*phi)})
      out <- list(y = t(ct), n = n, phi = phi, lambda = lambda)
    } 
  } else {
    # Neg binomial process model
    n <- rnbinom(n = nsite, size = sig_nb, mu = lambda)
    # binomial obs model 
    ct <- apply(cbind(n),1, function(x){rbinom(n = nvisit, size = x, prob = p)})  
    out <- list(y = t(ct), n = n, lambda = lambda, sig_nb = sig_nb)
  }
  
  message(paste("state model:", lhood))
  message(paste("obs model:", obs))
  return(out)
}

# prior pred nmix; "default priors" from Kery and Royle (2016)
# lambda ~ gamma(0.001, 0.001)
# p ~ beta(1,1)
# no priors in book for NB because they use RE to induce OD in their Bayesian analysis, 
# spAbundance uses, mu ~ gamma(0.001, 0.001) sigma ~ unif(0,100) 

## FUNCTION: prior_pred_y_nmix
##      args: hyper prior values for all potential priors 
##.           in typical nmix models std, negbinom, pois-pois
##      RETURNS: prior predicted data based on model specification
prior_pred_y_nmix <- function(a_lambda = 0.001, b_lambda = 0.001, 
                              a_p = 1, b_p = 1, 
                              nsite = 10, nvisit = 4, 
                              a_phi = 1, b_phi = 1, 
                              a_sig_nb = 0, b_sig_nb = 100,  
                              lhood = "pois", obs = "binom"){
  # generate pp data for pois-binom or pois-pois
  if (lhood == "pois"){ 
    # get draw from prior for pois samp model 
    lambda <- rgamma(n = 1, shape = a_lambda, rate = b_lambda)
    if (obs == "binom"){
      # pois binom prior
      p <- rbeta(n = 1, shape1 = a_p, shape2 = b_p)
      # generate pp data
      pp_dat <- gen_nmix(nsite = nsite, nvisit = nvisit, p = p, lambda = lambda, 
                         lhood = lhood, obs = obs)
      out <- pp_dat
    } else{
      # pois-pois priors
      phi <- rbeta(n = 1, shape1 = a_phi, shape2 = b_phi)  
      pp_dat <- gen_nmix(nsite = nsite, nvisit = nvisit, phi = phi , lambda = lambda, 
                         lhood = lhood, obs = obs)
      out <- pp_dat
    } 
    # generate pp data for nbinom-binom
  } else {
    # Neg binomial priors 
    lambda <- rgamma(n = 1, shape = a_lambda, rate = b_lambda)
    # prior for sigma used in spAbundance (doser 2023) is uniform 0, 100 
    sig_nb <- runif(1, a_sig_nb, b_sig_nb)
    # prior for p in detection process
    p <- rbeta(n = 1, shape1 = a_p, shape2 = b_p)
    # binomial obs model 
    pp_dat <- gen_nmix(nsite = nsite, nvisit = nvisit, p = p, lambda = lambda, 
                       sig_nb = sig_nb, 
                       lhood = lhood, obs = obs)
    out <- pp_dat
  }
  return(out)
}

## Generate many prior predictive dists
prior_pred_nmix <- function(nsim, 
                            a_lambda = 0.001, b_lambda = 0.001, 
                            a_p = 1, b_p = 1, 
                            a_phi = 1, b_phi = 1, 
                            a_sig_nb = NULL, b_sig_nb = NULL, 
                            nsite = 10, nvisit = 4, lhood = "pois", obs = "binom"){
  pp_n <- matrix(NA, nrow = nsim, ncol = nsite)
  pp_det <- rep(NA, nsim)
  plot_df <- tibble::tibble()
  for(i in 1:nsim){
    prior_gen <- suppressMessages(prior_pred_y_nmix(a_lambda = a_lambda, b_lambda = b_lambda,
                                                    a_p = a_p, b_p = b_p, nsite = nsite, nvisit = nvisit, 
                                                    a_phi = a_phi, b_phi = b_phi, 
                                                    a_sig_nb = a_sig_nb, b_sig_nb = b_sig_nb, 
                                                    lhood = lhood, obs = obs))
    
    ppy <- data.frame(prior_gen$y)
    names(ppy) <- paste0("V", 1:nvisit)
    pp_n[i, ] <- apply(prior_gen$y,1, mean) # prior_gen matrix will allow mean b/c homogenous class
    df_pp <- ppy |>
      dplyr::mutate(site = factor(paste0("S", 1:nsite), 
                                  levels = paste0("S", 1:nsite)), 
                    iter = i, 
                    n_state = prior_gen$n) |>
      tidyr::pivot_longer(cols = 1:nvisit, names_to = "visit", 
                          values_to = "y")
    plot_df <- rbind(plot_df, df_pp)
  }
  out <- list(plot_df = plot_df,
              pp_n = pp_n,
              pp_det = pp_det)
  message(paste("state model:", lhood))
  message(paste("obs model:", obs))
  return(out)
}

pp_plot_nmix <- function(df){
  df |> ggplot2::ggplot(ggplot2::aes(x = visit, y = site, fill = factor(y))) +
    ggplot2::geom_tile() + scale_fill_viridis_d() +
    labs(fill = "Count") + 
    facet_wrap(~iter)
}

# Recreate figure 3
# generate 20 pp datasets with n = 10, J = 4
# p ~ beta(1,1) and mu ~ gamma(0.001, 0.001) (default settings for args in prior_pred_occu)
set.seed(52827)
pp_nmix <- prior_pred_nmix(nsim = 20)
pp_nmix$plot_df |> 
  pp_plot_nmix()+ 
  theme_bw(base_size = 20) 

# Another set of 20 with default priors 
set.seed(52826)
test <- prior_pred_nmix(nsim = 20)
test$plot_df |> 
  pp_plot_nmix() + 
  theme_bw(base_size = 20) 


# Prior predictive check, max count
prior_pred_max_count <- function(plot_df, ...){
  df_max <- plot_df |> 
    group_by(iter) |> 
    summarise(max_y = max(y))
  return(df_max)
}

# prior check, 100 pp datasets Poisson-Nmix with default priors
set.seed(923)
pp_max <- prior_pred_nmix(nsim = 1000)

df_pp_max <- prior_pred_max_count(pp_max$plot_df)

df_pp_max |> 
  ggplot(aes(x = max_y)) + 
  geom_histogram()


# Negative binomial model for N
set.seed(73026)
test <- prior_pred_nmix(nsim = 20, lhood = "nbinom", a_sig_nb = 0, b_sig_nb = 100)

# get all 0s again..,
test$plot_df |> 
  pp_plot_nmix() + 
  theme_bw(base_size = 20) 

### Comparison of max count sin pp datasets across Poisson-Nmix and NBinom-Nmix 
# prior check, 100 pp datasets Poisson-Nmix with default priors
pp_max_nb <- prior_pred_nmix(nsim = 100, lhood = "nbinom", a_sig_nb = 0, b_sig_nb = 100)

pp_max_nb$plot_df |> 
  prior_pred_max_count() |>
  ggplot(aes(x = max_y)) + 
  geom_histogram()
