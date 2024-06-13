myldtnorm <- function(x, mean, sd, lower=-Inf, upper=Inf, log=TRUE){

  if(any(x<lower) || any(x>upper)){
    stop('x outside the range')
  }
  ld <- -0.5*log(2*pi)-log(sd)-(x-mean)^2/2/sd^2
  lp1 <- pnorm(upper, mean = mean, sd = sd, lower.tail = TRUE, log.p = TRUE)
  lp2 <- pnorm(lower, mean = mean, sd = sd, lower.tail = TRUE, log.p = TRUE)
  result <- ld-lp1-log(1-exp(lp2-lp1))

  return(result)

}

# truncated gamma pdf (log)
myldtgamma <- function(x, alpha, beta, lower=0, upper=Inf, log=TRUE){

  if(any(x<lower) || any(x>upper)){
    stop('x outside the range')
  }

  ld <- alpha*log(beta)-lgamma(alpha)+(alpha-1)*log(x)-beta*x
  lp1 <- pgamma(upper, shape = alpha, rate = beta, lower.tail = TRUE, log.p = TRUE)
  lp2 <- pgamma(lower, shape = alpha, rate = beta, lower.tail = TRUE, log.p = TRUE)
  result <- ld-lp1-log(1-exp(lp2-lp1))

  return(result)

}

# ------------------------------- MCMC Gibbs update ------------------------------

# ------------------- Update tau ----------------------
tau_log_prob <- function(tau_cc, u_obs_cc, s_obs_cc, mu_u, mu_s, sigma_u_2, sigma_s_2,
                         mu_tau_cc, var_tau_cc, tau_upper_cc){

  lp <- dtnorm(x = u_obs_cc, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE)+
    dtnorm(x = s_obs_cc, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE)+
    myldtgamma(x = tau_cc, alpha = mu_tau_cc^2/var_tau_cc, beta = mu_tau_cc/var_tau_cc, lower = 0,
               upper = tau_upper_cc, log = TRUE)

  if(is.infinite(lp) || is.na(lp)){
    lp <- myldtnorm(x = u_obs_cc, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE)+
      myldtnorm(x = s_obs_cc, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE)+
      myldtgamma(x = tau_cc, alpha = mu_tau_cc^2/var_tau_cc, beta = mu_tau_cc/var_tau_cc, lower = 0,
                 upper = tau_upper_cc, log = TRUE)
  }

  return(lp)
}

tau_update <- function(tau, u_obs, s_obs, k, alpha1_1, beta, gamma, q, sigma_u_2, sigma_s_2, t0, u01, compute_mean,
                       iter_num, sd_tau, mu_tau, var_tau, target_accept){
  # Number of cells
  N <- nrow(u_obs)

  # The upper bounds for tau (fixed value conditional on the rest)
  # vector of length N
  tau_upper <- ifelse(k==1, t0, Inf)

  # All are vector of length N
  tau_old <- tau; X_old <- -log(1/tau_old-1/tau_upper)
  sd_tau_old <- sd_tau

  # Save updated scaling factor
  sd_tau_new <- rep(NA, N)

  # Save updated tau
  tau_new <- rep(NA, N)

  accept_count <- 0
  n <- iter_num

  # Compute current (old) mean s and u for each cell
  mus_old <- compute_mean(tau = tau_old, t0 = t0, k = k, alpha1_1 = alpha1_1, beta = beta, gamma = gamma, q = q, u01 = u01)

  loop.result <- lapply(1:N, function(cc) {

    # Apply AMH
    X_new <- rnorm(n = 1, mean = X_old[cc], sd = sqrt(sd_tau_old[cc]))
    # Transform back to tau
    tau_new_cc <- 1/(exp(-X_new)+1/tau_upper[cc])

    # Compute new mean for s and u for this cell
    mu_new_cc <- compute_mean(tau = tau_new_cc, t0 = t0, k = k[cc], alpha1_1 = alpha1_1, beta = beta, gamma = gamma,
                              q = q, u01 = u01)

    # Compute log acceptance probability
    log_acceptance <- tau_log_prob(tau_cc = tau_new_cc, u_obs_cc = u_obs[cc], s_obs_cc = s_obs[cc],
                                   mu_u = mu_new_cc$u_mean, mu_s = mu_new_cc$s_mean, sigma_u_2 = sigma_u_2,
                                   sigma_s_2 = sigma_s_2, mu_tau_cc = mu_tau[k[cc]], var_tau_cc = var_tau[k[cc]],
                                   tau_upper_cc = tau_upper[cc])-
      tau_log_prob(tau_cc = tau_old[cc], u_obs[cc], s_obs[cc], mu_u = mus_old$u_mean[cc], mu_s = mus_old$s_mean[cc],
                   sigma_u_2, sigma_s_2, mu_tau[k[cc]], var_tau[k[cc]], tau_upper[cc])+
      ifelse(k[cc]==1, log(tau_new_cc)+log(tau_upper[cc]-tau_new_cc)-log(tau_old[cc])-log(tau_upper[cc]-tau_old[cc]),
             log(tau_new_cc)-log(tau_old[cc]))

    acceptance_tau <- exp(log_acceptance)
    acceptance_tau <- min(1, acceptance_tau)
    log_prob <- min(0,log_acceptance)

    # Decision
    outcome <- runif(1,min=0,max=1)
    if(log(outcome) > log_prob){
      X_new <- X_old[cc]
      tau_new_cc <- tau_old[cc]
      accept <- 0
    }else{
      accept <- 1
    }

    # Update scale parameter in the proposal distribution
    sd_tau_cc <- exp(log(sd_tau_old[cc])+n^(-0.7)*(acceptance_tau-target_accept))
    if(sd_tau_cc>exp(50)) {sd_tau_cc <- exp(50)}
    if(sd_tau_cc<exp(-50)) {sd_tau_cc <- exp(-50)}

    return(list(X_new=X_new, tau_new_cc=tau_new_cc,accept=accept, sd_tau_cc=sd_tau_cc))
  })# end lapply

  # Now move results to the new vectors
  tau_new <- sapply(loop.result, function(l) l$tau_new_cc)
  sd_tau_new <- sapply(loop.result, function(l) l$sd_tau_cc)
  accept_count <- sum(sapply(loop.result, function(l) l$accept))

  return(list(tau_new=tau_new,accept_count=accept_count, sd_tau_new=sd_tau_new))
}

# ------------------- Update prior params for tau: mu_tau, var_tau ----------------------
mu_tau_log_prob <- function(mu_tau, tau, var_tau, tau_upper, mu_star_base, sig_star_base){

  # tau is a vector containing tau in the corresponding state, others are singular values
  lp <- sum(myldtgamma(x = tau, alpha = mu_tau^2/var_tau, beta = mu_tau/var_tau,
                       lower = 0, upper = tau_upper, log = TRUE))+
    dtnorm(x = mu_tau, mean = mu_star_base, sd = sig_star_base, lower = 0, log = TRUE)

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtgamma(x = tau, alpha = mu_tau^2/var_tau, beta = mu_tau/var_tau,
                         lower = 0, upper = tau_upper, log = TRUE))+
      myldtnorm(x = mu_tau, mean = mu_star_base, sd = sig_star_base, lower = 0, log = TRUE)
  }

  return(lp)
}

mu_tau_update <- function(mu_tau, tau, var_tau, k, mu_star_base, sig_star_base, t0,
                          iter_num, sd_mu_tau, target_accept, MH_count){

  # Define inputs
  mu_tau_old <- mu_tau; X_old <- log(mu_tau_old)
  sd_mu_tau_old <- sd_mu_tau

  # Save updated scaling factor
  sd_mu_tau_new <- rep(NA, 2)

  # Save updated tau
  mu_tau_new <- rep(NA, 2)

  n <- iter_num
  accept_count <- rep(NA, 2)
  tau_upper <- c(t0, Inf)
  MH_used <- rep(NA, 2)

  for(j in 1:2){

    if(sum(k==j)==0){ #if the state is empty, draw from the prior
      mu_tau_new[j] <- rtnorm(n = 1, mean = mu_star_base[j], sd = sig_star_base[j], lower = 0, upper = Inf)

      # Transform to X
      X_new <- log(mu_tau_new[j])
      accept_count[j] <- 1

      # No update of the scale parameter in the proposal distribution
      sd_mu_tau_new[j] <- sd_mu_tau_old[j]

      MH_used[j] <- 0
    }else{ #apply adaptive MH

      X_new <- rnorm(n = 1, mean = X_old[j], sd = sqrt(sd_mu_tau_old[j]))

      mu_tau_new[j] <- exp(X_new)

      # Compute log acceptance probability
      log_acceptance <- mu_tau_log_prob(mu_tau = mu_tau_new[j], tau = tau[k==j], var_tau = var_tau[j], tau_upper = tau_upper[j],
                                        mu_star_base = mu_star_base[j], sig_star_base = sig_star_base[j])-
        mu_tau_log_prob(mu_tau = mu_tau_old[j], tau[k==j], var_tau[j], tau_upper[j], mu_star_base[j], sig_star_base[j])+
        log(mu_tau_new[j]) - log(mu_tau_old[j])

      acceptance_mu_tau <- exp(log_acceptance)
      acceptance_mu_tau <- min(1, acceptance_mu_tau)
      log_prob <- min(0,log_acceptance)

      # Decision
      outcome <- runif(1,min=0,max=1)
      if(log(outcome) > log_prob){
        X_new <- X_old[j]
        mu_tau_new[j] <- mu_tau_old[j]
        accept_count[j] <- 0
      }else{
        accept_count[j] <- 1
      }

      # Update scale parameter in the proposal distribution
      sd_mu_tau_new[j] <- exp(log(sd_mu_tau_old[j])+MH_count[j]^(-0.7)*(acceptance_mu_tau-target_accept))
      if(sd_mu_tau_new[j]>exp(50)) {sd_mu_tau_new[j] <- exp(50)}
      if(sd_mu_tau_new[j]<exp(-50)) {sd_mu_tau_new[j] <- exp(-50)}

      # MH step is used
      MH_used[j] <- 1
    } #end considering Gibbs or MH
  }

  return(list(mu_tau_new=mu_tau_new,accept_count=accept_count, sd_mu_tau_new=sd_mu_tau_new, MH_used=MH_used))
}


var_tau_log_prob <- function(var_tau, tau, mu_tau, tau_upper, eta_base, nu_base){

  # tau is a vector containing tau in the corresponding state, others are singular values
  lp <- sum(myldtgamma(x = tau, alpha = mu_tau^2/var_tau, beta = mu_tau/var_tau,
                       lower = 0, upper = tau_upper, log = TRUE))+
    dtnorm(x = var_tau, mean = eta_base, sd = nu_base, lower = 0, log = TRUE)

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtgamma(x = tau, alpha = mu_tau^2/var_tau, beta = mu_tau/var_tau,
                         lower = 0, upper = tau_upper, log = TRUE))+
      myldtnorm(x = var_tau, mean = eta_base, sd = nu_base, lower = 0, log = TRUE)
  }
  return(lp)
}

var_tau_update <- function(var_tau, tau, mu_tau, k, eta_base, nu_base, t0,
                           iter_num, sd_var_tau, target_accept, MH_count){

  # Define inputs
  var_tau_old <- var_tau; X_old <- log(var_tau_old)
  sd_var_tau_old <- sd_var_tau

  # Save updated scaling factor
  sd_var_tau_new <- rep(NA, 2)

  # Save updated tau
  var_tau_new <- rep(NA, 2)

  n <- iter_num
  accept_count <- rep(NA, 2)
  tau_upper <- c(t0, Inf)
  MH_used <- rep(NA, 2)

  for(j in 1:2){

    if(sum(k==j)==0){ #if the state is empty, draw from the prior
      var_tau_new[j] <- rtnorm(n = 1, mean = eta_base[j], sd = nu_base[j], lower = 0, upper = Inf)

      # Transform to X
      X_new <- log(var_tau_new[j])
      accept_count[j] <- 1

      # No update of the scale parameter in the proposal distribution
      sd_var_tau_new[j] <- sd_var_tau_old[j]

      MH_used[j] <- 0
    }else{ #apply adaptive MH

      X_new <- rnorm(n = 1, mean = X_old[j], sd = sqrt(sd_var_tau_old[j]))

      var_tau_new[j] <- exp(X_new)

      # Compute log acceptance probability
      log_acceptance <- var_tau_log_prob(var_tau = var_tau_new[j], tau = tau[k==j], mu_tau = mu_tau[j], tau_upper = tau_upper[j],
                                         eta_base = eta_base[j], nu_base = nu_base[j])-
        var_tau_log_prob(var_tau = var_tau_old[j], tau[k==j], mu_tau[j], tau_upper[j], eta_base[j], nu_base[j])+
        log(var_tau_new[j]) - log(var_tau_old[j])

      acceptance_var_tau <- exp(log_acceptance)
      acceptance_var_tau <- min(1, acceptance_var_tau)
      log_prob <- min(0,log_acceptance)

      # Decision
      outcome <- runif(1,min=0,max=1)
      if(log(outcome) > log_prob){
        X_new <- X_old[j]
        var_tau_new[j] <- var_tau_old[j]
        accept_count[j] <- 0
      }else{
        accept_count[j] <- 1
      }

      # Update scale parameter in the proposal distribution
      sd_var_tau_new[j] <- exp(log(sd_var_tau_old[j])+MH_count[j]^(-0.7)*(acceptance_var_tau-target_accept))
      if(sd_var_tau_new[j]>exp(50)) {sd_var_tau_new[j] <- exp(50)}
      if(sd_var_tau_new[j]<exp(-50)) {sd_var_tau_new[j] <- exp(-50)}

      # MH step is used
      MH_used[j] <- 1
    } #end considering Gibbs or MH

  }

  return(list(var_tau_new=var_tau_new,accept_count=accept_count, sd_var_tau_new=sd_var_tau_new, MH_used=MH_used))
}

# ------------------- Update t0 (or t0[2]-t0[1]) ----------------------
t0_log_prob <- function(t0, u_obs, s_obs, k, mu_u, mu_s, sigma_u_2, sigma_s_2, mu_0, var_0, mu_tau, var_tau){

  # From likelihood of s and u
  if(sum(k==2)==0){ #if no cell in k=2
    lp1 <- 0
  }else{
    lp1 <- sum(dtnorm(u_obs[k==2], mu_u[k==2], sqrt(sigma_u_2), lower = 0, log = TRUE))+
      sum(dtnorm(s_obs[k==2], mu_s[k==2], sqrt(sigma_s_2), lower=0, log = TRUE))

    if(is.infinite(lp1) || is.na(lp1)){
      lp1 <- sum(myldtnorm(u_obs[k==2], mu_u[k==2], sqrt(sigma_u_2), lower = 0, log = TRUE))+
        sum(myldtnorm(s_obs[k==2], mu_s[k==2], sqrt(sigma_s_2), lower=0, log = TRUE))
    }
  }

  # From prior for tau
  if(sum(k==1)==0){ #if no cell in k=1
    lp2 <- 0
  }else{
    # compute P(0 <= tau <= t0)
    la <- pgamma(q = t0, shape = mu_tau[1]^2/var_tau[1], rate = mu_tau[1]/var_tau[1], lower.tail = TRUE, log.p = TRUE)
    lp2 <- -sum(k==1)*la
  }

  # From prior for t0
  lp3 <- dgamma(x = t0, shape = mu_0^2/var_0, rate = mu_0/var_0, log = TRUE)

  return(lp1+lp2+lp3)
}

t0_update <- function(t0, u_obs, s_obs, k, alpha1_1, beta, gamma, q, sigma_u_2, sigma_s_2, tau, u01, compute_mean,
                      iter_num, sd_t0, mu_0, var_0, mu_tau, var_tau, target_accept){

  # Minimum value of t0 for non-zero posterior (fixed value conditional on the rest)
  t0_lower <- ifelse(any(k==1), max(tau[k==1]), 0)

  # Define inputs
  t0_old <- t0; X_old <- log(t0_old-t0_lower)
  sd_t0_old <- sd_t0

  n <- iter_num

  # Compute current (old) mean s and u for each cell
  mus_old <- compute_mean(tau, t0 = t0_old, k, alpha1_1, beta, gamma, q, u01)

  # Apply AMH
  X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(sd_t0_old))

  # Transform back to t0
  t0_new <- exp(X_new)+t0_lower
  # Compute new mean s and u for each cell
  mus_new <- compute_mean(tau, t0 = t0_new, k, alpha1_1, beta, gamma, q, u01)

  # Compute log acceptance probability
  log_acceptance <- t0_log_prob(t0 = t0_new, u_obs = u_obs, s_obs = s_obs, k = k, mu_u = mus_new$u_mean, mu_s = mus_new$s_mean,
                                sigma_u_2 = sigma_u_2, sigma_s_2 = sigma_s_2, mu_0 = mu_0, var_0 = var_0, mu_tau = mu_tau, var_tau = var_tau)-
    t0_log_prob(t0 = t0_old, u_obs, s_obs, k, mu_u = mus_old$u_mean, mu_s = mus_old$s_mean, sigma_u_2, sigma_s_2, mu_0,
                var_0, mu_tau, var_tau) +
    log(t0_new-t0_lower)-log(t0_old-t0_lower)

  acceptance_t0 <- exp(log_acceptance)
  acceptance_t0 <- min(1, acceptance_t0)
  log_prob <- min(0,log_acceptance)

  # Decision
  outcome <- runif(1,min=0,max=1)
  if(log(outcome) > log_prob){
    X_new <- X_old
    t0_new <- t0_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_t0_new <- exp(log(sd_t0_old)+n^(-0.7)*(acceptance_t0-target_accept))
  if(sd_t0_new>exp(50)) {sd_t0_new <- exp(50)}
  if(sd_t0_new<exp(-50)) {sd_t0_new <- exp(-50)}

  return(list(t0_new=t0_new, accept=accept, sd_t0_new=sd_t0_new))
}


# ------------------- Update alpha1_1 ----------------------
alpha1_1_log_prob <- function(alpha1_1, u_obs, s_obs, mu_u, mu_s, sigma_u_2, sigma_s_2, mu_alpha, sig_alpha,
                              alpha1_1_lower){

  lp <- sum(dtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))+
    sum(dtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
    dtnorm(x = alpha1_1, mean = mu_alpha, sd = sig_alpha, lower = alpha1_1_lower, log = TRUE)

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))+
      sum(myldtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
      myldtnorm(x = alpha1_1, mean = mu_alpha, sd = sig_alpha, lower = alpha1_1_lower, log = TRUE)
  }
  return(lp)
}

alpha1_1_update <- function(alpha1_1, u_obs, s_obs, t0, k, beta, gamma, q, sigma_u_2, sigma_s_2, tau, u01, compute_mean,
                            iter_num, sd_alpha1_1, mu_alpha, sig_alpha, target_accept){

  alpha1_1_lower <- beta*u01

  # Define inputs
  alpha1_1_old <- alpha1_1; X_old <- log(alpha1_1_old-alpha1_1_lower)
  sd_alpha1_1_old <- sd_alpha1_1

  n <- iter_num

  # Compute current (old) mean s and u for each cell
  mus_old <- compute_mean(tau, t0, k, alpha1_1=alpha1_1_old, beta, gamma, q, u01)

  # Apply AMH
  X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(sd_alpha1_1_old))

  # Transform back to alpha1_1
  alpha1_1_new <- exp(X_new)+alpha1_1_lower
  # Compute new mean s and u for each cell
  mus_new <- compute_mean(tau, t0, k, alpha1_1=alpha1_1_new, beta, gamma, q, u01)

  # Compute log acceptance probability
  log_acceptance <- alpha1_1_log_prob(alpha1_1 = alpha1_1_new, u_obs = u_obs, s_obs = s_obs, mu_u = mus_new$u_mean,
                                      mu_s = mus_new$s_mean, sigma_u_2 = sigma_u_2, sigma_s_2 = sigma_s_2,
                                      mu_alpha = mu_alpha, sig_alpha = sig_alpha, alpha1_1_lower = alpha1_1_lower) -
    alpha1_1_log_prob(alpha1_1 = alpha1_1_old, u_obs, s_obs, mu_u = mus_old$u_mean, mu_s = mus_old$s_mean, sigma_u_2,
                      sigma_s_2, mu_alpha, sig_alpha, alpha1_1_lower)+
    log(alpha1_1_new-alpha1_1_lower) - log(alpha1_1_old-alpha1_1_lower)

  acceptance_alpha <- exp(log_acceptance)
  acceptance_alpha <- min(1, acceptance_alpha)
  log_prob <- min(0,log_acceptance)

  # Decision
  outcome <- runif(1,min=0,max=1)
  if(log(outcome) > log_prob){
    X_new <- X_old
    alpha1_1_new <- alpha1_1_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_alpha1_1_new <- exp(log(sd_alpha1_1_old)+n^(-0.7)*(acceptance_alpha-target_accept))
  if(sd_alpha1_1_new>exp(50)) {sd_alpha1_1_new <- exp(50)}
  if(sd_alpha1_1_new<exp(-50)) {sd_alpha1_1_new <- exp(-50)}

  return(list(alpha1_1_new=alpha1_1_new, accept=accept, sd_alpha1_1_new=sd_alpha1_1_new))
}

# ------------------- Update gamma ----------------------
gamma_log_prob <- function(gamma, s_obs, mu_s, sigma_s_2, mu_gamma, sig_gamma){

  lp <- sum(dtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
    dtnorm(x = gamma, mean = mu_gamma, sd = sig_gamma, lower = 0, log = TRUE)

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
      myldtnorm(x = gamma, mean = mu_gamma, sd = sig_gamma, lower = 0, log = TRUE)
  }
  return(lp)
}

gamma_update <- function(gamma, s_obs, t0, k, alpha1_1, beta, q, sigma_s_2, tau, u01, compute_mean,
                         iter_num, sd_gamma, mu_gamma, sig_gamma, target_accept){

  # Define inputs
  gamma_old <- gamma; X_old <- log(gamma_old)
  sd_gamma_old <- sd_gamma

  n <- iter_num

  # Compute current (old) mean s and u for each cell
  mus_old <- compute_mean(tau, t0, k, alpha1_1, beta, gamma=gamma_old, q, u01)

  # Apply AMH
  X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(sd_gamma_old))

  # Transform back to alpha1_1
  gamma_new <- exp(X_new)
  # Compute new mean s and u for each cell
  mus_new <- compute_mean(tau, t0, k, alpha1_1, beta, gamma=gamma_new, q, u01)

  # Compute log acceptance probability
  log_acceptance <- gamma_log_prob(gamma = gamma_new, s_obs = s_obs, mu_s = mus_new$s_mean, sigma_s_2 = sigma_s_2,
                                   mu_gamma = mu_gamma, sig_gamma = sig_gamma)-
    gamma_log_prob(gamma = gamma_old, s_obs, mu_s = mus_old$s_mean, sigma_s_2, mu_gamma, sig_gamma)+
    log(gamma_new) - log(gamma_old)

  acceptance_gamma <- exp(log_acceptance)
  acceptance_gamma <- min(1, acceptance_gamma)
  log_prob <- min(0,log_acceptance)

  # Decision
  outcome <- runif(1,min=0,max=1)
  if(log(outcome) > log_prob){
    X_new <- X_old
    gamma_new <- gamma_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_gamma_new <- exp(log(sd_gamma_old)+n^(-0.7)*(acceptance_gamma-target_accept))
  if(sd_gamma_new>exp(50)) {sd_gamma_new <- exp(50)}
  if(sd_gamma_new<exp(-50)) {sd_gamma_new <- exp(-50)}

  return(list(gamma_new=gamma_new, accept=accept, sd_gamma_new=sd_gamma_new))

}

# ------------------- Update q ----------------------
q1_log_prob <- function(q1, u_obs, s_obs, mu_u, mu_s, sigma_u_2, sigma_s_2, q_lower, q_upper){

  lp <- sum(dtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))+
    sum(dtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
    dunif(x = q1, min = q_lower, max = q_upper, log = TRUE)

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))+
      sum(myldtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
      dunif(x = q1, min = q_lower, max = q_upper, log = TRUE)
  }

  return(lp)

}

q1_update <- function(q, u_obs, s_obs, alpha1_1, t0, k, beta, gamma, sigma_u_2, sigma_s_2, tau, u01, compute_mean,
                      iter_num, sd_q1, q_lower, q_upper, target_accept){


  # Define inputs
  q1_old <- q; X_old <- log(q1_old-q_lower)-log(q_upper-q1_old)
  sd_q1_old <- sd_q1

  n <- iter_num

  # Compute current (old) mean s and u for each cell
  mus_old <- compute_mean(tau, t0, k, alpha1_1, beta, gamma, q=q1_old, u01)

  # Apply AMH
  X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(sd_q1_old))

  # Transform back to q
  q1_new <- q_upper+(q_lower-q_upper)/(1+exp(X_new))

  # Compute new mean s and u for each cell
  mus_new <- compute_mean(tau, t0, k, alpha1_1, beta, gamma, q=q1_new, u01)

  # Compute log acceptance probability
  log_acceptance <- q1_log_prob(q1 = q1_new, u_obs = u_obs, s_obs = s_obs, mu_u = mus_new$u_mean, mu_s = mus_new$s_mean,
                                sigma_u_2 = sigma_u_2, sigma_s_2 = sigma_s_2, q_lower = q_lower, q_upper = q_upper)-
    q1_log_prob(q1 = q1_old, u_obs, s_obs, mu_u = mus_old$u_mean, mu_s = mus_old$s_mean, sigma_u_2, sigma_s_2,
                q_lower, q_upper)+
    log(q1_new-q_lower)+log(q_upper-q1_new)-
    log(q1_old-q_lower)-log(q_upper-q1_old)

  acceptance_q1 <- exp(log_acceptance)
  acceptance_q1 <- min(1, acceptance_q1)
  log_prob <- min(0,log_acceptance)

  # Decision
  outcome <- runif(1,min=0,max=1)
  if(log(outcome) > log_prob){
    X_new <- X_old
    q1_new <- q1_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_q1_new <- exp(log(sd_q1_old)+n^(-0.7)*(acceptance_q1-target_accept))
  if(sd_q1_new>exp(50)) {sd_q1_new <- exp(50)}
  if(sd_q1_new<exp(-50)) {sd_q1_new <- exp(-50)}


  return(list(q1_new=q1_new, accept=accept, sd_q1_new=sd_q1_new))

}

# ------------------- Update u01 ----------------------
u01_log_prob <- function(u01, u_obs, s_obs, mu_u, mu_s, sigma_u_2, sigma_s_2, u01_lower, u01_upper,
                         mu_alpha, sig_alpha, beta){

  lp <- sum(dtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))+
    sum(dtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
    dunif(x = u01, min = u01_lower, max = u01_upper, log = TRUE)-
    pnorm(q = beta*u01, mean = mu_alpha, sd = sig_alpha, lower.tail = FALSE, log.p = TRUE) #from prior for alpha

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))+
      sum(myldtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))+
      dunif(x = u01, min = u01_lower, max = u01_upper, log = TRUE)-
      pnorm(q = beta*u01, mean = mu_alpha, sd = sig_alpha, lower.tail = FALSE, log.p = TRUE)
  }

  return(lp)

}

u01_update <- function(u01, u_obs, s_obs, alpha1_1, t0, k, beta, gamma, q, sigma_u_2, sigma_s_2, tau, compute_mean,
                       iter_num, sd_u01, u01_lower, u01_upper, mu_alpha, sig_alpha, target_accept){


  # Define inputs
  u01_upper_ <- min(c(u01_upper, alpha1_1/beta))
  u01_old <- u01; X_old <- log(u01_old-u01_lower)-log(u01_upper_-u01_old)
  sd_u01_old <- sd_u01

  n <- iter_num

  # Compute current (old) mean s and u for each cell
  mus_old <- compute_mean(tau, t0, k, alpha1_1, beta, gamma, q, u01=u01_old)

  # Apply AMH
  X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(sd_u01_old))

  # Transform back to q
  u01_new <- u01_upper_+(u01_lower-u01_upper_)/(1+exp(X_new))

  # Compute new mean s and u for each cell
  mus_new <- compute_mean(tau, t0, k, alpha1_1, beta, gamma, q, u01=u01_new)

  # Compute log acceptance probability
  log_acceptance <- u01_log_prob(u01 = u01_new, u_obs = u_obs, s_obs = s_obs, mu_u = mus_new$u_mean, mu_s = mus_new$s_mean,
                                 sigma_u_2 = sigma_u_2, sigma_s_2 = sigma_s_2, u01_lower = u01_lower, u01_upper = u01_upper,
                                 mu_alpha = mu_alpha, sig_alpha = sig_alpha, beta = beta)-
    u01_log_prob(u01 = u01_old, u_obs, s_obs, mu_u = mus_old$u_mean, mu_s = mus_old$s_mean, sigma_u_2, sigma_s_2,
                 u01_lower, u01_upper, mu_alpha, sig_alpha, beta)+
    log(u01_new-u01_lower)+log(u01_upper_-u01_new)-
    log(u01_old-u01_lower)-log(u01_upper_-u01_old)

  acceptance_u01 <- exp(log_acceptance)
  acceptance_u01 <- min(1, acceptance_u01)
  log_prob <- min(0,log_acceptance)

  # Decision
  outcome <- runif(1,min=0,max=1)
  if(log(outcome) > log_prob){
    X_new <- X_old
    u01_new <- u01_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_u01_new <- exp(log(sd_u01_old)+n^(-0.7)*(acceptance_u01-target_accept))
  if(sd_u01_new>exp(50)) {sd_u01_new <- exp(50)}
  if(sd_u01_new<exp(-50)) {sd_u01_new <- exp(-50)}


  return(list(u01_new=u01_new, accept=accept, sd_u01_new=sd_u01_new))

}


# ------------------- Update variance parameters ----------------------

# ------------------ sigma_u_2 -------------
sigma_u_2_log_prob <- function(sigma_u_2, u_obs, mu_u, sigma_u_2_hat, invgamma_f){

  lp <- sum(dtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))-
    (invgamma_f+1)*log(sigma_u_2)-invgamma_f*sigma_u_2_hat/sigma_u_2

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtnorm(x = u_obs, mean = mu_u, sd = sqrt(sigma_u_2), lower = 0, log = TRUE))-
      (invgamma_f+1)*log(sigma_u_2)-invgamma_f*sigma_u_2_hat/sigma_u_2
  }

  return(lp)
}

sigma_u_2_update <- function(sigma_u_2, u_obs, mu_u_fixed, iter_num,
                             sd_sigma_u_2, sigma_u_2_hat, invgamma_f, target_accept){
  # Define inputs
  sigma_u_2_old <- sigma_u_2; X_old <- log(sigma_u_2_old)
  sd_sigma_u_2_old <- sd_sigma_u_2

  n <- iter_num

  # Apply AMH
  X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(sd_sigma_u_2_old))

  # Transform back to sigma_u_2
  sigma_u_2_new <- exp(X_new)

  # Compute log acceptance probability
  log_acceptance <- sigma_u_2_log_prob(sigma_u_2 = sigma_u_2_new, u_obs = u_obs, mu_u = mu_u_fixed,
                                       sigma_u_2_hat = sigma_u_2_hat, invgamma_f = invgamma_f) -
    sigma_u_2_log_prob(sigma_u_2_old, u_obs, mu_u_fixed, sigma_u_2_hat, invgamma_f)+
    log(sigma_u_2_new) - log(sigma_u_2_old)

  acceptance_sigma_u <- exp(log_acceptance)
  acceptance_sigma_u <- min(1, acceptance_sigma_u)
  log_prob <- min(0,log_acceptance)

  # Decision
  outcome <- runif(1,min=0,max=1)
  if(log(outcome) > log_prob){
    X_new <- X_old
    sigma_u_2_new <- sigma_u_2_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_sigma_u_2_new <- exp(log(sd_sigma_u_2_old)+n^(-0.7)*(acceptance_sigma_u-target_accept))
  if(sd_sigma_u_2_new>exp(50)) {sd_sigma_u_2_new <- exp(50)}
  if(sd_sigma_u_2_new<exp(-50)) {sd_sigma_u_2_new <- exp(-50)}

  return(list(sigma_u_2_new=sigma_u_2_new, accept=accept, sd_sigma_u_2_new=sd_sigma_u_2_new))
}

# ------------------ sigma_s_2 -------------
sigma_s_2_log_prob <- function(sigma_s_2, s_obs, mu_s, sigma_s_2_hat, invgamma_f){

  lp <- sum(dtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))-
    (invgamma_f+1)*log(sigma_s_2)-invgamma_f*sigma_s_2_hat/sigma_s_2

  if(is.infinite(lp) || is.na(lp)){
    lp <- sum(myldtnorm(x = s_obs, mean = mu_s, sd = sqrt(sigma_s_2), lower = 0, log = TRUE))-
      (invgamma_f+1)*log(sigma_s_2)-invgamma_f*sigma_s_2_hat/sigma_s_2
  }
  return(lp)
}

sigma_s_2_update <- function(sigma_s_2, s_obs, mu_s_fixed, iter_num,
                             sd_sigma_s_2, sigma_s_2_hat, invgamma_f, target_accept){
  # Define inputs
  sigma_s_2_old <- sigma_s_2; X_old <- log(sigma_s_2_old)
  sd_sigma_s_2_old <- sd_sigma_s_2

  n <- iter_num

  # Apply AMH
  X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(sd_sigma_s_2_old))

  # Transform back to sigma_u_2
  sigma_s_2_new <- exp(X_new)

  # Compute log acceptance probability
  log_acceptance <- sigma_s_2_log_prob(sigma_s_2 = sigma_s_2_new, s_obs = s_obs, mu_s = mu_s_fixed,
                                       sigma_s_2_hat = sigma_s_2_hat, invgamma_f = invgamma_f) -
    sigma_s_2_log_prob(sigma_s_2_old, s_obs, mu_s_fixed, sigma_s_2_hat, invgamma_f)+
    log(sigma_s_2_new) - log(sigma_s_2_old)

  acceptance_sigma_s <- exp(log_acceptance)
  acceptance_sigma_s <- min(1, acceptance_sigma_s)
  log_prob <- min(0,log_acceptance)

  # Decision
  outcome <- runif(1,min=0,max=1)
  if(log(outcome) > log_prob){
    X_new <- X_old
    sigma_s_2_new <- sigma_s_2_old
    accept <- 0
  }else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_sigma_s_2_new <- exp(log(sd_sigma_s_2_old)+n^(-0.7)*(acceptance_sigma_s-target_accept))
  if(sd_sigma_s_2_new>exp(50)) {sd_sigma_s_2_new <- exp(50)}
  if(sd_sigma_s_2_new<exp(-50)) {sd_sigma_s_2_new <- exp(-50)}

  return(list(sigma_s_2_new=sigma_s_2_new, accept=accept, sd_sigma_s_2_new=sd_sigma_s_2_new))
}

# ------------------- Update states k ----------------------
k_update <- function(u_obs, s_obs, alpha1_1, beta, gamma, q, sigma_u_2, sigma_s_2, tau, t0, u01, p,
                     mu_tau, var_tau, compute_mean){

  # number of cells
  N <- nrow(u_obs)

  # compute the mean (mu_u, mu_s) for each cell when k=1 and k=2, conditional on current params
  mean_by_k <- list(compute_mean(tau = tau, t0 = t0, k = rep(1,N), alpha1_1 = alpha1_1, beta = beta, gamma = gamma, q = q, u01 = u01),
                    compute_mean(tau = tau, t0 = t0, k = rep(2,N), alpha1_1 = alpha1_1, beta = beta, gamma = gamma, q = q, u01 = u01))

  tau_upper <- c(t0,Inf)

  prob <- c(p, 1-p)

  kk <- vapply(1:N, function(cc) {

    if(tau[cc]>t0){
      k_cc <- 2

      return(k_cc)
    }else{
      ## Set log probability
      LP <- vapply(1:2, function(j) {

        lp <- dtnorm(x = u_obs[cc], mean = mean_by_k[[j]]$u_mean[cc], sd = sqrt(sigma_u_2), lower = 0, log = TRUE)+
          dtnorm(x = s_obs[cc], mean = mean_by_k[[j]]$s_mean[cc], sd = sqrt(sigma_s_2), lower = 0, log = TRUE)+
          myldtgamma(x = tau[cc], alpha = mu_tau[j]^2/var_tau[j], beta = mu_tau[j]/var_tau[j], lower = 0,
                     upper = tau_upper[j], log = TRUE)+
          log(prob[j])

        if(is.infinite(lp) || is.na(lp)){
          lp <- myldtnorm(x = u_obs[cc], mean = mean_by_k[[j]]$u_mean[cc], sd = sqrt(sigma_u_2), lower = 0, log = TRUE)+
            myldtnorm(x = s_obs[cc], mean = mean_by_k[[j]]$s_mean[cc], sd = sqrt(sigma_s_2), lower = 0, log = TRUE)+
            myldtgamma(x = tau[cc], alpha = mu_tau[j]^2/var_tau[j], beta = mu_tau[j]/var_tau[j], lower = 0,
                       upper = tau_upper[j], log = TRUE)+
            log(prob[j])
        }

        return(lp)
      }, FUN.VALUE = numeric(1))

      # Compute the normalizing constant
      nc <- -max(LP)
      P <- exp(LP+nc)/sum(exp(LP+nc)) # LP is a vector of length 2
      k_cc <- sample(1:2, 1, prob=P)

      return(k_cc)
    }

  }, FUN.VALUE = numeric(1))

  return(kk)
}

# ------------------- Update probability p ----------------------
p_update <- function(k){

  N1 <- sum(k==1)
  N2 <- sum(k==2)

  p_new <- rbeta(1, N1+1, N2+1)

  return(p_new)
}




