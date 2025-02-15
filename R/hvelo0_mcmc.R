hvelo0_mcmc <- function(u_obs, s_obs, niter, burn_in=0, thinning=1,
                        mu_alpha, sig_alpha=1, mu_gamma, sig_gamma=0.4,
                        mu_star_base, sig_star_base, eta_base, nu_base, mu_0, var_0,
                        q_lower=0.3, q_upper=6, sigma_u_2_hat, sigma_s_2_hat, invgamma_f=1,
                        u01_lower, u01_upper, target_accept=0.44, verbose=TRUE,
                        k_initial=NULL, p_initial=NULL, alpha1_1_initial, gamma_initial,
                        tau_initial, mu_tau_initial, var_tau_initial, t0_initial,
                        q_initial, sigma_u_2_initial, sigma_s_2_initial, u01_initial,
                        k_fix=NULL){

  # we fix beta to be 1, and estimate q = lambda/beta
  # in the post-processing step, q does not need to be adjusted by beta
  # function to compute the mean for both s and u
  compute_mean <- function(tau, t0, k, alpha1_1, beta, gamma, q, u01){

    u_tilde <- exp(-beta*tau)
    u0_tilde <- exp(-beta*t0)
    q <- rep(q,2)

    # Define function to compute u(u_tilde), s(u_tilde)
    u_func <- function(u_tilde, u0, k, alpha1, alpha2, q, beta){
      u <- u0[k]*u_tilde+alpha1[k]/beta*(1-u_tilde)+alpha2[k]/(beta-beta*q[k])*(u_tilde-u_tilde^q[k])
      return(u)
    }

    s_func <- function(u_tilde, u0, s0, k, alpha1, alpha2, q, beta, gamma){
      s <- alpha1[k]/gamma+1/(gamma-beta)*(beta*u0[k]-alpha1[k]+alpha2[k]/(1-q[k]))*u_tilde+
        (s0[k]-beta/(gamma-beta)*(u0[k]-alpha1[k]/gamma+alpha2[k]/(gamma-beta*q[k])))*u_tilde^(gamma/beta)-
        alpha2[k]/(1-q[k])/(gamma-beta*q[k])*u_tilde^q[k]
      return(s)
    }

    # alpha
    alpha1 <- c(alpha1_1,beta*u01)
    alpha2 <- c(alpha1_1-beta*u01,(alpha1_1-beta*u01)*(u0_tilde^q[1]-1)) #alpha1[2]-alpha1[1]+alpha2[1]*e^{-lambda*t0}

    # initial conditions
    u0 <- c()
    u0[1] <- u01
    u0[2] <- u_func(u_tilde = u0_tilde, u0, k = 1, alpha1, alpha2, q, beta)

    s0 <- c()
    s0[1] <- beta*u01/gamma
    s0[2] <- s_func(u_tilde = u0_tilde, u0, s0, k = 1, alpha1, alpha2, q, beta, gamma)

    # mean
    u_mean <- u_func(u_tilde = u_tilde, u0 = u0, k = k, alpha1 = alpha1, alpha2 = alpha2, q = q, beta = beta)
    s_mean <- s_func(u_tilde = u_tilde, u0 = u0, s0 = s0, k = k, alpha1 = alpha1, alpha2 = alpha2, q = q, beta = beta, gamma = gamma)

    return(list(u_mean=u_mean, s_mean=s_mean))
  }

  # number of cells
  N <- nrow(u_obs)

  # beta fixed to be 1
  beta <- 1

  # Stop if both k_initial and k_fix are not given
  if(is.null(k_initial)&is.null(k_fix)){
    stop('Either initial values or true values should be given for states k!')
  }

  if(is.null(p_initial)&is.null(k_fix)){
    stop('Initial value for p should be given when k is not fixed!')
  }


  #------------------------ Step 0:Prepare for outputs -----------------
  k_output <- NULL
  p_output <- c()

  alpha1_1_output <- c()
  gamma_output <- c()

  tau_output <- NULL
  t0_output <- c()

  q_output <- c()

  # initial condition u0[1]
  u01_output <- c()

  sigma_u_2_output <- c()
  sigma_s_2_output <- c()

  # prior params for tau
  mu_tau_output <- NULL
  var_tau_output <- NULL

  # save cell means
  mus_output <- NULL

  # scale params in AMH
  sd_alpha1_1_output <- c()
  sd_gamma_output <- c()
  sd_tau_output <- NULL
  sd_t0_output <- c()
  sd_q_output <- c()
  sd_sigma_u_2_output <- c()
  sd_sigma_s_2_output <- c()

  sd_mu_tau_output <- NULL
  sd_var_tau_output <- NULL

  sd_u01_output <- c()

  # (Average) Acceptance probability
  acceptance_count_avg <- data.frame(alpha1_1_accept = rep(0,niter), gamma_accept = rep(0,niter),
                                     tau_accept = rep(0,niter), t0_accept = rep(0,niter),
                                     q_accept = rep(0,niter), u01_accept = rep(0, niter),
                                     sigma_u_2_accept = rep(0,niter), sigma_s_2_accept = rep(0,niter),
                                     mu_tau_1_accept = rep(0, niter), mu_tau_2_accept = rep(0, niter),
                                     var_tau_1_accept = rep(0, niter), var_tau_2_accept = rep(0, niter))

  # log-likelihood and posterior
  log_like <- rep(NA, niter)
  log_post <- rep(NA, niter)

  #----------------------- Step 1: Set the initial values as new values ----------------------------
  if(is.null(k_fix)){
    k_new <- k_initial
  }else{
    k_new <- k_fix
  }
  p_new <- p_initial

  alpha1_1_new <- alpha1_1_initial
  gamma_new <- gamma_initial

  tau_new <- tau_initial
  t0_new <- t0_initial

  q_new <- q_initial


  # initial condition u0[1]
  u01_new <- u01_initial

  sigma_u_2_new <- sigma_u_2_initial
  sigma_s_2_new <- sigma_s_2_initial

  # prior params for tau
  mu_tau_new <- mu_tau_initial; var_tau_new <- var_tau_initial

  # Total count of acceptance
  alpha1_1_count <- 0; gamma_count <- 0; tau_count <- 0; t0_count <- 0; q_count <- 0
  u01_count <- 0
  sigma_u_2_count <- 0; sigma_s_2_count <- 0; mu_tau_1_count <- 0; mu_tau_2_count <- 0; var_tau_1_count <- 0; var_tau_2_count <- 0

  if(max(tau_new[k_new==1])>=t0_new){
    stop('initial tau in k=1 greater than t0!!')
  }
  #----------------------- Step 2: Prepare for the adaptive covariance update -----------------------------

  # 1) For alpha1_1
  sd_alpha1_1_new <- 0.1

  # 2) For gamma
  sd_gamma_new <- 0.1

  # 3) For tau
  sd_tau_new <- rep(0.1,N)

  # 4) For t0
  sd_t0_new <- 0.1

  # 5) For q
  sd_q_new <- 0.1

  # 6) For sigma_u_2
  sd_sigma_u_2_new <- 0.1

  # 7) For sigma_s_2
  sd_sigma_s_2_new <- 0.1

  # 8) For parior params of tau
  sd_mu_tau_new <- rep(0.1,2)

  sd_var_tau_new <- rep(0.1,2)

  # 9) For initial condition u0[1]
  sd_u01_new <- 0.1

  #----------------------- Step 3: Updates -----------------------------
  output_index <- 0

  # count the number of iterations where MH is used to draw mu_tau and var_tau
  MH_count_mu_tau <- rep(2,2); MH_count_var_tau <- rep(2,2)

  start_time_mcmc <- Sys.time()
  if(verbose){
    # Show time
    print(paste('Start MCMC:', start_time_mcmc))
    cat('\n')

    pb <- txtProgressBar(min = 1, max = niter+1, style = 3)
  }

  # Iteration starts with iter = 2
  for(iter in 2:(niter+1)){

    if(verbose) {setTxtProgressBar(pb, iter)}
    # Starting value of the output index = 1
    # If the current iteration is greater than the burn in and divisible by the thinning index
    criterion <- (iter-1 > burn_in & (iter-1-burn_in)%%thinning == 0)

    if(criterion){
      output_index <- output_index + 1
      update <- TRUE
    }else{
      update <- FALSE
    }

    # 1) Update tau
    tau_output_sim <- tau_update(tau = tau_new, u_obs = u_obs, s_obs = s_obs, k = k_new, alpha1_1 = alpha1_1_new,
                                 beta = beta, gamma = gamma_new, q = q_new, sigma_u_2 = sigma_u_2_new, sigma_s_2 = sigma_s_2_new,
                                 t0 = t0_new, u01 = u01_new, compute_mean = compute_mean, iter_num = iter,
                                 sd_tau = sd_tau_new, mu_tau = mu_tau_new, var_tau = var_tau_new,
                                 target_accept)

    tau_new <- tau_output_sim$tau_new
    sd_tau_new <- tau_output_sim$sd_tau_new

    tau_count <- tau_count+tau_output_sim$accept_count
    acceptance_count_avg$tau_accept[iter-1] <- tau_count/((iter-1)*N)



    # 1-1) Update prior params for tau
    # a. mu_tau
    mu_tau_output_sim <- mu_tau_update(mu_tau = mu_tau_new, tau = tau_new, var_tau = var_tau_new, k = k_new, mu_star_base = mu_star_base,
                                       sig_star_base = sig_star_base, t0 = t0_new, iter_num = iter, sd_mu_tau = sd_mu_tau_new,
                                       target_accept = target_accept, MH_count = MH_count_mu_tau)

    mu_tau_new <- mu_tau_output_sim$mu_tau_new
    sd_mu_tau_new <- mu_tau_output_sim$sd_mu_tau_new

    mu_tau_1_count <- mu_tau_1_count+mu_tau_output_sim$accept_count[1]
    acceptance_count_avg$mu_tau_1_accept[iter-1] <- mu_tau_1_count/(iter-1)

    mu_tau_2_count <- mu_tau_2_count+mu_tau_output_sim$accept_count[2]
    acceptance_count_avg$mu_tau_2_accept[iter-1] <- mu_tau_2_count/(iter-1)

    MH_count_mu_tau <- MH_count_mu_tau+mu_tau_output_sim$MH_used

    # b. var_tau
    var_tau_output_sim <- var_tau_update(var_tau = var_tau_new, tau = tau_new, mu_tau = mu_tau_new, k = k_new, eta_base = eta_base,
                                         nu_base = nu_base, t0 = t0_new, iter_num = iter, sd_var_tau = sd_var_tau_new,
                                         target_accept = target_accept, MH_count = MH_count_var_tau)

    var_tau_new <- var_tau_output_sim$var_tau_new
    sd_var_tau_new <- var_tau_output_sim$sd_var_tau_new

    var_tau_1_count <- var_tau_1_count+var_tau_output_sim$accept_count[1]
    acceptance_count_avg$var_tau_1_accept[iter-1] <- var_tau_1_count/(iter-1)

    var_tau_2_count <- var_tau_2_count+var_tau_output_sim$accept_count[2]
    acceptance_count_avg$var_tau_2_accept[iter-1] <- var_tau_2_count/(iter-1)

    MH_count_var_tau <- MH_count_var_tau+var_tau_output_sim$MH_used



    # 2) Update t0
    t0_output_sim <- t0_update(t0 = t0_new, u_obs = u_obs, s_obs = s_obs, k = k_new, alpha1_1 = alpha1_1_new,
                               beta = beta, gamma = gamma_new, q = q_new, sigma_u_2 = sigma_u_2_new, sigma_s_2 = sigma_s_2_new,
                               tau = tau_new, u01 = u01_new, compute_mean = compute_mean, iter_num = iter,
                               sd_t0 = sd_t0_new, mu_0 = mu_0, var_0 = var_0, mu_tau = mu_tau_new, var_tau = var_tau_new, target_accept)

    t0_new <- t0_output_sim$t0_new
    sd_t0_new <- t0_output_sim$sd_t0_new

    t0_count <- t0_count+t0_output_sim$accept
    acceptance_count_avg$t0_accept[iter-1] <- t0_count/(iter-1)



    # 3) Update rate alpha1_1
    alpha1_1_output_sim <- alpha1_1_update(alpha1_1 = alpha1_1_new, u_obs = u_obs, s_obs = s_obs, t0 = t0_new, k = k_new,
                                           beta = beta, gamma = gamma_new, q = q_new, sigma_u_2 = sigma_u_2_new, sigma_s_2 = sigma_s_2_new,
                                           tau = tau_new, u01 = u01_new, compute_mean = compute_mean, iter_num = iter,
                                           sd_alpha1_1 = sd_alpha1_1_new, mu_alpha = mu_alpha, sig_alpha = sig_alpha,
                                           target_accept)

    alpha1_1_new <- alpha1_1_output_sim$alpha1_1_new
    sd_alpha1_1_new <- alpha1_1_output_sim$sd_alpha1_1_new

    alpha1_1_count <- alpha1_1_count+alpha1_1_output_sim$accept
    acceptance_count_avg$alpha1_1_accept[iter-1] <- alpha1_1_count/(iter-1)



    # 4) Update rate gamma
    gamma_output_sim <- gamma_update(gamma = gamma_new, s_obs = s_obs, t0 = t0_new, k = k_new, alpha1_1 = alpha1_1_new,
                                     beta = beta, q = q_new, sigma_s_2 = sigma_s_2_new, tau = tau_new, u01 = u01_new,
                                     compute_mean = compute_mean, iter_num = iter,
                                     sd_gamma = sd_gamma_new, mu_gamma = mu_gamma, sig_gamma = sig_gamma,
                                     target_accept)

    gamma_new <- gamma_output_sim$gamma_new
    sd_gamma_new <- gamma_output_sim$sd_gamma_new

    gamma_count <- gamma_count+gamma_output_sim$accept
    acceptance_count_avg$gamma_accept[iter-1] <- gamma_count/(iter-1)


    # 5) Update q (lambda/beta)
    q1_output_sim <- q1_update(q = q_new, u_obs = u_obs, s_obs = s_obs, alpha1_1 = alpha1_1_new, t0 = t0_new, k = k_new,
                               beta = beta, gamma = gamma_new, sigma_u_2 = sigma_u_2_new, sigma_s_2 = sigma_s_2_new,
                               tau = tau_new, u01 = u01_new, compute_mean = compute_mean, iter_num = iter,
                               sd_q1 = sd_q_new, q_lower = q_lower, q_upper = q_upper, target_accept)

    q_new <- q1_output_sim$q1_new
    sd_q_new <- q1_output_sim$sd_q1_new

    q_count <- q_count+q1_output_sim$accept
    acceptance_count_avg$q_accept[iter-1] <- q_count/(iter-1)



    # 6) Update initial condition u01
    u01_output_sim <- u01_update(u01 = u01_new, u_obs = u_obs, s_obs = s_obs, alpha1_1 = alpha1_1_new, t0 = t0_new, k = k_new,
                                 beta = beta, gamma = gamma_new, q = q_new, sigma_u_2 = sigma_u_2_new, sigma_s_2 = sigma_s_2_new,
                                 tau = tau_new, compute_mean = compute_mean, iter_num = iter,
                                 sd_u01 = sd_u01_new, u01_lower = u01_lower, u01_upper = u01_upper,
                                 mu_alpha = mu_alpha, sig_alpha = sig_alpha, target_accept)

    u01_new <- u01_output_sim$u01_new
    sd_u01_new <- u01_output_sim$sd_u01_new

    u01_count <- u01_count+u01_output_sim$accept
    acceptance_count_avg$u01_accept[iter-1] <- u01_count/(iter-1)



    # for updating variance parameters, the means for s and u are fixed
    mus_new <- compute_mean(tau_new, t0_new, k_new, alpha1_1_new, beta, gamma_new, q_new, u01_new)

    # 7) Update variance parameters
    # sigma_u_2
    sigma_u_2_output_sim <- sigma_u_2_update(sigma_u_2 = sigma_u_2_new, u_obs = u_obs, mu_u_fixed = mus_new$u_mean,
                                             iter_num = iter, sd_sigma_u_2 = sd_sigma_u_2_new,
                                             sigma_u_2_hat = sigma_u_2_hat, invgamma_f = invgamma_f, target_accept)

    sigma_u_2_new <- sigma_u_2_output_sim$sigma_u_2_new
    sd_sigma_u_2_new <- sigma_u_2_output_sim$sd_sigma_u_2_new

    sigma_u_2_count <- sigma_u_2_count+sigma_u_2_output_sim$accept
    acceptance_count_avg$sigma_u_2_accept[iter-1] <- sigma_u_2_count/(iter-1)


    # sigma_s_2
    sigma_s_2_output_sim <- sigma_s_2_update(sigma_s_2 = sigma_s_2_new, s_obs = s_obs, mu_s_fixed = mus_new$s_mean,
                                             iter_num = iter, sd_sigma_s_2 = sd_sigma_s_2_new,
                                             sigma_s_2_hat = sigma_s_2_hat, invgamma_f = invgamma_f, target_accept)

    sigma_s_2_new <- sigma_s_2_output_sim$sigma_s_2_new
    sd_sigma_s_2_new <- sigma_s_2_output_sim$sd_sigma_s_2_new

    sigma_s_2_count <- sigma_s_2_count+sigma_s_2_output_sim$accept
    acceptance_count_avg$sigma_s_2_accept[iter-1] <- sigma_s_2_count/(iter-1)



    # 8) Update states k
    if(is.null(k_fix)){
      k_new <- k_update(u_obs = u_obs, s_obs = s_obs, alpha1_1 = alpha1_1_new, beta = beta, gamma = gamma_new, q = q_new,
                        sigma_u_2 = sigma_u_2_new, sigma_s_2 = sigma_s_2_new, tau = tau_new, t0 = t0_new, u01 = u01_new,
                        p = p_new, mu_tau = mu_tau_new, var_tau = var_tau_new, compute_mean = compute_mean)
    }

    # 9) Update p
    if(is.null(k_fix)){
      p_new <- p_update(k = k_new)
    }


    # compute cell means
    mus_new <- compute_mean(tau_new, t0_new, k_new, alpha1_1_new, beta, gamma_new, q_new, u01_new)

    # compute log-likelihood and posterior
    log_like[iter-1] <- sum(dtnorm(x = u_obs, mean = mus_new$u_mean, sd = sqrt(sigma_u_2_new), lower = 0, log = TRUE))+
      sum(dtnorm(x = s_obs, mean = mus_new$s_mean, sd = sqrt(sigma_s_2_new), lower = 0, log = TRUE))
    if(is.infinite(log_like[iter-1]) || is.na(log_like[iter-1])){
      log_like[iter-1] <- sum(myldtnorm(x = u_obs, mean = mus_new$u_mean, sd = sqrt(sigma_u_2_new), lower = 0, log = TRUE))+
        sum(myldtnorm(x = s_obs, mean = mus_new$s_mean, sd = sqrt(sigma_s_2_new), lower = 0, log = TRUE))
    }

    log_prior <- dunif(x = u01_new, min = u01_lower, max = u01_upper, log = TRUE)+
      dtnorm(x = alpha1_1_new, mean = mu_alpha, sd = sig_alpha, lower = beta*u01_new, log = TRUE)+
      dtnorm(x = gamma_new, mean = mu_gamma, sd = sig_gamma, lower = 0, log = TRUE)+
      dunif(x = q_new, min = q_lower, max = q_upper, log = TRUE)+
      dgamma(x = t0_new, shape = mu_0^2/var_0, rate = mu_0/var_0, log = TRUE)+
      sum(myldtgamma(x = tau_new, alpha = mu_tau_new[k_new]^2/var_tau_new[k_new], beta = mu_tau_new[k_new]/var_tau_new[k_new], lower = 0,
                     upper = ifelse(k_new==1, t0_new, Inf), log = TRUE))+
      sum(dtnorm(x = mu_tau_new, mean = mu_star_base, sd = sig_star_base, lower = 0, log = TRUE))+
      sum(dtnorm(x = var_tau_new, mean = eta_base, sd = nu_base, lower = 0, log = TRUE))+
      ifelse(is.null(k_fix),sum(k_new==1)*log(p_new) + sum(k_new==2)*log(1-p_new), 0)-
      (invgamma_f+1)*log(sigma_u_2_new)-invgamma_f*sigma_u_2_hat/sigma_u_2_new-
      (invgamma_f+1)*log(sigma_s_2_new)-invgamma_f*sigma_s_2_hat/sigma_s_2_new

    if(is.infinite(log_prior) || is.na(log_prior)){
      log_prior <- dunif(x = u01_new, min = u01_lower, max = u01_upper, log = TRUE)+
        myldtnorm(x = alpha1_1_new, mean = mu_alpha, sd = sig_alpha, lower = beta*u01_new, log = TRUE)+
        myldtnorm(x = gamma_new, mean = mu_gamma, sd = sig_gamma, lower = 0, log = TRUE)+
        dunif(x = q_new, min = q_lower, max = q_upper, log = TRUE)+
        dgamma(x = t0_new, shape = mu_0^2/var_0, rate = mu_0/var_0, log = TRUE)+
        sum(myldtgamma(x = tau_new, alpha = mu_tau_new[k_new]^2/var_tau_new[k_new], beta = mu_tau_new[k_new]/var_tau_new[k_new], lower = 0,
                       upper = ifelse(k_new==1, t0_new, Inf), log = TRUE))+
        sum(myldtnorm(x = mu_tau_new, mean = mu_star_base, sd = sig_star_base, lower = 0, log = TRUE))+
        sum(myldtnorm(x = var_tau_new, mean = eta_base, sd = nu_base, lower = 0, log = TRUE))+
        ifelse(is.null(k_fix),sum(k_new==1)*log(p_new) + sum(k_new==2)*log(1-p_new), 0)-
        (invgamma_f+1)*log(sigma_u_2_new)-invgamma_f*sigma_u_2_hat/sigma_u_2_new-
        (invgamma_f+1)*log(sigma_s_2_new)-invgamma_f*sigma_s_2_hat/sigma_s_2_new
    }

    log_post[iter-1] <- log_like[iter-1]+log_prior
    #-------------------------- Step 4: Update simulated values ------------------------
    if(update == TRUE){
      k_output[[output_index]] <- k_new
      p_output[output_index] <- p_new

      alpha1_1_output[output_index] <- alpha1_1_new
      gamma_output[output_index] <- gamma_new

      tau_output[[output_index]] <- tau_new
      t0_output[output_index] <- t0_new

      q_output[output_index] <- q_new

      # initial condition u0[1]
      u01_output[output_index] <- u01_new

      sigma_u_2_output[output_index] <- sigma_u_2_new
      sigma_s_2_output[output_index] <- sigma_s_2_new

      # prior params for tau
      mu_tau_output[[output_index]] <- mu_tau_new
      var_tau_output[[output_index]] <- var_tau_new

      # save cell means
      mus_output[[output_index]] <- mus_new

      # scale params in AMH
      sd_alpha1_1_output[output_index] <- sd_alpha1_1_new
      sd_gamma_output[output_index] <- sd_gamma_new
      sd_tau_output[[output_index]] <- sd_tau_new
      sd_t0_output[output_index] <- sd_t0_new
      sd_q_output[output_index] <- sd_q_new
      sd_sigma_u_2_output[output_index] <- sd_sigma_u_2_new
      sd_sigma_s_2_output[output_index] <- sd_sigma_s_2_new
      sd_mu_tau_output[[output_index]] <- sd_mu_tau_new
      sd_var_tau_output[[output_index]] <- sd_var_tau_new
      sd_u01_output[output_index] <- sd_u01_new
    }


  }#end for loop

  # Return the list
  my_list <- list('k_output' = k_output, 'p_output' = p_output, 'alpha1_1_output' = alpha1_1_output, 'gamma_output' = gamma_output,
                  'tau_output' = tau_output, 't0_output' = t0_output, 'lambda_output' = q_output, 'u01_output' = u01_output,
                  'sigma_u_2_output' = sigma_u_2_output, 'sigma_s_2_output' = sigma_s_2_output,
                  'mu_tau_output' = mu_tau_output, 'var_tau_output' = var_tau_output, 'mus_output' = mus_output,
                  'acceptance_count_avg' = acceptance_count_avg, 'output_index' = output_index,
                  'log_like' = log_like, 'log_post' = log_post,
                  'mu_alpha' = mu_alpha, 'sig_alpha' = sig_alpha, 'mu_gamma' = mu_gamma, 'sig_gamma' = sig_gamma,
                  'mu_star_base' = mu_star_base, 'sig_star_base' = sig_star_base, 'eta_base' = eta_base, 'nu_base' = nu_base,
                  'mu_0' = mu_0, 'var_0' = var_0, 'lambda_lower' = q_lower, 'lambda_upper' = q_upper,
                  'u01_lower' = u01_lower, 'u01_upper' = u01_upper,
                  'sigma_u_2_hat' = sigma_u_2_hat, 'sigma_s_2_hat' = sigma_s_2_hat, 'invgamma_f' = invgamma_f)

  end_time <- Sys.time()
  diff_time <- difftime(end_time,start_time_mcmc)

  if(verbose) {
    close(pb)

    cat('\n')
    print(paste('End:',end_time))
    cat('\n')
    print(paste('MCMC running time:', round(diff_time, digits = 3),units(diff_time)))

  }

  return(my_list)
}


