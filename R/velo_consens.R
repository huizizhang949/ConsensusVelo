#' Implement consensus velocity for inference of a single gene
#'
#' @description
#' The function runs multiple chains in parallel to infer unknown parameters for a single gene,
#' including a preparation step to find good initial values where states are fixed to empirical values.
#'
#' @import msm
#' @import graphics
#' @importFrom stats dgamma pnorm pgamma rnorm runif dunif var prcomp density rbeta rgamma update
#' @import ggplot2
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param k.inits a matrix of initial values for the state. Each column corresponds to one set of initial values, with rows corresponding to cells.
#' Can use the output from \code{generate_k}.
#' @param empirical output from \code{get_empirical}.
#' @param epsilon a small positive value used to compute empirical \eqn{\tau}. Better to use the same value as used in \code{get_empirical}.
#' @param alpha1.sd standard deviation in the truncated normal prior for \eqn{\alpha_1^{(1)}}.
#' @param gamma.sd standard deviation in the truncated normal prior for \eqn{\gamma}.
#' @param t0.sd standard deviation in the Gamma prior for the switching point \eqn{t_0^{(2)}}.
#' @param lambda.lower,lambda.upper the lower and upper bound in the Uniform prior for \eqn{\lambda}.
#' @param seed a vector of values used as seed for each chain. Must be of the same length as the number of chains.
#' @param prep_niter,comp_niter integer. Total number of MCMC iterations for the preparation step and complete step.
#' @param prep_burn_in,comp_burn_in the length of burn-in period during MCMC for the preparation step and complete step.
#' @param prep_thinning,comp_thinning the thinning applied after burn-in period for the preparation step and complete step.
#' @param n_cores number of cores used for parallel computing.
#' @param n_chains number of chains to run.
#'
#' @usage velo_consens(u.obs, s.obs, k.inits, empirical, epsilon=1e-3, alpha1.sd, gamma.sd, t0.sd,
#'    lambda.lower=0, lambda.upper=8, seed, prep_niter, prep_burn_in=0, prep_thinning=1,
#'    comp_niter, comp_burn_in=0, comp_thinning=1, n_cores=30, n_chains=100)
#'
#' @return The output is a list of results for all chains. Within each list, it contains the following components for each chain:
#' \item{output_index}{total number of saved MCMC samples, taking into account of burn-in and thinning.}
#' \item{k_output}{a list of vectors for sampled states. The list is of length \code{output_index}.}
#' \item{p_output}{a vector of sampled hyper-parameters \eqn{p}. The vector is of length \code{output_index}.}
#' \item{alpha1_1_output}{a vector of sampled rates \eqn{\alpha_1^{(1)}}. The vector is of length \code{output_index}.}
#' \item{gamma_output}{samples for rates \eqn{\gamma}. Similar to \code{alpha1_1_output}.}
#' \item{tau_output}{samples for local time \eqn{\tau}. Similar to \code{k_output}.}
#' \item{t0_output}{samples for switching point \eqn{t_0^{(2)}}.}
#' \item{lambda_output}{samples for \eqn{\lambda}.}
#' \item{u01_output}{samples for initial condition \eqn{u_0^{(1)}}.}
#' \item{sigma_u_2_output, sigma_s_2_output}{samples for variance parameters for unspliced and spliced counts.}
#' \item{mu_tau_output}{a list of vectors for sampled hyper-parameters \eqn{\mu_{\tau,j}} for \eqn{j=1,2}. The list is of length \code{output_index}.}
#' \item{var_tau_output}{samples for hyper-parameters \eqn{\sigma^2_{\tau,j}} for \eqn{j=1,2}. Similar to \code{mu_tau_output}.}
#' \item{mus_output}{a list of length \code{output_index}. Each element of the list is a list of two vectors saving samples of \eqn{u(t), s(t)} from the ODE.}
#' \item{acceptance_count_avg}{average acceptance probabilities over iterations for parameters drawn using adaptive Metroplis-Hastings.}
#' \item{log_like, log_post}{log-likelihood and log-posterior at each iteration (ignoring burnin and thinning).}
#' \item{mu_alpha, sig_alpha}{prior mean and standard deviation for \eqn{\alpha_1^{(1)}}.}
#' \item{mu_gamma, sig_gamma}{prior mean and standard deviation for \eqn{\gamma}.}
#' \item{mu_star_base, sig_star_base}{prior mean and standard deviation for \eqn{\mu_{\tau,j}} for \eqn{j=1,2}.}
#' \item{eta_base, nu_base}{prior mean and standard deviation for \eqn{\sigma^2_{\tau,j}} for \eqn{j=1,2}.}
#' \item{mu_0, var_0}{prior mean and variance for \eqn{t_0^{(2)}}.}
#' \item{lambda_lower, lambda_upper}{the lower and upper bound in the Uniform prior for \eqn{\lambda}.}
#' \item{u01_lower, u01_upper}{the lower and upper bound in the Uniform prior for \eqn{u_0^{(1)}}.}
#' \item{sigma_u_2_hat, sigma_s_2_hat, invgamma_f}{prior parameters in the inverse-gamma prior for variance parameters.}
#'
#' @export
#'
#' @examples
#' set.seed(4343)
#' seed <- sample(1:1e5,size=Width)
#' consensus_result <- velo_consens(u.obs = u.obs, s.obs = s.obs,
#'   k.inits = k.inits, empirical = empirical, epsilon = 1e-3, alpha1.sd = 5,
#'   gamma.sd = 0.4, t0.sd = 3, lambda.lower = 0.3, lambda.upper = 8,
#'   seed = seed, prep_niter = 10, prep_burn_in = 0, prep_thinning = 1,
#'   comp_niter = Depth, comp_burn_in = 0, comp_thinning = 1, n_cores = 3, n_chains = Width)
velo_consens <- function(u.obs, s.obs, k.inits, empirical, epsilon=1e-3,
                         alpha1.sd, gamma.sd, t0.sd, lambda.lower=0, lambda.upper=8,
                         seed, prep_niter, prep_burn_in=0, prep_thinning=1,
                         comp_niter, comp_burn_in=0, comp_thinning=1, n_cores=30, n_chains=100){

  show_end <- function(i){
    print(paste('preparation for chain',i,'completed!'))
  }

  result <- pbapply::pblapply(1:n_chains, function(i){
  # result <- lapply(1:n_chains, function(i){

    # -------------- preparation step ------------

    k.random <- k.inits[,i]

    # alpha and gamma
    alpha1.hat <- empirical$params$alpha1.hat; gamma.hat <- empirical$params$gamma.hat

    # hyper parameters for mu_tau, var_tau
    mm <- tapply(empirical$tau.hat, empirical$k.hat, mean); vv <- tapply(empirical$tau.hat, empirical$k.hat, var)

    # t0
    m_t0 <- empirical$params$t0.hat; var_0 <- t0.sd^2

    # lambda
    q_lower <- lambda.lower; q_upper <- lambda.upper

    # variance
    sigma.hat.u <- empirical$params$sigma.hat.u; sigma.hat.s <- empirical$params$sigma.hat.s

    # u0
    u0min1 <- empirical$params$u0min1; u0min2 <- empirical$params$u0min2
    u0max1 <- empirical$params$u0max1; u0max2 <- empirical$params$u0max2

    # initial tau from random k
    if(empirical$alpha.type=='avg'){
      uu.random <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random])/(c(mean(c(u0min1,u0min2)),alpha1.hat)[k.random]-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random])
    }else{
      uu.random <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random])/(c(mean(c(u0min1,u0min2)),mean(c(u0max1,u0max2)))[k.random]-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random])
    }
    uu.random[uu.random<=0] <- epsilon
    uu.random[uu.random>=1] <- 1-epsilon

    tau.random <- -log(uu.random)

    # initial lambda
    set.seed(seed[i])
    q_initial <- runif(1, 0.5, 2)
    while(q_initial==1 || q_initial==gamma.hat){
      q_initial <- runif(1, 0.5, 2)
    }

    chain <- hvelo0_mcmc(u_obs=matrix(u.obs,ncol=1), s_obs=matrix(s.obs,ncol=1), niter = prep_niter,
                         burn_in = prep_burn_in, thinning = prep_thinning,
                         mu_alpha = alpha1.hat, sig_alpha = alpha1.sd, mu_gamma = gamma.hat, sig_gamma = gamma.sd,
                         mu_star_base = mm, sig_star_base = 2*mm, eta_base = vv, nu_base = 2*vv,
                         mu_0 = m_t0, var_0 = var_0, q_lower = lambda.lower, q_upper = lambda.upper,
                         sigma_u_2_hat = sigma.hat.u^2, sigma_s_2_hat = sigma.hat.s^2, invgamma_f = 1,
                         u01_lower = u0min1, u01_upper = u0min2,
                         target_accept = 0.44,
                         k_true = k.random, alpha1_1_initial = alpha1.hat, gamma_initial = gamma.hat, tau_initial = tau.random,
                         mu_tau_initial = mm, var_tau_initial = vv,
                         t0_initial = max(c(max(tau.random[k.random==1]),m_t0))+0.1, q_initial = q_initial,
                         sigma_u_2_initial = sigma.hat.u^2, sigma_s_2_initial = sigma.hat.s^2, u01_initial = mean(c(u0min1, u0min2)))

    show_end(i)
    # -------------- complete step ---------------

    ind <- seq(1,chain$output_index,by=1)
    n <- length(u.obs)

    # get initials
    alpha1_1_initial <- mean(chain$alpha1_1_output[ind]); gamma_initial <- mean(chain$gamma_output[ind])
    q_initial <- mean(sapply(chain$lambda_output[ind], function(l) l[1]))
    t0_initial <- mean(chain$t0_output[ind]); u01_initial <- mean(chain$u01_output[ind])
    tau_initial <- sapply(1:n, function(cc) mean(sapply(chain$tau_output[ind], function(l) l[cc])))
    sigma_u_2_initial <- mean(chain$sigma_u_2_output[ind]);sigma_s_2_initial <- mean(chain$sigma_s_2_output[ind])

    mu_tau_initial <- c(mean(sapply(chain$mu_tau_output[ind], function(l) l[1])),
                        mean(sapply(chain$mu_tau_output[ind], function(l) l[2])))
    var_tau_initial <- c(mean(sapply(chain$var_tau_output[ind], function(l) l[1])),
                         mean(sapply(chain$var_tau_output[ind], function(l) l[2])))

    if(max(tau_initial[k.random==1])>=t0_initial){
      t0_initial <- max(tau_initial[k.random==1])+0.1
    }

    set.seed(seed[i]-1)
    chain2 <- hvelo0_mcmc(u_obs=matrix(u.obs,ncol=1), s_obs=matrix(s.obs,ncol=1), niter = comp_niter,
                          burn_in = comp_burn_in, thinning = comp_thinning,
                          mu_alpha = alpha1.hat, sig_alpha = alpha1.sd, mu_gamma = gamma.hat, sig_gamma = gamma.sd,
                          mu_star_base = mm, sig_star_base = 2*mm, eta_base = vv, nu_base = 2*vv,
                          mu_0 = m_t0, var_0 = var_0, q_lower = lambda.lower, q_upper = lambda.upper,
                          sigma_u_2_hat = sigma.hat.u^2, sigma_s_2_hat = sigma.hat.s^2, invgamma_f = 1,
                          u01_lower = u0min1, u01_upper = u0min2,
                          target_accept = 0.44,
                          k_initial = k.random, alpha1_1_initial = alpha1_1_initial, gamma_initial = gamma_initial, tau_initial = tau_initial,
                          p_initial = mean(k.random==1), mu_tau_initial = mu_tau_initial, var_tau_initial = var_tau_initial,
                          t0_initial = t0_initial, q_initial = q_initial, sigma_u_2_initial = sigma_u_2_initial, sigma_s_2_initial = sigma_s_2_initial,
                          u01_initial = u01_initial)


  # })
  }, cl = n_cores)

  return(result)

}


#' Plot absolute difference in entropy
#'
#' @description
#' This function plot absolute difference in entropy to select a suitable value for \eqn{W} (number of chains) in consensus velocity.
#'
#' @import ggplot2
#'
#' @param Ws a vector of candidate values for widths (number of chains).
#' @param Ds either an integer or a vector of integers. Absolute differences will be computed for each value in \code{Ds}.
#' @param mcmc_k a list, each corresponding to the samples of k from a single chain.
#'
#' @return If \code{Ds} is an integer, a line plot showing absolute difference in entropy between different widths.
#' If \code{Ds} is a vector, for each value in \code{Ds} a line will be plotted.
#'
#' @export
#'
#' @examples
#' Ws <- c(1, 2, 3)
#' plot_entropy(Ws = Ws, Ds = 500, mcmc_k = mcmc_k_result)
plot_entropy <- function(Ws, Ds, mcmc_k){

  # Entropy
  # ---- a list of length = no of different D, each is a matrix of max(W)*n_obs ----
  # depend on values in D and max(W)
  k_mat_list_by_d <- lapply(Ds, function(d) {
    t(sapply(mcmc_k, function(k_output) {
      k_result <- k_output[[d]]
    }))
  })

  # ---- a list of list, the first level is of length = n_W, the second level is of length = n_D ----
  # within the list of list is a vector recording p(induction) for every cell

  # depend on the values in D and W

  cms_list <- NULL
  for (i in 1:length(Ws)) {
    cms_list[[i]] <- lapply(1:length(Ds), function(j) {
      if(Ws[i]==1) {
        vec <-  ifelse(k_mat_list_by_d[[j]][1,]==1, 1, 0)
      }else{
        vec <-  apply(k_mat_list_by_d[[j]][1:Ws[i],], 2, function(x) mean(x==1))
      }
      return(vec)
    })
  }

  entropy_vec_w <- data.frame(tidyr::expand_grid(Ws,Ds,entropy=NA))
  entropy_vec_w$entropy <- apply(entropy_vec_w, 1, function(vec) {
    i <- which(Ws==vec[1])
    j <- which(Ds==vec[2])

    p_c <- cms_list[[i]][[j]]

    entropy <- sapply(p_c, function(p) {
      if(p==1 || p==0){
        val <- 0
      }else{
        val <- -(p*log(p)+(1-p)*log(1-p))
      }

      return(val)
    })

    return(sum(entropy))
  })

  entropy_vec_w$AD_entropy <-  apply(entropy_vec_w, 1, function(vec) {
    i <- which(Ws==vec[1])
    j <- which(Ds==vec[2])

    if(i==1) {
      return(NA)
    }else{
      entropy_i_1 <- entropy_vec_w[entropy_vec_w$Ws==Ws[i-1] & entropy_vec_w$Ds==vec[2], 3]

      return(as.numeric(abs(vec[3]-entropy_i_1)))
    }

  })

  entropy_vec_w <- entropy_vec_w[!is.na(entropy_vec_w$AD_entropy),]

  if(length(Ds)==1){
    g1 <- ggplot(data = entropy_vec_w,aes(x=Ws,y=AD_entropy))+
      geom_line()+
      scale_colour_discrete(name = "D")+
      theme_bw()+theme(legend.position='bottom')+
      labs(title='Entropy',y='AD')
  }else{
    g1 <- ggplot(data = entropy_vec_w,aes(x=Ws,y=AD_entropy,col=factor(Ds)))+
      geom_line()+
      scale_colour_discrete(name = "D")+
      theme_bw()+theme(legend.position='bottom')+
      labs(title='Entropy',y='AD')
  }


  print(g1)
}



#' Plot absolute difference in moving average of log-posterior density
#'
#' @description
#' This function plot absolute difference in moving average of log-posterior density to select a suitable value for \eqn{D}
#' (number of iterations) in consensus velocity.
#'
#' @param Ws  either an integer or a vector of integers. Absolute differences will be computed for each value in \code{Ws}.
#' @param Ds  a vector of candidate values for depths (number of iterations).
#' @param mcmc_lp  a matrix of log-posterior density, where the rows correspond to chains and columns correspond to iterations.
#'
#' @return If \code{Ws} is an integer, a line plot showing absolute difference in moving average of log-posterior density between different depths.
#' If \code{Ws} is a vector, for each value in \code{Ws} a line will be plotted.
#'
#' @export
#'
#' @examples
#' Ds <- c(1, seq(100,1000,by=100))
#' plot_lp(Ws = 100, Ds = Ds, mcmc_lp = mcmc_lp_result)
plot_lp <- function(Ws, Ds, mcmc_lp){

  # compute moving average of log-posterior
  # cumlative sum for each chain
  cumsum_lp1 <- t(apply(mcmc_lp[,1:max(Ds)], 1, cumsum))
  # sum over rows to get cumsum_lp over 1,2,...,max(W) chains
  cumsum_lp2 <- apply(cumsum_lp1, 2, cumsum)
  # value in [i,j] is the moving average of log-posterior over all lp
  # before the j-th iteration from the first i chains
  moving_average_lp <- matrix(NA, nrow=max(Ws), ncol=max(Ds))
  for(i in 1:max(Ws)){
    moving_average_lp[i,] <- cumsum_lp2[i,]/c(1:max(Ds))/i
  }

  lp_mat <- data.frame(tidyr::expand_grid(Ws,Ds,lp=NA))
  lp_mat$lp <- apply(lp_mat,1,function(vec) {

    return(moving_average_lp[vec[1],vec[2]])

  })

  lp_mat$AD_lp <- apply(lp_mat,1,function(vec) {

    i <- which(Ws==vec[1])
    j <- which(Ds==vec[2])

    if(j==1) {
      return(NA)
    }else{
      lp_i_1 <- lp_mat[lp_mat$Ds==Ds[j-1] & lp_mat$Ws==vec[1], 3]

      return(as.numeric(abs(vec[3]-lp_i_1)))
    }

  })

  lp_mat <- lp_mat[!is.na(lp_mat$AD_lp),]

  if(length(Ws)==1){
    g1 <- ggplot(data = lp_mat,aes(x=Ds,y=AD_lp))+
      geom_line()+
      scale_colour_discrete(name = "W")+
      theme_bw()+theme(legend.position='bottom')+
      labs(title='Moving average of log-posterior',y='AD')
  }else{
    g1 <- ggplot(data = lp_mat,aes(x=Ds,y=AD_lp,col=factor(Ws)))+
      geom_line()+
      scale_colour_discrete(name = "W")+
      theme_bw()+theme(legend.position='bottom')+
      labs(title='Moving average of log-posterior',y='AD')
  }


  print(g1)
}








