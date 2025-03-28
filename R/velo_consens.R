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
#'
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param state_inits a matrix of initial values for the state. Each column corresponds to one set of initial values, with rows corresponding to cells.
#' Can use the output from \code{generate_k}.
#' @param empirical output from \code{get_empirical}.
#' @param mcmc a list of mcmc setup: :
#' \itemize{
#'  \item \code{prep_niter} number of iterations in the preparation step where the state is fixed.
#'  \item \code{prep_burn_in} number of iterations to discard as burn-in in the preparation step.
#'  \item \code{prep_thinning} after burn-in in the preparation step, save the sample at every \code{thinning} iterations.
#'  \item \code{comp_niter} number of iterations in the full algorithm.
#'  \item \code{comp_burn_in} number of iterations to discard as burn-in in the full algorithm.
#'  \item \code{comp_thinning} after burn-in in the full algorithm, save the sample at every \code{thinning} iterations.
#' }
#' @param epsilon a small positive value used to compute empirical \eqn{\tau}. Better to use the same value as used in \code{get_empirical}.
#' @param n_chains number of chains to run. Default to run the same number as the provided \code{state_inits}.
#' @param n_cores number of cores used for parallel computing.
#' @param cluster_type see \code{makeCluster} for parallel computing.
#' @param verbose if TRUE, print the running time for the preparation step and complete algorithm for the first chain.
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
#' mcmc_list <- list(prep_niter = 1000, prep_burn_in = 500, prep_thinning = 1,
#'     comp_niter = 1000, comp_burn_in = 500, comp_thinning = 1)
#' set.seed(1)
#' consensus_result <- velo_consens(u.obs = u.obs, s.obs = s.obs, state_inits = k.inits,
#'    empirical = empirical, mcmc = mcmc_list, epsilon = 1e-3, n_cores = 3, n_chains = 10)
velo_consens <- function(u.obs, s.obs, state_inits, empirical=list(), mcmc=list(), epsilon=1e-3,
                         n_chains=NULL, n_cores=NULL, cluster_type='FORK', verbose=TRUE){

  # check input
  if(!is.vector(u.obs) || !is.vector(s.obs)) {stop("s and u must be a vector")}
  if(!is.matrix(state_inits)) {stop("state_inits must be a matrix")}
  if(!is.list(empirical)) {stop("empirical must be a list")}
  if(!is.list(mcmc)) {stop("mcmc must be a list")}
  if(!is.null(n_chains) && (n_chains>ncol(state_inits))) {
    stop("the number of initial state is less than the required number of chains")
  }

  # mcmc settings (default)
  if(is.null(mcmc$prep_niter)) {prep_niter = 1000} else {prep_niter = mcmc$prep_niter}
  if(is.null(mcmc$prep_thinning)) {prep_thinning = 5} else {prep_thinning = mcmc$prep_thinning}
  if(is.null(mcmc$prep_burn_in)) {prep_burn_in = 0} else {prep_burn_in = mcmc$prep_burn_in}
  if(is.null(mcmc$comp_niter)) {comp_niter = 1000} else {comp_niter = mcmc$comp_niter}
  if(is.null(mcmc$comp_thinning)) {comp_thinning = 5} else {comp_thinning = mcmc$comp_thinning}
  if(is.null(mcmc$comp_burn_in)) {comp_burn_in = 0} else {comp_burn_in = mcmc$comp_burn_in}

  # parallel chains setup (default)
  if(is.null(n_cores)) {n_cores <- parallel::detectCores()-1}
  if(is.null(n_chains)) {n_chains <- 1}

  # --- for reproducibility ----
  # a seed for each chain
  seed <- sample(1:1e8,n_chains)

  # ---- priors ------
  df <- empirical$params
  # alpha and gamma
  alpha1.hat <- df$alpha1.hat; gamma.hat <- df$gamma.hat
  if(is.null(df$alpha1.sd)) {alpha1.sd <- alpha1.hat/6} else {alpha1.sd <- df$alpha1.sd}
  if(is.null(df$gamma.sd)) {gamma.sd <- gamma.hat/6} else {gamma.sd <- df$gamma.sd}

  # variance parameters
  if(is.null(df$invgamma.f)) {invgamma.f <- 1} else {invgamma.f <- df$invgamma.f}

  # hyper parameters for mu_tau, var_tau
  if(length(unique(empirical$k.hat))==1) {
    mm <- rep(mean(empirical$tau.hat),2); vv <- rep(var(empirical$tau.hat),2)
  } else {
    mm <- tapply(empirical$tau.hat, empirical$k.hat, mean)
    if(!any(table(empirical$k.hat)==1)) {
      vv <- tapply(empirical$tau.hat, empirical$k.hat, var)
    } else {
      vv <- rep(var(empirical$tau.hat),2)
    }
  }
  if(is.null(df$hyper.f)) {hyper.f <- 2} else {
    hyper.f <- df$hyper.f
  }

  # t0
  mu_t0 <- df$t0.hat
  if(is.null(df$t0.sd)) {t0.sd <- df$t0.hat/3} else {t0.sd <- df$t0.sd}
  var_t0 <- t0.sd^2

  # lambda
  if(is.null(df$lambda.lower)) {lambda.lower <- 0} else {lambda.lower <- df$lambda.lower}
  if(is.null(df$lambda.upper)) {lambda.upper <- 5} else {lambda.upper <- df$lambda.upper}

  # variance
  sigma.hat.u <- df$sigma.hat.u; sigma.hat.s <- df$sigma.hat.s

  # u0
  u0min1 <- df$u0min1; u0min2 <- df$u0min2
  u0max1 <- df$u0max1; u0max2 <- df$u0max2

  # util function
  show_end <- function(i,mode=c('preparation step','complete step'),time_elapse){
    system(sprintf('echo "\n%s"', paste(mode,'for chain',i,'completed:',time_elapse)))
  }

  # ---- run parallel chain ------
  cl <- parallel::makeCluster(n_cores,type=cluster_type)
  result <- pbapply::pblapply(1:n_chains, cl=cl, function(i){
  # result <- lapply(1:n_chains, function(i){

    tryCatch({
      # -------------- preparation step ------------

      k.random <- state_inits[,i]

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
      q_initial <- runif(1, 0, 2)
      while(q_initial==1 || q_initial==gamma.hat){
        q_initial <- runif(1, 0, 2)
      }

      start1 <- Sys.time()
      chain <- hvelo0_mcmc(u_obs=matrix(u.obs,ncol=1), s_obs=matrix(s.obs,ncol=1), niter = prep_niter,
                           burn_in = prep_burn_in, thinning = prep_thinning,
                           mu_alpha = alpha1.hat, sig_alpha = alpha1.sd, mu_gamma = gamma.hat, sig_gamma = gamma.sd,
                           mu_star_base = mm, sig_star_base = hyper.f*mm, eta_base = vv, nu_base = hyper.f*vv,
                           mu_0 = mu_t0, var_0 = var_t0, q_lower = lambda.lower, q_upper = lambda.upper,
                           sigma_u_2_hat = sigma.hat.u^2, sigma_s_2_hat = sigma.hat.s^2, invgamma_f = invgamma.f,
                           u01_lower = u0min1, u01_upper = u0min2,
                           target_accept = 0.44, verbose=FALSE,
                           k_fix = k.random, alpha1_1_initial = alpha1.hat, gamma_initial = gamma.hat, tau_initial = tau.random,
                           mu_tau_initial = mm, var_tau_initial = vv,
                           t0_initial = max(c(max(tau.random[k.random==1]),mu_t0))+0.1, q_initial = q_initial,
                           sigma_u_2_initial = sigma.hat.u^2, sigma_s_2_initial = sigma.hat.s^2, u01_initial = mean(c(u0min1, u0min2)))
      end1 <- Sys.time()

      if(verbose){
        if(i == 1){
          diff_time <- difftime(end1,start1)
          show_end(i,mode='preparation step',paste(round(diff_time, digits = 3),units(diff_time)))
        }
      }

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

      start2 <- Sys.time()

      set.seed(seed[i]-1)
      chain2 <- hvelo0_mcmc(u_obs=matrix(u.obs,ncol=1), s_obs=matrix(s.obs,ncol=1), niter = comp_niter,
                            burn_in = comp_burn_in, thinning = comp_thinning,
                            mu_alpha = alpha1.hat, sig_alpha = alpha1.sd, mu_gamma = gamma.hat, sig_gamma = gamma.sd,
                            mu_star_base = mm, sig_star_base = hyper.f*mm, eta_base = vv, nu_base = hyper.f*vv,
                            mu_0 = mu_t0, var_0 = var_t0, q_lower = lambda.lower, q_upper = lambda.upper,
                            sigma_u_2_hat = sigma.hat.u^2, sigma_s_2_hat = sigma.hat.s^2, invgamma_f = invgamma.f,
                            u01_lower = u0min1, u01_upper = u0min2,
                            target_accept = 0.44, verbose=FALSE,
                            k_initial = k.random, alpha1_1_initial = alpha1_1_initial, gamma_initial = gamma_initial, tau_initial = tau_initial,
                            p_initial = mean(k.random==1), mu_tau_initial = mu_tau_initial, var_tau_initial = var_tau_initial,
                            t0_initial = t0_initial, q_initial = q_initial, sigma_u_2_initial = sigma_u_2_initial, sigma_s_2_initial = sigma_s_2_initial,
                            u01_initial = u01_initial)

      end2 <- Sys.time()
      if(verbose){
        if(i == 1){
          diff_time <- difftime(end2,start2)
          show_end(i,mode='complete step',paste(round(diff_time, digits = 3),units(diff_time)))
        }
      }

      return(chain2)
    }, error = function(e){
      system(sprintf('echo "\n%s"', paste('Caught an error! Check the output for info')))
      return(paste('Error:',e))
    },
    warning = function(w){
      system(sprintf('echo "\n%s"', paste('Caught a warning! Check the output for info')))
      return(paste('Warning:',w))
    })

  }) #end lapply

  parallel::stopCluster(cl)

  return(result)

}


#' Plot absolute difference in entropy
#'
#' @description
#' This function plot absolute difference in entropy to select a suitable value for \eqn{W} (number of chains) in consensus velocity.
#'
#' @import ggplot2
#'
#' @param Ws a vector of candidate values for widths in an increasing order (number of chains).
#' @param Ds either an integer or a vector of integers (in an increasing order). Absolute differences will be computed for each value in \code{Ds}.
#' @param Ds_label labels of depths in the plot. In case there is thinning or burnin, \code{Ds} should be the index in the saved samples,
#' and \code{Ds_label} should be the original index in the whole sampling process. Should be of the same length as Ds.
#' @param mcmc_k a list, each corresponding to the (saved) samples of k from a single chain.
#'
#' @return If \code{Ds} is an integer, a line plot showing absolute difference in entropy between different widths.
#' If \code{Ds} is a vector, for each value in \code{Ds} a line will be plotted.
#'
#' @export
#'
#' @examples
#' Ws <- c(1, 2, 3)
#' plot_entropy(Ws = Ws, Ds = 500, mcmc_k = mcmc_k_result)
plot_entropy <- function(Ws, Ds, Ds_label=NULL, mcmc_k){

  if(!is.null(Ds_label) & (length(Ds)!=length(Ds_label))) {stop("Ds and Ds_label should be of the same length")}
  if(is.null(Ds_label)) {Ds_label=Ds}
  # Entropy
  # ---- a list of length = no of different D, each is a matrix of total_number_of_chain*n_obs ----
  # depend on values in D
  k_mat_list_by_d <- lapply(Ds, function(d) {
    t(sapply(mcmc_k, function(k_output) {
      k_result <- k_output[[d]]
    }))
  })

  # ---- a list of list, the first level is of length = n_W, the second level is of length = n_D ----
  # within the list of list is a vector recording p(induction) for every cell

  # depend on the values in D and W
  prob_list <- NULL
  for (i in 1:length(Ws)) {
    prob_list[[i]] <- lapply(1:length(Ds), function(j) {
      if(Ws[i]==1) {
        vec <-  ifelse(k_mat_list_by_d[[j]][1,]==1, 1, 0)
      }else{
        vec <-  colMeans(k_mat_list_by_d[[j]][1:Ws[i],]==1)
      }
      return(vec)
    })
  }

  entropy_vec_w <- data.frame(tidyr::expand_grid(Ws,Ds=Ds_label,entropy=NA))
  entropy_vec_w$entropy <- apply(entropy_vec_w, 1, function(vec) {
    i <- which(Ws==vec[1])
    j <- which(Ds_label==vec[2])

    p_c <- prob_list[[i]][[j]]

    entropy <- ifelse((p_c==0) | (p_c==1), 0, -(p_c*log(p_c)+(1-p_c)*log(1-p_c)))

    return(sum(entropy))
  })

  entropy_vec_w$AD_entropy <-  apply(entropy_vec_w, 1, function(vec) {
    i <- which(Ws==vec[1])
    j <- which(Ds_label==vec[2])

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








