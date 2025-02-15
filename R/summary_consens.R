#' Combine MCMC samples from multiple chains
#'
#' @description
#' This function combines MCMC samples for parameters from multiple chains into a clearer structure, based on selected \eqn{W} and \eqn{D}.
#'
#' @param consensus_result output from \code{velo_consensus}.
#' @param Width selected value for \eqn{W}. Only the first \code{Width} chains will be saved in the output.
#' @param ind sample index at which MCMC samples are saved in the output.
#' @param hyper if TRUE, samples for hyper-parameters \eqn{\mu_{\tau,j}} and \eqn{\sigma^2_{\tau,j}} are saved. These will be useful for posterior predictive checks.
#'
#' @return The output returns a list of results for first \code{Width} chains at selected sample index.
#' Within each list, it contains the following components for each chain:
#' \item{params}{a matrix of samples for parameters, including rates, \eqn{\lambda}, switching point,
#' initial condition, variance parameters and hyper-parameter \eqn{p}. Each column corresponds to samples for a variable.}
#' \item{tau}{a matrix of samples for local time \eqn{\tau}. Each column corresponds to the local time for a single cell.}
#' \item{mus}{a matrix of samples for \eqn{u(t),s(t)} from the ODE. Each column corresponds to a single cell.}
#' \item{v, t, k}{a matrix of samples for velocity \eqn{v} from the differential equation for \eqn{s(t)}, for latent time \eqn{t} and for states, respectively.}
#' \item{hyper}{if \code{hyper}=TRUE, a matrix of samples for hyper-parameters \eqn{\mu_{\tau,j}} and \eqn{\sigma^2_{\tau,j}}, otherwise it is NULL.}
#'
#' @export
#'
#' @examples
#' combined <- result_combine(consensus_result = consensus_result, Width = 3, ind = 10:100)
result_combine <- function(consensus_result, Width, ind=NULL, hyper=TRUE){

  params <- c('alpha1_1','gamma','lambda','t0','u01','sigma_u_2','sigma_s_2','p','tau','mus','mu_tau','var_tau','k')

  # the output is a list of length = no. of chains, each chain contains the output of the following:
  # c('alpha1_1','gamma','lambda','t0','u01','sigma_u_2','sigma_s_2','p','tau','mus','mu_tau','var_tau','k')
  # and velocity and t
  # ind is the sample index to be stored in each chain

  # the final output of this function is a list of length = no. of chains, within each list,
  # there are matrices for univariate parameters ('params'),
  # for 'tau', 'mus', 'v', 't', 'k', 'hyper'
  # nrow of every matrix = length(ind)

  if(is.null(ind)) {ind <- 1:consensus_result[[1]]$output_index}

  all_result <- lapply(1:Width, function(w) {

    my_list <- consensus_result[[w]]

    v_mcmc <- lapply(ind, function(i) {
      v_sample <- 1*my_list$mus_output[[i]]$u_mean-my_list$gamma_output[[i]]*my_list$mus_output[[i]]$s_mean
    })

    sample_output <- my_list[sapply(params,function(str) paste0(str,'_output'))]
    sample_output <- lapply(sample_output, function(l) l[ind])
    sample_output[[length(sample_output)+1]] <- v_mcmc
    names(sample_output)[length(sample_output)] <- 'v_output'

    t_mcmc <- lapply(ind, function(i) {
      t_sample <- my_list$tau_output[[i]]+ifelse(my_list$k_output[[i]]==1,0,my_list$t0_output[[i]])
      return(t_sample)
    })
    sample_output[[length(sample_output)+1]] <- t_mcmc
    names(sample_output)[length(sample_output)] <- 't_output'

    return(sample_output)
  })

  # prepare the matrix for univariate variables, each column is a variable
  params1 <- c('alpha1_1','gamma','lambda','t0','u01','sigma_u_2','sigma_s_2','p')
  params2 <- 'tau'; params3 <- 'mus'; params4 <- 'v'; params5 <- 't'; params6 <- 'k'; params7 <- c('mu_tau','var_tau')

  all_mat_list <- lapply(1:Width, function(w) {

    single_result <- all_result[[w]]

    # save univariate variables to a matrix
    mat1 <- sapply(params1, function(var_name) {
      single_result[[paste0(var_name,'_output')]]
    })
    colnames(mat1) <- params1

    # save tau samples to a matrix, each column is a variable (same order as the data)
    n <- length(single_result$tau_output[[1]])
    mat2 <- matrix(unlist(single_result$tau_output),nrow = length(ind),byrow=TRUE)
    if(ncol(mat2)!=n){
      stop('matrix for tau sample should have number of columns equal to data size!')
    }
    colnames(mat2) <- paste0('tau[',1:n,']')

    # save mean values u,s to a matrix, each column is a variable, first n cols are u, n+1 to 2n cols are s
    mat3_u <- matrix(unlist(lapply(single_result$mus_output, function(l) l$u_mean)),nrow=length(ind),byrow = TRUE)
    mat3_s <- matrix(unlist(lapply(single_result$mus_output, function(l) l$s_mean)),nrow=length(ind),byrow = TRUE)
    if((ncol(mat3_u)!=n) || (ncol(mat3_s)!=n)){
      stop('matrix for u_mean (or s_mean) sample should have number of columns equal to data size!')
    }
    colnames(mat3_u) <- paste0('u[',1:n,']')
    colnames(mat3_s) <- paste0('s[',1:n,']')

    mat3 <- cbind(mat3_u, mat3_s)


    # save v samples to a matrix, each column is a variable (same order as the data)
    n <- length(single_result$tau_output[[1]])
    mat4 <- matrix(unlist(single_result$v_output),nrow = length(ind),byrow=TRUE)
    if(ncol(mat4)!=n){
      stop('matrix for v sample should have number of columns equal to data size!')
    }
    colnames(mat4) <- paste0('v[',1:n,']')

    # save t samples to a matrix, each column is a variable (same order as the data)
    n <- length(single_result$tau_output[[1]])
    mat5 <- matrix(unlist(single_result$t_output),nrow = length(ind),byrow=TRUE)
    if(ncol(mat5)!=n){
      stop('matrix for t sample should have number of columns equal to data size!')
    }
    colnames(mat5) <- paste0('t[',1:n,']')

    # save k samples to a matrix, each column is a variable (same order as the data)
    n <- length(single_result$k_output[[1]])
    mat6 <- matrix(unlist(single_result$k_output),nrow = length(ind),byrow=TRUE)
    if(ncol(mat6)!=n){
      stop('matrix for k sample should have number of columns equal to data size!')
    }
    colnames(mat6) <- paste0('k[',1:n,']')

    # save prior params mu_tau, var_tau samples to a matrix, each column is a variable (order: mu_tau_1, mu_tau_2, var_tau_1, var_tau_2)
    if(hyper){

      mat7_mu <- matrix(unlist(single_result$mu_tau_output), nrow=length(ind), byrow = TRUE)
      mat7_var <- matrix(unlist(single_result$var_tau_output), nrow=length(ind), byrow = TRUE)
      if((ncol(mat7_mu)!=2) || (ncol(mat7_var)!=2)){
        stop('matrix for mu_tau (or var_tau) sample should have number of columns equal to 2!')
      }
      colnames(mat7_mu) <- paste0('mu_tau_',1:2)
      colnames(mat7_var) <- paste0('var_tau_',1:2)

      mat7 <- cbind(mat7_mu, mat7_var)
    }else{mat7 <- NULL}


    return(list(params=mat1, tau=mat2, mus=mat3, v=mat4, t=mat5, k=mat6, hyper=mat7))
  })

  return(all_mat_list)

}



#' Plot fitted phase portrait
#'
#' @description
#' This function plots the posterior samples for fitted phase portrait.
#'
#'
#' @param combined_result output from \code{result_combine}.
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param thinning an integer. Only plot for thinned samples.
#' @param title a character string for plot title.
#' @param cex see \code{par}.
#'
#' @return a phase portrait with samples for fitted lines in grey.
#'
#' @export
#'
#' @examples
#' plot_fit(combined_result = combined, u.obs = u.obs, s.obs = s.obs,
#'   thinning = 10, title = 'Example', cex = 1)
plot_fit <- function(combined_result, u.obs, s.obs, thinning=1, title='', cex=1){
  # combined_result is the output of function 'result_combine', a list of of length = no. of chains,
  # within each list, there are matrices for univariate parameters ('params'), for 'tau', 'mus', 'v', 't', 'k','hyper'

  # output: fitted phase portrait

  n <- ncol(combined_result[[1]]$tau)

  combined_result_t <- do.call(rbind, lapply(combined_result, function(l) l$t))
  combined_result_tau <- do.call(rbind, lapply(combined_result, function(l) l$tau))
  combined_result_means <- do.call(rbind, lapply(combined_result, function(l) l$mus))
  combined_result_u <- combined_result_means[,1:n]; combined_result_s <- combined_result_means[,(n+1):(2*n)]
  combined_result_k <- do.call(rbind, lapply(combined_result, function(l) l$k))

  # apply thinning
  plot_ind <- seq(thinning,nrow(combined_result_t),by=thinning)
  plot(s.obs,u.obs,pch=20,cex=cex,xlim=c(0,max(s.obs)),ylim=c(0,max(u.obs)),main=title,xlab='s',ylab='u')
  for (l in plot_ind) {
    for (kk in 1:2) {
      ix <- which(combined_result_k[l,]==kk)
      o <- order(combined_result_tau[l,][ix])
      lines(combined_result_s[l,][ix][o],combined_result_u[l,][ix][o],col='grey')
    }
  }
}


#' Plot predicted velocity in the phase portrait
#'
#' @description
#' This function computes and plots samples for predicted velocities, including posterior mean.
#'
#'
#' @param combined_result output from \code{result_combine}.
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param delta a time increment to compute cell future position (default to 1).
#' @param obs_ind a vector of cell indices for which predicted velocities will be computed.
#' @param thinning compute predicted velocity for thinned samples. The posterior mean is also based on thinned samples.
#' @param cex see \code{par}.
#' @param transparency level of transparency added to the observed values (between 0 and 1).
#'
#' @return a phase portrait overlaid by posterior samples for predicted velocities for selected cells (grey), and posterior mean is shown in red arrows.
#'
#' @export
#'
#' @examples
#' plot_predicted_velocity(combined_result = combined, u.obs = u.obs, s.obs = s.obs, delta = 1,
#'     obs_ind = c(1,10,15,20,25,30,35,40,80,100,105,115,120),
#'     thinning = 10, cex = 0.2, transparency = 0.5)
plot_predicted_velocity <- function(combined_result, u.obs, s.obs, delta=1, obs_ind=NULL,
                                    thinning=1, cex=1, transparency=1){

  n <- length(u.obs)
  if(is.null(obs_ind)) {obs_ind <- sample(1:n,size=min(c(5,n)))}

  combined_result_params <- do.call(rbind, lapply(combined, function(l) l$params))
  combined_result_means <- do.call(rbind, lapply(combined, function(l) l$mus))
  combined_result_u <- combined_result_means[,1:n]; combined_result_s <- combined_result_means[,(n+1):(2*n)]
  rm(combined_result_means)
  combined_result_t <- do.call(rbind, lapply(combined, function(l) l$t))

  compute_mean_vec <- function(tau, t0, k, alpha1_1, beta, gamma, q, u01){

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

    return(c(u_mean, s_mean))
  }
  compute_mean_vec <- Vectorize(compute_mean_vec)


  ix <- seq(thinning,nrow(combined_result_params),by=thinning)
  ccs <- obs_ind

  # a list of length = no. of selected cells, first item is a vector of length 2, first entry for E(u(t+1)-u(t)), second for E(s(t+1)-s(t))
  # second item is a matrix of dim = no. of iters *2, recording all the samples of u(t+1)-u(t), and s(t+1)-s(t)

  delta_mus <- lapply(ccs, function(cc) {

    # dim = 2 * no. of iters
    t_future <- combined_result_t[ix,cc]+delta
    k_future <- ifelse(t_future<=combined_result_params[ix,'t0'],1,2)
    tau_future <- ifelse(k_future==1,t_future,t_future-combined_result_params[ix,'t0'])

    # 2 * no. of iter, first row=all samples of u for a cell, second row for s
    mus <- compute_mean_vec(tau_future,combined_result_params[ix,'t0'],k_future,
                            combined_result_params[ix,'alpha1_1'],
                            beta=1, combined_result_params[ix,'gamma'],
                            combined_result_params[ix,'lambda'],
                            combined_result_params[ix,'u01'])

    temp <- rbind(mus[1,]-combined_result_u[ix,cc],mus[2,]-combined_result_s[ix,cc])

    return(list(avg=rowMeans(unname(temp)),samples=t(unname(temp))))

  })

  # a matrix of dim = no. of selected cells * 2, first column for E(u(t+1)-u(t)), second for E(s(t+1)-s(t))
  avg_delta <- do.call(rbind,lapply(delta_mus,function(l) l$avg))
  # a matrix of dim = no. of selected cells * no. of iter, recording all samples for u(t+1)-u(t)
  delta_u_by_iter <- t(do.call(cbind,lapply(delta_mus,function(l) l$samples[,1])))
  delta_s_by_iter <- t(do.call(cbind,lapply(delta_mus,function(l) l$samples[,2])))

  plot(s.obs,u.obs,type='n',xlab='s',ylab='u')
  points(s.obs,u.obs,pch=20,cex=cex,col=scales::alpha('black',transparency))

  a <- lapply(1:length(ix), function(ii) {
    lapply(1:length(ccs),function(jj) {
      lines(x=c(s.obs[ccs[jj]],s.obs[ccs[jj]]+delta_s_by_iter[jj,ii]),
            y=c(u.obs[ccs[jj]],u.obs[ccs[jj]]+delta_u_by_iter[jj,ii]),col='grey',lwd=0.5)
    })
  })
  shape::Arrows(x0=s.obs[ccs],y0=u.obs[ccs],
                x1=s.obs[ccs]+avg_delta[,2],y1=u.obs[ccs]+avg_delta[,1],
                arr.type = 'curved',col='red',arr.length = 0.2)
}


#' Plot joint posterior distributions for two parameters
#'
#' @description
#' This function visualizes the joint posterior distributions for a pair of parameters to check correlations.
#'
#'
#' @param combined_result output from \code{result_combine}.
#' @param param1,param2 names of the parameters. Must be one of the following:
#' 1) gene-specific parameters 'alpha1_1','gamma','lambda','t0','u01','sigma_u_2','sigma_s_2','p'.
#' 2) latent time 't[1]', 't[2]',...
#' 3) local time 'tau[1]','tau[2]', ...
#' 4) velocity 'v[1]', 'v[2]', ...
#'
#' @return a 2D kernel density estimation for joint posterior distributions.
#' @export
#'
#' @examples
#' plot_distribution_2d(combined_result = combined, param1 = 't0', param2 = 'gamma')
plot_distribution_2d <- function(combined_result, param1, param2){

  combined_result_params <- do.call(rbind, lapply(combined, function(l) l$params))
  combined_result_t <- do.call(rbind, lapply(combined, function(l) l$t))
  combined_result_tau <- do.call(rbind, lapply(combined, function(l) l$tau))
  combined_result_v <- do.call(rbind, lapply(combined, function(l) l$v))

  all_mat <- cbind(combined_result_params, combined_result_t, combined_result_tau, combined_result_v)

  g1 <- ggplot(data=data.frame(x=all_mat[,param1],y=all_mat[,param2]), aes(x=x, y=y) ) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette=4, direction=1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw()+
    theme(legend.position = 'none',text = element_text(size = 15))+
    labs(x=param1, y=param2)

  print(g1)
}




















