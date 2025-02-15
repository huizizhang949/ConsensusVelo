#' Implement a post-processing step to infer a gene-shared latent time
#'
#' @description
#' This function uses posterior summaries from fitting individual genes to estimate a shared latent time and splicing rate.
#'
#'
#' @import rjags
#'
#' @param x a cell-by-gene matrix. Each entry is the posterior mean of log gene-wise latent time, conditional on the state with a larger posterior probability \code{(k.mode)}.
#' @param k.mode a cell-by-gene matrix. Each entry is the state (induction: 1, repression: 2) with a larger posterior probability.
#' @param var.logt an array \code{[cell, gene, state]} of posterior variance of log gene-wise latent time, for each state.
#' @param p.hat an array \code{[cell, gene, state]} of posterior probability of belonging to each state for each cell in each gene.
#' @param n_chains number of chains to run.
#' @param niter integer. Total number of MCMC iterations.
#' @param burn_in the length of burn-in period during MCMC.
#' @param thinning the thinning applied after burn-in period.
#'
#' @return The function returns a list of results for all chains. A \code{mcmc.list} object. Within each list, it contains a matrix of samples for each chain.
#' The columns correspond to samples for the following parameters: \eqn{\mu_{c,j}, \  \tilde{\sigma}^2_{c,j}, \   \rho_{c,j}, \  -\log(\beta_g)},
#' for \eqn{j=1,2}, for each cell \eqn{c} and gene \eqn{g}.
#' @export
#'
#' @examples
#' post_result <- mcmc_shared_time(x = y, k.mode = k.mode, var.logt = var.logt,
#'   p.hat = p.hat, n_chains = 2, niter = 200, burn_in = 100, thinning = 1)
mcmc_shared_time <- function(x, k.mode, var.logt, p.hat, n_chains=1, niter=1000, burn_in=0, thinning=1){

  y <- x
  n <- nrow(y); G <- ncol(y)

  # empirical b
  # b=-log(beta)
  b <- c()
  b[1] <- 0
  for(g in 2:G){
    dd <- lapply(1:2,function(j) {
      ix <- which(k.mode[,1]==j & k.mode[,g]==j)
      y[ix,1]-y[ix,g]
    })
    b[g] <- mean(unlist(dd))
  }

  # priors and initials
  # mu_c
  temp <- lapply(1:n, function(cc) {
    if(length(unique(k.mode[cc,]))>1){
      unlist(tapply(y[cc,]+b, k.mode[cc,],mean))
    }else{
      rep(mean(y[cc,]+b), 2)
    }
  })
  mu_c_initial <- do.call(rbind,temp)

  # sigma_2_c
  temp <- lapply(1:n, function(cc) {
    if(length(unique(k.mode[cc,]))>1){
      val <- unlist(tapply(y[cc,]+b, k.mode[cc,],var))
    }else{
      val <- rep(var(y[cc,]+b), 2)
    }
    if(any(is.na(val))){
      val1 <- val[!is.na(val)]
      val <- rep(val1,2)
    }
    return(val)
  })
  sigma_2_c_initial <- do.call(rbind,temp)

  # hyper parameters rho to indicate probabilities of each state
  rho.hat <- t(sapply(1:n, function(cc) {
    val <- mean(p.hat[cc,,1])
    if(val==1){
      return(c(0.99,0.01))
    }else if(val==0){
      return(c(0.01,0.99))
    }
    return(c(val,1-val))
  }))

  model.data <- list(n=n, G=G,  y=y, mu=mu_c_initial, sig_2=mu_c_initial^2, phi=sigma_2_c_initial,
                     rho_hat=rho.hat, sigma_2_hat=var.logt, k=k.mode)

  # the first chain is initialized with all parameters set to their empirical values
  model.inits <- list( list(mu_c=mu_c_initial,sigma_2_c_inv=1/sigma_2_c_initial,
                              rho=rho.hat,
                              b=c(NA,b[2:G]),
                              ".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 66))

  # the rest chains will initialize rho and b randomly
  if(n_chains>1){

    more.inits <- lapply(1:(n_chains-1), function(i){

      set.seed(i)
      rho_temp <- matrix(runif(2*n),nrow=n,ncol=2)
      rho_temp <- rho_temp/rowSums(rho_temp)
      list(mu_c=mu_c_initial,sigma_2_c_inv=1/sigma_2_c_initial,
           rho=rho_temp,
           b=c(NA,rnorm(G-1)),
           ".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = i+1)
    })

    model.inits <- c(model.inits, more.inits)
  }


  model_file <-
    " model {

  for(i in 1:n)
  {
  for(g in 1:G)
  {
    y[i,g] ~ dnorm(mu_c[i,k[i,g]]-b[g], 1/(sigma_2_hat[i,g,k[i,g]]+sigma_2_c[i,k[i,g]]))
  }
  }

  for(i in 1:n){
  for(j in 1:2){
  mu_c[i,j] ~ dnorm(mu[i,j],1/sig_2[i,j])
  sigma_2_c_inv[i,j] ~ dgamma(1,phi[i,j])
  sigma_2_c[i,j] <- 1/sigma_2_c_inv[i,j]
  }
  }

  for(i in 1:n){
    rho[i,1:2] ~ ddirch(rho_hat[i,1:2])
  for(g in 1:G){
  k[i,g] ~ dcat(rho[i,1:2])
  }
  }

  b[1] <- 0
  for(g in 2:G){
  b[g] ~ dnorm(0,0.01)
  }

 }
"


  model.run <- jags.model(file=textConnection(model_file),
                            data=model.data,
                            inits=model.inits,
                            n.chains=n_chains)

  print('Complete compiling')

  update(model.run, n.iter=burn_in)
  print('Complete burnin')

  n.iter = niter-burn_in
  model.coda <- coda.samples(model.run, variable.names=c('mu_c','sigma_2_c','rho','b'),
                               n.iter=n.iter, thin = thinning)

  return(model.coda)

}


#' Compute gene-shared latent for each cell
#'
#' @description
#' This function calculates gene-shared latent time for each cell, based on posterior samples from the post-processing step.
#'
#' @param post_result output from \code{mcmc_shared_time}.
#' @param n number of cells.
#' @param G number of genes.
#'
#' @return a list of results for all chains. A \code{mcmc.list} object. Within each list, it contains a matrix of samples for shared latent time.
#' Each column corresponds to the posterior samples for a single cell.
#' @export
#'
#' @examples
#' shared_time_samples <- get_shared_time(post_result = post_result, n = 400, G = 26)
get_shared_time <- function(post_result, n, G){

  # get index for every param
  b_ind <- 1:G; #-log(beta)
  mu_c_ind <- (b_ind[G]+1):(b_ind[G]+n*2)
  rho_ind <- (mu_c_ind[n*2]+1):(mu_c_ind[n*2]+n*2)
  sigma_2_c_ind <- (rho_ind[n*2]+1):(rho_ind[n*2]+n*2)

  # number of chains
  J <- length(post_result)
  # number of iterations within each chain
  N <- nrow(post_result[[1]])

  mu_c <- coda::mcmc.list(lapply(1:J,function(j) {
    post_result[[j]][,mu_c_ind]
  }))

  rho <- coda::mcmc.list(lapply(1:J,function(j) {
    post_result[[j]][,rho_ind]
  }))


  shared_t <- coda::mcmc.list(lapply(1:J,function(j) {
    temp <- t(sapply(1:N,function(ii) {
      val1 <- exp(mu_c[[j]][ii,1:n])
      val2 <- exp(mu_c[[j]][ii,1:n+n])

      val <- ifelse(rho[[j]][ii,1:n]>0.5,val1,val2)
      return(val)
    }))

    colnames(temp) <- paste0('shared_t[',1:n,']')

    return(coda::mcmc(temp))
  }))

  return(shared_t)

}


#' Calculate posterior probability of one cell before another
#'
#' @description
#' This function computes posterior probability of one cell before another, based on MCMC samples for gene-shared latent time.
#'
#'
#' @param shared_time_samples output from \code{get_shared_time}.
#'
#' @return a matrix of probabilities, of dimension \eqn{n \times n} where \eqn{n} is the number of cells. The cells
#' are ordered by posterior median of shared time from top to bottom, and left to right.
#' @export
#'
#' @examples
#' p_order <- order_probability(shared_time_samples = shared_time_samples)
order_probability <- function(shared_time_samples){

  shared_time_median <- summary(shared_time_samples)$quantiles[,3]

  shared_time_samples <- do.call(rbind,shared_time_samples)

  # number of observations
  n <- ncol(shared_time_samples)
  # total number of iterations across chains
  N <- nrow(shared_time_samples)

  p_order <- matrix(0,nrow=n,ncol=n)
  for (i in 1:N) {
    p_order <- p_order+outer(shared_time_samples[i,],shared_time_samples[i,],'<')
  }
  p_order <- p_order/N

  # order cells by posterior median
  o <- order(shared_time_median,decreasing = TRUE)

  return(p_order[o,o])
}




#' Visualize shared time on PCA space
#'
#' @description
#' This function applies PCA on concatenated matrices for observed spliced and unspliced counts, and then visualize
#' the point estimate of shared latent time after normalizing it on \eqn{[0,1]}.
#'
#'
#' @param u.obs a cell-by-gene matrix of unspliced counts.
#' @param s.obs a cell-by-gene matrix of spliced counts.
#' @param shared_time a vector of (point estimate) of gene-shared latent time.
#'
#' @return a plot with PC2 against PC1, and cells are colored by normalized shared latent time.
#' @export
#'
#' @examples
#' shared_time_pca(u.obs = u.obs, s.obs = s.obs, shared_time = shared_time_median)
shared_time_pca <- function(u.obs, s.obs, shared_time){

  normx <- function(x){
    return(((x-min(x))/(max(x)-min(x))))
  }

  y_all <- cbind(s.obs,u.obs)
  pca <- prcomp(y_all,center=TRUE, scale. = TRUE)

  my_theme <- theme(plot.title = element_text(size = 20, face = "bold"),
                    legend.title=element_text(size=20), axis.text = element_text(size=15),
                    legend.text=element_text(size=15), axis.title=element_text(size=20),
                    legend.position = 'right')

  g1 <- ggplot(data = data.frame(pc1=pca$x[,1],pc2=pca$x[,2],t=normx(shared_time)))+
    geom_point(aes(x=pc1,y=pc2,col=t))+
    theme_bw()+
    scale_color_viridis_c()+
    labs(x='PC1',y='PC2')+
    my_theme

  print(g1)
}



#' Visualize predicted velocity on PCA space
#'
#' @description
#' This function computes the posterior mean for the inverse of the splicing rate, which is then used to adjust all
#' the parameters across genes. After rescaling, the predicted velocities across genes are computed and then
#' visualized on the PCA space.
#'
#' @param combined_all_genes a list of results for all genes. Each item is the output from \code{result_combine} for a single gene.
#' @param post_result output from \code{mcmc_shared_time}.
#' @param u.obs a cell-by-gene matrix of unspliced counts.
#' @param s.obs a cell-by-gene matrix of spliced counts.
#' @param groups a vector of strings for colors to show different groups of cells, if available.
#' @param delta a time increment to compute cell future position (default to 1).
#' @param obs_ind a vector of cell indices for which predicted velocities will be computed. If not provided, 5 cells will be randomly selected.
#' @param thinning compute predicted velocity for thinned samples. The posterior mean is also based on thinned samples.
#' @param n_cores number of cores used for parallel computing.
#' @param cluster_type see \code{makeCluster} for parallel omputing.
#' @param cex see \code{par}.
#' @param transparency level of transparency added to the observed values (between 0 and 1).
#'
#' @return a plot with PC2 against PC1, overlaid by posterior samples for predicted velocities for selected cells (grey),
#' and posterior mean is shown in red arrows. If \code{groups} is provided, cells are colored based on the groups.
#' @export
#'
#' @examples
#' velocity_pca(combined_all_genes = combined_all_genes, post_result = post_result,
#'     u.obs = u.obs, s.obs = s.obs, n_cores=2,]
#'     obs_ind = seq(40,400,by=40), thinning = 1, cex = 0.2)
velocity_pca <- function(combined_all_genes, post_result, u.obs, s.obs, groups = NULL,
                         delta = 1, obs_ind = NULL, thinning = 1,
                         n_cores=NULL, cluster_type='FORK',
                         cex = 1, transparency = 1){

  if(is.null(n_cores)) {n_cores <- parallel::detectCores()-1}

  # combine samples into a single matrix
  post_result <- do.call(rbind, post_result)

  n <- nrow(u.obs); G <- ncol(s.obs); if(is.null(obs_ind)) {obs_ind <- sample(1:n,size=min(c(5,n)))}
  ccs <- obs_ind
  # total number of iterations across chains for single-gene results
  N <- nrow(do.call(rbind, lapply(combined_all_genes[[1]], function(l) l$params)))

  # b=-log(beta)
  b_ind <- 1:G
  b_samples <- post_result[,b_ind]

  # a=1/beta
  a_pmean <- colMeans(exp(b_samples))

  ix <- seq(thinning,N,by=thinning)

  avg_delta_u_adj <- matrix(NA,nrow=length(ccs),ncol=G);avg_delta_s_adj <- matrix(NA,nrow=length(ccs),ncol=G);
  delta_u_by_iter_adj <- array(NA,dim = c(length(ccs),G,length(ix))); delta_s_by_iter_adj <- array(NA,dim = c(length(ccs),G,length(ix)))

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

  # for (g in 1:G) {
  cl <- parallel::makeCluster(n_cores,type=cluster_type)
  delta_mus_all <- pbapply::pblapply(1:G, function(g) {

    combined_result_params <- do.call(rbind, lapply(combined_all_genes[[g]], function(l) l$params))
    combined_result_means <- do.call(rbind, lapply(combined_all_genes[[g]], function(l) l$mus))
    combined_result_u <- combined_result_means[,1:n]; combined_result_s <- combined_result_means[,(n+1):(2*n)]
    rm(combined_result_means)
    combined_result_t <- do.call(rbind, lapply(combined_all_genes[[g]], function(l) l$t))

    # ----  compute E(u(t+1)-u(t)) and E(s(t+1)-s(t)) for every sample, for some selected iterations ix ----
    # a list of length = no. of cells, first item is a vector of length 2, first entry for E(u(t+1)-u(t)), second for E(s(t+1)-s(t))
    # second item is a matrix of dim = no. of iters *2, recording all the samples of u(t+1)-u(t), and s(t+1)-s(t)
    delta_mus <- lapply(ccs, function(cc) {

      # dim = 2 * no. of iters
      t_future <- combined_result_t[ix,cc]*a_pmean[g]+delta
      k_future <- ifelse(t_future<=combined_result_params[ix,'t0']*a_pmean[g],1,2)
      tau_future <- ifelse(k_future==1,t_future,t_future-combined_result_params[ix,'t0']*a_pmean[g])

      # 2 * no. of iters, first row=all samples of u for a cell, second row for s
      mus <- compute_mean_vec(tau_future,combined_result_params[ix,'t0']*a_pmean[g],k_future,
                              combined_result_params[ix,'alpha1_1']/a_pmean[g],
                              beta=1/a_pmean[g], combined_result_params[ix,'gamma']/a_pmean[g],
                              combined_result_params[ix,'lambda'],
                              combined_result_params[ix,'u01'])

      temp <- rbind(mus[1,]-combined_result_u[ix,cc],mus[2,]-combined_result_s[ix,cc])

      return(list(avg=rowMeans(unname(temp)),samples=t(unname(temp))))

    })

  })
  parallel::stopCluster(cl)

  for(g in 1:G){
    delta_mus <- delta_mus_all[[g]]
    # a matrix of dim = no. of cells * 2, first column for E(u(t+1)-u(t)), second for E(s(t+1)-s(t))
    avg_delta <- do.call(rbind,lapply(delta_mus,function(l) l$avg))
    avg_delta_u_adj[,g] <- avg_delta[,1]; avg_delta_s_adj[,g] <- avg_delta[,2]
    # a matrix of dim = no. of cells * no. of iter, recording all samples for u(t+1)-u(t)
    delta_u_by_iter_adj[,g,] <- t(do.call(cbind,lapply(delta_mus,function(l) l$samples[,1])))
    delta_s_by_iter_adj[,g,] <- t(do.call(cbind,lapply(delta_mus,function(l) l$samples[,2])))

  }
  # ------ pca by concatenating u and s-------

  y_all <- cbind(s.obs,u.obs)
  pca <- prcomp(y_all,center=TRUE, scale. = TRUE)

  avg_future_pca_adj <- scale(cbind(s.obs[ccs,]+avg_delta_s_adj,u.obs[ccs,]+avg_delta_u_adj),
                              center = pca$center,scale = pca$scale)%*%pca$rotation

  plot(pca$x[,1],pca$x[,2],type='n',main=paste0(''),xlab='PC1',ylab='PC2')
  if(!is.null(groups)){
    points(pca$x[,1],pca$x[,2],pch=20,cex=cex,col=scales::alpha(groups,transparency))
  }else{
    points(pca$x[,1],pca$x[,2],pch=20,cex=cex,col=scales::alpha('black',transparency))
  }
  a <- lapply(1:length(ix), function(ii) {

    future_pca <- scale(cbind(s.obs[ccs,]+delta_s_by_iter_adj[,,ii],u.obs[ccs,]+delta_u_by_iter_adj[,,ii]),
                        center = pca$center,scale = pca$scale)%*%pca$rotation

    lapply(1:length(ccs),function(jj) {
      lines(x=c(pca$x[,1][ccs[jj]],future_pca[jj,1]),
            y=c(pca$x[,2][ccs[jj]],future_pca[jj,2]),col='grey',lwd=0.5)
    })
  })
  shape::Arrows(x0=pca$x[,1][ccs],y0=pca$x[,2][ccs],
         x1=avg_future_pca_adj[,1],y1=avg_future_pca_adj[,2],
         arr.type = 'curved',col='red',arr.length = 0.2)




}






