#' Posterior predictive checks based on one single replicate
#'
#' @description
#' This function generates one replicated dataset and compares that with the observed data, for a single gene.
#'
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param combined_result output from \code{result_combine}.
#'
#' @return two plots comparing generated observations with observed data, for unspliced and spliced counts, respectively.
#' @export
#'
#' @examples
#' set.seed(4)
#' ppc_single(u.obs = u.obs, s.obs = s.obs, combined_result = combined)
ppc_single <- function(u.obs, s.obs, combined_result){

  n <- length(u.obs)

  combined_result_means <- do.call(rbind, lapply(combined_result, function(l) l$mus))
  combined_result_u <- combined_result_means[,1:n]; combined_result_s <- combined_result_means[,(n+1):(2*n)]
  rm(combined_result_means)
  combined_result_params <- do.call(rbind, lapply(combined_result, function(l) l$params))

  # total number of MCMC samples
  n_sample <- nrow(combined_result_u)

  # random select a MCMC sample
  ii <- sample(1:n_sample,1)

  u.rep <- truncnorm::rtruncnorm(n,a=0,b=Inf,mean=combined_result_u[ii,],sd = sqrt(combined_result_params[ii,'sigma_u_2']))
  s.rep <- truncnorm::rtruncnorm(n,a=0,b=Inf,mean=combined_result_s[ii,],sd = sqrt(combined_result_params[ii,'sigma_s_2']))

  par(mfrow=c(1,2))
  plot(u.obs,u.rep,ylab='replicate',xlab='observed',main='u',pch=20,cex=0.6)
  abline(0,1,col='red')
  plot(s.obs,s.rep,ylab='replicate',xlab='observed',main='s',pch=20,cex=0.6)
  abline(0,1,col='red')

}



#' Posterior predictive checks based on multiple replicates
#'
#' @description
#' The function generates multiple replicated datasets, and show observed data that is not covered by
#' the credible intervals of the replicates. It then compares empirical \eqn{\gamma} and quantiles between
#' replicates and observed data
#'
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param combined_result output from \code{result_combine}.
#' @param n_replicate number of replicates.
#' @param prob the target probability of the credible interval.
#' @param quantiles a vector of quantiles to compare between the replicates and true data.
#' @param u.quant,s.quant each is a vector of two values for lower and upper quantiles to find
#' empirical repression and induction steady states for unspliced and spliced counts. Should be the same as used in \code{get_empirical}.
#' @param gamma.hat.obs empirical \eqn{\gamma} from the observed data.
#'
#' @return If there are cells whose observed values are not covered by the credible intervals, the first plot will show three panels:
#' 1) a phase portrait where these cells are colored in red, 2) and 3) the observed values and credible intervals will be shown for these cells,
#' for \eqn{u} and \eqn{s}, separately. If all cells are covered, the first plot will not be shown.
#'
#' The second plot with three panels will always be shown: 1) comparing empirical \eqn{\gamma} which is represented as the slope of the lines,
#' 2) and 3) comparing quantiles from the replicated datasets against those from the observed data, for \eqn{u} and \eqn{s}, separately. In each panel,
#' the replicates are shown as grey and the observed data is shown as red.
#' @export
#'
#' @examples
#' set.seed(4)
#' ppc_multiple(u.obs = u.obs, s.obs = s.obs, combined_result = combined,
#'   n_replicate = 200, prob = 0.95, quantiles=seq(0.05,0.95,by=0.05),
#'   u.quant = c(0.35,0.8), s.quant = c(0.35,0.8), gamma.hat.obs = empirical$params$gamma.hat)
ppc_multiple <- function(u.obs, s.obs, combined_result, n_replicate, prob = 0.95,
                         quantiles = seq(0.05,0.95,by=0.05),
                         u.quant = c(0.35,0.8), s.quant = c(0.35,0.8), gamma.hat.obs){

  n <- length(u.obs)

  combined_result_means <- do.call(rbind, lapply(combined_result, function(l) l$mus))
  combined_result_u <- combined_result_means[,1:n]; combined_result_s <- combined_result_means[,(n+1):(2*n)]
  rm(combined_result_means)
  combined_result_params <- do.call(rbind, lapply(combined_result, function(l) l$params))

  # total number of MCMC samples
  n_sample <- nrow(combined_result_u)

  # choose random MCMC samples
  ix <- sample(1:n_sample,n_replicate)

  # n_replicate * n_cell
  u.rep <- t(sapply(ix,function(ii) truncnorm::rtruncnorm(n,a=0,b=Inf,mean=combined_result_u[ii,],sd = sqrt(combined_result_params[ii,'sigma_u_2']))))
  s.rep <- t(sapply(ix,function(ii) truncnorm::rtruncnorm(n,a=0,b=Inf,mean=combined_result_s[ii,],sd = sqrt(combined_result_params[ii,'sigma_s_2']))))

  # HPD from the replicates
  u.rep.quant <- t(coda::HPDinterval(coda::mcmc(u.rep), prob = prob))
  s.rep.quant <- t(coda::HPDinterval(coda::mcmc(s.rep), prob = prob))

  # find cells whose quantiles of replicates do not cover observed values
  # and plot
  ccs <- which(u.obs>u.rep.quant[2,]|u.obs<u.rep.quant[1,]|s.obs>s.rep.quant[2,]|s.obs<s.rep.quant[1,])
  if(length(ccs)>0){

    par(mfrow=c(1,3),cex=1.2)

    plot(s.obs,u.obs,xlim=c(0,max(s.obs)),ylim=c(0,max(u.obs)),pch=20,cex=0.6,xlab='s',ylab='u')
    points(s.obs[ccs],u.obs[ccs],col='red',pch=20,cex=0.6)

    # for u
    ccs <- which(u.obs>u.rep.quant[2,]|u.obs<u.rep.quant[1,])
    plot(1:length(ccs),u.obs[ccs],ylim=c(0,max(u.obs[ccs])*1.2),xlab='Index',ylab='u',pch=20,cex=0.7,main='u')
    for(i in 1:length(ccs)){
      lines(rep(i,2),coda::HPDinterval(coda::mcmc(u.rep[,ccs[i]]),prob=prob),col='grey')
    }
    # for s
    ccs <- which(s.obs>s.rep.quant[2,]|s.obs<s.rep.quant[1,])
    plot(1:length(ccs),s.obs[ccs],ylim=c(0,max(s.obs[ccs])*1.2),xlab='Index',ylab='s',pch=20,cex=0.7,main='s')
    for(i in 1:length(ccs)){
      lines(rep(i,2),coda::HPDinterval(coda::mcmc(s.rep[,ccs[i]]),prob=prob),col='grey')
    }

  }

  # ---- empirical gamma and quantiles from replicates ----
  gamma.rep <- sapply(1:n_replicate, function(ii) {
    k.hat <- rep(1,n)
    k.hat[s.rep[ii,]<=quantile(s.rep[ii,],p=s.quant[1])&
            u.rep[ii,]<=quantile(u.rep[ii,],p=u.quant[1])] <- 3
    k.hat[s.rep[ii,]>=quantile(s.rep[ii,],p=s.quant[2])&
            u.rep[ii,]>=quantile(u.rep[ii,],p=u.quant[2])] <- 4
    c(u.rep[ii,][k.hat>2]%*%s.rep[ii,][k.hat>2]/(s.rep[ii,][k.hat>2]%*%s.rep[ii,][k.hat>2]))

  })

  # quantiles
  quant.obs.u <- quantile(u.obs,p=quantiles)
  quant.rep.u <- lapply(1:n_replicate, function(ii) {
    quantile(u.rep[ii,],p=quantiles)
  })
  quant.obs.s <- quantile(s.obs,p=quantiles)
  quant.rep.s <- lapply(1:n_replicate, function(ii) {
    quantile(s.rep[ii,],p=quantiles)
  })


  par(mfrow=c(1,3))
  plot(0:1,gamma.hat.obs*(0:1),type='l',xlab='',ylab='',main='Empirical gamma')
  for(i in 1:n_replicate){
    lines(0:1,gamma.rep[i]*(0:1),col='grey')
  }
  lines(0:1,gamma.hat.obs*(0:1),col='red',type = 'l')
  legend('topleft',legend=c('replicate','observed'),col=c('grey','red'),lty=1)

  plot(quant.obs.u,quant.rep.u[[1]],type='n',main='Quantiles from u',,xlab='observed',ylab='replicate')
  for(i in 1:n_replicate){
    lines(quant.obs.u,quant.rep.u[[i]],col='grey')
  }
  abline(0,1,col='red')

  plot(quant.obs.s,quant.rep.s[[1]],type='n',main='Quantiles from s',xlab='observed',ylab='replicate')
  for(i in 1:n_replicate){
    lines(quant.obs.s,quant.rep.s[[i]],col='grey')
  }
  abline(0,1,col='red')
}




#' Compute posterior predictive p-values based on mixed predictive distribution
#'
#' @description
#' This function computes posterior predictive p-values based on multiple replicates, by comparing
#' generated counts with observed counts for \eqn{s} and \eqn{u}, separately. The comparison is done
#' for each state separately, and based only on cells with a large posterior probability of belonging to that state.
#'
#'
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param combined_result output from \code{result_combine}.
#' @param n_replicate number of replicates.
#' @param prob a numeric value between 0 and 1. Cells with posterior probabilities of being in a state greater than \code{prob}
#' will be used to compare their observed counts with generated counts.
#'
#' @return two histograms showing p-values for each state separately.
#' Also outputs the p-values for \eqn{u} and \eqn{s} for each state, separately.
#' @export
#'
#' @examples
#' set.seed(134)
#' ppp_values <- ppp_mixed(u.obs = u.obs, s.obs = s.obs,
#'   combined_result = combined, n_replicate = 1000, prob = 0.9)
ppp_mixed <- function(u.obs, s.obs, combined_result, n_replicate, prob=0.9){

  combined_result_params <- do.call(rbind, lapply(combined_result, function(l) l$params))
  combined_result_k <- do.call(rbind, lapply(combined_result, function(l) l$k))
  combined_result_hyper <- do.call(rbind, lapply(combined_result, function(l) l$hyper))

  p_ind <- apply(combined_result_k,2,function(x) mean(x==1))
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

  # induction
  ccs <- which(p_ind>prob)
  kk=1

  # a list of length = n_ccs
  us.replicate.ind <- lapply(ccs,function(cc) {

    i1 <- which(combined_result_k[,cc]==kk)
    i2 <- sample(i1,n_replicate)

    us.sim <- t(sapply(i2,function(ii) {

      mu_tau <- combined_result_hyper[ii,paste0('mu_tau_',kk)]
      var_tau <- combined_result_hyper[ii,paste0('var_tau_',kk)]

      if(kk==1){
        tau.sim <- truncdist::rtrunc(1,'gamma',a=0,b=combined_result_params[ii,'t0'],shape=mu_tau^2/var_tau,scale=var_tau/mu_tau)
      }else{
        tau.sim <- rgamma(1,shape=mu_tau^2/var_tau,rate=mu_tau/var_tau)
      }

      temp <- compute_mean(tau.sim, combined_result_params[ii,'t0'],
                           kk, combined_result_params[ii,'alpha1_1'], beta=1,
                           combined_result_params[ii,'gamma'], combined_result_params[ii,'lambda'],
                           combined_result_params[ii,'u01'])


      val <- truncnorm::rtruncnorm(2,a=0,b = Inf,mean = unlist(temp),sd = sqrt(combined_result_params[ii,c('sigma_u_2','sigma_s_2')]))

      return(val)
    }))

  })

  # each row corresponds to one sample, columns correspond to cells
  u.replicate.ind <- do.call(cbind,lapply(us.replicate.ind,function(l) l[,1]))
  s.replicate.ind <- do.call(cbind,lapply(us.replicate.ind,function(l) l[,2]))
  # ppp values
  ppp.u.ind <- sapply(1:length(ccs),function(i) {
    mean(u.replicate.ind[,i]>u.obs[ccs[i]])
  })
  ppp.s.ind <- sapply(1:length(ccs),function(i) {
    mean(s.replicate.ind[,i]>s.obs[ccs[i]])
  })


  # repression
  ccs <- which(p_ind<(1-prob))
  kk=2

  # a list of length = n_ccs
  us.replicate.repres <- lapply(ccs,function(cc) {

    i1 <- which(combined_result_k[,cc]==kk)
    i2 <- sample(i1,n_replicate)

    us.sim <- t(sapply(i2,function(ii) {

      mu_tau <- combined_result_hyper[ii,paste0('mu_tau_',kk)]
      var_tau <- combined_result_hyper[ii,paste0('var_tau_',kk)]

      if(kk==1){
        tau.sim <- truncdist::rtrunc(1,'gamma',a=0,b=combined_result_params[ii,'t0'],shape=mu_tau^2/var_tau,scale=var_tau/mu_tau)
      }else{
        tau.sim <- rgamma(1,shape=mu_tau^2/var_tau,rate=mu_tau/var_tau)
      }

      temp <- compute_mean(tau.sim, combined_result_params[ii,'t0'],
                           kk, combined_result_params[ii,'alpha1_1'], beta=1,
                           combined_result_params[ii,'gamma'], combined_result_params[ii,'lambda'],
                           combined_result_params[ii,'u01'])


      val <- truncnorm::rtruncnorm(2,a=0,b = Inf,mean = unlist(temp),sd = sqrt(combined_result_params[ii,c('sigma_u_2','sigma_s_2')]))

      return(val)
    }))

  })

  # each row corresponds to one sample, columns correspond to cells
  u.replicate.repres <- do.call(cbind,lapply(us.replicate.repres,function(l) l[,1]))
  s.replicate.repres <- do.call(cbind,lapply(us.replicate.repres,function(l) l[,2]))

  # ppp values
  ppp.u.repres <- sapply(1:length(ccs),function(i) {
    mean(u.replicate.repres[,i]>u.obs[ccs[i]])
  })
  ppp.s.repres <- sapply(1:length(ccs),function(i) {
    mean(s.replicate.repres[,i]>s.obs[ccs[i]])
  })


  par(mfrow=c(1,2))
  hist(c(ppp.u.ind,ppp.s.ind),breaks = 10,xlab='ppp',main='Induction')
  hist(c(ppp.u.repres,ppp.s.repres),breaks = 10,xlab='ppp',main='Repression')

  return(list(ppp.u.ind=ppp.u.ind,ppp.s.ind=ppp.s.ind,ppp.u.repres=ppp.u.repres,ppp.s.repres=ppp.s.repres))
}







