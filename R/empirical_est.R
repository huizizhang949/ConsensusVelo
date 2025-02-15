#' Derive empirical estimates for parameters
#'
#' @description
#' The funciton calculates empirical estimates for all parameters (excluding \eqn{\lambda}) based on the model with a constant transcription rate.
#'
#' @importFrom stats quantile dist mad
#'
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param u.quant a vector of two values for lower and upper quantiles to find empirical repression and induction steady states for unspliced counts.
#' @param s.quant a vector of two values for lower and upper quantiles to find empirical repression and induction steady states for spliced counts.
#' @param alpha a character string indicating which method is used to calculate empirical alpha. One of 'max' or 'avg'.
#' If 'max', empirical value is based on the maximum \eqn{u} or \eqn{s}, otherwise it is based on the average.
#' @param epsilon a small positive value used to compute empirical \eqn{\tau}.
#' @param plot if TRUE, plot empirical \eqn{\gamma} and empirical states.
#'
#'
#' @return The output contains the following items:
#' \item{tau.hat}{a vector of empirical \eqn{\tau} for all cells.}
#' \item{k.hat}{a vector of empirical states \eqn{k} for all cells, consisting of two states (1 for induction, 2 for repression).}
#' \item{k.hat.4}{a vector of empirical states \eqn{k} for all cells, consisting of four states (3 for repression steady state, 4 for induction steady state).}
#' \item{params}{a matrix of empirical estimates for rates, switching point, initial condition (including its prior hyper-parameters) and variance parameters.}
#' \item{alpha.type}{a character string indicating which method is used to calculate empirical alpha. One of 'max' or 'avg'.}
#'
#' @export
#'
#' @examples
#' empirical <- get_empirical(u.obs = u.obs, s.obs = s.obs,
#'   u.quant = c(0.35,0.8), s.quant = c(0.35,0.8), alpha='max')
get_empirical <- function(u.obs, s.obs, u.quant, s.quant, alpha='max', epsilon=1e-3, plot=TRUE){

  C <- length(u.obs)

  k.hat <- rep(1,C)
  k.hat[s.obs<=quantile(s.obs,p=s.quant[1])&
          u.obs<=quantile(u.obs,p=u.quant[1])] <- 3
  k.hat[s.obs>=quantile(s.obs,p=s.quant[2])&
          u.obs>=quantile(u.obs,p=u.quant[2])] <- 4

  gamma.hat <- c(u.obs[k.hat>2]%*%s.obs[k.hat>2]/(s.obs[k.hat>2]%*%s.obs[k.hat>2]))

  k.hat[u.obs-gamma.hat*s.obs>0&k.hat<3] <- 1
  k.hat[u.obs-gamma.hat*s.obs<=0&k.hat<3] <- 2
  if(plot){
    plot(s.obs,u.obs,col=k.hat)
    abline(a=0,b=gamma.hat,col='red')
  }
  k.hat.4 <- k.hat

  sigma.hat.u <- mean(c(mad(u.obs[k.hat==3]),mad(u.obs[k.hat==4])))
  sigma.hat.s <- mean(c(mad(s.obs[k.hat==3]),mad(s.obs[k.hat==4])))

  if(alpha=='avg'){
    alpha1.hat <- max(mean(u.obs[k.hat.4==4]), mean(s.obs[k.hat.4==4]*gamma.hat))
  }else{
    alpha1.hat <- max(u.obs,s.obs*gamma.hat)
  }

  u0min1 = min(u.obs[k.hat == 3])*0.9
  u0min2 = max(u.obs[k.hat == 3])*1.1

  u0max1 = min(u.obs[k.hat == 4])*0.9
  u0max2 = max(u.obs[k.hat == 4])*1.1

  k.hat[u.obs-gamma.hat*s.obs>0] <- 1
  k.hat[u.obs-gamma.hat*s.obs<=0] <- 2

  if(alpha=='avg'){
    uu.hat <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.hat])/(c(mean(c(u0min1,u0min2)),alpha1.hat)[k.hat]-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.hat])
  }else{
    uu.hat <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.hat])/(c(mean(c(u0min1,u0min2)),mean(c(u0max1,u0max2)))[k.hat]-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.hat])
  }

  uu.hat[uu.hat<=0] <- epsilon
  uu.hat[uu.hat>=1] <- 1-epsilon
  tau.hat <- -log(uu.hat)
  t0.hat <- max(tau.hat[k.hat==1])

  params <- data.frame(alpha1.hat=alpha1.hat, gamma.hat=gamma.hat, t0.hat=t0.hat,
                       u01.hat=mean(c(u0min1,u0min2)), u0min1=u0min1, u0min2=u0min2, u0max1=u0max1, u0max2=u0max2,
                       sigma.hat.u=sigma.hat.u, sigma.hat.s=sigma.hat.s)

  empirical <- list(tau.hat=tau.hat, k.hat=k.hat, k.hat.4=k.hat.4, params=params, alpha.type=alpha)

  return(empirical)
}


#' Generate initial values for states for uncertain cells
#'
#' @description
#' The function generates initial values for states for uncertains that will be used when running multiple chains, based on two different ways.
#'
#' @param u.obs a vector of observed unspliced counts.
#' @param s.obs a vector of observed spliced counts.
#' @param empirical output from \code{get_empirical}.
#' @param n_inits number of initializations (chains) to generate for.
#' @param type a character string indicating which method is used to generate states. One of 'centre' or 'min'.
#' If 'centre', the probability to generate states is based on the distance to the centre points, otherwise it is based on the points closest to the steady states.
#' @param plot if TRUE, plot the centre points or closest points.
#'
#' @return a matrix of initial values. Each column corresponds to one set of initial values, with rows corresponding to cells.
#'
#' @export
#'
#' @examples
#' k.inits1 <- generate_k(u.obs = u.obs, s.obs = s.obs, empirical = empirical,
#'   n_inits = 50, type = 'centre')
generate_k <- function(u.obs, s.obs, empirical, n_inits=50, type='min', plot=TRUE){

  k.hat.4 <- empirical$k.hat.4
  k.hat <- empirical$k.hat
  ix_steady <- which(k.hat.4>2)
  dat <- data.frame(u=u.obs,s=s.obs)
  dat_scale <- data.frame(scale(dat))

  if(type=='centre'){

    on_centre <- c(mean(dat_scale$s[k.hat.4==1]),mean(dat_scale$u[k.hat.4==1]))
    off_centre <- c(mean(dat_scale$s[k.hat.4==2]),mean(dat_scale$u[k.hat.4==2]))

    if(plot){
      plot(dat_scale$s,dat_scale$u,col=k.hat.4,main='scaled data', xlab='s', ylab='u')
      points(on_centre[1],on_centre[2],col='yellow2',pch=16,cex=1.5)
      points(off_centre[1],off_centre[2],col='yellow2',pch=16,cex=1.5)
    }


    prob1 <- t(sapply(ix_steady, function(ii) {
      obs <- c(dat_scale$s[ii],dat_scale$u[ii])
      LP <- c(-sum((obs-on_centre)^2),-sum((obs-off_centre)^2))
      nc <- -max(LP)
      P <- exp(LP+nc)/sum(exp(LP+nc))

      return(P)
    }))

    k.inits <- replicate(n_inits, {
      k.steady <- sapply(1:nrow(prob1), function(ii){
        s <- sample(1:2, size=1, prob=prob1[ii,])
      })
      k.random <- k.hat
      k.random[ix_steady] <- k.steady

      return(k.random)
    })

  }else{
    ix_off <- which(k.hat.4==3)
    ix_on <- which(k.hat.4==4)

    distance <- as.matrix(stats::dist(dat_scale,diag = FALSE,upper = FALSE))

    distance_off1 <- distance[(k.hat.4==1) & dat_scale$u>max(dat_scale$u[k.hat.4==3]), which(k.hat.4==3)]
    close_off1_ix<- as.numeric(names(which.min(apply(distance_off1,1,mean))))

    distance_off2 <- distance[(k.hat.4==2) & dat_scale$s>max(dat_scale$s[k.hat.4==3]), which(k.hat.4==3)]
    close_off2_ix<- as.numeric(names(which.min(apply(distance_off2,1,mean))))

    distance_on1 <- distance[(k.hat.4==1) & dat_scale$s<min(dat_scale$s[k.hat.4==4]), which(k.hat.4==4)]
    close_on1_ix<- as.numeric(names(which.min(apply(distance_on1,1,mean))))

    distance_on2 <- distance[(k.hat.4==2) & dat_scale$u<min(dat_scale$u[k.hat.4==4]), which(k.hat.4==4)]
    close_on2_ix<- as.numeric(names(which.min(apply(distance_on2,1,mean))))

    variance_off <- var(c(distance_off1[which.min(apply(distance_off1,1,mean)),],
                          distance_off2[which.min(apply(distance_off2,1,mean)),]))

    variance_on <- var(c(distance_on1[which.min(apply(distance_on1,1,mean)),],
                         distance_on2[which.min(apply(distance_on2,1,mean)),]))


    if(plot){
      plot(dat_scale$s,dat_scale$u,col='grey', main='scaled data', xlab='s', ylab='u')
      points(dat_scale$s[k.hat.4==3 & k.hat==1],dat_scale$u[k.hat.4==3 & k.hat==1])
      points(dat_scale$s[k.hat.4==3 & k.hat==2],dat_scale$u[k.hat.4==3 & k.hat==2],col='red')
      points(dat_scale[close_off1_ix,2:1],col='yellow2',pch=16,cex=1.5)
      points(dat_scale[close_off2_ix,2:1],col='yellow2',pch=16,cex=1.5)

      points(dat_scale$s[k.hat.4==4 & k.hat==1],dat_scale$u[k.hat.4==4 & k.hat==1])
      points(dat_scale$s[k.hat.4==4 & k.hat==2],dat_scale$u[k.hat.4==4 & k.hat==2],col='red')
      points(dat_scale[close_on1_ix,2:1],col='yellow2',pch=16,cex=1.5)
      points(dat_scale[close_on2_ix,2:1],col='yellow2',pch=16,cex=1.5)

    }


    prob2 <- t(sapply(ix_steady, function(ii) {
      obs <- c(dat_scale$s[ii],dat_scale$u[ii])
      if(k.hat.4[ii]==3){
        LP <- c(-sum((obs-dat_scale[close_off1_ix,2:1])^2)/2/variance_off,-sum((obs-dat_scale[close_off2_ix,2:1])^2)/2/variance_off)
      }else{
        LP <- c(-sum((obs-dat_scale[close_on1_ix,2:1])^2)/2/variance_on,-sum((obs-dat_scale[close_on2_ix,2:1])^2)/2/variance_on)
      }

      nc <- -max(LP)
      P <- exp(LP+nc)/sum(exp(LP+nc))

      return(P)
    }))

    k.inits <- replicate(n_inits, {
      k.steady <- sapply(1:nrow(prob2), function(ii){
        s <- sample(1:2, size=1, prob=prob2[ii,])
      })
      k.random <- k.hat
      k.random[ix_steady] <- k.steady

      return(k.random)
    })

  }

  return(k.inits)
}




