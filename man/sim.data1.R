library(ggplot2)
library(gridExtra)
# variance = 1(u), 1(s), 200 obs,

library(truncnorm)
#function for u(t) and s(t)
# time-variable alpha in both states
# correct explicit solution
u_func3 <- function(u0,t,t0,k,alpha1,alpha2,lambda,beta){
  tau_k <- t-t0[k]

  # alpha1 <- c(alpha1,u0[1]*beta)
  # alpha2 <- c(alpha2,alpha1[2]-alpha1[1]+alpha2[1]*exp(-lambda[1]*(t0[2]-t0[1])))
  u <- u0[k]*exp(-beta*tau_k)+alpha1[k]/beta*(1-exp(-beta*tau_k))+alpha2[k]/(beta-lambda[k])*(exp(-beta*tau_k)-exp(-lambda[k]*tau_k))
  return(u)
}

s_func3 <- function(u0,s0,t,t0,k,alpha1,alpha2,lambda,beta,gamma){
  tau_k <- t-t0[k]

  # alpha1 <- c(alpha1,u0[1]*beta)
  # alpha2 <- c(alpha2,alpha1[2]-alpha1[1]+alpha2[1]*exp(-lambda[1]*(t0[2]-t0[1])))
  # s <- s0[k]*exp(-gamma*tau_k)+alpha1[k]/gamma*(1-exp(-gamma*tau_k))+(alpha1[k]-beta*u0[k])/(gamma-beta)*(exp(-gamma*tau_k)-exp(-beta*tau_k))+
  #   beta*alpha2[k]*(exp(-beta*tau_k)/(beta-lambda[k])/(gamma-beta)-exp(-lambda[k]*tau_k)/(beta-lambda[k])/(gamma-lambda[k])-
  #                     exp(-gamma*tau_k)/(gamma-beta)/(gamma-lambda[k]))
  s <- alpha1[k]/gamma+1/(gamma-beta)*(beta*u0[k]-alpha1[k]+beta*alpha2[k]/(beta-lambda[k]))*exp(-beta*tau_k)+
    (s0[k]-beta/(gamma-beta)*(u0[k]-alpha1[k]/gamma+alpha2[k]/(gamma-lambda[k])))*exp(-gamma*tau_k)-
    beta*alpha2[k]/(beta-lambda[k])/(gamma-lambda[k])*exp(-lambda[k]*tau_k)
  return(s)
}


# t <- c(seq(0,4,length.out = 60),seq(4,8,length.out=40),seq(8,12,length.out = 70),seq(12,20,length.out = 30))
# t <- c(seq(0,8,length.out=100), seq(8.1,20, length.out=100))
t <- seq(0,20,length.out=200)
t0 <- c(0,8)

state <- ifelse(t<=t0[2],1,2)

#tau
tau <- t-ifelse(state==1,t0[1],t0[2])

#compute initial u0 and s0
beta <- 1
gamma <- 1.2

u0 <- c()
u0[1] <- 30
s0 <- c()
s0[1] <- beta*u0[1]/gamma
# alphas for induction
alpha1=100;alpha2=alpha1-beta*u0[1];lambda=rep(0.8,2)

q <- lambda/beta

uu0 <- exp(-beta*(t0[2]-t0[1]))
uu <- exp(-beta*tau)


u0[2] <- u_func3(u0,t=t0[2],t0=t0,k=1,alpha1,alpha2,lambda,beta)
u0

s0[2] <- s_func3(u0,s0,t=t0[2],t0=t0,k=1,alpha1,alpha2,lambda,beta,gamma)

s0

# alphas for repression
alpha1 <- c(alpha1,u0[1]*beta)
alpha2 <- c(alpha2,alpha1[2]-alpha1[1]+alpha2[1]*exp(-lambda[1]*(t0[2]-t0[1])))


alpha_tau <- function(tau,t0,k,u0,alpha1,alpha2,lambda){
  # alpha22 <- -alpha1+alpha2*exp(-lambda*t0[2])
  # alpha1 <- c(alpha1,u0[1]*beta)
  # alpha2 <- c(alpha2,alpha1[2]-alpha1[1]+alpha2[1]*exp(-lambda[1]*(t0[2]-t0[1])))

  alpha1[k]-alpha2[k]*exp(-lambda[k]*tau)
  # ifelse(k==1,alpha1-alpha2*exp(-lambda*tau),alpha1*exp(-lambda*tau))
}

alpha_obs <- alpha_tau(tau,t0,u0,k=state,alpha1,alpha2,lambda)
par(mfrow=c(1,2))
plot(t,alpha_obs,col='blue')
abline(h=100)
abline(h=u0[1])
abline(v=t0[2])




### Mean counst for each gene
s.mean <- s_func3(u0,s0,t=t,t0=t0,k=state,alpha1,alpha2,lambda,beta=beta,gamma=gamma)

u.mean <- u_func3(u0,t=t,t0=t0,k=state,alpha1,alpha2,lambda,beta=beta)

##plot the phase protrait
# par(mfrow=c(1,2))
plot(s.mean,u.mean,col='red',xlim=c(0,100),ylim=c(0,110),type='l')
# lines(s.mean,u.mean)
# lines(s.mean.base,u.mean.base,col='black')
# points(s.mean[210],u.mean[210])
# points(s.mean.base[210],u.mean.base[210])

# lines(s.mean.wrong,u.mean.wrong,col='blue')
# lines(s.mean.wrong1,u.mean.wrong1,col='green')

# abline(a=0,b=1,col='grey')

#slope
# apply(s.mean, 2, max)/apply(u.mean, 2, max)

#u versus t
par(mfrow=c(1,2))
plot(t,u.mean,type='l',col='red')
# lines(t,u.mean.base,col='black')
abline(v=t0[2],col='grey')
# lines(t,u.mean.wrong,col='blue')
# lines(t,u.mean.wrong1,col='green')

plot(t,s.mean,type='l',col='red')
# lines(t,s.mean.base,col='black')
abline(v=t0[2],col='grey')
# lines(t,s.mean.wrong,col='blue')
# lines(t,s.mean.wrong1,col='green')

# plot(t,u.mean,type='l',ylim=c(0,110))
# lines(t,s.mean, col='blue')

par(mfrow=c(1,1))

# truncated normal
sigma_u <- 1
sigma_s <- 1

set.seed(937)
u.obs <- rtruncnorm(length(u.mean),a = 0,mean = u.mean,sd = sigma_u)
s.obs <- rtruncnorm(length(s.mean),a = 0,mean = s.mean,sd = sigma_s)

plot(s.obs,u.obs,main=paste0('Truncated-normal\nsd=(',sigma_u,', ', sigma_s,')'),
     xlim=c(0,max(s.obs)),ylim=c(0,max(u.obs)))
lines(s.mean,u.mean,col='red')

plot(s.obs,u.obs,main=paste0('Truncated-normal\nsd=(',sigma_u,', ', sigma_s,')'),
     xlim=c(0,max(s.obs)),ylim=c(0,max(u.obs)),
     col=state)


# k.hat for initialization of gamma, alpha, prior mean of sigma, prior params for u0,s0
k.hat <- rep(1,length(s.obs))
k.hat[s.obs<=quantile(s.obs,p=0.35)&u.obs<=quantile(u.obs,p=0.35)] <- 3
k.hat[s.obs>=quantile(s.obs,p=0.8)&u.obs>=quantile(u.obs,p=0.8)] <- 4
plot(s.obs,u.obs,col=k.hat)

# sted <- u.obs[k.hat==4]
# mean(sted)
# mad(sted)
# thre <- mean(sted)+1.9*mad(sted)
# topu <- u.obs[u.obs>thre]
# mad(topu)
# plot(density(rtruncnorm(10000,a=0,mean=alpha1.hat,sd=mad(topu))))

# sum(u.obs>thre)


gamma.hat <- c(u.obs[k.hat>2]%*%s.obs[k.hat>2]/(s.obs[k.hat>2]%*%s.obs[k.hat>2]))

k.hat[u.obs-gamma.hat*s.obs>0&k.hat<3] <- 1
k.hat[u.obs-gamma.hat*s.obs<=0&k.hat<3] <- 2
plot(s.obs,u.obs,col=k.hat)
k.hat.4 <- k.hat

sigma.hat.u <- mean(c(mad(u.obs[k.hat==3]),mad(u.obs[k.hat==4])))
sigma.hat.s <- mean(c(mad(s.obs[k.hat==3]),mad(s.obs[k.hat==4])))

alpha1.hat <- max(u.obs, s.obs*gamma.hat)

u0min1 = min(u.obs[k.hat == 3])*0.9
u0min2 = max(u.obs[k.hat == 3])*1.1

u0max1 = min(u.obs[k.hat == 4])*0.9
u0max2 = max(u.obs[k.hat == 4])*1.1

# s0max1 = min(s.obs[k.hat == 4])*0.9
# s0max2 = max(s.obs[k.hat == 4])*1.1

# ------- k.hat for initialization of state  (first chain) ---------
# for (i in 1:200) {
#   if(k.hat[i]>2){
#     k.hat[i] <- sample(1:2,1)
#   }
# }
k.hat[u.obs-gamma.hat*s.obs>0] <- 1
k.hat[u.obs-gamma.hat*s.obs<=0] <- 2
plot(s.obs,u.obs,col=k.hat)

# ------- uu.hat for initialization of uu  (first chain, inverse u(t) using the original explicit solution) ---------
uu.hat <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.hat])/(c(mean(c(u0min1,u0min2)),mean(c(u0max1,u0max2)))[k.hat]-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.hat])
uu.hat[uu.hat<=0] <- 0.001
uu.hat[uu.hat>=1] <- 0.999

# based on true k
uu.hat.ktrue <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[state])/(c(mean(c(u0min1,u0min2)),mean(c(u0max1,u0max2)))[state]-c(alpha1.hat,mean(c(u0min1,u0min2)))[state])
uu.hat.ktrue[uu.hat.ktrue<=0] <- 0.001
uu.hat.ktrue[uu.hat.ktrue>=1] <- 0.999

print(1-mean(k.hat!=state))
v <- u.mean-gamma*s.mean

(u.obs[120]-c(alpha1.hat,0)[2])/(c(mean(c(u0min1,u0min2)),mean(c(u0max1,u0max2)))[2]-c(alpha1.hat,0)[2])

tau.hat <- -log(uu.hat)
tau.hat.true <- -log(uu.hat.ktrue)
plot(tau.hat.true,tau,col=state)
abline(0,1)


theta <- data.frame(alpha1=alpha1[1],gamma=gamma,lambda=lambda[1],delta_t=t0[2]-t0[1],
                    sigma_u=sigma_u, sigma_s=sigma_s, u0=u0[1], s0=s0[1])




# ----------- second chain: initialize with different k and tau ------------
# randomly sample cells in empirical steady state (1 or 2). The probability is proportional to
# their distance to observed centre in induction and repression (k.hat.4==1 or 2)
# so we have different initial k and tau now (only differ in empirical steady state)
# all priors are not changed (especially for tau!!!), and only initials for tau and k change

# ix_off <- which(k.hat.4==3)
# ix_on <- which(k.hat.4==4)
ix_steady <- which(k.hat.4>2)

# on_centre <- c(mean(s.obs[k.hat.4==1]),mean(u.obs[k.hat.4==1]))
# off_centre <- c(mean(s.obs[k.hat.4==2]),mean(u.obs[k.hat.4==2]))

# plot(s.obs,u.obs,col=k.hat.4)
# points(on_centre[1],on_centre[2],col='magenta',pch=16,cex=1.5)
# points(off_centre[1],off_centre[2],col='magenta',pch=16,cex=1.5)

dat <- data.frame(u=u.obs,s=s.obs)
dat_scale <- data.frame(scale(dat))

on_centre <- c(mean(dat_scale$s[k.hat.4==1]),mean(dat_scale$u[k.hat.4==1]))
off_centre <- c(mean(dat_scale$s[k.hat.4==2]),mean(dat_scale$u[k.hat.4==2]))

plot(dat_scale$s,dat_scale$u,col=k.hat.4)
points(on_centre[1],on_centre[2],col='magenta',pch=16,cex=1.5)
points(off_centre[1],off_centre[2],col='magenta',pch=16,cex=1.5)

# prob \propto exp(-d^2)

prob1 <- t(sapply(ix_steady, function(ii) {
  obs <- c(dat_scale$s[ii],dat_scale$u[ii])
  LP <- c(-sum((obs-on_centre)^2),-sum((obs-off_centre)^2))
  nc <- -max(LP)
  P <- exp(LP+nc)/sum(exp(LP+nc))

  return(P)
}))


p1 <- rep(NA,200)
p1[ix_steady] <- prob1[,1]
dat <- data.frame(u=u.obs,s=s.obs,p=p1,steady=ifelse(k.hat.4>2,'y','n'),
                  state=state)


set.seed(4542)
k.steady1 <- sapply(1:nrow(prob1), function(ii){
  s <- sample(1:2, size=1, prob=prob1[ii,])
})

# CCR
mean(k.hat[-ix_steady]==state[-ix_steady])
# 1
mean(k.hat[ix_steady]==state[ix_steady])
# 0.6632653


k.random1 <- k.hat
k.random1[ix_steady] <- k.steady1
mean(k.random1==state)
# 0.92
mean(k.random1[-ix_steady]==state[-ix_steady])
# 1
mean(k.random1[ix_steady]==state[ix_steady])
# 0.8367347


# ------- uu.random1 for initialization of uu  (second chain, inverse u(t) using the original explicit solution) ---------
uu.random1 <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random1])/(c(mean(c(u0min1,u0min2)),mean(c(u0max1,u0max2)))[k.random1]-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random1])
uu.random1[uu.random1<=0] <- 0.001
uu.random1[uu.random1>=1] <- 0.999

print(mean(k.random1==state))
# 0.92

tau.random1 <- -log(uu.random1)

which(tau.hat!=tau.random1)==which(k.hat!=k.random1)




# ----------- third chain: a different way to compute k ------------
# randomly sample cells in empirical steady state (1 or 2). The prob2ability is proportional to
# their distance to observed centre in induction and repression (k.hat.4==1 or 2)
# so we have different initial k and tau now (only differ in empirical steady state)
# all priors are not changed (especially for tau!!!), and only initials for tau and k change

ix_off <- which(k.hat.4==3)
ix_on <- which(k.hat.4==4)
ix_steady <- which(k.hat.4>2)

# on_centre <- c(mean(s.obs[k.hat.4==1]),mean(u.obs[k.hat.4==1]))
# off_centre <- c(mean(s.obs[k.hat.4==2]),mean(u.obs[k.hat.4==2]))

# plot(s.obs,u.obs,col=k.hat.4)
# points(on_centre[1],on_centre[2],col='magenta',pch=16,cex=1.5)
# points(off_centre[1],off_centre[2],col='magenta',pch=16,cex=1.5)

dat <- data.frame(u=u.obs,s=s.obs)
dat_scale <- data.frame(scale(dat))


distance <- as.matrix(dist(dat_scale,diag = FALSE,upper = FALSE))

distance_off1 <- distance[which(k.hat.4==1), which(k.hat.4==3)]
close_off1_ix<- as.numeric(names(which.min(apply(distance_off1,1,mean))))

distance_off2 <- distance[which(k.hat.4==2), which(k.hat.4==3)]
close_off2_ix<- as.numeric(names(which.min(apply(distance_off2,1,mean))))

distance_on1 <- distance[which(k.hat.4==1), which(k.hat.4==4)]
close_on1_ix<- as.numeric(names(which.min(apply(distance_on1,1,mean))))

distance_on2 <- distance[which(k.hat.4==2), which(k.hat.4==4)]
close_on2_ix<- as.numeric(names(which.min(apply(distance_on2,1,mean))))

# pink points: point in k=1 closet (smallest average distance) to black
# p_{i1} \prop to exp(-d{i,pink1}^2)
# p_{i2} \prop to exp(-d{i,pink2}^2)
#
# length(distance_on1[which.min(apply(distance_on1,1,mean)),])
# sum(k.hat.4==3)

# v_on1 <- var(distance_on1[which.min(apply(distance_on1,1,mean)),])
# v_on2 <- var(distance_on2[which.min(apply(distance_on2,1,mean)),])

# variance <- mean(c(v_on2, v_on1))

variance_off <- var(c(distance_off1[which.min(apply(distance_off1,1,mean)),],
                      distance_off2[which.min(apply(distance_off2,1,mean)),]))

variance_on <- var(c(distance_on1[which.min(apply(distance_on1,1,mean)),],
                     distance_on2[which.min(apply(distance_on2,1,mean)),]))


plot(dat_scale$s,dat_scale$u,col='grey')
points(dat_scale$s[k.hat.4==3 & k.hat==1],dat_scale$u[k.hat.4==3 & k.hat==1])
points(dat_scale$s[k.hat.4==3 & k.hat==2],dat_scale$u[k.hat.4==3 & k.hat==2],col='red')
points(dat_scale[close_off1_ix,2:1],col='magenta',pch=16,cex=1)
points(dat_scale[close_off2_ix,2:1],col='magenta',pch=16,cex=1)

points(dat_scale$s[k.hat.4==4 & k.hat==1],dat_scale$u[k.hat.4==4 & k.hat==1])
points(dat_scale$s[k.hat.4==4 & k.hat==2],dat_scale$u[k.hat.4==4 & k.hat==2],col='red')
points(dat_scale[close_on1_ix,2:1],col='magenta',pch=16,cex=1)
points(dat_scale[close_on2_ix,2:1],col='magenta',pch=16,cex=1)

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

hist(prob2[1,],breaks = 20)

set.seed(523)
k.steady2 <- sapply(1:nrow(prob2), function(ii){
  s <- sample(1:2, size=1, prob=prob2[ii,])
})
# CCR
mean(k.hat[-ix_steady]==state[-ix_steady])
# 1
mean(k.hat[ix_steady]==state[ix_steady])
# 0.6632653

table(k.steady2, k.hat[ix_steady])
table(k.steady2, state[ix_steady])
par(mfrow=c(2,2))
plot(s.obs,u.obs,col='grey',main='empirical k (no randomization)')
points(s.obs[ix_steady],u.obs[ix_steady],col=k.hat[ix_steady])
plot(s.obs,u.obs,col='grey',main='truth')
points(s.obs[ix_steady],u.obs[ix_steady],col=state[ix_steady])
# plot(s.obs,u.obs,col='grey',main='initial k1')
# points(s.obs[ix_steady],u.obs[ix_steady],col=k.steady1)
plot(s.obs,u.obs,col='grey',main='initial k2')
points(s.obs[ix_steady],u.obs[ix_steady],col=k.steady2)

k.random2 <- k.hat
k.random2[ix_steady] <- k.steady2
mean(k.random2==state)
# 0.79
mean(k.random2[-ix_steady]==state[-ix_steady])
# 1
mean(k.random2[ix_steady]==state[ix_steady])
# 0.5714286


# ------- uu.random2 for initialization of uu  (second chain, inverse u(t) using the original explicit solution) ---------
uu.random2 <- (u.obs-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random2])/(c(mean(c(u0min1,u0min2)),mean(c(u0max1,u0max2)))[k.random2]-c(alpha1.hat,mean(c(u0min1,u0min2)))[k.random2])
uu.random2[uu.random2<=0] <- 0.001
uu.random2[uu.random2>=1] <- 0.999

print(mean(k.random2==state))
# 0.79

tau.random2 <- -log(uu.random2)
which(tau.hat!=tau.random2)==which(k.hat!=k.random2)

# set.seed(923)
# rep.k.steady <- replicate(10000,{
#   val <- sapply(1:nrow(prob2), function(ii){
#     s <- sample(1:2, size=1, prob=prob2[ii,])
#   })
#   return(val)
# })
# dim(rep.k.steady)
# # 98 10000
# ccr <- apply(rep.k.steady,2,function(x) mean(x==state[ix_steady]))
# boxplot(ccr)
# hist(ccr)
# plot(apply(rep.k.steady,1,function(x) mean(x==1)),prob2[,1])
# abline(0,1)

# library(mclust)
# mc1 <- Mclust(data=dat_scale[k.hat.4==1,2:1],G = 1)
# plot(mc1,what = c("classification"))
# summary(mc1)
# mc1
# mc2 <- Mclust(data=dat_scale[k.hat.4==2,2:1],G = 1)
# plot(mc2,what = c("classification"))
#
# on_centre_mc <- mc1$parameters$mean; off_centre_mc <- mc2$parameters$mean
#
# on_centre
#
# plot(dat_scale$s,dat_scale$u,col=k.hat.4)
# points(on_centre_w[1],on_centre_w[2],col='orange',pch=16,cex=1.5)
# points(off_centre_w[1],off_centre_w[2],col='orange',pch=16,cex=1.5)
#
# points(on_centre[1],on_centre[2],col='magenta',pch=16,cex=1.5)
# points(off_centre[1],off_centre[2],col='magenta',pch=16,cex=1.5)
#
# points(on_centre_mc[1],on_centre_mc[2],col='brown',pch=16,cex=1.5)
# points(off_centre_mc[1],off_centre_mc[2],col='brown',pch=16,cex
