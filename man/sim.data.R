library(truncnorm)
# function for u(t) and s(t)
# time-dependent alpha in both states
u_func3 <- function(u0,t,t0,k,alpha1,alpha2,lambda,beta){
  tau_k <- t-t0[k]
  u <- u0[k]*exp(-beta*tau_k)+alpha1[k]/beta*(1-exp(-beta*tau_k))+alpha2[k]/(beta-lambda[k])*(exp(-beta*tau_k)-exp(-lambda[k]*tau_k))
  return(u)
}

s_func3 <- function(u0,s0,t,t0,k,alpha1,alpha2,lambda,beta,gamma){
  tau_k <- t-t0[k]
  s <- alpha1[k]/gamma+1/(gamma-beta)*(beta*u0[k]-alpha1[k]+beta*alpha2[k]/(beta-lambda[k]))*exp(-beta*tau_k)+
    (s0[k]-beta/(gamma-beta)*(u0[k]-alpha1[k]/gamma+alpha2[k]/(gamma-lambda[k])))*exp(-gamma*tau_k)-
    beta*alpha2[k]/(beta-lambda[k])/(gamma-lambda[k])*exp(-lambda[k]*tau_k)
  return(s)
}

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

u0[2] <- u_func3(u0,t=t0[2],t0=t0,k=1,alpha1,alpha2,lambda,beta)
u0

s0[2] <- s_func3(u0,s0,t=t0[2],t0=t0,k=1,alpha1,alpha2,lambda,beta,gamma)

s0

# alphas for repression
alpha1 <- c(alpha1,u0[1]*beta)
alpha2 <- c(alpha2,alpha1[2]-alpha1[1]+alpha2[1]*exp(-lambda[1]*(t0[2]-t0[1])))

### Mean counst for each gene
s.mean <- s_func3(u0,s0,t=t,t0=t0,k=state,alpha1,alpha2,lambda,beta=beta,gamma=gamma)

u.mean <- u_func3(u0,t=t,t0=t0,k=state,alpha1,alpha2,lambda,beta=beta)

#plot the phase protrait
par(mfrow=c(1,3))
plot(s.mean,u.mean,col='red',type='l')

#u versus t
plot(t,u.mean,type='l',col='red')
abline(v=t0[2],col='grey')
#s versus t
plot(t,s.mean,type='l',col='red')
abline(v=t0[2],col='grey')

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

