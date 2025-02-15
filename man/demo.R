library(ggplot2)
library(pheatmap)
# ------- Fit a single gene ----------

# get empirical estimates for parameters from the data
empirical <- get_empirical(u.obs = u.obs, s.obs = s.obs, u.quant = c(0.35,0.8), s.quant = c(0.35,0.8), alpha='max')

# generate multiple initial values for state k
# first type (based on the distance to the centre points of two states)
set.seed(231)
k.inits1 <- generate_k(u.obs = u.obs, s.obs = s.obs, empirical = empirical, n_inits = 50, type = 'centre')

# second type (based on the distance to the closest points of two states)
set.seed(7439)
k.inits2 <- generate_k(u.obs = u.obs, s.obs = s.obs, empirical = empirical, n_inits = 50, type = 'min')

# combine and shuffle
set.seed(234)
rand_ind <- sample(1:100,size=100)
k.inits <- cbind(k.inits1,k.inits2)[,rand_ind]

# consensus velocity
# number of chains
Width <- 3
Depth <- 1000

# define other prior parameters, if not provided, default values will be used
# standard deviations of alpha1, gamma and switching point
empirical$params$alpha1.sd=5
empirical$params$gamma.sd=0.4
empirical$params$t0.sd=3
# upper and lower bound of lambda
empirical$params$lambda.lower=0; empirical$params$lambda.upper=5
# the parameter in the inverse-gamma prior for the two variance parameters: IG(invgamma.f, invgamma.f * empirical)
empirical$params$invgamma.f=1
# the parameter in the prior for mu_tau (also for var_tau), control the variance: TruncN(mean, (hyper.f*mean)^2)
empirical$params$hyper.f=2


# mcmc setup (if not provided, default values will be used)
mcmc_list <- list(prep_niter = 1000, prep_burn_in = 500, prep_thinning = 1,
                  comp_niter = Depth, comp_burn_in = 0, comp_thinning = 1)

set.seed(1)
consensus_result <- velo_consens(u.obs = u.obs, s.obs = s.obs, state_inits = k.inits, empirical = empirical,
                                 mcmc = mcmc_list, epsilon = 1e-3, n_cores = 3, n_chains = Width)


set.seed(234)
rand_ind <- sample(1:100,size=100)
k_lp_result <- lapply((1:100)[rand_ind], function(w) {

  load(paste0('~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/code/velocity_MCMC/sim_4_13_1_consensus/run_result/complete/complete_',w,'.RData'),)

  return(list(k=my_list$k_output,lp=my_list$log_post))
})

# get posterior samples for k and log-posterior
k_lp_result <- lapply(consensus_result, function(l) {

  # load(paste0('~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/code/velocity_MCMC/sim_4_13_1_consensus/run_result/complete/complete_',w,'.RData'),)

  return(list(k=l$k_output,lp=l$log_post))
})
# a list, each corresponding to the samples of k from a single chain
mcmc_k_result <- lapply(k_lp_result,function(l) l$k)
# a matrix: n_chains * n_iterations
mcmc_lp_result <- t(sapply(k_lp_result,function(l) l$lp))
rm(k_lp_result)


# candidate values for W
Ws <- c(1, 2, 3)
Ws <- c(1, 5, seq(10,100,by=5))
Ws <- c(1,seq(2,10,by=2))
# plot absolute difference in entropy, conditional on a fixed D (or several Ds).
plot_entropy(Ws = Ws, Ds = 500, mcmc_k = mcmc_k_result)

# candidate values for D
Ds <- c(1, seq(100,1000,by=100))

# plot absolute difference in moving average of log-posterior, conditional on a fixed W (or several Ws).
plot_lp(Ws = 1, Ds = Ds, mcmc_lp = mcmc_lp_result)


# suppose we decide to use all the chains and the last 100 (saved) samples in each chain
# below we clean up the results
load("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/code/velocity_MCMC/sim_4_13_1_consensus/combined_W100_D9k_10k_thin10.RData")

# the final output of the function below is a list of length = no. of chains. Within each list,
# there are matrices for univariate parameters ('params'),
# for 'tau', 'mus', 'v', 't', 'k', 'hyper'
# the rows for of every matrix correspond to samples, with nrow = length(ind)
combined <- result_combine(consensus_result = consensus_result, Width = 10, ind = 900:1000)

# to show fitted phase portrait (with some thinning)
plot_fit(combined_result = combined, u.obs = u.obs, s.obs = s.obs, thinning = 10, title = 'Example', cex = 1)

# to compute posterior probability of induction for every cell
combined_result_k <- do.call(rbind, lapply(combined, function(l) l$k))
p_ind <- apply(combined_result_k,2,function(x) mean(x==1))

# plot the probability in phase portrait
ggplot()+
  geom_point(aes(x=s.obs,y=u.obs,col=p_ind))+
  labs(x='s',y='u')+
  scale_color_gradient2(midpoint=0.5, low="blue", mid="yellow",
                        high="red", space ="Lab", name='p (induction)')+
  theme_bw()+theme(legend.position='bottom')


# to show predicted velocity for a subset of the cells
plot_predicted_velocity(combined_result = combined, u.obs = u.obs, s.obs = s.obs, delta = 1,
                        obs_ind = c(1,10,15,20,25,30,35,40,80,100,105,115,120),
                        thinning = 10, cex = 0.2, transparency = 0.5)


# to check the joint posterior distribution between parameters and their correlation
# below we plot 2d density plots for selected parameters
# the parameter names can be one of the following
# 'alpha1_1','gamma','lambda','t0','u01','sigma_u_2','sigma_s_2','p'
# 't[1]', 't[2]', .....
# 'tau[1]','tau[2], .....
plot_distribution_2d(combined_result = combined, param1 = 't0', param2 = 'gamma')
plot_distribution_2d(combined_result = combined, param1 = 't0', param2 = 't[100]')


# ------ Posterior predictive checks -------------
# --- a single replicate ---
# compare generated observations
set.seed(4)
ppc_single(u.obs = u.obs, s.obs = s.obs, combined_result = combined)

# --- multiple replicates -----
set.seed(4)
ppc_multiple(u.obs = u.obs, s.obs = s.obs, combined_result = combined, n_replicate = 200, prob = 0.95,
             quantiles=seq(0.05,0.95,by=0.05),
             u.quant = c(0.35,0.8), s.quant = c(0.35,0.8), gamma.hat.obs = empirical$params$gamma.hat)

# the function plots histograms for ppp values and also outputs them
set.seed(134)
ppp_values <- ppp_mixed(u.obs = u.obs, s.obs = s.obs, combined_result = combined, n_replicate = 100, prob = 0.9)

hist(c(ppp_values$ppp.u.ind,ppp_values$ppp.s.ind),breaks = 20)
hist(c(ppp_values$ppp.u.repres,ppp_values$ppp.s.repres),breaks = 20)

# ------- Post-processing step -----------------
load('man/y.RData')
load('man/k.mode.RData')
load('man/p.hat.RData')
load('man/var.logt.RData')
load('man/u.obs.RData'); load('man/s.obs.RData')
# the post-processing step is implemented using rjags
# the output is a list of length = n_chains, within each list is a matrix
# the row corresponds to each iteration, columns correspond to variables (mu_c, sigma_2_c, rho_c, b=-log(beta))
post_result <- mcmc_shared_time(x = x, k.mode = k.mode, var.logt = var.logt, p.hat = p.hat, n_chains = 3,
                                niter = 200, burn_in = 100, thinning = 1)

par(mfrow=c(4,2))
plot(post_result, ask=T)

# to obtain posterior samples for shared time
shared_time_samples <- get_shared_time(post_result = post_result, n = 400, G = 26)
plot(shared_time_samples, ask=T)

# compute posterior probability of each cell having a smaller latent time than another cell and order them
# by posterior median
p_order <- order_probability(shared_time_samples = shared_time_samples)

# visualize the probabilities
my_palette <- colorRampPalette(c("blue", 'yellow','red'))(n = 400)
pheatmap(p_order,scale = "none",fontsize = 14,show_rownames = F,show_colnames = F,
         color = my_palette,border_color=NA,cluster_rows = FALSE, cluster_cols = FALSE,na_col = 'white')

# visualize posterior median of shared time (normalized on [0,1]) on PCA space
shared_time_median <- summary(shared_time_samples)$quantiles[,3]
shared_time_pca(u.obs = u.obs, s.obs = s.obs, shared_time = shared_time_median)

# visualize posterior median of shared time for a specific gene in the phase portrait
g=1
ggplot(data = data.frame(u=u.obs[,g],s=s.obs[,g],t=shared_time_median))+
  geom_point(aes(x=s,y=u,col=t))+
  theme_bw()+
  labs(title=paste0("Gene-shared latent time"))+
  scale_color_viridis_c()+
  theme(legend.position = 'right')



# compute velocity vectors, after adjusting parameters by beta
combined_all_genes <- lapply(1:26,function(g) {

  file_name <- list.files(path=paste0('/Users/zhanghuizi/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/code/velocity_MCMC/t_shared_simulation/',g,'gene'),pattern="combined*", full.names="TRUE")
  load(file_name)

  for(i in 1:60){
    colnames(combined[[i]]$params)[3] <- 'lambda'
  }

  return(combined)

})
par(mfrow=c(1,1))

velocity_pca(combined_all_genes = combined_all_genes, post_result = post_result,
             groups = rep(c('blue','yellow'),each=200),
             u.obs = u.obs, s.obs = s.obs,
             obs_ind = seq(40,400,by=40), thinning = 1, cex = 0.2)

# ------ save data --------
sim.data <- list(u.obs = u.obs, s.obs = s.obs, u.mean = u.mean, s.mean = s.mean,
                 t = t, tau = tau, state = state, t0 = t0[2], alpha1_1 = alpha1[1],
                 beta = 1, gamma = gamma, lambda = lambda[1], u01 = u0[1],
                 sigma_s = sigma_s, sigma_u = sigma_u)
save(sim.data, file='data/sim.data.rda')

post.summary <- list(u.obs = u.obs, s.obs = s.obs, x = y, k.mode = k.mode, var.logt = var.logt, p.hat = p.hat)
save(post.summary, file='data/post.summary.rda')
