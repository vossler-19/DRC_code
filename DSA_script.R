### load stan library 
setwd("/Users/Harley/Documents/DRC Research/2D_Model")

require(dplyr)

require(rstan);
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

tau.prob=function(beta,gamma,rho,B=100){tau=1; for(i in 1:B){
tau=1-exp(-(beta/gamma)*(tau+rho))}; return(tau) 
}

### stochastic SIR simulation  ###
# fix random numbers for reproducibility
datao <- read.csv("drc_1a.csv")
set.seed(124); 	
data<-read.csv("drc_1a.csv");
k=400                                 	# number of datapoints 
kk<-length(data[,1])	
idx=sample(1:kk,kk)							#model thinning  
data=data[idx,-1]								# down to k datapoints 
#data=data[idx,]
data=data[order(data[,1]),]
T=max(data[,1])
not.censored=ifelse(data[,2]<=T,1,0)
data[,2]=pmin(T,data[,2])
data_SIR<-list(k=kk,t0=0,ti=data[,1],Ti=data[,2],event=not.censored)
### MCMC via STAN  ###
fit <- stan(
   file = "DSA_greg2.stan",  # Stan program
   data = data_SIR,       # named list of data
   chains = 2,             # number of Markov chains
   warmup = 1000,          # number of warmup iterations per chain
   iter = 3000,            # total number of iterations per chain
   cores = 2,              # number of cores
   refresh = 1000,         # show progress every 'refresh' iterations
   control=list(adapt_delta=.9) )

print(summary(fit))

#print(plot(fit, show_density=T,pars=c("N","lp__","rho",),include=FALSE,fill_color="green"))
print(plot(fit, show_density=T,pars=c("beta","gamma"),include=TRUE,fill_color="green"))

#str(extract(fit)) list of samples 
#extract(fit$N) e.g., posterior sample of N  

beta_m = mean(extract(fit)$beta)
gamma_m = mean(extract(fit)$gamma)
rho_m = mean(extract(fit)$rho) 
R0 = mean(extract(fit)$R0)
N= mean(extract(fit)$N)
tau=tau.prob(beta=beta_m,gamma=gamma_m,rho=rho_m);
#cat('beta=',beta_m,'gamma=',gamma_m,'R0=',R0,'rho=',rho_m,'tau=',tau,'N=',N,'\n')
cat('beta=',beta_m,'gamma=',gamma_m,'rho=',rho_m,'tau=',tau)
print(traceplot(fit, pars=c("beta", "gamma", "rho")))
print(source('DSA_plot.R'))


