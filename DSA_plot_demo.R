require(roperators)
require(matrixStats)

## Determine step size and T=final time of data
step.size=0.01  # numerical step
T = max(data[,1]) + 2*step.size
gr=seq(0,T,by=step.size)  # grid for some plots 


# one dim euler 
euler.one <- function (t = T, dt = step.size, fun = f, ic=1)
{ 
p <- length(ic)
 n <- t/dt 
 xmat <- matrix(0,ncol=p,nrow=n)
 x <- ic 
 xmat[1,] <- x
 for (i in 2:n) { 
 	x <- x + fun(x)*dt
	xmat[i,] <- x
 } 
ts(xmat,start=0,deltat=dt) 
} 


# general euler 
euler <- function (t = T, dt = step.size, fun = f, ic=c(1,rho,0)) 
{ 
p <- length(ic)
 n <- t/dt 
 xmat <- matrix(0,ncol=p,nrow=n)
 x <- ic 
 xmat[1,] <- x
 for (i in 2:n) { 
 	x <- x + fun(x)*dt
	xmat[i,] <- x
 } 
ts(xmat,start=0,deltat=dt) 
} 

# survivalone dim ode (epidemic curve)  
ode.ds <- function(x, a=gamma, b=beta, c=rho*beta)  -(b*x*(1-x) +a*x*log(x)+c*x)

# full survival model including recovery (sir ode)
ode.sir<-function(x,k1=beta,k2=gamma)
{
	c(-k1*x[1]*x[2], k1*x[1]*x[2]-k2*x[2], k2*x[2])
}


### Sample 500 times from posteriors of beta, gamma, rho
set.seed(124)
S = sample(3000, 500, replace = FALSE)
N=500

betasamp <- as.data.frame(extract(fit)$beta[S])
gammasamp <- as.data.frame(extract(fit)$gamma[S])
rhosamp <- as.data.frame(extract(fit)$rho[S])



### Calculate confidence bounds 
##### Onset 
prob = c(0.025, 0.975)
beta = median(extract(fit)$beta)
gamma=median(extract(fit)$gamma)
rho=median(extract(fit)$rho)

res1=euler.one(fun=ode.ds) # surv fcn
res2=euler(fun=ode.sir,ic=c(1,rho,0))
#hist(data[,1],prob=T,main="Onset")
onset <- -ode.ds(res1)/max(1-res1)
hosp <- gamma*(res2[,2]-rho*exp(-gamma*gr[-1]))/max(1-res1)

# Calculate the density of infection
sdot.t <- data.frame(values=seq(1:nrow(res1)))
for (i in 1:N) {
  beta = betasamp[i,]
  gamma = gammasamp[i,]
  rho = rhosamp[i,]
  
  res1_o=euler.one(fun=ode.ds) # surv fcn
  
  s <- -ode.ds(res1_o)/max(1-res1_o) #point-wise probabilities 
  sdot.t[,paste0("s",i)] <- s
}
# calculate standard deviation and percentile values for this distribution of density values
sdot.t <- sdot.t[,-1]
stdevs <- rowSds(as.matrix(sdot.t))

percentile_o <- apply(as.matrix(sdot.t), MARGIN = 1, FUN = quantile, prob=prob)
  

###############################################################################

### hospitalization
# repeat steps from above but for density of recovery
idot.t <- data.frame(values=seq(1:nrow(res1)))
for (i in 1:N) {
  beta = betasamp[i,]
  gamma = gammasamp[i,]
  rho = rhosamp[i,]
  
  res2_h=euler(fun=ode.sir,ic=c(1,rho,0))
  
  h <- gamma*(res2[,2]-rho*exp(-gamma*gr[-1]))/max(1-res1)
  idot.t[,paste0("s",i)] <- h
}
idot.t <- idot.t[,-1]
stdevs <- rowSds(as.matrix(idot.t))

percentile_h <- apply(as.matrix(idot.t), MARGIN = 1, FUN = quantile, prob=prob)


###################### Density Plots ##############################################

onset_plot <- ggplot() + #geom_density(aes(x=onset), data = datao, adjust=2) + 
  geom_histogram(aes(x=onset, y=..density..), data = datao, bins = 16, fill="grey93", 
                 color="black") +
  theme_bw() +
  geom_line(aes(x=gr[-1], y=onset), col="red", size=1.5) + 
  geom_line(aes(x=gr[-1], y=percentile_o[1,]),linetype="dotted", col="blue", size=0.8) +
  geom_line(aes(x=gr[-1], y=percentile_o[2,]), col="blue", linetype="dotted", size=0.8) +
  geom_ribbon(aes(ymin=percentile_o[1,], ymax=percentile_o[2,], x=gr[-1]), 
              fill="lightblue", alpha=0.6) + 
  ylab("Density") + xlab("Day of Wave 1") +
  theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.5))) +
  theme(axis.text.x = element_text(size=rel(1.4))) +
  theme(axis.text.y = element_text(size=rel(1.4)))
print(onset_plot)

datao_h <- datao %>%
  filter(hosp <= max(onset))

hosp_plot <- ggplot() + #geom_density(aes(x=hosp), data = datao_h, adjust=2) +
  geom_histogram(aes(x=hosp, y=..density..), data = datao_h, bins = 16, fill="grey93", 
                 color="black") +
  theme_bw() +
  geom_line(aes(x=gr[-1], y=hosp), col="red", size=1.5) + 
  geom_line(aes(x=gr[-1], y=percentile_h[1,]),linetype="dotted", col="blue", size=0.8) +
  geom_line(aes(x=gr[-1], y=percentile_h[2,]), linetype="dotted", col="blue", size=0.8) +
  geom_ribbon(aes(ymin=percentile_h[1,], ymax=percentile_h[2,], x=gr[-1]), 
              fill="lightblue", alpha=0.6) +
  ylab("Density") + xlab("Day of Wave 1") +
  theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.5))) +
  theme(axis.text.x = element_text(size=rel(1.4))) +
  theme(axis.text.y = element_text(size=rel(1.4)))
print(hosp_plot)

################################################################################
# Draws out the values of density according to time points in the empirical data
# these values are used for calculating derived parameters in following steps
r1 <- approx(time(res2),res2[,1],data[,1])$y
r2 <- approx(time(res2),res2[,2],data[,1])$y

# Produces value for negative log-likelihood
ll=0
for (i in 1:k){
  ll = ll + (log(r1[i])+log(r2[i])+not.censored[i]*log(gamma)-
    gamma*(data_SIR$Ti[i]-data_SIR$ti[i]))
}
ll = ll + k*(log(beta)-log(1-r1[k]))
ll


# Produce distirbutions of n and kinfty to calculate mean, standard deviation
n <- c()
kinfty <- c()
for (i in 1:N) {
  beta = betasamp[i,]
  gamma = gammasamp[i,]
  rho = rhosamp[i,]
  
  res1=euler.one(fun=ode.ds) # surv fcn
  res2=euler(fun=ode.sir,ic=c(1,rho,0))
  r1 <- approx(time(res2),res2[,1],data[,1])$y
  r2 <- approx(time(res2),res2[,2],data[,1])$y
  
  kT=nrow(datao)
  n[i] <- kT/(1-min(res1))
  
  kinfty[i] <- (tau*kT)/(1-min(res1))
}

print(list(mean(n), sd(n), names=c("Mean n", "St. Dev. n")))
print(list(mean(kinfty), sd(kinfty), names=c("Mean kinfty", "St. Dev. kinfty")))

# Histograms of derived parameters n, R0, kinfty
hist(n, prob=T, main="Effective Population Size")
hist(kinfty, prob=T, main="Estimated Final Cases")
hist(extract(fit)$R0, main="Distribution of R0", xlab = "R0")

# Create a table of model parameters 
model_params1 <- matrix(c(median(extract(fit)$beta), median(extract(fit)$gamma),
                          median(extract(fit)$rho),median(extract(fit)$R0),median(n),
                          median(kinfty),beta_m, gamma_m,  rho_m, mean(extract(fit)$R0), mean(n), mean(kinfty), sd(extract(fit)$beta), 
                         sd(extract(fit)$gamma), sd(extract(fit)$rho), sd(extract(fit)$R0), sd(n), 
                         sd(kinfty)), ncol = 3)
colnames(model_params1) <- c("Median","Mean", "Standard Deviation")
rownames(model_params1) <- c("Beta", "Gamma", "Rho", "R0", "n", "kinfty")
#write.csv(model_params1, "model_params1.csv")



