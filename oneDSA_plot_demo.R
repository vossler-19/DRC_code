#euler.prob
delta =0 

require(roperators)
require(matrixStats)

nn=1  #(165 in thousands) 
rate= 1#.6

a=gamma_m
b=beta_m
c=beta_m*rho_m

step.size=0.01  # numerical step
T=max(data) + step.size*2 # time horizon (in days )

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

# survival ode with delta 
ode.sur <- function(x, a=gamma, b=beta, c=rho*beta)  -(b*x*(1-x) +a*x*log(x)+c*x)

# without delta 
ode.sur0 <- function(x, a=gamma, b=beta, c=rho*beta)  -(b*x*(1-x) +a*x*log(x)+c*x)


# Plot a histogram of the empirical density of infection
hist(data,prob=T,ylim=c(0,.01))


################ New Code ####################################################
# Randomly sample 500 observations of each parameter from the posterior distirbutions
set.seed(126)
S = sample(3000, 500, replace = FALSE)
N=500

betasamp <- as.data.frame(extract(fit)$beta[S])
gammasamp <- as.data.frame(extract(fit)$gamma[S])
rhosamp <- as.data.frame(extract(fit)$rho[S])


##### Onset 
# Select the median of the posteriors for beta, gamma, rho for density plots
prob = c(0.025, 0.975)
beta = median(extract(fit)$beta)
gamma=median(extract(fit)$gamma)
rho=median(extract(fit)$rho)

res1=euler.one(fun=ode.sur0) # surv fcn

onset <- -ode.sur0(res1)/max(1-res1)

## For each observation of our sample, derive the probability of being susceptible
sdot.t <- data.frame(values=seq(1:nrow(res1)))
for (i in 1:N) {
  beta = betasamp[i,]
  gamma = gammasamp[i,]
  rho = rhosamp[i,]
  
  res1_o=euler.one(fun=ode.sur0) # surv fcn
  
  s <- -ode.sur0(res1_o)/max(1-res1_o)
  sdot.t[,paste0("s",i)] <- s
}
sdot.t <- sdot.t[,-1]
stdevs <- rowSds(as.matrix(sdot.t))
percentile_o <- apply(as.matrix(sdot.t), MARGIN = 1, FUN = quantile, prob=prob)

## Plot the predicted density of infection against the empirical density (hist)
onset_plot <- ggplot() + 
  geom_histogram(aes(x=onset, y=..density..), data = datao, bins = 16, fill="grey93", 
                 color="black") +
  theme_bw() +
  geom_line(aes(x=gr[-1], y=onset), col="red", size=1.5) + 
  geom_line(aes(x=gr[-1], y=percentile_o[1,]),linetype="dotted", col="blue", size=0.8) +
  geom_line(aes(x=gr[-1], y=percentile_o[2,]), col="blue", linetype="dotted", size=0.8) +
  geom_ribbon(aes(ymin=percentile_o[1,], ymax=percentile_o[2,], x=gr[-1]), 
              fill="lightblue", alpha=0.6) + 
  ylab("Density") + xlab("Day of Epidemic") +
  theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.5))) +
  theme(axis.text.x = element_text(size=rel(1.4))) +
  theme(axis.text.y = element_text(size=rel(1.4))) +
  ggtitle("Predicted Density Infection Curve")
print(onset_plot)

## Calculate distribution of n and kinfty
n <- c()
kinfty <- c()
for (i in 1:N) {
  beta = betasamp[i,]
  gamma = gammasamp[i,]
  rho = rhosamp[i,]
  
  res1=euler.one(fun=ode.sur0) # surv fcn
  
  kT=nrow(datao)
  n[i] <- kT/(1-min(res1))
  
  kinfty[i] <- (tau*kT)/(1-min(res1))
}

print(list(mean(n), sd(n), names=c("Mean n", "St. Dev. n")))
print(list(mean(kinfty), sd(kinfty), names=c("Mean kinfty", "St. Dev. kinfty")))

hist(n, main="Effective Population Size")

hist(kinfty, main="Estimated Final Cases")

hist(extract(fit)$R0, main="Distribution of R0", xlab = "R0")

#Save the estimated parameters in a table
model_params <- matrix(c(median(extract(fit)$beta), median(extract(fit)$gamma), 
                         median(extract(fit)$rho), median(extract(fit)$R0), median(n), 
                         median(kinfty),beta_m, gamma_m,  rho_m, mean(extract(fit)$R0), 
                         mean(n), mean(kinfty), sd(extract(fit)$beta), sd(extract(fit)$gamma), 
                         sd(extract(fit)$rho), sd(extract(fit)$R0), sd(n), sd(kinfty)), ncol = 3)
colnames(model_params) <- c("Median", "Mean", "Standard Deviation")
rownames(model_params) <- c("Beta", "Gamma", "Rho", "R0", "n", "kinfty")
#write.csv(model_params, "model_params3_norec.csv")

# Calculates the negative log-likelihood of the model
r1 <- approx(time(res1),res1[,1],datao[,1])$y
ll=0
for (i in 1:k){
  ll = ll + log(r1[i])+log(beta*(1+rho-r1[i])+gamma*log(r1[i]))
}
ll = ll + -k*log(1-r1[k])
print(ll)

