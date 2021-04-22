
functions {

real[] oneSIR(real t , real[] y , real[] parms , real[] rdata , int[] idata){

real beta = parms[1];
real gamma = parms[2];
real rho = parms[3];
real dydt[1];

dydt[1] = -beta*y[1]*(1+rho-y[1])-gamma*y[1]*log(y[1]);

return dydt;
 }
}

data {

  int<lower=0> k;   //number of infected individuals 
  real<lower=0> t0;     // initial time 
  real<lower=0> ti[k];  // times of infection  
    
}

transformed data{

real x_r[0];
int x_i[0];

}

parameters {

  real<lower=0.15> beta;               
  real<lower=0, upper=beta> gamma;               
  real<lower=0.001, upper=0.02> rho;   
   
}

transformed parameters{
real R0 = beta/gamma;

}



model {
real temp[k,1];
real parms[3];
real init[1];

parms[1] = beta; 
parms[2]= gamma;
parms[3]= rho;
init[1] = 1;


temp = integrate_ode_rk45(oneSIR,init,t0,ti,parms,x_r,x_i,1.0E-7,1.0E-7,1.0E7); 

for(i in 1:k){
target += log(temp[i,1])+log(beta*(1+rho-temp[i,1])+gamma*log(temp[i,1]));
}
target += -k*log(1-temp[k,1])+gamma_lpdf(beta|.002,.002)+gamma_lpdf(gamma|.002,.002);  

}

