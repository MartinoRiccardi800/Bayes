//########################################################
//###  NORMAL REGRESSION 	       ########################
//#######################################################

///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 0> N;       // number of data
	//int<lower = 0> p_fix;   // number of covariates, fixed effects
	int<lower = 0> ngr;	// number of groups
	//int<lower = 0> M[N];    // numerosities in each genus: array
	
	real<lower = 0> Y[N];  	// response vector
	vector[N] X;
	vector[ngr] D;
	//matrix[N, p_fix] X;     // design matrix (fixed effects)
	matrix[117, 3] G;     	// gropus (dummy) allocation  2925 75
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[ngr] beta;        	// regression coefficients 
	//vector[ngr] theta;      	// (group specific) random effects
	vector[ngr] sigma2_beta;	// variances for the prior on beta
	//vector[ngr] sigma2_theta;	// variances for the prior on theta
	real<lower = 0> tau2;
	real<lower = 0> theta;
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters 
{
	vector[N] mu;
	vector[ngr] mu2;
	for(w in 0:((ngr/3)-1)){
	for(i in 1:117){
	  mu[i+ 117*w] = (row(G, i) * (segment(beta, 3*w + 1 , 3)) ) * X[i+ 117*w] ;
	}
	}
	for(p in 1:ngr){
	  mu2[p] = theta * D[p];
	}
}

////////////////// MODEL ////////////////////////
model {
  
	// Likelihood     
	for (s in 1:N)
	{
		// print(" s = ", s);
		Y[s] ~ normal(mu[s], pow(tau2, 0.5));  
	} 

	for (j in 1:ngr) 
	{
	 	beta[j] ~ normal(mu2[j], pow(sigma2_beta[j], 0.5));
		sigma2_beta[j] ~ inv_gamma(2., 20000.);
	}
  
  theta ~ normal(1,1);
  
	tau2 ~ inv_gamma(2., 20000.);
	

}

////////////////// GENERATED QUANTITIES ////////////////////////
generated quantities 
{
  	vector[N] log_lik;
  	for (j in 1:N){
    		log_lik[j] = normal_lpdf(Y[j] | mu[j], pow(tau2, 0.5));
  	}
}
