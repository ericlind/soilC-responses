/* Following combination of Stan manual (v2.14.0) p 147-148
* and the Stan users group thread
* "Dealing with divergent transitions in example model: From Wishart to LKJ priors for correlation."
* 
* Adds additional offset matrix z to help hyperprior correlation matrix
*/

data{
    // indexes
    int<lower=1> N;		// number plots
    int<lower=1> J;		// number sites
    int<lower=1> L;		// number site-level predictors
    int<lower=1> K;		// number of plot-level predictors
    int<lower=1, upper=J> site[N]; // parent site of plot

    // variables
    vector[N] y;		// response
    matrix[N, K] x;		// plot-level predictor matrix including intercept
    matrix[J, L] u;	// site-level predictors including site intercept
}
parameters {
    matrix[K, J] z;		//     cholesky_factor_corr[K] L_Omega;    vector<lower=0,upper=pi()/2>[K] tau_unif; 	// reparameterization
    matrix[L, K] gamma;		// group coeffs    real<lower=0> sigma;		// prediction error scale
}
transformed parameters {
    matrix[J, K] beta;		 // indiv coeffs by group
    vector<lower=0>[K] tau; 	// prior scale   
    for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]); // reparameterization
    beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * z)'; // slopes are deterministic plus random offset
}
model {    to_vector(z) ~ normal(0, 1);
    L_Omega ~ lkj_corr_cholesky(2);    to_vector(gamma) ~ normal(0, 5);// Likelihood    
      y ~ normal(rows_dot_product(x, beta[site]), sigma);}
generated quantities {
    vector[N] log_lik;
    vector[K] betasite;
    for (n in 1:N) {
        betasite =  to_vector(beta[site[n]]); // vectorize for the lpdf to work
	log_lik[n] = normal_lpdf(y[n] | (x[n] * betasite), sigma); // likelihood of plot observation under model
    }
}



