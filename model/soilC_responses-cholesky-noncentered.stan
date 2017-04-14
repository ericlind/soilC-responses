/* Following combination of Stan manual p 147-148
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
    matrix[K, J] z;		// 
    matrix[L, K] gamma;		// group coeffs
}
transformed parameters {
    matrix[J, K] beta;		 // indiv coeffs by group
    vector<lower=0>[K] tau; 	// prior scale   
    for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]); // reparameterization
    beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * z)'; // slopes are deterministic plus random offset
}
model {
    L_Omega ~ lkj_corr_cholesky(2);
      y ~ normal(rows_dot_product(x, beta[site]), sigma);



