functions{
  // function that describes the changes in CAR+ counts in FO B cells
  real CAR_positive_FOB(real time){
    real F0 = exp(11.722278); real B0 = exp(4.475064); real n = 4.781548 ; real X = 6.943644 ; real q = 5;
    real value = F0 + (B0 * time^n) * (1 - ((time^q)/((X^q) + (time^q))));
    return value;
   }

   real CAR_negative_MZB(real time){
    real M0 = exp(14.06); real nu = 0.0033; real b0 = 20.58;
    real value = M0 * (1 + exp(-nu * (time - b0)^2));
    return value;
   }

   real Total_FoB(real time){
    real M0 = exp(16.7); real nu = 0.004; real b0 = 20;
    real value = M0 * (1 + exp(-nu * (time - b0)^2));
    return value;
   }

   // function that contains the ODE equations to be used in ODE solver
   array[] real ODE_sys(real time,  array[] real y, array[] real parms, array[] real rdata, array[] int idata) {
     // the array of parameters invloved in ode equations
     real alpha = parms[1];
     real beta = parms[2];
     real mu = parms[3];
     real delta = parms[4];
     real lambda = parms[5];
     real nu = parms[6];

     real t0 = 4.0;
     real alpha_tau = alpha/(1 + exp(nu * (time-t0)^2));
     real mu_tau = mu/(1 + exp(nu * (time-t0)^2));
     real beta_tau = beta/(1 + exp(nu * (time-t0)^2));

     // the system of ODEs
     array[3] real dydt;
     // CAR positive GCB cells in WT
     dydt[1] = alpha_tau * Total_FoB(time) - delta * y[1];
     // CAR positive MZB cells in WT
     dydt[2] = mu_tau * Total_FoB(time) + beta_tau * CAR_negative_MZB(time) - lambda * y[2];

     // CAR positive MZB cells in N2KO
     dydt[3] = beta_tau * CAR_negative_MZB(time) - lambda * y[3];
     return dydt;
   }

   array[,] real solve_ODE_sys(array[] real solve_time, array[] real init_cond, array[] real parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     array[numdim, 3] real y_sol;
     y_sol = integrate_ode_rk45(ODE_sys, init_cond, 4.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }
}
data{
  int<lower  = 1> numObs1;                           // number of observations for donor fractions may be different than cell counts
  int<lower  = 1> numPred;
  array[numObs1] real solve_time;                    // time points for donor fractions
  array[numObs1] real<lower=0> CAR_MZ_counts;
  array[numObs1] real<lower=0> CAR_GC_counts;
  array[numObs1] real<lower=0> CAR_MZN2_counts;
  array[numObs1] real<lower=0> CAR_GCN2_counts;
  array[numPred] real ts_pred;
}

parameters{
  // parameters to sample with boundary conditions
  real<lower = 0> alpha;
  real<lower = 0> beta;
  real<lower = 0> mu;
  real<lower = 0> delta;
  real<lower = 0> lambda;
  real<lower = 0> nu;
  real<lower = 0> M0N2;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  }


transformed parameters{
  array[numObs1, 3] real y_hat;     // declaring the array for ODE solution

  array[numObs1] real CAR_GCcounts_mean;
  array[numObs1] real CAR_MZcounts_mean;
  array[numObs1] real CAR_GCN2counts_mean;
  array[numObs1] real CAR_MZN2counts_mean;

  array[6] real parms;                  // declaring the array for parameters
  array[3] real init_cond;              // declaring the array for state variables

  real CAR_GC0 = exp(11.5);              // transformed parameters for better/faster sampling
  real CAR_MZ0 = exp(10.8);              // transformed parameters for better/faster sampling
  real CAR_MZ0N2k0 = exp(M0N2);          // transformed parameters for better/faster sampling

  // initial conditions and parameters
  init_cond[1] = CAR_GC0;
  init_cond[2] = CAR_MZ0;
  init_cond[3] = CAR_MZ0N2k0;

  parms[1] = alpha;
  parms[2] = beta;
  parms[3] = mu;
  parms[4] = delta;
  parms[5] = lambda;
  parms[6] = nu;

  y_hat[1] = init_cond;
  // solution of the system of ODEs for the predictor values
  y_hat[2:] = solve_ODE_sys(solve_time[2:], init_cond, parms);

  for (i in 1:numObs1){
    CAR_GCcounts_mean[i] = y_hat[i, 1];
    CAR_MZcounts_mean[i] = y_hat[i, 2];
  }

  for (i in 1:numObs1){
    CAR_GCN2counts_mean[i] = y_hat[i, 1];
    CAR_MZN2counts_mean[i] = y_hat[i, 3];
  }
}

model{
  // prior distribution for model parameters
  alpha ~ normal(0.01, 0.5);
  beta ~ normal(0.01, 0.5);
  mu ~ normal(0.01, 0.5);
  nu ~ normal(0.01, 0.5);
  delta ~ normal(0.04, 0.3);
  lambda ~ normal(0.1, 0.3);
  M0N2 ~ normal(8, 1.5);

  sigma1 ~ normal(0, 2.5);
  sigma2 ~ normal(0, 2.5);
  sigma3 ~ normal(0, 2.5);

  // model fitting on to data
  log(CAR_GC_counts) ~ normal(log(CAR_GCcounts_mean), sigma1);
  log(CAR_MZ_counts) ~ normal(log(CAR_MZcounts_mean), sigma2);
  log(CAR_GCN2_counts) ~ normal(log(CAR_GCN2counts_mean), sigma1);
  log(CAR_MZN2_counts) ~ normal(log(CAR_MZN2counts_mean), sigma3);
}

generated quantities{
  // ODE predictions
  array[numPred, 3] real y_hat_pred;
  // variables for model predictions
  array[numPred] real y1_mean_pred, y2_mean_pred, y3_mean_pred, y4_mean_pred;
  array[numPred] real FOtoCARMZ_pred, MZtoCARMZ_pred, FOtoCARGC_pred;
  array[numPred] real mu_pred, alpha_pred, beta_pred;
  // variables for model predictions with stdev
  array[numPred] real CAR_MZcounts_pred, CAR_GCcounts_pred;
  array[numPred] real CAR_MZN2counts_pred, CAR_GCN2counts_pred;
  // Residuals
  vector[numObs1] resid_d1, resid_d2; vector[numObs1] resid_d3, resid_d4;
  // log likelihoods
  vector[numObs1] log_lik1, log_lik2; vector[numObs1] log_lik3, log_lik4;

   //ODE solution
   y_hat_pred[1] = init_cond;
   y_hat_pred[2:] = solve_ODE_sys(ts_pred[2:], init_cond, parms);

   // model predictions with stdev
   for (i in 1:numPred){
     //CAR GC
     y1_mean_pred[i] = y_hat_pred[i, 1];
     CAR_GCcounts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma1));
     //CAR MZ
     y2_mean_pred[i] = y_hat_pred[i, 2];
     CAR_MZcounts_pred[i] = exp(normal_rng(log(y2_mean_pred[i]), sigma2));

     //CAR GC in N2KO
     y3_mean_pred[i] = y_hat_pred[i, 1];
     CAR_GCN2counts_pred[i] = exp(normal_rng(log(y3_mean_pred[i]), sigma1));
     //CAR MZ in N2KO
     y4_mean_pred[i] = y_hat_pred[i, 3];
     CAR_MZN2counts_pred[i] = exp(normal_rng(log(y4_mean_pred[i]), sigma3));

     // Influx into CAR MZ
     mu_pred[i] = mu/(1 + exp(nu *(ts_pred[i] - 4.0)^2));
     beta_pred[i] = beta/(1 + exp(nu *(ts_pred[i] - 4.0)^2));
     FOtoCARMZ_pred[i] = mu_pred[i] * Total_FoB(ts_pred[i])/y2_mean_pred[i];
     MZtoCARMZ_pred[i] = beta_pred[i] * CAR_negative_MZB(ts_pred[i])/y2_mean_pred[i];
     // Influx into CAR GC
     alpha_pred[i] = alpha/(1 + exp(nu *(ts_pred[i] - 4.0)^2));
     FOtoCARGC_pred[i] = alpha_pred[i] * Total_FoB(ts_pred[i])/y1_mean_pred[i];
   }

   // calculating the log predictive accuracy for each point
   for (n in 1:numObs1) {
     resid_d1[n] = log(CAR_GC_counts[n]) - log(CAR_GCcounts_mean[n]);
     resid_d2[n] = log(CAR_MZ_counts[n]) - log(CAR_MZcounts_mean[n]);

     log_lik1[n] = normal_lpdf(log(CAR_GC_counts[n]) | log(CAR_GCcounts_mean[n]), sigma1);
     log_lik2[n] = normal_lpdf(log(CAR_MZ_counts[n]) | log(CAR_MZcounts_mean[n]), sigma2);
   }

   // calculating the log predictive accuracy for each point
   for (n in 1:numObs1) {
     resid_d3[n] = log(CAR_GCN2_counts[n]) - log(CAR_GCN2counts_mean[n]);
     resid_d4[n] = log(CAR_MZN2_counts[n]) - log(CAR_MZN2counts_mean[n]);

     log_lik3[n] = normal_lpdf(log(CAR_GCN2_counts[n]) | log(CAR_GCN2counts_mean[n]), sigma1);
     log_lik4[n] = normal_lpdf(log(CAR_MZN2_counts[n]) | log(CAR_MZN2counts_mean[n]), sigma3);
   }
}
