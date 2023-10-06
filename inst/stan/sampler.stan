data {
  int<lower=0> N;
  vector[N] completeness;
  real betap1;
  real betap2;
  real core_threshold;
  real rare_threshold;
}

transformed data {
  real cdf_core = beta_cdf(core_threshold | betap1, betap2);
  real cdf_rare = beta_cdf(rare_threshold | betap1, betap2);
}

generated quantities {
  real prior_sample = beta_rng(betap1, betap2);

  real u_core = uniform_rng(cdf_core, 1);
  real u_notcore = uniform_rng(0, cdf_core);
  real prior_core = inv_inc_beta(betap1, betap2, u_core);
  real prior_notcore = inv_inc_beta(betap1, betap2, u_notcore);

  real u_rare = uniform_rng(0, cdf_rare);
  real u_notrare = uniform_rng(cdf_rare, 1);
  real prior_rare = inv_inc_beta(betap1, betap2, u_rare);
  real prior_notrare = inv_inc_beta(betap1, betap2, u_notrare);

  vector[N] obs_vec_core;
  vector[N] obs_vec_notcore;
  vector[N] obs_vec_rare;
  vector[N] obs_vec_notrare;

  for(i in 1:N){
    obs_vec_rare[i] = bernoulli_rng(prior_rare * completeness[i]);
    obs_vec_notrare[i] = bernoulli_rng(prior_notrare * completeness[i]);
    obs_vec_core[i] = bernoulli_rng(prior_core * completeness[i]);
    obs_vec_notcore[i] = bernoulli_rng(prior_notcore * completeness[i]);
  }

  real obs_core = sum(obs_vec_core);
  real obs_notcore = sum(obs_vec_notcore);
  real obs_rare = sum(obs_vec_rare);
  real obs_notrare = sum(obs_vec_notrare);
}
