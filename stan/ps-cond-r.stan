data {
  int<lower=0> N_r;
  int<lower=0> N_z;
  int<lower=0> N_U;
  int<lower=0> N_S;
  int<lower=0> N_a;
  // ordering of counts: Y_0,S_0 Y_0,S_1 Y_1,S_1
  /* array[N_r, N_z, N_a, 4] int<lower=0> xtab;  */
  array[N_r, N_z, N_a * 4] int<lower=0> obs;
  array[N_U, N_z] int<lower=0,upper=1> S1_cond_U_Z;
  // Column 1 is the number of S,
  // Column 2 is the number of groups corresponding to the S
  array[N_S,2] int N_p_S;
  array[N_z] int n_Y_p_z;
  array[N_z] int n_Y0_p_z;
  /* real<lower=0, upper=1> sn_S; */
  /* real<lower=0, upper=1> sn_Y; */
  /* real<lower=0, upper=1> sp_S; */
  /* real<lower=0, upper=1> sp_Y; */
  vector[N_r] prop_r;
  vector[N_U] prior_counts_theta;
}
transformed data {
  int N_Y_1_z = max(n_Y_p_z);
  array[N_z, N_Y_1_z] int idx_1;
  array[N_z, max(n_Y0_p_z)] int idx_0;
  int N_Y_1 = 0;

  for (z in 1:N_z) {
    int i_1 = 1;
    int i_0 = 1;
    for (u in 1:N_U) {
      if (S1_cond_U_Z[u,z] == 1) {
        idx_1[z,i_1] = u;
        i_1 += 1;
        N_Y_1 += 1;
      } else {
        idx_0[z,i_0] = u;
        i_0 += 1;
      }
    }
  }
  print(idx_0);
  print(idx_1);
}
parameters {
  real<lower=0, upper=1> sn_S_raw;
  real<lower=0, upper=1> sn_Y_raw;
  real<lower=0, upper=1> sp_S_raw;
  real<lower=0, upper=1> sp_Y_raw;
  array[N_r] simplex[N_U] theta;
  array[N_U] vector[N_a - 1] p_a_raw;
  array[N_r - 1] vector[N_a - 1] beta_a;
  array[N_a] vector<lower=0, upper=1>[N_Y_1] p_y_1;
}
transformed parameters {
  // ordering of counts: Y_0,S_0 Y_0,S_1 Y_1,S_0 Y_1,S_1
  real sn_S = sn_S_raw * 0.5 + 0.5;
  real sn_Y = sn_Y_raw * 0.5 + 0.5;
  real sp_Y = sp_Y_raw * 0.5 + 0.5;
  real sp_S = sp_S_raw * 0.5 + 0.5;
  array[N_r, N_z, N_a] vector[4] obs_pars;
  array[N_U, N_r] vector[N_a] p_a;
  for (u in 1:N_U)
    for (r in 1:N_r)
      if (r == 1) {
        p_a[u,r] = softmax(append_row(p_a_raw[u],[0]'));
      } else {
        p_a[u,r] = softmax(append_row(p_a_raw[u]
                                        + beta_a[r-1],[0]'));
      }
  for (r in 1:N_r) {
    for (z in 1:N_z) {
      int idx_start = z == 1 ? 1 : sum(n_Y_p_z[1:(z-1)]) + 1;
      int idx_end = sum(n_Y_p_z[1:z]);
      for (a in 1:N_a) {
        real p_y_0_s_0;
        real p_y_0_s_1;
        real p_y_1_s_1;
        vector[n_Y_p_z[z]] p_y_1_z = p_y_1[a,idx_start:idx_end];
        vector[n_Y0_p_z[z]] pa_0;
        vector[n_Y_p_z[z]] pa_1;
        for (i in 1:n_Y0_p_z[z])
          pa_0[i] = p_a[idx_0[z][i]][r,a];
        for (i in 1:n_Y_p_z[z])
          pa_1[i] = p_a[idx_1[z][i]][r,a];
        p_y_0_s_0 =
          dot_product(theta[r,idx_0[z][1:n_Y0_p_z[z]]], pa_0);
        p_y_0_s_1 =
          dot_product(theta[r,idx_1[z][1:n_Y_p_z[z]]] .* pa_1,
                      1 - p_y_1_z);
        p_y_1_s_1  =
          dot_product(theta[r,idx_1[z][1:n_Y_p_z[z]]] .* pa_1,
                      p_y_1_z);
        obs_pars[r,z,a,1] = (1 - sn_S) * (1 - sn_Y) * p_y_1_s_1
          + (1 - sn_S) * sp_Y * p_y_0_s_1
          + sp_S * sp_Y * p_y_0_s_0;
        obs_pars[r,z,a,2] = sn_S * (1 - sn_Y) * p_y_1_s_1
          + sn_S * sp_Y * p_y_0_s_1
          + (1 - sp_S) * sp_Y * p_y_0_s_0;
        obs_pars[r,z,a,3] = (1 - sn_S) * sn_Y * p_y_1_s_1
          + (1 - sn_S) * (1 - sp_Y) * p_y_0_s_1
          + sp_S * (1 - sp_Y) * p_y_0_s_0;
        obs_pars[r,z,a,4] = sn_S * sn_Y * p_y_1_s_1
          + sn_S * (1 - sp_Y) * p_y_0_s_1
          + (1 - sp_S) * (1 - sp_Y) * p_y_0_s_0;
      }
    }
  }
}
model {
  for (r in 1:N_r) {
    for (z in 1:N_z) {
      /* array[N_a * 4] int obs; */
      vector[N_a * 4] par_vec;
      for (a in 1:N_a) {
        /* obs[((a - 1) * 4 + 1):(a * 4)] = xtab[r,z,a]; */
        par_vec[((a - 1) * 4 + 1):(a * 4)] = obs_pars[r,z,a];
      }
      obs[r,z] ~ multinomial(par_vec);
      
    }
  }
  sp_S_raw ~ beta(10,2);
  sn_S_raw ~ beta(4,2);
  sn_Y_raw ~ beta(5,2);
  sp_Y_raw ~ beta(4,2);
  for (r in 1:N_r)
    theta[r] ~ dirichlet(prior_counts_theta);
  for (u in 1:N_U)
    p_a_raw[u] ~ normal(0, 1.7);
  for (r in 1:(N_r - 1))
    beta_a[r] ~ normal(0, 1.7);

}
generated quantities {
  vector[N_z - 1] marg_ve;
  matrix[N_a, N_z - 1] ve_a;
  matrix[N_r,N_z-1] ve_r;

  for (z in 1:(N_z - 1)) {
    for (a in 1:N_a) {
      ve_a[a,z] = 1 - p_y_1[a,sum(n_Y_p_z[1:(z + 1)])] / p_y_1[a,n_Y_p_z[1]];
    }
  }

  for (z in 1:(N_z - 1)) {
    vector[N_r] num_ve;
    vector[N_r] denom_ve;
    for (r in 1:N_r) {
      num_ve[r] = dot_product(p_y_1[,sum(n_Y_p_z[1:(z+1)])], to_array_1d(p_a[N_U,r]));
      denom_ve[r] =  dot_product(p_y_1[,n_Y_p_z[1]], to_array_1d(p_a[N_U,r]));
      ve_r[r,z] = 1 - num_ve[r] / denom_ve[r];
    }
    marg_ve[z] = 1 - dot_product(num_ve, prop_r) / dot_product(denom_ve, prop_r);
  }
}
