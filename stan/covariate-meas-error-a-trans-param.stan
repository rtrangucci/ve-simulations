functions {
  real integral(real l, real sp) {
    return log(
      (-1 + l) / 420 *(60 + 60*l^6 + l^5*(60 - 70*sp) - 1029*sp + 3654*sp^2 - 
                               5985*sp^3 + 5180*sp^4 - 2310*sp^5 + 420*sp^6 + 
                               2*l^4*(30 - 77*sp + 42*sp^2) + 
                               l^3*(60 - 259*sp + 294*sp^2 - 105*sp^3) + 
                               l^2*(60 - 399*sp + 714*sp^2 - 525*sp^3 + 140*sp^4) + 
                               l*(60 - 609*sp + 1554*sp^2 - 1785*sp^3 + 980*sp^4 - 210*sp^5)) + (-1 + sp)^6 * sp*log(sp) - (-1 + sp)^6*sp*log(-1 + l + sp)
               );
  }
  real integrand(real x, real xc,
                 array[] real theta,
                 array[] real x_r,
                 array[] int x_i) {
    real lp = 0;
    /* for (i in 1:x_i[1]) { */
    /*   lp += (x_r[2] - 1) * log(x - theta[i]); */
    /* } */
    lp += -x_i[1]
      * (x_r[1] + x_r[2] - 1)
      * log(theta[x_i[1]+1] - 1 + x)
      + (x_r[3] - 1) * log(x)
      + (x_r[4] - 1) * log1m(x);
    return exp(lp);
  }
  real vec_max(array[] real arr) {
    real max_el = arr[1];
    for (i in 2:size(arr)) {
      if (arr[i] > max_el)
        max_el = arr[i];
    }
    return max_el;
  }
}
data {
  int<lower=0> N_r;
  int<lower=0> N_z;
  int<lower=0> N_U;
  int<lower=0> N_S;
  int<lower=0> N_a;
  int<lower=0> N_x;
  // ordering of counts: Y_0 Y_1,S_0 Y_1,S_1
  array[N_r, N_z, N_a, N_x, 4] int<lower=0> xtab; 
  array[N_U, N_z] int<lower=0,upper=1> S1_cond_U_Z;
  // Column 1 is the number of S,
  // Column 2 is the number of groups corresponding to the S
  array[N_S,2] int N_p_S;
  array[N_z] int n_Y_p_z;
  array[N_z] int n_Y0_p_z;
  real sn_S_shape_1;
  real sn_S_shape_2;
  real sp_S_shape_1;
  real sp_S_shape_2;
  real sn_Y_shape_1;
  real sn_Y_shape_2;
  real sp_Y_shape_1;
  real sp_Y_shape_2;
  vector[N_U-1] prior_means;
  vector[N_U-1] prior_sds;
  array[N_x] vector[N_r] p_r_given_x;
  array[N_a] vector[N_x] p_x_given_a;
  array[N_r] vector[N_x] p_x_given_r;
  matrix[N_r,N_x] prop_r_x;
  vector[N_x] prop_x;
  array[N_U] int N_post;
  int N_gt_1;
  /* real sn_S; */
  /* real sn_Y; */
  /* real sp_Y; */
  /* real sp_S; */
  array[N_r] vector[N_U - 1] theta_raw_metric;
  array[N_x - 1] vector[N_U - 1] beta_raw_metric;
  array[N_U] vector[N_a - 1] p_a_raw_metric;
  array[N_x - 1] vector[N_a - 1] beta_a_metric;
  /* array[N_a] vector<lower=0, upper=1>[N_Y_1] p_y_1; */
  real sn_S_raw_metric;
  /* real<lower=0, upper=1> sn_Y_raw; */
  real sp_S_raw_metric;
  real sp_Y_raw_metric;
  array[N_a, sum(n_Y_p_z)] real p_y_1_tilde_raw_metric;
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
  array[N_r] vector[N_U - 1] theta_raw;
  array[N_x - 1] vector[N_U - 1] beta_raw;
  array[N_U] vector[N_a - 1] p_a_raw;
  array[N_x - 1] vector[N_a - 1] beta_a;
  /* array[N_a] vector<lower=0, upper=1>[N_Y_1] p_y_1; */
  real sn_S_raw;
  /* real<lower=0, upper=1> sn_Y_raw; */
  real sp_S_raw;
  real sp_Y_raw;
  array[N_a, N_Y_1] real p_y_1_tilde_raw;
}
transformed parameters {
  real sn_S = 0.5 + 0.5 * inv_logit(sn_S_raw * sn_S_raw_metric);
  /* real sn_Y = sn_Y_raw * 0.5 + 0.5; */
  real sp_Y = 0.5 + 0.5 * inv_logit(sp_Y_raw * sp_Y_raw_metric);
  real sp_S = 0.5 + 0.5 * inv_logit(sp_S_raw * sp_S_raw_metric);
  // ordering of counts: Y_0,S_0; Y_0,S_1; Y_1,S_0; Y_1,S_1
  array[N_r, N_z, N_a, N_x] vector[4] obs_pars;
  array[N_r, N_x] vector[N_U] theta;
  array[N_U, N_x] vector[N_a] p_a;
  array[N_a, N_Y_1] real  p_y_1_tilde;
  /* array[N_a] vector[N_Y_1] p_y_1; */
  for (a in 1:N_a)
    for (n in 1:N_Y_1)
      p_y_1_tilde[a,n] = 1 - sp_Y
        + sp_Y * inv_logit(p_y_1_tilde_raw[a,n] * p_y_1_tilde_raw_metric[a,n]);
  
  for (u in 1:N_U)
    for (x in 1:N_x)
      if (x == 1) {
        p_a[u,x] = softmax(append_row(p_a_raw[u] .* p_a_raw_metric[u],[0]'));
      } else {
        p_a[u,x] = softmax(append_row(p_a_raw[u] .* p_a_raw_metric[u]
                                      + beta_a[x-1] .* beta_a_metric[x-1],[0]'));
      }
  for (r in 1:N_r) {
    for (x in 1:N_x) {
      if (x == 1) {
        theta[r,x] = softmax(append_row(theta_raw[r] .* theta_raw_metric[r],[0]'));
      } else {
        theta[r,x] = softmax(append_row(theta_raw[r] .* theta_raw_metric[r]
                                        + beta_raw[x-1] .* beta_raw_metric[x-1],[0]'));
      }
      for (z in 1:N_z) {
        int idx_start = z == 1 ? 1 : sum(n_Y_p_z[1:(z-1)]) + 1;
        int idx_end = sum(n_Y_p_z[1:z]);
        for (a in 1:N_a) {
          real p_y_0_s_0;
          real p_y_0_s_1;
          real p_y_1_s_1;
          vector[n_Y_p_z[z]] p_y_1_z = to_vector(p_y_1_tilde[a,idx_start:idx_end]);
          vector[n_Y0_p_z[z]] pa_0;
          vector[n_Y_p_z[z]] pa_1;
          for (i in 1:n_Y0_p_z[z])
            pa_0[i] = p_a[idx_0[z][i]][x,a];
          for (i in 1:n_Y_p_z[z])
            pa_1[i] = p_a[idx_1[z][i]][x,a];
          p_y_0_s_0 =
            dot_product(theta[r,x,idx_0[z][1:n_Y0_p_z[z]]], pa_0);
          p_y_0_s_1 =
            dot_product(theta[r,x,idx_1[z][1:n_Y_p_z[z]]] .* pa_1,
                        1 - p_y_1_z);
          p_y_1_s_1  =
            dot_product(theta[r,x,idx_1[z][1:n_Y_p_z[z]]] .* pa_1,
                        p_y_1_z);
          // S = 0, Y = 0
          obs_pars[r,z,a,x,1] = (1 - sn_S) * p_y_0_s_1
             + sp_S * sp_Y * p_y_0_s_0;
          // S = 1, Y = 0
          obs_pars[r,z,a,x,2] = sn_S * p_y_0_s_1
            + (1 - sp_S) * sp_Y * p_y_0_s_0;
          // S = 0, Y = 1
          obs_pars[r,z,a,x,3] = (1 - sn_S) * p_y_1_s_1
            + sp_S * (1 - sp_Y) * p_y_0_s_0;
          // S = 1, Y = 1
          obs_pars[r,z,a,x,4] = sn_S * p_y_1_s_1
            + (1 - sp_S) * (1 - sp_Y) * p_y_0_s_0;
        }
      }
    }
  }
}
model {
  for (r in 1:N_r) {
    for (z in 1:N_z) {
      for (x in 1:N_x) {
        array[N_a * 4] int obs;
        vector[N_a * 4] par_vec;
        for (a in 1:N_a) {
          obs[((a - 1) * 4 + 1):(a * 4)] = xtab[r,z,a,x];
          par_vec[((a - 1) * 4 + 1):(a * 4)] = obs_pars[r,z,a,x];
        }
        obs ~ multinomial(par_vec);
      }
    }
  }

  /* array[N_r] vector[N_U - 1] theta_raw; */
  /* array[N_x - 1] vector[N_U - 1] beta_raw; */
  /* array[N_U] vector[N_a - 1] p_a_raw; */
  /* array[N_x - 1] vector[N_a - 1] beta_a; */
  /* /\* array[N_a] vector<lower=0, upper=1>[N_Y_1] p_y_1; *\/ */
  /* real sn_S_raw; */
  /* /\* real<lower=0, upper=1> sn_Y_raw; *\/ */
  /* real sp_S_raw; */
  /* real sp_Y_raw; */
  /* array[N_a, N_Y_1] real p_y_1_tilde_raw; */
  for (r in 1:N_r)
    theta_raw[r] ~ normal(prior_means, prior_sds ./ theta_raw_metric[r]);
  for (x in 1:(N_x-1))
    beta_raw[x] ~ normal(0, rep_vector(1.7, N_U - 1) ./ beta_raw_metric[x]);
  for (u in 1:N_U)
    p_a_raw[u] ~ normal(0, rep_vector(1.7, N_a - 1) ./ p_a_raw_metric[u]);
  for (x in 1:(N_x-1))
    beta_a[x] ~ normal(0, rep_vector(1.7, N_a - 1) ./ beta_a_metric[x]);

  sn_S_raw ~ normal(0, 1 / sn_S_raw_metric);
  sp_S_raw ~ normal(0, 5 / sp_S_raw_metric);
  /* sn_Y_raw ~ beta(sn_Y_shape_1,sn_Y_shape_2); */
  sp_Y_raw ~ normal(0, 1 / sp_Y_raw_metric);
  to_array_1d(p_y_1_tilde_raw) ~ normal(0,rep_vector(1.7, N_a * N_Y_1) ./ to_vector(to_array_1d(p_y_1_tilde_raw_metric)));
  /* target += integral(vec_max(to_array_1d(p_y_1_tilde)), sp_Y); */
  /* target += log(integrate_1d(integrand, vec_max(to_array_1d(p_y_1_tilde)),1.0, */
  /*                            append_array(to_array_1d(p_y_1_tilde), */
  /*                                         {sp_Y}), */
  /*                            {1, 1, 7, 3}, {N_x * N_a * N_Y_1}, */
  /*                            rel_error */
  /*                           )); */
}
generated quantities {
  /*
    P_S_given_U_Z (rows are S_P_0,
    cols are treatment)
    Can only have Y \in {0,1}
    when S == 1
          [,1] [,2]
    [1,]    0    0
    [2,]    1    0
    [3,]    0    1
    [4,]    1    1
  VE indices for p_y_1 given from
  column-first traversal of matrix above
  Index 1 is p(Y(0) = 1 | S = (1,0))
  Index 2 is p(Y(0) = 1 | S = (1,1))
  Index 3 is p(Y(1) = 1 | S = (0,1))
  Index 4 is p(Y(1) = 1 | S = (1,1))
   */
  matrix[N_a,N_x] ve;
  matrix[N_a, N_z - 1] ve_a;
  matrix[N_r,N_z-1] ve_S;
  matrix[N_x,N_z-1] marg_ve_S_by_x;
  matrix[N_x,N_U] p_S;
  vector[N_z - 1] marg_ve_S;
  vector[N_z - 1] marg_ve_a_x;
  matrix[N_x, N_z - 1] marg_ve;
  array[N_r,N_x,N_z] real p_s_1;
  {
    array[N_a, N_x, N_Y_1] real p_y_1;
    for (a in 1:N_a)
      for (x in 1:N_x)
        for (u in 1:N_Y_1)
          p_y_1[a,x,u] = p_y_1_tilde[a,u] - (1 - sp_Y);

    for (r in 1:N_r)
      for (x in 1:N_x)
        for (z in 1:N_z)
          p_s_1[r,x,z] = sum(theta[r,x,idx_1[z]]);

    for (z in 1:(N_z - 1)) {
      for (a in 1:N_a) {
        for (x in 1:N_x)
          ve[a,x] = 1 - p_y_1[a,x,sum(n_Y_p_z)] / p_y_1[a,x,n_Y_p_z[1]];
        ve_a[a,z] =
          1 -
          dot_product(to_vector(p_y_1[a,,sum(n_Y_p_z[1:(z + 1)])]),
                      p_x_given_a[a])
          / dot_product(to_vector(p_y_1[a,,n_Y_p_z[1]]),
                        p_x_given_a[a]);
      }
    }
    for (r in 1:N_r)
      for (z in 1:(N_z-1))
        ve_S[r,z] = 1 - dot_product(to_vector(p_s_1[r,,z+1]),p_x_given_r[r])
          / dot_product(to_vector(p_s_1[r,,1]),
                        p_x_given_r[r]);
    for (x in 1:N_x) {
      for (z in 1:(N_z - 1))
        marg_ve_S_by_x[x,z] = 1 - dot_product(to_vector(p_s_1[,x,z+1]),
                                          p_r_given_x[x])
        / dot_product(to_vector(p_s_1[,x,1]),p_r_given_x[x]);
    }
    for (x in 1:N_x) {
      for (u in 1:N_U)
        p_S[x,u] = dot_product(to_vector(theta[,x,u]),p_r_given_x[x]);
    }
    for (z in 1:(N_z - 1))
      marg_ve_S[z] = 1 - sum(to_matrix(p_s_1[,,z+1]) .* prop_r_x) / 
      sum(to_matrix(p_s_1[,,1]) .* prop_r_x);

    for (z in 1:(N_z - 1)) {
      vector[N_x] num_ve;
      vector[N_x] denom_ve;
      for (x in 1:N_x) {
        num_ve[x] = dot_product(p_y_1[,x,sum(n_Y_p_z[1:(z+1)])], to_array_1d(p_a[N_U,x]));
        denom_ve[x] =  dot_product(p_y_1[,x,n_Y_p_z[1]], to_array_1d(p_a[N_U,x]));
        marg_ve[x,z] = 1 - num_ve[x] / denom_ve[x];
      }
      marg_ve_a_x[z] = 1 - dot_product(num_ve, prop_x) / dot_product(denom_ve, prop_x);
    }
  }
}
