
data {
  int N;
  int N_ex;
  int N_each;
  vector[N] ex;
  matrix[N,N_each] Y;
  real obs;
  array[N_each] int start;
  int boundary1;
  int boundary2;
}


parameters {
  vector[N] w;
  vector[N_ex] b_ex;
  real<lower=0> s_w;
  real<lower=0> s_b_ex;
}


transformed parameters {
  matrix[N_ex,N_each] b_ex_each;
  matrix[N,N_each] alpha_each;

  // Insert 0 to b_ex_each
  for (n in 1:N_each) {

    if (start[n] >= 3) {
      for (t in 1:(start[n]-2)) {
        b_ex_each[t,n] = 0;
      }
      for (n in 1:N_each) {
        for (t in (start[n]-1):N_ex) {
          b_ex_each[t,n] = b_ex[t-(start[n]-1)+1];
        }
      }
    }//if (start[n] >= 3) {

    if (start[n] == 1 | start[n] == 2){
      for (t in 1:N_ex) {
        b_ex_each[t,n] = b_ex[t];
      }
    }//if (start[n] == 1 | start[n] == 2){

  }//for (n in 1:N_each)

  // True state of velocity of each response variable
  for (n in 1:N_each) {
    for (t in 1:boundary1) {
      alpha_each[t,n] = w[t];
    }
    if (boundary2 <= N) {
      for (t in boundary2:N) {
        alpha_each[t,n] = w[t];
      }
    }
    for (t in (boundary1+1):(boundary2-1)) {
      alpha_each[t,n] = w[t] + b_ex_each[t-boundary1,n] * ex[t];
    }
  }

}


model {
  // State equation of w
  for (t in 1:N) {
    w[t] ~ normal(0, s_w);
  }

  // Coefficient of ex (common)
  b_ex[1] ~ normal(0, s_b_ex);
  for (t in 2:N_ex) {
    b_ex[t] ~ normal(b_ex[t-1], s_b_ex);
  }

  // Observation equation of Y
  for (n in 1:N_each) {
    for (t in 1:N) {
      Y[t,n] ~ normal(alpha_each[t,n], obs);
    }
  }
}


generated quantities {
  vector[N] alpha;
  vector[N] dist;

  // Common state of velocity
  for (t in 1:boundary1) {
    alpha[t] = w[t];
  }
  if (boundary2 <= N) {
    for (t in boundary2:N) {
      alpha[t] = w[t];
    }
  }
  for (t in (boundary1+1):(boundary2-1)) {
    alpha[t] = w[t] + b_ex[t-boundary1] * ex[t];
  }

  // Common state of distance
  for (t in 1:N) {
    dist[t] = sum(alpha[1:t]);
  }
  dist = dist - sum(alpha[1:boundary1]);
}
