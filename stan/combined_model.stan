data {
  int<lower=0> d;                  // number of days
  int<lower=0> r;                  // number of regions
  int<lower=0> w;                  // size of the convolution window
  array[d, r] int<lower=0> H;      // hospitalisations
  array[d, r] int<lower=0> D;      // deaths
  
  array[r] real<lower=0> I0;       // PCR-positive cases on 01-08-2020
  
  array[w] real Pi;                // infection-to-hospitalisation delay, Pi[1] corresponds to t=1
  array[w] real Tau;               // infection-to-death delay, Pi[1] corresponds to t=1
  array[w] real g;                 // generation interval, g[1] corresponds to t=1
}

parameters {
  real<lower=0> phi_h;             //neg-binomial over-dispersion (hospitalisations)
  real<lower=0> phi_d;             //neg-binomial over-dispersion (deaths)
  real<lower=0, upper=1> alpha;    //infection-to-hospitalisation probability
  real<lower=0, upper=1> ifr;      //infection-to-death probability
  real<lower=0> sigma;             //random walk typical step size
  
  array[r] real<lower=0> i0;
  array[r] real<lower=0> R0;
  
  array[d, r] real rdm_walk;
}

transformed parameters {
  array[d+w, r] real R;               // effective reproduction number
  array[d+w, r] real incidence;       // incidence, incidence[1] corresponds to t=1
  array[d, r] real h;               // E[hospitalisations]
  array[d, r] real deaths;          // E[deaths]

  for(i in 1:d+w) {
    for(j in 1:r) {
      if (i <= w) {
        R[i,j] = R0[j];
        incidence[i, j] = i0[j];
      } else {
        R[i,j] = R0[j] * exp(rdm_walk[i-w, j]);
        incidence[i, j] = R[i, j] * (dot_product(incidence[i-w:i-1, j], reverse(g)));
        h[i-w, j] = alpha * (dot_product(incidence[i-w:i-1, j], reverse(Pi)));
        deaths[i-w, j] = ifr * (dot_product(incidence[i-w:i-1, j], reverse(Tau)));
      }
    }
  }
}

model {
  phi_h ~ gamma(6.25, 0.125);
  phi_d ~ normal(0., 5.);
  alpha ~ normal(0.028, 0.002);
  ifr   ~ normal(0.0054, 0.002);
  sigma ~ normal(0.3, 0.02);
  
  for(i in 1:r){
    R0[i] ~ normal(1, 0.1);
    i0[i] ~ exponential(1./(3.0*I0[i]));
  }
  
  for(i in 1:d)
    for(j in 1:r) {
      if (i==1)
        rdm_walk[i,j] ~ normal(0., sigma);
      else
        rdm_walk[i,j] ~ normal(rdm_walk[i-1, j], sigma);
      H[i, j] ~ neg_binomial_2(h[i, j], phi_h);
      D[i, j] ~ neg_binomial_2(deaths[i, j], phi_d);
  }
}