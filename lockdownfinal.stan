data {
  int <lower=1> M; // number of Rt prediction models
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  real<lower=0> x[N2]; // index of days (starting at 1)
  int deaths[N2]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  vector[N2] f; // h * s
  #int EpidemicStart[M];
  real SI[N2]; // fixed pre-calculated SI using emprical data from Neil
}

#transformed data {
 # real delta = 1e-5;
#}

parameters {
  real<lower=0> r0; //  initial R
  real<lower=0> rt[9]; // r varying with time after lockdown
 # real<lower=0> kappa;
  real<lower=0> y;
  real<lower=0> phi;
  real<lower=0> tau;
}

transformed parameters {
   # int q;
    real convolution;
    vector[N2] prediction = rep_vector(0,N2);
    vector[N2] E_deaths  = rep_vector(0,N2);
    vector[N2] Rt = rep_vector(0,N2);
   # vector[N2] Rt = rep_vec(0,N2);
      prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days

      for (i in 1:36) {
      #  q=i%10+1;
        Rt[i]=r0;
      }

      for (i in 37:43) {
        Rt[i]=rt[1];
      }
     
     for (i in 44:50) {
        Rt[i]=rt[2];
      }

   for (i in 51:57) {
        Rt[i]=rt[3];
      }

   for (i in 58:64) {
        Rt[i]=rt[4];
      }

    for (i in 65:71) {
        Rt[i]=rt[5];
      }

    for (i in 72:78) {
        Rt[i]=rt[6];
      }

     for (i in 79:85) {
        Rt[i]=rt[7];
      }

      for (i in 86:94) {
        Rt[i]=rt[8];
      }    
           for (i in 95:N2) {
        Rt[i]=rt[9];
      }

      for (i in (N0+1):N2) {
        convolution=0;
        for(j in 1:(i-1)) {
          convolution += prediction[j]*SI[i-j]; // Correctd 22nd March
        }
        prediction[i] = Rt[i] * convolution;
      }
      
      E_deaths[1]= 1e-9;
      for (i in 2:N2){
        E_deaths[i]= 0;
        for(j in 1:(i-1)){
          E_deaths[i] += prediction[j]*f[i-j];
        }
      }
}


model {
  tau ~ exponential(0.03);
      y ~ exponential(1.0/tau);
  phi ~ normal(0,5);
  #kappa ~ normal(0,0.5);
  r0 ~ normal(2.4, 0.5); // citation needed 
  rt ~ normal(0.8, 0.25);

    for(i in 1:N){
       deaths[i] ~ neg_binomial_2(E_deaths[i],phi); 
    }
}

generated quantities {
    vector[N] lp0 = rep_vector(1000,N); // log-probability for LOO for the counterfactual model
    vector[N] lp1 = rep_vector(1000,N); // log-probability for LOO for the main model
    vector[M] diffcsdeaths81 = rep_vector(0,M);
    vector[M] diffcscases81 = rep_vector(0,M);
    vector[M] diffcsdeaths82 = rep_vector(0,M);
    vector[M] diffcscases82 = rep_vector(0,M);
    vector[M] diffcsdeaths83 = rep_vector(0,M);
    vector[M] diffcscases83 = rep_vector(0,M);
    vector[M] diffcsdeaths84 = rep_vector(0,M);
    vector[M] diffcscases84 = rep_vector(0,M);
    vector[M] cscases = rep_vector(0,M);
    vector[M] csdeaths = rep_vector(0,M);
    vector[2] q=rep_vector(0,2);
    real convolution0;
    matrix[N2,M] prediction0 = rep_matrix(0,N2,M);
    matrix[N2,M] E_deaths0  = rep_matrix(0,N2,M);
    matrix[N2,M] Rt0 = rep_matrix(0,N2,M);
   real r;
   int m;
   int t;
   int u;

    r=0.60 ;
    q[1]=0;
    q[2]=0.05;
    m=1;
    t=1;
    u=1;
while (r<1.15){
 for (n in q) {
  for (p in q) {
   for (i in 1:113) {
     Rt0[i,m]=Rt[i];
     }
     for (i in 114:127) {
     Rt0[i,m]=r;
      }

      for (i in 128:145) {
       Rt0[i,m]=r+n;
      }
      
    for (i in 146:N2) {
    Rt0[i,m]=r+n+p;
    }
  if (t==1) { 
  m=m+1;
  }
  t=t+1;
}
    t=1;
  if (u==1) { 
  m=m+1;
  }
  u=u+1;
}
  u=1;
  m=m+1;
  r=r+0.05;
}



Rt0[1:113,45] = Rt[1:113];
Rt0[114:N2, 45]=rep_vector(0.81,N2-113);
Rt0[1:113,46] = Rt[1:113];
Rt0[114:N2, 46]=rep_vector(0.82,N2-113);
Rt0[1:113,47] = Rt[1:113];
Rt0[114:N2, 47]=rep_vector(0.83,N2-113);
Rt0[1:113,48] = Rt[1:113];
Rt0[114:N2, 48]=rep_vector(0.84,N2-113);




   for (v in 1:M) {
      prediction0[1:N0,v] = rep_vector(y,N0); // learn the number of cases in the first N0 days
      for (i in (N0+1):N2) {
        convolution0=0;
        for(j in 1:(i-1)) {
          convolution0 += prediction0[j,v]*SI[i-j]; // Correctd 22nd March
        }
        prediction0[i,v] = Rt0[i,v] * convolution0;
      }
      
      E_deaths0[1,v]= 1e-9;
      for (i in 2:N2){
        E_deaths0[i,v]= 0;
        for(j in 1:(i-1)){
          E_deaths0[i,v] += prediction0[j,v]*f[i-j];
        }
      }

}

      for(i in 1:N){
        lp0[i] = neg_binomial_2_lpmf(deaths[i] | E_deaths[i],phi); 
        lp1[i] = neg_binomial_2_lpmf(deaths[i] | E_deaths0[i,1],phi); 
      }

    for (i in 1:M) {
      csdeaths[i]=sum(E_deaths0[,i]);
      cscases[i]=sum(prediction0[,i]);

      diffcsdeaths81[i]=sum(E_deaths0[,i])-sum(E_deaths0[,45]);
      diffcscases81[i]=sum(prediction0[,i])-sum(prediction0[,45]);
    diffcsdeaths82[i]=sum(E_deaths0[,i])-sum(E_deaths0[,46]);
      diffcscases82[i]=sum(prediction0[,i])-sum(prediction0[,46]);
      diffcsdeaths83[i]=sum(E_deaths0[,i])-sum(E_deaths0[,47]);
      diffcscases83[i]=sum(prediction0[,i])-sum(prediction0[,47]);
      diffcsdeaths84[i]=sum(E_deaths0[,i])-sum(E_deaths0[,48]);
      diffcscases84[i]=sum(prediction0[,i])-sum(prediction0[,48]);


    }

}


