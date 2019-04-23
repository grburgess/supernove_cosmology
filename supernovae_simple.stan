
functions {

  real H(real z, real Om, real w, real h0){
    
    real Ode;
    real zp;
    Ode =  1-Om;
    zp  = 1+z;
   
    return h0* sqrt(pow(zp,3) *Om + Ode )  ;}

  real integrand(real z, real Om, real Ode, real wp1){


    real zp1 = z+1;
    real Hsquare = Om * pow(zp1,3) + Ode *pow(zp1, 3*wp1);
    
    return  inv_sqrt(Hsquare);


  }


  real distance_modulus(real z, real Om, real Ode, real wp1, real h0){

    int N_itr = 10; 
    real N_div = 10.;
    real z0 = 0.; 
    
    real h= (z-z0)/N_div; 
    //vector[N_itr+2] s;
    real s_sum=0.5 * integrand(z0, Om, Ode, wp1);
    real dl;
    //    s[1] = 0.5 * integrand(z0, Om, Ode, wp1);
    
    
    for(n in 1:N_itr)
      {
        s_sum +=  integrand(z0 +n*h, Om, Ode, wp1);
      }
    
    s_sum += 0.5*integrand(z, Om, Ode, wp1);

    s_sum *= h;
    
    dl = (1+z) * (299792458.0*1E-3/ h0) * s_sum;

    return 5*log10(dl)+25;

  }


}



data{
    
  
  int<lower=0> N_sn;           // Number of supernovae
  real<lower=0> h0;
  vector[N_sn] zcmb;           //  redshifts
  vector[N_sn] m_obs;             // Observerd m
  vector[N_sn] m_sigma;             // Observerd m err


  int N_model;
  vector[N_model] z_model;

  
}    
   

parameters {

    
  real<lower=0, upper=1> Om; // yes, I do not have to explain
  real<lower=-20.,upper=-18.> M0; // intrinsic brightness
  real<lower=-3., upper=-0.> tau_log;
  real<lower=-2, upper=-0.0> w;  
    
        
  vector<lower=-10, upper=10>[N_sn] x1_latent;      // 
  vector<lower=-5, upper=5>[N_sn] c_latent;      // 
    
  vector[N_sn] mb;
    
   
  
}



transformed parameters {
    

  real tau;
  real Ode = 1-Om;
  real wp1 = w+1;

  vector[N_sn] dist_mods_latent;
  vector[N_sn] m_latent;
  
  tau = pow(10.,tau_log);

  for(n in 1:N_sn){

    dist_mods_latent[n] = distance_modulus(zcmb[n], Om, Ode, wp1, h0);
 
  }


  m_latent = M0 +  dist_mods_latent;
    
}




model {
    
  mb ~ normal( m_latent, tau);
  m_obs ~ normal(mb, m_sigma);

} 

generated quantities {


  vector[N_model] mb_curve;
  vector[N_model] mu_curve;

  for (n in 1:N_model) {

    mu_curve[n] = distance_modulus(z_model[n], Om, Ode, wp1, h0);


  }

  mb_curve = M0 + mu_curve;


}



