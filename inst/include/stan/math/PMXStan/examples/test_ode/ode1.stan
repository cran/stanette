data{
    int<lower=0> NOBS;                 // number of observations for each patient
    int<lower=0> NEVTS;                // number of observations for each patient
    vector[NOBS] conc;                 // observations of concentration
    vector[2] inits;                   // 
    vector<lower=0>[NEVTS] obs_time;   // observation time for all patients
    vector<lower=0>[NEVTS] evid;       // 
    vector<lower=0>[NEVTS] amt;        // 
}

parameters{
    vector<lower=-3.0, upper=3.0>[3] theta; //fixed effect
    real<lower=0> sigma2;
}

transformed parameters{
    real<lower=0> ka;
    real<lower=0> ke;
    real<lower=0> V;
    real<lower=0> sigma;
    vector[NOBS] y_pred;

    ka <- exp( theta[1] );
    ke <- exp( theta[2] );
    V  <- exp( theta[3] );
    sigma <- sqrt(sigma2);
    
    {
        vector[NOBS] g;
        vector[3] params;
        params[1] <- ka;
        params[2] <- ke;
        params[3] <- V;
        g <- generic_ode_interface(
            params,
            inits,
            obs_time,
            evid,
            amt,  
            1E-4,
            1E-4,
            NOBS,
            2);
        for(j in 1:NOBS)
            y_pred[j] <- g[j]/V;
    } //end of local variable
}//end of transformed parameters block

model{
    for(k in 1:3){
        theta[k] ~ normal(0.,1000.);
    }
    sigma2 ~ inv_gamma(.01,.01);
    conc ~ normal(y_pred, sigma);
}
