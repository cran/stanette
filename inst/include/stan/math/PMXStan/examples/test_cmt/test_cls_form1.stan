data{
    int<lower=0> NSUB;                      // number of patients
    int<lower=0> NOBS[NSUB];                // number of observations for each patient
    int<lower=0> NDOSE[NSUB];               // number of doses for each patient
    vector[sum(NOBS)] conc;                 // observations of concentration
    vector<lower=0>[sum(NOBS)] obs_time;    // observation time for all patients
    vector<lower=0>[sum(NDOSE)] dose_time;  // dose time for all patients
    vector<lower=0>[sum(NDOSE)] dose_amt;   // dose amounts for all patients 
    vector<lower=0>[sum(NDOSE)] inf_time;   // infusion time for all doses
}
 
parameters{
    vector<lower=-3.0, upper=3.0>[4] theta;       //fixed effect
    vector<lower=-3.0, upper=3.0>[4] eta[NSUB];   //between-subject random effect
    vector<lower=0, upper=2>[4] sigma2_eta;       //variance of between-subject random effect
    real<lower=0> sigma2;                         //variance of intra-individual random effect
}
 
transformed parameters{
    vector[NSUB] CL;
    vector[NSUB] V;
    vector[NSUB] Q;
    vector[NSUB] V2;

    vector<lower=0>[4] sigma_eta;  
    real<lower=0> sigma;
    vector[sum(NOBS)] y_pred;
 
    sigma_eta[1] <- sqrt(sigma2_eta[1]);
    sigma_eta[2] <- sqrt(sigma2_eta[2]);
    sigma_eta[3] <- sqrt(sigma2_eta[3]);
    sigma_eta[4] <- sqrt(sigma2_eta[4]);

    sigma <- sqrt(sigma2);

    for(n in 1:NSUB){
        CL[n] <- exp( theta[1] + sigma_eta[1] * eta[n,1] );
        V[n] <- exp( theta[2] + sigma_eta[2] * eta[n,2] );
        Q[n] <- exp( theta[3] + sigma_eta[3] * eta[n,3] );
        V2[n] <- exp( theta[4] + sigma_eta[4] * eta[n,4] );

    }
    
    {
        int y_index;
        int dose_index;
        y_index <- 1;
        dose_index <- 1;

        for(i in 1:NSUB){
            vector[NOBS[i]] g;
            vector[4] params;           
            params[1] <- CL[i];
            params[2] <- V[i];
            params[3] <- Q[i];
            params[4] <- V2[i];

            g <- linear_cmpt_iv_infusion(
                     segment(obs_time, y_index, NOBS[i]),
                     segment(dose_time, dose_index, NDOSE[i]),
                     segment(dose_amt, dose_index, NDOSE[i]), 
                     segment(inf_time, dose_index, NDOSE[i]),
                     params, 
                     2,    // number of compartment(s) 
                     1);   // parameterization option: 1(CL_V), 2(micro_rate)

            for(j in 1:NOBS[i])
                y_pred[y_index + j - 1] <- g[j];
            y_index <- y_index + NOBS[i];
            dose_index <- dose_index + NDOSE[i];

        }// end of for loop
    } //end of local variable
}//end of transformed parameters block
 
model{
    for(k in 1:4){
        for(i in 1:NSUB)
            eta[i,k] ~ normal(0.,1.);
        theta[k] ~ normal(0.,1000.);
        sigma2_eta[k] ~ inv_gamma(.01,.01);
    }
    sigma2 ~ inv_gamma(.01,.01);
    conc ~ normal(y_pred, sigma);
}
