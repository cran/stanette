################################################################################################
# A collection of functions to generate model specific Stan codes based on pre-defined templates
# 
# Code last updated: 6Apr2016 by Yuan Xiong

###########################################################
### template code for popPK model with closed form solution
template_str_pk_cls <- "data{
    int<lower=0> NSUB;                      // number of patients
    int<lower=0> NOBS[NSUB];                // number of observations for each patient
    int<lower=0> NDOSE[NSUB];               // number of doses for each patient
    vector[sum(NOBS)] conc;                 // observations of concentration
    vector<lower=0>[sum(NOBS)] obs_time;    // observation time for all patients
    vector<lower=0>[sum(NDOSE)] dose_time;  // dose time for all patients
    vector<lower=0>[sum(NDOSE)] dose_amt;   // dose amounts for all patients <%= tinf_decl_str %>
}
 
parameters{
    vector<lower=-5, upper=5>[<%= npars %>] theta;       //fixed effect
    vector[<%= npars %>] eta[NSUB];              //inter-individual random effect
    vector<lower=0, upper=2>[<%= npars %>] sigma_eta;     //variance of inter-individual random effect
    real<lower=0> sigma;              //variance of intra-individual random effect
}
 
transformed parameters{
 <%= theta_def_str %>
    vector[sum(NOBS)] y_pred;

    for(n in 1:NSUB){
 <%= theta_calc_str %>
    }
    
    {
        int y_index;
        int dose_index;
        y_index = 1;
        dose_index = 1;

        for(i in 1:NSUB){
            vector[NOBS[i]] g;
            vector[<%= npars %>] params;           
 <%= param_def_str %>
            g = linear_cmpt_<%= rte_str %>(
                     segment(obs_time, y_index, NOBS[i]),
                     segment(dose_time, dose_index, NDOSE[i]),
                     segment(dose_amt, dose_index, NDOSE[i]), <%= tinf_input_str %>
                     params, 
                     <%= ncmpts %>,    // number of compartment(s) 
                     <%= idxpar %>);   // parameterization option: 1(CL_V), 2(micro_rate)

            for(j in 1:NOBS[i])
                y_pred[y_index + j - 1] = g[j];
            y_index = y_index + NOBS[i];
            dose_index = dose_index + NDOSE[i];

        }// end of for loop
    } //end of local variable
}//end of transformed parameters block
 
model{
    for(k in 1:<%= npars %>){
        for(i in 1:NSUB)
            eta[i,k] ~ normal(0.,1.);
        theta[k] ~ normal(0.,1000.);
        sigma_eta[k] ~ normal(0.,1000.);
    }
    sigma ~ normal(0.,1000.);
    conc ~ normal(y_pred, sigma);
}

generated quantities{
    vector[sum(NOBS)] log_lik;
    for(n in 1:sum(NOBS)){
        log_lik[n] = normal_lpdf(conc[n] | y_pred[n], sigma);
    }
}"

###########################################################
### template code for popPK model with ODE solution
template_str_pk_ode <- "data{
    int<lower=0> NSUB;                      // number of patients
    int<lower=0> NOBS[NSUB];                // number of observations for each patient
    int<lower=0> NEVTS[NSUB];               // number of events (dosing and observation) for each patient
    vector[sum(NOBS)] conc;                 // observations of concentration
    vector<lower=0>[sum(NEVTS)] evt_time;   // events (dosing and observation) time for all patients
    vector<lower=0>[sum(NEVTS)] evid;       // events definition
    <%= amt_rate_decl_str %>
    vector[<%= ncmpts %>*NSUB] inits;                   // initial states for ODE solver
}
 
parameters{
    vector<lower=-3.0, upper=3.0>[<%= npars %>] theta;       //fixed effect
    vector<lower=-3.0, upper=3.0>[<%= npars %>] eta[NSUB];   //between-subject random effect
    vector<lower=0, upper=2>[<%= npars %>] sigma2_eta;       //variance of between-subject random effect
    real<lower=0> sigma2;                         //variance of intra-individual random effect
}
 
transformed parameters{
 <%= theta_def_str %>
    vector<lower=0>[<%= npars %>] sigma_eta;  
    real<lower=0> sigma;
    vector[sum(NOBS)] y_pred;
 
 <%= eta_calc_str %>
    sigma = sqrt(sigma2);

    for(n in 1:NSUB){
 <%= theta_calc_str %>}
    
    {
        int y_index;
        int evt_index;
        int init_index;
        y_index = 1;
        evt_index = 1;
        init_index = 1;
        
        for(i in 1:NSUB){
            vector[NOBS[i]] g;
            vector[<%= npars %>] params;           
 
 <%= param_def_str %>         
            g = generic_ode_interface(
	                     params,
	                     segment(inits, init_index, <%= ncmpts %>),
	                     segment(evt_time, evt_index, NEVTS[i]),
	                     segment(evid, evt_index, NEVTS[i]),
	                     segment(<%= dosing_input_str %>, evt_index, NEVTS[i]),  
	                     1E-4,
	                     1E-4,
	                     NOBS[i],
                     1 + <%= ncmpts.aux %>);  // index of state variable to fit

            for(j in 1:NOBS[i])
                y_pred[y_index + j - 1] = g[j]/V[i];
            y_index = y_index + NOBS[i];
            evt_index = evt_index + NEVTS[i];
            init_index = init_index + <%= ncmpts %>;  
        }// end of for loop
    } //end of local variable
}//end of transformed parameters block
 
model{
    for(k in 1:<%= npars %>){
        for(i in 1:NSUB)
            eta[i,k] ~ normal(0.,1.);
        theta[k] ~ normal(0.,1000.);
        sigma2_eta[k] ~ inv_gamma(.01,.01);
    }
    sigma2 ~ inv_gamma(.01,.01);
    conc ~ normal(y_pred, sigma);
}

generated quantities{
    vector[sum(NOBS)] log_lik;
    for(n in 1:sum(NOBS)){
        log_lik[n] = normal_lpdf(conc[n] | y_pred[n], sigma);
    }
}"

###########################################################
### list of PK parameters
pkpar.all <- list(
                 clr = c("CL","V","Q","V2","Q2","V3"), 
                 mcc = c("k0","V","k12","k21","k13","k31"),
                 aux.cls = c("ka","tlag"),
                 aux.ode = "ka"
                )

###########################################################
### template for PK model ODEs
template_pk_ode <- 
    list(
         onecmpt = "k0 = CL/V;
<%= absp_depot_str %>d/dt(centr) = <%= input_def_str %> - k0*centr;
",
         twocmpt = "k0 = CL/V;
k12 = Q/V;
k21 = Q/V2;
<%= absp_depot_str %>d/dt(centr) = <%= input_def_str %> - (k0+k12)*centr + k21*peri;
d/dt(peri)  = k12*centr - k21*peri;
",
         thrcmpt = "k0 = CL/V; 
k12 = Q/V;
k21 = Q/V2;
k13 = Q2/V;
k31 = Q2/V3;
<%= absp_depot_str %>d/dt(centr) = <%= input_def_str %> - (k0+k12+k13)*centr + k21*peri1 + k31*peri2;
d/dt(peri1) = k12*centr - k21*peri1;
d/dt(peri2) = k13*centr - k31*peri2;
"
        )

#########################################################
### Generate Stan model code dynamically for a 
### population PK model 
generateStanCodePK <- function(model.specs)
{
# Input
#    model.specs: model specifications
# Output
#    A Stan source code string according to model specification
#
# Modified 15Jul2015
#   - take model.specs as input
# Modified 5May2015
#   - code generation for PK ODE solution added
#   - for ODE solution, generate ODEs for the PK model automatically
# Created 21Apr2015 by Yuan Xiong for PK closed form solution


m.path <- model.specs$path
m.pk.struct <- model.specs$pk.struct
m.route <- model.specs$route
m.pk.solver <- model.specs$solver
m.pk.param <- model.specs$pk.param

# check/create model path
if (!file.exists(m.path)) {
    dir.create(m.path)
}

# model solver
if(m.pk.solver %in% c("closed_form","ODE")){
  switch(m.pk.solver,
               "closed_form" = {
                   mpksolver <- "cls"
                   template_str <- template_str_pk_cls
               },
               "ODE" = {
                   mpksolver <- "ode"
                   template_str <- template_str_pk_ode
                   generatePKODEs(m.path, m.pk.struct, m.route)
                   ode_str <- paste(readLines(file.path(m.path,"ode.txt")),collapse="\n")
                   instant.stan.extension(ode_str)
                   if(m.route!= "1st_order_abs")
                   { # regenerate ODEs for users
                     generatePKODEs(m.path, m.pk.struct, m.route,1)
                   }
               }
  )
} else stop("Unrecognized PK model solver.")

# PK model structure
if(m.pk.struct %in% c("1-cmpt","2-cmpt","3-cmpt")){
  switch(m.pk.struct,
         "1-cmpt" = {
             ncmpts.base <- 1
             npars.base <- 2 
             mpkstruct <- "1cmpt"
         },
         "2-cmpt" = {
             ncmpts.base <- 2
             npars.base <- 4 
             mpkstruct <- "2cmpt"
         },
         "3-cmpt" = {
             ncmpts.base <- 3
             npars.base <- 6 
             mpkstruct <- "3cmpt"
         }
  )
} else stop("Unrecognized PK model structure.")

# drug administration
if(m.route %in% c("1st_order_abs","IV_bolus","IV_infusion")){
  switch(m.route,
         "1st_order_abs" = {
             switch(m.pk.solver,
                    "closed_form" = {
                        ncmpts.aux <- 0
                        npars.aux <- 2
                        par.aux <- pkpar.all$aux.cls
                    },
                    "ODE" = {
                        ncmpts.aux <- 1
                        npars.aux <- 1
                        par.aux <- pkpar.all$aux.ode
                    }
             )
             rte_str <- "1order_absor"      # cls
             tinf_decl_str <- NULL          # cls
             tinf_input_str <- NULL         # cls
             amt_rate_decl_str <- "vector<lower=0>[sum(NEVTS)] dose_amt;   // dose amounts for all patients "      # ODE
             dosing_input_str <- "dose_amt" # ODE
             mpkadmin <- "1stabsorp"
         },
         "IV_bolus" = {
             ncmpts.aux <- 0
             npars.aux <- 0
             par.aux <- NULL
             rte_str <- "iv_bolus"          # cls
             tinf_decl_str <- NULL          # cls
             tinf_input_str <- NULL         # cls
             amt_rate_decl_str <- "vector<lower=0>[sum(NEVTS)] dose_amt;   // dose amounts for all patients "      # ODE
             dosing_input_str <- "dose_amt" # ODE
             mpkadmin <- "ivbolus"
         },
         "IV_infusion" = {
             ncmpts.aux <- 0
             npars.aux <- 0
             par.aux <- NULL
             rte_str <- "iv_infusion"										 # cls
             tinf_decl_str <- "\n    vector<lower=0>[sum(NDOSE)] inf_time;   // infusion time for all doses"	 # cls
             tinf_input_str <-"\n                     segment(inf_time, dose_index, NDOSE[i])," 			 # cls
             amt_rate_decl_str <- "vector[sum(NEVTS)] dose_rate;           // infusion rates for all patients"	 # ODE
             dosing_input_str <- "dose_rate"									 # ODE
             mpkadmin <- "ivinfs"
         }
  )
} else stop("Unrecognized drug administration.")

# model parameterization
if(m.pk.param %in% c("CL_V","micro_rate")){
  switch(m.pk.param,
         "CL_V" = {
             idxpar <- 1
             mpkparam <- "clearance"
         },
         "micro_rate" = {
             idxpar <- 2
             mpkparam <- "microconst"
         }
  )
} else stop("Unrecognized PK model parameterization.")

# get model parameter vector
ncmpts <- ncmpts.base + ncmpts.aux
npars <- npars.base + npars.aux
par.vec <- c(pkpar.all[[idxpar]][1:npars.base], par.aux)  

# get strings to be defined in the template
theta_def_lst <- lapply(1:npars, function(ix)
                    sprintf("   vector[NSUB] %s;\n", par.vec[ix])
                )
theta_def_str <- do.call("paste", theta_def_lst)

eta_calc_lst <- lapply(1:npars, function(ix)
                   sprintf("   sigma_eta[%d] = sqrt(sigma2_eta[%d]);\n", ix, ix)
               )
eta_calc_str <- do.call("paste", eta_calc_lst)

theta_calc_lst <- lapply(1:npars, function(ix)
                     sprintf("       %s[n] = exp( theta[%d] + sigma_eta[%d] * eta[n,%d] );\n", par.vec[ix], ix, ix, ix)
                 )
theta_calc_str <- do.call("paste", theta_calc_lst)

param_def_lst <- lapply(1:npars, function(ix)
	            sprintf("           params[%d] = %s[i];\n", ix, par.vec[ix])
                )
param_def_str <- do.call("paste", param_def_lst)

# save stan code
stanfilename <- paste("popPK_",mpkstruct,"_",mpkadmin,"_",mpkparam,"_",mpksolver,".stan",sep="")
brew(text = template_str, output = file.path(m.path, stanfilename))

list(stanfilename = stanfilename, ntheta = npars, neta = npars)
}

 
###############################################################
### Generate ODE system automatically for a population PK model 
generatePKODEs <- function(m.path, m.pk.struct, m.route, ode.mode=0)
{
# Input
#    as defined in the S3 class PMXStanModel
#    ode.mode: whether or not (1/0) to include infusion_rate of dose_amt into the ODEs
# Output
#    A string of a set of ODEs according to PK model specification
#
# modified 14May2015
#   - generate 2 sets of ODEs for the IV_bolus and IV_infusion cases
# Created 5May2015 by Yuan Xiong

# check/create model path
if (!file.exists(m.path)) {
    dir.create(m.path)
}

switch(m.pk.struct,
       "1-cmpt" = {
           template_str <- template_pk_ode$onecmpt
       },
       "2-cmpt" = {
           template_str <- template_pk_ode$twocmpt
       },
       "3-cmpt" = {
           template_str <- template_pk_ode$thrcmpt
       }
)

switch(m.route,
       "1st_order_abs" = {
           absp_depot_str <- "d/dt(depot) =-ka*depot;\n"
           input_def_str <- "ka*depot"
       },
       "IV_bolus" = {
           absp_depot_str <- NULL
           input_def_str <- NULL
           if(ode.mode>0) input_def_str <- "DoseAmount"
       },
       "IV_infusion" = {
           absp_depot_str <- NULL
           input_def_str <- NULL
           if(ode.mode>0) input_def_str <- "InfusionRate"
       }
)

brew(text = template_str, output = file.path(m.path, "ode.txt"))

}

#########################################################
# PKPD models in ODE form

template_str_pkpd <- "data{
    int<lower=0> NSUB;                      // number of patients
    int<lower=0> NOBS[NSUB];                // number of observations for each patient
    int<lower=0> NEVTS[NSUB];               // number of events (dosing&observation) for each patient
    vector[sum(NOBS)] y;                    // observations of response
    vector<lower=0>[sum(NEVTS)] evt_time;   // events (dosing&observation) time for all patients
    vector<lower=0>[sum(NEVTS)] evid;       // events definition
    <%= amt_rate_decl_str %>
    vector[<%= nstates %>*NSUB] inits;                   // initial states for ODE solver
}

parameters{
    vector<lower=-5, upper=5>[<%= ntheta %>] theta;       
    <%= rnd_eff_decl_str %>
    real<lower=0> sigma;                         //variance of intra-individual random effect
}

transformed parameters{
    <%= pars_decl_str %>

    vector[sum(NOBS)] y_pred;

    for(n in 1:NSUB){
        <%= pars_calc_str %>
    }
    
    {
        int y_index;
        int evt_index;
        int init_index;
        y_index    = 1;
        evt_index  = 1;
        init_index = 1;

        for(i in 1:NSUB){
            vector[NOBS[i]] g;
            vector[<%= npars %>] params;
            <%= pars_spec_str %>
 
            g = generic_ode_interface(
                params,
                segment(inits, init_index, <%= nstates %>),
                segment(evt_time, evt_index, NEVTS[i]),
                segment(evid, evt_index, NEVTS[i]),
                segment(<%= dosing_input_str %>, evt_index, NEVTS[i]),  
                1E-4,
                1E-4,
                NOBS[i],
                <%= m.obs.state %>);  // index of state variable to fit
            for(j in 1:NOBS[i])
                y_pred[y_index + j - 1] = g[j];
            y_index = y_index + NOBS[i];
            evt_index = evt_index + NEVTS[i];
            init_index = init_index + <%= nstates %>;
        } //end of for loop
    } //end of local variable
}//end of transformed parameters block

model{
    for(k in 1:<%= ntheta %>){
        theta[k] ~ normal(0.,1000.);
    }
    <%= rnd_eff_model_str %>
    sigma ~ normal(0.,1000.);
    y ~ normal(y_pred, sigma);
}

generated quantities{
    vector[sum(NOBS)] log_lik;
    for(n in 1:sum(NOBS)){
        log_lik[n] = normal_lpdf(y[n] | y_pred[n], sigma);
    }
}"

#########################################################
### Generate Stan model code dynamically for a PKPD model 
generateStanCodePKPD <- function(model.specs)
{
# Input
#    model.specs: model specifications
# Output
#    A Stan source code string according to model specification
#
# Modified 20Jan2016
# Created 8Jul2015 by Yuan Xiong 

m.path <- model.specs$path
m.route <- model.specs$route
m.obs.state <- model.specs$obs.state
m.theta <- model.specs$theta
m.eta <- model.specs$eta
m.const <- model.specs$fixed

# check/create model path
if (!file.exists(m.path)) {
    dir.create(m.path)
}

### get input
pars <- scan(sprintf("%s/ODE_PARS.txt",tempdir()), character(0), quiet = T)
# total number of parameters in the ODE
npars <- length(pars)

# number of parameters to be estimated
ntheta <- length(m.theta) 
# number of parameters with between-subject variability
neta <- length(m.eta)

# get theta
idx.theta <- sapply(m.theta, function(theta) match(theta, pars))

# get constant
if (!is.null(m.const)) {
  const.par.name <- names(m.const)
}

# get eta
flag.eta <- sapply(m.theta, function(theta) as.numeric(theta%in%m.eta))
count.eta <- cumsum(flag.eta)

### dynamic strings in STAN code
# number of state variables
states <- scan(sprintf("%s/STATE_VARS.txt",tempdir()), character(0), quiet = T)
nstates <- length(states)

# amt_rate_decl_str and dosing_input_str
# drug administration
if(m.route %in% c("1st_order_abs","IV_bolus","IV_infusion")){
  if (m.route == "IV_infusion") {
    amt_rate_decl_str <- "vector[sum(NEVTS)] dose_rate;           // infusion rates for all patients"	 # ODE
    dosing_input_str <- "dose_rate"
  } else {
    amt_rate_decl_str <- "vector<lower=0>[sum(NEVTS)] dose_amt;   // dose amounts for all patients "      # ODE
    dosing_input_str <- "dose_amt" 
  }
} else stop("Unrecognized drug administration.")

# pars_decl_str
pars_decl_lst <- lapply(1:ntheta, function(ix)
                    if (ix==1) sprintf("vector[NSUB] %s; ", names(idx.theta)[ix])
                    else sprintf("\n    vector[NSUB] %s; ", names(idx.theta)[ix])
                )
pars_decl_str <- do.call("paste", pars_decl_lst)

# pars_calc_str
pars_calc_lst <- lapply(1:ntheta, function(ix) {
                    if (flag.eta[ix]) {
                      if (ix==1) sprintf("%s[n] = exp( theta[%d] + sigma_eta[%d] * eta[n,%d] ); ", names(idx.theta)[ix], ix, count.eta[ix], count.eta[ix]) 
                      else sprintf("\n        %s[n] = exp( theta[%d] + sigma_eta[%d] * eta[n,%d] ); ", names(idx.theta)[ix], ix, count.eta[ix], count.eta[ix]) 
                    } else {
                      if (ix==1) sprintf("%s[n] = exp( theta[%d] );", names(idx.theta)[ix], ix)
                      else sprintf("\n        %s[n] = exp( theta[%d] ); ", names(idx.theta)[ix], ix)
                    }
                })
pars_calc_str <- do.call("paste", pars_calc_lst)

# pars_spec_str
pars_spec_lst <- lapply(1:npars, function(ix) {
                    if (ix%in%idx.theta) {  # to be estimated
                      if (ix==1) sprintf("params[%d] = %s[i]; ", ix, pars[ix]) 
                      else sprintf("\n            params[%d] = %s[i]; ", ix, pars[ix]) 
                    } else { # constant
                      if (ix==1) sprintf("params[%d] = %g; ", ix, m.const[const.par.name==pars[ix]])
                      else sprintf("\n            params[%d] = %g; ", ix, m.const[const.par.name==pars[ix]])
                    }
                })
pars_spec_str <- do.call("paste", pars_spec_lst)

# rnd_eff_decl_str, sigma_eta_decl_str, sigma_eta_calc_str, rnd_eff_model_str
if (neta>0){  # at least one eta
  rnd_eff_decl_str <- sprintf("vector[%d] eta[NSUB]; \n    vector<lower=0, upper=2>[%d] sigma_eta; ", neta, neta)
  rnd_eff_model_str <- sprintf("for(k in 1:%d){ \n        for(i in 1:NSUB) \n            eta[i,k] ~ normal(0.,1.); \n        sigma_eta[k] ~ normal(0.,1000.); \n    } ", neta)

} else { # no eta
  rnd_eff_decl_str <- NULL
  rnd_eff_model_str <- NULL
}

# save stan code
stanfilename <- "PKPD.stan"
brew(text = template_str_pkpd, output = file.path(m.path, stanfilename))

list(stanfilename = stanfilename, ntheta = ntheta, neta = neta)

}
