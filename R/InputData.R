#' @title Transformation of a NONMEM-readable dataset to a Stan-readable list
#'
#' @description 
#' Prepares input data (a list) for the model-specific Stan code from a conventional PK/PD dataset for NONMEM 
#' (usually a table in .txt or .csv format). 
#'
#' @details
#' \code{prepareInputData} rewrites a table consisting of relevant columns according to NONMEM conventions to a \code{list} 
#' consisting of named vectors compatible to the auto-generated model-specific Stan code. It also allows users to specify 
#' column names for initial values of the observed state variable for an ODE system, and column names for relevant 
#' covariates that will be explored more through modeling or post-processing.
#'
#' Currently, \code{prepareInputData} takes two different formats (.txt or .csv) of a NONMEM data file. The minimally 
#' required columns in this data file are: \code{TIME}, \code{DV}, \code{AMT},and \code{EVID}. If \code{data.type} is set 
#' to "population", a column \code{ID} is also required. If the route of drug administration is set to "IV_infusion" 
#' (through argument \code{route} for \code{PMXStanModel}), a column \code{RATE} is also required.
#' 
#' For more details on the data structure of the returned list, please see \emph{Value} section.
#' 
#' @param data.type a string to specify the type of the data ("population" or "individual").
#' @param data.source the name of an R data frame that serves as the source data. If not provided, the next argument
#'                    \code{data.file} will be checked to locate the data in the system.
#' @param data.file a string to specify the path of the data file. It will be ignored if the previous argument 
#'                  \code{data.source} is already provided; but it will be required if \code{data.source} is not provided.
#' @param model a \code{PMXStanModel} object that will be fitted to the current dataset.
#' @param inits a string providing the colume name of initial values of the observed state variable for each individual 
#'              in the current dataset, or a numerical vector of initial values for each individual.
#' @param covar a string or a vector of strings providing colume name(s) of covariate(s) in the current dataset, 
#'              which needs to be used either in model building or post-processing (e.g. goodness-of-fit plotting
#'              \code{\link{obs.vs.pred}} and \code{\link{rsd.vs.pred}} for a \code{\link{PMXStanFit}} object), or both.
#'
#' @return 
#' The returned list contains four elements:
#' \item{data.type}{a string indicating the data type for future reference; the same as the argument \code{data.type}.}
#' \item{standata}{a list by itself to be fed into the model-specific Stan code.}
#' \item{ID}{a vector of IDs whose observations have been included in \code{standata} above. Note that a proof checking
#'     step has been performed to make sure that every subject included in Stan input data has both valid dosing and
#'     observation records; therefore, this set of \code{ID}'s might be the same as those in the original data table,
#'     but it might also be a subset only. The purpose of this step is to avoid errors at the future sampling step due
#'     to lack of dosing records or observations.
#' }
#' \item{d.cov}{a data frame with columns of subject IDs and values for each selected covariates defined by the \code{covar} 
#      argument. Note that for now it only handles \emph{time-invariant} covariates. 
#' }
#'
#' A special consideration should be noted when handling data from a single subject, i.e. when \code{data.type} is set to 
#' "individual". Due to the strict requirement with data type specification in Stan grammar and the underlying mechanisms 
#' how C++ handles different data types, to avoid the comlexity caused by data from a single subject, a "dummy" subject
#' has been attached in these cases, with a single dose of 0 and the initial observation set the same as the real subject.
#'
#' @seealso 
#' \code{\link[base]{list}} for general information on \code{R} lists; 
#' \code{\link{PMXStanModel}} for the initialization of a \code{PMXStanModel} object;
#' \code{\link{PMXStanFit}} for the generation of a \code{PMXStanFit} object, running a model-specific Stan executable 
#' and linking to an input data list to generate posterior samples of parameters.
#' 
#' @author Yuan Xiong and Wenping Wang
#'
#' @examples
#' m <- PMXStanModel(path = "pk_m1")
#' data("examples_data")
#' dat <- prepareInputData(data.source = d1_nm_poppk, model = m)
#' str(dat)
#
# Documentation last updated: 6Apr2016
#
# Code last updated: 2Jun2016 by Yuan Xiong
# Modified 18Jan2016
#   - included user-specified initis and covariates
# Modified 8Jul2015
#   - handled non-zero initial values for ODE models
# Modified 23Apr2015 
#   - excluded patients with ndose=0 or nobs=0
#   - handled individual data input
# Created 2Apr2015 by Yuan Xiong

utils::globalVariables(c("EVID", "sel.ds", "sel.obs"))
prepareInputData <- function(data.type = "population",
                             data.source = NULL, 
                             data.file = NULL,
                             model = NULL,
                             inits = NULL,
                             covar = NULL
                            ) 
{
###########################
# Check input arguments

# check data type
if(!hasArg(data.type)) 
  message("Data type is set to \"population\" by default.")
if(!data.type %in% c("individual","population"))
  stop("Unrecognized data type.")

# check data source/file
if(!hasArg(data.source)){
  if(!hasArg(data.file)) 
    stop("Either source data (data.source) or the location of the source data (data.file) needs to be provided.")
} else {
  if(hasArg(data.file)) 
    warning("Since source data (data.source) has been provided, the location of the data (data.file) will be ignored.")
}

# check model object
if(!hasArg(model)) 
  stop("A PMXStanModel object must be provided.")
  
# check inits for observations
model.specs <- model$get.model.specs()

# get parameters  
modeltype <- model.specs$type
m.route <- model.specs$route
m.obs.idx <- model.specs$obs.state
m.pk.solver <- model.specs$solver

# double check model specifications
if(modeltype == "PKPD" & m.pk.solver == "closed_form")
{
  m.pk.solver <- "ODE"   # reset for consistency for all PKPD models
  warning("The system will be solved by ODE solver, since the model is set as PKPD type.")
}

if(hasArg(data.source)) {
  d0 = data.source
} else { # data.source not provided but data.file provided
  ### check data format
  namelength <- nchar(data.file)
  datext <- substr(data.file, namelength-3, namelength)  
  if(datext %in% c(".txt",".csv")) {
    switch(datext,
           ".txt" = {
               d0 <- read.table(data.file, header = T)
           },
           ".csv" = {
               d0 <- read.csv(data.file, header = T)
           }
    )
  } else stop("Input data need to be in .txt or .csv format.")
}

### check minimally required columns
# ID(for population data),TIME,DV,AMT,EVID,RATE(for IV infusion only)
datcols <- colnames(d0)
if(!any(datcols == "ID"))
  if(data.type == "population"){ 
    #for population data, need ID
    stop("Population data need to include ID.")
  } else{ 
    # for individual data, add an ID column if missing
    d0$ID <- 1
  }
if(!any(datcols == "TIME")) stop("Input data need to include TIME.")
if(!any(datcols == "AMT")) stop("Input data need to include AMT.")
if(!any(datcols == "EVID")) stop("Input data need to include EVID.")
if(!any(datcols == "DV")) {
  if(any(datcols == "LIDV")){
    d0$DV <- d0$LIDV
  } else{ # no DV, no LIDV
    if(any(datcols == "Cp")) {
      d0$DV <- d0$Cp
    } else{ # no DV, no LIDV, no Cp
      stop("Input data need to include DV or LIDV or Cp.")
    }
  }
}
## for infusion, also need RATE
if(m.route == "IV_infusion")
  if(!any(datcols == "RATE")) stop("Input data need to include RATE for IV infusion dosing.")

### handle data
# sort each individual's data by time
d0 <- d0[order(d0$ID, d0$TIME),]

# check original data
pts.id <- unique(d0$ID)
# number of patients
nsub0 <- length(pts.id)

# check ndose and label patients with ndose=0
d0.ds <- subset(d0, EVID==1)
# number of doses for each patient
ndose0 <- t(sapply(1:nsub0, function(idx) c(pts.id[idx], sum(d0.ds$ID==pts.id[idx]))))
ndose0 <- data.frame(ID = ndose0[,1], sel.ds = (ndose0[,2]>0))

# check nobs and label patients with nobs=0
d0.dv <- subset(d0, EVID==0)
# number of Cp observations for each patient
nobs0 <- t(sapply(1:nsub0, function(idx) c(pts.id[idx], sum(d0.dv$ID==pts.id[idx]))))
nobs0 <- data.frame(ID = nobs0[,1], sel.obs = (nobs0[,2]>0))

# select patients with both nonzero numbers of doses and Cp
d1 <- merge(d0, ndose0, by="ID")
d1 <- merge(d1, nobs0, by="ID")
d2 <- subset(d1, sel.ds=="TRUE" & sel.obs=="TRUE")
d2 <- d2[,1:(ncol(d2)-2)]
d2 <- d2[order(d2$ID, d2$TIME, -d2$EVID),]  # dosing ahead of observation

# get patients
pts.id <- unique(d2$ID)
# number of patients
nsub <- length(pts.id)

# double check available data
if(nsub == 0) stop("No subjects have both valid dosing and valid observation information.")

# double check data type
if(nsub > 1 & data.type == "individual")
  stop("Data contain more than one individual, but data type was specified as individual.")
if(nsub == 1 & data.type == "population")
  stop("Data contain valid information from only one individual, but data type was specified as population.")

# get inits
if(!hasArg(inits)) {
  m.obs.init <- rep(0, nsub)
  message("Initial values of observations were not provided by user. They are set to zero(s) by default.")
} else {
  if (is.numeric(inits)) {
    if (length(inits)==nsub) {
      m.obs.init <- inits
    } else {
      cat("IDs of subjections with both valid dosing and valid observation information:\n")
      print(pts.id)
      stop("Number of eligible subjects (printed above) does not match length of user-specified initial observations.")
    }
  } else { 
    if (class(inits)=="character") {
      if(!any(datcols == inits)) {
        stop("User-specified column for initial observations cannot be found in data file.")
      } else {
        d.init <- d2[, c("ID",inits)]
        d.init <- d.init[!duplicated(d.init$ID),]
        m.obs.init <- d.init[, inits]
      }
    } else {
      stop("Input \"inits\" needs to be of either numeric or character class.")
    }
  }
}

# get covariates
ncov <- 0
d.cov <- NULL
if (!is.null(covar)) {
  d.cov <- data.frame(ID = pts.id)
  for (cvr in covar) {
    if (!any(datcols == cvr)) {
      warning(sprintf("User-specified covariate %s cannot be found in data file.", cvr))
    } else {
      d.cvr <- d2[, c("ID",cvr)]
      d.cvr <- d.cvr[!duplicated(d.cvr$ID),]
      d.cov <- cbind(d.cov, d.cvr[, 2])
      ncov <- ncov + 1
      colnames(d.cov)[ncov+1] <- cvr
    }
  }
  if (ncov==0) {
    d.cov <- NULL
    warning("None of the user-specified covariate(s) can be found in data file.")
  }
}

# subset for observations of dependent variable
d.dv <- subset(d2, EVID==0)
# number of observations for each patient
nobs <- sapply(1:nsub, function(idx) sum(d.dv$ID==pts.id[idx]))
# observations time for each Cp value
obs_time <- d.dv$TIME
# concentrations
dv <- d.dv$DV

# subset for dosing
d.ds <- subset(d2, EVID==1)
# dose amount
dose_amt <- d.ds$AMT
# dose time
dose_time <- d.ds$TIME

switch(m.pk.solver,
       "closed_form" = {
           # get concentration (same as DV)
           conc <- dv
           
           # get dosing information
	   ndose <- sapply(1:nsub, function(idx) sum(d.ds$ID==pts.id[idx]))
	   
	   if(m.route == "IV_infusion") {
	     inf_time <- d.ds$AMT/d.ds$RATE
	   }
	   
	   # handle individual data
	   if(nsub==1){
	     nsub <- 2
	     nobs <- c(nobs,1)
	     obs_time <- c(obs_time,0)
	     conc <- c(conc,0)
	     ndose <- c(ndose,1)
	     dose_amt <- c(dose_amt,0)
	     dose_time <- c(dose_time,0)
	     if(m.route == "IV_infusion") inf_time <- c(inf_time,inf_time[length(inf_time)])
	   }
	   
	   # output
	   if(m.route == "IV_infusion"){
	     standat <- list(NSUB=nsub, NOBS=nobs, obs_time=obs_time, conc=conc,
	          NDOSE=ndose, dose_amt=dose_amt, dose_time=dose_time,
	          inf_time=inf_time
	         )
	   } else{
	     standat <- list(NSUB=nsub, NOBS=nobs, obs_time=obs_time, conc=conc,
	          NDOSE=ndose, dose_amt=dose_amt, dose_time=dose_time
	         )
           }
           
       },
       "ODE" = {
           d3 <- d2
           
           # handle individual data
           # create a fake individual with 0 dose
           if(nsub==1){
             nsub <- 2
             nobs <- c(nobs,1)
             ID.fake <- pts.id + 1
             pts.id <- c(pts.id, ID.fake)
             # create a fake observation record
             d.dv.fake <- d.dv[1,]
             d.dv.fake$ID <- ID.fake
             dv <- c(dv, d.dv.fake$DV)
             # create a fake dosing record
             d.ds.fake <- d.ds[1,]
             d.ds.fake$ID <- ID.fake
             d.ds.fake$AMT <- 0
             # data for the fake individual
             d.fake <- rbind(d.dv.fake, d.ds.fake)
             # attach to the original data
             d3 <- rbind(d3, d.fake)
             d3 <- d3[order(d3$ID, d3$TIME, -d3$EVID),] 
             # handle inits
             m.obs.init <- c(m.obs.init, m.obs.init)
           }
           
           # handle infusion case
           if(m.route == "IV_infusion") 
           {   # need to add infusion-off events
	       inf_duration <- d.ds$AMT/d.ds$RATE
	       inf_off_time <- dose_time + inf_duration
	       inf_off_rate <- -d.ds$RATE
	       d.inf.off <- d.ds
	       d.inf.off$TIME <- inf_off_time
	       d.inf.off$RATE <- inf_off_rate
	       d.ds.exp <- rbind(d.ds, d.inf.off)
	       d.ds.exp <- d.ds.exp[order(d.ds.exp$ID, d.ds.exp$TIME),]
	       # combine with observation events
	       d3 <- rbind(d.dv, d.ds.exp)
	       d3 <- d3[order(d3$ID, d3$TIME, -d3$EVID),]
	   } 
           
           # nevts
	   nevts <- sapply(1:nsub, function(idx) sum(d3$ID==pts.id[idx]))
           # evt_time
           evt_time <- d3$TIME
           # evid
           evid <- d3$EVID
           # inits
           states <- scan("STATE_VARS.txt", character(0), quiet=T)
           n.state.var <- length(states)
           ## set to all 0's for PK model
           inits <- rep(0, n.state.var*nsub)
           ## need to plug in initial values of response state for PKPD model
           if(!is.null(m.obs.init))
           {
             if (length(m.obs.init)==nsub) {
               inits[m.obs.idx*seq(0,(nsub-1))+m.obs.idx] <- m.obs.init
             } else {  # number of initial values inconsistent with nsub
               stop("Mis-match found between initial values of number of subject. Please check pts.id and change m.obs.init accordingly.")
             }
           }
           
           # output
	   if(m.route == "IV_infusion"){
	     if(modeltype == "PK"){
	       standat <- list(NSUB=nsub, NOBS=nobs, NEVTS=nevts,
	                      conc=dv, evt_time=evt_time, evid=10101*evid,
                              dose_rate=d3$RATE, inits=inits
	                     )
	     } else{ # PKPD: may need to change evid according to ODE inputs
	       standat <- list(NSUB=nsub, NOBS=nobs, NEVTS=nevts,
	       	              y=dv, evt_time=evt_time, evid=10101*evid,
	                      dose_rate=d3$RATE, inits=inits
	                     )
	     }
	   } else{
	     if(modeltype == "PK"){
	       standat <- list(NSUB=nsub, NOBS=nobs, NEVTS=nevts,
	                      conc=dv, evt_time=evt_time, evid=101*evid,
                              dose_amt=d3$AMT, inits=inits
	                     )
	     } else{ # PKPD: may need to change evid according to ODE inputs
	       standat <- list(NSUB=nsub, NOBS=nobs, NEVTS=nevts,
	       	              y=dv, evt_time=evt_time, evid=101*evid,
	                      dose_amt=d3$AMT, inits=inits
	                     )
	       }
             }    
           }
)# end of switch

if (data.type == "individual") {
  list(data.type = data.type, standata = standat, ID = pts.id[1], d.cov = d.cov)
} else { # population
  list(data.type = data.type, standata = standat, ID = pts.id, d.cov = d.cov)
}

}
