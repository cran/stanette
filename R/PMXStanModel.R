#' @title Creation of a \code{PMXStanModel} object 
#'
#' @description 
#' Initializes an object of class \code{PMXStanModel} with methods for generating and compiling Stan code and 
#' querying model specifications. 
#'
#' @details
#' \code{PMXStanModel} serves as an interface for practical PK/PD modeling using Stan under a Bayesian framework. 
#' The first step of building such a model is for a user to provide model specifications (for more details on 
#' specification arguments and default values, please refer to the \emph{Arguments} section). With a proper set of 
#' specifications, a model-specific Stan source code is then generated based on a generic template code for the 
#' associated model type. Users then can choose to compile the auto-generated Stan source code directly 
#' to a self-contained platform-specific executable, or to modify the Stan source code according to their own modeling 
#' strategies before compiling it to an executable. 
#'
#' The template Stan code serves at least two purposes. First, it presents a grammar-corrected ready-to-be-used code 
#' that already takes care of most technical details necessary to build such a model in Stan; secondly, it makes sure 
#' that everything needed for running Stan sampling is internally consistent, from preparing the compatible data list, 
#' interpreting various dosing events and schedule, to calling the appropriate solver. In a template Stan code, model 
#' parameters to be estimated are defined as \code{theta}'s and inter-individual variabilities for parameters are defined 
#' as \code{eta}'s, to follow the conventions in pharmacometrics practice. All priors are set as non-informative. 
#' It needs to be kept in mind that the code is fully accessible and modifiable by a user; therefore, a user can add 
#' his/her own customized code by changing certain part of the auto-generated code conveniently. The modified Stan code  
#' can be re-compiled at any time as long as the changes made by users comply with standard Stan grammar. It is suggested  
#' that in these cases, a user begins with making small changes and testing the compilation (as well as compatibility with 
#' the input data if necessary) gradually before going to more significant modifications.
#' 
#' @param type a string to specify the type of model: "PK" (default) or "PKPD".
#' @param path a string to specify the path that will be used to store the model, Stan input data, and model fitting 
#'             and diagnostic results.
#' @param route a string to specify mode of drug administration: "1st_order_abs" (default), "IV_bolus", or "IV_infusion".
#' @param solver a string to specify which solver to be used for a PK model: "closed_form" (default) or "ODE";
#'               ignored for PKPD models in ODE form.
#' @param ode a string to specify equations of the ODE system for a PKPD model; ignored for PK models.
#' @param pk.struct a string to specify the PK model structure: "1-cmpt", "2-cmpt" (default), or "3-cmpt"; 
#'                  ignored for PKPD models in ODE form.
#' @param pk.param a string to specify the method for PK model parameterization: "CL_V" (default), or "micro_rate";
#'                 ignored for PKPD models in ODE form.
#' @param obs.state an integer to specify the index of the state variable corresponding to observed data for a PKPD model in ODE form;
#'                  ignored for PK model with closed-form solution, but required for PKPD models.
#' @param theta a string or vector of strings to specify which parameters will be estimated for a PKPD model in ODE form; 
#'              ignored for PK model with closed-form solution, but cannot be \code{NULL} for PKPD models.
#' @param eta a string or vector of strings to specify which parameters will have inter-individual variability for a PKPD model 
#'            in ODE form; ignored for PK model with closed-form solution, and can be \code{NULL} for PKPD models.
#' @param fixed a vector of strings to specify which parameter values will be fixed at constants for a PKPD model in ODE 
#'              form; ignored for PK model with closed-form solution, and can be \code{NULL} for PKPD models.
#' @param compile a logical variable indicating whether to compile the generated Stan code during the initialization 
#'                process (\code{TRUE}) or not (\code{FALSE}, as default). Note that a \code{PMXStanModel} object can be 
#'                compiled at any time after initialization by calling the \code{compile} method (see \emph{Value} 
#'                and \emph{Details} sections).
#'
#' @return 
#' A \code{PMXStanModel} object, with the following list of methods:
#' \item{get.model.specs}{returns model specification.}
#' \item{generate.stancode}{generates Stan code according to model specification.}
#' \item{get.state.var}{returns names of state variables in the ODE system. }
#' \item{get.ode.par}{returns names of parameters in the ODE system.}
#' \item{get.ntheta}{returns number of parameters to be estimated for an ODE-based PKPD model.}
#' \item{get.neta}{returns number of inter-individual variabilities set by user for an ODE-based PKPD model.}
#' \item{get.stan.file}{returns path of the auto-generated Stan file.}
#' \item{compile.stanmodel}{compiles the Stan file, either auto-generated (default) or with a path specified by user.   
#'                          This function has a generic form \code{\link{compile}}. 
#'                         }
#' \item{get.compile.status}{returns a logical indicator whether the model-associated Stan file has been compiled or not.}
#' \item{retrieve.stanmodel}{returns the compiled Stan model.}
#' \item{print.model}{returns the compiling status and the path of the model-associated Stan file, as well as model specifications.
#'                    This function has a generic form \code{\link{print}}.   
#'                   }
#'
#' @seealso 
#' \code{\link{prepareInputData}} for transformation of a NONMEM-readable dataset to a list compatiable to auto-generated 
#' model-specific Stan code;
#' \code{\link{PMXStanFit}} for how to run the compiled Stan executable for the model, link it with the input data and 
#' generate posterior samples of parameters, and perform model diagnostics;
#' \code{rstan-stanmodel} for a class specifically referred to the compiled model.
#'
#' @references 
#' The Stan Development Team. \emph{Stan Modeling Language User's Guide and Reference Manual}. \url{http://mc-stan.org}. 
#'
#' @author Yuan Xiong and Wenping Wang
#'
#' @examples
#' \dontrun{
#' ### A population PK model
#' m1 <- PMXStanModel(path = "pk_m1", compile = TRUE)
#' print(m1)
#'
#' ### A population PKPD model
#' ode <- "
#'   C2 = centr/V;
#'   d/dt(depot) =-ka*depot;
#'   d/dt(centr) = ka*depot - ke*centr;
#'   d/dt(eff) = (1+Emax*C2/(C2+EC50))*Kin - Kout*eff;
#' "
#' instant.stan.extension(ode)
#'
#' m5 <- PMXStanModel(type = "PKPD", 
#'                    path = "pkpd_m1",
#'                    ode = ode, 
#'                    theta= c("Emax","EC50"),
#'                    eta = c("Emax","EC50"),
#'                    fixed = c(V=1, ka=0.5, ke=0.4, Kin=0.5, Kout=0.5),
#'                    obs.state = 3
#'                   )
#' compile(m5)
#' }
#
# Documentation last updated: 6Apr2016
#
# Code last updated: 4Apr2016 by Yuan Xiong
# Created 28Oct2015 by Yuan Xiong

PMXStanModel <- function(type = "PK", 
                        path = "model_temp",
                        route = "1st_order_abs",
                        solver = NULL,
                        ode = NULL,
                        pk.struct = "2-cmpt",
                        pk.param = "CL_V",
                        obs.state = NULL,
                        theta = NULL,
                        eta = NULL,
                        fixed = NULL,
                        compile = FALSE
                       )
{ 
  ########################################################
  # Validate model specification inputs
  
  # check input arguments
  if(!hasArg(type)) message("Model type is set to \"PK\" by default.")
  if(!hasArg(path)) message("Model path is set to \"model_temp\" by default.")
  if(!hasArg(route)) message("Drug administration is set to \"1st_order_abs\" by default.")
  
  # initialization
  .model.specs <- NULL
  .state.var <- NULL
  .ode.par <- NULL
  .ntheta <- NULL
  .neta <- NULL
  
  if(type %in% c("PK","PKPD")){
  
    switch(type,
               "PK" = {
                 if(!hasArg(pk.struct)) 
                   message("PK model struct is set to \"2-cmpt\" by default.")
                 if(!(pk.struct %in% c("1-cmpt","2-cmpt","3-cmpt")))
		   stop("Unrecognized PK model structure. Please select from \"1-cmpt\", \"2-cmpt\", or \"3-cmpt\".")

                 if (!hasArg(solver)) {
                   solver <- "closed_form"
                   message("Solver for the PK model is set to \"closed_form\" by default.")
                 }
                 if(!(solver %in% c("closed_form","ODE")))
		   stop("Unrecognized model solver. Please select from \"closed_form\" or \"ODE\".")
                 if (solver=="closed_form" & !is.null(ode)) {
                   message("Input ODE(s) will be ignored since closed-form solution was chosen. \nIf user-specific ODE 
                   system is preferred, please change model type to \"PKPD\".")
                 }
                 if (solver=="ODE" & !is.null(ode)) {
                   message("ODE(s) for PK models will be generated automatically from model specification. \n
                   If user-specific ODE system is preferred, please change model type to \"PKPD\".")
                 }
                 
                 if(!hasArg(pk.param)) 
		   message("PK model parameterization is set to \"CL_V\" by default.")
                 if(!(pk.param %in% c("CL_V","micro_rate")))
		   stop("Unrecognized PK model parameterization. Please select from \"CL_V\" or \"micro_const\".")
                 
               },
               "PKPD" = {
                 if (is.null(ode) ) {
                   stop("Please provide ODE system equation(s).")
                 } else {
                   .state.var <- scan("STATE_VARS.txt", character(0), quiet=T)
                   .ode.par <- scan("ODE_PARS.txt", character(0), quiet=T)

                 }
                 .npars <- length(.ode.par)
                 
                 # NOTE: theta cannot be NULL for PKPD models
                 if (is.null(theta)){
                   stop("Please provide system parameter(s) to be estimated.")
                 }
                 if (FALSE%in%(theta%in%.ode.par)) {
                   stop("Mis-specified parameter(s) in theta found.")
                 }
                 idx.theta <- sapply(theta, function(idvtheta) match(idvtheta, .ode.par))
                 .ntheta <- length(theta) 
                 
                 # NOTE: eta can be NULL for PKPD models
                 if (!is.null(eta)){ 
                   if (FALSE%in%(eta%in%theta)) stop("Mis-specified parameter(s) in eta found.")
                 }
                 flag.eta <- sapply(theta, function(idvtheta) as.numeric(idvtheta%in%eta))
                 count.eta <- cumsum(flag.eta)
                 .neta <- length(eta)
                 
                 # NOTE: argument "fixed" can be NULL for PKPD models
                 if(hasArg(fixed)) { # constant parameters provided
                   const.par.name <- names(fixed)
		   if(FALSE%in%(const.par.name%in%.ode.par)) { # parameter name not found
		     stop("Mis-specified parameter(s) in argument 'fixed' found.")
		   } else { # parameter names found
		     ### get index for constant parameters
		     idx.const <- sapply(const.par.name, function(par.const) match(par.const, .ode.par))
		     idx.all <- sort(c(idx.theta, idx.const))
		     if (length(idx.all)!=.npars) { # length doesn't match
		       stop("All parameters need to be specified either in arguments 'theta' or 'fixed', but not in both.")
		     } else { # length matches
		       if (FALSE%in%(idx.all==seq(.npars))) stop("Mis-specified parameters in arguments 'theta' and/or 'fixed' found.")
		     }
		   }
                 } else { # constant parameters not provided
                   if (.ntheta<.npars) stop("All parameters need to be specified either in arguments 'theta' or 'fixed'.")
                 }
                 
                 if (is.null(obs.state)){
                   n.state.var <- length(.state.var)
                   obs.state <- n.state.var
                   message("Index of observed state variable is set to be the index of the last state variable by 
                   default.")
                 }
                 if (is.null(solver)) solver <- "ODE"
                 pk.struct <- NULL
                 pk.param <- NULL
               }
          )     
    } else stop("Please provide a proper model type; needs to be \"PK\" or \"PKPD\".")
  
  if(!(route %in% c("1st_order_abs","IV_bolus","IV_infusion")))
    stop("Unrecognized drug administration. Please select from \"1st_order_abs\", \"IV_bolus\", or \"IV_infusion\".")
  
  # validated specifications
  .model.specs <- list(type = type, 
  		      path = path,
  		      route = route,
  		      solver = solver,
  		      ode = ode,
  		      pk.struct = pk.struct,
  		      pk.param = pk.param,
  		      obs.state = obs.state,
  		      theta = theta,
  		      eta = eta,
  		      fixed = fixed
                     )
   
  #######################################################
  # Generate Stan code  
  .stanfilename <- NULL
  .stanfilepath <- NULL
  
  generate.stancode <- function() 
  {
    # Modified 24Nov2015 by Yuan Xiong
    #   - return numbers of theta's and eta's for the convenience of post-processing
    # Modified 11Nov2015 from v7 by Yuan Xiong 
    #   - a unified wrapper function for both PK and PKPD model code generation
    
    if (.model.specs$type == "PK") {
      stanfilereturn <- generateStanCodePK(.model.specs)
    } else {
      stanfilereturn <- generateStanCodePKPD(.model.specs)
    }
    .stanfilename <<- stanfilereturn$stanfilename
    .ntheta <<- stanfilereturn$ntheta
    .neta <<- stanfilereturn$neta
    
    .stanfilepath <<- file.path(.model.specs$path, .stanfilename)
      
  }
  ### force to run
  generate.stancode()

  #####################################################
  # Compile Stan model
  .stanmodel <- NULL
  .compile.status <- FALSE
  
  compile.stanmodel <- function(stanfilepath=NULL) 
  {
    if (is.null(stanfilepath)) stanfilepath <- .stanfilepath
    .stanmodel <<- stan_model(stanfilepath)
    .compile.status <<- TRUE
  }
  ### force to run or not dependent on input
  if(compile) compile.stanmodel()
  
  #####################################################
  # Retrieve a compiled Stan model
  retrieve.stanmodel <- function() 
  {
    if (is.null(.stanmodel)) message("The model has not been compiled.")
    .stanmodel 
  }
  
  ######################################################
  # Print key attributes of the model at the screen
  print.model <- function() 
  {
    cat(
          sprintf(" Compiling status: %s \n\n", .compile.status),
          sprintf("Path of Stan file: %s \n\n", .stanfilepath)
    )
    
    cat(" Model specifications: \n")
    print(.model.specs)      
  }
  
  #####################################################
  # Wrap up
  out <- list(get.model.specs = function() .model.specs,      
             generate.stancode = generate.stancode,
             get.state.var = function() .state.var,
             get.ode.par = function() .ode.par,
             get.ntheta = function() .ntheta,
             get.neta = function() .neta,
             get.stan.file = function() .stanfilepath,
             compile.stanmodel = compile.stanmodel,
             get.compile.status = function() .compile.status,
             retrieve.stanmodel = retrieve.stanmodel,
             print.model = print.model
            )            
  class(out) <- "PMXStanModel"
  out
  
} # end of PMXStanModel

#########################################################
### generic functions

#' @title Compilation of Stan code
#'
#' @description
#' Compiles a piece of model-specific Stan code, as part of a \code{PMXStanModel} object, into an executable.
#'
#' @details
#' This is a generic version of the method \code{compile.stanmodel()} for the \code{\link{PMXStanModel}} class. 
#' The compilation step can also be performed simultaneous during the initialization process of a 
#' \code{PMXStanModel} object, by setting the argument \code{compile = TRUE}.
#'
#' @param model a \code{PMXStanModel} object.
#' @param stanfilepath a string for user to specify the path of the Stan code. The default is NULL, and the
#'                      location of the Stan file to be compiled will be decided by \code{model} above. 
#'
#' @return No explicit return; a successful compilation will generate an executable code as part of the current
#'         \code{PMXStanModel} object.
#' 
#' @seealso
#' \code{\link{PMXStanModel}} for the method \code{compile.stanmodel()} and the argument \code{compile}.
#'
#' @examples
#' \dontrun{
#' m <- PMXStanModel(path = "pk_m1")
#' compile(m)
#' 
#' m$compile.stanmodel()
#'
#' m <- PMXStanModel(path = "pk_m1", compile = T)
#' }

compile <- function(model, stanfilepath = NULL)
{
    UseMethod("compile", model)
}

compile.PMXStanModel <- function(model, stanfilepath = NULL)
{   
    model$compile.stanmodel(stanfilepath)
    
    if(!is.null(stanfilepath)){
      if(stanfilepath!=model$get.stan.file()) {
        warning("Path of compiled Stan file is different from the path of auto-generated Stan file.")
      }
    }
}

# print model
print.PMXStanModel <- function(x, ...)
{
  x$print.model()
}

#' @title Replication of an existing \code{PMXStanModel} object
#'
#' @description
#' Makes a copy of an existing \code{PMXStanModel} object into a new path.
#'
#' @details
#' The \code{copy()} function for a \code{PMXStanModel} object aims to provide a convenient way for the
#' following scenarios:
#'
#' - Fitting the same model to different datasets. With each copy of the model, a new input dataset  
#'   can be used for generating samples, and each will result in a \code{PMXStanFit} objects. These results
#'   can then be compared and summarized.
#' 
#' - An exploratory, stage-wise modeling exercise, commonly used when current knowledge is not sufficient
#'   to choose among options through the model building process. One can start from a copy of the same
#'   "parent" version of the model, make specific changes one at a time, and later make a choice based on
#'   certain criteria to compare different fittings.
#'
#' @param model a \code{PMXStanModel} object.
#' @param newpath a string for user to specify the path to store the new model object. 
#'                If not provided, the new path will be named by adding a suffix "_copy" to the path 
#'                of the model being copied.
#' @param compile a logical variable indicating whether to compile the Stan code in the new \code{PMXStanModel}  
#'                object during the initialization process (\code{TRUE}) or not (\code{FALSE}, as default). 
#'                
#' @return A replicated version of the specified existing \code{PMXStanModel} object.
#' 
#' @seealso
#' \code{\link{PMXStanModel}} for the creation of a new \code{PMXStanModel} object.
#'
#' @author Yuan Xiong and Wenping Wang
#'
#' @examples
#' \dontrun{
#' m1 <- PMXStanModel(path = "pk_m1")
#' m2 <- copy(m1, newpath = "pk_m2", compile = TRUE)
#' }
copy <- function(model, newpath = NULL, compile = FALSE)
{
    UseMethod("copy", model)
}

copy.PMXStanModel <- function(model, newpath = NULL, compile = FALSE)
{
  # get model specification
  model.specs <- model$get.model.specs()
  
  # set model path if not provided
  if (!hasArg(newpath)) {
    newpath <- paste(model.specs$path, "_copy", sep="")
  }
  
  # generate a new PMXStanModel object
  model.new <- PMXStanModel(type = model.specs$type, 
  			   path = newpath,
  			   route = model.specs$route,
  			   solver = model.specs$solver,
  			   ode = model.specs$ode,
  			   pk.struct = model.specs$pk.struct,
  			   pk.param = model.specs$pk.param,
  			   obs.state = model.specs$obs.state,
  			   theta = model.specs$theta,
  			   eta = model.specs$eta,
  			   fixed = model.specs$fixed,
  			   compile = compile
  			  )
  model.new
}
