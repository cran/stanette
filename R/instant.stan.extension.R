#' generate a customized NUTS-compatible ODE solver
#'
#' generate a customized NUTS-compatible ODE solver
#' 
#' @param ode_str ODE string
#' @param covar a vector of character with names of covariates
#' @param multi_state a logical if multi-state 
#' @return NULL
#' @details
#'    generate a customized NUTS-compatible ODE solver
#' 
#' @author Wenping Wang
#' @examples
#' \dontrun{
#' ode <- "
#' d/dt(depot) =-KA*depot;
#' d/dt(centr) = KA*depot - KE*centr;
#' "
#' instant.stan.extension(ode)
#' }
instant.stan.extension = function(ode_str=NULL, covar=NULL, multi_state=FALSE)
{
    if (is.null(ode_str)) {
      stop("please provide ODE string")
    }
    
    if (!is.loaded("parse_ode")) {
      dyn.load(system.file("libs/odeparser.so", package = "stanette"))
    }
    
    .dir = sprintf("%s/",tempdir())
    fname = function(s) {
		sprintf("%s/%s", .dir, s)
	}
	.mod = fname("model.txt")
    
    if (is.null(covar)) {
		opts = c(
		  "include/generic_ode_interface_template.txt",
		  "include/generic_ode_interface_template_multi_state.txt"
		)
		.tmpl = system.file(opts[multi_state+1], package = "stanette")
		cat(ode_str, file=.mod)
		nvar = 0
	} else {
		.tmpl = system.file("include/generic_ode_interface_template_cov.txt", package = "stanette")
		covar = strsplit(covar, "[,| \t]+")[[1]]
		nvar = length(covar)
		if (prod(nchar(covar))*nvar == 0)
			stop("unrecoganized covar string")
		paste(covar, "= 9999.999 + -9999.999;")
		cat(paste(covar, "= 9999.999 + -9999.999;"), file=.mod, sep="\n")
		cat(ode_str, file="model.txt", append=TRUE)
	}
	.extn = system.file("include/stan/math/PMXStan/", package = "stanette")
	x = .C("parse_ode", .tmpl, .mod, sprintf("%s/generic_ode_interface.hpp", .extn), as.character(nvar), .dir)

    pars = scan(fname("ODE_PARS.txt"), what="", quiet=TRUE)
    cat("A new ODE extension for Stan has been created.\n")
    cat(sprintf("System parameters are: %s\n", paste(pars, collapse = " ")))
    invisible()
}
