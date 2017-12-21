####
#
# Penalized Ensemble Predictor ("PEP") S4 -class/object 
# Package 'ePCR'
# Teemu Daniel Laajala, teelaa@utu.fi
# Description: Collects a set of 'psp's to an ensemble that can be used for novel predictions
#
####

# Make sure few formal classes are defined to avoid warnings at new class definition
setOldClass("glmnet")
setOldClass("coxnet")
setOldClass("Surv")
setOldClass("impute")

#' Penalized Ensemble Predictor (PEP) S4-class ensemble consisting of individual PSP-members
#'
#' This class constructs an ensemble of individual Penalized Ensemble Predictor (PSP-class) members. Each member contributes to the model output equally, and ensemble-level functions wrap up individual predictions into an averaged ensemble prediction. The user may define an arbitrary number of PSPs and tailor them to suit the particular needs, and then provide them as a list to the PEP-constructor. As such, constructing well tailored individual ensemble members (of PSP-class) in order to produce a powerful ensemble (of PEP-class) is important on both levels.
#'
#' @slot PSPs List of PSP-objects that will be treated as equal members of the ensemble
#' @slot description A character string describing the structure or purpose of the ensemble
#' @slot features A character list of variable/feature names
#' @slot dictionary A named list of above variables/features and their more precise description
#' @slot predens A function for compiling all predictions from the PSPs into consensus prediction
#' @slot prednorm A function for normalizing the predictions e.g. to risk scores in [0,1]
#' @name PEP-class
#' @rdname PEP-class
#' @examples
#' \dontrun{
#' # The PEP-construction is wrapped in NOT RUN, because cross-validating multiple PSPs
#' # is very time consuming especially if a tight grid of alpha/lambda is to be explored.
#' # The simulated data from Turku University Hospital (TYKS) is used as an example:
#' data(TYKSSIMU)
#' 
#' # Two cohorts and corresponding data matrices:
#' head(xMEDISIMU)
#' head(xTEXTSIMU)
#' # Two survival responses:
#' head(yMEDISIMU)
#' head(xTEXTSIMU)
#'
#' # Search L1/L2 norm alpha-grid with 10 values between [0,1]
#' aseq <- seq(from=0, to=1, by=0.1)
#' # Lambda sequence penalization is of 100 length conditional for each alpha
#' nlamb <- 100
#' 
#' library(survival)
#' # Create three ensemble members; one for MEDI cohort, one for TEXT cohort,
#' # and finally one member that combines both cohorts simultaneously in a coxnet
#' psp1 <- new("PSP", x = rbind(xMEDISIMU, xTEXTSIMU), 
#'	y = Surv(rbind(yMEDISIMU, yTEXTSIMU)[,"surv"]),
#'	plot = TRUE, alphaseq = aseq, scorefunc = score.cindex, seed = 1,
#'	folds = 10, nlambda = nlamb)
#' psp2 <- new("PSP", x = xMEDISIMU, 
#'	y = Surv(yMEDISIMU[,"surv"]),
#'	plot = TRUE, alphaseq = aseq, scorefunc = score.cindex, seed = 1,
#'	folds = 10, nlambda = nlamb)
#' psp3 <- new("PSP", x = xTEXTSIMU, 
#' 	y = Surv(yTEXTSIMU[,"surv"]),
#'	plot = TRUE, alphaseq = aseq, scorefunc = score.cindex, seed = 1,
#'	folds = 10, nlambda = nlamb)
#' par(mfrow=c(1,3))
#' plot(psp1); plot(psp2); plot(psp3); # Inspect the alpha/lambda surfaces
#'
#' # Create an ensemble of the above 3 members
#' simuens <- new("PEP", PSPs = list(psp1, psp2, psp3))
#' simuens
#' # Ready PEP-object can be used for novel predictions etc
#'
#' }
#' 
#' # Run example predictions from a previously optimized PEP-model
#' data(ePCRmodels)
#' data(TYKSSIMU)
#' 
#' # Perform risk predictions from the joint cohort ensemble member as an example
#' MEDIpred <- predict(TYKS@PSPs[[1]]@fit, s=TYKS@PSPs[[1]]@optimum["Lambda"], 
#'	newx = conforminput(TYKS@PSPs[[1]], xMEDISIMU))[,1]
#' TEXTpred <- predict(TYKS@PSPs[[1]]@fit, s=TYKS@PSPs[[1]]@optimum["Lambda"], 
#'	newx = conforminput(TYKS@PSPs[[1]], xTEXTSIMU))[,1]
#'
#' # Risk scores obtained for the new patients (arbitrary unit as per Cox regression)
#' head(MEDIpred)
#' head(TEXTpred)
#'
#' @exportClass PEP
setClass(
	Class="PEP",
	# variables of the object; slots:
	representation = representation(
		PSPs = "list", # A list of PSP objects
		description = "character", # Preferably a rather lengthy (pre-formatted) text of what is the purpose/aim of this particular PEP
		features = "character", # A character vector of the features that are required to perform ensemble predictions with this PEP
		dictionary = "list", # A functional dictionary for mapping alternate feature names to 'features' vector
		predens = "function", # A function that compiles predictions from all PSPs and compiles them into a single risk score
		prednorm = "function" # A function that normalizes the predictions
	),
	prototype = prototype( # Empty prototype with no PSPs but proposals for ensemble construction and normalization
		PSPs = list(), 
		description = "An extended description here of what the PEP contains and preferably PSP-specific explanations.",
		features = "",
		predens = meanrank,
		prednorm = normriskrank
	)
)

###
#	Initializer - Will collect invidual PSPs to a PEP
###
setMethod("initialize", "PEP",
	function(.Object,
		PSPs
){
	if(!class(PSPs)=="list"){
		stop("PSPs should be a list that consists of PSP-objects")	
	}else if(!all(unlist(lapply(PSPs, FUN=function(z) class(z)=="PSP")))){
		stop("The list PSPs should consist solely of PSP-objects")
	}
	.Object@PSPs = PSPs
	.Object@features = unique(unlist(lapply(.Object@PSPs, FUN=function(z) z@features)))
	
	# Construct 'features' vector based on what all variables PSPs require (unique columns in data matrix)
	
	return(.Object)
})

###
#	Functions of the PEP class (ensemble of PSPs)
###

#' PEP-methods
#' @name PEP-methods
NULL

# Show object that "cat"s object to R terminal if asked
setMethod("show", "PEP",
	function(object){
		cat("Penalized Ensemble Predictor\n")
		cat("Count of PSPs: ", length(object@PSPs))
		cat("\n")
	}
)
#' Print a PEP object
#' @rdname PEP-methods
#' @param x Generic x
#' @param ... Additional custom parameters passed on
#' @export
setMethod("print", "PEP",
	function(x, ...){
		cat("Penalized Ensemble Predictor\n")
		cat("Count of PSPs: ", length(x@PSPs))
		cat("\n")
	}
)

## TODO
# By default the mean CV surface in terms of alpha/lambda is plotted using hamlet-package's hmap-function
# @rdname PEP-methods
# @param y Generic y
# @param ... Additional custom parameters passed on
# @export plot
##setMethod("plot", "PEP",
##	function(x, y, ...){
##		stop("No plot function available yet for the whole PEP")
##	}
##)

#' predict.PEP: Predict for a novel patient from current PEP-ensemble
#' @rdname PEP-methods
#' @param object PEP-ensemble model object
#' @param type Type of prediction; either "response" or "ensemble"
#' @param newx New data matrix
#' @export
setMethod("predict", "PEP",
	function(object, type="response", newx){
		if(!missing(newx)){
			preds <- lapply(object@PSPs, FUN=function(z){
				glmnet::predict.coxnet(z@fit, newx=as.matrix(z@x.expand(newx)), s=z@optimum["Lambda"])
			})		
		}else{
			preds <- lapply(object@PSPs, FUN=function(z){
				glmnet::predict.coxnet(z@fit, newx=as.matrix(z@x.expand(z@x)), s=z@optimum["Lambda"])
			})		
		}
		if(type=="ensemble"){
			if(missing(newx)) stop("For a novel ensemble prediction a 'newx' should be provided")
			preds <- object@prednorm(object@predens(preds))
			names(preds) <- rownames(newx)
			preds
		}else{
			preds
		}
	}
	# Ensemble prediction
)

###
#	Other interesting utility functions for PEPs, such as KM
###

## TODO
# Kaplan-Meier with division at a given cutoff point within [0,1]
#
# @rdname PEP-methods
# @param ... Additional custom parameters passed on
# @param cutoff Cutoff point for division
# @export
###setGeneric("PEP.KM", function(object, cutoff=0.5) { standardGeneric("PEP.KM") })
# @rdname PSP-methods
# @exportMethod PEP.KM
###setMethod("PEP.KM", "PEP",
###	function(object, cutoff=0.5){
###		print(paste(cutoff, "\n", object@PSPl))
###	}
###)
