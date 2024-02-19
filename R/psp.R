####
#
# Penalized SINGLE Predictor ("PSP") S4 -class/object 
# Package 'ePCR'
# Teemu Daniel Laajala, teelaa@utu.fi
# Fits and saves the relevant parameters for an alpha/lambda grid optimization for a single penalized Cox regression model
#
####

# Make sure few formal classes are defined to avoid warnings at new class definition
setOldClass("glmnet")
setOldClass("coxnet")
setOldClass("Surv")
setOldClass("impute")

#' Penalized Single Predictor (PSP) S4-class as a member of PEP-ensembles
#'
#' PSP is a single penalized Cox regression model, where an alpha/lambda grid has been optimized using cross-validation and a chosen prediction metric. PSPs are single entities that will compile together into PEPs, the ensemble objects that will average over multiple PSPs to generate an ensemble prediction. Typically a single PSP models a part of the data, such as a cohort strata.
#'
#' @name PSP-class
#' @rdname PSP-class
#' @import methods
#' @importFrom impute impute.knn
#' @slot description A general user-provided string describing the PSP
#' @slot features A character vector indicating feature names
#' @slot strata Information whether data matrix x included substrata (will be used in plotting functions etc)
#' @slot alphaseq The sequence of alpha values to test, ranging between [0,1]; alpha = 0 being ridge regression, 0 < alpha < 1 being elastic net and alpha = 1 being LASSO
#' @slot cvfolds The number of cross-validation folds to utilize; by default 10
#' @slot nlambda The amount of lambda values utilized in each regularization path; by default 100 as in glmnet-package
#' @slot cvmean A matrix indicating the mean CV performance in alpha/lambda grid (preferred over median)
#' @slot cvmedian A matrix indicating the median CV performance in alpha/lambda grid
#' @slot cvstdev A matrix indicating the standard deviation in CV performance over the folds in the alpha/lambda grid
#' @slot cvmin A matrix indicating minimum CV performance in alpha/lambda grid
#' @slot cvmax A matrix indicating maximum CV performance in alpha/lambda grid
#' @slot score The scoring function, user-defined or one provided by ePCR package such as score.cindex or score.iAUC
#' @slot cvrepeat Number of cross-validation procedures to run multiple times and then average over, in order to reduce the effect of binning samples
#' @slot impute The imputation function used if provided matrix 'x' includes missing values; by default the impute.knn-function from BioConductor package 'impute'
#' @slot optimum The optimum in alpha/lambda grid, with optimal alpha and similarly for lambda
#' @slot seed The initial random seed used for cross-validation
#' @slot x The input data matrix
#' @slot x.expand A function that allows expansion of matrix 'x' to include interactions between variables; if no such are desired, this should be an identity function
#' @slot y The Surv-object as in survival-package, which serves as the response y
#' @slot fit The glmnet coxnet-object obtained with optimal alpha
#' @slot criterion The optimizing criterion; by default "min" for minimizing CV-error
#' @slot dictionary A list of discriptions for each variable
#' @slot regAUC A numeric vector for the AUC under regularization curve as computed by integrateRegCurve-function
#' @examples
#' # As an example, illustrate a naive PSP built on the small medication cohort
#' data(TYKSSIMU)
#' library(survival)
#' # Minimal example with much fewer patients and variables
#' psp_ex <- new("PSP", alphaseq=c(0.2, 0.8), nlambda=20, folds=3,
#' 	x = xMEDISIMU[1:80,c(1:20,40:50)], y = yMEDISIMU[1:80,"surv"],
#'	seeds = 1, score=score.cindex)
#'
#' plot(psp_ex) # Optimization surface of alpha/lambda
#' 
#' # Illustrate the use of some PSP-methods:
#' PSP.KM(psp_ex, cutoff = 0.5) # Kaplan-Meier
#' PSP.PCA(psp_ex) # PCA plot of training data
#' PSP.BOX(psp_ex) # Boxplots, here for the first training variable
#' PSP.CSP(psp_ex) # Cumulative survival probabilities for the training data
#' invisible(PSP.NA(psp_ex)) # Time-to-event Nelson-Aalen heuristic algorithm
#'
#' \dontrun{
#' # Computationally intensive novel PSP-fitting is omitted from the test runs
#' # Functions for readily fitted PSP-objects are illustrated above
#' data(TYKSSIMU)
#' library(survival)
#' psp_meditext <- new("PSP", x = rbind(xMEDISIMU, xTEXTSIMU), 
#'	y = Surv(rbind(yMEDISIMU, yTEXTSIMU)[,"surv"]),
#'	plot = TRUE, alphaseq = seq(0, 1, by=.01), scorefunc = score.cindex, 
#'	seed = 1, folds = 10, nlambda = 100)
#' plot(psp_meditext)
#' }
#' @exportClass PSP
setClass(
	Class="PSP",
	# variables of the object; slots:
	representation = representation(
		description = "character", # Preferably a lengthy description of what the PSP is (e.g. what batch)
		features = "character", # A character vector of the features that are required to perform ensemble predictions with this PSP
		strata = "factor", # If data was somehow stratified / in batches, each unique character instance is a batch of observations
		alphaseq = "numeric", # Each L1/L2 norm alpha \in [0,1] tested
		cvfolds = "numeric", # Number of folds in the CV
		nlambda = "numeric", # each row, f(alpha), has a different lambda sequence
		cvmean = "matrix", # Numeric matrix of CV-score means
		cvmedian = "matrix", # Numeric matrix of CV-score medians
		cvstdev = "matrix", # Numeric matrix of CV-score standard deviations
		cvmin = "matrix", # Numeric matrix of CV-score minimums
		cvmax = "matrix", # Numeric matrix of CV-score maximums
		score = "function", # Function fo score results in CV, e.g. C-index, iAUC etc
		cvrepeat = "numeric", # If CV was run multiple times and averaged over
		impute = "function", # Imputation function to use in case data is missing
		optimum = "numeric", # The location in cvmat[i,j] where the optimum is found
		seed = "numeric", # Random seed(s) used to initialize cross-validations
		x = "matrix", # Training input matrix X
		x.expand = "function", # Function that is used to expand x into xFull; since many new variables are extracted, it's not feasible to save the whole x in each PSP
		y = "Surv", # Training response vector/matrix Y as a survival object
		fit = "coxnet", # A fitted glmnet object for the optimal parameters (or lambda is a sequence but contains the optimum)
		criterion = "character", # What optimum criterion to use; by standard minimum mean CV error (iAUC or c-index)
		dictionary = "list", # A list of synonyms and explanations for variable names that may help map new predictions to the current fitted model
		regAUC = "numeric" # Values of running integrateRegCurves-function for the regularized model
	),
	# Proposed default values for some fields
	prototype = prototype(
		alphaseq = seq(from=0, to=1, by=.1),
		nlambda = 100,
		cvfolds = 10,
		score = score.cindex,
		criterion = "min",
		dictionary = list(),
		impute = function(x) { impute::impute.knn(as.matrix(x))$data } , # Notice: impute.knn is from Bioconductor, not CRAN
		x.expand = function(x) {
			cbind(
				x, # All original variables
				interact.all(x[,colnames(x)[which(apply(x, MARGIN=2, FUN=function(z) length(table(z))>10 | !all(z == round(z,0))))]]), # All interactions between numeric variables
				interact.part(input=x, first=colnames(x[,colnames(x)[which(apply(x, MARGIN=2, FUN=function(z) length(unique(z))>6 | !all(z == round(z,0))))]]), second=colnames(x[,colnames(x)[which(apply(x, MARGIN=2, FUN=function(z)  length(unique(z))<=6 & all(z == round(z,0))))]])) # All interactions where one is a numeric and the other is a binary variable
			)
		}
	)
)

###
#	Functions of the PSP class
###

###
#	Initializer! Will also run the CV, most likely heavy to call
###
setMethod("initialize", "PSP",
	function(.Object, 
		alphaseq = .Object@alphaseq, # Sequence of alphas
		nlambda = .Object@nlambda,  # How many lambda-penalization values will be tested per alpha
		folds = .Object@cvfolds, # Number of cross-validation folds
		x, y, # Data matrix X and the (survival) response Y
		seeds, # Should a fixed seed be utilized for reproducibility
		cvrepeat = 1, # Should cross-validation be repeated multiple times
		scorefunc = score.iAUC, # Scoring function -- changed default to iAUC
		plot = FALSE, # Should the CV Error curves be plotted per each alpha
		verb = 0, # Level of verbosity
		criterion = "min",
		strata = as.factor(rep(1, times=nrow(x))), # Which criterion to use to select optimal parameters; 1 = maximal mean CV score (c-index or iAUC), 2 = ...
		dictionary, # A named list object where within each unique feature name can be two objects: 'description' and 'synonyms'
		x.expand
		# , ... # removing ... during bugfixing to see if any parameters leak through
	){
		if(verb>-1) cat("--- Initializing new PSP object ---\n\n")
		# Run parameters saved
		.Object@alphaseq = alphaseq
		.Object@nlambda = nlambda
		.Object@cvfolds = folds
		.Object@score = scorefunc
		.Object@strata = strata
		# Custom interaction / data matrix x expansion function
		if(!missing(x.expand)) .Object@x.expand = x.expand
		if(!sum(is.na(x))==0){
			if(verb>-1) cat("--- Missing entries detected in x, running defined imputation function ---\n\n")
			x <- .Object@impute(x)
		}
		# Save data matrix and response vector
		.Object@x = as.matrix(x)
		.Object@y = y

		tmps <- list()
		# Inspecting the cv repeats & provided RNG seeds
		# If user has not provided cvrepeat parameter, it is assumed to be the sae as the number of provided seeds
		if(!missing(seeds)){
			if(!cvrepeat==length(seeds)) cvrepeat <- length(seeds)
		}
		for(z in 1:cvrepeat){	
			if(verb>-1) cat("--- Cross-validation (", folds, "-folds) repeat run ", z, " of ", cvrepeat, " ---\n\n")
			# Run the cross-validation in the grid
			# Give multiple seeds
			if(!missing(seeds)){
				tmps[[length(tmps)+1]] <- cv.grid(alphaseq = alphaseq, x = .Object@x.expand(x), y = y, folds = folds, nlamb=nlambda, scorefunc=scorefunc, plot=plot, verb=verb, seed=seeds[z])
			# Using random seeds and repeating cv according to desired cvrepeats
			}else{
				tmps[[length(tmps)+1]]  <- cv.grid(alphaseq = alphaseq, x = .Object@x.expand(x), y = y, folds = folds, nlamb=nlambda, scorefunc=scorefunc, plot=plot, verb=verb)
			}
		}
		.Object@cvrepeat = cvrepeat
		.Object@seed = seeds


		# If multiple CV runs were done, the key CV statistics are averaged over the binning repeats
		.Object@cvmean = Reduce("+", lapply(tmps, FUN=function(z) z[["mean"]]))/length(tmps)
		.Object@cvmedian = Reduce("+", lapply(tmps, FUN=function(z) z[["median"]]))/length(tmps)
		.Object@cvstdev = Reduce("+", lapply(tmps, FUN=function(z) z[["stdev"]]))/length(tmps)
		.Object@cvmin = Reduce("+", lapply(tmps, FUN=function(z) z[["min"]]))/length(tmps)
		.Object@cvmax = Reduce("+", lapply(tmps, FUN=function(z) z[["max"]]))/length(tmps)

		# Optimal fit according to the chosen criterion
		if(criterion %in% c("min", "min.mean", "max", "max.mean")){
			opt <- which(.Object@cvmean == max(.Object@cvmean, na.rm=T), arr.ind=TRUE)
			if(verb>0){
				cat("Content of detected 'opt' (optimum in alpha/lambda grid):\n")
				print(opt)
			}
			# Single optimum detected
			if(dim(opt)[1]==1){
				# Single unique optimum according to the criterion
				alphaopt <- as.numeric(rownames(opt)[1])
				lambdaopt <- opt[1,2] # The only optimum, second column is the lambda value
			# Multiple optima
			}else if(nrow(opt)==0){
				warning("Unable to detect optima in the CV grid; please inspect the CV results manually.\nThis could result for example due to a too high CV fold-count in relation to the N.")
				alphaopt <- NA
				lambdaopt <- NA
			}else{
				warning("Multiple optima detected in the CV grid; choosing the one with the highest alpha and then highest lambda, but manual inspection is highly encouraged.")
				# Favor high alpha, and then high lambda
				opt <- opt[order(opt[,1], rev(opt[,2])),] # Notice rev-addition to lambda (2nd column) - lambda values are in descending order
				opt <- opt[nrow(opt),,drop=F]			
				# Pick these this alpha/lambda optimum
				alphaopt <- as.numeric(rownames(opt)[1])
				lambdaopt <- opt[1,2] # The only optimum, second column is the lambda value
			}
		}else if(criterion %in% c("lambda.1se", "l.1se", "l1se")){ # The normal 1se definition from 'glmnet'
			opt <- which(.Object@cvmean == max(.Object@cvmean), arr.ind=TRUE)
			alphaopt <- as.numeric(rownames(opt)[1])
			lambdaopt <- opt[1,2] # The only optimum, second column is the lambda value
			# Start walking so as long as cvmean is within 1se of the minimal cvmean error (i.e. maximum iAUC or c-index)
			# This will provide a more penalized and conservative model and is generally advocated by 'glmnet' as cv minimum error will still overfit
			# Standard error
			se <- .Object@cvstdev[which(alphaseq==alphaopt), lambdaopt]/sqrt(folds)
			# Step from current minimum until we reach intercept model or min-|1se| < min-error criterion
			for(i in lambdaopt:.Object@nlambda){
				if(abs(max(.Object@cvmean) - .Object@cvmean[which(alphaseq==alphaopt),i]) > se){
					lambdaopt <- i-1
					break;
				}else if(i==.Object@nlambda){
					lambdaopt = .Object@nlambda
				}
			}
		}else if(criterion %in% c("alpha.1se", "a.1se", "a1se")){ # Instead of using conservative lambda, we'll use conservative alpha moving towards LASSO within 1 standard error of minimum
			opt <- which(.Object@cvmean == max(.Object@cvmean), arr.ind=TRUE)
			alphaopt <- as.numeric(rownames(opt)[1])
			lambdaopt <- opt[1,2] # The only optimum, second column is the lambda value
			# Start walking so as long as cvmean is within 1se of the minimal cvmean error (i.e. maximum iAUC or c-index)
			# This will converge towards a LASSO-like model as long as we stay within 1se of the CV optimum
			# Standard error
			se <- .Object@cvstdev[which(alphaseq==alphaopt), lambdaopt]/sqrt(folds)
			# Step from current minimum until we reach intercept model or min-|1se| < min-error criterion
			for(i in which(alphaseq==alphaopt):length(.Object@alphaseq)){
				if(abs(max(.Object@cvmean) - .Object@cvmean[i,lambdaopt]) > se){
					alphaopt <- i-1
					break;
				}else if(i==length(.Object@alphaseq)){
					alphaopt = alphaseq[length(.Object@alphaseq)]
				}
			}
		}else{
			stop("Illegal 'criterion' parameter, should be one of: min, lambda.1se, alpha.1se")
		}
		# Try to fit model object; run through the ePCR model fitting procedure even if optimums are not detected
		# If errors occur this allows the user to still customize this vector based on the cv-slots
		if(all(!is.na(c(alphaopt, lambdaopt)))){
			# Fit the actual glmnet/coxnet object based on the obtained cross-validation results
			.Object@fit = glmnet::glmnet(x = as.matrix(.Object@x.expand(x)), y = y, family = "cox", 
				nlambda = nlambda, 
				alpha = alphaopt)
			if(verb>-1) cat("--- Computing AUCs for regularization curves for coefficients --- \n\n")
			# Run integrateRegCurve for the final model fit
			.Object@regAUC = integrateRegCurve(.Object@fit)
			# In some cases the lambda sequence is shorter than indicated by nlambda; in that case pick the last lambda index that was informative
			if(length(.Object@fit$lambda)<lambdaopt) lambdaopt <- length(.Object@fit$lambda)
			# Identified optimum parameters according to the criterion
			.Object@optimum <- c(Alpha = alphaopt, AlphaIndex = which(alphaseq==alphaopt), Lambda = .Object@fit$lambda[lambdaopt], LambdaIndex = lambdaopt)
		}else{
			warning("Could not fit the coxnet object to slot @fit after running the CV - problematic entries in the alpha or lambda optimum (slot @optimum)")
		}
		
		if(verb>-1) cat("--- Generating feature list and dictionary --- \n\n")
		.Object@features = unique(colnames(x))
		.Object@dictionary = vector("list", length(.Object@features))
		names(.Object@dictionary) <- .Object@features
		if(!missing(dictionary)){ # If user has provided a feasible dictionary (by default one is provided alongside ePCR)
			# TODO
		}
		if(verb>-1) cat("--- New PSP object successfully created --- \n\n")
		return(.Object)
	}
)

# Show object that "cat"s object to R terminal if asked
setMethod("show", "PSP",
	function(object){
		cat("PSP ePCR object\n")
		cat("N observations: ", dim(object@x)[1], "\n")
		cat("Optimal alpha: ", object@optimum["Alpha"], "\n")
		cat("Optimal lambda: ", object@optimum["Lambda"], "\n")
		cat("Optimal lambda index: ", object@optimum["LambdaIndex"], "\n")
	}
)

#' PSP-methods
#' @note Please refer to the PSP-class examples for applying these PSP-methods
#' @name PSP-methods
NULL

#' print.PSP: Print general information of PSPs contents to the terminal
#' @rdname PSP-methods
#' @param x Generic x
#' @param ... Additional custom parameters passed on
#' @export
setMethod("print", "PSP",
	function(x, ...){
		cat("PSP ePCR object\n")
		cat("N observations: ", dim(x@x)[1], "\n")
		cat("Optimal alpha: ", x@optimum["Alpha"], "\n")
		cat("Optimal lambda: ", x@optimum["Lambda"], "\n")
		cat("Optimal lambda index: ", x@optimum["LambdaIndex"], "\n")
	}
)

#' plot.PSP: By default the mean CV surface in terms of alpha/lambda is plotted using hamlet-package's hmap-function
#' @rdname PSP-methods
#' @param y Generic y
#' @param bias Bias for skewing the color in heatmap key plotting
#' @export
setGeneric("plot", function(x, y, ...) { standardGeneric("plot") })
#' @rdname PSP-methods
#' @aliases plot,PSP,ANY-method
setMethod("plot", "PSP",
	function(x, y, bias=0.1, ...){
		heatcv(x, bias=bias, ...)
	}
)

#' coef.PSP: Default PSP coef-function extracts only the optimum parameters, not whole lambda-range
#' @rdname PSP-methods
#' @export
setMethod("coef", "PSP",
	function(object){
		#predict.coxnet(object@fit, s = object@optimum["Lambda"], type = "coefficients")
		glmnet::predict.glmnet(object@fit, s = object@optimum["Lambda"], type = "coefficients")
	}
)
#' predict.PSP: Predict for a novel patient from current PSP
#' @rdname PSP-methods
#' @param verb Level of verbosity
#' @export
setMethod("predict", "PSP",
	function(object, type="response", newx, verb=0){
		if(!missing(newx)){
			if(sum(is.na(newx))>0){
				if(verb>-1) cat("--- Missing entries detected in newx, running defined imputation function ---\n\n")
				newx <- impute::impute.knn(newx)$data
			}
			# Novel data prediction
			predict(object@fit, newx=as.matrix(object@x.expand(newx)), type=type, s=object@optimum["Lambda"])
		}else{
			# Demonstrate using the data matrix in stored in ePCR object
			predict(object@fit, newx=as.matrix(object@x.expand(object@x)), type=type, s=object@optimum["Lambda"])
		}
	}
)

###
#	Other interesting utility functions for PSPs, such as plotting
###

#' PSP.KM: Kaplan-Meier with division at a given cutoff point within [0,1]
#'
#' @rdname PSP-methods
#' @param cutoff Cutoff point for division
#' @export
setGeneric("PSP.KM", function(object, ...) { standardGeneric("PSP.KM") })
#' @rdname PSP-methods
#' @aliases PSP.KM,PSP,ANY-method
setMethod("PSP.KM", "PSP",
	function(object, cutoff=0.5){
		class <- as.factor(c("Low", "High")[(normriskrank(predict(object@fit, newx=as.matrix(object@x.expand(object@x)), type="response", s=object@optimum["Lambda"]))>=cutoff)+1])
		surv <- survival::survfit(object@y ~ class)
		plot(surv, col=1:2)
		legend("bottomleft", col=1:2, legend=c("Lower cutoff", "Higher cutoff"), bty="n", lwd=1)
	}
)

#' PSP.PCA: Principal Component Plot of a single PSP, showing 2 principal axes with a colouring if strata have been indicated; newx can also be plotted in relation to fitting data
#'
#' @param type Types of variables to include; recognizes (int)eger, (bin)ary and (num)eric
#' @param shuffle Shuffle plotting order
#' @param z Should data centering and scaling should be conducted
#' @param cex Zooming multiplier
#' @param col Vector of color numbers or names to use for strata
#' @param pch Point type to use (refer to par-function pch-parameter)
#' @docType methods
#' @rdname PSP-methods
#' @export
setGeneric("PSP.PCA", function(object, ...) { standardGeneric("PSP.PCA") })
#' @rdname PSP-methods
#' @aliases PSP.PCA,PSP,ANY-method
setMethod("PSP.PCA", "PSP",
	function(object, # PSP object
		newx, # Whether to plot new data as well
		expanded=TRUE, # Whether to plot x.expand data (with interactions) or just raw x
		type = "all", # Type of PCA variables to plot; "all" includes all, "bin"/"binary" includes 0/1, "int"/"integer" includes ...-1,0,1,2... and "num"/"numeric" includes all non-integer
		shuffle=TRUE, # Whether to shuffle the drawing order of points; helps with very overlapping groups
		z=TRUE, # Should data be centered and scaled; if both TRUE, so called z-trans
		cex=1, col=c("aquamarine", "coral", "royalblue", "black"), pch=16){ # Additional tunable plotting parameters
		if(!missing(newx)){ # Plot data along wth new data
			if(expanded){
				dat <- object@x.expand(rbind(object@x, newx[,colnames(object@x)])) # Include expanded data (i.e. interactions)
			}else{
				dat <- rbind(object@x, newx[,colnames(object@x)]) # Just raw data without x.expand
			}
			# Collect strata information for colouring
			stratas <- c(as.character(object@strata),rep("New data", times=nrow(newx))) # newx is a new strata
		}else{ # Plot only existing data PCA with strata
			if(expanded){
				dat <- object@x.expand(rbind(object@x)) # Include expanded data (i.e. interactions)
			}else{
				dat <- object@x # Just raw data without x.expand
			}
			# Collect strata information for colouring
			stratas <- as.character(object@strata)
		}
		# Plot PCA for variable subsets
		if(type == "all"){ # Include all variables
			included <- 1:ncol(dat)
		}else if (type %in% c("bin", "binary")){ # Only binary variables
			included <- which(unlist(apply(dat, MARGIN=2, FUN=function(z) all(z %in% c(0,1)))))
		}else if (type %in% c("int", "integer")){ # Only integer variables (includes binary)
			included <- which(unlist(apply(dat, MARGIN=2, FUN=function(z) all(round(z,0) == z))))
		}else if (type %in% c("num", "numeric")){ # Only non-integer variables (e.g. doubles etc)
			included <- which(unlist(apply(dat, MARGIN=2, FUN=function(z) !all(round(z,0) == z))))
		}else{ # Otherwise fall back to all included
			included <- 1:ncol(dat)
		}
		# Z-score normalization
		if(z) dat <- zt(dat, addz=FALSE)
		# Perform PCA
		pca <- stats::prcomp(dat[,included])
		# Shuffle plotting order of patients and plot 2 main principal components
		if(shuffle){
			ord <- sample(1:nrow(pca$x))
			plot(pca$x[ord,1:2], pch=pch, cex=cex, col=col[as.factor(stratas)][ord], xlab="PC1", ylab="PC2")
		# Use original (most likely batch-ordered) ordering
		}else{
			ord <- 1:nrow(pca$x)
			plot(pca$x[,1:2], pch=pch, cex=cex, col=col[as.factor(stratas)], xlab="PC1", ylab="PC2")
		}
		# Plot legend
		legend("topright", pch=16, col=col[1:length(unique(stratas))][as.numeric(unique(as.factor(stratas)[ord]))], unique(stratas))
	}
)
#' PSP.BOX: Boxplot of a single variable in a PSP in respect to strata, for outlier detection and observing variable distributions
#'
#' @rdname PSP-methods
#' @param var Name of variable to plot
#' @param expanded Should data matrix expansion through interactions be included
#' @export
setGeneric("PSP.BOX", function(object, ...) { standardGeneric("PSP.BOX") })
#' @rdname PSP-methods
#' @aliases PSP.BOX,PSP,ANY-method
setMethod("PSP.BOX", "PSP",
	function(object, newx, var=colnames(object@x)[1], expanded=FALSE){
		if(!missing(newx)){ # Plot data along wth new data
			if(expanded){
				dat <- object@x.expand(rbind(object@x, newx[,colnames(object@x)])) # Include expanded data (i.e. interactions)
			}else{
				dat <- rbind(object@x, newx[,colnames(object@x)]) # Just raw data without x.expand
			}
			# Collect strata information for colouring
			stratas <- c(as.character(object@strata),rep("New data", times=nrow(newx))) # newx is a new strata
		}else{ # Plot only existing data PCA with strata
			if(expanded){
				dat <- object@x.expand(rbind(object@x)) # Include expanded data (i.e. interactions)
			}else{
				dat <- object@x # Just raw data without x.expand
			}
			# Collect strata information for colouring
			stratas <- as.character(object@strata)
		}
		# Boxplot
		boxplot(as.numeric(dat[,var]) ~ as.factor(stratas), range=0)
	}
)

#' PSP.CSP: Cumulative survival probabilities
#'
#' @docType methods
#' @rdname PSP-methods
#' @param object PSP-object
#' @param newx New data matrix
#' @param t Sequence of time points to evaluate cumulative survival probabilities at
#' @param plot Plot the corresponding functionality
#' @export
setGeneric("PSP.CSP", function(object, ...) { standardGeneric("PSP.CSP") })
#' @rdname PSP-methods
#' @aliases PSP.CSP,PSP,ANY-method
setMethod("PSP.CSP", "PSP",
	function(object, newx, t=seq(from=1, to=36*30.5, by=1), plot=FALSE){
		if(!missing(newx)){
			# Predict for the new provided data matrix x
			x <- object@x.expand(newx)
		}else{
			# Predict for the model matrix x
			x <- object@x.expand(object@x)
		}
		csp = TimeSurvProb(object@fit, time=object@y[,"time"], event=object@y[,"status"], s=object@optimum["Lambda"], olddata=object@x.expand(object@x), newdata=x, times=t, plot=plot)
		csp
	}
)

#' PSP.NA: Nelson-Aalen with time-to-event prediction at point t = F^-1(0.5)
#'
#' @docType methods
#' @rdname PSP-methods
#' @export
setGeneric("PSP.NA", function(object, ...) { standardGeneric("PSP.NA") })
#' @rdname PSP-methods
#' @aliases PSP.NA,PSP,ANY-method
setMethod("PSP.NA", "PSP",
	function(object, newx, plot=TRUE){
		if(!missing(newx)){
			# Predict for the new provided data matrix x
			x <- object@x.expand(newx)
		}else{
			# Predict for the model matrix x
			x <- object@x.expand(object@x)
		}
		beta <- predict(object@fit, s=as.numeric(object@optimum["Lambda"]), type="coefficients")[,1]
		na = NelsonAalen(b=beta, time=object@y[,"time"], events=object@y[,"status"], Xold=object@x.expand(object@x), Xnew=x, tpred=seq(from=0, to=max(object@y[,"time"], na.rm=T), length.out=1000), plot=plot)
		na
	}
)
