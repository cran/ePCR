####
#
# Functions intended mainly for internal use within the ePCR-package
# Teemu Daniel Laajala, teelaa@utu.fi
#
####

###
#
#	GENERAL HELPER FUNCTIONS
#
###

#' Integrate the area over/under the regularization path of a penalized regression model
#'
#' This function evaluates the overall significance of a regularized regression coefficient in a penalized Cox model. It takes into account the whole range of lambda-penalization parameter, and computes the area over or under the regularization curve. This gives more insight into the importance of a regression coefficient over the whole range of lambda, instead of evaluating it at a single optimal lambda point determined typically using cross-validation.
#' 
#' @param fit A regularized regression model fited using glmnet
#' @param weighted Should the regularization curve be weighted by the corresponding lambda (as higher lambda pushes coefficients to zero)
#'
#' @return Integrated area over or under a regularization curve using the trapezoid method from the pracma-package
#'
#' @examples
#' # Exemplify one PSP of the readily fitted ensembles
#' data(ePCRmodels)
#' RegAUC <- cbind(
#'	integrateRegCurve(fit = DREAM@PSPs[[1]]@fit),
#'	integrateRegCurve(fit = DREAM@PSPs[[2]]@fit),
#'	integrateRegCurve(fit = DREAM@PSPs[[3]]@fit)
#' )
#' SortRegAUC <- RegAUC[order(apply(RegAUC, MARGIN=1, 
#' 	FUN=function(z) abs(mean(z)) ), decreasing=TRUE),]
#' colnames(SortRegAUC) <- c(DREAM@PSPs[[1]]@description, 
#'	DREAM@PSPs[[2]]@description,
#'	DREAM@PSPs[[3]]@description)
#' SortRegAUC[1:10,] # Top 10 coefficients according to (absolute) regularization curve auc
#' @importFrom pracma trapz
#' @export
# Integrate the area under a regularization curve (coefficient estimates as a function of lambda) in a glmnet object
integrateRegCurve <- function(fit, weighted=FALSE){
	coefsmat = as.matrix(fit$beta) # Coefficient (rows) matrices as a function of lambda (cols)
	lambdas = fit$lambda # Tested penalization values (lambda)
	# Integrate each coefficient (row) separately
	apply(coefsmat, MARGIN=1, FUN=function(z){
		# Add polygon points at the x=0 line at the left-most and right-most point
		xs <- c(lambdas[1], lambdas, lambdas[length(lambdas)])
		if(weighted){
			# Weight by lambda when integrating
			ys <- c(0, z*lambdas, 0)
		}else{
			ys <- c(0, z, 0)
		}
		# From polyarea:  Areas to the left of the vertices are positive, those to the right are counted negative.
		# flip x/y axes so the interpretation is in respect to x=0 axis
		# -> switch to trapezoidal approximation
		-pracma::trapz(xs, ys) # Take negative; AUC should be positive for curves above x>0, negative for curves x<0
	})
}

#' Bootstrapped testing of regression coefficients in a penalized model
#'
#' The purpose of this function is to evaluate a p-value-like statistic for penalized regression coefficients. A fixed number of bootstrapped datasets are generated, and the model coefficients are fitted to these bootstrapped datasets using the pre-determined lambda.
#' 
#' @param fit A regularized regression model fit as provided by the glmnet-package
#' @param lambda The pre-fixed corresponding optimal lambda value, typically determined using cross-validation (e.g. cv.glmnet$lambda.1se or cv.glmnet$lambda.min in glmnet)
#' @param boot The number of bootstrapped datasets to generate
#' @param epsilon The tolerance around beta = 0 to still count estimates as zero
#' 
#' @return Significance values for regression coefficients, defined as the proportion of bootstrapped model fits where coefficient did not shrink within epsilon of zero or where it did not flip sign.
#' @note Notice that this is a highly experimental function, and that many statisticians argue that computing p-values does not make sense for penalized models. The null hypothesis is not well defined, as the bias (regularization) pushes the regression coefficients towards zero. Therefore the null hypothesis is not known and the interpretation is not the conventional regression coefficient p-value.
#' @examples
#' \dontrun{
#' # Computationally too intensive to run bootstrapped fits <5s
#' data(TYKSSIMU)
#' library(survival)
#' x <- as.matrix(xMEDISIMU)
#' y <- yMEDISIMU[,"surv"]
#' nlambda <- 30
#' psp1 <- new("PSP", alphaseq=c(0, 0.5, 1), nlambda = nlambda, folds = 3, x = x, y = y, seeds = 1)
#' .Object <- psp1
#' alphaopt <- psp1@optimum["Alpha"]
#' bs <- bootstrapRegCoefs(fit = psp1@fit, lambda = psp1@optimum["Lambda"], boot = 100)
#' # Histogram of bootstrapped ps
#' hist(bs$ps, breaks=100)
#' }
#' @import glmnet
#' @export
bootstrapRegCoefs <- function(fit, lambda, boot=1000, epsilon = 10^-6){
	# Original data objects need to be present in the environment, as they're not stored in the fitted object
	x = eval(expr=as.list(fit$call)$x)
	y = eval(expr=as.list(fit$call)$y)
	betas = glmnet::predict.coxnet(fit, type="coefficients", s=lambda)[,1]
	# Construct a bootstrap model fit call
	callboot = fit$call
	callboot$x = quote(xboot)
	callboot$y = quote(yboot)
	# Run boot-count of bootstrapped datasets, columns are bootstrap runs and rows are the coefficients at given lambda
	boots = do.call("cbind", lapply(1:boot, FUN=function(z){
		r = sample(1:nrow(x), replace=T) # Bootstrap indices
		xboot = x[r,] # Bootstrap the data matrix X
		if(is.vector(y)){ # Bootstrap the response, with either a vector or some other format (such as 2-column surv object)
			yboot = y[r] # Conventional response vector
		}else{
			yboot = y[r,] # Possibly a survival object or such
		}
		rfit = eval(expr=callboot)
		glmnet::predict.coxnet(rfit, type="coefficients", s=lambda)[,1]
	}))
	# Compute p-value like statistics
	ps = unlist(lapply(1:nrow(boots), FUN=function(z){
		sum(abs(boots[z,]) < epsilon | (betas[z]<0) != (boots[z,]<0))/boot # Proportion of coefs within epsilon of zero or flipped sign
	}))
	names(ps) = names(betas)
	list(boots = boots, ps = ps)
}

###
#
#	FUNCTIONS INTENDED MAINLY FOR INTERNAL USE BUT STILL EXPORTED
#
###

#' Extended function for z-transformation, filling non-finite values and changes column names at will
#'
#' An extended function of the standard z-score standardization of a vector in R (i.e. function 'scale'). Supports filling in non-finite values as well as re-naming variables to distinguish them from non-standardized variables.
#'
#' @param x A data matrix for which the columns are to be standardized
#' @param fillfinite The value to fill non-finite values with, by default zero.
#' @param addz Boolean indicating whether letter 'z' should be appended to the variable names to indicate the standardization
#' @param saveattr Boolean for if an 'attr' should be attached to the standardized vector, similar to how the R default function 'scale' conserves the centering and scaling values
#' @return z-score standardized values (zero mean and unit variation), with non-finite values imputed by zero by default.
#' @examples
#' somedata <- cbind(rnorm(100), runif(100))
#' normdata <- zt(somedata)
#' head(normdata)
#' apply(normdata, MARGIN=2, FUN=mean)
#' apply(normdata, MARGIN=2, FUN=sd)
#' @importFrom stats sd
#' @export
zt <- function(x, 
	fillfinite = 0, # If there are missing or non-finite values in the matrix, should a certain value be used to impute these
	addz = T, # Whether "z" should be added to the column name to indicate it's been z-transformed
	saveattr = T # Whether mean and sd should be saved as an attribute
){
	z <- as.matrix(apply(x, MARGIN=2, FUN=function(z){
		(z - mean(z, na.rm=T))/stats::sd(z, na.rm=T)
	}))
	if(!is.na(fillfinite) & !is.null(fillfinite)) z[!is.finite(z)] <- fillfinite
	if(addz & !is.null(colnames(z))) colnames(z) <- paste(paste("z", colnames(z), sep=""))
	z
}

#' Compute all pairwise interactions between the columns of a data matrix
#'
#' The function multiplies the columns (variables) of a matrix or a data.frame with each other, and produces a new matrix where all pairwise interactions are present. This also includes multiplying a column with its self, thus effectively returning a squared column.
#' 
#' @param input A data matrix (of class matrix or data.frame) for which all column-wise multiplications are to be computed
#' @return A matrix where columns of the original data matrix have been multiplied, indicating column names coupled with a colon in-between
#' @examples
#' set.seed(1)
#' somedata <- data.frame(a = rnorm(10), b = rnorm(10), c = runif(10), d = runif(10))
#' somedata
#' allinteract <- interact.all(somedata)
#' allinteract
#' @export
interact.all <- function(input){
	output <- do.call("cbind", lapply(1:ncol(input), FUN=function(z){ 
		do.call("cbind", lapply(z:ncol(input), FUN=function(x){
			tmp <- data.frame(input[,z] * input[,x])
			colnames(tmp)[1] <- paste(colnames(input)[z], ":", colnames(input)[x], sep="") # Multiplication symbol changed to : as per R style instead of x to avoid confusion with variables containing x
			tmp
		}))
	}))
	output
}

#' Compute a chosen set of pairwise interactions between two sets of columns in a data matrix
#'
#' Similar to interact.all-function, but here user provides two sets of variables, and each pairwise combination between these two sets is multiplied. These pairwise interactions are then returned as a new data matrix, with a colon indicating which variables were multiplied.
#'
#' @param input The input data matrix, of either class matrix or data.frame
#' @param first The first set of columns to combine with each of the members of the second set, as either integers or column names
#' @param second The second set of columns to combine with each of the members of the first set, as either integers or column names
#' @return A data matrix with multiplied columns as indicated using the sets 'first' and 'second'
#' @examples
#' set.seed(1)
#' somedata <- data.frame(a = rnorm(10), b = rnorm(10), c = runif(10), d = runif(10))
#' somedata
#' someinteract <- interact.part(somedata, first = c("a", "b"), second = c("c", "d"))
#' someinteract
#' @export
interact.part <- function(input, first, second){
	output <- do.call("cbind", lapply(first, FUN=function(z){ 
		do.call("cbind", lapply(second, FUN=function(x){
			tmp <- data.frame(input[,z] * input[,x])
			colnames(tmp)[1] <- paste(z, ":", x, sep="") # Multiplication symbol changed to : as per R style instead of x to avoid confusion with variables containing x
			tmp
		}))
	}))
	output
}



#' Scoring function for evaluating survival prediction through concordance index (c-index)
#' 
#' C-index (Concordance index) of the predicted vs. true answer, i.e. proportion of pairs that go in correct direction over all pairwise comparisons
#' 
#' @param pred Numeric risk score for each event
#' @param time A vector of event or censoring times
#' @param event A binary valued vector that indicates either death (1) or right-censoring (0)
#' @param real A previously constructed Surv-object instead of providing time and event
#' @return Corcordance index (c-index) of the prediction
#
# Pred should be a numeric risk score
# 'time' and indicates censoring or event time, while 'event' indicates either death (1) or right-censoring (0)
#' @importFrom survival Surv concordance
#' @examples
#' # A random prediction ought to be near 0.5
#' # c-index is not sensitive to time scale, as it tests pairwise prediction accuracy
#' set.seed(1); prediction <- sample(1:20)
#' time <- seq(from=1000, to=50, by=-50)
#' event <- c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1)
#' library(survival)
#' score.cindex(pred = prediction, real = Surv(time=time, event=event))
#' @export
score.cindex = function(pred, time, event, real){ 
	# Survival object for the observed data
	if(!missing(real) & !class(real)=="Surv"){
		surv <- survival::Surv(time, event)
	}else{
		surv <- real
	}
	# Compute c-index object
	cindex <- survival::concordance(surv ~ pred)$concordance
	# Return concordance index (c-index)
	cindex	
}

#' Scoring function for evaluating survival prediction by time-wise integrated AUC
#' 
#' Time-wise integrated prediction for survival is performed by this scoring function using the timeROC-package. 
#' It's offered as an alternative to the score.cindex-function with the difference that time-wise integrated AUC is sensitive to the choice of time-window. 
#' By default (as similar to DREAM 9.5 mCRPC challenge), the AUCs are determined at 6 to 30 months, and the AUC is then normalized to a score within [0,1]. Notice that for studies shorter or longer than this proposed time window, the scoring function should be adjusted accordingly.
#' 
#' @param pred Numeric risk score for each event
#' @param time A vector of event or censoring times
#' @param event A binary valued vector that indicates either death (1) or right-censoring (0)
#' @param real A previously constructed Surv-object instead of providing time and event
#' @param times Time-points at which to evaluate the iAUC
#' @return The integrated area under the ROC-curve over time
#' 
#' @importFrom timeROC timeROC
#' @importFrom Bolstad2 sintegral
#' @importFrom survival Surv
# 24-month time-integrated AUC (iAUC) as in DREAM competition, courtesy of challenge organized and modified by TDL
# Integrated between 6 to 30 months thus 24 month period
# Pred should be a numeric risk score
# 'time' and indicates censoring or event time, while 'event' indicates either death (1) or right-censoring (0)
#' @examples
#' # A random prediction ought to be near 0.5 
#' # iAUC is sensitive to the choice of time points to test AUC at
#' set.seed(1); prediction <- sample(1:20)
#' time <- seq(from=1000, to=50, by=-50)
#' event <- c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1)
#' library(survival)
#' score.iAUC(pred = prediction, real = Surv(time=time, event=event))
#' @export
score.iAUC = function(pred, time, event, real, 
	# compute iAUC from 6 to 30 months; sequence of months
	times = seq(6,30,by=1) * 30.5){
	# Survival object for the observed data
	if(!missing(real) & !class(real)=="Surv"){
		surv <- survival::Surv(time, event)
	}else{
		surv <- real
	}
	# Compute timeROC AUCs at time points and then integrate area under the curve
	if(class(real)=="Surv"){ # Open up a Surv-object if such was provided
		event = real[,"status"]
		time = real[,"time"]
	}
	aucs <- timeROC::timeROC(T=time,
		  delta=event,
		  marker=pred,
		  cause=1,
		  weighting="marginal",
		  times=times,
		  iid=FALSE)$AUC
	# Simpsons rules for integrating under curve
	iAUC <- Bolstad2::sintegral(times, aucs)$int / (max(times) - min(times))
	# Return integrated AUC between the month-sequence in interval 6 months to 30 months by default
	iAUC
}

#' Function that creates customized cross-validation folds
#'
#' The function creates two matrices and returns them as members of a list. 'train' holds training sample indices, columns are cv-folds and rows are sample indices to hold as training samples in each fold. 'test' holds test/validation sample indices, columns are cv-folds and rows are sample indices to hold as test samples in each fold
#' 
#' @param x The original data matrix or data.frame
#' @param fold Number of desired cross-validation folds; preferably between 3 (3-fold CV) and nrow(x) (LOO-CV)
#' @param strata Indicator if some strata should be balanced over the bins; if no balancing is required, the vector should consist of a single value with length equal to rows in x. Otherwise each strata/batch should be indicated as a unique member in the vector.
#' @examples 
#' data(TYKSSIMU)
#' cvfolds <- cv(x = xMEDISIMU, fold = 3)
#' cvfolds$train
#' cvfolds$test
#' @param shuffle Whether the indices for the data matrix should be shuffled prior to assigning them to train/test bins
#' @param seed A random seed for reproducibility
#' @export
cv <- function(
	# Original x data frame
	x,
	# Number of CV-folds
	fold = 10,
	# Should some strata be balanced over the bins? By default no, should be a vector equal to the number of rows in x
	strata = rep(1, times=nrow(x)),
	# Should data be randomized when creating cv-folds in addition to allocating to bins
	shuffle = TRUE,
	# Seed number for reproducibility
	# If NULL then seed is not set
	seed = NULL
){
	# Seed number
	if(!is.null(seed)) set.seed(seed)

	# Allocate folds
	uniqs <- unique(strata)
	folds <- rep(1:fold, times=ceiling(nrow(x)/fold))[1:nrow(x)]
	ifelse(shuffle,
		ord <- unlist(lapply(uniqs, FUN=function(z) sample(which(strata==z)))),
		ord <- unlist(lapply(uniqs, FUN=function(z) which(strata==z))))
	
	whichfold <- vector(length=nrow(x))
	whichfold[ord] <- folds

	# Construct test and train sets based on the folds
	test <- lapply(1:fold, FUN=function(z) which(z==whichfold))
	train <- lapply(1:fold, FUN=function(z) do.call("c", test[-z]))
	
	# Return cv sample indices per each fold for train and test sets respectively
	list(train = train, test = test)
}


#' Cross-validation runs for risk predition at a single value of alpha
#'
#' Run n-fold cross-validation for a chosen prediction metric at a single value of the L1/L2 norm alpha. A suitable lambda sequence is determined by glmnet, and the cross-validation returns a prediction matrix over the folds over various lambda. This function is mostly called by the higher hierarchy functions, such as cv.grid, which allows varying also the alpha-parameter.
#'
#' @param x The data matrix to use for predictions
#' @param y The response for coxnet; preferably a preconstructed Surv-object
#' @param folds Number of cross-validation folds
#' @param alpha Chosen L1/L2 norm parameter lambda
#' @param nlamb Number of lambda values 
#' @param verb Integer indicating level of verbosity, where 0 is silent and 1 provides additional information
#' @param scorefunc Chosen scoring function, e.g. score.cindex or score.iAUC
#' @param plot Should a CV-performance curve be plotted as a function of lambda, indicating min/max/mean/median of CV performance over the folds
#' @return A matrix of cross-validation scores, where rows correspond to CV folds and columns to various lambda values chosen by glmnet
#' @examples 
#' data(TYKSSIMU)
#' library(survival)
#' ydat <- Surv(event = yMEDISIMU[,"DEATH"], time = yMEDISIMU[,"LKADT_P"])
#' set.seed(1)
#' cvs <- cv.alpha(x = xMEDISIMU, y = ydat, alpha = 0.5, folds = 5, 
#' 	nlamb = 50, verb = 1, scorefunc = score.cindex, plot = TRUE)
#' cvs
#' @importFrom stats median
#' @export
cv.alpha <- function(
	# Data matrix (no missing values, and readily determined features to try with glmnet)
	x,
	# Response Y vector
	y,
	# Count of cross-validation folds
	folds = 10,
	# A single alpha value
	alpha = 0.5,
	# Number of lambda values to vary within each alpha; lets glmnet determine this suitable lambda sequence by itself
	nlamb = 100,
	# Level if verbosity, verb>=1 indicates debugging
	verb = 0,
	# Which score function to use (readily provided funcs: score.cindex, score.iAUC)
	scorefunc,
	# By default no plotting at a single lambda value is not done, although this may be very interesting in many cases
	plot = FALSE
){
	# Indexing of cross validation, build the inner loops
	crossval <- cv(x=x, fold=folds) 
	# Change factory settings in glmnet to make sure whole lambda sequence is run
	#glmnet.control(fdev = 0, devmax = 1)

	if(verb>0) print("Determining suitable lambda-sequence...")
	
	# Determine a suitable lambda sequence for all data
	if(class(y)=="Surv"){
		lamb <- glmnet::glmnet(y=y, x=as.matrix(x), family="cox", alpha=alpha, nlambda=nlamb)$lambda
	}else{
		lamb <- glmnet::glmnet(y=Surv(time=y[,"LKADT_P"], event=y[,"DEATH"]), x=as.matrix(x), family="cox", alpha=alpha, nlambda=nlamb)$lambda
	}

	# The inner loop is used for conventional 10-fold CV (or however many were specified by the user)
	res <- do.call("rbind", lapply(1:folds, FUN=function(fold){
		# Build the inner loop data matrices
		train.x <- x[crossval$train[[fold]],]
		test.x <- x[crossval$test[[fold]],]
		# Inner test
		train.y <- y[crossval$train[[fold]],]
		test.y <- y[crossval$test[[fold]],]

		if(verb>0) print("Running cross-validation loop instance...")

		## GLMNET FIT	
		if(class(y)=="Surv"){
			fit <- glmnet::glmnet(y=train.y, x=as.matrix(train.x), family="cox", alpha=alpha, lambda=lamb)
		}else{
			fit <- glmnet::glmnet(y=Surv(time=train.y[,"LKADT_P"], event=train.y[,"DEATH"]), x=as.matrix(train.x), family="cox", alpha=alpha, lambda=lamb)		
		}
		
		if(verb>0) print("Predicting from the CV fit...")

		## INNER FOLD PREDICTION
		pred <- glmnet::predict.coxnet(fit, newx=as.matrix(test.x), s=lamb, type="response")
		if(verb>0) print("Obtaining score for the CV fit in predicting the left-out part...")

		# Compute prediction scoring per each lambda prediction
		apply(pred, MARGIN=2, FUN=function(pred.y){
			# Compute iAUC / c-index / itc according to the scoring script provided by the user
			score <- scorefunc(
				pred = pred.y,
				real = test.y
			)
		})
	}))
	colnames(res) <- paste("L_", 1:length(lamb),sep="")
	rownames(res) <- paste("CV_",1:folds,sep="")
	if(verb>0) print(paste("Alpha", alpha, "run, plotting min/mean/median/max CV score curve"))
	# Reset glmnet settings to factory
	#glmnet.control(factory = TRUE)

	# Should a mean curve over the folds be plotted
	# Allows min/max as well
	if(plot){ 
		# Plot mean curve
		plot(x=log(lamb), 
			y=apply(res, MARGIN=2, FUN=function(z) mean(unlist(z))), 
			xlab="log(Lambda)", ylab="CV Error", type="l", col="red", ylim=c(0, 1), lwd=2)
		# Median over CVs
		points(x=log(lamb), y=apply(res, MARGIN=2, FUN=function(z) stats::median(unlist(z))), type="l", lwd=2, col="orange")
		# Max over CVs
		points(x=log(lamb), y=apply(res, MARGIN=2, FUN=function(z) max(unlist(z))), type="l", lwd=1, col="blue")
		# Min over CVs
		points(x=log(lamb), y=apply(res, MARGIN=2, FUN=function(z) min(unlist(z))), type="l", lwd=1, col="blue")
		legend("bottom", horiz=T, cex=0.8, col=c("blue", "orange", "red"), lwd=c(1,2,2), legend=c("Min/Max", "Median", "Mean"), bty="n")
	}	
	res
}

#' Cross-validation runs for risk predition for a grid of predetermined alpha values and their conditional lambda values
#'
#' Expanded Cross-Validation function to run the whole CV in the lambda/alpha grid instead of just lambda-sequence with a pre-specified alpha
#'
#' @param alphaseq Sequence of alpha values to test, which should be within [0,1] (with alpha = 0 being ridge regression, 0 < alpha < 1 being elastic net, and alpha = 1 being LASSO)
#' @param seed Random number generation seed for reproducibility
#' @param x Data matrix x
#' @param y The Surv-object response y
#' @param folds Number of folds in the cross-validation
#' @param nlamb Number of lambda values to test in each alpha; notice that these lambda values vary conditional to alpha
#' @param verb Level of verbosity, with 0 as silent and 1 with additional output
#' @param scorefunc Chosen scoring function, e.g. score.cindex or score.iAUC
#' @param plot Whether a performance should be plotted at each varying alpha-value similar to cv.alpha-plots
#' @return List of matrices of cross-validation performance values over the alpha/lambda grid for mean/median/min/max/stdev of the chosen performance metric, with rows indicating various alpha-values and columns indicating lambda-values.
#' @examples
#' data(TYKSSIMU)
#' library(survival)
#' ydat <- Surv(event = yMEDISIMU[,"DEATH"], time = yMEDISIMU[,"LKADT_P"])
#' cvs <- cv.grid(x = xMEDISIMU, y = ydat, folds = 3, nlamb = 30, alphaseq = seq(0, 1, by=5), 
#' 	scorefunc = score.iAUC, plot = TRUE, seed = 1)
#' cvs
#' @importFrom stats median
#' @export
cv.grid <- function(
	# Sequence of alpha values
	alphaseq = seq(from=0, to=1, by=.1),
	# Should a seed number be set, highly recommended albeit not done by default
	seed,
	## Majority of the rest are just propagated as-is to the cv.alpha-function 
	# Data matrix (no missing values, and readily determined features to try with glmnet)
	x,
	# Response Y vector
	y,
	# Count of cross-validation folds
	folds = 10,
	# Number of lambda values to vary within each alpha; lets glmnet determine this suitable lambda sequence by itself
	nlamb = 100,
	# Level if verbosity, verb>=1 indicates debugging
	verb = 0,
	# Which score function to use (readily provided funcs: score.cindex, score.iAUC)
	scorefunc,
	# By default no plotting at a single lambda value is not done, although this may be very interesting in many cases
	plot = FALSE
){
	cv.list <- list()
	for(a in alphaseq){
		#if(!missing(seed)) set.seed(seed) # Use the same binning for each alpha if seed is desired
		set.seed(seed) # Use the same binning for each alpha if seed is desired
		if(verb>=0) print(paste("alpha", a))
		cva <- cv.alpha(x=x, y=y, alpha=a, 
			folds = folds, nlamb = nlamb, 
			scorefunc = scorefunc, verb = verb, plot = plot)
		# Sometimes lambda sequence runs short, repeat the last lambda with non-zero parameters to get a full compatible matrix
		while(ncol(cva)<nlamb){
			cva <- cbind(cva, REP = cva[,ncol(cva)])
		}
		cv.list[[length(cv.list)+1]] <- cva

		names(cv.list)[length(cv.list)] <- a
	}
	# Return matrices of key statistics at different places such as alpha/lambda CV mean, median, min, and max
	list(mean = .cvbind(cv.list, utilfunc=mean), median = .cvbind(cv.list, utilfunc=stats::median), min = .cvbind(cv.list, utilfunc=min), max = .cvbind(cv.list, utilfunc=max), stdev = .cvbind(cv.list, utilfunc=sd))
}

###
#
#	STRICTLY INTERNAL FUNCTIONS
#
###

# c060 internal base survival function
#
# Function below is taken from within the c060 R-package since it is not exported but its functionality was required.
# Rights to below code belongs to the c060-package maintainers and contributors.
# Concordantly, it is kept strictly internal function also within the ePCR-package.
.basesurv <- function (response, lp, times.eval = NULL, centered = FALSE) 
{
    if (is.null(times.eval)) 
        times.eval <- sort(unique(response[, 1]))
    t.unique <- sort(unique(response[, 1][response[, 2] == 1]))
    alpha <- length(t.unique)
    for (i in 1:length(t.unique)) {
        alpha[i] <- sum(response[, 1][response[, 2] == 1] == 
            t.unique[i])/sum(exp(lp[response[, 1] >= t.unique[i]]))
    }
    obj <- stats::approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, 
        rule = 2)
    if (centered) 
        obj$y <- obj$y * exp(mean(lp))
    obj$z <- exp(-obj$y)
    names(obj) <- c("times", "cumBaseHaz", "BaseSurv")
    return(obj)
}

# Bind a list of different alpha-parameter runs
.cvbind <- function(list, utilfunc=mean){
	tmp <- do.call("rbind", lapply(list, FUN=function(z){
		apply(z, MARGIN=2, FUN=function(x) utilfunc(unlist(x)))
	}))
	rownames(tmp) <- names(list)
	tmp
}

# Trim 'fat' of glmnet/coxnet objects to conserve space; some functionality may not be any more usable as-is after this
.trim.coxnet <- function(
	# A glmnet-object, preferably coxnet
	object
	){

	# Return trimmed object
	object
}

# Trim 'fat' of PEP/PSP objects to conserve space; some functionality may not be any more usable as-is after this
.trim.ePCR <- function(
	# PSP or PEP
	object,
	# Should cross-validation mean be trimmed; heatmaps cannot be plotted if so
	rm.cvmean = F
	){
	
	# Function for specificly trimming psp
	trim <- function(psp){
	
	}
	
	if(class(object)=="PSP"){
		object <- trim(object)
	}else if(class(object)=="PEP"){
	
	}else{
		stop(paste("object input should be either a PSP or a PEP object, current object class:", class(object)))
	}

	# Return trimmed object
	object
}



