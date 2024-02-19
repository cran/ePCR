####
#
#	Functions intended mainly for external user-based use of the ePCR-package
#	Teemu Daniel Laajala, teelaa@utu.fi
#
####

###
#
#	Functions for ensemble-level operations	
#
###

#' Compute mean of predicted risk ranks for an ePCR ensemble
#'
#' @param x A list or a matrix of risk scores given per each ensemble member (each column or list member is considered an equal member of the ensemble)
#'
#' @return An averaged predicted risk rank over all the ensemble members
#' 
#' @note Extensively called by the 'predict'-function for PEP-objects when risk predictions are performed over the ensemble
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @export
meanrank = function(x){ # x is a list or matrix of risk scores per ensemble member (column = ensemble member)
	if("list" %in% class(x)){
		apply(do.call("cbind", lapply(x, FUN=rank)), MARGIN=1, FUN=mean)
	}else if("matrix" %in% class(x)){
		apply(apply(x, MARGIN=2, FUN=rank), MARGIN=1, FUN=mean)
	}else{
		stop(paste("Error in mean rank computation for x, invalid class:", paste(class(x), collapse=" ")))
	}
}

#' Normalize ensemble risk scores to ranks and then to uniform range
#'
#' @param x A list or a matrix of risk scores given per each ensemble member (each column or list member is considered an equal member of the ensemble 
#'
#' @return An averaged predicted risk rank over all the ensemble members that has been normalized to the range [0,1] based on: (x - min(x)) / (max(x) - min(x)) -> [0,1]
#' 
#' @note Normalizes 'predict'-function calls for PEP-objects after calling 'meanrank'-function
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @export
normriskrank <- function(x){
	x <- meanrank(x)
	(x - min(x)) / (max(x) - min(x))
}

#' Conform the dimensions of a new input data matrix to a readily fitted PEP or PSP object
#'
#' @param object A readily fitted PSP or PEP object
#' @param newx A data matrix or a data.frame which to expand
#' @return An expanded data matrix for which the dimensions conform to the regression coefficients in the PSP or PEP
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @export
conforminput <- function(object, newx){
	if(inherits(object,"PSP")){
		#feats <- rownames(glmnet::predict.coxnet(object@fit, s=object@optimum["Lambda"], type="coefficients"))
		feats <- rownames(predict(object@fit, s=object@optimum["Lambda"], type="coefficients"))
		expandx <- object@x.expand(as.matrix(newx))
		comformx <- as.data.frame(lapply(feats, FUN=function(z){
			if(z %in% colnames(expandx)){
				expandx[,z]
			}else{
				rep(0, times=nrow(newx))
			}
		}))
		colnames(comformx) <- feats
		as.matrix(comformx[,feats])
	}else if(inherits(object,"PEP")){
		#feats <- rownames(glmnet::predict.coxnet(object@PSPs[[1]]@fit, s=object@PSPs[[1]]@optimum["Lambda"], type="coefficients"))
		feats <- rownames(predict(object@PSPs[[1]]@fit, s=object@PSPs[[1]]@optimum["Lambda"], type="coefficients"))
		expandx <- object@PSPs[[1]]@x.expand(as.matrix(newx))
		comformx <- as.data.frame(lapply(feats, FUN=function(z){
			if(z %in% colnames(expandx)){
				expandx[,z]
			}else{
				rep(0, times=nrow(newx))
			}
		}))
		colnames(comformx) <- feats
		as.matrix(comformx[,feats])
	}else{
		stop("'object' should be either a PSP or a PEP S4-object")
	}
}


###
#
#	Functions for prediction and relevant operations
#
###

#' Predict cumulative survival probabilities for new data at given time points
#'
#' Given a readily fitted regularized Cox regression model, this function predicts the cumulative survival probabilities for new data at time points determined by the user. The function uses c060-package's functionality for computing base hazard, and then performs linear predictions for new observations using the fitted regularized Cox regression model.
#'
#' @param fit A single regularized Cox regression model fitted using glmnet
#' @param time Time to events for the training data
#' @param event Event indicators for the training data (0 censored, 1 event)
#' @param olddata The old data matrix used to fit the original 'fit' glmnet-object
#' @param newdata The new data matrix for which to predict time-to-event prediction (should comform to the old data matrix)
#' @param s The optimal lambda parameter as used in the glmnet-package for its fit objects
#' @param times The time points at which to estimate the cumulative survival probabilities (by default in days)
#' @param plot Should the cumulative survival probabilities be plotted as a function of time
#'
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @return Cumulative survival probabilities at the chosen time points
#' @import survival
#' @import glmnet
#' @export
TimeSurvProb <- function(fit, time, event, olddata, newdata, s, times=c(1:36)*30.5, plot=FALSE){
	# Survival response as 2-column matrix
	sur <- as.matrix(survival::Surv(time=time, event=event))
	# Linear predictor
	#lp1 <- as.numeric(glmnet::predict.coxnet(fit, newx = as.matrix(olddata), s = s, type = "link"))
	lp1 <- as.numeric(predict(fit, newx = as.matrix(olddata), s = s, type = "link"))
	#lp2 <- as.numeric(glmnet::predict.coxnet(fit, newx = as.matrix(newdata), s = s, type = "link"))
	lp2 <- as.numeric(predict(fit, newx = as.matrix(newdata), s = s, type = "link"))
	# Use base cumulative base hazard computed from package c060
	# c060 does not export basesurv-function; code is from c060
	#basesur <- c060:::basesurv(response=sur, lp=lp1, times.eval=times)
	basesur <- .basesurv(response=sur, lp=lp1, times.eval=times)
	# Compute survival probabilities for newdata in the desired time points
	probs <- exp(exp(lp2) %*% -t(basesur$cumBaseHaz))    
	# Should the predicted time-to-event cases be plotted
	if(plot){
		plot.new(); plot.window(xlim=c(0,max(times, na.r=T)), ylim=c(0,1))
		axis(1); axis(2); box(); title(xlab="Time (days)", ylab="Cumulative survival probability")
		apply(probs, MARGIN=1, FUN=function(z){
			points(times, z, type="l", lwd=1)
		})
		abline(v=max(time[event==1], na.rm=T), col="grey", lwd=2)
		legend("bottomleft", col=c("black", "grey"), lwd=c(1,2), legend=c("Individuals", "Max obs event time"))
	}
	# Format return
	colnames(probs) = times
	rownames(probs) = rownames(newdata)
	ranks = rank(probs[,ncol(probs)])
	names(ranks) = rownames(newdata)
	list(ranks = ranks, probs = probs)
}

#' Cox-Oakes extension of the Nelson-Aalen estimates for a Cox model
#'
#' Implementing the heuristic Cox and Oakes extension of the Nelson-Aalen estimate for Cox model to extract individual-specific survival. Time-to-event predictions are then given at the first time point at which an individual reaches an event probability of 50%.
#' 
#' @param b Beta coefficients at optimal glmnet coxnet model (lambda, alpha)
#' @param Xold Data matrix with rows as individuals (i.e. training data)
#' @param Xnew Possible new prediction data matrix (if omitted the training data is used, or if it is a new dataset the columns should comform to the training data and new individuals be provided as rows)
#' @param events Deaths or right-censoring per each individual (1 death 0 alive censored) for Xold
#' @param time Times to event or censoring for Xold
#' @param tpred The predicted time points; more tight grid gives a smoother curve
#' @param plot Should an individualized plot be plotted to show how the cumulative survival curves behave
#'
#' @examples
#' data(TYKSSIMU)
#' library(survival)
#' xdat <- as.matrix(xMEDISIMU)
#' ydat <- yMEDISIMU[,"surv"]
#' @return Predicted times-to-event predictions either for the training data or an optional provided novel dataset
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @note See section 3.6 at http://data.princeton.edu/pop509/NonParametricSurvival.pdf for the reference
#' @export
NelsonAalen <- function(
	b, # Beta coefficients at optimal glmnet coxnet model (lambda, alpha)
	Xold, # Data matrix with rows as individuals (i.e. teaching data)
	Xnew, # Possible new prediction covariates (can be the same teaching data or a new dataset with same columns and new individuals as rows)
	events, # Deaths or not per each individual (1 death 0 alive censored) for Xold
	time, # Times to event or censoring for Xold
	tpred = 0:round(max(time,na.rm=T),0), # The predicted time points; more tight grid gives a smoother curve
	plot=FALSE # Should an individualized plot be plotted to show how the cumulative survival curves behave
){
	# Need for %*%, could be a data.frame otherwise
	Xold <- as.matrix(Xold)
	if(!missing(Xnew)){
		Xnew <- as.matrix(Xnew)
	}else{
		Xnew <- as.matrix(Xold)
	}
	
	# The exponent term exp(x^'_k \hat{\beta}) - notation mixes a lot of i:s and j:s here, and t used here instead for transpose...
	# Pre-compute exponent terms, always pick k:th from vector here-on
	expterms <- unlist(lapply(1:nrow(Xold), FUN=function(k) as.numeric(exp(Xold[k,,drop=F] %*% b))))
	# Deaths AT the exact discrete moment t_i
	# Precompute deaths at t_i:s -> d_i terms
	### JAN '17; handle dis as a list to allow possibility for non-discrete event times
	### Double check with discrete timed data, should produce same result as previously
	dis <- vector("list", length=length(unique(time)))
	names(dis) <- as.character(unique(time))
	
	for(i in names(dis)){
		dis[i] <- sum(events[which(as.character(time) == i)] == 1)
	}
	h0tis <- unlist(lapply(1:nrow(Xold), FUN=function(i){
		###nom <- dis[[time[i]]]
		nom <- as.numeric(dis[as.character(time[i])])
		denom <- sum(expterms[which(time[i] < time)])
		nom/denom
	}))
	# If new predictions desired, re-compute exponential terms
	if(!missing(Xnew)) expterms <- unlist(lapply(1:nrow(Xnew), FUN=function(k) as.numeric(exp(Xnew[k,,drop=F] %*% b))))

	preds <- lapply(tpred, FUN=function(tz){
		# Cumulative hazard function \hat{V}_0 (t)
		V0t <- sum(h0tis[which(time <= tz)])
		# Cumulative survival
		S0t <- exp(-V0t*expterms)
		list(hazards = V0t*expterms, survival = S0t)
	})
	hazards <- do.call("rbind", (lapply(preds, FUN=function(z) z[[1]])))
	survival <- do.call("rbind", (lapply(preds, FUN=function(z) z[[2]])))
	colnames(hazards) <- colnames(survival) <- rownames(Xnew)
	rownames(hazards) <- rownames(survival) <- paste("t_", tpred, sep="")

	# F50s are time points for which the cumulative survival probabilities F(t) = 0.5
	F50s <- apply(survival, MARGIN=2, FUN=function(z){
		tpred[which(z<0.5)[1]]
	})
	F50cens <- F50s
	F50s[is.na(F50s)] <- max(time[events==1], na.rm=T)
	F50s <- survival::Surv(time=F50s, event=!is.na(F50cens))
	
	res <- list(hazards = t(hazards), F50s = F50s, survival = t(survival))
	# Individual survival plots with possibly plot>=2 to indicate predicted 0.5 cumulative survival days
	if(plot){
		plot.new()
		plot.window(xlim=c(0, max(time, na.rm=T)), ylim=c(0,1))
		box(); axis(1); axis(2)
		apply(res$survival, MARGIN=1, FUN=function(z){
			points(tpred, z, type="l", col="black")
		})
		title(xlab="Time", ylab="Cumulative Survival Probability")
		abline(v=max(time[events==1], na.rm=T), col="grey", lwd=2)
		if(as.numeric(plot)>=2){
			abline(h=0.5, lwd=2, col="red")
			lapply(F50s[,"time"], FUN=function(z) points(x=c(z, z), y=c(0, 0.5), col="red", lwd=1, type="b", pch=16))
		}
	}
	# Return list of predicted time-to-events
	res
}


###
#
#	Functions for visualizations
#
###

#' Plot a heatmap of the prediction performance statistic as a function of lambda and alpha combinations
#'
#' This function plots a heatmap of cross-validation results by varying the penalization/regularization parameter (lambda, x-axis), together with the corresponding L1/L2 norm parameter alpha (i.e. LASSO, elastic net, ridge regression). The optimal spot in the parameter grid gives insight into the behavior of the regularization in respect to the norms, but note that the lambda-parameter on x-axis is not constant given a conditional alpha-parameter; rather it is a suitable vector chosen by the glmnet-package.
#'
#' @param psp An S4-class PSP-object to plot, as built using the ePCR-package
#' @param bias Bias in color palette (skews it to favor distinguishing high values better by default)
#' @param by.rownames Show every n:th row name (helps for dense axis labels)
#' @param by.colnames Show every n:th column name (helps for dense axis labels)
#' @param paletcol Names for colours to include in the heatmap palette
#' @param paletncol Number of colours on the color key
#' @param xlab Label for the x-axis (typically log-lambda penalization parameter)
#' @param ylab Label for the y-axis (typically alpha-value indicating LASSO, elastic net or ridge regression)
#' @param main Main label on top of the heatmap
#' @param plot.opt Should the best (highest) performance statistic be indicated as a large dot on the heatmap
#' @param plot.1sd Should boundaries of the optimal performance statistic area be outlined as within 1 standard deviation of the optimal spot (note: experimental). This attempts to mimic the 1sd-optimum suggested in the glmnet-package for cross-validation for a constant alpha parameter but for 2 dimensions.
#' @param ... additional parameters passed on to the hmap-function of hamlet-package
#'
#' @examples
#' data(ePCRmodels)
#' par(mfrow=c(1,3))
#' heatcv(DREAM@PSPs[[1]], main=DREAM@PSPs[[1]]@description, by.rownames=10, by.colnames=10)
#' heatcv(DREAM@PSPs[[2]], main=DREAM@PSPs[[2]]@description, by.rownames=10, by.colnames=10)
#' heatcv(DREAM@PSPs[[3]], main=DREAM@PSPs[[3]]@description, by.rownames=10, by.colnames=10)
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @note The heatmap plotting is compatible with the default plot-region in a R graphic canvas. The function hmap from the same author's hmap-package can be highly customized to fit more specific needs.
#' @import grDevices graphics 
#' @importFrom hamlet hmap hmap.key
#' @export
heatcv <- function(
	psp, # PSP-object to plot
	bias=0.1, # Bias in color palette (skews it to favor distinguishing high values better by default)
	by.rownames = 1, # Show every n:th rowname
	by.colnames = 1, # Show every n:th colname
	paletcol = c("cyan","blue","black","red","orange"), # Colors utilized in the heatmap palette
	paletncol = 1000, # Number of colours for binning in the cv heatmap
	xlab = "Alpha-dependent log-Lambda", # X-axis label
	ylab = "Alpha", # Y-axis label
	main = "", # Main title
	plot.opt = TRUE, # Should optimum be plotted as a point? By default a purple dot
	plot.1sd = FALSE, # Should the boundaries of 1 standard error within optimum be plotted in heatmap
	... # Additional parameters
){
	# The base information is the mean of performation statistic obtained averaged over the cross-validation folds
	mat <- psp@cvmean
	# Heatmap 'hmap' function of the package hamlet
	palet <- colorRampPalette(paletcol, bias=bias)(paletncol)
	# Only show every 5th row/column name in the heatmap
	rownames(mat)[-seq(from=1, to=nrow(mat), by=by.rownames)] <- NA
	colnames(mat)[-seq(from=1, to=ncol(mat), by=by.colnames)] <- NA
	
	# Furthermore, the variation in the cross-validation folds can be utilized to supplement the mean
	# EXPERIMENTAL; the area can be very non-convex, thus good intuitive visualization can be hard to achieve
	if(plot.1sd){	
		mat1sd <- abs(psp@cvmean - psp@cvmean[which(psp@cvmean==max(psp@cvmean), arr.ind=T)])<psp@cvstdev[which(psp@cvmean==max(psp@cvmean), arr.ind=T)]
		mat1sd[which(!mat1sd, arr.ind=T)] <- NA
		mat1sd[which(mat1sd, arr.ind=T)] <- "purple"
		for(i in 2:(nrow(mat1sd)-1)){
			for(j in 2:(nrow(mat1sd)-1)){
				npurple <- 0 + 
					as.numeric(!is.na(mat1sd[i-1,j-1])) + 
					as.numeric(!is.na(mat1sd[i-1,j])) + 
					as.numeric(!is.na(mat1sd[i,j-1])) + 
					as.numeric(!is.na(mat1sd[i+1,j])) + 
					as.numeric(!is.na(mat1sd[i+1,j+1])) + 
					as.numeric(!is.na(mat1sd[i,j+1])) + 
					as.numeric(!is.na(mat1sd[i+1,j-1])) + 
					as.numeric(!is.na(mat1sd[i-1,j+1]))
				if(npurple>4){
					mat1sd[i,j] <- NA	
				}
			}
		}
		border <- mat1sd
	}else{
		border <- matrix(NA, nrow=nrow(psp@cvstdev), ncol=ncol(psp@cvstdev))
	}
	# Pinpoint the optimal performance statistic
	if(plot.opt){
		h <- hamlet::hmap(as.matrix(mat), Rowv=NA, Colv=NA, col=palet, border=border)
		# Point of optimum
		points(
			x=mean(c(h$xmatseq[psp@optimum["LambdaIndex"]+1], h$xmatseq[psp@optimum["LambdaIndex"]])), 
			y=mean(c(rev(h$ymatseq)[which(psp@alphaseq==psp@optimum["Alpha"])+1], rev(h$ymatseq)[which(psp@alphaseq==psp@optimum["Alpha"])])), 
			col="pink", cex=1.5, pch=16)
		legend("topright", bty="n", col="pink", pch=16, "CV Optimum")
	}else{
		h <- hamlet::hmap(as.matrix(mat), Rowv=NA, Colv=NA, col=palet, border=border)
	}
	title(xlab=xlab, ylab=ylab, main=main)
	# Add color key
	hamlet::hmap.key(h)
}
