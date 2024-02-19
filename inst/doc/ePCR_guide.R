## -----------------------------------------------------------------------------

library(ePCR)
# Kernel density simulated patients from Turku University Hospital (TYKS)
# Data consists of TEXT cohort (text-search found patients) 
# and MEDI (patients identified using medication and few keywords)
data(TYKSSIMU)
# The following data matrices x and survival responses y become available
head(xTEXTSIMU); head(yTEXTSIMU) 
head(xMEDISIMU); head(yMEDISIMU)

library(survival)


## ----message=FALSE, warning=FALSE---------------------------------------------

testset <- 1:30
# Medication cohort fit
# Leaving out patients into a separate test set using negative indices
psp_medi <- new("PSP", 
 	# Input data matrix x (example data loaded previously)
 	x = xMEDISIMU[-testset,],
 	# Response vector, 'surv'-object
 	y = yMEDISIMU[-testset,"surv"],
 	# Seeds for reproducibility
 	seeds = c(1,2),
 	# If user wishes to run the CV binning multiple times,
 	# this is possible by averaging over them for smoother CV heatmap.
 	cvrepeat = 2,
 	# Using the concordance-index as prediction accuracy in CV
 	score = score.cindex,
 	# Alpha sequence
 	alphaseq = seq(from=0, to=1, length.out=6),
 	# Using glmnet's default nlambda of 100
 	nlambda = 100,
 	# Running the nominal 10-fold cross-validation
 	folds = 10,
 	# x.expand slot is a function that would allow interaction terms
 	# For the sake of the simplicity we will consider identity function
 	x.expand = function(x) { as.matrix(x) }
)


## ----message=FALSE, warning=FALSE---------------------------------------------

# Text run similar to above
# Leaving out patients into a separate test set using negative indices
psp_text <- new("PSP", 
 	x = xTEXTSIMU[-testset,],
 	y = yTEXTSIMU[-testset,"surv"],
 	seeds = c(3,4),
 	cvrepeat = 2,
 	score = score.cindex,
 	alphaseq = seq(from=0, to=1, length.out=6),
 	nlambda = 100,
 	folds = 10,
 	x.expand = function(x) { as.matrix(x) }
)


## -----------------------------------------------------------------------------

# Taking a look on the show-method for PSP:
psp_medi


## ----fig1, fig.height = 7, fig.width = 7, fig.align = "center"----------------

# Plot the CV-surface of the fitted PSP:
plot(psp_medi, 
 	# Showing only every 10th row and column name (propagated to heatcv-function)
 	by.rownames=10, by.colnames=10, 
 	# Adjust main title and tilt the bias of the color key legend (see ?heatcv)
 	main="C-index CV for psp_medi", bias=0.2)


## ----fig2, fig.height = 7, fig.width = 7, fig.align = "center"----------------

plot(psp_text, 
 	# Showing only every 10th row and column name (propagated to heatcv-function)
 	by.rownames=10, by.colnames=10, 
 	# Adjust main title and tilt the bias of the color key legend (see ?heatcv)
 	main="C-index CV for psp_text", bias=0.2)	


## -----------------------------------------------------------------------------
psp_medi@optimum
psp_text@optimum
slotNames(psp_medi)

## -----------------------------------------------------------------------------

pep_tyks <- new("PEP",
 	# The main input is the list of PSP objects
 	PSPs = list(psp_medi, psp_text)
)
# These PSPs were constructed using the example code above.
pep_tyks


## -----------------------------------------------------------------------------

# Conduct naive test set evaluation
xtest <- rbind(xMEDISIMU[testset,], xTEXTSIMU[testset,])
ytest <- rbind(yMEDISIMU[testset,], yTEXTSIMU[testset,])
# Perform survival prediction based on the PEP-ensemble we've created
xpred <- predict(pep_tyks, newx=as.matrix(xtest), type="ensemble")
# Construct a survival object using the Surv-class
ytrue <- Surv(time = ytest[,"surv"][,"time"], event = ytest[,"surv"][,"status"])
# Test c-index between our constructed ensemble prediction and true response
tyksscore <- score.cindex(pred = xpred, real = ytrue)
print(paste("TYKS example c-index:", round(tyksscore, 4)))


## -----------------------------------------------------------------------------
data(ePCRmodels)
class(DREAM)
class(TYKS)

## -----------------------------------------------------------------------------

# Create a DREAM-matching data input matrix from our xtest and the full data matrix
xtemp <- conforminput(DREAM, xtest)
# Predict survival for our hospital registry example dataset 
dreampred <- predict(DREAM, 
 	# Providing full new data and average prediction over the ensemble members
 	newx=xtemp, type="ensemble",
 	# Defining that we don't want any further data matrix feature extraction
 	# The call to conforminput above already formatted the input data
 	x.expand = as.matrix
)


## -----------------------------------------------------------------------------
# Test c-index between the DREAM ensemble prediction and TYKS true response
dreamscore <- score.cindex(pred = dreampred, real = ytrue)
print(paste("DREAM example c-index:", round(dreamscore, 4)))

## -----------------------------------------------------------------------------
sessionInfo()

