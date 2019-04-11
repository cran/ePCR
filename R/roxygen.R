####
#
#	This file exists solely for the purpose of roxygen2-documenting the data and package etc prior to wrapping the documentation into .Rd
#	Teemu Daniel Laajala, teelaa@utu.fi
#
####

#' Ensemble Penalized Cox Regression Modeling for Overall Survival and Time-to-Event Prediction in Advanced Prostate Cancer
#'
#' @name ePCR
#' @references Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Bioinformatics. 2018 Jun 15. doi: 10.1093/bioinformatics/bty477.
#' @references Guinney J, Wang T, Laajala TD, et al. Prediction of overall survival for patients with metastatic castration-resistant prostate cancer: development of a prognostic model through a crowdsourced challenge with open clinical trial data. Lancet Oncol 2017; 18: 132-142.
#' @references Laajala TD, Guinney J, Costello JC. Community mining of open clinical trial data. Oncotarget 2017; 8: 81721-81722. doi: 10.18632/oncotarget.20853.
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @docType package
NULL

#' FIMM-UTU DREAM winning implementation of an ensemble of Penalized Cox Regression models for mCPRC research (ePCR)
#'
#' @name DREAM
#' @docType data
#' @references Guinney J, Wang T, Laajala TD, et al. Prediction of overall survival for patients with metastatic castration-resistant prostate cancer: development of a prognostic model through a crowdsourced challenge with open clinical trial data. Lancet Oncol 2017; 18: 132-142.
#' @note Notice that in order to save space, some slots in the S4 object have been set to null.
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @usage data(ePCRmodels)
"DREAM"

#' ePCR model fitted to the Turku University Hospital cohorts (all features)
#'
#' @name TYKS
#' @docType data
#' @references Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Bioinformatics. 2018 Jun 15. doi: 10.1093/bioinformatics/bty477.
#' @note Notice that in order to save space, some slots in the S4 object have been set to null.
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @usage data(ePCRmodels)
"TYKS"

#' ePCR model fitted to the Turku University Hospital cohorts (features derived from text mining only)
#'
#' @name TYKS_reduced
#' @docType data
#' @references Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Bioinformatics. 2018 Jun 15. doi: 10.1093/bioinformatics/bty477.
#' @note Notice that in order to save space, some slots in the S4 object have been set to null.
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}
#' @usage data(ePCRmodels)
"TYKS_reduced"

#' TYKSSIMU - simulated data matrices and survival responses from Turku University Hospital
#'
#' @name TYKSSIMU
#' @docType data
#' @usage data(TYKSSIMU)
#' @references Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to a real-world patient cohort. Bioinformatics. 2018 Jun 15. doi: 10.1093/bioinformatics/bty477.
#' @examples
#' data(TYKSSIMU)
#' head(xTEXTSIMU)
#' head(xMEDISIMU)
#' head(yTEXTSIMU)
#' head(yMEDISIMU)
#' dim(xTEXTSIMU)
#' dim(xMEDISIMU)
#' @author Teemu Daniel Laajala \email{teelaa@@utu.fi}, Mika Murtojärvi \email{mianmu2@hotmail.com}
#' @format xTEXTSIMU: 
#' @format xMEDISIMU:
#' @format yTEXTSIMU: 
#' @format yMEDISIMU: 
NULL

#' xMEDISIMU: Simulated prostate cancer data from Turku University Hospital (data matrix x, Medication-cohort)
#'
#' @rdname TYKSSIMU
"xMEDISIMU"

#' xTEXTSIMU: Simulated prostate cancer data from Turku University Hospital (data matrix x, Text-cohort)
#'
#' @rdname TYKSSIMU
"xTEXTSIMU"

#' yMEDISIMU: Simulated prostate cancer data from Turku University Hospital (survival response y, Medication-cohort)
#'
#' @rdname TYKSSIMU
"yMEDISIMU"

#' yTEXTSIMU: Simulated prostate cancer data from Turku University Hospital (survival response y, Text-cohort)
#'
#' @rdname TYKSSIMU
"yTEXTSIMU"
