# ePCR
CRAN R package - ePCR: Ensemble Penalized Cox Regression for Survival Prediction

## Description

The top-performing ensemble-based Penalized Cox Regression (ePCR) framework developed during the 
DREAM 9.5 mCRPC Prostate Cancer Challenge <https://www.synapse.org/ProstateCancerChallenge> presented 
in Guinney J, Wang T, Laajala TD, et al. (2017) <doi:10.1016/S1470-2045(16)30560-5> is provided here-in, 
together with the corresponding follow-up work. While initially aimed at modeling the most advanced stage 
of prostate cancer, metastatic Castration-Resistant Prostate Cancer (mCRPC), the modeling framework has 
subsequently been extended to cover also the non-metastatic form of advanced prostate cancer (CRPC). Readily 
fitted ensemble-based model S4-objects are provided, and a simulated example dataset based on a real-life 
cohort is provided from the Turku University Hospital, to illustrate the use of the package. Functionality 
of the ePCR methodology relies on constructing ensembles of strata in patient cohorts and averaging over them, 
with each ensemble member consisting of a highly optimized penalized/regularized Cox regression model. 
Various cross-validation and other modeling schema are provided for constructing novel model objects.

## Citation

Methodology:

Guinney J*, Wang T*, Laajala TD*, Winner KK, Bare JC, Neto EC, Khan SA, Peddinti G, Airola A, Pahikkala T, Mirtti T, Yu T, Bot BM, Shen L, Abdallah K, Norman T, Friend S, Stolovitzky G, Soule H, Sweeney CJ, Ryan CJ, Scher HI, Sartor O, Xie Y, Aittokallio T, Zhou FL, Costello JC; Prostate Cancer Challenge DREAM Community. _Prediction of overall survival for patients with metastatic castration-resistant prostate cancer: development of a prognostic model through a crowdsourced challenge with open clinical trial data._ Lancet Oncol. 2017 Jan;18(1):132-142. doi: 10.1016/S1470-2045(16)30560-5. Epub 2016 Nov 16. PMID: 27864015; PMCID: PMC5217180.

ePCR R-package:

Laajala TD, Murtojärvi M, Virkki A, Aittokallio T. _ePCR: an R-package for survival and time-to-event prediction in advanced prostate cancer, applied to real-world patient cohorts._ Bioinformatics. 2018 Nov 15;34(22):3957-3959. doi: 10.1093/bioinformatics/bty477. PMID: 29912284; PMCID: PMC6223370.
