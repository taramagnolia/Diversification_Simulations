########  Read in matrix list of occurrence data for preservation scenarios
paMatrixList <- readRDS("paMatrixList_constantLam_random0.rds")

######  30R: Constant preservation rate ###### 
	constantP <- 0.3
	sampMatrixList <- vector("list", length = length(paMatrixList))
	sampVec <- numeric(length(paMatrixList))
	for(i in 1:length(paMatrixList)){
		sampMatrixList[[i]] <- get.sampMatrix.constant(paMatrixList[[i]], constantP)
		sampVec[i] <- sum(sampMatrixList[[i]])/sum(paMatrixList[[i]])
	}
###### 30R

###### IncR: Linear increasing preservation rate ###### 
	startP <- 0.1
	endP <- 0.5
	sampMatrixList <- vector("list", length = length(paMatrixList))
	sampVec <- numeric(length(paMatrixList))
	for(i in 1:length(paMatrixList)){
		sampMatrixList[[i]] <- get.sampMatrix.linearInc(paMatrixList[[i]], startP, endP)
		sampVec[i] <- sum(sampMatrixList[[i]])/sum(paMatrixList[[i]])
	}
###### IncR

###### PulseR: Shifting preservation rate ######
	shiftingP <- c(rep(0.1, times =16), 0.35, 0.6, 0.527, 0.463, 0.406, 0.357, 0.313, rep(0.3, times = 17))
	sampMatrixList <- vector("list", length = length(paMatrixList))
	sampVec <- numeric(length(paMatrixList))
	for(i in 1:length(paMatrixList)){
		sampMatrixList[[i]] <- get.sampMatrix.shifting(paMatrixList[[i]], shiftingP)
		sampVec[i] <- sum(sampMatrixList[[i]])/sum(paMatrixList[[i]])
	}
###### PulseR

###### FreqR: Fossil Occurrences Preservation Rate ######
	presFossOcc <- read.csv("~/Dropbox/Research_Div Models/GeneraFADLADforTom/PresRateFossilOcc.csv", header = T, stringsAsFactor = F)
	noNA <- presFossOcc$R_red50[!is.na(presFossOcc$R_red50)]	
	sampMatrixList <- vector("list", length = length(paMatrixList))
	sampVec <- numeric(length(paMatrixList))
	for(i in 1:length(paMatrixList)){
		shiftingP <- c(sample(presFossOcc$R_red50, size = 5, replace = T), presFossOcc$R_red50, sample(presFossOcc$R_red50, size = 5, replace = T))
		for(k in 1:length(shiftingP)){
			if(is.na(shiftingP[k])){
				shiftingP[k] <- sample(noNA, size = 1)
			}
		}
		shiftingP <- rev(shiftingP)
		sampMatrixList[[i]] <- get.sampMatrix.shifting(paMatrixList[[i]], shiftingP)
		sampVec[i] <- sum(sampMatrixList[[i]])/sum(paMatrixList[[i]])
	}
###### FreqR

###### StratR: Macrostrat Preservation Rate ######
	macrostrat <- read.csv("~/Dropbox/Research_Div Models/macrostrat_rateParameter.csv", header = T, stringsAsFactor = F)
	length(macrostrat$prop_red50)
	shiftingP <- c(rep(macrostrat$prop_red50[1], times = 5), macrostrat$prop_red50, rep(macrostrat$prop_red50[30], times = 5))
	length(shiftingP)
	shiftingP <- rev(shiftingP)
	sampMatrixList <- vector("list", length = length(paMatrixList))
	sampVec <- numeric(length(paMatrixList))
	for(i in 1:length(paMatrixList)){
		sampMatrixList[[i]] <- get.sampMatrix.shifting(paMatrixList[[i]], shiftingP)
		sampVec[i] <- sum(sampMatrixList[[i]])/sum(paMatrixList[[i]])
	}
###### StratR


