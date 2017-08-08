## Calculating Three-Timer (Alroy 2008, 2010) diversification metrics for Basin and Range rodents
# TM Smiley
# July 2017

rm(list = ls())

setwd("~/gitCode/Diversification_Simulations")
source('Liow2014functions_corrected2017.R')
master.data <- read.csv("Miomap_MasterData_faunalANDlocality_withRegion_speciesOnly.csv", header = T, stringsAsFactors = F)

# BR
timeVec <- 35:5 
data <- master.data
data <- data[data$Region == 3,]
data <- data[which(data$Taxon.Order == "Rodentia"),]


species <- unique(data$Taxon.Name)
diversityMatrixMaxAge <- matrix(0, nrow = length(species), ncol = length(timeVec));
diversityMatrixMinAge <- matrix(0, nrow = length(species), ncol = length(timeVec));
diversityMatrixRandAge <- matrix(0, nrow = length(species), ncol = length(timeVec));
colnames(diversityMatrixMaxAge) <- 35:5
colnames(diversityMatrixMinAge) <- 35:5
colnames(diversityMatrixRandAge) <- 35:5

for(i in 1:length(species)){
	tempDS <- data[which(data$Taxon.Name == species[i]),]
	maxVec <- ceiling(tempDS$Maximum.Age)
	minVec <- ceiling(tempDS$Minimum.Age)
	randVec <- NULL
	for(m in 1:nrow(tempDS)){
		randVec <- c(randVec, sample(c(ceiling(tempDS$Maximum.Age[m]): ceiling(tempDS$Minimum.Age[m])), 1))

	}
	diversityMatrixMaxAge[i,] <- as.numeric(35:5 %in% as.character(maxVec))
	diversityMatrixMinAge[i,] <- as.numeric(35:5 %in% as.character(minVec))
	diversityMatrixRandAge[i,] <- as.numeric(35:5 %in% as.character(randVec))
}

#Alroy metrics
alroyMax <- alroyrates(diversityMatrixMaxAge, timeVec)
alroyMin <- alroyrates(diversityMatrixMinAge, timeVec)
alroyRand <- alroyrates(diversityMatrixRandAge, timeVec)	

