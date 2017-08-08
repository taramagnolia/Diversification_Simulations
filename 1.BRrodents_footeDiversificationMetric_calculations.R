## Calculating per-capita (Foote 2000) diversification metrics for Basin and Range rodents
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

## Foote metrics
### (start with 'data' that has been subset by desired region)
# Find unique species
uniqueSpecies <- unique(data$Taxon.Name);
uniqueSpecies <- uniqueSpecies[!is.na(uniqueSpecies)]

# Assign first and last appearance age (and Family) for each unique species
maxAge <- numeric(length(uniqueSpecies));
minAge <- numeric(length(uniqueSpecies));
familyNames <- character(length(uniqueSpecies));
orderNames <- character(length(uniqueSpecies));
for(i in 1:length(uniqueSpecies)){
	record <- grep(pattern = uniqueSpecies[i], x = data$Taxon.Name);
	tempDS <- data[record,];
	maxAge[i] <- max(tempDS$Maximum.Age);
	minAge[i] <- min(tempDS$Minimum.Age);
	familyNames[i] <- tempDS$Taxon.Family[1];
	orderNames[i] <- tempDS$Taxon.Order[1]
}

#Create and save a new dataframe with Family, Species, Min Age, Max Age
brRecords <- data.frame("Order" = orderNames, "Family" = familyNames, "Species" = uniqueSpecies, "MinAge" = minAge, "MaxAge" = maxAge, stringsAsFactors = FALSE)
FLdata <- brRecords

## for summing diversity, raw and no singleton
timeVec <- 35:5

# Diversity Matrix: presence/absence of each species for each time interval
diversityMatrix <- matrix(0, nrow = nrow(FLdata), ncol = length(timeVec));
for(k in 1:nrow(FLdata)){
	ctr <- 1;
	for(i in timeVec){
		if(FLdata$MinAge[k] <= i & FLdata$MaxAge[k] > (i-1)){
			diversityMatrix[k, ctr] <- 1;
		}
	ctr <- ctr +1;
	}
}

colnames(diversityMatrix) <- timeVec
dataFoote <- diversityMatrix

#Foote diversification metrics
FA=function(x){min(which(x==1))}
LA=function(x){max(which(x==1))}
FAD <- apply(dataFoote,1, FA)
LAD <- apply(dataFoote,1, LA)
FA <- timeVec[FAD]
LA <- timeVec[LAD]
	
foote <- footerates(FA, LA, (timeVec))
foote <- as.data.frame(foote)
rownames(foote) <- timeVec[1:30]