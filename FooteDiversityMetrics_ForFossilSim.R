

## FUNCTION TO GET PRESENCE/ABSENCE MATRIX
# arguments:
# lineageRanges = lineageRange data from get.LineageRange Function
# age = age parameter from running your tree

# Modified for Silvestro birthdeath_simulator output
get.paMatrix <- function(lineageRanges, age){
	dataTable <- lineageRanges

	# Fill in diversity matrix with 0/1 
	diversityMatrix <- matrix(0, nrow = nrow(dataTable), ncol = age)	
	for(k in 1:nrow(dataTable)){
		ctr <- 1
		for(i in age:1){
			if(dataTable$te[k] <= i & dataTable$ts[k] > (i - 1) ){
				diversityMatrix[k, ctr] <- 1
			}
			ctr <- ctr + 1	
		}
	}	
	return(diversityMatrix)
	
}




#############
## FUNCTION TO IMPOSE PRESERVATION SCENARIO
#paMatrix <- paMatrixList[[1]]

# Constant sampling prob
#sampProb <- rep(0.5, age)
# Linear increasing to present sampling function
#sampProb <- seq(from = 0.1, to = 0.5, by = ((0.5-0.1)/(age-1)))

# Constant
get.sampMatrix.constant <- function(paMatrix, constantP){
	for(i in 1:ncol(paMatrix)){
		index <- which(paMatrix[,i] == 1)
		keep <- sample(1:length(index), constantP*length(index), replace = FALSE)
		paMatrix[index[-keep],i] <- 0
	}
	return(paMatrix)
}

#Linear Increasing
get.sampMatrix.linearInc <- function(paMatrix, startP, endP){
	incP <- seq(from = startP, to = endP, by = ((endP-startP)/(ncol(paMatrix)-1)))
	for(i in 1:ncol(paMatrix)){
		index <- which(paMatrix[,i] == 1)
		keep <- sample(1:length(index), incP[i]*length(index), replace = FALSE)
		paMatrix[index[-keep],i] <- 0
	}
	return(paMatrix)	
	
}

#Shifting - function has to be defined by user
#shiftingProb <- c(0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.5)
get.sampMatrix.shifting <- function(paMatrix, shiftingP, age){
	for(i in 1:ncol(paMatrix)){
		index <- which(paMatrix[,i] == 1)
		keep <- sample(1:length(index), shiftingP[i]*length(index), replace = FALSE)
		paMatrix[index[-keep],i] <- 0
	}
	return(paMatrix)		
}

#sum(paMatrixList[[1]])
#prob.30 <- get.presMatrix.constant(paMatrixList[[1]], 0.3)
#sum(prob.30)/sum(paMatrixList[[1]])
#prob.50 <- get.presMatrix.constant(paMatrixList[[1]], 0.5)
#sum(prob.50)/sum(paMatrixList[[1]])
#prob.1inc5 <- get.presMatrix.linearInc(paMatrixList[[1]], 0.1, 0.5)
#sum(prob.1inc5)/sum(paMatrixList[[1]])


###### FUNCTION TO RANGE THROUGH PA_MATRIX AFTER SAMPLING (AND REMOVE EMPTY ROWS?)
#tmp <- sampMatrixList[[1]]
#test <- apply(tmp, 1, function(x) paste(x, collapse=''))
#sapply(test, function(x) grepl("10+1", x))
#require(gsubfn)

#sampMatrix <- sampMatrixList[[1]]

get.paSampledMatrix <- function(sampMatrix){
	tmp <- apply(sampMatrix, 1, function(x) paste(x, collapse=''))
	for(i in 1:length(tmp)){
		while(grepl("10+1", tmp[i])){
			tmp[i] <- gsubfn("10+1", function (x) paste(replicate(nchar(x), "1"), collapse = ""), 					tmp[i])
		}
	}
	tmp <- do.call('rbind', strsplit(tmp, split = "")) 
	mode(tmp) <- 'numeric'
	return(tmp)
}

#############
## FUNCTION TO CALCULATE FOOTE DIVERSIFICATION METRICS
# arguments:
# paMatrix = presence/absence matrix from get.paMatrix Function

get.diversityMetrics <- function(paMatrix){
	diversityMatrix <- paMatrix

	# Sum raw diversity per interval
	diversityRaw <- numeric(ncol(diversityMatrix));
	for(i in 1:ncol(diversityMatrix)){
		diversityRaw[i] <- sum(diversityMatrix[,i]);
	}

	# Diversity No Singletons Matrix: presence/absence of each non-singleton species for each time interval
	divNoSMatrix <- diversityMatrix;
	for(k in 1:nrow(divNoSMatrix)){
		if(sum(diversityMatrix[k,]) <= 1 ){
			divNoSMatrix[k,] <- 0;
		}
	}

	# Sum no singleton diversity per interval
	diversityNoS <- numeric(ncol(diversityMatrix));
	for(i in 1:ncol(divNoSMatrix)){
		diversityNoS[i] <- sum(divNoSMatrix[,i]); 
	}

	# Nbt Matrix: Range through taxa matrix
	NbtMatrix <- matrix(0, nrow = nrow(diversityMatrix), ncol = ncol(diversityMatrix));
	for(j in 1:nrow(divNoSMatrix)){
		for(k in 2:(ncol(divNoSMatrix)-1)){
			if(divNoSMatrix[j,k] == 1 & divNoSMatrix[j,(k-1)] == 1 & divNoSMatrix[j,(k+1)] == 1){
				NbtMatrix[j,k] <- 1;
			}
		}
	}

	# Count range-through taxa per interval
	Nbt <- numeric(ncol(diversityMatrix));
	for(i in 1:ncol(NbtMatrix)){
		Nbt[i] <- sum(NbtMatrix[,i]); 
	}

	# Nft Matrix: First top taxa matrix (Origination events)
	NftMatrix <- matrix(0, nrow = nrow(diversityMatrix), ncol = ncol(diversityMatrix));
	for(j in 1:nrow(divNoSMatrix)){
		for(k in 2:ncol(divNoSMatrix)){
			if(divNoSMatrix[j,k] == 1 & divNoSMatrix[j,(k-1)] == 0){
				NftMatrix[j,k] <- 1;
			}
		}
	}

	# Sum Nft taxa per interval
	Nft <- numeric(ncol(diversityMatrix));
	for(i in 1:ncol(NftMatrix)){
		Nft[i] <- sum(NftMatrix[,i]); 
	}

	# Nbl Matrix: Last bottom taxa matrix (Extinction events)
	NblMatrix <- matrix(0, nrow = nrow(diversityMatrix), ncol = ncol(diversityMatrix));
	for(j in 1:nrow(divNoSMatrix)){
		for(k in 1:(ncol(divNoSMatrix)-1)){
			if(divNoSMatrix[j,k] == 1 & divNoSMatrix[j,(k+1)] == 0){
				NblMatrix[j,k] <- 1;
			}
		}
	}
		
	# Sum Nbl taxa per interval
	Nbl <- numeric(ncol(diversityMatrix));
	for(i in 1:ncol(NblMatrix)){
		Nbl[i] <- sum(NblMatrix[,i]); 
	}


	#################################
	# Rates Calculations:
	# p_i -> origination rate
	# q_i -> extinction rate
	# d_i -> diversification rate
	# t_i -> turnover rate

	p_i <- numeric(ncol(diversityMatrix));
	q_i <- numeric(ncol(diversityMatrix));
	d_i <- numeric(ncol(diversityMatrix));
	t_i <- numeric(ncol(diversityMatrix));
	for(i in 1:ncol(diversityMatrix)){
		p_i[i] <- log((Nbt[i]+Nft[i])/Nbt[i]);
		q_i[i] <- log((Nbt[i] + Nbl[i])/Nbt[i]);
		d_i[i] <- (Nft[i] - Nbl[i])/Nbt[i];
		t_i[i] <- (Nft[i] + Nbl[i])/Nbt[i];
	}

	# Create and save a new dataframe with interval boundaries and diversity metrics
	divMetricsDF <- data.frame("Interval" = 0:(ncol(diversityMatrix)-1), "DiversityRaw" = diversityRaw, "DivNoSingletons" = diversityNoS, "Nbt" = Nbt, 	"Nft_O" = Nft, "Nbl_E" = Nbl, "OrigRate" = p_i, "ExtRate" = q_i, "DivRate" = d_i, "TurnoverRate" = t_i, stringsAsFactors = FALSE);
	return(divMetricsDF)
}


#############
## FUNCTION TO GET MAX, MIN AND MEAN VALUES FROM A LIST OF DIVERSITY METRICS

## Functions to get max and min of simulation values:
#table = rbind table of metric from each simulation (each list element)
# 2 in apply function refers to the margin --> apply function to columns of table
colMax <- function(table) {
	if (any(!is.finite(table))) {
		table[which(!is.finite(table), arr.ind = TRUE)] <- NA
	}
	ind <- apply(table, 2, function(x) all(is.na(x)))
	res <- apply(table[,which(ind == FALSE)], 2, max, na.rm = TRUE)
	if (length(ind) > 0) {
		vec <- vector(length = ncol(table))
		vec[which(ind == TRUE)] <- NA
		vec[which(ind == FALSE)] <- res
		res <- vec
	}
	return(res)
}


colMin <- function(table) {
	if (any(!is.finite(table))) {
		table[which(!is.finite(table), arr.ind = TRUE)] <- NA
	}
	ind <- apply(table, 2, function(x) all(is.na(x)))
	res <- apply(table[,which(ind == FALSE)], 2, min, na.rm = TRUE)
	if (length(ind) > 0) {
		vec <- vector(length = ncol(table))
		vec[which(ind == TRUE)] <- NA
		vec[which(ind == FALSE)] <- res
		res <- vec
	}	
	return(res)
}

colMeans2 <- function(table) {
	if (any(!is.finite(table))) {
		table[which(!is.finite(table), arr.ind = TRUE)] <- NA
	}
	ind <- apply(table, 2, function(x) all(is.na(x)))
	res <- apply(table[,which(ind == FALSE)], 2, mean, na.rm = TRUE)
	if (length(ind) > 0) {
		vec <- vector(length = ncol(table))
		vec[which(ind == TRUE)] <- NA
		vec[which(ind == FALSE)] <- res
		res <- vec
	}
	return(res)
}



## Function to get full stats table
get.metricStats <- function(metricList){
	
	#paDiv.mean <- colMeans2(do.call('rbind', lapply(metricList, function(x) x$DiversityRaw)))
	#paDiv.min <- colMin(do.call('rbind', lapply(metricList, function(x) x$DiversityRaw)))
	#paDiv.max <- colMax(do.call('rbind', lapply(metricList, function(x) x$DiversityRaw)))

	#nsDiv.mean <- colMeans2(do.call('rbind', lapply(metricList, function(x) x$DivNoSingletons)))
	#nsDiv.min <- colMin(do.call('rbind', lapply(metricList, function(x) x$DivNoSingletons)))
	#nsDiv.max <- colMax(do.call('rbind', lapply(metricList, function(x) x$DivNoSingletons)))

	origRate.mean <- colMeans2(do.call('rbind', lapply(metricList, function(x) x$OrigRate)))
	origRate.min <- colMin(do.call('rbind', lapply(metricList, function(x) x$OrigRate)))
	origRate.max <- colMax(do.call('rbind', lapply(metricList, function(x) x$OrigRate)))

	extRate.mean <- colMeans2(do.call('rbind', lapply(metricList, function(x) x$ExtRate)))
	extRate.min <- colMin(do.call('rbind', lapply(metricList, function(x) x$ExtRate)))
	extRate.max <- colMax(do.call('rbind', lapply(metricList, function(x) x$ExtRate)))

	divRate.mean <- colMeans2(do.call('rbind', lapply(metricList, function(x) x$DivRate)))
	divRate.min <- colMin(do.call('rbind', lapply(metricList, function(x) x$DivRate)))
	divRate.max <- colMax(do.call('rbind', lapply(metricList, function(x) x$DivRate)))


	#toRate.mean <- colMeans2(do.call('rbind', lapply(metricList, function(x) x$TurnoverRate)))
	#toRate.min <- colMin(do.call('rbind', lapply(metricList, function(x) x$TurnoverRate)))
	#toRate.max <- colMax(do.call('rbind', lapply(metricList, function(x) x$TurnoverRate)))

	#df <- data.frame(paDiv.mean, paDiv.min, paDiv.max, nsDiv.mean, nsDiv.min, nsDiv.max, origRate.mean, origRate.min, origRate.max, extRate.mean, extRate.min, extRate.max, divRate.mean, divRate.min, divRate.max, toRate.mean, toRate.min, toRate.max)
	df <- data.frame(origRate.mean, origRate.min, origRate.max, extRate.mean, extRate.min, extRate.max, divRate.mean, divRate.min, divRate.max)
	return(df)
}



#####################
######################
# FUNCTION TO PLOT TREE, LINEAGE RANGES, DIVERSITY, AND FOOTE DIVERSITY METRICS FOR ONE TREE
# arguments:
# tree = simulated tree
# lineageRange = lineageRange output from getLineageRange Function
# divMetrics = diversification metrics output from diversityMetrics Function 

# LTT from phytools
#require(phytools)

plot.treeToRates <- function(tree, lineageRanges, divMetrics){
	
	dataframe <- divMetrics
	
	#Plot parameters
	quartz.options(height = 12, width = 12, dpi =72)
	par(mfcol = c(3,2))
	par(oma = c(1,1,1,3), mar = c(3,3,0,3))

	# Plot phylogeny and durations
	plot.window(xlim = c(0,max(lineageRanges$lineageEnd)), ylim = c(0,length(tree$tip.label)))
	plot(tree)
	axis(side = 1)
	
	plot.new()
	plot.window(xlim = c(0,max(lineageRanges$lineageEnd+0.5)), ylim = c(0,length(tree$tip.label)))
	segments(x0 = lineageRanges$lineageStart, y0 = 1:length(tree$tip.label), x1 = lineageRanges$lineageEnd)
	#text(rownames(lineageRanges), x = 6.2, y = 1:length(z[[1]]$tip.label), )
	axis(side = 1)
	axis(side = 2)

	# Plot diversity
	col1 <- "#66a61e";
	col2 <-  "#e7298a"
	lty1 <- 1;
	lty2 <- 2;
	plot.new();
	plot.window(xlim = c(0,max(lineageRanges$lineageEnd+0.5)), ylim = c(0, max(dataframe$DiversityRaw)));

	lines(x = dataframe$Interval+0.5, y = dataframe$DiversityRaw, lwd=3, col=col1);
	lines(x = dataframe$Interval+0.5, y = dataframe$DivNoSingletons, lwd=3, lty = lty2, col=col2);

	# assign axes, labels and legend
	axis(side=1)
	axis(side=2, at=seq(0, max(dataframe$DiversityRaw), by=2), cex.axis=1.2);
	#mtext(side=1, text="Millions of years before present", line=3, cex=1);
	mtext(side=2, text="Diversity", line=3, cex=1);
	#mtext(side=3, text = title, line=1, cex=1.2);
	legend("topleft", inset = 0.05, c("PA Diversity", "No Singleton Diversity"), cex = 1, col = c(col1, col2), lty = c(lty1,lty2), lwd = 2.5, bty = "n");


	# Plot sp/ext rates
	col1 <- "#d95f02";
	col2 <- "#7570b3";
	lty1 <- 1;
	lty2 <- 1;

	plot.new();
	plot.window(xlim = c(0,max(lineageRanges$lineageEnd+0.5)), ylim = c(0, max(c(dataframe$OrigRate[which(is.finite(dataframe$OrigRate))], dataframe$ExtRate[which(is.finite(dataframe$ExtRate))]), na.rm = T)));

	## plot lines
	# lines for origination rate through time
	lines(x = dataframe$Interval+0.5, y = dataframe$OrigRate, lwd=3, col=col1);
	# lines for confidence intervals
	#for(i in 1:length(timeKeep)){
		#lines(x = c(timeKeep[i],timeKeep[i],timeKeep[i]), y = c(p_CIlow[i], metricsKeep$OrigRate[i], p_CIhigh[i]), lwd = 1, col = colOCI);
	#}
	# points for significant intervals
	#points(x = sigPintervals, y = origSig, pch = 8, col = col1, cex = 1.5);

	# lines for extinction rate through time
	lines(x = dataframe$Interval+0.5, y = dataframe$ExtRate, lwd=3, lty = lty2, col=col2);
	# lines for confidence intervals
	#for(i in 1:length(timeKeep)){
	#	lines(x = c(timeKeep[i],timeKeep[i],timeKeep[i]), y = c(q_CIlow[i], metricsKeep$ExtRate[i], q_CIhigh[i]), lwd = 1, col = colECI);
	#}
	# points for significant intervals
	#points(x = sigQintervals, y = extSig, pch = 8, col = col2, cex = 1.5);

	axis(side = 1)
	axis(side = 2)
	# assign axes, labels and legend
	#axis(side=1, at=seq(xmax, xmin, by = -2), cex.axis=1.2);
	#axis(side=2, at=seq(0, ymax, by=0.5), cex.axis=1.2);
	#mtext(side=1, text="Millions of years before present", line=3, cex=1);
	#mtext(side=2, text="Per Capita Rate", line=3, cex=1);
	#mtext(side=3, text = title, line=1, cex=1.2)
	#legend("topleft", inset = 0.05, c("Origination Rate", "Extinction Rate", "Significant Rate"), cex = 1, col = c(col1, col2, "black"), lty = c(lty1,lty2, NA), pch = c(NA,NA,8), lwd = 2.5);
	legend("topleft", inset = 0.05, c("Origination Rate", "Extinction Rate"), cex = 1, col = c(col1, col2), lty = c(lty1,lty2), lwd = 2.5, bty = "n")


	## Plot div/to rates
	col1 <- "#1b9e77";
	col2 <- "#e6ab02";
	lty1 <- 1;
	lty2 <- 1;

	# Plot parameters
	plot.new();
	plot.window(xlim = c(0,max(lineageRanges$lineageEnd+0.5)), ylim = c(min(c(dataframe$DivRate[which(is.finite(dataframe$DivRate))], dataframe$TurnoverRate[which(is.finite(dataframe$TurnoverRate))]), na.rm = T), max(c(dataframe$DivRate[which(is.finite(dataframe$DivRate))], dataframe$TurnoverRate[which(is.finite(dataframe$Turnover))]), na.rm = T)));

	abline(h = 0, col = "gray", lty = 2, lwd = 0.7)
	## plot lines
	# lines for origination rate through time
	lines(x = dataframe$Interval+0.5, y = dataframe$DivRate, lwd=3, col=col1);
	# lines for confidence intervals
	#for(i in 1:length(timeKeep)){
		#lines(x = c(timeKeep[i],timeKeep[i],timeKeep[i]), y = c(p_CIlow[i], metricsKeep$OrigRate[i], p_CIhigh[i]), lwd = 1, col = colOCI);
	#}
	# points for significant intervals
	#points(x = sigPintervals, y = origSig, pch = 8, col = col1, cex = 1.5);

	# lines for extinction rate through time
	lines(x = dataframe$Interval+0.5, y = dataframe$TurnoverRate, lwd=3, lty = lty2, col=col2);
	# lines for confidence intervals
	#for(i in 1:length(timeKeep)){
	#	lines(x = c(timeKeep[i],timeKeep[i],timeKeep[i]), y = c(q_CIlow[i], metricsKeep$ExtRate[i], q_CIhigh[i]), lwd = 1, col = colECI);
	#}
	# points for significant intervals
	#points(x = sigQintervals, y = extSig, pch = 8, col = col2, cex = 1.5);

	axis(side = 1)
	axis(side = 2)
	# assign axes, labels and legend
	#axis(side=1, at=seq(xmax, xmin, by = -2), cex.axis=1.2);
	#axis(side=2, at=seq(0, ymax, by=0.5), cex.axis=1.2);
	#mtext(side=1, text="Millions of years before present", line=3, cex=1);
	#mtext(side=2, text="Per Capita Rate", line=3, cex=1);
	#mtext(side=3, text = title, line=1, cex=1.2)
	#legend("topleft", inset = 0.05, c("Diversification Rate", "Turnover Rate", "Significant Rate"), cex = 1, col = c(col1, col2, "black"), lty = c(lty1,lty2, NA), pch = c(NA,NA,8), lwd = 2.5);
	legend("topleft", inset = 0.05, c("Diversification Rate", "Turnover Rate"), cex = 1, col = c(col1, col2), lty = c(lty1,lty2), lwd = 2.5, bty = "n");
	
	ltt(tree, log.lineages = FALSE, gamma = FALSE)
	legend("topleft", inset = 0.05, c("LTT plot"), cex = 2, bty = "n")
}
	


#####################
######################
# FUNCTION TO PLOT DIVERSITY AND FOOTE DIVERSITY METRICS FOR A TREE LIST -> DIVERSITY METRICS FOR SIMULATIONS
# arguments:
# tree = simulated tree
# lineageRange = lineageRange output from getLineageRange Function
# divMetrics = diversification metrics output from diversityMetrics Function 

plot.divMetrics <- function(metricStats, lwd = 3, alpha = 0.5, cex.axis = 1.5, cex = 1){	

#Plot parameters
	quartz.options(height = 12, width = 6, dpi =72)
	par(mfcol = c(3,1))
	par(oma = c(3,2,1,1), mar = c(3,3,1,1))


## Plot diversity
	col1 <- "#66a61e";
	col2 <-  "#e7298a"
	col1Range <- as.vector(col2rgb(col1))/255
	col2Range <- as.vector(col2rgb(col2))/255
	plot.new();
	plot.window(xlim = c(0,nrow(metricStats)), ylim = c(0, max(metricStats$paDiv.max)));
	#pa diversity
	polygon(x = c(1:nrow(metricStats), nrow(metricStats):1), y = c(metricStats$paDiv.min, rev(metricStats$paDiv.max)), col = rgb(red = col1Range[1], green = col1Range[2], blue = col1Range[3], alpha = alpha), border = NA)
	lines(x = 1:nrow(metricStats), y = metricStats$paDiv.mean, col = col1, lwd = lwd)
	#no singleton diversity
	polygon(x = c(1:nrow(metricStats), nrow(metricStats):1), y = c(metricStats$nsDiv.min, rev(metricStats$nsDiv.max)), col = rgb(red = col2Range[1], green = col2Range[2], blue = col2Range[3], alpha = alpha), border = NA)
	lines(x = 1:nrow(metricStats), y = metricStats$nsDiv.mean, col = col2, lwd = lwd)
	# assign axes, labels and legend
	axis(side=1, cex.axis = cex.axis)
	axis(side = 2, cex.axis = cex.axis)
	#axis(side=2, at=seq(0, max(metricStats$paDiv.max), by=2), cex.axis=cex);
	#mtext(side=1, text="Interval", line=3, cex=1);
	mtext(side=2, text="Diversity", line=3, cex=cex);
	legend("topleft", inset = 0.05, c("PA Diversity", "No Singleton Diversity"), cex = 1, col = c(col1, col2), lty = c(1,1), lwd = 2.5, bty = "n");
	

## Plot orig/ext rates
	col1 <- "#d95f02";
	col2 <- "#7570b3";
	col1Range <- as.vector(col2rgb(col1))/255
	col2Range <- as.vector(col2rgb(col2))/255
	plot.new();
	plot.window(xlim = c(0,nrow(metricStats)), ylim = c(min(c(metricStats$origRate.min, metricStats$extRate.min), na.rm = T), max(c(metricStats$origRate.max, metricStats$extRate.max), na.rm = T)));
	# Origination Rates
	x <- c(1:nrow(metricStats), nrow(metricStats):1)
	y <- c(metricStats$origRate.min, rev(metricStats$origRate.max))
	x <- x[which(!is.na(y))]
	y <- y[which(!is.na(y))]
	polygon(x, y, col = rgb(red = col1Range[1], green = col1Range[2], blue = col1Range[3], alpha = alpha), border = NA)
	lines(x = 1:nrow(metricStats), y = metricStats$origRate.mean, col = col1, lwd = lwd)
	# Extinction Rates
	x <- c(1:nrow(metricStats), nrow(metricStats):1)
	y <- c(metricStats$extRate.min, rev(metricStats$extRate.max))
	x <- x[which(!is.na(y))]
	y <- y[which(!is.na(y))]
	polygon(x , y , col = rgb(red = col2Range[1], green = col2Range[2], blue = col2Range[3], alpha = alpha), border = NA)
	lines(x = 1:nrow(metricStats), y = metricStats$extRate.mean, col = col2, lwd = lwd)
	# assign axes, labels and legend
	axis(side=1, cex.axis = cex.axis)
	axis(side=2, cex.axis = cex.axis)
	#mtext(side=1, text="Interval", line=3, cex=cex);
	mtext(side=2, text="Per Capita Rate", line=3, cex=cex);
	legend("topright", inset = 0.05, c("Origination Rate", "Extinction Rate"), cex = 1, col = c(col1, col2), lty = c(1,1), lwd = 2.5, bty = "n")
	
## Plot div/to rates
	col1 <- "#1b9e77";
	col2 <- "#e6ab02";
	col1Range <- as.vector(col2rgb(col1))/255
	col2Range <- as.vector(col2rgb(col2))/255
	plot.new();
	plot.window(xlim = c(0,nrow(metricStats)), ylim = c(min(c(metricStats$divRate.min, metricStats$toRate.min), na.rm = T), max(c(metricStats$divRate.max, metricStats$toRate.max), na.rm = T)));
	# Diversification Rates
	x <- c(1:nrow(metricStats), nrow(metricStats):1)
	y <- c(metricStats$divRate.min, rev(metricStats$divRate.max))
	x <- x[which(!is.na(y))]
	y <- y[which(!is.na(y))]
	polygon(x , y, col = rgb(red = col1Range[1], green = col1Range[2], blue = col1Range[3], alpha = alpha), border = NA)
	lines(x = 1:nrow(metricStats), y = metricStats$divRate.mean, col = col1, lwd = lwd)
	# Turnover Rates
	x <- c(1:nrow(metricStats), nrow(metricStats):1)
	y <- c(metricStats$toRate.min, rev(metricStats$toRate.max))
	x <- x[which(!is.na(y))]
	y <- y[which(!is.na(y))]
	polygon(x , y, col = rgb(red = col2Range[1], green = col2Range[2], blue = col2Range[3], alpha = alpha), border = NA)
	lines(x = 1:nrow(metricStats), y = metricStats$toRate.mean, col = col2, lwd = lwd)
	# assign axes, labels and legend
	axis(side=1, cex.axis = cex.axis)
	axis(side=2, cex.axis = cex.axis)
	mtext(side=1, text="Interval", line=3, cex=cex);
	mtext(side=2, text="Per Capita Rate", line=3, cex=cex);
	legend("topright", inset = 0.05, c("Diversification Rate", "Turnover Rate"), cex = 1, col = c(col1, col2), lty = c(1,1), lwd = 2.5, bty = "n")
}
