## Calculating capture-mark-recapture diversification metrics for Basin and Range rodents
# TM Smiley 
# July 2017


rm(list = ls())

setwd("~/gitCode/Diversification_Simulations")
paMatrixList <- readRDS("paMatrixList_toUseForCMR.rds") # presence-absence matrices (min, max, rand age estimate) for North American rodents

# reading in files and processing the data with MARK
names <- c("Max", "Min", "Rand")
for(k in 1:length(paMatrixList)){
  filename <- paste0("diversityMatrix", names[k], "_CMR.inp")
  for(i in 1:nrow(paMatrixList[[k]])){
    paMatrix <- paMatrixList[[k]]
    vec <- paste(c(paMatrix[i,], " 1;"), collapse = "")
    write(vec, filename, append = ifelse(i == 1, F, T))
  }
}

files <- list.files(, pattern = "diversityMatrix")
neo.pradsenList <- list()
neo.pradLamdaList <- list()
for(i in 1:length(files)){
  neo<- convert.inp(files[i])
  neo.pradsenList[[i]] <- process.data(neo, model= "Pradsen")
  neo.pradLamdaList[[i]] <- process.data(neo, model = "Pradlambda")
}

# Model specifications
#time-varying
Phi.t = list(formula=~time)
Gamma.t = list(formula=~time)
p.t = list(formula=~time)
L.t = list(formula=~time)

#run time-varying model ONLY (two different parameterizations to get phi, gamma and lambda)
timeList <- list()
#timeDList <- list()
for(i in 1:length(files)){
  timeList[[i]] = mark(neo.pradsenList[[i]],  model.parameters = list(Phi = Phi.t, p = p.t, Gamma = Gamma.t), delete = T, brief = T)
  #	timeDList[[i]]=mark(neo.pradLamda,  model.parameters=list(Phi=Phi.t, p=p.t, Lambda=L.t))
}

dim(time$results$real)
cbind(1:91, time$results$real)

#### 
# store results into a metric List for plotting
metricList <- list()
for(i in 1:length(files)){
  time <- timeList[[i]]
  #timeD <- timeDList[[i]]
  phi = time$results$real$estimate[1:30]
  p = time$results$real$estimate[31:61]
  gam = time$results$real$estimate[62:91]
  #div=timeD$results$real$estimate[X:X]
  pcs = -log(gam) #instantaenous rate for comparison
  pce = -log(phi) #instantaenousrate for comparison
  pcd = pcs - pce
  df <- data.frame("phi" = phi, "gam" = gam, "ExtRate" = pce, "OrigRate" = pcs, "DivRate" = pcd)
  metricList[[i]] <- df
}
