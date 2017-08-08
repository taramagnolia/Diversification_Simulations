# Code for a basic Pradel model
#30 05 2013 Lee Hsiang Liow

### Modified from Liow and Finarelli 2014 (TMS) ###


# MARK, written for Windows, must be downloaded before RMark can be used.Find MARK at http://warnercnr.colostate.edu/~gwhite/mark/mark.htm 
library(RMark) 

paMatrixList <- readRDS("paMatrixList_toUseForCMR.rds") # presence-absence matrices (min, max, rand age estimate) for North American rodents

# reading in files, generating input file for processing in MARK
names <- c("Max", "Min", "Rand")
for(k in 1:length(paMatrixList)){
  filename <- paste0("diversityMatrix", names[k], "_CMR.inp")
  for(i in 1:nrow(paMatrixList[[k]])){
    paMatrix <- paMatrixList[[k]]
    vec <- paste(c(paMatrix[i,], " 1;"), collapse = "")
    write(vec, filename, append = ifelse(i == 1, F, T))
  }
}

# make Pradel models
files <- list.files(, pattern = "diversityMatrix")
neo.pradsenList <- list()
# neo.pradLamdaList <- list()
for(i in 1:length(files)){
  neo <- convert.inp(files[i])
  neo.pradsenList[[i]] <- process.data(neo, model= "Pradsen") # Pradel's survival and senority
  # neo.pradLamdaList[[i]] <- process.data(neo, model = "Pradlambda") # Pradel's survival and "growth"
}


############model specifications####################
#time-varying
Phi.t = list(formula=~time)
Gamma.t = list(formula=~time)
p.t = list(formula=~time)
L.t = list(formula=~time)

#run time-varying model ONLY (two different parameterizations to get phi, gamma and lambda)
timeList <- list()
#timeDList <- list()
for(i in 1:length(files)){
  timeList[[i]] = mark(neo.pradsenList[[i]],  model.parameters=list(Phi = Phi.t, p = p.t, Gamma = Gamma.t), delete = T, brief = T)
  #	timeDList[[i]]=mark(neo.pradLamda,  model.parameters=list(Phi=Phi.t, p=p.t, Lambda=L.t))
}

dim(timeList[[1]]$results$real)
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

