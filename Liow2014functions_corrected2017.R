#Code for a basic Pradel model, Foote rates from (Foote 2000 Palaoebiology), Alroy three timer rates from (Alroy 2009 PNAS)
#30 05 2013 Lee Hsiang Liow

## CORRECTED 2017 (see Liow and Finarelli 2017)

#####################################################
#Foote function with corrections (Liow and Finarelli 2017)
footerates <- function(FA, LA, timeint){
	da = cbind(FA, LA)
	NFt = NA # number of top crosser (this was previously misnamed re Foote 2000)
	Nbt = NA # number of double crosser (this was previously misnamed re Foote 2000)
	NbL = NA # number of bottom crosser
	q = NA # per capita extinction rate
	p = NA # per capita origination rate
	for (i in 1:(length(timeint)-1)){
		Ft = da[which(da[,1] <= timeint[i] & da[,1] > timeint[i+1] & da[,2] <= timeint[i+1]),] #top crosser
		bt = da[which(da[,1] > timeint[i] & da[,2] <= timeint[i+1]),] #double crosser
		bL = da[which(da[,1] > timeint[i] & da[,2] <= timeint[i] & da[,2] > timeint[i+1]),] #bottom crosser
		#Ft=da[which(da[,1]>timeint[i] & da[,2]>timeint[i+1] & da[,2]<=timeint[i]),] #old code
		#bt= da[which(da[,1]>timeint[i] & da[,2]<timeint[i+1]),] #old code
		#bL=da[which(da[,1]<timeint[i] & da[,1]>=timeint[i+1] & da[,2]<timeint[i+1]),] #old code
		
		# If First-top is not empty, insert 1 to cell
		if (!is.null(nrow(Ft))){
			NFt[i] = nrow(Ft)
		}else{
			NFt[i] = 1
		}
		
		# If Bottom and Top is not empty, insert 1 to cell
		if (!is.null(nrow(bt))){
			Nbt[i] = nrow(bt)
		}else{
			Nbt[i] = 1
		}
		
		# If Bottom-last is not empty, insert 1 to cell
		if (!is.null(nrow(bL))){
			NbL[i] = nrow(bL)
		}else{
			NbL[i] = 1
		}
	}
		
	footeq = -log(Nbt / (NbL + Nbt)) #extinction
	footep = -log(Nbt / (NFt + Nbt)) #origination
	res = cbind(footeq, footep)
  	colnames(res) = c("ExtRate", "OrigRate")
  	rownames(res) = timeint[1:(length(timeint) - 1)]
  	#return(res)
  	return(list(NFt = NFt, NbL = NbL, Nbt = Nbt, Nb = Nbt + NbL, Nt = Nbt + NFt, OrigRate = footep, 	ExtRate = footeq))

}


##################### 
#Alroy Three Timer Rates


alroyrates <- function(data, times){
	timeint <- length(times)

	#timeint = number of time bins total	
	#Five Classes of  "occurrences"
  	SIB = apply(data, 2, sum) # one timers
  	twoTi = NULL          # two timers (before and within)
  	twoTipo = NULL        # two timers (within and after)
  	threeT = NULL         # three timers
  	partT = NULL          # part timers
    
  	for (i in 2:(timeint - 1)){ # second time interval to second to last only 
 		twoTi[i] = length(which(data[,i] == 1 & data[,i-1] == 1))
     	twoTipo[i] = length(which(data[,i] == 1 & data[,i+1] == 1))
     	threeT[i] = length(which(data[,i] == 1 & data[,i-1] == 1 & data[,i+1] == 1))
     	partT[i] = length(which(data[,i] == 0 & data[,i-1] == 1 & data[,i+1] == 1))
  	}

	twoTi = twoTi[2:(timeint - 1)]
	twoTipo = twoTipo[2:(timeint - 1)]
	threeT = threeT[2:(timeint - 1)]
	partT = partT[2:(timeint - 1)]

	Ap = threeT / (threeT + partT) # preservation probability

	mu = log(twoTi / threeT) # this is from the second time interval
	lam = log(twoTipo / threeT)

	Cmu = NULL # corrected extinction rate
	Clam = NULL # corrected origination rate

### *** Code modified from Liow and Finarelli 2014 (by TMS) to include ln(preservation rate) according to Alroy 2010 *** ###
	# for (i in 1:timeint-2) # old code
	# {Cmu[i]=mu[i]+Ap[i+1] # old code #values from second time interval but origination not available for second interval for lam
	# Clam[i]=lam[i]+Ap[i-1]} # old code


	for (i in 1:timeint - 2){
		Cmu[i] = mu[i] + log(Ap[i + 1]) # ln(preservation rate) 
		#values from second time interval but origination not available for second interval for lam
		Clam[i] = lam[i] + log(Ap[i - 1]) # ln(preservation rate)
		}
	
	res = cbind(mu, lam, Cmu, Clam)
  	colnames(res) = c("mu", "lam", "ExtRate", "OrigRate")
  	rownames(res) = times[2:(length(times) - 1)]
  	return(res)
}

