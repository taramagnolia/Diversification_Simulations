#Code for a basic Pradel model, Foote rates from (Foote 2000 Palaoebiology), Alroy three timer rates from (Alroy 2009 PNAS)
#30 05 2013 Lee Hsiang Liow

#####################################################
#Foote function


footerates <- function(FA, LA, timeint){
  
  da=cbind(FA, LA)
  NFt=NA            #bottom crosser
  Nbt=NA            #double crosser
  NbL=NA              #upper crosser
  q =NA
  p=NA
  
  for (i in 1:(length(timeint)-1))
    
  {
    Ft=da[which(da[,1]>timeint[i] & da[,2]>timeint[i+1] & da[,2]<=timeint[i]),]          #bottom crosser
    bt= da[which(da[,1]>timeint[i] & da[,2]<timeint[i+1]),]                             #double crosser
    bL=da[which(da[,1]<timeint[i] & da[,1]>=timeint[i+1] & da[,2]<timeint[i+1]),]         #upper crosser
    
    if (!is.null(nrow(Ft)))  
    {NFt[i]=nrow(Ft)
    } else{
      NFt[i]=1}
    
    if (!is.null(nrow(bt)))  
    {Nbt[i]=nrow(bt)
    }else{
      Nbt[i]=1}
    
    if (!is.null(nrow(bL)))  
    {NbL[i]=nrow(bL)
    }else{
      NbL[i]=1}
      }
  
  footeq=-log( Nbt/(NFt+ Nbt))    
  footep=-log( Nbt/(NbL+ Nbt))
  res=cbind(footeq, footep)
  colnames(res)=c("ExtRate","OrigRate")
  rownames(res)=timeint[1:(length(timeint)-1)]
  return(res)
}





##################### 
#Alroy Three Timer Rates


alroyrates <- function(data, times){
	timeint <- length(times)

	#timeint = number of time bins total	
	#Five Classes of  "occurrences"
  	SIB=apply(data,2,sum) #one timers
  	twoTi=NULL          #two timers (before and within)
  	twoTipo=NULL          #two timers (within and after)
  	threeT=NULL            #three timers
  	partT=NULL              #part timers
    
  	for (i in 2:(timeint-1))#second time interval to second to last only 
  
  		{  twoTi[i]=length(which(data[,i]==1 & data[,i-1]==1))
     		twoTipo[i]=length(which(data[,i]==1 & data[,i+1]==1))
     		threeT[i]=length(which(data[,i]==1 & data[,i-1]==1 & data[,i+1]==1))
     		partT[i]=length(which(data[,i]==0 & data[,i-1]==1 & data[,i+1]==1))
  		}

	twoTi=twoTi[2:(timeint-1)]
	twoTipo=twoTipo[2:(timeint-1)]
	threeT=threeT[2:(timeint-1)]
	partT=partT[2:(timeint-1)]

	Ap=threeT/(threeT+partT) #preservation probability

	mu=log(twoTi/threeT) #this is from the second time interval
	lam=log(twoTipo/threeT)

	Cmu=NULL #corrected extinction rate
	Clam= NULL #corrected origination rate

	for (i in 1:timeint-2)
		{Cmu[i]=mu[i]+log(Ap[i+1]) #values from second time interval but origination not available for second interval for lam
		Clam[i]=lam[i]+log(Ap[i-1])}
	
	res=cbind(mu, lam, Cmu, Clam)
  	colnames(res)=c("mu","lam", "ExtRate", "OrigRate")
  	rownames(res)=times[2:(length(times)-1)]
  	return(res)
}

