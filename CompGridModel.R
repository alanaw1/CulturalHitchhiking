#### Computational Grid Model #### 
### Alan J. Aw, Marcus W. Feldman, Tian Chen Zeng ###
### This Version: 1/30/2018 ####

######### THIS CODE PRODUCES A FOUR-PLOT GRAPH ###########

### Set seed for reproducibility ###
set.seed(123) 

### Variables ###
R = 100 # number of reps
T = 60 # number of generations, T 
C = 100  # number of cultural groups, |C|
H = 2500  # number of haplotypes, |H|
N.total = 2e4 # N_total
rho = 18e-5 # the average death rate, rho 
          # Note: always multiply this by (C-1) to compute mean decline in a cultural group size
          # after one generation.
alpha.c = 1.00 + 1e-2*((1:C)-1) # differential fitness of cultural groups (1e-2 change to 0)
mu = 2e-3  # haplogroup mutation rate, mu  

### Helper Functions ###
#library("entropy")
## Normalized Shannon index (see p. 7 of Supplementary Information)
normalized.shannonIndex = function(vec) {
  freqs = vec/sum(vec)
  vec.entropy=entropy::entropy(freqs)
  n.of.types=0
  for (i in 1:length(vec)) {
    if (vec[i] != 0) {
      n.of.types=n.of.types+1
    }
  }
  if (n.of.types == 1) { 
    return(0)
  } else {
    return(vec.entropy / log(n.of.types))
  }
}

## Richness (counts number of distinct types whose sizes are positive)
richness = function(vec) {
  n.types=0
  for (i in 1:length(vec)) {
    if (vec[i] != 0) {
      n.types=n.types+1
    }
  }
  return(n.types)
}

### Population Size and Positive Part ###
N=array(0,dim=c(C,H,T,R),dimnames=list(culture=1:C,haplogroup=1:H,generation=1:T,rep=1:R))

## For patrilineal, uncomment next three lines and comment the other option.
for (c in 1:C) {
  N[c,c,1,] = N.total / C
}

## For non-patrilineal, uncomment code below and comment code above. 
#for (c in 1:C) { 
#  for (h in 1:C) {
#    N[c,h,1,] = N.total / (C*C)
#  }
#}

## Function to obtain positive part of a float ##
posPart = function(num) { max(num,0) }

### Tracks relevant parameteric trajectories of male populations ###
recPopSize = array(0,dim=c(T,R))
recEntropy = array(0,dim=c(T,R))
recRichness = array(0,dim=c(T,R))
recCultRichness = array(0,dim=c(T,R))
recCultEntropy = array(0,dim=c(T,R))

### Implement simulation and save results ###
## Repeat This ##
for(r in 1:R) {
  
  ## Loop through time ##
  for(t in 2:T) {
    
    ## Killing ## 
    for(c in 1:C) {
      deathRate = alpha.c[c] * rho * sum(N[-c,,t-1,r]) * sqrt(sum(N[c,,t-1,r])) # sqrt of cultural group
      nDeaths = rpois(n=1,lambda=deathRate)
      if (sum(N[c,,t-1,r])>0) { # in case the cultural group is empty already
        newDeaths = rmultinom(n=1,size=nDeaths,prob=N[c,,t-1,r]) 
        #newHaplotypeSizes = round(NewBirths[c,]-newDeaths)
        N[c,,t,r] = sapply((N[c,,t-1,r]-newDeaths),posPart)
      } else {
        N[c,,t,r] = N[c,,t-1,r]
      }
    }
    ## Mutation ## 
    for(c in 1:C) {
      argmax = 0
      maxhaplotype = 0
      for(h in 1:H) {
        if (N[c,h,t,r] > maxhaplotype) {
          maxhaplotype = N[c,h,t,r]
          argmax = h
        }
      }
      if (argmax != 0){
        move = round(mu*N[c,argmax,t,r])
        N[c,argmax,t,r] = N[c,argmax,t,r]-move
        allExceptArgmax = rep(1,H)
        allExceptArgmax[argmax] = 0
        mutants = rmultinom(n=1,move,prob=as.vector(allExceptArgmax))
        for (h in 1:H) {
          N[c,h,t,r] = N[c,h,t,r] + mutants[h]
        }
      }
    }
    
    ## Weed off small cultural groups ##
    fusionCounter = 0 # Should we do a fusion counter? Alan on 1/30/2018
    for (c in 1:C) {
      if (sum(N[c,,t,r]) < 20) { 
        ## to merge with the neighbouring cultural group (no geography). Added by Alan on 1/30/2018
        u = 2*rbinom(1,1,0.5) # pick +/- 1 with prob 0.5 each
        cPrime = c + u 
        if (cPrime > C) cPrime = 1 # since c must have been C and u = +1
        if (cPrime == 0) cPrime = C # since c must have been 1 and u = -1
        for (h in 1:H) { 
          N[cPrime,h,t,r] = N[cPrime,h,t,r] + N[c,h,t,r] # shift those from cultural group c to c'
          N[c,h,t,r] = 0 # empty cultural group c
        }
      }
    }
    
    ## Scaling / Normalisation ## 
    rescale = N.total / sum(N[,,t,r])
    for (c in 1:C) {
      for (h in 1:H) {
        N[c,h,t,r] = N[c,h,t,r] * rescale
      }
    }
    
    ## Replenish empty cultural groups ## (Added by Alan Aw on 9/26/2017)
    for (c in 1:C) {
      cg.size = sum(N[c,,t,r])
      if (cg.size == 0) {
        largest.c = 0
        for (o in 1:C) {
          if (sum(N[o,,t,r]) > largest.c) {
            largest.c = o
          }
        }
        N[c,,t,r] = 0.5 * N[largest.c,,t,r]
        N[largest.c,,t,r] = 0.5 * N[largest.c,,t,r]
      }
    }
  } ## End of time loop
  
  ## Save results ##
  totVsTime = apply(N[,,,r],3,sum)
  Hap = array(dim=c(H,T),dimnames=list(haplogroup=1:H,generation=1:T)) # haplogroup size array
  for(t in 1:T) {
    for(h in 1:H) {
      Hap[h,t] = sum(N[,h,t,r])
    }
  }
  Cult = array(dim=c(C,T),dimnames=list(culturalgroup=1:C,generation=1:T)) # cultural group size array
  for(t in 1:T) {
    for(c in 1:C) {
      Cult[c,t] = sum(N[c,,t,r])
    }
  }
  shan = apply(Hap[,],2,normalized.shannonIndex) # haplogroup diversity over time 
  hap.rich = apply(Hap[,],2,richness) # haplotype richness over time
  cult.rich = apply(Cult[,],2,richness) # cultural group richness over time
  cult.shan = apply(Cult[,],2,normalized.shannonIndex) # cultural group diversity over time
  recPopSize[,r] = totVsTime
  recEntropy[,r] = shan
  recRichness[,r] = hap.rich
  recCultRichness[,r] = cult.rich
  recCultEntropy[,r] = cult.shan
}

### Plot results ###
avePopSize = c()
for (t in 1:T) {
  avePopSize = c(avePopSize,mean(recPopSize[t,])) 
}
aveEntropy = c()
for (t in 1:T) {
  aveEntropy = c(aveEntropy,mean(recEntropy[t,])) 
}
upperLim.avePopSize = c()
for (t in 1:T) {
  upperLim.avePopSize=c(upperLim.avePopSize, mean(recPopSize[t,]) + 2 * sd(recPopSize[t,]))
}
lowerLim.avePopSize = c()
for (t in 1:T) {
  lowerLim.avePopSize=c(lowerLim.avePopSize, mean(recPopSize[t,]) - 2 * sd(recPopSize[t,]))
}
upperLim.aveEntropy = c()
for (t in 1:T) {
  upperLim.aveEntropy=c(upperLim.aveEntropy, mean(recEntropy[t,]) + 2 * sd(recEntropy[t,]))
}
lowerLim.aveEntropy = c()
for (t in 1:T) {
  lowerLim.aveEntropy=c(lowerLim.aveEntropy, max(0,mean(recEntropy[t,]) - 2 * sd(recEntropy[t,])))
}

aveRichness = c()
for (t in 1:T) {
  aveRichness = c(aveRichness,mean(recRichness[t,])) 
}
upperLim.aveRichness = c()
for (t in 1:T) {
  upperLim.aveRichness=c(upperLim.aveRichness, mean(recRichness[t,]) + 2 * sd(recRichness[t,]))
}
lowerLim.aveRichness = c()
for (t in 1:T) {
  lowerLim.aveRichness=c(lowerLim.aveRichness, mean(recRichness[t,]) - 2 * sd(recRichness[t,]))
}

aveCultRichness = c()
for (t in 1:T) {
  aveCultRichness = c(aveCultRichness,mean(recCultRichness[t,])) 
}
upperLim.aveCultRichness = c()
for (t in 1:T) {
  upperLim.aveCultRichness=c(upperLim.aveCultRichness, min(C,mean(recCultRichness[t,]) + 2 * sd(recCultRichness[t,])))
}
lowerLim.aveCultRichness = c()
for (t in 1:T) {
  lowerLim.aveCultRichness=c(lowerLim.aveCultRichness, max(0,mean(recCultRichness[t,]) - 2 * sd(recCultRichness[t,])))
}

aveCultEntropy = c()
for (t in 1:T) {
  aveCultEntropy = c(aveCultEntropy,mean(recCultEntropy[t,])) 
}
upperLim.aveCultEntropy = c()
for (t in 1:T) {
  upperLim.aveCultEntropy=c(upperLim.aveCultEntropy, mean(recCultEntropy[t,]) + 2 * sd(recCultEntropy[t,]))
}
lowerLim.aveCultEntropy = c()
for (t in 1:T) {
  lowerLim.aveCultEntropy=c(lowerLim.aveCultEntropy, max(0,mean(recCultEntropy[t,]) - 2 * sd(recCultEntropy[t,])))
}


par(mfrow=c(2,2),mar=c(4,4,1.6,0.5), oma=c(1.2,2,1,1))
# Plot population size trajectory
#plot(1:T,avePopSize,type="l",log="y",ylim=c(1,1e8),xlab="Number of Generations",ylab="Population Size",col="black",lwd=2,main="Number of Individuals through Time")
#lines(1:T,upperLim.avePopSize,type="l", lty=2,col="black")
#lines(1:T,lowerLim.avePopSize,type="l",lty=2,col="black")
# Plot haplogroup richness trajectory
plot(1:T,aveRichness,type="n",ylim=c(0,max(recRichness)),xlab="",cex.axis=1.15,cex.lab=1.3,cex.main=1.2,
     ylab="No. of Distinct Haplogroups",main="A")
for(r in 1:R) {
  lines(1:T,recRichness[,r],lwd=1,col="gold") # lightpink
}
lines(1:T,upperLim.aveRichness,type="l", lty=2,col="black")
lines(1:T,lowerLim.aveRichness,type="l", lty=2,col="black")
lines(1:T,aveRichness,type="l",lwd=2,col="darkgoldenrod4") # red

# Plot haplotype diversity trajectory
plot(1:T,aveEntropy,type="n",ylim=c(0,max(recEntropy)),xlab="",cex.axis=1.15,cex.lab=1.3,cex.main=1.2,
     ylab="Haplogroup Diversity",main="B")
for(r in 1:R) {
  lines(1:T,recEntropy[,r],lwd=1,col="gold") # lightpink
}
lines(1:T,upperLim.aveEntropy,type="l", lty=2,col="black")
lines(1:T,lowerLim.aveEntropy,type="l",lty=2,col="black")
lines(1:T,aveEntropy,type="l",lwd=2,col="darkgoldenrod4") # red

# Plot cultural group richness trajectory
plot(1:T,aveCultRichness,type="n",ylim=c(0,max(recCultRichness)),xlab="Number of Generations",cex.axis=1.15,cex.lab=1.3,cex.main=1.2,
     ylab="No. of Distinct Cultural Groups",main="C")
for(r in 1:R) {
  lines(1:T,recCultRichness[,r],lwd=1,col="gold") # lightpink
}
lines(1:T,upperLim.aveCultRichness,type="l", lty=2,col="black")
lines(1:T,lowerLim.aveCultRichness,type="l", lty=2,col="black")
lines(1:T,aveCultRichness,type="l",lwd=2,col="darkgoldenrod4") # red

# Plot cultural group diversity trajectory
plot(1:T,aveCultEntropy,type="n",ylim=c(0,max(recCultEntropy)),xlab="Number of Generations",cex.axis=1.15,
     cex.lab=1.3,cex.main=1.2,ylab="Cultural Group Diversity",main="D")
for(r in 1:R) {
  lines(1:T,recCultEntropy[,r],lwd=1,col="gold") # lightpink
}
lines(1:T,upperLim.aveCultEntropy,type="l", lty=2,col="black")
lines(1:T,lowerLim.aveCultEntropy,type="l",lty=2,col="black")
lines(1:T,aveCultEntropy,type="l",lwd=2,col="darkgoldenrod4")  # red

###########%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##############
##### Code for Plotting the Haplogroup Lineage Trajectories ##### 
#### RUN ONLY AFTER SAVING THE FOUR-PLOT GRAPH ####

dev.off()
par(mfrow=c(1,1),mar=c(4,4,1.6,0.5), oma=c(0.5,0.5,0.5,0.5),family="Helvetica")
K=H
ave.Hap = array(dim=c(K,T),dimnames=list(haplogroup=1:K,generation=1:T)) # k is ranked haplogroup
## Sorting creates the matrix of ranks of haplogroup sizes according to who's largest by T generations
rankingMatrix = array(0,dim=c(R,H)) # to track the ranking of haplogroup sizes each rep
tempRankingVector = array(0,dim=c(H)) # to collect haplogroup sizes at generation T
for(r in 1:R) {
  for(h in 1:H){
    tempRankingVector[h]=sum(N[,h,T,r])
  }
  rankingMatrix[r,] = sort(tempRankingVector,index.return=TRUE)$ix
}

## Computes average ranked haplogroup sizes across time
for (t in 1:T) { 
  for (k in 1:K) {
    localSum = 0
    for (r in 1:R) { 
      localSum = localSum + sum(N[,rankingMatrix[r,k],t,r]) 
    }
    ave.Hap[k,t] = localSum / R
  }
}
plot(1:T,ave.Hap[1,],type="n",ylim=c(0,max(ave.Hap)),xlab="Number of Generations",cex.lab=1.3,
     cex.axis=1.2,ylab="Number of Individuals")
for(h in 1:H) {
  lines(1:T,ave.Hap[h,],lwd=1,col="blue")
}
#library(ggplot2)
#qplot()

### Plotting side-by-side (one rep vs average of 100 reps) ###
par(mfrow=c(1,2),mar=c(2,4,1.6,0.5), oma=c(1.2,2,1,1),pty="s") # remember to save 6.8 x 4.2 portrait
## Average of 100 reps 
plot(1:T,ave.Hap[1,],type="n",ylim=c(0,max(ave.Hap)+3000),xlab="Number of Generations",cex.lab=1.3,
     cex.axis=1.2,ylab="Number of Individuals",main="E")
for(h in 1:H) {
  lines(1:T,ave.Hap[h,],lwd=1,col="deeppink2")
}
## First rep trajectories only
plot(1:T,ave.Hap[1,],type="n",ylim=c(0,max(ave.Hap)+3000),xlab="Number of Generations",cex.lab=1.3,
     cex.axis=1.2,ylab="",main="F")
the.hap = array(dim=c(H,T))
for (h in 1:H) {
  for (t in 1:T) {
    the.hap[h,t] = sum(N[,h,t,1])
  }
  lines(1:T,the.hap[h,],col="deeppink2")
}


