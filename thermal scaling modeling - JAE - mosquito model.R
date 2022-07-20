require(deSolve) 
require(scales)
require(ggplot2)
require(dplyr)
require(RColorBrewer)
require(viridis)
require(rcartocolor)
require(gridExtra)


temp <- seq(0,45,0.1) #temperature range

ind.matrix <- matrix(nrow=1000,ncol=length(temp)) # sets length

set.seed(1234)


#### POPULATION-LEVEL MODEL ####

# Biting rate: a, Briere

a.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
a.c <- runif(dim(a.matrix)[1], min=0.5,max=1.3)
a.tmin <- runif(dim(a.matrix)[1], min=0,max=10)
a.tmax <- a.tmin + runif(dim(a.matrix)[1], min=15, max=35)

# Divide a by 1000 to get approximate scale for contact rate #

for(j in 1:dim(a.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<a.tmin[j] || temp[i]>a.tmax[j],      
           a.matrix[j,i]<-0,
           a.matrix[j,i]<- (a.c[j]*temp[i]*(temp[i]-a.tmin[j])*((a.tmax[j]-temp[i])^(1/2)) / 3000) ) # Divide by 1000 to get proper scale for contact rate
  }
}


# Vector competence: bc

bc.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
bc.c <- runif(dim(bc.matrix)[1], min=0.5,max=1.3)
bc.tmin <- runif(dim(bc.matrix)[1], min=0,max=10)
bc.tmax <- bc.tmin + runif(dim(bc.matrix)[1], min=15, max=35)

# Divide a by 1000 to get approximate scale for contact rate #

for(j in 1:dim(bc.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<bc.tmin[j] || temp[i]>bc.tmax[j],      
           bc.matrix[j,i]<-0,
           bc.matrix[j,i] <- -bc.c[j]*(temp[i]-bc.tmin[j])*(temp[i]-bc.tmax[j]) / 500)  # Divide by 1000 to get proper scale for contact rate
  }
}


# Lifespan #

lifespan.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
lifespan.c <- runif(dim(lifespan.matrix)[1], min=0.5,max=1.3)
lifespan.tmin <- runif(dim(lifespan.matrix)[1], min=0,max=10)
lifespan.tmax <- lifespan.tmin + runif(dim(lifespan.matrix)[1], min=15, max=35)


# Divide lifespan by 1,000,000 to get approximate scale for probability of infection #

for(j in 1:dim(lifespan.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<lifespan.tmin[j] || temp[i]>lifespan.tmax[j],      
           lifespan.matrix[j,i]<-0,
           lifespan.matrix[j,i]<- -lifespan.c[j]*(temp[i]-lifespan.tmin[j])*(temp[i]-lifespan.tmax[j]) / 2)
  }
}

# mortality rate mu
temp.mu.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

temp.mu.matrix[,] <- 1/lifespan.matrix[,]

mu.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# set any value greater than 5 to 5
for(i in 1:length(temp.mu.matrix)){
  ifelse(temp.mu.matrix[i] < 5,
         mu.matrix[i] <- temp.mu.matrix[i],
         mu.matrix[i] <- 5)
}



# Host PDR: PDR # 

PDR.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
PDR.c <- runif(dim(PDR.matrix)[1], min=0.5,max=1.3)
PDR.tmin <- runif(dim(PDR.matrix)[1], min=0,max=10)
PDR.tmax <- PDR.tmin + runif(dim(PDR.matrix)[1], min=15, max=35)

for(j in 1:dim(PDR.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<PDR.tmin[j] || temp[i]>PDR.tmax[j],      
           PDR.matrix[j,i]<-0,
           PDR.matrix[j,i]<- (PDR.c[j]*temp[i]*(temp[i]-PDR.tmin[j])*((PDR.tmax[j]-temp[i])^(1/2)) ) / 5000 )
  }
}



# EIP 
EIP.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

EIP.matrix[,] <- 1/PDR.matrix[,]



# pea #

pea.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
pea.c <- runif(dim(pea.matrix)[1], min=0.5,max=1.3)
pea.tmin <- runif(dim(pea.matrix)[1], min=0,max=10)
pea.tmax <- pea.tmin + runif(dim(pea.matrix)[1], min=15, max=35)


# Divide pea by 1,000,000 to get approximate scale for probability of infection #

for(j in 1:dim(pea.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<pea.tmin[j] || temp[i]>pea.tmax[j],      
           pea.matrix[j,i]<-0,
           pea.matrix[j,i]<- -pea.c[j]*(temp[i]-pea.tmin[j])*(temp[i]-pea.tmax[j]) / 500)
  }
}



# Host EFD: EFD # 

EFD.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
EFD.c <- runif(dim(EFD.matrix)[1], min=0.5,max=1.3)
EFD.tmin <- runif(dim(EFD.matrix)[1], min=0,max=10)
EFD.tmax <- EFD.tmin + runif(dim(EFD.matrix)[1], min=15, max=35)

for(j in 1:dim(EFD.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<EFD.tmin[j] || temp[i]>EFD.tmax[j],      
           EFD.matrix[j,i]<-0,
           EFD.matrix[j,i]<- (EFD.c[j]*temp[i]*(temp[i]-EFD.tmin[j])*((EFD.tmax[j]-temp[i])^(1/2)) ) / 500 )
  }
}


# Host MDR: MDR # 

MDR.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
MDR.c <- runif(dim(MDR.matrix)[1], min=0.5,max=1.3)
MDR.tmin <- runif(dim(MDR.matrix)[1], min=0,max=10)
MDR.tmax <- MDR.tmin + runif(dim(MDR.matrix)[1], min=15, max=35)

for(j in 1:dim(MDR.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<MDR.tmin[j] || temp[i]>MDR.tmax[j],      
           MDR.matrix[j,i]<-0,
           MDR.matrix[j,i]<- (MDR.c[j]*temp[i]*(temp[i]-MDR.tmin[j])*((MDR.tmax[j]-temp[i])^(1/2)) ) / 7000 )
  }
}



N <- 10000 # constant human population
r <- 0.01 # constant human recovery rate


# R0 calculations

R0.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

for(j in 1:dim(R0.matrix)[1]){
  for(i in 1:length(temp)){
    
    R0.matrix[j,i] = (((a.matrix[j,i]^2)*bc.matrix[j,i]*exp((-mu.matrix[j,i])/PDR.matrix[j,i])*EFD.matrix[j,i]*pea.matrix[j,i]*MDR.matrix[j,i])/
    (N*r*(mu.matrix[j,i]^3)))^(1/2)
      

  }
}



# calculate infected days for our individual-level parasitism metric

###### Logistic equation to generate proportion of infectious mosquitos over time (b) (Ohm et al. 2018)

# b = bmax / (1 + exp(-k* (t-tm)) )

# where bmax is the upper asymptote (or max transmission prevalence i,.e. vector competence, or inf. or trans. prob.)
# tm = t50, the time at which 50% of infected vectors have become infectious (EIP)
# and k is a rate for fitted logistic models.

TT <- seq (0,365,1) # 1 year time length, same as survival curves in other script #


k = 1 # Setting k = 1. Unimodal temperature dependence in Ohm et al. 2018

survival<-function(t,y,p){ # function to simulate exponential survival curve (constant hazard)
  mu<-p[1];
  N<-y[1];
  
  dN <- -mu*N
  
  list(c(dN))
}

N0 <- 1

infected.days <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# need to run through for all 1000 sims 
for(j in 1:dim(ind.matrix)[1]){

  print(j)
  # This gets proportion infected over time at each temperature
  temp.prop.infected.matrix <- matrix(nrow=length(TT),ncol=length(temp)) # temporary for each simulation
  
  bc <- bc.matrix[j,] # set vector competence for each sim
  EIP <- EIP.matrix[j,] # set vector competence for each sim
  
  
for(i in 1:length(temp)){
  for(t in 1:length(TT)){
    temp.prop.infected.matrix[t,i]<- (bc[i] / (1 + exp(-k * (TT[t]-EIP[i]))))
  }
}
  
  # Then survival curve over time at each temperature
  temp.surv.curves.matrix <- matrix(nrow=length(TT),ncol=length(temp)) 
  
  mu <- mu.matrix[j,]
  
  for(i in 1:length(temp)){
    parms <- mu[i]
    
    
    predictions <-lsoda(N0,TT,survival,parms)	
    
    # If probability of being alive is < 0.00001, set to 0 so that R doesn't run into problems close to 0.
    surv.curve <- c()
    for(m in 1:length(TT)){
      surv.curve[m] <- ifelse(predictions[m,2]>0.00001,surv.curve[m]<-predictions[m,2],surv.curve[m]<-0)
    }
    
    temp.surv.curves.matrix[,i] <- surv.curve
  }
  
  inf.days <- c()

  for(i in 1:length(temp)){
    weighted.prop.inf <- temp.surv.curves.matrix[,i]*temp.prop.infected.matrix[,i]
    temp.infected.days <- sum(weighted.prop.inf)
    inf.days[i] <- temp.infected.days
  }
  
  infected.days[j,] <- inf.days
  
}


R0.topt <- c()
ind.topt <- c()
vec <- c() # 0 and 1s for if the epidemic spreads (If R0 = 0 or not)

for(j in 1:dim(R0.matrix)[1]){
  R0.topt[j] <- temp[which.max((R0.matrix[j,]))]
  ind.topt[j] <- temp[which.max((infected.days[j,]))]
  vec[j] <- isTRUE(R0.topt[j]>0)
}


data.mat <- data.frame(ind.topt,R0.topt,vec)


# filter out ones where R0 was 0 
data.mat <- data.mat %>% 
  filter(vec==TRUE)

# Save the thermal optima data in a CSV to plot after #
write.csv(data.mat, file="~/mosquito_simulations.csv")
