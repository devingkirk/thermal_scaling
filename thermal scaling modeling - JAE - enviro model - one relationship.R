require(deSolve) 
require(scales)
require(ggplot2)
require(dplyr)
require(RColorBrewer)
require(viridis)
require(rcartocolor)
require(gridExtra)

# Environmentally-transmitted systems #

##### Describe a curve for individual-level parasitism (i.e., average parasite load)

temp <- seq(0,45,0.1) #temperature range

ind.matrix <- matrix(nrow=1000,ncol=length(temp))

set.seed(1222)

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax

ind.c <- runif(dim(ind.matrix)[1], min=0.5,max=1.3)
ind.tmin <- runif(dim(ind.matrix)[1], min=0,max=10)
ind.tmax <- ind.tmin + runif(dim(ind.matrix)[1], min=15, max=35)

for(j in 1:dim(ind.matrix)[1]){
for(i in 1:length(temp)){
  ifelse(temp[i]<ind.tmin[j] || temp[i]>ind.tmax[j],      
         ind.matrix[j,i]<-0,
         ind.matrix[j,i]<- (-ind.c[j]*(temp[i]-ind.tmin[j])*(temp[i]-ind.tmax[j]) ) / 2 )
}
}


#### POPULATION-LEVEL MODEL ####

# Contact rate: chi #

chi.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
chi.c <- runif(dim(chi.matrix)[1], min=0.5,max=1.3)
chi.tmin <- runif(dim(chi.matrix)[1], min=0,max=10)
chi.tmax <- chi.tmin + runif(dim(chi.matrix)[1], min=15, max=35)

# Divide chi by 1000 to get approximate scale for contact rate #

for(j in 1:dim(chi.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<chi.tmin[j] || temp[i]>chi.tmax[j],      
           chi.matrix[j,i]<-0,
           chi.matrix[j,i]<- (chi.c[j]*temp[i]*(temp[i]-chi.tmin[j])*((chi.tmax[j]-temp[i])^(1/2)) / 2000) ) # Divide by 1000 to get proper scale for contact rate
  }
}


# Probability of infection: sigma #

sigma.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
sigma.c <- runif(dim(sigma.matrix)[1], min=0.5,max=1.3)
sigma.tmin <- runif(dim(sigma.matrix)[1], min=0,max=10)
sigma.tmax <- sigma.tmin + runif(dim(sigma.matrix)[1], min=15, max=35)


# Divide sigma by 1,000,000 to get approximate scale for probability of infection #

for(j in 1:dim(sigma.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<sigma.tmin[j] || temp[i]>sigma.tmax[j],      
           sigma.matrix[j,i]<-0,
           sigma.matrix[j,i]<- -sigma.c[j]*(temp[i]-sigma.tmin[j])*(temp[i]-sigma.tmax[j])  / 10000000)
  }
}



# Parasite-induced mortality: alpha # 

# Make this proportional to load 
alpha.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

alpha.matrix[,] <- ind.matrix[,]/1000


# Parasites released after host death

# Make this proportional to load 
omega.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

omega.matrix[,] <- ind.matrix[,]


# Parasites released over time

lambda.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax

lambda.c <- runif(dim(lambda.matrix)[1], min=0.5,max=1.3)
lambda.tmin <- runif(dim(lambda.matrix)[1], min=0,max=10)
lambda.tmax <- lambda.tmin + runif(dim(lambda.matrix)[1], min=15, max=35)

for(j in 1:dim(lambda.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<lambda.tmin[j] || temp[i]>lambda.tmax[j],      
           lambda.matrix[j,i]<-0,
           lambda.matrix[j,i]<- (-lambda.c[j]*(temp[i]-lambda.tmin[j])*(temp[i]-lambda.tmax[j]) ) / 40 )
  }
}



# Host birth.rate: birth.rate # 

birth.rate.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
birth.rate.c <- runif(dim(birth.rate.matrix)[1], min=0.5,max=1.3)
birth.rate.tmin <- runif(dim(birth.rate.matrix)[1], min=0,max=10)
birth.rate.tmax <- birth.rate.tmin + runif(dim(birth.rate.matrix)[1], min=15, max=35)

for(j in 1:dim(birth.rate.matrix)[1]){
  for(i in 1:length(temp)){
    ifelse(temp[i]<birth.rate.tmin[j] || temp[i]>birth.rate.tmax[j],      
           birth.rate.matrix[j,i]<-0,
           birth.rate.matrix[j,i]<- (birth.rate.c[j]*temp[i]*(temp[i]-birth.rate.tmin[j])*((birth.rate.tmax[j]-temp[i])^(1/2)) ) / 100 )
  }
}


# For host background mortality rate, assuming mortality follows an inverted quadratic

mu.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
mu.inter <- runif(dim(mu.matrix)[1], min=.1,max=1) # this parameter also known as c
mu.n.slope <- runif(dim(mu.matrix)[1], min=.02,max=0.03) # this parameter also known as b
mu.qd <-  runif(dim(mu.matrix)[1], min=0.0008, max=0.0009) #this parameter also known as a


for(j in 1:dim(mu.matrix)[1]){
  for(i in 1:length(temp)){

temp.value <- c()
temp.value[i] <- (mu.qd[j]*(temp[i])^2)-(mu.n.slope[j]*temp[i])+mu.inter[j]

ifelse(temp.value/5<0.00667, #if 1/5 value is less than 0.00667 (corresponding to a lifespan of 150 days), set to 0.00667, else set to 1/5 the value
       mu.matrix[j,i]<-0.00667,
       mu.matrix[j,i]<-temp.value[i]/5)
}
}



# For parasite background mortality rate, assuming mortality follows an inverted quadratic

theta.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

# pull 1000 different tmin, c, and add somewhere between 10 and 25 to get tmax
theta.inter <- runif(dim(theta.matrix)[1], min=.1,max=1) # this parameter also known as c
theta.n.slope <- runif(dim(theta.matrix)[1], min=.02,max=0.03) # this parameter also known as b
theta.qd <-  runif(dim(theta.matrix)[1], min=0.0008, max=0.0009) #this parameter also known as a


for(j in 1:dim(theta.matrix)[1]){
  for(i in 1:length(temp)){
    
    temp.value <- c()
    temp.value[i] <- (theta.qd[j]*(temp[i])^2)-(theta.n.slope[j]*temp[i])+theta.inter[j]
    
    ifelse(temp.value/5<0.0667, #if 1/5 value is less than 0.005 (corresponding to a lifespan of 200 days), set to 0.005, else set to 1/5 the value
           theta.matrix[j,i]<-0.0667,
           theta.matrix[j,i]<-temp.value[i]/5)
  }
}




# approximating density as birth rate / death rate


density.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

density.matrix[,] <- birth.rate.matrix[,]/(mu.matrix[,])




# R0 calculations for different scenarios #

R0.matrix <- matrix(nrow=dim(ind.matrix)[1],ncol=length(temp))

for(j in 1:dim(R0.matrix)[1]){
  for(i in 1:length(temp)){
    
    R0.matrix[j,i] = (chi.matrix[j,i]*sigma.matrix[j,i]*density.matrix[j,i] / theta.matrix[j,i]) * 
  ((lambda.matrix[j,i] / mu.matrix[j,i] + alpha.matrix[j,i]) + omega.matrix[j,i])

  }
}


R0.topt <- c()
ind.topt <- c()
vec <- c() # 0 and 1s for if the epidemic spreads (If R0 = 0 or not)


for(j in 1:dim(R0.matrix)[1]){
  R0.topt[j] <- temp[which.max((R0.matrix[j,]))]
  ind.topt[j] <- temp[which.max((ind.matrix[j,]))]
  vec[j] <- isTRUE(R0.topt[j]>0)
}


data.mat <- data.frame(ind.topt,R0.topt,vec)


# filter out ones where R0 was 0 
data.mat <- data.mat %>% 
  filter(vec==TRUE)

# Save the thermal optima data in a CSV to plot after #
write.csv(data.mat, file="~/enviro_simulations_only_omega_prop.csv")

