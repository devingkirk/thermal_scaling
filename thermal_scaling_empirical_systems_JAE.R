# Load packages #
require(scales)
require(dplyr)
require(wesanderson)
require(plotrix)

# Set working directory #

# Read in the data #
data <- read.csv("~/empirical_thermal_optima_data_JAE.csv")

head(data)

# Data has 24 systems: 13 VBD systems, 11 environmentally-transmitted systems #

# For each mosquito system, we have three different measures of individual-level parasitism: infected days, vector competence, parasite development rate
# We will use infected days for the overall analysis, though the data also has the other mosquito metrics to compare #

# systems is the one measure for each system #
systems <- data[1:24,]

# separate out vector-borne systems versus environmentally-transmitted #
vector.systems <- systems %>%
  filter(Transmission.type == "vector")

enviro.systems <- systems %>%
  filter(Transmission.type == "environment")


# Test for correlations in each group #
vector.sys.cor <- cor.test(x=vector.systems$ind_parasitism_topt, y=vector.systems$pop_parasitism_topt, method="pearson"); vector.sys.cor

# For environmental systems, first test for correlation with all 11 systems #
full.enviro.sys.cor <- cor.test(x=enviro.systems$ind_parasitism_topt, y=enviro.systems$pop_parasitism_topt, method="pearson"); full.enviro.sys.cor

# and then test for correlation with the 4 environmental systems where both Topt peak at intermediate temps, not the end of the range examined 
enviro.systems.subset <- enviro.systems %>%
  filter(exclude_for_TPC == FALSE)

subset.enviro.sys.cor <- cor.test(x=enviro.systems.subset$ind_parasitism_topt, y=enviro.systems.subset$pop_parasitism_topt, method="pearson"); subset.enviro.sys.cor



# Now test two further subsets #

# Looking for correlation in the systems where trait-based models predict population-level parasitism
# i.e. the two measures are not independent but temperature effects are isolated #

non.ind.systems <- systems %>%
  filter(measures_nonindependent == TRUE)

non.ind.cor <- cor.test(x=non.ind.systems$ind_parasitism_topt, y=non.ind.systems$pop_parasitism_topt, method="pearson"); non.ind.cor



# Looking for correlation in the systems where individual- and population-level measures  are obtained separately 
# i.e. the two measures are independent but the data is not necessarily isolating temperature effects #

ind.systems <- systems %>%
  filter(measures_nonindependent == FALSE)

ind.cor <- cor.test(x=ind.systems$ind_parasitism_topt, y=ind.systems$pop_parasitism_topt, method="pearson"); ind.cor



# Overall correlation in 24 systems:

overall.sys.cor <- cor.test(x=systems$ind_parasitism_topt, y=systems$pop_parasitism_topt, method="pearson"); overall.sys.cor





## Plot figures ##
# One color for VBD, one color for environmental with all intermediate TPC peaks, one color for env. with peak at end of range
  
col.pal <- wes_palette("FantasticFox1",n=3,type="discrete")


col.vec <- c(rep(col.pal[1],13),
             col.pal[2],col.pal[2],col.pal[2],col.pal[3],     
             col.pal[2],  col.pal[3],  col.pal[3],  col.pal[3],  
             col.pal[2],  col.pal[2],  col.pal[2])
             


png("~/thermal_scaling_JAE_fig_1.png", width=3000, height=4000, res=400, family="Arial")
par(mfrow=c(2,1))
par(oma=c(1,1,1,1))
par(mar=c(5,4,2,2))

plot(vector.systems$ind_parasitism_topt, vector.systems$pop_parasitism_topt, ylim=c(18,30), xlim=c(18,30), xlab=NA,ylab=NA,type="n")
abline(a=0,b=1,lty=2, lwd=1.5)
points(vector.systems$ind_parasitism_topt, vector.systems$pop_parasitism_topt, bg=alpha(col.vec[1:13],0.9), pch=23,cex=2)
mtext(expression(paste("Population-level parasitism T"[opt]," ("~degree~"C)")),side=2,las=0.5,line=3,outer=FALSE,cex=1.25)
mtext(expression(paste("Individual-level parasitism T"[opt]," ("~degree~"C)")),side=1,las=1,line=3,outer=FALSE,cex=1.25)
text(x=18,y=29.5,labels="a)", font=1, cex=1.5)
text(x=20,y=29.5,labels="r = 0.82; p < 0.001", font=1)
legend(x=18,y=29.5,legend=c("1:1 line"), lty=c(2), bty="n")
title(main="Empirical vector-borne systems")

plot(enviro.systems$ind_parasitism_topt, enviro.systems$pop_parasitism_topt, ylim=c(10,27), xlim=c(10,27), xlab=NA,ylab=NA,type="n")
abline(a=0,b=1,lty=2, lwd=1.5)
points(enviro.systems$ind_parasitism_topt, enviro.systems$pop_parasitism_topt, bg=alpha(col.vec[14:24],0.9), pch=21,cex=2)
mtext(expression(paste("Population-level parasitism T"[opt]," ("~degree~"C)")),side=2,las=0.5,line=3,outer=FALSE,cex=1.25)
mtext(expression(paste("Individual-level parasitism T"[opt]," ("~degree~"C)")),side=1,las=1,line=3,outer=FALSE,cex=1.25)
text(x=10,y=26.5,labels="b)",font=1, cex=1.5)
text(x=14.6,y=26.5,labels="r = 0.74 with all systems; p = 0.01 ", font=1)
text(x=14,y=25.25,labels="r = 0.65 with subset; p = 0.35", font=1)
title(main="Empirical environmentally-transmitted systems")
legend(x=10.15,y=25.2,legend=c(expression(paste("T"[opt]," intermediate")),expression(paste("T"[opt]," at end of examined range")),
                             "1:1 line"),
       pt.bg=c(alpha(col.pal[3],0.9),alpha(col.pal[2]),0.9,NA), pch=c(21,21,NA), lty=c(NA,NA,2), bty="n")

dev.off()




## Figure 2
systems.2 <- systems %>%
  filter(!is.na(host_topt_ind))

col.vec.2 <- c(rep(col.pal[1],9),
             col.pal[2],col.pal[2],col.pal[3],col.pal[3],     
             col.pal[2],  col.pal[3],  col.pal[3],  col.pal[2])


sys.pch <- c(rep(23,9),rep(21,8))


# Figure 2 #
png("~/thermal_scaling_JAE_fig2.png", width=2400, height=2400, res=400, family="Arial")
par(mfrow=c(1,1))
par(oma=c(1,1,1,1))
par(mar=c(4,4,2,2))
plot(systems.2$ind_parasitism_rel_host, systems.2$popparasitism_rel_host, ylim=c(-15,15), xlim=c(-15,15), cex=1.25,lwd=1, type="n",xlab=NA,ylab=NA)
abline(a=0,b=0,lty=1, col="black")
abline(v=0,b=0,lty=1,col="black")
abline(a=0,b=1,lty=2)

draw.circle(x=0,y=0,radius=2.5,lwd=1, border="black")

points(systems.2$ind_parasitism_rel_host, systems.2$pop_parasitism_rel_host,  cex=1.5,pch=sys.pch,lwd=1.2,
       bg=alpha(col.vec.2,0.9),col="black")
mtext(expression(paste("Individual-level parasitism T"[opt]," - host performance T"[opt]," ("~degree~"C)")),side=1,las=1,line=3,outer=FALSE,cex=1)
mtext(expression(paste("Population-level parasitism T"[opt]," - host performance T"[opt]," ("~degree~"C)")),side=2,las=0.5,line=3,outer=FALSE,cex=1)


dev.off()

