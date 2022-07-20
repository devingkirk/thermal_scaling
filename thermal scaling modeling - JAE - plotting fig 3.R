
# This script is to plot model simulation data #
require(scales)


# Set working directory that contains csv files of the four different simulations #

# Bring in CSV files of thermal optima data from simulations  #
mosq.data.mat <- read.csv("~/Desktop/mosquito_simulations.csv")
enviro.data.both.mat <- read.csv("~/Desktop/enviro_simulations_both_prop.csv")
enviro.data.omega.mat <- read.csv("~/Desktop/enviro_simulations_only_omega_prop.csv")
enviro.data.neither.mat <- read.csv("~/Desktop/enviro_simulations_neither_proportional.csv")

mosq.cor <- cor.test(x=mosq.data.mat$ind.topt, y=mosq.data.mat$R0.topt, method="pearson"); mosq.cor
enviro.both.cor <- cor.test(x=enviro.data.both.mat$ind.topt, y=enviro.data.both.mat$R0.topt, method="pearson"); enviro.both.cor
enviro.omega.cor <- cor.test(x=enviro.data.omega.mat$ind.topt, y=enviro.data.omega.mat$R0.topt, method="pearson"); enviro.omega.cor
enviro.neither.cor <- cor.test(x=enviro.data.neither.mat$ind.topt, y=enviro.data.neither.mat$R0.topt, method="pearson"); enviro.neither.cor



png("~/Desktop/thermal_scaling_JAE_fig_3.png", width=3960, height=3960, res=400, family="Arial")
par(mfrow=c(2,2))
par(oma=c(1,1,1,1))
par(mar=c(5,4,4,2))

plot(mosq.data.mat$ind.topt, mosq.data.mat$R0.topt, ylim=c(7,28), xlim=c(7,28), xlab=NA,ylab=NA,type="n")
abline(a=0,b=1,lty=2, lwd=1.5)
points(mosq.data.mat$ind.topt, mosq.data.mat$R0.topt,pch=19,lwd=1.2,cex=1,
       bg=alpha("black",0.5), col=alpha("black",0.5))
mtext(expression(paste("Population-level parasitism T"[opt]," ("~degree~"C)")),side=2,las=0.5,line=3,outer=FALSE,cex=1.25)
mtext(expression(paste("Individual-level parasitism T"[opt]," ("~degree~"C)")),side=1,las=1,line=3,outer=FALSE,cex=1.25)
text(x=10,y=27,labels="a) r = 0.68", font=1, cex=1.5)
title(font.main=1,cex.main=1.1,main="Simulated vector-borne systems")


plot(enviro.data.both.mat$ind.topt, enviro.data.both.mat$R0.topt, ylim=c(7,28), xlim=c(7,28), xlab=NA,ylab=NA,type="n")
abline(a=0,b=1,lty=2, lwd=1.5)
points(enviro.data.both.mat$ind.topt, enviro.data.both.mat$R0.topt,pch=19,lwd=1.2,cex=1,
       bg=alpha("black",0.5), col=alpha("black",0.5))
mtext(expression(paste("Population-level parasitism T"[opt]," ("~degree~"C)")),side=2,las=0.5,line=3,outer=FALSE,cex=1.25)
mtext(expression(paste("Individual-level parasitism T"[opt]," ("~degree~"C)")),side=1,las=1,line=3,outer=FALSE,cex=1.25)
text(x=10,y=27,labels="b) r = 0.41", font=1, cex=1.5)
title(font.main=1,cex.main=1.1,main="Simulated environmentally-transmitted systems: \n two parasite-shedding parameters scale with parasite load")

plot(enviro.data.omega.mat$ind.topt, enviro.data.omega.mat$R0.topt, ylim=c(7,28), xlim=c(7,28), xlab=NA,ylab=NA,type="n")
abline(a=0,b=1,lty=2, lwd=1.5)
points(enviro.data.omega.mat$ind.topt, enviro.data.omega.mat$R0.topt,pch=19,lwd=1.2,cex=1,
       bg=alpha("black",0.5), col=alpha("black",0.5))
mtext(expression(paste("Population-level parasitism T"[opt]," ("~degree~"C)")),side=2,las=0.5,line=3,outer=FALSE,cex=1.25)
mtext(expression(paste("Individual-level parasitism T"[opt]," ("~degree~"C)")),side=1,las=1,line=3,outer=FALSE,cex=1.25)
text(x=10,y=27,labels="c) r = 0.26", font=1, cex=1.5)
title(font.main=1,cex.main=1.1,main="Simulated environmentally-transmitted systems: \n one parasite-shedding parameter scales with parasite load")

plot(enviro.data.neither.mat$ind.topt, enviro.data.neither.mat$R0.topt, ylim=c(7,28), xlim=c(7,28), xlab=NA,ylab=NA,type="n")
abline(a=0,b=1,lty=2, lwd=1.5)
points(enviro.data.neither.mat$ind.topt, enviro.data.neither.mat$R0.topt,pch=19,lwd=1.2,cex=1,
       bg=alpha("black",0.5), col=alpha("black",0.5))
mtext(expression(paste("Population-level parasitism T"[opt]," ("~degree~"C)")),side=2,las=0.5,line=3,outer=FALSE,cex=1.25)
mtext(expression(paste("Individual-level parasitism T"[opt]," ("~degree~"C)")),side=1,las=1,line=3,outer=FALSE,cex=1.25)
text(x=10,y=27,labels="d) r = 0.00", font=1, cex=1.5)
title(font.main=1,cex.main=1.1,main="Simulated environmentally-transmitted systems: \n no parasite-shedding parameters scale with parasite load")

dev.off()

