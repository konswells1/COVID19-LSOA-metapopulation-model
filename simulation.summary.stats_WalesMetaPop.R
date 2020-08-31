
Drive <- "F"

dir.library <- paste(Drive, ":/Kons_backup/R Library", sep = "")
.libPaths(dir.library)

dir.project <- paste(Drive, ":/Kons_backup/2020.06_Covid19_MetaPopWales/", sep = "")
dir_CoV.data <- paste(Drive,  ":/Kons_backup/CoV_Data/", sep = "")
dir_UK.geo <- paste(Drive,  ":/Kons_backup/CoV_Data/Data_UK.Geography/", sep = "")

## Folder name
SIM = "sims_WalesMetaPop_200811/MetaWales"
dir.sim.cloud <- paste(dir.project, "/", SIM, sep="")
dir.analysis <- paste(dir.project, "/sims_WalesMetaPop_200811/", sep="")

#######################################
##
##   COVID-19 M4/Wales  metapopulation model   
##
#######################################


library(raster)
library(data.table)
library(lhs)
library(sf)
library(sp)
library(rgeos)
library(dismo)
library(ggplot2)
library(GGally)
library("ggpubr")
library(ggpubr)
library(beeswarm)
library(RColorBrewer)  
library(colorspace)

crs.WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
clgr.1 <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)
clgr.2 <- colorRampPalette(rev(brewer.pal(11, 'BrBG')))(100)
cl.1 <- brewer.pal(3, name = "Set1")
cl.6 <- brewer.pal(6, name = "Dark2")
cl.4 <- brewer.pal(4, name = "Set2")

setwd(dir.analysis)
load("WSP_MetaWales.RData")

dir.sim = dir.sim.cloud
setwd(dir.sim)

# Vectors of parameter values from hypercube
scenario <- Hypercube$scenario
beta <- Hypercube$beta
p.asympt <- Hypercube$p.asympt
inf.asympt <- Hypercube$inf.asympt
disp.rate <- Hypercube$disp.rate
p.resist <- Hypercube$p.resist
delta <- Hypercube$delta
gravD <- Hypercube$gavD
p.trace.EAI <- Hypercube$p.trace.EAI
p.trace.I <- Hypercube$p.trace.I
lockdown.treshold <- Hypercube$lockdown.treshold
lockdown.time <- Hypercube$lockdown.time
lockdown.dist <- Hypercube$lockdown.dist
lockdown.factor <- 1 -Hypercube$lockdown.factor    # Reverse direction of modelled factor for easier interepretation

out_names <- c("nS",  "nE",  "nA", "nI",  "nR", "A.new", "I.new", "time.Lockdown")

sample.demogr <- sort(rep(1:NSample_demogr, 4))

parameter_names <- names(Hypercube)[-1]
parameter_names[which(parameter_names=="gavD")] <- "gravD"

NPar <- length(parameter_names)

## Check file existence
exist.out <- exist.gravi <- rep(NA, NSim)
for(x in 1:NSim){
	file.i <- paste("Out_ind_", x, ".RData", sep="") 
	if(file.exists(file.i)){
		exist.out[x] <- 1
	} 	
	file.j <- paste("pop.grav.index_", x, ".RData", sep="") 
	if(file.exists(file.j)){
		exist.gravi[x] <- 1
	} 
}
table(exist.out); table(exist.gravi)

######
## Load and summarize output into summary objects

####
## Population-level number of I and A ('epidemic size': total number of individuals over modelled time period)
setwd(dir.sim)
popOut_nI <- array(NA, dim=c(NSim, NPop))

for(x in 1:NSim){
	file.i <- paste("Out_ind_", x, ".RData", sep="") 
	if(file.exists(file.i)){
 		load(file.i)
		for(m in 1:NPop){
			popOut_nI[x,m] <- sum(Out_ind[m,, which(out_names=="I.new")], na.rm=T)
		}
	}
}


Mat_pop.grav.index <- matrix(NA, nrow=NSample_demogr,  ncol=NPop)
for(x in 1:NSample_demogr){
	simno <- which(scenario==1)[x]
	file.j <- paste("pop.grav.index_", simno, ".RData", sep="") 
	if(file.exists(file.j)){
		load(file.j) 
		Mat_pop.grav.index[x,] <- pop.grav.index
	}
}


# Epidemic size for each simulation
sim_episize <- rep(NA, NSim)
for(x in 1:NSim){
	sim_episize[x] <- sum(popOut_nI[x,])
}


## Relative epidemic size for each management scenario in relation to control
local_episize_basel <- 
local_rel.episize_traceEAI <- 
local_rel.episize_traceI <- 
local_rel.episize_lockd <- array(NA, dim=c(NSample_demogr, NPop))

metapop_rel.episize_traceEAI <- 
metapop_rel.episize_traceI <- 
metapop_rel.episize_lockd <- rep(NA, NSample_demogr)

metapop_var.episize_basel <- 
metapop_var.episize_traceEAI <- 
metapop_var.episize_traceI <- 
metapop_var.episize_lockd <- rep(NA, NSample_demogr)

for(x in 1:NSample_demogr){
	local_episize_basel[x,] <- popOut_nI[which(sample.demogr==x&scenario==1),]
	local_rel.episize_traceEAI[x,] <- round(popOut_nI[which(sample.demogr==x&scenario==2),]/(popOut_nI[which(sample.demogr==x&scenario==1),]+1),3)
	local_rel.episize_traceI[x,] <- round(popOut_nI[which(sample.demogr==x&scenario==3),]/(popOut_nI[which(sample.demogr==x&scenario==1),]+1),3)
	local_rel.episize_lockd[x,] <- round(popOut_nI[which(sample.demogr==x&scenario==4),]/(popOut_nI[which(sample.demogr==x&scenario==1),]+1),3)

	metapop_rel.episize_traceEAI[x] <- round(sum(popOut_nI[which(sample.demogr==x&scenario==2),])/sum(popOut_nI[which(sample.demogr==x&scenario==1),]+1),3)
	metapop_rel.episize_traceI[x] <- round(sum(popOut_nI[which(sample.demogr==x&scenario==3),])/sum(popOut_nI[which(sample.demogr==x&scenario==1),]+1),3)
	metapop_rel.episize_lockd[x] <- round(sum(popOut_nI[which(sample.demogr==x&scenario==4),])/sum(popOut_nI[which(sample.demogr==x&scenario==1),]+1),3)

	if(!is.na(popOut_nI[which(sample.demogr==x),1])){
	metapop_var.episize_basel[x] <- var(popOut_nI[which(sample.demogr==x&scenario==1),]/ max(popOut_nI[which(sample.demogr==x&scenario==1),], na.rm=T))
	metapop_var.episize_traceEAI[x] <- var(popOut_nI[which(sample.demogr==x&scenario==2),]/ max(popOut_nI[which(sample.demogr==x&scenario==2),], na.rm=T))
	metapop_var.episize_traceI[x] <- var(popOut_nI[which(sample.demogr==x&scenario==3),]/ max(popOut_nI[which(sample.demogr==x&scenario==3),], na.rm=T))
	metapop_var.episize_lockd[x] <- var(popOut_nI[which(sample.demogr==x&scenario==4),]/ max(popOut_nI[which(sample.demogr==x&scenario==4),], na.rm=T))
	}else{}
}


## "Urban-rural gradient": Correlation between population-level relative epidemic size with managements and population gravity/density
# Positive correlation: increase in relative epidemic size with increased connectivity/density -> less efficient control in 'urban area' 

metapop_gradient_basel <- rep(NA, NSample_demogr)
metapop_gradient_traceEAI <- rep(NA, NSample_demogr)
metapop_gradient_traceI <- rep(NA, NSample_demogr)
metapop_gradient_lockd <- rep(NA, NSample_demogr)

for(x in 1:NSample_demogr){
	# Compute gravity distance matrix and population-level gravity index for sample=specific parameter gravD
	gravD <- Hypercube$gavD  [which(Hypercube$scenario==1)] [x]

	if(!is.na(local_rel.episize_traceEAI[x,1]) & !is.na(Mat_pop.grav.index[x,1])){
	metapop_gradient_basel[x] <- cor.test(Mat_pop.grav.index[x,], local_episize_basel[x,], method='spearman')$estimate
	metapop_gradient_traceEAI[x] <- cor.test(Mat_pop.grav.index[x,], local_rel.episize_traceEAI[x,], method='spearman')$estimate
	}else{}
	if(!is.na(local_rel.episize_traceI[x,1])& !is.na(Mat_pop.grav.index[x,1])){
	metapop_gradient_traceI[x] <- cor.test(Mat_pop.grav.index[x,], local_rel.episize_traceI[x,], method='spearman')$estimate
	}else{}
	if(!is.na(local_rel.episize_lockd[x,1])& !is.na(Mat_pop.grav.index[x,1])){
	metapop_gradient_lockd[x] <- cor.test(Mat_pop.grav.index[x,], local_rel.episize_lockd[x,], method='spearman')$estimate
	}else{}
}


# Dataframe of parameters and output
out_sim <- data.frame(
	sim_episize = sim_episize,
	Hypercube
) 

  Out_metapop <- data.frame(
	rel.episize_traceEAI = metapop_rel.episize_traceEAI,
	rel.episize_traceI = metapop_rel.episize_traceI,
	rel.episize_lockd = metapop_rel.episize_lockd,

	gradient_basel = metapop_gradient_basel,
	gradient_traceEAI = metapop_gradient_traceEAI,
	gradient_traceI = metapop_gradient_traceI,
	gradient_lockd = metapop_gradient_lockd,

	beta =Hypercube$beta [which(Hypercube$scenario==1)],
	p.asympt = Hypercube$p.asympt [which(Hypercube$scenario==1)],
	inf.asympt = Hypercube$inf.asympt [which(Hypercube$scenario==1)],
	disp.rate = Hypercube$disp.rate [which(Hypercube$scenario==1)],
	p.resist = Hypercube$p.resist [which(Hypercube$scenario==1)],
	delta = Hypercube$delta [which(Hypercube$scenario==1)],
	gravD = Hypercube$gavD  [which(Hypercube$scenario==1)],
	p.trace.EAI = Hypercube$p.trace.EAI [which(Hypercube$scenario==2)],
	p.trace.I = Hypercube$p.trace.I [which(Hypercube$scenario==3)],
	lockdown.treshold = Hypercube$lockdown.treshold [which(Hypercube$scenario==4)],
	lockdown.time = Hypercube$lockdown.time [which(Hypercube$scenario==4)],
	lockdown.dist = Hypercube$lockdown.dist [which(Hypercube$scenario==4)],
	lockdown.factor = 1 - Hypercube$lockdown.factor [which(Hypercube$scenario==4)]  # Reverse direction of modelled factor for easier interepretation
)
Out_metapop <- Out_metapop[complete.cases(Out_metapop),]


###################
# Analysis: sensitivity analysis of relative epidemic size of baseline scenarios

form_episize.basel <- gradient_basel ~ scale(beta) + scale(p.asympt) +  scale(inf.asympt) + scale(disp.rate) + scale(p.resist) + scale(delta) + scale(gravD)
glm_episize.basel <- glm(form_episize.basel, data= Out_metapop)
summary(glm_episize.basel)

# Run BRT with cross-validation method using gbm.step()
gbm_episize.basel <- gbm.step(data=Out_metapop, gbm.x = match(c("beta", "p.asympt", "inf.asympt", "disp.rate", "p.resist", "delta", "gravD"), names(Out_metapop)), 
		gbm.y = which( names(Out_metapop)=="gradient_basel"), family = "gaussian", n.trees = 10000, tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)
summary(gbm_episize.basel)
gbm_episize.basel$contributions

###################
# Analysis: sensitivity analysis of relative epidemic size with different control strategies

form_episize_traceEAI <- rel.episize_traceEAI ~ scale(beta) + scale(p.asympt) +  scale(inf.asympt) + scale(disp.rate) + scale(p.resist) + scale(delta) + scale(gravD) + scale(p.trace.EAI)
glm_episize_traceEAI <- glm(form_episize_traceEAI, data= Out_metapop)
summary(glm_episize_traceEAI)

form_episize_traceI <- rel.episize_traceI ~ scale(beta) + scale(p.asympt) +  scale(inf.asympt) + scale(disp.rate) + scale(p.resist) + scale(delta) + scale(gravD) + scale(p.trace.I) 
glm_episize_traceI <- glm(form_episize_traceI, data= Out_metapop)
summary(glm_episize_traceI)

form_episize_lockd <- rel.episize_lockd ~ scale(beta) + scale(p.asympt) +  scale(inf.asympt) + scale(disp.rate) + scale(p.resist) + scale(delta) + scale(gravD) + scale(lockdown.treshold) + scale(lockdown.time) + scale(lockdown.dist) + scale(lockdown.factor)
glm_episize_lockd <- glm(form_episize_lockd, data= Out_metapop)
summary(glm_episize_lockd)


# Run BRT with cross-validation method using gbm.step()
gbm_episize_traceEAI <- gbm.step(data=Out_metapop, gbm.x = match(c("beta", "p.asympt", "inf.asympt", "disp.rate", "p.resist", "delta", "gravD", "p.trace.EAI"), names(Out_metapop)), 
		gbm.y = which( names(Out_metapop)=="rel.episize_traceEAI"), family = "gaussian", n.trees = 10000, tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)
summary(gbm_episize_traceEAI)
gbm_episize_traceEAI$contributions

gbm_episize_traceI <- gbm.step(data=Out_metapop, gbm.x = match(c("beta", "p.asympt", "inf.asympt", "disp.rate", "p.resist", "delta", "gravD", "p.trace.I"), names(Out_metapop)), 
		gbm.y = which( names(Out_metapop)=="rel.episize_traceI"), family = "gaussian", n.trees = 10000, tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)
summary(gbm_episize_traceI)
gbm_episize_traceI$contributions

gbm_episize_lockd <- gbm.step(data=Out_metapop, gbm.x = match(c("beta", "p.asympt", "inf.asympt", "disp.rate", "p.resist", "delta", "gravD", "lockdown.treshold", "lockdown.time", "lockdown.dist", "lockdown.factor"), names(Out_metapop)), 
		gbm.y = which( names(Out_metapop)=="rel.episize_lockd"), family = "gaussian", n.trees = 10000, tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)
summary(gbm_episize_lockd)
gbm_episize_lockd$contributions

###################
# Analysis: sensitivity analysis the correlation in population centrality and relative epidemic size across urban-rural gradients

form_gradient_traceEAI <- gradient_traceEAI ~ scale(beta) + scale(p.asympt) +  scale(inf.asympt) + scale(disp.rate) + scale(p.resist) + scale(delta) + scale(gravD) + scale(p.trace.EAI)
glm_gradient_traceEAI <- glm(form_gradient_traceEAI, data= Out_metapop)
summary(glm_gradient_traceEAI)

form_gradient_traceI <- gradient_traceI ~ scale(beta) + scale(p.asympt) +  scale(inf.asympt) + scale(disp.rate) + scale(p.resist) + scale(delta) + scale(gravD) + scale(p.trace.I)
glm_gradient_traceI <- glm(form_gradient_traceI, data= Out_metapop)
summary(glm_gradient_traceI)

form_gradient_lockd <- gradient_lockd ~ scale(beta) + scale(p.asympt) +  scale(inf.asympt) + scale(disp.rate) + scale(p.resist) + scale(delta) + scale(gravD) + scale(lockdown.treshold) + scale(lockdown.time) + scale(lockdown.dist) + scale(lockdown.factor)
glm_gradient_lockd <- glm(form_gradient_lockd, data= Out_metapop)
summary(glm_gradient_lockd)


gbm_gradient_traceEAI <- gbm.step(data=Out_metapop, gbm.x = match(c("beta", "p.asympt", "inf.asympt", "disp.rate", "p.resist", "delta", "gravD", "p.trace.EAI"), names(Out_metapop)), 
		gbm.y = which( names(Out_metapop)=="gradient_traceEAI"), family = "gaussian", n.trees = 10000, tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)
summary(gbm_gradient_traceEAI)
gbm_gradient_traceEAI$contributions


gbm_gradient_traceI <- gbm.step(data=Out_metapop, gbm.x = match(c("beta", "p.asympt", "inf.asympt", "disp.rate", "p.resist", "delta", "gravD", "p.trace.I"), names(Out_metapop)), 
		gbm.y = which( names(Out_metapop)=="gradient_traceI"), family = "gaussian", n.trees = 10000, tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)
summary(gbm_gradient_traceI)
gbm_gradient_traceI$contributions

gbm_gradient_lockd <- gbm.step(data=Out_metapop, gbm.x = match(c("beta", "p.asympt", "inf.asympt", "disp.rate", "p.resist", "delta", "gravD", "lockdown.treshold", "lockdown.time", "lockdown.dist", "lockdown.factor"), names(Out_metapop)), 
		gbm.y = which( names(Out_metapop)=="gradient_lockd"), family = "gaussian", n.trees = 10000, tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)
summary(gbm_gradient_lockd)
gbm_gradient_lockd$contributions

##############
# Table of parameter relative importances

predictor.all <- paste(c(paste("scale(", parameter_names[-NPar], ")+", sep=""), paste("scale(", parameter_names[NPar], ")", sep="")), collapse='')
predictor.all <- gsub("gavD", "gravD", predictor.all)


Dat.summary <- data.frame(Par=parameter_names)
Dat.summary$RelInfl_EpiSize_traceEAI <- gbm_episize_traceEAI$contributions$rel.inf[match(Dat.summary$Par, rownames(gbm_episize_traceEAI$contributions))]
Dat.summary$RelInfl_EpiSize_traceI <- gbm_episize_traceI$contributions$rel.inf[match(Dat.summary$Par, rownames(gbm_episize_traceI$contributions))]
Dat.summary$RelInfl_EpiSize_lockd <- gbm_episize_lockd$contributions$rel.inf[match(Dat.summary$Par, rownames(gbm_episize_lockd$contributions))]
Dat.summary$RelInfl_gradient_traceEAI <- gbm_gradient_traceEAI$contributions$rel.inf[match(Dat.summary$Par, rownames(gbm_gradient_traceEAI$contributions))]
Dat.summary$RelInfl_gradient_traceI <- gbm_gradient_traceI$contributions$rel.inf[match(Dat.summary$Par, rownames(gbm_gradient_traceI$contributions))]
Dat.summary$RelInfl_gradient_lockd <- gbm_gradient_lockd$contributions$rel.inf[match(Dat.summary$Par, rownames(gbm_gradient_lockd$contributions))]

Dat.summary$hjust_EpiSize_traceEAI <- ifelse(glm(as.formula(paste("rel.episize_traceEAI ~ ", predictor.all)), data= Out_metapop)$coefficients[-1]>0, 1, -1)
Dat.summary$hjust_EpiSize_traceI <- ifelse(glm(as.formula(paste("rel.episize_traceI ~ ", predictor.all)), data= Out_metapop)$coefficients[-1]>0, 1, -1)
Dat.summary$hjust_EpiSize_lockd <- ifelse(glm(as.formula(paste("rel.episize_lockd ~ ", predictor.all)), data= Out_metapop)$coefficients[-1]>0, 1, -1)
Dat.summary$hjust_gradient_traceEAI <- ifelse(glm(as.formula(paste("gradient_traceEAI ~ ", predictor.all)), data= Out_metapop)$coefficients[-1]>0, 1, -1)
Dat.summary$hjust_gradient_traceI <- ifelse(glm(as.formula(paste("gradient_traceI ~ ", predictor.all)), data= Out_metapop)$coefficients[-1]>0, 1, -1)
Dat.summary$hjust_gradient_lockd <- ifelse(glm(as.formula(paste("gradient_lockd ~ ", predictor.all)), data= Out_metapop)$coefficients[-1]>0, 1, -1)
Dat.summary[is.na(Dat.summary)] <- 0

Dat.EpiSize.long <- data.frame(Par=rep(names(Hypercube)[-1],3))
Dat.EpiSize.long$Type <- c(rep("3) Trace all", NPar), rep("1) Trace symptomatic only", NPar), rep("2) Regional lockdown", NPar))
Dat.EpiSize.long$RelInfl <- c(Dat.summary$RelInfl_EpiSize_traceEAI, Dat.summary$RelInfl_EpiSize_traceI, Dat.summary$RelInfl_EpiSize_lockd)
Dat.EpiSize.long$hjust <- c(Dat.summary$hjust_EpiSize_traceEAI, Dat.summary$hjust_EpiSize_traceI, Dat.summary$hjust_EpiSize_lockd)

Dat.Gradient.long <- data.frame(Par=rep(names(Hypercube)[-1],3))
Dat.Gradient.long$Type <- c(rep("3) Trace all", NPar), rep("1) Trace symptomatic only", NPar), rep("2) Regional lockdown", NPar))
Dat.Gradient.long$RelInfl <- c(Dat.summary$RelInfl_gradient_traceEAI, Dat.summary$RelInfl_gradient_traceI, Dat.summary$RelInfl_gradient_lockd)
Dat.Gradient.long$hjust <- c(Dat.summary$hjust_gradient_traceEAI, Dat.summary$hjust_gradient_traceI, Dat.summary$hjust_gradient_lockd)

# Parameter names as used in  plot
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="beta")] <- "Transmission rate"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="p.asympt")] <- "Asymptomatics proportion"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="inf.asympt")] <- "Asymptomatics infectiousness"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="disp.rate")] <- "Travel frequency"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="p.resist")] <- "Initally recovered"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="delta")] <- "Density dependence"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="gavD")] <- "Distance weight in gravity model"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="p.trace.EAI")] <- "Isolation all infected ind."
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="p.trace.I")] <- "Isolation symptomatic ind."
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="lockdown.treshold")] <- "Lockdown threshold"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="lockdown.time")] <- "Lockdown duration"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="lockdown.dist")] <- "Lockdown travel limit"
Dat.EpiSize.long$Par[which(Dat.EpiSize.long$Par=="lockdown.factor")] <- "Lockdown stringency"

Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="beta")] <- "Transmission rate"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="p.asympt")] <- "Asymptomatics proportion"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="inf.asympt")] <- "Asymptomatics infectiousness"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="disp.rate")] <- "Travel frequency"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="p.resist")] <- "Initally recovered"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="delta")] <- "Density dependence"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="gavD")] <- "Distance weight in gravity model"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="p.trace.EAI")] <- "Isolation all infected ind."
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="p.trace.I")] <- "Isolation symptomatic ind."
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="lockdown.treshold")] <- "Lockdown threshold"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="lockdown.time")] <- "Lockdown duration"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="lockdown.dist")] <- "Lockdown travel limit"
Dat.Gradient.long$Par[which(Dat.Gradient.long$Par=="lockdown.factor")] <- "Lockdown stringency"


###############
# Plots

plot_EpiSize.RelInfl <- ggplot(Dat.EpiSize.long, aes(x=Par, y=RelInfl))+
	coord_flip() + 
	geom_bar(stat='identity', show.legend = F, aes(fill=ifelse(hjust==1,  "#377EB8", "#E41A1C")))+
	facet_wrap(~Type) +
	ggtitle("Relative epidemic size") +
	xlab("Parameter") +
	ylab("Relative influence [%]") +
	theme(text = element_text(size=10),
          title = element_text(size=12),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=15),
	    axis.text.y.left = element_text(size=12))

plot_EpiSize.RelInfl


plot_Gradient.RelInfl <- ggplot(Dat.Gradient.long, aes(x=Par, y=RelInfl))+
	coord_flip() + 
	geom_bar(stat='identity', show.legend = FALSE, aes(fill=ifelse(hjust==1,  "red", "blue"))) +
	facet_wrap(~Type) +
	ggtitle("Urban-rural epidemic gradient (correlation)") +
	xlab("Parameter") +
	ylab("Relative influence [%]" ) +
	theme(text = element_text(size=10),
          title = element_text(size=12),
          strip.text.x = element_text(size=12),
          strip.text.y = element_text(size=15),
	    axis.text.y.left = element_text(size=12))

plot_Gradient.RelInfl 

plot_beeswarm <- function(){ 
	par(mfrow=c(1,1), oma=c(4,4,1,1))
	cl.4 <- brewer.pal(4, name = "Set2")
	group.order <- ifelse(scenario==1, "g.1", "g.x")
	group.order <- ifelse(scenario==2, "g.4", group.order)
	group.order <- ifelse(scenario==3, "g.2", group.order)
	group.order <- ifelse(scenario==4, "g.3", group.order)
	group.labels <- c("Baseline", "Trace symptomatic only", "2) Local lockdown", "3) Trace all")
	group.labels <- c("Baseline", "Trace sympt.", "Lockdown", "Trace all")
	group.cols <- cl.4[as.numeric(as.factor(group.order))]
	# Random selection of 5000 points
	sel <- sample(1:NSim, 10000, replace=FALSE)

	beeswarm(log10(sim_episize+1)[sel] ~ group.order[sel], data = out_sim,
	    pch = 19, cex=0.5,  pwcol = group.cols[sel],
	    xlab = "", ylab = "",
	    labels = group.labels, cex.lab =1.4, cex.axis=1.4, las=1)
	    ##legend("topright", legend = c("", ""), title = "")
	mtext(expression("log"[10]*"(Epidemic size)"), 2, outer = F, font = 1, line = 2.9, cex = 1.5, las = 3)
}
plot_beeswarm ()	


## Plot correlations in relative epidemic size and gradients
Out_metapop$rel.episize_traceI[which(Out_metapop$rel.episize_traceI >5)] <- 1
plot_corr.traceI <- ggplot(data.frame(x=Out_metapop$gradient_traceI, y= Out_metapop$rel.episize_traceI), aes(x=x, y=y)) + 
geom_point(size=1.9, aes(color = Out_metapop$beta)) +
labs(title = "Trace symptomatic only") +
ylab("Relative epidemic size") +
xlab("") + 
scale_color_gradient(low = "blue", high = "red", name="beta") + theme(legend.position = "none")

Out_metapop$rel.episize_lockd[which(Out_metapop$rel.episize_lockd >5)] <- 1
plot_corr.Lockd <- ggplot(data.frame(x=Out_metapop$gradient_lockd, y= Out_metapop$rel.episize_lockd), aes(x=x, y=y)) + 
geom_point(size=1.9, aes(color = Out_metapop$beta)) +
labs(title = "Local lockdown") + xlab("Urban-rural gradient (correlation strength)") + scale_color_gradient(low = "blue", high = "red", name="beta") + theme(legend.position = "none") +
ylab("")


plot_corr.traceEAI <- ggplot(data.frame(x= Out_metapop$gradient_traceEAI, y= Out_metapop$rel.episize_traceEAI), aes(x=x, y=y)) + 
geom_point(size=1.9, aes(color = Out_metapop$beta)) +
labs(title = "Trace all") +
xlab("" ) + scale_color_gradient(low = "blue", high = "red", name="beta") +
ylab("")

plot_Corr.Episize.Gradient <- ggarrange(plot_corr.traceI, plot_corr.Lockd, plot_corr.traceEAI, ncol=3, nrow=1, common.legend = TRUE)
plot_Corr.Episize.Gradient


#####
# Plot maps

dat.grav <- data.frame(east = pop_eastings[sel_pop], north= pop_northings[sel_pop], 
			grav=scale(Mat_pop.grav.index[which(Hypercube$gavD[which(Hypercube$scenario==1)]==1)[1],]))
plot_map.Wales.gravity <- ggplot(data = map_LSOA.Wales) +  geom_sf() + ##theme_classic() +
    xlim(range(pop_eastings[sel_pop])) + ylim(range(pop_northings[sel_pop])) +
    geom_point(data=dat.grav, aes(x = east, y = north, color= grav)) +     
    scale_color_gradient(low="darkblue", high="yellow", name= "Connectivity") +
    xlab("Easting") + ylab("Northting")    
plot_map.Wales.gravity



pop_col <- ifelse(pop_is.sel==1, "#E7298A", "#1B9E77")
plot_map.UK.LSOA <- ggplot(data = map_LSOA.Wales) + geom_sf() + theme_classic() +
     xlim(range(pop_eastings)) + ylim(range(pop_northings)) +
    geom_point(data=LSOA_Pop, aes(x = Eastings, y = Northings), color= pop_col) +
    xlab("Easting") + ylab("Northting")    
plot_map.UK.LSOA

ggarrange(plot_map.UK.LSOA, plot_map.Wales.gravity, ncol=2, nrow=1, common.legend=F)



###############
# Summary for results

# p.trace.EAI at which no scenario has rel. epidemic size > 0.05
max(round(Hypercube$p.trace.EAI [which(Hypercube$scenario==2)][which(metapop_rel.episize_traceEAI>0.05)],2), na.rm=T)

# 
targetval.1 = 0.05
max(round(Hypercube$p.asympt [which(Hypercube$scenario==3)][which(metapop_rel.episize_traceI<targetval.1)],2), na.rm=T)
min(round(Hypercube$p.resist [which(Hypercube$scenario==3)][which(metapop_rel.episize_traceI<targetval.1)],2))
max(round(Hypercube$delta [which(Hypercube$scenario==3)][which(metapop_rel.episize_traceI<targetval.1)],2))



# Regional lockdown to lead relative epidemic sizes < 0.5
targetval = 0.05
min(round(Hypercube$lockdown.treshold [which(Hypercube$scenario==4)][which(metapop_rel.episize_lockd>targetval)],2), na.rm=T)
min(round(Hypercube$lockdown.factor [which(Hypercube$scenario==4)][which(metapop_rel.episize_lockd>targetval)],2), na.rm=T)
max(round(Hypercube$lockdown.time [which(Hypercube$scenario==4)][which(metapop_rel.episize_lockd>targetval)],2), na.rm=T)
max(round(Hypercube$lockdown.dist [which(Hypercube$scenario==4)][which(metapop_rel.episize_lockd>targetval)],2), na.rm=T)

plot(Hypercube$lockdown.treshold [which(Hypercube$scenario==4)], local_rel.episize_lockd, ylim=c(0,1), xlab= "% infectious (I) individuals isolated", ylab="", pch = 19, cex.lab=1.4, col =  clgr.1[round(100*Hypercube$p.asympt[which(Hypercube$scenario==3)])])


`
