######################################
######################################
####    R script to accompany:
####    'A Bayesian hierarchical model with integrated covariate selection 
####    and misclassification matrices to estimate neonatal causes of death'
####    Mulick et al. (2021)
####    7 January 2021
######################################
######################################

# This code reproduces the analysis from the final model in the published paper.
# All the data are available in the 'Data' subfolder and can be read with this R code.

rm(list=ls())     # Clear workspace



#######   PACKAGES  ##########
#install.packages(c("R2jags", "doParallel", "plyr", "tidyverse", "data.table"))
require(R2jags)
require(doParallel)
require(plyr)
require(tidyverse)
require(data.table)



#######   LOAD DATA AND FUNCTIONS  ##########
my.data  <- "Bayesian_COD_data.RData"
load(file=paste0("Data/", my.data))
my.func  <- "Bayesian_COD_functions.R"
source(paste0("Data/", my.func))



#######   DEFINE VALUES FOR PARAMETERS IN FUNCTIONS   ##########
my.modl  <- "Data/Bayesian_COD_model.txt"    # Text file with Bayesian model
my.lamb  <- 30     # lambda parameter for LASSO
my.rsdl  <- 0.21   # random effects standard deviation upper limit
my.vxr   <- vxr    # random effects variable name
my.vxf   <- vxf    # fixed effects variable names
my.vdt   <- vdt    # true causes of death names
my.quad  <- FALSE  # quadratic terms required?
my.iter  <- 10000  # number of mcmc SAMPLES
my.burn  <- 2000   # number of draws used for burn-in
my.chas  <- 4      # number of chains
my.thin  <- 1      # thin factor 
my.prin  <- 1      # don't print results while calculating?
my.samc  <- T      # Save JAGS object?
my.sade  <- F      # Save matrix of deaths?
my.savx  <- F      # Save matrix of explanatory variables?



#######   RUN MODEL   ##########
nam <- paste("va",my.lamb,my.rsdl,sep="_")
assign(nam,f.e1(STUD=studies, DEAT=deaths, LAMB=my.lamb, RSDL=my.rsdl, NAME=nam))



#######   MAKE OUT OF SAMPLE PREDICTIONS FROM MODEL   ##########
mcva <- pva <- list()

# matrix of model coefficients
mcva <- f.par(get(nam))

# define relevant random effect terms for each country
RVA=sapply(data.predict$isocode, function(x) c(filter(studies, isocode==x, natrep==T)$studyid))

  # do India manually - one nationally representative study (studyid==1), use for all states
  RVA[nchar(names(RVA))==4] <- as.integer(1)

# make predictions (overall neonatal period: per.early==per.late==0)
pva <- f.pr2(mcva, PD=mutate(data.predict, per.early=0, per.late=0, premvslbw=1, id=paste(isocode,year,sep=".")), PE="overall", RT=RVA)

# Do India manually: collapse state predictions to entire country
India <- pva
India$PF <- India$PF[which(nchar(rownames(India$PF)) == 9), , ] 
India$PR <- India$PR[which(nchar(rownames(India$PR)) == 9), , ] 

  # Years of interest
  years <- as.numeric(unique(substr(rownames(India$PF), 6, 9)))
  
  # Create new objects for predicted distributions
  PFNew <- array(NA, dim = c(length(years), dim(India$PF)[2:3]),
                 dimnames = list(paste('IND', years, sep = '.'),
                                 vdt, NULL))
  PRNew <- array(NA, dim = c(length(years), dim(India$PR)[2:3]),
                 dimnames = list(paste('IND', years, sep = '.'),
                                 vdt, NULL))
  
  # Loop through MCMC samples and create national COD distribution
  for (j in 1:dim(India$PF)[3]) {
    
    for (k in 1:2) {
      
      # Fixed effects
      if (k == 1) x <- as.data.frame(India$PF[, , j])
      
      # Random effects 
      if (k == 2) x <- as.data.frame(India$PR[, , j])
      
      # Recover envelopes
      x$isocode <- substr(rownames(x), 1, 4)
      x$year <- substr(rownames(x), 6, 9)
      x <- merge(x, data.predict[, c('isocode', 'year', 'envelope')],
                 by = c('isocode', 'year'), all.x = T)
      
      # Collapse Indian states
      x[, paste(vdt)] <- x[, paste(vdt)] * x$envelope
      x <- x[, !names(x) %in% 'isocode']    
      x <- aggregate(x[, !names(x) %in% 'year'], by = list(x$year), FUN = sum, na.rm = T)
      
      # Re-calculate fractions
      x[, paste(vdt)] <- x[, paste(vdt)] / rowSums(x[, paste(vdt)], na.rm = T)
      x <- x[order(x$Group.1), ]
      
      # Save national fractions in new object
      if (k == 1) PFNew[, , j] <- as.matrix(x[, paste(vdt)])
      if (k == 2) PRNew[, , j] <- as.matrix(x[, paste(vdt)])
      
    }
    
  }
  
  # Update Indian object with National data
  India$PF <- PFNew
  India$PR <- PRNew
  rm(PFNew, PRNew)
  
  

#######   CALCULATE CI FOR PREDICTIONS   ##########

# get 95% CIs
  #India
  pci <- list()
  pci[[1]] <- f.pci2(India)
  #rest
  pci[[2]] <- f.pci2(pva)

pci <- do.call(rbind,pci)
pci <- pci[nchar(pci$iso)==3, ] #exclude Indian states

