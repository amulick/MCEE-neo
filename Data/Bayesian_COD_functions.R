#######################################################################
# Functions script accompanying 'A Bayesian hierarchical model with integrated 
# covariate selection and misclassification matrices to estimate 
# neonatal and child causes of death' 
#
# by Mulick AR, Oza S, Prieto-Merino D, Villavicencio F, Cousens S, Perin J
#
# BAYESIAN MODEL ESTIMATION
# f.e1:  accepts study (death) input with misclassification matrices and covariates
#
# PREDICT COD DISTRIBUTION (OUT-OF-SAMPLE)
# f.par:  gather model coefficients from MCMC array and other data returned by f.e1
# f.pr2:  Predict COD distributions from covariates (supplied in function) and data returned by f.par
# f.pci2: Calculate credible intervals from data returned by f.pr2
#
# 11 December 2020
#######################################################################



##################################################Z
####
####   Model estimation
####  
f.e1 <- function(PARA=1, STUD=my.stud, DEAT=my.deat, VDT=my.vdt, 
                 VXF=my.vxf, VXR=my.vxr, MODL=my.modl, LAMB=my.lamb, 
                 RSDL=my.rsdl, QUAD=my.quad, CHAS=my.chas, ITER=my.iter, 
                 BURN=my.burn, THIN=my.thin, SEED=123, PRIN=my.prin, 
                 SAMC=my.samc, SAVX=my.savx, SADE=my.sade, NAME=NULL) {
  ##
  ##  PARA = 1 to run CHAINS in parallel and 0 to run CHAINS in serial
  ##  ------------------------- #
  ##  STUD = dataset of studies
  ##  DEAT = dataset of reported causes of death by study with missclassification matrix
  ##  VDT  = names of true causes of death (string vector)
  ##  VXF  = names of fixed effects variables (string vector)
  ##  VXR  = name of random effects variable (string vector)
  ##  -------------- #
  ##  MODL = txt file with BUGS model
  ##  LAMB = Lambda value for Poisson
  ##  RSDL = limit for the random effects SD
  ##  QUAD = NULL if no quadratic terms desired (default); !is.null(QUAD) adds quadratic terms
  ##  -------------- #
  ##  CHAS = number of chains
  ##  ITER = number of mcmc SAMPLES
  ##  BURN = number of draws used for burn-in (currently arbitrarily set to 500)
  ##  THIN = thin factor 
  ##  SEED = random seed for jags
  ##  -------------- #
  ##  PRIN = Print stage of computation
  ##  SAMC = Save MCMC samples? TRUE/FALSE
  ##  SAVX = T, Save explanatory variables
  ##  SADE = T, Save matrix of deaths
  ##  NAME = Optional name for the model
  ##
  ptime <- proc.time()
  
  ## Prepare DATA
  if (QUAD) {
    STUD[, paste0(VXF, ".q")] <- STUD[, VXF]^2
    VXF <- c(VXF, paste0(VXF, ".q"))
  }
  # Add numeric variable for Re group to data
  STUD$rG <- as.numeric(as.factor(STUD[,VXR]))
  # drop deaths from studies not in estimation sample
  DEAT <- droplevels(DEAT[DEAT$sid %in% STUD$sid, ])
  # Save DEAT if required for return
  if(SADE) sde <- DEAT else sde <- NULL
  
  ## Prepare JAGS DATA
  # Constants for model
  K  <- length(VXF)+1            # Number of fixed-effects parameter
  S  <- dim(STUD)[1]             # Number of studies 
  C  <- length(VDT)              # Number of true causes of DEAT
  rN <- length(unique(STUD$rG))  # number of RE groups
  # First reported death
  fi <- as.numeric(tapply(1:nrow(DEAT), DEAT$sid, min))
  # Last reported death
  la <- as.numeric(tapply(1:nrow(DEAT), DEAT$sid, max))
  # Number of deaths in the study
  N <- as.numeric(tapply(DEAT$n, DEAT$sid, sum))
  # Create GM matrix
  GM <- as.matrix(DEAT[, c(VDT,"n")])
  ## Scale & centre continuous covariates (including squared terms)
  XM <- cbind(CONS=rep(1,S), as.matrix(STUD[,VXF]))
  XM[, -1] <- apply(XM[, -1], 2, scale)
  # List with JAGS data
  jd.r1 <- list(S = S, C = C, K = K, rN = rN, lambda = LAMB, 
                rsdlim = RSDL, fi = fi, la = la, N = N,
                rG = STUD$rG, XM = XM, GM = GM)
  
  ## Parameters to follow
  jp.r1 <- c("B","re","sd")
  ## Initial values, randomly generated from distributions with large variance
  ## except 'sd' parameters, which have fixed starting values
  ji.r1 <- function() list(B=matrix(c(rep(NA,K),rnorm((K*(C-1)),0,8)),K,C),
                           re=matrix(c(rep(NA,rN),rnorm((rN*(C-1)),0,8)),rN,C),
                           sd=c(NA,rep(0.03,(C-1))))
  
  # Remove all unnecessary objects before calling JAGS
  rm(DEAT, XM, GM, la, fi, N)

  ## SAMPLES
  if(PARA==1){
    ## Set parameters in the Global environment
    assign("my.iter", ITER, envir=.GlobalEnv)
    assign("my.burn", BURN, envir=.GlobalEnv)
    assign("my.thin", THIN, envir=.GlobalEnv)
    assign("my.modl", MODL, envir=.GlobalEnv)
    jm.r1 <- jags.parallel(data=jd.r1, parameters.to.save=jp.r1, 
                           inits=ji.r1, model.file=paste(my.modl), 
                           jags.seed = SEED, n.chains=CHAS, 
                           n.iter=my.iter, n.burnin=my.burn, 
                           n.thin=my.thin, 
                           export_obj_names=c("my.iter","my.burn","my.thin","my.modl")) # Parallel chains
  }
  if(PARA==0){
    set.seed(SEED)
    jm.r1 <- jags(data=jd.r1, parameters.to.save=jp.r1, inits=ji.r1, 
                  model.file=paste(MODL), n.chains=CHAS, n.iter=ITER,
                  n.burnin=BURN, n.thin=THIN) # NON-parallel chains
  }
  
  ###  Prepare Return list  ###
  # Define parameters of the model
  para <- list(name = NAME, time = ptime, PARA = PARA, 
               VDT = VDT, VXF = VXF, VXR = VXR, MODL = MODL, 
               LAMB = LAMB, RSDL = RSDL, QUAD = QUAD)
  # Save means and SD of predictors for standardisation:
  dv <- data.frame(xvar=VXF, mean = colMeans(STUD[, VXF], na.rm = T), 
                   sd = apply(STUD[, VXF], 2, sd, na.rm=T), 
                   stringsAsFactors = F, row.names = NULL)
  # decide what prediction variables to save
  if(SAVX==T) svx <- c("sid",VXR,"rG",VXF) else svx <- c("sid",VXR,"rG")
  # Remove MCMC samples if not wanted
  if(SAMC==F){
    jm.r1$BUGSoutput$sims.array <- NULL
    jm.r1$BUGSoutput$sims.list <- NULL # don't save these
    jm.r1$BUGSoutput$sims.matrix <- NULL # don't save these
  }
  ptime <- proc.time() - ptime
  if(PRIN==1) print(paste("time = ",ptime[3]))
  return(list(param = para, Vars = dv, Studies = STUD[,svx], 
              Deaths = sde, output=jm.r1))
}
###
###  END OF FUNCTION
###
################################################ #



##################################################Z
####
####   Gather model data from f.e1 for prediction
####
f.par <- function(MO, NP=500){
  
  ## MO    JAGS object with posterior MCMC coefficient distribution
  ## NP    Number of coefficient sets from which to estimate credible intervals
  
  ## Cubes of fixed and random effects coefficients
  
  # extract relevant components from the prediction object
  VXF  <- MO$Vars$xvar # array of means of covariates
  VDT <- MO$param$VDT
  SA  <- MO$output$BUGSoutput$sims.array  # array of estimations
  # Recover simulation parameters from BUGS output
  K   <- length(VXF) # number of explanatory variables (excluding intercept)
  C   <- length(VDT) # number of underlying causes of death
  I   <- dim(SA)[1]  # number of iterations
  H   <- dim(SA)[2] # number of chains
  cB  <- which(substr(dimnames(SA)[[3]], 1, 2) == 'B[') # betas
  cR  <- which(substr(dimnames(SA)[[3]], 1, 3) == 're[') # random effects
  R   <- length(cR)/(C-1) # number of random effects
  print(paste("K=",K,", C=",C, ", I=",I, ", H=",H, ", R=",R))
  # Samples to keep from the simulations
  rc <- sample(H, NP, replace=T)  # random selection of chains
  ri <- sample(I, NP, replace=F) # random iterations
  # Prepare arrays
  BM <- array(NA, dim = c(K+1, C, NP))  # array of fixed betas
  dimnames(BM) <- list(c("constant",VXF), VDT, paste(rc,ri,sep="."))
  RM <- array(NA, dim = c(R, C, NP))  # array of random effects
  dimnames(RM) <- list(MO$Studies[match(c(1:R), MO$Studies$rG), 2], VDT, paste(rc,ri,sep="."))
  # populate arrays
  for (i in 1:NP){
    # Matrix of coefficients
    BM[,,i] <- cbind(rep(0,K+1), matrix(SA[ri[i],rc[i],cB], nrow=(K+1) )) #betas
    RM[,,i] <- cbind(rep(0,R),   matrix(SA[ri[i],rc[i],cR], nrow=R)) #random effects
  }
  return(list(Model=MO$param$name, VX=MO$Vars, VDT=VDT, RE=MO$Studies, Summary=MO$output$BUGSoutput$summary, BM=BM, RM=RM))
}
###
###  END OF FUNCTION
###
################################################ #



################################################ #
####
####   Predictions with an array of coefficients from f.par()
####   
f.pr2 <- function(PA, PD, RT, ID="id", PE=""){
  
  ## PA   object produced with function f.par()
  ## PD   data set with covariates to be used in the prediction
  ## RT   list of applicable random terms for each ID in PD
  ## ID   variable in PD that uniquely identify observations
  ## PE   Period label ("early","late","any")
  
  ## This function makes predictions with fixed effects only and then adds 
  ## ONE random effect term selected at random among some candidates in RT in each MCMC iteration.
  
  # Prepare prediction data
  VXF <- PA$VX$xvar  # list of prediction covariates
  S   <- dim(PD)[1] 
  K   <- dim(PA$BM)[1]-1 # number of covariates excluding constant
  C   <- dim(PA$BM)[2]   # number of causes of death
  N   <- dim(PA$BM)[3]   # number of simulations
  if(sum(!(VXF %in% names(PD)))>0){
    stop(paste('Variables', paste(VXF[!(VXF %in% names(PD))], collapse=" , "), ' not in the data'))
  }
  # Prepare raw variables dataset from prediction sample
  for(k in 1:K) PD[,VXF[k]] <- (PD[,VXF[k]] - PA$VX$mean[k])/PA$VX$sd[k]
  DX <- cbind(cons=rep(1,S), as.matrix(PD[,VXF]))
  # Prepare matrix of selected random terms for each ID
  RE <- array(NA, dim = c(S, N))
  for(i in 1:S) {
    if(length(RT[[i]]>0)) {
      RE[i,] <- sample(as.character(RT[[i]]), N, replace=T)
    } else {
      RE[i,]  <- sample(dimnames(PA$RM)[[1]], N, replace=T)
    }
  }
  # Prepare matrix to store predictions
  LF <- array(NA, dim = c(S, C, N))  # matrix of fixed Predictions
  dimnames(LF) <- list(PD[,ID], dimnames(PA$BM)[[2]], dimnames(PA$BM)[[3]])
  LR <- LF  # matrix of random Predictions
  
  # PREDICTION for each iteration
  for (i in 1:N){
    LF[,,i] <- DX %*% PA$BM[,,i] # logodds with fixed effects
    LR[,,i] <- LF[,,i] + PA$RM[RE[,i],,i] # logodds of random effects
  }
  PF <- aperm(apply(LF, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
  PR <- aperm(apply(LR, c(1,3), function(x) exp(x)/sum(exp(x))), perm=c(2,1,3))
  return(list(Model=PA$Model, Period=PE, VX=PA$VX, PF=PF, PR=PR) )
}
###
###  END OF FUNCTION
###
################################################ #



############################################## #
###
###   CALCULATE CIs from PREDICTIONS returned by f.pr2
###
f.pci2 <- function(PM, CI=95){
  
  ## PM   object from function f.pr2()
  ## CI   Credible interval coverage
  
  q <- c((100-CI)/200, (100+CI)/200)
  
  for(c in 1:(dim(PM$PF)[2])){
    xf <- data.frame(t(apply(PM$PF[,c,], 1, function(x) c(mean(x, na.rm=T), sd(x, na.rm=T), quantile(x, q[1], na.rm=T), quantile(x, q[2], na.rm=T)))))
    xf$cause <- dimnames(PM$PF)[[2]][c]
    xf$type <- "fixed"
    xf$iso  <- sapply(rownames(xf), function(x)  strsplit(x, ".",fixed=T)[[1]][1])
    xf$year <- as.numeric(sapply(rownames(xf), function(x)  strsplit(x, ".",fixed=T)[[1]][2]))
    xr <- data.frame(t(apply(PM$PR[,c,], 1, function(x) c(mean(x, na.rm=T), sd(x, na.rm=T), quantile(x, q[1], na.rm=T), quantile(x, q[2], na.rm=T)))))
    xr$cause <- dimnames(PM$PR)[[2]][c]
    xr$type <- "random"
    xr$iso  <- sapply(rownames(xr), function(x)  strsplit(x, ".",fixed=T)[[1]][1])
    xr$year <- as.numeric(sapply(rownames(xr), function(x)  strsplit(x, ".",fixed=T)[[1]][2]))
    if(c==1) PP <- rbind(xf,xr) else PP <- rbind(PP,xf,xr)
  }
  names(PP)[1:4] <- c("me","sd","lo","hi")
  PP$model <- PM$Model
  PP$period <- PM$Period
  return(PP[,c(7:10,6,5,1:4)])
}
###
###  END OF FUNCTION
###
################################################ #



