# Literature:

  * Tanaka, E. Simple outlier detection for a multi‐environmental field trial. Biometrics. 2020; 1– 9. https://doi.org/10.1111/biom.13216

# R code:

args <- commandArgs(trailingOnly=TRUE)

istart <- as.integer(args[1])
iend <- as.integer(args[2])

library(asreml)
library(myf)
library(doParallel)
library(foreach)

no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)

load("1-caige2016.Rdata")

nsim <- 1000
true.final1 <- c(m1$gammas, coef(m1)$fixed[, 1])
allest.final1 <- matrix(NA, nrow=nsim, ncol=length(true.final1))
colnames(allest.final1) <- names(true.final1)
true.final0 <- c(m0$gammas, coef(m0)$fixed[, 1])
allest.final0 <- matrix(NA, nrow=nsim, ncol=length(true.final0))
colnames(allest.final0) <- names(true.final0)
ff <- coef(m1)$fixed[,1]

result.array <- array(NA, dim=c(22, nsim, nobs), dimnames=list("Statistics"=c("pi0", "pti0", "qti0", "qi0", "ti0", "ri0", "CDi0", "LRi0",
                                                                              "pi1", "qi1", "pti1", "qti1", "ti1", "ri1", "CDi1", "LRi1",
                                                                              "MSi0", "WTi0", "WTi1", "LRti0", "MSi1", "LRti1"),
                                                               "Simulation"=1:nsim, "Observation"=1:nobs))
# contamination to be introduced
nout.mat <- matrix(c(5, 5, 5, 3, 3, 3, 1, 1, 1), ncol=3, nrow=3)
nout.all <- apply(nout.mat, 2, sum)
nout <- c(5, 3, 1)

obsOut1 <- matrix(NA, nrow=nsim, ncol=nout.all[1])
obsOut2 <- matrix(NA, nrow=nsim, ncol=nout.all[2])
obsOut3 <- matrix(NA, nrow=nsim, ncol=nout.all[3])

errorOut1 <- matrix(NA, nrow=nsim, ncol=nout.all[1])
errorOut2 <- matrix(NA, nrow=nsim, ncol=nout.all[2])
errorOut3 <- matrix(NA, nrow=nsim, ncol=nout.all[3])

true.sim <- list()
est.sim <- list()

update2.asreml <- function(object, maxiter=500, text=".") {
  iter <- 1
  while(!object$conv) {
    object <- update(object, trace=F, stepsize=0.0001)
    iter <- iter + 1
    if(iter > maxiter) {
      line <- paste0("CAIGE\tLow\t", aseed, "\t", isim, "\t", text)
      write(line, file=paste0("logCAIGElow", istart, "_", iend, ".txt"), append=TRUE)
      break
    }
  }  
  object
}

# # simulation --------------------------------------------------------------
for(isim in istart:iend) {
  set.seed(20180119 + isim)
  # choose 3 sites to contaminate 
  site.cont <- sample(enams, ncol(nout.mat))
  
  ########## data simulation
  simdat.df <- dat
  simdat.df$Yield <- NA
  simdat.df$Yield <- ff["(Intercept)"] + ff[paste0("LocSite_", as.character(simdat.df$LocSite))] 
  # first generate the genetic effects
  ugxe <- GEmatchol %*% rnorm(m * t)
  nn <- rownames(ugxe); ugxe <- setNames(as.vector(ugxe), nn)
  
  true.sim[[isim]] <- list()
  true.sim[[isim]][["ugxe"]] <- ugxe
  est.sim[[isim]] <- list()
  
  for(aexpt in enams) {
    atrial.df <- droplevels(subset(dat, LocSite==aexpt))
    nrow <- nlevels(atrial.df$Row)
    ncol <- nlevels(atrial.df$Column)
    nrep <- nlevels(atrial.df$Block)
    urep.mod <- rnorm(nrep, 0, sqrt(est.final[aexpt, 'block']))
    ucol.mod <- rnorm(ncol, 0, sqrt(est.final[aexpt, 'ranc']))
    urow.mod <- rnorm(nrow, 0, sqrt(est.final[aexpt, 'ranr']))
    gxe.name <- paste0(aexpt, ":", as.character(atrial.df$Geno))
    error.mod <- rnorm(nrow*ncol, 0, sd=sqrt(est.final[aexpt, 'sigm']))
    simdat.df[simdat.df$LocSite==aexpt, "Yield"] <-  simdat.df[simdat.df$LocSite==aexpt, "Yield"] + 
      ugxe[gxe.name] +
      urep.mod[as.numeric(atrial.df$Block)] +
      ucol.mod[as.numeric(atrial.df$Column)] + 
      urow.mod[as.numeric(atrial.df$Row)] + error.mod
    
    true.sim[[isim]][[aexpt]][["urep"]] <- urep.mod
    true.sim[[isim]][[aexpt]][["ucol"]] <- ucol.mod
    true.sim[[isim]][[aexpt]][["urow"]] <- urow.mod
    true.sim[[isim]][[aexpt]][["error"]] <- error.mod
    
    # add contamination 
    if(aexpt==site.cont[1]) {
      obsOut1[isim, ] <- sample(atrial.df$Obs, nout.all[1])
      errorOut1[isim, ] <- c(rnorm(nout[1], 4 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])),
                             rnorm(nout[1], 7 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])),
                             rnorm(nout[1], 10 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])))
      simdat.df[obsOut1[isim, ], "Yield"] <- simdat.df[obsOut1[isim, ], "Yield"] + errorOut1[isim, ]
    } else {
      if(aexpt==site.cont[2]) {
        obsOut2[isim, ] <- sample(atrial.df$Obs, nout.all[2])
        errorOut2[isim, ] <- c(rnorm(nout[2], 4 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])),
                               rnorm(nout[2], 7 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])),
                               rnorm(nout[2], 10 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])))
        simdat.df[obsOut2[isim, ], "Yield"] <- simdat.df[obsOut2[isim, ], "Yield"] + errorOut2[isim, ]
      } else {
        if(aexpt==site.cont[3]) {
          obsOut3[isim, ] <- sample(atrial.df$Obs, nout.all[3])
          errorOut3[isim, ] <- c(rnorm(nout[3], 4 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])),
                                 rnorm(nout[3], 7 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])),
                                 rnorm(nout[3], 10 * sqrt(est.final[aexpt,'sigm']), sqrt(est.final[aexpt,'sigm'])))
          simdat.df[obsOut3[isim, ], "Yield"] <- simdat.df[obsOut3[isim, ], "Yield"] + errorOut3[isim, ]
        }
      }
    }
  }
  
  ######## DIAG model 
  m0sim <- m0
  m0sim <- update(m0sim, data=simdat.df, trace=F)
  m0sim <- update2.asreml(m0sim, text="m0sim")
  m0sim <- update(m0sim, aom=T, trace=F)
  allest.final0[isim, ] <- c(m0sim$gammas, coef(m0sim)$fixed[, 1])
  
  # studentised conditional residual
  result.array["ti0", isim, ] <- ti0 <- m0sim$aom$R[, "stdCondRes"]
  result.array["ri0", isim, ] <- m0sim$aom$R[, "rinvE"]
  
  # p-value 
  result.array["pi0", isim, ] <- p0 <- 2 * (1 - pnorm(abs(ti0)))
  
  for(aexpt in enams) {
    ii <- simdat.df$Obs[simdat.df$LocSite==aexpt]
    ntmp <- length(ii)
    simdat2.df <- droplevels(simdat.df[ii,])
    # q-value
    result.array["qi0", isim, ii] <- p.adjust(p0[ii], method="holm")
    
    # likelihood ratio analytic from fixing gamma
    sub_m0sim <- asreml(Yield ~ 1, random=~Geno + Block + Row + Column, data=droplevels(simdat.df[ii,]), trace=F)
    pp <- 1
    result.array["LRi0", isim, ii] <- ((ntmp - pp - 1) * log((ntmp - pp - 1) / (ntmp - pp -ti0[ii]^2)) - log(ti0[ii]^2)) 
    
    # GCD
    isFixed <- asreml:::asreml.guzpfx(sub_m0sim$gammas.con) %in% c("Fixed", "Boundary")
    H <- solve(asreml:::asreml.ltri2mat(sub_m0sim$ai)[!isFixed, !isFixed])
    theta <- sub_m0sim$gammas[!isFixed] * sub_m0sim$sigma2
    result.array["CDi0", isim, ii] <- foreach(iobs=1:ntmp, .combine=c) %dopar% {
      Mtmp2 <- sub_m0sim
      simdat2.df <- droplevels(simdat.df[ii,])
      simdat2.df$Yield[iobs] <- NA
      Mtmp2 <- update(Mtmp2, data=simdat2.df, trace=F, maxiter=1)
      theta_i <- Mtmp2$gammas[!isFixed] * Mtmp2$sigma2
      t(theta - theta_i) %*% H %*% (theta - theta_i)
    }
    
    # MSOM & VSOM & LR -- Exact
    result.array[c("MSi0", "WTi0", "pti0", "LRti0"), isim, ii] <- foreach(iobs=1:ntmp, .combine=cbind) %dopar% {
      simdat2.df[["out1"]] <- 0
      simdat2.df[["out1"]][iobs] <- 1
      m0simMSOM <- update(sub_m0sim, data=simdat2.df, fixed.=~ . + out1, trace=F)
      m0simMSOM <- update2.asreml(m0simMSOM, text=paste0("m0simMSOM_iobs", iobs))
      WTi0 <- m0simMSOM$yssqu[["out1"]] / m0simMSOM$sigma2
      m0simVSOM <- update(sub_m0sim, data=simdat2.df, random.=~ . + out1, trace=F)
      m0simVSOM <- update2.asreml(m0simVSOM, text=paste0("m0simVSOM_iobs", iobs))
      c(coef(m0simMSOM, list=T)$out1[,1], WTi0, 1 - pchisq(WTi0, 1), abs(2 * (m0simVSOM$loglik - sub_m0sim$loglik)))
    }
    result.array["qti0", isim, ii] <- p.adjust(result.array["pti0", isim, ii], method="holm")
  }
  
  ######## US model
  m1sim <- m1
  m1sim <- update(m1sim, data=simdat.df, trace=F)
  m1sim <- update2.asreml(m1sim, text="m1sim")
  m1sim <- update(m1sim, aom=T, trace=F)
  allest.final1[isim, ] <- c(m1sim$gammas, coef(m1sim)$fixed[, 1])
  
  est.sim[[isim]][["ugxe"]] <- coef(m1sim, pattern="LocSite:Geno")[, 1]
  
  # studentised conditional residual
  result.array["ti1", isim, ] <- ti <- m1sim$aom$R[, "stdCondRes"]
  result.array["ri1", isim, ] <- m1sim$aom$R[, "rinvE"]
  
  # p-value & q-value
  result.array["pi1", isim, ] <- p1 <- 2 * (1 - pnorm(abs(ti)))
  result.array["qi1", isim, ] <- p.adjust(p1, method="holm")
  
  # likelihood ratio analytic form fixing gamma
  pp <- 7
  result.array["LRi1", isim, ] <- ((nobs - pp - 1) * log((nobs - pp - 1) / (nobs - pp -ti^2)) - log(ti^2)) 
  
  # GCD
  isFixed <- asreml:::asreml.guzpfx(m1sim$gammas.con) %in% c("Fixed", "Boundary")
  H <- solve(asreml:::asreml.ltri2mat(m1sim$ai)[!isFixed, !isFixed])
  theta <- m1sim$gammas[!isFixed]
  result.array["CDi1", isim, ] <- foreach(iobs=1:nobs, .combine=c) %dopar% {
    Mtmp2 <- m1sim
    simdat2.df <- simdat.df
    simdat2.df$Yield[iobs] <- NA
    Mtmp2 <- update(Mtmp2, data=simdat2.df, maxiter=1, trace=F)
    theta_i <- Mtmp2$gammas[!isFixed]
    t(theta - theta_i) %*% H %*% (theta - theta_i)
  }
  
  # MSOM & LR -- Exact
  result.array[c("MSi1", "WTi1", "pti1", "LRti1"), isim, ] <- foreach(iobs=1:nobs, .combine=cbind) %dopar% {
    simdat2.df <- simdat.df
    simdat2.df[["out1"]] <- 0
    simdat2.df[["out1"]][iobs] <- 1
    m1simMSOM <- update(m1sim, data=simdat2.df, fixed.=~ . + out1, trace=F)
    m1simMSOM <- update2.asreml(m1simMSOM, text=paste0("m1simMSOM_iobs", iobs))
    WTi1 <- m1simMSOM$yssqu[["out1"]] 
    m1simVSOM <- update(m1sim, data=simdat2.df, random.=~ . + out1, trace=F)
    m1simVSOM <- update2.asreml(m1simVSOM, text=paste0("m1simVSOM_iobs", iobs))
    c(coef(m1simMSOM, list=T)$out1[,1], WTi1, 1 - pchisq(WTi1, 1), abs(2 * (m1simVSOM$loglik - m1sim$loglik)))
  }
  result.array["qti1", isim, ] <- p.adjust(result.array["pti1", isim, ], method="holm")
  
  ############ Robust regression
  ## DIAG model
  obsOutlier0 <- which(result.array["qi0", isim, ] < 0.05)
  if(length(obsOutlier0) > 0) {
    outnames0 <- c()
    for(aout in obsOutlier0) {
      simdat.df[[paste0("out", aout)]] <- 0
      simdat.df[[paste0("out", aout)]][aout] <- 1
      outnames0 <- c(outnames0, paste0("out", aout))
    }
    
    # VSOM 
    m1simVSOM0 <- update(m1sim, data=simdat.df, random.=as.formula(paste("~ . + ", paste(outnames0, collapse=" + "))))
    m1simVSOM0 <- update2.asreml(m1simVSOM0, text="m1simVSOM0_robust_out0")
    est.sim[[isim]][["ugxe_VSOM0"]] <- coef(m1simVSOM0, pattern="LocSite:Geno")[, 1]
    
    # MSOM
    m1simMSOM0 <- update(m1sim, data=simdat.df, fixed.=as.formula(paste("~ . + ", paste(outnames0, collapse=" + "))))
    m1simMSOM0 <- update.asreml(m1simMSOM0, text="m1simMSOM0_robust_out0")
    est.sim[[isim]][["ugxe_MSOM0"]] <- coef(m1simMSOM0, pattern="LocSite:Geno")[, 1]
    
    # Case Deletion
    simdat2.df <- simdat.df
    simdat2.df$Yield[obsOutlier0] <- NA
    m1simCD0 <- update(m1sim, data=simdat2.df)
    m1simCD0 <- update2.asreml(m1simCD0, text="m1simCD0_robust_out0")
    est.sim[[isim]][["ugxe_CD0"]] <- coef(m1simCD0, pattern="LocSite:Geno")[, 1]
    
  } else {
    est.sim[[isim]][["ugxe_VSOM0"]] <- est.sim[[isim]][["ugxe_MSOM0"]] <- est.sim[[isim]][["ugxe_CD0"]] <- rep(NA, m * t)
  }
  
  ## US model
  obsOutlier1 <- which(result.array["qi1", isim,] < 0.05)
  if(length(obsOutlier1) > 0) {
    outnames1 <- c()
    for(aout in obsOutlier1) {
      simdat.df[[paste0("out", aout)]] <- 0
      simdat.df[[paste0("out", aout)]][aout] <- 1
      outnames1 <- c(outnames1, paste0("out", aout))
    }
    
    # VSOM 
    m1simVSOM <- update(m1sim, data=simdat.df, random.=as.formula(paste("~ . + ", paste(outnames1, collapse=" + "))))
    m1simVSOM <- update2.asreml(m1simVSOM, text="m1simVSOM_robust_out1")
    est.sim[[isim]][["ugxe_VSOM1"]] <- coef(m1simVSOM, pattern="LocSite:Geno")[, 1]
    
    # MSOM
    m1simMSOM <- update(m1sim, data=simdat.df, fixed.=as.formula(paste("~ . + ", paste(outnames1, collapse=" + "))))
    m1simMSOM <- update2.asreml(m1simMSOM, text="m1simMSOM_robust_out1")
    est.sim[[isim]][["ugxe_MSOM1"]] <- coef(m1simMSOM, pattern="LocSite:Geno")[, 1]
    
    # Case Deletion
    simdat2.df <- simdat.df
    simdat2.df$Yield[obsOutlier1] <- NA
    m1simCD <- update(m1sim, data=simdat2.df)
    m1simCD <- update2.asreml(m1simCD, text="m1simCD_robust_out1")
    est.sim[[isim]][["ugxe_CD1"]] <- coef(m1simCD, pattern="LocSite:Geno")[, 1]
    
  } else {
    est.sim[[isim]][["ugxe_VSOM1"]] <- est.sim[[isim]][["ugxe_MSOM1"]] <- est.sim[[isim]][["ugxe_CD1"]] <- rep(NA, m * t)
  }
  
  ## US model - MS
  obsOutlier1t <- which(result.array["qti1", isim,] < 0.05)
  if(length(obsOutlier1t) > 0) {
    outnames1t <- c()
    for(aout in obsOutlier1t) {
      simdat.df[[paste0("out", aout)]] <- 0
      simdat.df[[paste0("out", aout)]][aout] <- 1
      outnames1t <- c(outnames1t, paste0("out", aout))
    }
    
    # VSOM 
    m1simVSOMt <- update(m1sim, data=simdat.df, random.=as.formula(paste("~ . + ", paste(outnames1t, collapse=" + "))))
    m1simVSOMt <- update2.asreml(m1simVSOMt, text="m1simVSOMt_robust")
    est.sim[[isim]][["ugxe_VSOM1t"]] <- coef(m1simVSOMt, pattern="LocSite:Geno")[, 1]
    
    # MSOM
    m1simMSOMt <- update(m1sim, data=simdat.df, fixed.=as.formula(paste("~ . + ", paste(outnames1t, collapse=" + "))))
    m1simMSOMt <- update2.asreml(m1simMSOMt, text="m1simMSOMt_robust")
    est.sim[[isim]][["ugxe_MSOM1t"]] <- coef(m1simMSOMt, pattern="LocSite:Geno")[, 1]
    
    # Case Deletion
    simdat2.df <- simdat.df
    simdat2.df$Yield[obsOutlier1t] <- NA
    m1simCDt <- update(m1sim, data=simdat2.df)
    m1simCDt <- update2.asreml(m1simCDt, text="m1simCDt_robust")
    est.sim[[isim]][["ugxe_CD1t"]] <- coef(m1simCDt, pattern="LocSite:Geno")[, 1]
    
  } else {
    est.sim[[isim]][["ugxe_VSOM1t"]] <- est.sim[[isim]][["ugxe_MSOM1t"]] <- est.sim[[isim]][["ugxe_CD1t"]] <- rep(NA, m * t)
  }
  
  ## US model - MS
  obsOutlier1true <- c(obsOut1[isim, ], obsOut2[isim, ], obsOut3[isim, ])
  outnames1true <- c()
  for(aout in obsOutlier1true) {
    simdat.df[[paste0("out", aout)]] <- 0
    simdat.df[[paste0("out", aout)]][aout] <- 1
    outnames1true <- c(outnames1true, paste0("out", aout))
  }
  
  # VSOM 
  m1simVSOMtrue <- update(m1sim, data=simdat.df, random.=as.formula(paste("~ . + ", paste(outnames1true, collapse=" + "))))
  m1simVSOMtrue <- update2.asreml(m1simVSOMtrue, text="m1simVSOMtrue_robust")
  est.sim[[isim]][["ugxe_VSOM1true"]] <- coef(m1simVSOMtrue, pattern="LocSite:Geno")[, 1]
  
  # MSOM
  m1simMSOMtrue <- update(m1sim, data=simdat.df, fixed.=as.formula(paste("~ . + ", paste(outnames1true, collapse=" + "))))
  m1simMSOMtrue <- update2.asreml(m1simMSOMtrue, text="m1simMSOMtrue_robust")
  est.sim[[isim]][["ugxe_MSOM1true"]] <- coef(m1simMSOMtrue, pattern="LocSite:Geno")[, 1]
  
  # Case Deletion
  simdat2.df <- simdat.df
  simdat2.df$Yield[obsOutlier1true] <- NA
  m1simCDtrue <- update(m1sim, data=simdat2.df)
  m1simCDtrue <- update2.asreml(m1simCDtrue, text="m1simCDtrue_robust")
  est.sim[[isim]][["ugxe_CD1true"]] <- coef(m1simCDtrue, pattern="LocSite:Geno")[, 1]
  
  }
