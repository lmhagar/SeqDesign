## this code file contains the functions used to implement Algorithm 1 for
## the PLATINUM-CAN example, as well as the functions to run the confirmation
## simulations and construct bootstrap confidence intervals

## load the necessary packages
require(rjags)
require(coda)
require(ggplot2)
require(foreach)
require(doParallel)
require(doSNOW)

## define expit and logit functions
expit <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(x){
  return(log(x) - log(1-x))
}

## this function is used to simulate data, implement censoring, approximate
## posteriors, implement the marginalization process, and get complementary
## posterior probabilities for the sequential design involving survival data
getCompProbs <- function(n, c_vec, cutpoints, hazards, beta, probs, censors,
                         deltaL, R = 10000, jags_model, prec_beta, prec_logh, 
                         burnin = 1000, nchains = 1, nthin = 2, ndraws = 5000, 
                         seed = 1){
  
  ## the inputs are described as follows:
  ## n: the sample size for the first interim analysis
  ## c_vec: the vector of constants c to reflect how much more data are 
  ##        collected in the subsequent stages
  ## cutpoints: the vector of interval boundaries for the piecewise
  ##            exponential function
  ## beta: a 3 x 1 vector of coefficients for the hazard function
  ## probs: a 2 x 1 vector for the probability of each binary covariate being 1
  ## censors: a vector of probabilities indicating what proportion of patients
  ##          are censored at the beginning of each interval
  ## deltaL: the lower interval endpoint for the hazard ratio hypothesis
  ## R: the number of simulation repetitions used to estimate each sampling distribution
  ## jags_model: the name of the JAGS file (must be in the working directory)
  ## prec_beta: the prior precisions for the regression coefficients
  ## prec_logh: the prior precisions for the log hazards
  ## burnin: the number of burnin iterations for MCMC
  ## nchains: the number of MCMC chains
  ## nthin: the thinning parameter for MCMC
  ## ndraws: the number of MCMC draws that should be retained
  ## seed: a numeric seed used for random number generation
  
  ## this subfunction generates data before accounting for dropout
  dataGen <- function(n, cutpoints, hazards, beta, probs) {
    
    ## Generate covariates x1 and x2
    x <- cbind(rbinom(n, 1, probs[1]), rbinom(n, 1, probs[2]))
    
    ## Compute interval widths
    intervals <- c(0, cutpoints)
    widths <- diff(intervals)
    
    ## Identify unique covariate combinations
    uc <- unique(x)
    
    ## Precompute adjusted hazards for each unique combination
    adjusted_hazards_list <- lapply(1:nrow(uc), function(i) {
      linpred <- sum(beta * c(uc[i, ], uc[i,1]*uc[i,2]))
      hazards * exp(linpred)
    })
    
    ## Initialize results and generate pseudorandom sequence for
    ## CDF inversion
    y <- rep(0, n)
    event <- rep(0, n)
    u <- runif(n)
    
    for (i in seq_along(adjusted_hazards_list)) {
      ## Select patients belonging to subtreatment group
      group_idx <- which(x[,1] == uc[i,1] & x[,2] == uc[i,2])
      
      ## Retrieve adjusted hazards for this group
      adjusted_hazards <- adjusted_hazards_list[[i]]
      
      ## Simulate survival times for this group
      for (j in group_idx) {
        cumulative_hazard <- 0
        for (k in seq_along(widths)) {
          ## Increment cumulative hazard
          incremental_hazard <- adjusted_hazards[k] * widths[k]
          cumulative_hazard <- cumulative_hazard + incremental_hazard
          
          ## Check if event occurs in this interval based on the CDF and u realization
          if (u[j] > exp(-cumulative_hazard)) {
            y[j] <- intervals[k] - (log(u[j]) + cumulative_hazard - incremental_hazard)/adjusted_hazards[k]
            event[j] <- 1
            break
          }
        }
        
        ## If no event experienced, survival time exceeds last interval
        if (event[j] == 0) {
          y[j] <- max(intervals)
        }
      }
    }
    
    return(data.frame(y = y, event = event, x1 = x[, 1], x2 = x[, 2]))
  }
  
  ## this subfunction censors data to account for dropout
  dataDrop <- function(dat, censors, n, cutpoints){
    ## determine time of dropout (if applicable)
    ## dropout indicators 
    drop.mat <- cbind(matrix(rbinom(n*length(cutpoints), 1, rep(censors, each = n)), 
                             byrow = FALSE, nrow = n), rep(1, n))
    
    drop.vec <- apply(drop.mat, 1, function(x){which(x == 1)[1]})
    
    intervals <- c(0, cutpoints)
    for (k in 1:length(censors)){
      dat$event <- ifelse(drop.vec != k, dat$event, 
                          ifelse(dat$y < intervals[k], dat$event, 0))
      dat$y <- ifelse(drop.vec != k, dat$y, 
                    ifelse(dat$y < intervals[k], dat$y, intervals[k]))
    }
    
    ## round survival times up to the nearest day
    dat$y <- ceiling(dat$y)
    
    return(dat)
  }
  
  ## this subfunction processes the data for the JAGS model
  processData <- function(dat, cutpoints){
    
    ## obtain summary statistics for the data
    ## Compute interval widths
    intervals <- c(0, cutpoints)
    widths <- diff(intervals)
    
    ## get the total number of events by interval segment
    ev <- dat$event
    tm <- dat$y
    e_seg <- sum((tm<=cutpoints[1])*ev)
    for (k in 2:length(cutpoints)){
      e_seg <- c(e_seg, sum((tm > cutpoints[k-1] & tm<=cutpoints[k])*ev))
    }
    
    ## get the total events by covariate combination
    ev_x1 <- sum((dat$x1)*ev)
    ev_x2 <- sum((dat$x2)*ev)
    ev_x1x2 <- sum((dat$x1)*ev*dat$x2)
    e_sub <- c(ev_x1, ev_x2, ev_x1x2)
    
    ## get the total time at risk for x1 = 0 and x2 = 0
    tm00 <- subset(dat, dat$x1 == 0 & dat$x2 == 0)$y
    e_00 <- sum(pmin(tm00, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_00 <- c(e_00, sum(pmin(pmax(tm00 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## get the total time at risk for x1 = 1 and x2 = 0
    tm10 <- subset(dat, dat$x1 == 1 & dat$x2 == 0)$y
    e_10 <- sum(pmin(tm10, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_10 <- c(e_10, sum(pmin(pmax(tm10 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## get the total time at risk for x1 = 0 and x2 = 1
    tm01 <- subset(dat, dat$x1 == 0 & dat$x2 == 1)$y
    e_01 <- sum(pmin(tm01, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_01 <- c(e_01, sum(pmin(pmax(tm01 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## get the total time at risk for x1 = 1 and x2 = 1
    tm11 <- subset(dat, dat$x1 == 1 & dat$x2 == 1)$y
    e_11 <- sum(pmin(tm11, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_11 <- c(e_11, sum(pmin(pmax(tm11 - cutpoints[k-1], 0), widths[k])))
    } 
    
    return(list(e_seg, e_sub, e_00, e_10, e_01, e_11))
  }
  
  ## this subfunction returns conditional posterior draws to be input into the
  ## marginalization process
  
  getConditional <- function(jags_model, dat_list, prec_beta, prec_logh, burnin,
                             nchains, nthin, ndraws){
    
    ## initialize model
    model.fit <- jags.model(file=jags_model,
                            data=list(e_seg = dat_list[[1]], e_sub = dat_list[[2]], 
                                      e_00 = dat_list[[3]], e_10 = dat_list[[4]], 
                                      e_01 = dat_list[[5]], e_11 = dat_list[[6]], K = length(dat_list[[1]]), 
                                      prec_beta = prec_beta, prec_logh = prec_logh, zero = 0), 
                            n.chains = nchains, quiet = TRUE)
    
    ## update models and extract posterior draws
    update(model.fit, burnin, progress.bar = "none")
    model.samples <- coda.samples(model.fit, c("beta[1]", "beta[2]", "beta[3]"), 
                                  n.iter=ndraws, thin=nthin, progress.bar = "none")
    
    as.matrix(model.samples)
  }
  
  ## this subfunction implements the marginalization process for the hazard ratio
  ## and gets the complementary posterior probabilities for the hypothesis of interest
  getCompProbs <- function(mod.samples, dat_vec, deltaL){
    mm <- nrow(mod.samples)
    nn <- length(dat_vec)
    
    ## implement marginalization process
    haz.marg <- NULL
    for (k in 1:mm){
      mu0 <- exp(mod.samples[k, 2]*dat_vec)
      w0 <- rexp(nn, 1)
      w0 <- w0/sum(w0)
      
      mu1 <- exp(mod.samples[k, 1] + mod.samples[k, 3]*dat_vec)*mu0
      w1 <- rexp(nn, 1)
      w1 <- w1/sum(w1)
      
      ## get sample from log-HR
      haz.marg[k] <- log(sum(w1*mu1)) - log(sum(w0*mu0)) 
    }
    
    ## get complementary posterior probability
    kd.haz <- density(na.omit(haz.marg))
    np.prob <- mean(pnorm(log(deltaL), ifelse(is.finite(haz.marg), haz.marg, 
                                          max(na.omit(haz.marg)) + 1), kd.haz$bw, lower.tail = TRUE))
    
    ## if complementary probability is very small, use a single normal approximation to the 
    ## posterior for numerical stability
    if (np.prob < 0.00004){
      pnorm(log(deltaL), mean(na.omit(haz.marg)), sd(na.omit(haz.marg)), lower.tail = TRUE)
    } else {
      np.prob
    }
  }
  
  ## we now implement this in parallel across all analyses
  comp.probs <- foreach(i=1:R, .packages=c('rjags', 'coda'), .combine=rbind,
                     .options.snow=opts) %dopar% {
                       
                       ## use previously defined functions to get complementary probabilities
                       ## for the first analysis
                       set.seed(seed + i)
                       dat.acc <- dataGen(n = n, cutpoints = cutpoints, 
                                          hazards = hazards, beta = beta, probs = probs)
                       dat.acc <- dataDrop(dat.acc, censors = censors, n = n, cutpoints = cutpoints)
                       dat.proc <- processData(dat = dat.acc, cutpoints = cutpoints)
                       jags.path <- paste(getwd(), '/', jags_model, sep='')
                       post.cond <- getConditional(jags_model = jags.path, 
                                                   dat_list = dat.proc, prec_beta = prec_beta, 
                                                   prec_logh = prec_logh, burnin = burnin, nchains = nchains,
                                                   nthin = nthin, ndraws = ndraws*nthin)
                       comps <- getCompProbs(mod.samples = post.cond, dat_vec = dat.acc[, 4], deltaL = deltaL)
                       
                       ## augment the data set by repeating this process to add data from later stages
                       for (j in 2:length(c_vec)){
                         dat.new <- dataGen(n = (c_vec[j]-c_vec[j-1])*n, cutpoints = cutpoints, 
                                            hazards = hazards, beta = beta, probs = probs)
                         dat.new <- dataDrop(dat.new, censors = censors, n = (c_vec[j]-c_vec[j-1])*n, cutpoints = cutpoints)
                         dat.acc <- rbind(dat.acc, dat.new)
                         dat.proc <- processData(dat = dat.acc, cutpoints = cutpoints)
                         post.cond <- getConditional(jags_model = jags.path, 
                                                     dat_list = dat.proc, prec_beta = prec_beta, 
                                                     prec_logh = prec_logh, burnin = burnin, nchains = nchains,
                                                     nthin = nthin, ndraws = ndraws*nthin)
                         comps.temp <- getCompProbs(mod.samples = post.cond, dat_vec = dat.acc[,4], deltaL = deltaL)
                         comps <- cbind(comps, comps.temp)
                       }
                       
                       ## return vector of complementary probabilities for each analysis
                       comps
                     }
  
  return(comp.probs)
}

## this function is very similar to the previous one, but it is used to get the
## complementary probabilities used for confirmation purposes
confirmCompProbs <- function(n, c_vec, cutpoints, hazards, beta, probs, censors,
                             deltaL, R = 10000, jags_model, prec_beta, prec_logh, 
                             burnin = 1000, nchains = 1, nthin = 2, ndraws = 5000, 
                             seed = 1){
  
  ## the inputs are described as follows:
  ## n: the sample size for the first interim analysis
  ## c_vec: the vector of constants c to reflect how much more data are 
  ##        collected in the subsequent stages
  ## cutpoints: the vector of interval boundaries for the piecewise
  ##            exponential function
  ## beta: a 3 x 1 vector of coefficients for the hazard function
  ## probs: a 2 x 1 vector for the probability of each binary covariate being 1
  ## censors: a vector of probabilities indicating what proportion of patients
  ##          are censored at the beginning of each interval
  ## deltaL: the lower interval endpoint for the hazard ratio hypothesis
  ## R: the number of simulation repetitions used to estimate each sampling distribution
  ## jags_model: the name of the JAGS file (must be in the working directory)
  ## prec_beta: the prior precisions for the regression coefficients
  ## prec_logh: the prior precisions for the log hazards
  ## burnin: the number of burnin iterations for MCMC
  ## nchains: the number of MCMC chains
  ## nthin: the thinning parameter for MCMC
  ## ndraws: the number of MCMC draws that should be retained
  ## seed: a numeric seed used for random number generation
  
  ## this subfunction generates data before accounting for dropout
  dataGen <- function(n, cutpoints, hazards, beta, probs) {
    
    ## Generate covariates x1 and x2
    x <- cbind(rbinom(n, 1, probs[1]), rbinom(n, 1, probs[2]))
    
    ## Compute interval widths
    intervals <- c(0, cutpoints)
    widths <- diff(intervals)
    
    ## Identify unique covariate combinations
    uc <- unique(x)
    
    ## Precompute adjusted hazards for each unique combination
    adjusted_hazards_list <- lapply(1:nrow(uc), function(i) {
      linpred <- sum(beta * c(uc[i, ], uc[i,1]*uc[i,2]))
      hazards * exp(linpred)
    })
    
    ## Initialize results and generate pseudorandom sequence for
    ## CDF inversion
    y <- rep(0, n)
    event <- rep(0, n)
    u <- runif(n)
    
    for (i in seq_along(adjusted_hazards_list)) {
      ## Select patients belonging to subtreatment group
      group_idx <- which(x[,1] == uc[i,1] & x[,2] == uc[i,2])
      
      ## Retrieve adjusted hazards for this group
      adjusted_hazards <- adjusted_hazards_list[[i]]
      
      ## Simulate survival times for this group
      for (j in group_idx) {
        cumulative_hazard <- 0
        for (k in seq_along(widths)) {
          ## Increment cumulative hazard
          incremental_hazard <- adjusted_hazards[k] * widths[k]
          cumulative_hazard <- cumulative_hazard + incremental_hazard
          
          ## Check if event occurs in this interval based on the CDF and u realization
          if (u[j] > exp(-cumulative_hazard)) {
            y[j] <- intervals[k] - (log(u[j]) + cumulative_hazard - incremental_hazard)/adjusted_hazards[k]
            event[j] <- 1
            break
          }
        }
        
        ## If no event experienced, survival time exceeds last interval
        if (event[j] == 0) {
          y[j] <- max(intervals)
        }
      }
    }
    
    return(data.frame(y = y, event = event, x1 = x[, 1], x2 = x[, 2]))
  }
  
  ## this subfunction censors data to account for dropout
  dataDrop <- function(dat, censors, n, cutpoints){
    ## determine time of dropout (if applicable)
    ## dropout indicators 
    drop.mat <- cbind(matrix(rbinom(n*length(cutpoints), 1, rep(censors, each = n)), 
                             byrow = FALSE, nrow = n), rep(1, n))
    
    drop.vec <- apply(drop.mat, 1, function(x){which(x == 1)[1]})
    
    intervals <- c(0, cutpoints)
    for (k in 1:length(censors)){
      dat$event <- ifelse(drop.vec != k, dat$event, 
                          ifelse(dat$y < intervals[k], dat$event, 0))
      dat$y <- ifelse(drop.vec != k, dat$y, 
                      ifelse(dat$y < intervals[k], dat$y, intervals[k]))
    }
    
    ## round survival times up to the nearest day
    dat$y <- ceiling(dat$y)
    
    return(dat)
  }
  
  ## this subfunction processes the data for the JAGS model
  processData <- function(dat, cutpoints){
    
    ## obtain summary statistics for the data
    ## Compute interval widths
    intervals <- c(0, cutpoints)
    widths <- diff(intervals)
    
    ## get the total number of events by interval segment
    ev <- dat$event
    tm <- dat$y
    e_seg <- sum((tm<=cutpoints[1])*ev)
    for (k in 2:length(cutpoints)){
      e_seg <- c(e_seg, sum((tm > cutpoints[k-1] & tm<=cutpoints[k])*ev))
    }
    
    ## get the total events by covariate combination
    ev_x1 <- sum((dat$x1)*ev)
    ev_x2 <- sum((dat$x2)*ev)
    ev_x1x2 <- sum((dat$x1)*ev*dat$x2)
    e_sub <- c(ev_x1, ev_x2, ev_x1x2)
    
    ## get the total time at risk for x1 = 0 and x2 = 0
    tm00 <- subset(dat, dat$x1 == 0 & dat$x2 == 0)$y
    e_00 <- sum(pmin(tm00, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_00 <- c(e_00, sum(pmin(pmax(tm00 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## get the total time at risk for x1 = 1 and x2 = 0
    tm10 <- subset(dat, dat$x1 == 1 & dat$x2 == 0)$y
    e_10 <- sum(pmin(tm10, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_10 <- c(e_10, sum(pmin(pmax(tm10 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## get the total time at risk for x1 = 0 and x2 = 1
    tm01 <- subset(dat, dat$x1 == 0 & dat$x2 == 1)$y
    e_01 <- sum(pmin(tm01, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_01 <- c(e_01, sum(pmin(pmax(tm01 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## get the total time at risk for x1 = 1 and x2 = 1
    tm11 <- subset(dat, dat$x1 == 1 & dat$x2 == 1)$y
    e_11 <- sum(pmin(tm11, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_11 <- c(e_11, sum(pmin(pmax(tm11 - cutpoints[k-1], 0), widths[k])))
    } 
    
    return(list(e_seg, e_sub, e_00, e_10, e_01, e_11))
  }
  
  ## this subfunction returns conditional posterior draws to be input into the
  ## marginalization process
  
  getConditional <- function(jags_model, dat_list, prec_beta, prec_logh, burnin,
                             nchains, nthin, ndraws){
    
    ## initialize model
    model.fit <- jags.model(file=jags_model,
                            data=list(e_seg = dat_list[[1]], e_sub = dat_list[[2]], 
                                      e_00 = dat_list[[3]], e_10 = dat_list[[4]], 
                                      e_01 = dat_list[[5]], e_11 = dat_list[[6]], K = length(dat_list[[1]]), 
                                      prec_beta = prec_beta, prec_logh = prec_logh, zero = 0), 
                            n.chains = nchains, quiet = TRUE)
    
    ## update models and extract posterior draws
    update(model.fit, burnin, progress.bar = "none")
    model.samples <- coda.samples(model.fit, c("beta[1]", "beta[2]", "beta[3]"), 
                                  n.iter=ndraws, thin=nthin, progress.bar = "none")
    
    as.matrix(model.samples)
  }
  
  ## this subfunction implements the marginalization process for the hazard ratio
  ## and gets the complementary posterior probabilities for the hypothesis of interest
  getCompProbs <- function(mod.samples, dat_vec, deltaL){
    mm <- nrow(mod.samples)
    nn <- length(dat_vec)
    
    ## implement marginalization process
    haz.marg <- NULL
    for (k in 1:mm){
      mu0 <- exp(mod.samples[k, 2]*dat_vec)
      w0 <- rexp(nn, 1)
      w0 <- w0/sum(w0)
      
      mu1 <- exp(mod.samples[k, 1] + mod.samples[k, 3]*dat_vec)*mu0
      w1 <- rexp(nn, 1)
      w1 <- w1/sum(w1)
      
      ## get sample from log-HR
      haz.marg[k] <- log(sum(w1*mu1)) - log(sum(w0*mu0)) 
    }
    
    ## get complementary posterior probability based on posterior samples
    mean(haz.marg <= log(deltaL))
  }
  
  ## we now implement this in parallel across all analyses
  comp.probs <- foreach(i=1:R, .packages=c('rjags', 'coda'), .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## this function adjusts for fractional sample sizes
                          nss <- ceiling(n*c_vec)
                          n.diff <- diff(nss)
                          
                          set.seed(seed + i)
                          dat.acc <- dataGen(n = nss[1], cutpoints = cutpoints, 
                                             hazards = hazards, beta = beta, probs = probs)
                          dat.acc <- dataDrop(dat.acc, censors = censors, n = nss[1], cutpoints = cutpoints)
                          dat.proc <- processData(dat = dat.acc, cutpoints = cutpoints)
                          jags.path <- paste(getwd(), '/', jags_model, sep='')
                          post.cond <- getConditional(jags_model = jags.path, 
                                                      dat_list = dat.proc, prec_beta = prec_beta, 
                                                      prec_logh = prec_logh, burnin = burnin, nchains = nchains,
                                                      nthin = nthin, ndraws = ndraws*nthin)
                          comps <- getCompProbs(mod.samples = post.cond, dat_vec = dat.acc[, 4], deltaL = deltaL)
                          
                          for (j in 2:length(c_vec)){
                            dat.new <- dataGen(n = n.diff[j-1], cutpoints = cutpoints, 
                                               hazards = hazards, beta = beta, probs = probs)
                            dat.new <- dataDrop(dat.new, censors = censors, n = n.diff[j-1], cutpoints = cutpoints)
                            dat.acc <- rbind(dat.acc, dat.new)
                            dat.proc <- processData(dat = dat.acc, cutpoints = cutpoints)
                            post.cond <- getConditional(jags_model = jags.path, 
                                                        dat_list = dat.proc, prec_beta = prec_beta, 
                                                        prec_logh = prec_logh, burnin = burnin, nchains = nchains,
                                                        nthin = nthin, ndraws = ndraws*nthin)
                            comps.temp <- getCompProbs(mod.samples = post.cond, dat_vec = dat.acc[,4], deltaL = deltaL)
                            comps <- cbind(comps, comps.temp)
                          }
                          
                          ## return vector of complementary probabilities for each analysis
                          comps
                        }
  
  return(comp.probs)
}

## this function creates a matrix of stopping probabilities for a range of sample
## sizes based on our linear approximations. The output from this function is used 
## in the paper to construct plots and bootstrap confidence intervals
stop_mat <- function(m1, m2, n1s, n2s, lb, ub, by, gam){
  # m1 is matrix of sequential probs (first sampling distribution estimate)
  # m2 is the matrix of sequential probs (second sampling distribution estimate)
  # n1s are the sample sizes for the first sampling distribution estimate
  # n2s are the sample sizes for the second sampling distribution estimate
  # lb is lower bound for the first analysis
  # ub is upper bound for the first analysis
  # by is the increment for the first analysis
  # gam are the thresholds
  
  ## need negation of logit since we work with complementary probabilities
  l1 <- function(x){
    -1*logit(x)
  }
  
  # get logits for the sequential probs in both sampling distribution estimates
  ls <- apply(m1, 2, l1)
  li <- apply(m2, 2, l1)
  
  # adjust an infinite logits
  for (j in 1:ncol(ls)){
    ls[,j] <- ifelse(ls[,j] == -Inf, min(subset(ls, is.finite(ls[,j]))[,j]) - 1, ls[,j])
    ls[,j] <- ifelse(ls[,j] == Inf, max(subset(ls, is.finite(ls[,j]))[,j]) + 1, ls[,j])
  }
  
  for (j in 1:ncol(li)){
    li[,j] <- ifelse(li[,j] == -Inf, min(subset(li, is.finite(li[,j]))[,j]) - 1, li[,j])
    li[,j] <- ifelse(li[,j] == Inf, max(subset(li, is.finite(li[,j]))[,j]) + 1, li[,j])
  }
  
  # get indexes to combine individual logits later
  ls <- cbind(ls, seq(1, nrow(ls), 1))
  li <- cbind(li, seq(1, nrow(li), 1))
  
  slopes <- NULL
  ints <- NULL
  
  ## construct the slopes separately for each analysis
  for (j in 1:ncol(m1)){
    ls_s <- ls[order(ls[,j]),j]
    li_s <- li[order(li[,j]),j]
    
    l_slope <- (li_s - ls_s)/(n2s[j]-n1s[j])
    l_int <- ls_s - l_slope*n1s[j]
    
    # reorder according to smaller sample size
    l_slope[ls[order(ls[,j]),ncol(ls)]] <- l_slope 
    l_int[ls[order(ls[,j]),ncol(ls)]] <- l_int 
    
    slopes <- cbind(slopes, l_slope)
    ints <- cbind(ints, l_int)
  }
  
  ## create matrix to calculate cumulative stopping probabilities
  samps <- seq(lb,ub,by)
  ratios <- n1s/n1s[1]
  res.mat <- matrix(0, ncol = ncol(m1), nrow = length(samps))
  for (i in 1:length(samps)){
    ## check which logits are larger than first threshold based on
    ## the linear approximations
    j <- 1
    assign(paste0("d", j), ints[,j] + samps[i]*ratios[j]*slopes[,j])
    stop.temp <- get(paste0("d", j)) >= logit(gam[1])
    res.mat[i, 1] <- mean(stop.temp)
    ## repeat process for the other analyses
    for (j in 2:ncol(m1)){
      assign(paste0("d", j), ints[,j] + samps[i]*ratios[j]*slopes[,j])
      
      ## ensure that previously stopped repetitions are still stopped
      stop.temp <- ifelse(stop.temp, 1, 
                          ifelse(get(paste0("d", j)) >= logit(gam[j]), 1, 0))
      
      res.mat[i, j] <- mean(stop.temp)
    }
  }
  
  ## output matrix where rows correspond to sample sizes and columns
  ## correspond to the stopping probability for a given analysis
  lst <- NULL
  for (j in 1:ncol(m1)){
    lst <- c(lst, res.mat[,j])
  }
  
  return(lst)
  
}