## this code file contains the functions used to implement Algorithm 1 for
## the beta example, as well as the functions to run the confirmation
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

## this function is used to get values for the two shape parameters
## given a uniform design prior on the q-quantile
getShapes <- function(R, low = 0.0225, up = 0.0275, 
                      q = 0.99, total = 400, seed = 1){
  
  ## the inputs are described as follows:
  ## R: the number of Monte Carlo iterations
  ## low: the lower endpoint for the design prior
  ## up: the upper endpoint for the design prior
  ## q: the quantile that defines the target of inference
  ## total: the sum of the two shape parameters
  ## seed: a numeric seed used for random number generation
  
  ## generate quantile values from a uniform distribution
  dd <- runif(R)*(up-low) + low
  
  ## use uniroot to get the alpha value
  find_alpha <- function(alpha, total, q, p){
    qbeta(q, alpha, total - alpha) - p
  }
  uuu <- function(lower, upper, total, q, p){
    uniroot(find_alpha, lower = lower, upper = upper, total = total, q = q, p = p)$root
  }
  
  ## implement for each quantile
  alphas <- sapply(dd, uuu, lower = 0.01, upper = 100, total = total, q = q)
  
  ## sort data frame by effect size and return
  res <- cbind(delta = dd, alpha = alphas, beta = 400 - alphas)
  return(res[order(res[,1]), ])
}

## extract the functions from inside getShapes() to obtain
## shape parameters under the conditional approach. The inputs are
## the same as in getShapes()
uuu <- function(lower, upper, total, q, p){
  uniroot(find_alpha, lower = lower, upper = upper, total = total, q = q, p = p)$root
}

find_alpha <- function(alpha, total, q, p){
  qbeta(q, alpha, total - alpha) - p
}

## this function is used to get the sufficient statistics for the beta
## samples; these results are used to parallelize computing the
## sampling distributions of posterior predictive probabilities
getBetaSummaries <- function(n, c_vec, shapes, seed = 1){
  
  ## the inputs are described as follows:
  ## n: sample size for the first analysis
  ## c_vec: the vector dictating the spacing of the analyses
  ## shapes: a matrix of shape parameters, returned by getShapes()
  ## seed: a numeric seed used for random number generation
  
  ## extract information from shapes matrix
  R <- nrow(shapes)
  deltas <- shapes[,1]
  alphas <- shapes[,2]
  betas <- shapes[,3]
  
  ## get sample sizes for analyses
  ns <- ceiling(n*c_vec)
  
  ## we now implement this in parallel across all analyses
  comp.probs <- foreach(i=1:R, .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## generate data and compute sufficient statistics for first analysis
                          set.seed(seed + i)
                          dat <- rbeta(ns[1], alphas[i], betas[i])
                          dat.proc <- c(sum(log(dat)), sum(log(1-dat)), length(dat))
                          dat.save <- c(deltas[i], dat.proc)
                          
                          ## repeat across all analyses
                          for (j in 2:length(c_vec)){
                            dat.new <- rbeta(n = ns[j]-ns[j-1], alphas[i], betas[i])
                            dat.new.proc <- c(sum(log(dat.new)), sum(log(1-dat.new)), length(dat.new))
                            dat.proc <- dat.proc + dat.new.proc
                            dat.save <- c(dat.save, dat.proc)
                          }
                          
                          ## return vector of complementary probabilities for each analysis
                          dat.save
                        }
}

## this function is used to approximate posteriors and get complementary
## posterior and posterior predictive probabilities for the beta example at a particular
## analysis
getCompPred <- function(ni, nf, summaries, deltaU, jags_model, hyper, 
                        burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000, M = 1000){
  
  ## the inputs are described as follows:
  ## ni: the sample size at the current analysis
  ## ni: the sample size at the final analysis
  ## summaries: a matrix corresponding to sufficient statistics for the current analysis
  ## deltaU: the upper interval endpoint for the hypothesis H1
  ## jags_model: the name of the JAGS file (must be in the working directory)
  ## hyper: the upper interval endpoint for the uniform priors on the shape parameters
  ## burnin: the number of burnin iterations for MCMC
  ## nchains: the number of MCMC chains
  ## nthin: the thinning parameter for MCMC
  ## ndraws: the number of MCMC draws that should be retained
  
  ## this function gets posterior draws given the beta data summaries
  getPostBeta <- function(jags_model, dat_vec, hyper, burnin = 1000, 
                          nchains = 1, nthin = 1, ndraws = 5000){
    
    ## extract the sufficient statistics
    sum_logy <- dat_vec[1]
    sum_log1my <- dat_vec[2]
    n <- dat_vec[3]
    
    ## initialize and run the jags model
    model.fit <- jags.model(file=jags_model,
                            data=list(n=n, sum_logy = sum_logy, sum_log1my = sum_log1my,
                                      upper0 = hyper, zero = 0), 
                            n.chains = nchains, quiet = TRUE)
    
    update(model.fit, burnin, progress.bar = "none")
    model.samples <- coda.samples(model.fit, c("alpha", "beta", "phi"), 
                                  n.iter=ndraws, thin=nthin, progress.bar = "none")
    
    ## return the results
    as.matrix(model.samples)
  }
  
  ## this function gets takes the posterior draws and obtains the complement
  ## of the posterior probability that H1 is true
  getCompProbs <- function(mod.samples, deltaU, q0 = 0.99){
    
    qs <- qbeta(q0, mod.samples[,1], mod.samples[,2])
    qs <- log(qs) - log(1 - qs)
    
    ## get complementary posterior probability using kernel density estimation on the logit scale
    kd.qs <- density(na.omit(qs))
    np.prob <- mean(pnorm(log(deltaU) - log(1-deltaU), ifelse(is.finite(qs), qs,
                                                              ifelse(qs > 0,  max(na.omit(qs)) + 1,
                                                                     min(na.omit(qs)) - 1)),
                          kd.qs$bw, lower.tail = FALSE))
    
    ## use a normal approximation if the complementary probability is extremely small
    if (np.prob < 0.00004){
      pnorm(log(deltaU) - log(1-deltaU), mean(na.omit(qs)), sd(na.omit(qs)), lower.tail = FALSE)
    } else {
      np.prob
    }
  }
  
  R <- nrow(summaries)
  
  ## we now implement this in parallel across all Monte Carlo iterations
  comp.probs <- foreach(i=1:R, .packages=c('rjags', 'coda'), .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## extract data summary for current analysis
                          dat.init.proc <- c(summaries[i,1], summaries[i,2], ni)
                          jags.path <- paste(getwd(), '/', jags_model, sep='')
                          post <- getPostBeta(jags_model = jags.path, dat_vec = dat.init.proc, 
                                              hyper = hyper, burnin = burnin, nchains = nchains,
                                              nthin = nthin, ndraws = ndraws*nthin)
                          comps <- getCompProbs(mod.samples = post, deltaU = deltaU)
                          
                          ## save M posterior draws to compute posterior predictive probabilities
                          post.keep <- post[sample(1:nrow(post), M, replace = FALSE),]
                          alpha.sim <- post.keep[,1]
                          beta.sim <- post.keep[,2]
                          
                          ## get comp probs used to construct posterior predictive probabilities
                          for (k in 1:M){
                            dat.new <- rbeta(n = nf - ni, alpha.sim[k], beta.sim[k])
                            dat.new.proc <- c(sum(log(dat.new)), sum(log(1-dat.new)), nf - ni)
                            dat.proc <- dat.init.proc + dat.new.proc
                            post <- getPostBeta(jags_model = jags.path, dat_vec = dat.proc, 
                                                hyper = hyper, burnin = burnin, nchains = nchains,
                                                nthin = nthin, ndraws = ndraws*nthin)
                            comps.temp <- getCompProbs(mod.samples = post, deltaU = deltaU)
                            comps <- c(comps, comps.temp)
                          }
                          
                          ## return vector of complementary probabilities for this analysis
                          ## the first one is the posterior probability for the current analysis;
                          ## the remaining ones will be used to compute posterior predictive probabilities
                          comps
                        }
}

## this function is used to approximate posteriors and get complementary
## posterior probabilities for the beta example at the final analysis (i.e., 
## we don't need to compute posterior predictive probabilities here). The inputs
## are the same as for the getCompPred() function.
getCompPost <- function(n, summaries, deltaU, jags_model, hyper, 
                        burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000, M = 1000){
  
  getPostBeta <- function(jags_model, dat_vec, hyper, burnin = 1000, 
                          nchains = 1, nthin = 1, ndraws = 5000){
    
    sum_logy <- dat_vec[1]
    sum_log1my <- dat_vec[2]
    n <- dat_vec[3]
    
    model.fit <- jags.model(file=jags_model,
                            data=list(n=n, sum_logy = sum_logy, sum_log1my = sum_log1my,
                                      upper0 = hyper, zero = 0), 
                            n.chains = nchains, quiet = TRUE)
    
    update(model.fit, burnin, progress.bar = "none")
    model.samples <- coda.samples(model.fit, c("alpha", "beta", "phi"), 
                                  n.iter=ndraws, thin=nthin, progress.bar = "none")
    
    as.matrix(model.samples)
  }
  
  getCompProbs <- function(mod.samples, deltaU, q0 = 0.99){
    
    qs <- qbeta(q0, mod.samples[,1], mod.samples[,2])
    qs <- log(qs) - log(1 - qs)
    
    ## get complementary posterior probability
    kd.qs <- density(na.omit(qs))
    np.prob <- mean(pnorm(log(deltaU) - log(1-deltaU), ifelse(is.finite(qs), qs,
                                                              ifelse(qs > 0,  max(na.omit(qs)) + 1,
                                                                     min(na.omit(qs)) - 1)),
                          kd.qs$bw, lower.tail = FALSE))
    if (np.prob < 0.00004){
      pnorm(log(deltaU) - log(1-deltaU), mean(na.omit(qs)), sd(na.omit(qs)), lower.tail = FALSE)
    } else {
      np.prob
    }
  }
  
  R <- nrow(summaries)
  
  ## we now implement this in parallel across all analyses
  comp.probs <- foreach(i=1:R, .packages=c('rjags', 'coda'), .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## initial analysis
                          dat.init.proc <- c(summaries[i,1], summaries[i,2], n)
                          jags.path <- paste(getwd(), '/', jags_model, sep='')
                          post <- getPostBeta(jags_model = jags.path, dat_vec = dat.init.proc, 
                                              hyper = hyper, burnin = burnin, nchains = nchains,
                                              nthin = nthin, ndraws = ndraws*nthin)
                          comps <- getCompProbs(mod.samples = post, deltaU = deltaU)
                          
                          ## return vector of complementary probabilities for this analysis
                          comps
                        }
}

## this function is used to get the estimate posterior predictive probabilities given the
## posterior probabilities returned by getCompPred().
ppp.dens.logit <- function(x, gam){
  
  ## the inputs are described as follows:
  ## x: a vector of posterior probabilities from the posterior predictive distribution
  ## gam: the success threshold for the final analysis
  
  ## approximate this probability using kernel density approximation
  kd <- density(na.omit(logit(as.numeric(x))))
  np.prob <- mean(pnorm(log(1 - gam) - log(gam), ifelse(is.finite(logit(as.numeric(x))), logit(as.numeric(x)),
                                                        ifelse(logit(as.numeric(x)) > 0,  max(na.omit(logit(as.numeric(x)))) + 1,
                                                               min(na.omit(logit(as.numeric(x)))) - 1)),
                        kd$bw, lower.tail = FALSE))
  
  ## use normal approximation for stability if the complementary probability is extremely small
  if (np.prob < 0.00004){
    pnorm(log(1 - gam) - log(gam), mean(na.omit(logit(x))), sd(na.omit(logit(x))), lower.tail = FALSE)
  } else {
    np.prob
  }
}

## this function is used to get the estimate posterior predictive probabilities given the
## posterior probabilities returned by getCompPred(). This function is used in the
## confirmatory simulations where the probabilities are estimated using empirical averages.
## The inputs are the same as ppp.dens.logit().
ppp.prob <- function(x, gam){
  mean(x >= 1 - gam)
}

## this function is used to estimate the joint sampling distribution of posterior and 
## posterior predictive probabilities across all analyses
jointSamp <- function(n, c_vec, ll, gam, J = 10){
  
  ## the inputs are described as follows:
  ## n: sample size for the first analysis
  ## c_vec: the vector dictating the spacing of the analyses
  ## ll: a string used to index the .csv files (i.e., "pred" or "cond")
  ## gam: the success threshold for the final analysis
  ## J: the number of .csv files that the results are split into
  
  res.full <- NULL
  ## repeat for each analysis
  for (i in 1:length(c_vec)){
    res.t <- NULL
    ## repeat for each file
    for (j in 1:J){
      res.read <- read.csv(paste0(ll, "_n_",n,"_t_", i, "_batch_", j, ".csv"))
      ## extract vector of posterior probabilities (for observed data)
      res.temp <- matrix(res.read[, 1], ncol = 1)
      ## for all analyses but the final one, compute posterior predictive probability
      ## using the remaining probabilities
      if (i < length(c_vec)){
        res.temp <- cbind(res.temp, apply(res.read[, -1], 1, ppp.dens.logit, gam = gam))
      }
      res.t <- rbind(res.t, res.temp)
    }
    res.full <- cbind(res.full, res.t)
  }
  return(res.full)
}

## this function is used to create linear approximations to posterior and posterior predictive 
## probabilities under the conditional approach
getLines <- function(m1, m2, n1s, n2s){
  
  ## the inputs are described as follows:
  # m1 is matrix of sequential probs (first joint sampling distribution estimate)
  # m2 is the matrix of sequential probs (second joint sampling distribution estimate)
  # n1s are the sample size for the first joint sampling distribution estimate
  # n2s are the sample size for the second joint sampling distribution estimate
  
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
  
  return(cbind(ints, slopes))
}

## this function is used to create linear approximations to posterior and posterior predictive 
## probabilities under the predictive approach. It uses getLines() to create linear approximation
## after the posterior summaries have been binned into subgroups based on the true value of the 
## target of inference used to generate the data. Remember that the posterior summaries are
## pre-sorted by the true value of the target of inference. The inputs for this function are
## the same as getLines(), except for J, which denotes the number of subgroups.
getLinesPred <- function(m1, m2, n1s, n2s, J){
  lines.res <- NULL
  n.inc <- nrow(m1)/J
  
  ## contstruct linear approximations separately for each subgroup
  for (j in 1:J){
    mat1.temp <- m1[seq((j-1)*n.inc+1, j*n.inc, 1),]
    mat2.temp <- m2[seq((j-1)*n.inc+1, j*n.inc, 1),]
    lines.temp <- getLines(mat1.temp, mat2.temp, n1s, n2s)
    lines.res <- rbind(lines.res, lines.temp)
  }
  return(lines.res)
}

## this function uses linear approximations to create a matrix of stopping probabilities across 
## a range of sample sizes that can be plotted or used to find optimal sample sizes;
## the output is slightly different for the bootstrapping procedure as described below
stopMat <- function(lines, c_vec, lb, ub, by, gam, xi, boot = FALSE){
  
  ## the inputs are described as follows:
  ## lines: a matrix of intercepts and slopes returned by getLines()
  ## c_vec: the vector dictating the spacing of the analyses
  ## lb: the smallest sample size (n1) for which the stopping probs are estimated
  ## ub: the largest sample size (n1) for which the stopping probs are estimated
  ## by: the increments between lb and ub at which the stopping probs are estimated
  ## gam: the vector of sucess thresholds for all analyses
  ## xi: the vector of failure thresholds for the first T-1 analyses
  ## boot: a binary indicator (FALSE if output is not being used for bootstrapping)
  
  ## get vector of sample sizes and ratios
  samps <- seq(lb,ub,by)
  ratios <- c_vec/c_vec[1]
  
  ## extract the number of decisions
  tps <- 0.5*ncol(lines)

  ## extract intercepts and slopes from the lines matrix
  res.mat <- matrix(0, ncol = tps, nrow = length(samps))
  ints <- lines[, seq(1, tps, 1)]
  slopes <- lines[, seq(tps + 1, 2*tps, 1)]
  
  ## repeat process for each sample size
  for (i in 1:length(samps)){
    
    ## check which logits are larger than first success threshold based on
    ## the linear approximations
    j <- 1
    stop.p <- ints[,j] + samps[i]*ratios[j]*slopes[,j] >= logit(gam[j])
    res.mat[i, j] <- mean(stop.p)
    
    ## repeat process for failure thresholds at first analysis
    stop.f <- ints[,j + 1] + samps[i]*ratios[j + 1]*slopes[,j + 1] < logit(xi[j])
    ## make sure we didn't just stop for success
    stop.f <- ifelse(stop.p, 0, stop.f)
    res.mat[i, j+1] <- mean(stop.f)
    
    ## repeat this process for the other interim analyses
    for (j in 2:(0.5*length(c_vec) - 0.5)){
      ## check success thresholds
      stop.p.temp <- ints[,2*j-1] + samps[i]*ratios[2*j-1]*slopes[,2*j-1] >= logit(gam[j])
      stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
      res.mat[i, 2*j-1] <- mean(stop.p)
      
      ## repeat process for failure thresholds
      stop.f.temp <- ints[,2*j] + samps[i]*ratios[2*j]*slopes[,2*j] < logit(xi[j])
      stop.f <- ifelse(stop.p, 0, ifelse(stop.f, 1, stop.f.temp))
      res.mat[i, 2*j] <- mean(stop.f)
    }
    ## implement process for final analysis (success thresholds only)
    j <- 0.5*length(c_vec) + 0.5
    stop.p.temp <- ints[,2*j-1] + samps[i]*ratios[2*j-1]*slopes[,2*j-1] >= logit(gam[j])
    stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
    res.mat[i, 2*j-1] <- mean(stop.p)
    
  }
  
  ## if not bootstrapping, return a matrix where rows correspond to n1 sample sizes
  ## and columns correspond to decisions
  if (boot == FALSE){
    return(res.mat)
  } else {
    
    ## otherwise, return a matrix that is collapsed into a single row
    lst <- NULL
    for (j in 1:tps){
      lst <- c(lst, res.mat[,j])
    }
    return(lst)
  }
}

## this function estimates the stopping probabilities for a single sample
## size (n1) using an estimate of the joint sampling distribution (i.e., no
## linear approximations)
stopData <- function(joint, gam, xi){
  
  ## the inputs are described as follows:
  ## joint: a matrix of complementary probabilities from the joint sampling dist
  ## gam: the vector of sucess thresholds for all analyses
  ## xi: the vector of failure thresholds for the first T-1 analyses
  
  ## get the number of decisions
  tps <- 0.5*ncol(joint) + 0.5
  res <- NULL
  
  ## convert the complementary probabilities to desired probabilities
  joint <- 1 - joint
  
  ## check which logits are larger than first threshold based on data
  j <- 1
  stop.p <- joint[,j] >= gam[j]
  res[j] <- mean(stop.p)
  ## repeat process for failure thresholds
  stop.f <- joint[, j+1] < xi[j]
  stop.f <- ifelse(stop.p, 0, stop.f)
  res[j+1] <- mean(stop.f)
  
  ## repeat process for other interim analyses
  for (j in 2:(tps - 1)){
    stop.p.temp <- joint[,2*j-1] >= gam[j]
    stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
    res[2*j-1] <- mean(stop.p)
    ## repeat process for failure thresholds
    stop.f.temp <- joint[,2*j] < xi[j]
    stop.f <- ifelse(stop.p, 0, ifelse(stop.f, 1, stop.f.temp))
    res[2*j] <- mean(stop.f)
  }
  
  ## implement process for final analysis (stopping only for success)
  j <- tps
  stop.p.temp <- joint[,2*j-1] >= gam[j]
  stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
  res[2*j-1] <- mean(stop.p)
  
  return(res)
}

## this function is used to estimate the joint sampling distribution of posterior and 
## posterior predictive probabilities across all analyses. The difference between this 
## function and jointSamp() is that ppp.prob() is used to compute posterior predictive 
## probabilities instead of ppp.dens.logit(). The inputs are the same as for jointSamp()
jointSampDat <- function(n, c_vec, ll, gam){
  res.full <- NULL
  for (i in 1:length(c_vec)){
    res.t <- NULL
    for (j in 1:10){
      res.read <- read.csv(paste0(ll, "_n_",n,"_t_", i, "_batch_", j, ".csv"))
      res.temp <- matrix(res.read[, 1], ncol = 1)
      if (i < length(c_vec)){
        res.temp <- cbind(res.temp, apply(res.read[, -1], 1, ppp.prob, gam = gam))
      }
      res.t <- rbind(res.t, res.temp)
    }
    res.full <- cbind(res.full, res.t)
  }
  return(res.full)
}

## this function is the same as getCompPred() except that the posterior probabilities are 
## estimated using empirical averages (used for confirmation purposes). The inputs are the 
## same as getCompPred().
getCompPredMean <- function(ni, nf, summaries, deltaU, jags_model, hyper, 
                        burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000, M = 1000){
  
  getPostBeta <- function(jags_model, dat_vec, hyper, burnin = 1000, 
                          nchains = 1, nthin = 1, ndraws = 5000){
    
    sum_logy <- dat_vec[1]
    sum_log1my <- dat_vec[2]
    n <- dat_vec[3]
    
    model.fit <- jags.model(file=jags_model,
                            data=list(n=n, sum_logy = sum_logy, sum_log1my = sum_log1my,
                                      upper0 = hyper, zero = 0), 
                            n.chains = nchains, quiet = TRUE)
    
    update(model.fit, burnin, progress.bar = "none")
    model.samples <- coda.samples(model.fit, c("alpha", "beta", "phi"), 
                                  n.iter=ndraws, thin=nthin, progress.bar = "none")
    
    as.matrix(model.samples)
  }
  
  getCompProbs <- function(mod.samples, deltaU, q0 = 0.99){
    
    qs <- qbeta(q0, mod.samples[,1], mod.samples[,2])
    mean(qs > deltaU)
  }
  
  R <- nrow(summaries)
  
  ## we now implement this in parallel across all analyses
  comp.probs <- foreach(i=1:R, .packages=c('rjags', 'coda'), .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## initial analysis
                          dat.init.proc <- c(summaries[i,1], summaries[i,2], ni)
                          jags.path <- paste(getwd(), '/', jags_model, sep='')
                          post <- getPostBeta(jags_model = jags.path, dat_vec = dat.init.proc, 
                                              hyper = hyper, burnin = burnin, nchains = nchains,
                                              nthin = nthin, ndraws = ndraws*nthin)
                          comps <- getCompProbs(mod.samples = post, deltaU = deltaU)
                          
                          post.keep <- post[sample(1:nrow(post), M, replace = FALSE),]
                          alpha.sim <- post.keep[,1]
                          beta.sim <- post.keep[,2]
                          
                          ## get comp probs used to construct posterior predictive probs
                          for (k in 1:M){
                            dat.new <- rbeta(n = nf - ni, alpha.sim[k], beta.sim[k])
                            dat.new.proc <- c(sum(log(dat.new)), sum(log(1-dat.new)), nf - ni)
                            dat.proc <- dat.init.proc + dat.new.proc
                            post <- getPostBeta(jags_model = jags.path, dat_vec = dat.proc, 
                                                hyper = hyper, burnin = burnin, nchains = nchains,
                                                nthin = nthin, ndraws = ndraws*nthin)
                            comps.temp <- getCompProbs(mod.samples = post, deltaU = deltaU)
                            comps <- c(comps, comps.temp)
                          }
                          
                          ## return vector of complementary probabilities for this analysis
                          comps
                        }
}

## this function is the same as getCompPost() except that the posterior probabilities are 
## estimated using empirical averages (used for confirmation purposes). The inputs are the 
## same as getCompPost(). This function is again used for the final analysis.
getCompPostMean <- function(n, summaries, deltaU, jags_model, hyper, 
                        burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000, M = 1000){
  
  getPostBeta <- function(jags_model, dat_vec, hyper, burnin = 1000, 
                          nchains = 1, nthin = 1, ndraws = 5000){
    
    sum_logy <- dat_vec[1]
    sum_log1my <- dat_vec[2]
    n <- dat_vec[3]
    
    model.fit <- jags.model(file=jags_model,
                            data=list(n=n, sum_logy = sum_logy, sum_log1my = sum_log1my,
                                      upper0 = hyper, zero = 0), 
                            n.chains = nchains, quiet = TRUE)
    
    update(model.fit, burnin, progress.bar = "none")
    model.samples <- coda.samples(model.fit, c("alpha", "beta", "phi"), 
                                  n.iter=ndraws, thin=nthin, progress.bar = "none")
    
    as.matrix(model.samples)
  }
  
  getCompProbs <- function(mod.samples, deltaU, q0 = 0.99){
    
    qs <- qbeta(q0, mod.samples[,1], mod.samples[,2])
    mean(qs > deltaU)
  }
  
  R <- nrow(summaries)
  
  ## we now implement this in parallel across all analyses
  comp.probs <- foreach(i=1:R, .packages=c('rjags', 'coda'), .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## initial analysis
                          dat.init.proc <- c(summaries[i,1], summaries[i,2], n)
                          jags.path <- paste(getwd(), '/', jags_model, sep='')
                          post <- getPostBeta(jags_model = jags.path, dat_vec = dat.init.proc, 
                                              hyper = hyper, burnin = burnin, nchains = nchains,
                                              nthin = nthin, ndraws = ndraws*nthin)
                          comps <- getCompProbs(mod.samples = post, deltaU = deltaU)
                          
                          ## return vector of complementary probabilities for this analysis
                          comps
                        }
}