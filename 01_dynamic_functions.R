## this code file contains the functions used to implement Algorithm 1 for
## the dynamic borrowing example in Section 6.1 of the main text

## define expit and logit functions
expit <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(x){
  return(log(x) - log(1-x))
}

## this function is used to obtain the posterior weights for the mixture distribution
## the inputs are as follows:
## y: the number of successes (in one intervention)
## n: the number of observations (in one intervention)
## a: the first shape parameter for the informative prior
## b: the second shape parameter for the informative prior
## a0: the first shape parameter for the diffuse prior
## b0: the second shape paramter for the diffuse prior
wpost = function(y, n, a, b, a0 = 1, b0 = 1) {
  lp = lbeta(a + y, b + n - y) - lbeta(a, b)
  lp0 = lbeta(a0 + y, b0 + n - y) - lbeta(a0, b0)
  lp_max = max(lp, lp0)
  p = exp(lp - lp_max)
  p0 = exp(lp0 - lp_max)
  wp = p/(p + p0)
  return(wp)
}

## this is the function to sample from the mixture posterior distribution for
## one intervention; the inputs are as follows:
## w: the posterior weight for the informative prior (calculated using wpost() function)
## a: the first shape parameter for the informative prior
## b: the second shape parameter for the informative prior 
## y: the number of successes (in one intervention)
## n: the number of observations (in one intervention)
## M: number of observations simulated from the mixture prior
betaFMsamp <- function(w, a, b, y, n, M = 1000) {
  idx <- rbinom(M, 1, w) + 1  # now 1 or 2
  as <- c(a, 1)[idx] + y
  bs <- c(b, 1)[idx] + n - y
  rbeta(M, as, bs)
}

## this function is used to get the sufficient statistics for the binomial
## samples; these results are used to parallelize computing the
## sampling distributions of posterior predictive probabilities
getBinSummaries <- function(n, c_vec, ratios, probs, R, seed = 1){
  
  ## the inputs are described as follows:
  ## n: sample size for the first analysis
  ## c_vec: the vector dictating the spacing of the analyses
  ## ratios: proportion allocation to treatment in each stage
  ## probs: a vector of success probabilities
  ## R: number of simulation repetitions
  ## seed: a numeric seed used for random number generation
  
  ## get sample sizes for analyses
  ns <- ceiling(n*c_vec)
  ndiff <- c(ns[1], diff(ns))
  
  ## we now implement this in parallel across all analyses
  dat.summ <- foreach(i=1:R, .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## generate data and compute sufficient statistics for first analysis
                          set.seed(seed + i)
                          
                          n1diff <- rbinom(length(c_vec), ndiff, 1 - ratios)
                          n2diff <- ndiff - n1diff
                          
                          dat.proc <- c(n1diff[1], rbinom(1, n1diff[1], probs[1]), n2diff[1], rbinom(1, n2diff[1], probs[2]))
                          dat.save <- dat.proc
                          
                          ## repeat across all analyses
                          for (j in 2:length(c_vec)){
                            dat.new.proc <- c(n1diff[j], rbinom(1, n1diff[j], probs[1]), n2diff[j], rbinom(1, n2diff[j], probs[2]))
                            dat.proc <- dat.proc + dat.new.proc
                            dat.save <- c(dat.save, dat.proc)
                          }
                          
                          ## return vector of complementary probabilities for each analysis
                          dat.save
                        }
  
  return(dat.summ)
}

## same function as getBinSummaries() used under the predictive approach 
## in design scenario 4
getBinSummariesPred <- function(n, c_vec, ratios, probs, R, seed = 1){
  
  ## the inputs are described as follows:
  ## n: sample size for the first analysis
  ## c_vec: the vector dictating the spacing of the analyses
  ## ratios: proportion allocation to treatment in each stage
  ## probs: a vector of success probabilities
  ## R: number of simulation repetitions
  ## seed: a numeric seed used for random number generation
  
  ## get sample sizes for analyses
  ns <- ceiling(n*c_vec)
  ndiff <- c(ns[1], diff(ns))
  
  ## we now implement this in parallel across all analyses
  dat.summ <- foreach(i=1:R, .combine=rbind,
                      .options.snow=opts) %dopar% {
                        
                        ## generate data and compute sufficient statistics for first analysis
                        set.seed(seed + i)
                        
                        n1diff <- rbinom(length(c_vec), ndiff, 1 - ratios)
                        n2diff <- ndiff - n1diff
                        
                        dat.proc <- c(n1diff[1], rbinom(1, n1diff[1], probs[i,1]), n2diff[1], rbinom(1, n2diff[1], probs[i,2]))
                        dat.save <- dat.proc
                        
                        ## repeat across all analyses
                        for (j in 2:length(c_vec)){
                          dat.new.proc <- c(n1diff[j], rbinom(1, n1diff[j], probs[i,1]), n2diff[j], rbinom(1, n2diff[j], probs[i,2]))
                          dat.proc <- dat.proc + dat.new.proc
                          dat.save <- c(dat.save, dat.proc)
                        }
                        
                        ## return vector of complementary probabilities for each analysis
                        dat.save
                      }
  
  return(dat.summ)
}

## this function is used to approximate posteriors and get complementary
## posterior and posterior predictive probabilities for the binomial example at a particular
## analysis
getCompPred <- function(ni, nf, summaries, deltaL = 0, mappar, ratio = 0.5,
                        ndraws = 5000, M = 1000){
  
  ## the inputs are described as follows:
  ## ni: the sample size at the current analysis
  ## ni: the sample size at the final analysis
  ## summaries: a matrix corresponding to sufficient statistics for the current analysis
  ## deltaL: the upper interval endpoint for the hypothesis H1
  ## ratio: treatment allocation ratio for future stages
  ## mappar: (a, b) parameters for MAP prior
  ## ndraws: the number of MCMC draws that should be retained
  ## M: number of predictive samples for posterior predictive probabilities
  
  wpost = function(y, n, a, b, a0 = 1, b0 = 1) {
    lp = lbeta(a + y, b + n - y) - lbeta(a, b)
    lp0 = lbeta(a0 + y, b0 + n - y) - lbeta(a0, b0)
    lp_max = max(lp, lp0)
    p = exp(lp - lp_max)
    p0 = exp(lp0 - lp_max)
    wp = p/(p + p0)
    return(wp)
  }
  
  betaFMsamp <- function(w, a, b, y, n, M = 1000) {
    idx <- rbinom(M, 1, w) + 1  # now 1 or 2
    as <- c(a, 1)[idx] + y
    bs <- c(b, 1)[idx] + n - y
    rbeta(M, as, bs)
  }
  
  ## this function gets posterior draws given the binomial data summaries
  getPostBin <- function(dat_vec, mappar, ndraws = 5000){
    
    ## extract the sufficient statistics and priors
    a1 <- as.numeric(mappar[1])
    b1 <- as.numeric(mappar[2])
    
    y0 <- as.numeric(dat_vec[2])
    n0 <- as.numeric(dat_vec[1])
    y1 <- as.numeric(dat_vec[4])
    n1 <- as.numeric(dat_vec[3])
    
    w0 <- wpost(y0, n0, a1, b1)
    s0 <- betaFMsamp(w0, a1, b1, y0, n0, M = ndraws)
    
    s1 <- rbeta(ndraws, 1 + y1, 1 + n1 - y1)
    
    ## return the results
    cbind(s0, s1)
  }
  
  ## this commented out function was used for the diffuse priors in scenario 3
  # getPostBin <- function(dat_vec, mappar, ndraws = 5000){
  #   
  #   ## extract the sufficient statistics and priors
  #   a1 <- as.numeric(mappar[1])
  #   b1 <- as.numeric(mappar[2])
  #   
  #   y0 <- as.numeric(dat_vec[2])
  #   n0 <- as.numeric(dat_vec[1])
  #   y1 <- as.numeric(dat_vec[4])
  #   n1 <- as.numeric(dat_vec[3])
  # 
  #   s0 <- rbeta(ndraws, 1 + y0, 1 + n0 - y0)
  #   s1 <- rbeta(ndraws, 1 + y1, 1 + n1 - y1)
  #   
  #   ## return the results
  #   cbind(s0, s1)
  # }
  
  ## this function gets takes the posterior draws and obtains the complement
  ## of the posterior probability that H1 is true
  estimate_pp <- function(treatment, control, lmargin) {
    ldiff <- log(treatment - control + 1) - log(1 - (treatment - control))
    kd <- density(na.omit(ldiff))
    prob <- mean(pnorm(lmargin, ifelse(is.finite(ldiff), ldiff, max(na.omit(ldiff)) + 1), kd$bw, lower.tail = TRUE))
    if (prob < 0.00004) prob <- pnorm(lmargin, mean(na.omit(ldiff)), sd(na.omit(ldiff)), lower.tail = TRUE)
    return(prob)
  }
  
  R <- nrow(summaries)
  
  ## we now implement this in parallel across all Monte Carlo iterations
  comp.probs <- foreach(i=1:R, .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## extract data summary for current analysis
                          dat.init.proc <- summaries[i,]
                          post <- getPostBin(dat_vec = dat.init.proc, mappar = mappar, ndraws = ndraws)
                          comps <- estimate_pp(treatment = post[,2], control = post[,1], lmargin = deltaL)
                          
                          ## save M posterior draws to compute posterior predictive probabilities
                          post.keep <- post[sample(1:nrow(post), M, replace = FALSE),]
                          p0.sim <- post.keep[,1]
                          p1.sim <- post.keep[,2]
                          
                          n.new <- ceiling(nf - ni)
                          n1diff.new <- rbinom(M, n.new, 1 - ratio)
                          n2diff.new <- n.new - n1diff.new
                          ## get comp probs used to construct posterior predictive probabilities
                          for (k in 1:M){
                            dat.new.proc <- c(n1diff.new[k], rbinom(1, n1diff.new[k], p0.sim[k]), 
                                              n2diff.new[k], rbinom(1, n2diff.new[k], p1.sim[k]))
                            dat.proc <- dat.init.proc + dat.new.proc
                            post <- getPostBin(dat_vec = dat.proc, mappar = mappar, ndraws = ndraws)
                            comps.temp <- estimate_pp(treatment = post[,2], control = post[,1], lmargin = deltaL)
                            comps <- c(comps, comps.temp)
                          }
                          
                          ## return vector of complementary probabilities for this analysis
                          ## the first one is the posterior probability for the current analysis;
                          ## the remaining ones will be used to compute posterior predictive probabilities
                          comps
                        }
}

## this function is used to approximate posteriors and get complementary
## posterior probabilities for the binomial example at the final analysis (i.e., 
## we don't need to compute posterior predictive probabilities here). The inputs
## are the same as for the getCompPred() function.
getCompPost <- function(summaries, deltaL = 0, mappar, ratio = 0.5,
                        ndraws = 5000, M = 1000){
  
  ## the inputs are described as follows:
  ## summaries: a matrix corresponding to sufficient statistics for the current analysis
  ## deltaL: the upper interval endpoint for the hypothesis H1
  ## ratio: treatment allocation ratio for future stages
  ## ndraws: the number of MCMC draws that should be retained
  ## M: number of predictive samples for posterior predictive probabilities
  
  wpost = function(y, n, a, b, a0 = 1, b0 = 1) {
    lp = lbeta(a + y, b + n - y) - lbeta(a, b)
    lp0 = lbeta(a0 + y, b0 + n - y) - lbeta(a0, b0)
    lp_max = max(lp, lp0)
    p = exp(lp - lp_max)
    p0 = exp(lp0 - lp_max)
    wp = p/(p + p0)
    return(wp)
  }
  
  betaFMsamp <- function(w, a, b, y, n, M = 1000) {
    idx <- rbinom(M, 1, w) + 1  # now 1 or 2
    as <- c(a, 1)[idx] + y
    bs <- c(b, 1)[idx] + n - y
    rbeta(M, as, bs)
  }
  
  ## this function gets posterior draws given the binomial data summaries
  getPostBin <- function(dat_vec, mappar, ndraws = 5000){
    
    ## extract the sufficient statistics and priors
    a1 <- as.numeric(mappar[1])
    b1 <- as.numeric(mappar[2])
    
    y0 <- as.numeric(dat_vec[2])
    n0 <- as.numeric(dat_vec[1])
    y1 <- as.numeric(dat_vec[4])
    n1 <- as.numeric(dat_vec[3])
    
    w0 <- wpost(y0, n0, a1, b1)
    s0 <- betaFMsamp(w0, a1, b1, y0, n0, M = ndraws)
    
    s1 <- rbeta(ndraws, 1 + y1, 1 + n1 - y1)
    
    ## return the results
    cbind(s0, s1)
  }
  
  ## this commented out function was used for the diffuse priors in scenario 3
  # getPostBin <- function(dat_vec, mappar, ndraws = 5000){
  #   
  #   ## extract the sufficient statistics and priors
  #   a1 <- as.numeric(mappar[1])
  #   b1 <- as.numeric(mappar[2])
  #   
  #   y0 <- as.numeric(dat_vec[2])
  #   n0 <- as.numeric(dat_vec[1])
  #   y1 <- as.numeric(dat_vec[4])
  #   n1 <- as.numeric(dat_vec[3])
  # 
  #   s0 <- rbeta(ndraws, 1 + y0, 1 + n0 - y0)
  #   s1 <- rbeta(ndraws, 1 + y1, 1 + n1 - y1)
  #   
  #   ## return the results
  #   cbind(s0, s1)
  # }
  
  ## this function gets takes the posterior draws and obtains the complement
  ## of the posterior probability that H1 is true
  estimate_pp <- function(treatment, control, lmargin) {
    ldiff <- log(treatment - control + 1) - log(1 - (treatment - control))
    kd <- density(na.omit(ldiff))
    prob <- mean(pnorm(lmargin, ifelse(is.finite(ldiff), ldiff, max(na.omit(ldiff)) + 1), kd$bw, lower.tail = TRUE))
    if (prob < 0.00004) prob <- pnorm(lmargin, mean(na.omit(ldiff)), sd(na.omit(ldiff)), lower.tail = TRUE)
    return(prob)
  }
  
  R <- nrow(summaries)
  
  ## we now implement this in parallel across all Monte Carlo iterations
  comp.probs <- foreach(i=1:R, .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ## extract data summary for current analysis
                          dat.init.proc <- summaries[i,]
                          post <- getPostBin(dat_vec = dat.init.proc, mappar = mappar, ndraws = ndraws)
                          comps <- estimate_pp(treatment = post[,2], control = post[,1], lmargin = deltaL)
                          
                          ## return vector of complementary probabilities for this analysis
                          ## the first one is the posterior probability for the current analysis;
                          ## the remaining ones will be used to compute posterior predictive probabilities
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

## this function is used to tune the success thresholds before accounting for failure
## note: this function can also be used to tune failure thresholds (just switch 
## m1 and m10 as well as m2 and m20 as inputs, and let alp be (beta_1, ..., beta_{T-1})
## from a beta-spending function) 
stop_mat <- function(m1, m2, m10, m20, 
                     n1s, n2s, lb, ub, by, alp){
  # m1 is matrix of sequential probs (first sampling distribution estimate)
  # m2 is the matrix of sequential probs (second sampling distribution estimate)
  # m10 is matrix of sequential probs under Psi0 (first sampling distribution estimate)
  # m20 is the matrix of sequential probs under Psi0 (second sampling distribution estimate)
  # n1s are the sample sizes for the first sampling distribution estimate
  # n2s are the sample sizes for the second sampling distribution estimate
  # lb is lower bound for the first analysis
  # ub is upper bound for the first analysis
  # by is the increment for the first analysis
  # alp is the alpha-spending function (vector)
  
  ## need negation of logit since we work with complementary probabilities
  l1 <- function(x){
    -1*logit(x)
  }
  
  ## for hypothesis H1
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
  
  ## for hypothesis H0
  # get logits for the sequential probs in both sampling distribution estimates
  ls0 <- apply(m10, 2, l1)
  li0 <- apply(m20, 2, l1)
  
  # adjust an infinite logits
  for (j in 1:ncol(ls0)){
    ls0[,j] <- ifelse(ls0[,j] == -Inf, min(subset(ls0, is.finite(ls0[,j]))[,j]) - 1, ls0[,j])
    ls0[,j] <- ifelse(ls0[,j] == Inf, max(subset(ls0, is.finite(ls0[,j]))[,j]) + 1, ls0[,j])
  }
  
  for (j in 1:ncol(li0)){
    li0[,j] <- ifelse(li0[,j] == -Inf, min(subset(li0, is.finite(li0[,j]))[,j]) - 1, li0[,j])
    li0[,j] <- ifelse(li0[,j] == Inf, max(subset(li0, is.finite(li0[,j]))[,j]) + 1, li0[,j])
  }
  
  # get indexes to combine individual logits later
  ls0 <- cbind(ls0, seq(1, nrow(ls0), 1))
  li0 <- cbind(li0, seq(1, nrow(li0), 1))
  
  slopes0 <- NULL
  ints0 <- NULL
  
  ## construct the slopes separately for each analysis
  for (j in 1:ncol(m10)){
    ls_s0 <- ls0[order(ls0[,j]),j]
    li_s0 <- li0[order(li0[,j]),j]
    
    l_slope0 <- (li_s0 - ls_s0)/(n2s[j]-n1s[j])
    l_int0 <- ls_s0 - l_slope0*n1s[j]
    
    # reorder according to smaller sample size
    l_slope0[ls0[order(ls0[,j]),ncol(ls0)]] <- l_slope0 
    l_int0[ls0[order(ls0[,j]),ncol(ls0)]] <- l_int0 
    
    slopes0 <- cbind(slopes0, l_slope0)
    ints0 <- cbind(ints0, l_int0)
  }
  
  ## create matrix to calculate cumulative stopping probabilities
  samps <- seq(lb,ub,by)
  ratios <- n1s/n1s[1]
  res.mat <- matrix(0, ncol = ncol(m1), nrow = length(samps))
  res.mat0 <- matrix(0, ncol = ncol(m10), nrow = length(samps))
  gam.mat <- matrix(0, ncol = length(alp), nrow = length(samps))
  for (i in 1:length(samps)){
    print(samps[i])
    ## find appropriate decision threshold (analysis 1)
    j <- 1
    assign(paste0("d0", j), ints0[,j] + samps[i]*ratios[j]*slopes0[,j])
    
    gam.temp.lb <- 0
    gam.temp.ub <- 1
    while(gam.temp.ub - gam.temp.lb > 0.00010001){
      gam.temp <- round(0.5*(gam.temp.lb + gam.temp.ub), 4)
      stop.bs <- mean(get(paste0("d0", j)) >= logit(gam.temp))
      if (stop.bs > alp[1]){
        gam.temp.lb <- gam.temp
      } else {
        gam.temp.ub <- gam.temp
      }
    }
    
    assign(paste0("gam", j), gam.temp.ub)
    gam.mat[i, 1] <- get(paste0("gam", j))
    stop.temp0 <- get(paste0("d0", j)) >= logit(get(paste0("gam", j)))
    res.mat0[i, 1] <- mean(stop.temp0)
    
    ## get stopping probability under H1 for first analysis
    assign(paste0("d", j), ints[,j] + samps[i]*ratios[j]*slopes[,j])
    stop.temp <- get(paste0("d", j)) >= logit(get(paste0("gam", j)))
    res.mat[i, 1] <- mean(stop.temp)
    ## repeat process for the other analyses
    for (j in 2:ncol(m1)){
      assign(paste0("d0", j), ints0[,j] + samps[i]*ratios[j]*slopes0[,j])
      
      gam.temp.lb <- 0
      gam.temp.ub <- 1
      while(gam.temp.ub - gam.temp.lb > 0.00010001){
        gam.temp <- round(0.5*(gam.temp.lb + gam.temp.ub), 4)
        stop.bs <- mean(ifelse(stop.temp0, 1, 
                               ifelse(get(paste0("d0", j)) >= logit(gam.temp), 1, 0)))
        if (stop.bs > alp[j]){
          gam.temp.lb <- gam.temp
        } else {
          gam.temp.ub <- gam.temp
        }
      }
      
      assign(paste0("gam", j), gam.temp.ub)
      gam.mat[i, j] <- get(paste0("gam", j))
      stop.temp0 <- ifelse(stop.temp0, 1, 
                           ifelse(get(paste0("d0", j)) >= logit(get(paste0("gam", j))), 1, 0))
      res.mat0[i, j] <- mean(stop.temp0)
      
      assign(paste0("d", j), ints[,j] + samps[i]*ratios[j]*slopes[,j])
      
      ## ensure that previously stopped repetitions are still stopped
      stop.temp <- ifelse(stop.temp, 1, 
                          ifelse(get(paste0("d", j)) >= logit(get(paste0("gam", j))), 1, 0))
      
      res.mat[i, j] <- mean(stop.temp)
    }
  }
  
  ## output matrix where rows correspond to sample sizes and columns
  ## correspond to the stopping probability for a given analysis
  lst <- list(res.mat, res.mat0, gam.mat)
  
  return(lst)
  
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

## this function uses linear approximations to create a matrix of stopping probabilities across 
## a range of sample sizes that can be plotted or used to find optimal sample sizes
stopMatBoth <- function(lines, c_vec, lb, ub, by, gam, xi){
  
  ## the inputs are described as follows:
  ## lines: a matrix of intercepts and slopes returned by getLines()
  ## c_vec: the vector dictating the spacing of the analyses
  ## lb: the smallest sample size (n1) for which the stopping probs are estimated
  ## ub: the largest sample size (n1) for which the stopping probs are estimated
  ## by: the increments between lb and ub at which the stopping probs are estimated
  ## gam: the vector of sucess thresholds for all analyses
  ## xi: the vector of failure thresholds for the first T-1 analyses
  
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
  return(res.mat)
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

## matrix for stopping from data for design 5 (only posterior probabilities)
stop_dat <- function(m1, gam){
  # m1 is matrix of sequential probs (first sampling distribution estimate)
  # m2 is the matrix of sequential probs (second sampling distribution estimate)
  # n1s are the sample sizes for the first sampling distribution estimate
  # n2s are the sample sizes for the second sampling distribution estimate
  # lb is lower bound for the first analysis
  # ub is upper bound for the first analysis
  # by is the increment for the first analysis
  # gam are the thresholds
  
  ## need negation of logit since we work with complementary probabilities
  m1 <- 1 - m1
  
  res.vec <- NULL
  
  j <- 1
  assign(paste0("d", j), m1[,j])
  stop.temp <- get(paste0("d", j)) >= gam[j]
  res.vec[j] <- mean(stop.temp)
  ## repeat process for the other analyses
  for (j in 2:ncol(m1)){
    assign(paste0("d", j), m1[,j])
    
    ## ensure that previously stopped repetitions are still stopped
    stop.temp <- ifelse(stop.temp, 1, 
                        ifelse(get(paste0("d", j)) >= gam[j], 1, 0))
    
    res.vec[j] <- mean(stop.temp)
  }
  
  return(res.vec)
  
}