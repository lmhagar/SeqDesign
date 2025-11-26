## this code file contains the code to reproduce the numerical results and
## and plots in Section 6.1 of the text. Please run "01_ex_functions.R" first
## to ensure that the necessary packages and functions are loaded

## load the more packages
require(cowplot)
require(ggpubr)
require(mvtnorm)
library(gridExtra)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get matrix of beta coefficients were each row corresponds to a scenario
## scenario 1 is the only one used in Section 6.1 of the paper
betas.mat <- rbind(c(log(1.25), 0, log(1.4)-log(1.25)),
                   c(log(1.3), 0, log(1.6)-log(1.3)), 
                   c(log(1.35), 0, log(1.8)-log(1.35)))
## get probabilities for censoring vector
censors.v1 <- c(0.0250, 0.0378, 0.0768, 0.1097)

## estimate sampling distributions for equally spaced setting
## we estimate many sampling distributions here for plotting
ns <- seq(50, 250, 25)
j <- 3
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- getCompProbs(n = ns[i], c_vec = c_vec = c(1, 2, 3, 4, 5), cutpoints = c(7, 14, 21, 28),
                             hazards = c(0.055, 0.095, 0.04, 0.0200),
                             beta = as.numeric(betas.mat[j,]),
                             probs = c(0.5, 1/3), censors = censors.v1,
                             deltaL = 1, R = 10000, jags_model = "jags_surv.txt", 
                             prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                             burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000, 
                             seed = 100 + 10000)
  write.csv(probs.temp, paste0("comp_probs_T5_equal_scen", j,"_", ns[i], ".csv"), row.names = FALSE)
}

## now we get the sampling distribution estimates at na = 100 under H0 (to tune decision thresholds)
ns <- 100
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- getCompProbs(n = ns[i], c_vec = c_vec = c(1, 2, 3, 4, 5), cutpoints = c(7, 14, 21, 28),
                             hazards = c(0.055, 0.095, 0.04, 0.0200),
                             beta = rep(0,3),
                             probs = c(0.5, 1/3), censors = censors.v1,
                             deltaL = 1, R = 10000, jags_model = "jags_surv.txt", 
                             prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                             burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000, 
                             seed = 100 + 30000)
  write.csv(probs.temp, paste0("comp_probs_H0_T5_equal_scen", j,"_", ns[i], ".csv"), row.names = FALSE)
}

## retain the original vector of sample sizes for plotting
ns <- seq(50, 250, 25)

## now we fine tune the decision thresholds where failure thresholds
## are fixed
tuneSuccessDTs <- function(samp, fails, t1E){
  gam.low <- 0.5
  gam.up <- 1
  while (gam.up - gam.low > 0.0001){
    gam.next <- 0.5*(gam.low + gam.up)
    treshs <- rep(gam.next, ncol(samp))
    samp.success <- NULL
    samp.fail <- NULL
    for (i in 1:ncol(samp)){
      samp.success <- cbind(samp.success, samp[,i] >= treshs[i])
      if (i < ncol(samp)){
        samp.fail <- cbind(samp.fail, samp[,i] < fails[i])
      }
    }
    
    samp.bool <- samp.success[,1]
    for (i in 2:ncol(samp)){
      samp.temp <- ifelse(samp.fail[,i-1],
                          0, samp.success[,i])
      samp.bool <- cbind(samp.bool, samp.temp)
    }
    
    t1E.temp <- mean(as.numeric(rowSums(samp.bool)) > 0)
    
    ## increase gamma.next until type I error rate is small enough
    if (t1E.temp > t1E){
      gam.low <- gam.next
    } else {
      gam.up <- gam.next
    }
  }
  return(gam.up)
}

## load in sampling distribution estimate from Algorithm 2
sampH0.eq <- read.csv("comp_probs_H0_T5_equal_scen3_100.csv")
sampH0.eq <- 1 - sampH0.eq 

## get decision thresholds for each design
for (k in 2:5){
  assign(paste0("gamma.eq", k),
         rep(tuneSuccessDTs(sampH0.eq[, 1:k], rep(0.5, k-1), 0.025), k))
  
}

## this function is used to get our linear approximations given the sampling
## distribution estimates contained in the .csv files
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
    if (ncol(ints) > 3){
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

## create linear approximations for the appendix figure

## load in the sampling distribution estimates
probs.eq100 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 100, ".csv"))
probs.eq100 <- probs.eq100[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
probs.eq200 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 200, ".csv"))
probs.eq200 <- probs.eq200[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]

## get linear approximations for each design
for (k in 2:5){
  assign(paste0("lin.prep.eq", k),
         getLines(m1 = probs.eq100[, 1:(2*k - 1)], m2 = probs.eq200[, 1:(2*k - 1)], 100*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                  200*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)]))
  
}

## use linear approximations to obtain matrices of stopping probabilities
for (k in 2:5){
  assign(paste0("lin.eq", k),
         stopMat(lines = get(paste0("lin.prep.eq", k)), c_vec = c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                 50, 250, by = 1, 
                 gam = get(paste0("gamma.eq", k)), xi = rep(0.5, k-1)))
  
}

## use this function to quickly get confirmation probabilities for all analyses
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
  if (tps > 2){
    for (j in 2:(tps - 1)){
      stop.p.temp <- joint[,2*j-1] >= gam[j]
      stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
      res[2*j-1] <- mean(stop.p)
      ## repeat process for failure thresholds
      stop.f.temp <- joint[,2*j] < xi[j]
      stop.f <- ifelse(stop.p, 0, ifelse(stop.f, 1, stop.f.temp))
      res[2*j] <- mean(stop.f)
    }
  }
  
  ## implement process for final analysis (stopping only for success)
  j <- tps
  stop.p.temp <- joint[,2*j-1] >= gam[j]
  stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
  res[2*j-1] <- mean(stop.p)
  
  return(res)
}

## for all designs

ns <- seq(50, 250, 25)
## repeat with the equally spaced analyses
for (k in 2:5){
  assign(paste0("stop.eq",k),  NULL)
  for (i in 1:length(ns)){
    print(ns[i])
    data.temp <- read.csv(paste0("comp_probs_T5_equal_scen3_", ns[i], ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
    stop.temp <- stopData(data.temp[, 1:(2*k-1)], get(paste0("gamma.eq",k)), rep(0.5, k-1))
    
    assign(paste0("stop.eq",k), rbind(get(paste0("stop.eq",k)), stop.temp))
  }
}

## now get the bootstrap confidence intervals for each T value

## start with T = 2
k <- 2
## read in initial samples
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 100, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp1 <- samp1[, 1:(2*k - 1)]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 200, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp2 <- samp2[, 1:(2*k - 1)]

registerDoSNOW(cl)
pb <- txtProgressBar(max = MM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

logit <- function(x){log(x) - log(1-x)}

## get stopping probabilities based on bootstrap sampling distribution estimates
boot.eq2 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      lines.temp <- getLines(m1 = samp1.temp, m2 = samp2.temp, 100*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                                             200*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)])
                      
                      stopMat(lines = lines.temp, c_vec = c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                              50, 250, by = 1, 
                              gam = gamma.eq2, xi = rep(0.5, k-1), boot = TRUE)
                    }

## get quantiles of stopping probabilities for the percentile bootstrap method
lb.eq2 <- as.numeric(apply(boot.eq2, 2, quantile, probs = 0.025))
ub.eq2 <- as.numeric(apply(boot.eq2, 2, quantile, probs = 0.975))

## save results
write.csv(lb.eq2, "lb_eq2_appd.csv", row.names = FALSE)
write.csv(ub.eq2, "ub_eq2_appd.csv", row.names = FALSE)

## repeat with T = 3
k <- 3
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 100, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp1 <- samp1[, 1:(2*k - 1)]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 200, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp2 <- samp2[, 1:(2*k - 1)]

boot.eq3 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      lines.temp <- getLines(m1 = samp1.temp, m2 = samp2.temp, 100*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                                             200*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)])
                      
                      stopMat(lines = lines.temp, c_vec = c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                              50, 250, by = 1, 
                              gam = gamma.eq3, xi = rep(0.5, k-1), boot = TRUE)
                    }

lb.eq3 <- as.numeric(apply(boot.eq3, 2, quantile, probs = 0.025))
ub.eq3 <- as.numeric(apply(boot.eq3, 2, quantile, probs = 0.975))

write.csv(lb.eq3, "lb_eq3_appd.csv", row.names = FALSE)
write.csv(ub.eq3, "ub_eq3_appd.csv", row.names = FALSE)

## repeat with T = 4
k <- 4
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 100, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp1 <- samp1[, 1:(2*k - 1)]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 200, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp2 <- samp2[, 1:(2*k - 1)]

boot.eq4 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      lines.temp <- getLines(m1 = samp1.temp, m2 = samp2.temp, 100*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                                             200*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)])
                      
                      stopMat(lines = lines.temp, c_vec = c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                              50, 250, by = 1, 
                              gam = gamma.eq4, xi = rep(0.5, k-1), boot = TRUE)
                    }

lb.eq4 <- as.numeric(apply(boot.eq4, 2, quantile, probs = 0.025))
ub.eq4 <- as.numeric(apply(boot.eq4, 2, quantile, probs = 0.975))

write.csv(lb.eq4, "lb_eq4_appd.csv", row.names = FALSE)
write.csv(ub.eq4, "ub_eq4_appd.csv", row.names = FALSE)

## repeat with T = 5
k <- 5
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 100, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp1 <- samp1[, 1:(2*k - 1)]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen3_", 200, ".csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)]
samp2 <- samp2[, 1:(2*k - 1)]

boot.eq5 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      lines.temp <- getLines(m1 = samp1.temp, m2 = samp2.temp, 100*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                                             200*c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)])
                      
                      stopMat(lines = lines.temp, c_vec = c(1, 1, 2, 2, 3, 3, 4, 4, 5)[1:(2*k - 1)], 
                              50, 250, by = 1, 
                              gam = gamma.eq5, xi = rep(0.5, k-1), boot = TRUE)
                    }

lb.eq5 <- as.numeric(apply(boot.eq5, 2, quantile, probs = 0.025))
ub.eq5 <- as.numeric(apply(boot.eq5, 2, quantile, probs = 0.975))

write.csv(lb.eq5, "lb_eq5_appd.csv", row.names = FALSE)
write.csv(ub.eq5, "ub_eq5_appd.csv", row.names = FALSE)

## get colour palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## create the figure in Appendix D.2
## create data frames for the first analysis
## start with the components of the data frames based on our linear approximations
n.cols <- 201
j <- 1

for (k in 2:5){
  
  assign(paste0("df1.t", k), data.frame(n = c(seq(50, 250, 1), ns, rep(seq(50, 250, 1), 2)),
                                        prob = c(get(paste0("lin.eq", k))[,1], get(paste0("stop.eq", k))[,1], 
                                                 read.csv(paste0("lb_eq", k, "_appd.csv"))$x[1:n.cols], read.csv(paste0("ub_eq", k, "_appd.csv"))$x[1:n.cols]),
                                        type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
}

## get first subplot for each design
for (k in 2:4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(''), y= bquote('Stopping Probability')) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1)) 
}

## axis labels are different for the final row
for (k in 5){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(italic(n)[1]), y= bquote('Stopping Probability')) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1)) 
}

## repeat with the subplots for the second analysis
j <- 2
for (k in 2:5){
  
  assign(paste0("df", j, ".t", k), data.frame(n = 2*c(seq(50, 250, 1), ns, rep(seq(50, 250, 1), 2)),
                                              prob = c(get(paste0("lin.eq", k))[,2*j-1], get(paste0("stop.eq", k))[,2*j-1], 
                                                       read.csv(paste0("lb_eq", k, "_appd.csv"))$x[(j*n.cols + 1):((j+1)*n.cols)], read.csv(paste0("ub_eq", k, "_appd.csv"))$x[(j*n.cols + 1):((j+1)*n.cols)]),
                                              type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
}

for (k in 2:4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(''), y= bquote('')) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1)) 
}

for (k in 5){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(italic(n)[2]), y= bquote('')) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1)) 
}

## repeat with the subplots for the third analysis
j <- 3
for (k in 3:5){
  
  assign(paste0("df", j, ".t", k), data.frame(n = 3*c(seq(50, 250, 1), ns, rep(seq(50, 250, 1), 2)),
                                              prob = c(get(paste0("lin.eq", k))[,2*j-1], get(paste0("stop.eq", k))[,2*j-1], 
                                                       read.csv(paste0("lb_eq", k, "_appd.csv"))$x[((j+1)*n.cols + 1):((j+2)*n.cols)], read.csv(paste0("ub_eq", k, "_appd.csv"))$x[((j+1)*n.cols + 1):((j+2)*n.cols)]),
                                              type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
}

## these subplots are only relevant for the designs with T >= 3
for (k in 3:4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(''), y= bquote('')) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1))  
}

for (k in 5){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(italic(n)[3]), y= bquote('')) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1))  
}

## repeat process to get subplots for the fourth analysis
j <- 4
for (k in 4:5){
  assign(paste0("df", j, ".t", k), data.frame(n = 4*c(seq(50, 250, 1), ns, rep(seq(50, 250, 1), 2)),
                                              prob = c(get(paste0("lin.eq", k))[,2*j-1], get(paste0("stop.eq", k))[,2*j-1], 
                                                       read.csv(paste0("lb_eq", k, "_appd.csv"))$x[((j+2)*n.cols + 1):((j+3)*n.cols)], read.csv(paste0("ub_eq", k, "_appd.csv"))$x[((j+2)*n.cols + 1):((j+3)*n.cols)]),
                                              type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
}

for (k in 4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(''), y= bquote('')) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1))
}

for (k in 5){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(italic(n)[4]), y= bquote('')) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1))
}

## repeat process to get subplots for the fifth analysis
j <- 5
for (k in 5){
  
  assign(paste0("df", j, ".t", k), data.frame(n = 5*c(seq(50, 250, 1), ns, rep(seq(50, 250, 1), 2)),
                                              prob = c(get(paste0("lin.eq", k))[,2*j-1], get(paste0("stop.eq", k))[,2*j-1], 
                                                       read.csv(paste0("lb_eq", k, "_appd.csv"))$x[((j+3)*n.cols + 1):((j+4)*n.cols)], read.csv(paste0("ub_eq", k, "_appd.csv"))$x[((j+3)*n.cols + 1):((j+4)*n.cols)]),
                                              type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
}

for (k in 5){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 1)) +
           labs(color  = "", linetype = "") +
           labs(x= bquote(italic(n)[5]), y= bquote('')) +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
           theme(plot.title = element_text(size=20,face="bold",
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18)) +
           theme(legend.text=element_text(size=18)) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
           theme(legend.position="none") +
           ylim(0,1)) 
}

## get matrix for plotting grid

m <- matrix(NA, 4, 5)
m[1, 1:2] <- 1:2
m[2, 1:3] <- 3:5
m[3, 1:4] <- 6:9
m[4, 1:5] <- 10:14

## get simplified legend for plotting
df5.t5.leg <- subset(df5.t5, !(df5.t5$type %in% c(4, 8)))

j <- 5; k <- 5
assign(paste0("plot", j, ".t", k, ".leg"), 
       ggplot(get(paste0("df",j,".t", k, ".leg")), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Equal)  ", "Simulated (Equal)  ", 
                                                  "95% Bootstrap CI (Equal)  "),
                            values = cbb[c(4, 1, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Equal)  ", "Simulated (Equal)  ", 
                                                     "95% Bootstrap CI (Equal)  "),
                               values = rep(c("solid", "longdash", "dashed"), 1)) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(n)[5]), y= bquote('')) +
         theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=14),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="bottom",
               legend.key.size = unit(0.925, "cm")) +
         ylim(0,1)) 


fig.row1 <- grid.arrange(grobs = list(plot1.t2 + theme(plot.margin=unit(c(0,0.25,-0.35,0.1),"cm")), 
                                      plot2.t2 + theme(plot.margin=unit(c(0,0.25,-0.35,-0.65),"cm")),
                                      plot1.t3 + theme(plot.margin=unit(c(0,0.25,-0.35,0.1),"cm")), 
                                      plot2.t3 + theme(plot.margin=unit(c(0,0.25,-0.35,-0.65),"cm")), 
                                      plot3.t3 + theme(plot.margin=unit(c(0,0.25,-0.35,-0.65),"cm")),
                                      plot1.t4 + theme(plot.margin=unit(c(0,0.25,-0.35,0.1),"cm")), 
                                      plot2.t4 + theme(plot.margin=unit(c(0,0.25,-0.35,-0.65),"cm")), 
                                      plot3.t4 + theme(plot.margin=unit(c(0,0.25,-0.35,-0.65),"cm")),
                                      plot4.t4 + theme(plot.margin=unit(c(0,0.25,-0.35,-0.65),"cm")),
                                      plot1.t5 + theme(plot.margin=unit(c(0,0.25,0.1,0.1),"cm")), 
                                      plot2.t5 + theme(plot.margin=unit(c(0,0.25,0.1,-0.65),"cm")), 
                                      plot3.t5 + theme(plot.margin=unit(c(0,0.25,0.1,-0.65),"cm")), 
                                      plot4.t5 + theme(plot.margin=unit(c(0,0.25,0.1,-0.65),"cm")), 
                                      plot5.t5 + theme(plot.margin=unit(c(0,0.25,0.1,-0.65),"cm"))), 
                         layout_matrix = m, widths = c(1.2475, rep(1,4)), heights = c(rep(1, 3), 1.1475))

fig_final <- plot_grid(fig.row1, get_legend(get(paste0("plot5.t5.leg"))), ncol = 1, rel_heights = c(1, .05))

# output as .pdf file for the article
pdf(file = paste0("Fig_AppD2.pdf"),   # The directory you want to save the file in
    width = 1.25*12, # The width of the plot in inches
    height = 1.25*10) # The height of the plot in inches

fig_final

dev.off()

## get the sample sizes for the confirmatory sampling distribution estimates

## T = 2
mat.eq2 <- cbind(seq(50, 250, 1), lin.eq2)
ind.eq2 <- min(which(mat.eq2[,4] >= 0.8))
n1.eq2 <- mat.eq2[ind.eq2,1]

## Expected total sample size
n1.eq2*(mat.eq2[ind.eq2, 2] + mat.eq2[ind.eq2, 3]) + 2*n1.eq2*(1 - mat.eq2[ind.eq2, 2] - mat.eq2[ind.eq2, 3])

## T = 3
mat.eq3 <- cbind(seq(50, 250, 1), lin.eq3)
ind.eq3 <- min(which(mat.eq3[,6] >= 0.8))
n1.eq3 <- mat.eq3[ind.eq3,1]

## Expected total sample size
n1.eq3*(mat.eq3[ind.eq3, 2] + mat.eq3[ind.eq3, 3]) + 2*n1.eq3*(mat.eq3[ind.eq3, 4] - mat.eq3[ind.eq3, 2] + 
                                                                 mat.eq3[ind.eq3, 5] - mat.eq3[ind.eq3, 3]) +
  3*n1.eq3*(1- mat.eq3[ind.eq3, 5] - mat.eq3[ind.eq3, 4])

## T = 4
mat.eq4 <- cbind(seq(50, 250, 1), lin.eq4)
ind.eq4 <- min(which(mat.eq4[,8] >= 0.8))
n1.eq4 <- mat.eq4[ind.eq4,1]

## Expected total sample size
n1.eq4*(mat.eq4[ind.eq4, 2] + mat.eq4[ind.eq4, 3]) + 2*n1.eq4*(mat.eq4[ind.eq4, 4] - mat.eq4[ind.eq4, 2] +
                                                                 mat.eq4[ind.eq4, 5] - mat.eq4[ind.eq4, 3]) + 
  3*n1.eq4*(mat.eq4[ind.eq4, 6] - mat.eq4[ind.eq4, 4] +
              mat.eq4[ind.eq4, 7] - mat.eq4[ind.eq4, 5]) +
  4*n1.eq4*(1- mat.eq4[ind.eq4, 7] - mat.eq4[ind.eq4, 6])

## T = 5
mat.eq5 <- cbind(seq(50, 250, 1), lin.eq5)
ind.eq5 <- min(which(mat.eq5[,10] >= 0.8))
n1.eq5 <- mat.eq5[ind.eq5,1]

## Expected total sample size
n1.eq5*(mat.eq5[ind.eq5, 2] + mat.eq5[ind.eq5, 3]) + 2*n1.eq5*(mat.eq5[ind.eq5, 4] - mat.eq5[ind.eq5, 2] +
                                                                 mat.eq5[ind.eq5, 5] - mat.eq5[ind.eq5, 3]) + 
  3*n1.eq5*(mat.eq5[ind.eq5, 6] - mat.eq5[ind.eq5, 4] +
              mat.eq5[ind.eq5, 7] - mat.eq5[ind.eq5, 5]) +
  4*n1.eq5*(mat.eq5[ind.eq5, 8] - mat.eq5[ind.eq5, 6] +
              mat.eq5[ind.eq5, 9] - mat.eq5[ind.eq5, 7]) +
  5*n1.eq5*(1- mat.eq5[ind.eq5, 9] - mat.eq5[ind.eq5, 8])

## now get the confirmatory power and type I error rate estimate for the table
## in the main text

## implement for T = 2

ns <- 143
j <- 3
## confirm power
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 60)*100)
  write.csv(probs.temp, paste0("confirm_probs_T2_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

## confirm type I error rate
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 80)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T2_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

## implement for T = 3

ns <- 103
j <- 3
## confirm power
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2, 3), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 60)*100)
  write.csv(probs.temp, paste0("confirm_probs_T3_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2, 3), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 80)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T3_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

## implement for T = 4
ns <- 80
j <- 3
## confirm power
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2, 3, 4), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 60)*100)
  write.csv(probs.temp, paste0("confirm_probs_T4_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2, 3, 4), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 80)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T4_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

## implement for T = 5

j <- 3
ns <- 70
## confirm power
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2, 3, 4, 5), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 60)*100)
  write.csv(probs.temp, paste0("confirm_probs_T5_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 2, 3, 4, 5), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 80)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T5_equal_scen3_", ns[i], ".csv"), row.names = FALSE)
}

## compute the confirmatory estimates of the stopping probabilities

## implement for T = 2

## type I error rate
stopData(read.csv(paste0("confirm_probs_H0_T2_equal_scen3_143.csv"))[, c(1, 1, 2)], 
         gamma.eq2, rep(0.5, 1))

## power
stopData(read.csv(paste0("confirm_probs_T2_equal_scen3_143.csv"))[, c(1, 1, 2)], 
         gamma.eq2, rep(0.5, 1))

## implement for T = 3
stopData(read.csv(paste0("confirm_probs_H0_T3_equal_scen3_103.csv"))[, c(1, 1, 2, 2, 3)], 
         gamma.eq3, rep(0.5, 2))

stopData(read.csv(paste0("confirm_probs_T3_equal_scen3_103.csv"))[, c(1, 1, 2, 2, 3)], 
         gamma.eq3, rep(0.5, 2))

## implement for T = 4
stopData(read.csv(paste0("confirm_probs_H0_T4_equal_scen3_80.csv"))[, c(1, 1, 2, 2, 3, 3, 4)], 
         gamma.eq4, rep(0.5, 3))

stopData(read.csv(paste0("confirm_probs_T4_equal_scen3_80.csv"))[, c(1, 1, 2, 2, 3, 3, 4)], 
         gamma.eq4, rep(0.5, 3))

## implement for T = 5
stopData(read.csv(paste0("confirm_probs_H0_T5_equal_scen3_70.csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)], 
         gamma.eq5, rep(0.5, 4))

stopData(read.csv(paste0("confirm_probs_T5_equal_scen3_70.csv"))[, c(1, 1, 2, 2, 3, 3, 4, 4, 5)], 
         gamma.eq5, rep(0.5, 4))