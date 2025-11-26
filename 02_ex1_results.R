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

## estimate sampling distributions for backloaded setting
## we estimate many sampling distributions here for plotting
j <- 1
ns <- seq(100, 500, 50)
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- getCompProbs(n = ns[i], c_vec = c_vec = c(1, 1.5, 2, 2.5, 3), cutpoints = c(7, 14, 21, 28),
                             hazards = c(0.055, 0.095, 0.04, 0.0200),
                             beta = as.numeric(betas.mat[j,]),
                             probs = c(0.5, 1/3), censors = censors.v1,
                             deltaL = 1, R = 10000, jags_model = "jags_surv.txt", 
                             prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                             burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000, 
                             seed = i*100)
  write.csv(probs.temp, paste0("comp_probs_T5_back_scen", j,"_", ns[i], ".csv"), row.names = FALSE)
}

## repeat estimation process for equally spaced analyses
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

## now we get the sampling distribution estimates at na = 200 under H0 (to tune decision thresholds)
ns <- 200
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- getCompProbs(n = ns[i], c_vec = c_vec = c(1, 1.5, 2, 2.5, 3), cutpoints = c(7, 14, 21, 28),
                             hazards = c(0.055, 0.095, 0.04, 0.0200),
                             beta = rep(0,3),
                             probs = c(0.5, 1/3), censors = censors.v1,
                             deltaL = 1, R = 10000, jags_model = "jags_surv.txt", 
                             prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                             burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000, 
                             seed = i*100 + 20000)
  write.csv(probs.temp, paste0("comp_probs_H0_T5_back_scen", j,"_", ns[i], ".csv"), row.names = FALSE)
}

## repeat estimation process for equally spaced analyses
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
ns <- seq(100, 500, 50)

## now we fine tune the decision thresholds
## function to get get the decision threshold for second analysis
getGamma2 <- function(gamma.prev, alpha.cum, info.times){
  cor.mat <- outer(info.times, info.times, FUN = function(x, y) sqrt(pmin(x, y) / pmax(x, y)))
  error_fn <- function(gamma.prev, cor.mat, gamma.new, target.err){
    pmvnorm(lower = c(-Inf, qnorm(gamma.new)), upper = c(qnorm(gamma.prev[1]), Inf), mean = rep(0, 2), corr = cor.mat) - target.err
  }
  return(uniroot(error_fn, lower = 0.5, upper = 0.999999, gamma.prev = gamma.prev, cor.mat = cor.mat, target.err = alpha.cum[2] - alpha.cum[1]))
}

## function to get the decision threshold for third analysis
getGamma3 <- function(gamma.prev, alpha.cum, info.times){
  cor.mat <- outer(info.times, info.times, FUN = function(x, y) sqrt(pmin(x, y) / pmax(x, y)))
  error_fn <- function(gamma.prev, cor.mat, gamma.new, target.err){
    pmvnorm(lower = c(-Inf, -Inf, qnorm(gamma.new)), upper = c(qnorm(gamma.prev[1]), qnorm(gamma.prev[2]), Inf), mean = rep(0, 3), corr = cor.mat) - target.err
  }
  return(uniroot(error_fn, lower = 0.5, upper = 0.999999, gamma.prev = gamma.prev, cor.mat = cor.mat, target.err = alpha.cum[3] - alpha.cum[2]))
}

## function to get the decision threshold for fourth analysis
getGamma4 <- function(gamma.prev, alpha.cum, info.times){
  cor.mat <- outer(info.times, info.times, FUN = function(x, y) sqrt(pmin(x, y) / pmax(x, y)))
  error_fn <- function(gamma.prev, cor.mat, gamma.new, target.err){
    pmvnorm(lower = c(-Inf, -Inf, -Inf, qnorm(gamma.new)), 
            upper = c(qnorm(gamma.prev[1]), qnorm(gamma.prev[2]), qnorm(gamma.prev[3]), Inf), 
            mean = rep(0, 4), corr = cor.mat) - target.err
  }
  return(uniroot(error_fn, lower = 0.5, upper = 0.999999, gamma.prev = gamma.prev, cor.mat = cor.mat, target.err = alpha.cum[4] - alpha.cum[3]))
}

## function to get the decision threshold for fifth analysis
getGamma5 <- function(gamma.prev, alpha.cum, info.times){
  cor.mat <- outer(info.times, info.times, FUN = function(x, y) sqrt(pmin(x, y) / pmax(x, y)))
  error_fn <- function(gamma.prev, cor.mat, gamma.new, target.err){
    pmvnorm(lower = c(-Inf, -Inf, -Inf, -Inf, qnorm(gamma.new)), 
            upper = c(qnorm(gamma.prev[1]), qnorm(gamma.prev[2]), qnorm(gamma.prev[3]), qnorm(gamma.prev[4]), Inf), 
            mean = rep(0, 5), corr = cor.mat) - target.err
  }
  return(uniroot(error_fn, lower = 0.5, upper = 0.999999, gamma.prev = gamma.prev, cor.mat = cor.mat, target.err = alpha.cum[5] - alpha.cum[4]))
}

## Define information fractions
info.back2 <- c(1, 1.5)/1.5
info.back3 <- c(1, 1.5, 2)/2
info.back4 <- c(1, 1.5, 2, 2.5)/2.5
info.back5 <- c(1, 1.5, 2, 2.5, 3)/3

info.eq2 <- c(1, 2)/2
info.eq3 <- c(1, 2, 3)/3
info.eq4 <- c(1, 2, 3, 4)/4
info.eq5 <- c(1, 2, 3, 4, 5)/5

## get alpha vectors using alpha-spending function for each setting
a.back2 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.back2)))
a.back3 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.back3)))
a.back4 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.back4)))
a.back5 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.back5)))

a.eq2 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.eq2)))
a.eq3 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.eq3)))
a.eq4 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.eq4)))
a.eq5 <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(info.eq5)))

## get initial decision threshold values for the backloaded case
for (k in 2:5){
  assign(paste0("gamma.back", k), NULL)
  assign(paste0("gamma.back", k),
         c(get(paste0("gamma.back", k)), 1 - get(paste0("a.back", k))[1]))
  assign(paste0("gamma.back", k),
         c(get(paste0("gamma.back", k)), 
           getGamma2(1 - get(paste0("a.back", k))[1], get(paste0("a.back", k))[1:2], get(paste0("info.back", k))[1:2])$root))
  if (k > 2){
    assign(paste0("gamma.back", k),
           c(get(paste0("gamma.back", k)), 
             getGamma3(get(paste0("gamma.back", k)),
                       get(paste0("a.back", k))[1:3], get(paste0("info.back", k))[1:3])$root))
  }
  if (k > 3){
    assign(paste0("gamma.back", k),
           c(get(paste0("gamma.back", k)), 
             getGamma4(get(paste0("gamma.back", k)),
                       get(paste0("a.back", k))[1:4], get(paste0("info.back", k))[1:4])$root))
  }
  if (k > 4){
    assign(paste0("gamma.back", k),
           c(get(paste0("gamma.back", k)), 
             getGamma5(get(paste0("gamma.back", k)),
                       get(paste0("a.back", k))[1:5], get(paste0("info.back", k))[1:5])$root))
  }
}

## get initial decision threshold values for the equally spaced case
for (k in 2:5){
  assign(paste0("gamma.eq", k), NULL)
  assign(paste0("gamma.eq", k),
         c(get(paste0("gamma.eq", k)), 1 - get(paste0("a.eq", k))[1]))
  assign(paste0("gamma.eq", k),
         c(get(paste0("gamma.eq", k)), 
           getGamma2(1 - get(paste0("a.eq", k))[1], get(paste0("a.eq", k))[1:2], get(paste0("info.eq", k))[1:2])$root))
  if (k > 2){
    assign(paste0("gamma.eq", k),
           c(get(paste0("gamma.eq", k)), 
             getGamma3(get(paste0("gamma.eq", k)),
                       get(paste0("a.eq", k))[1:3], get(paste0("info.eq", k))[1:3])$root))
  }
  if (k > 3){
    assign(paste0("gamma.eq", k),
           c(get(paste0("gamma.eq", k)), 
             getGamma4(get(paste0("gamma.eq", k)),
                       get(paste0("a.eq", k))[1:4], get(paste0("info.eq", k))[1:4])$root))
  }
  if (k > 4){
    assign(paste0("gamma.eq", k),
           c(get(paste0("gamma.eq", k)), 
             getGamma5(get(paste0("gamma.eq", k)),
                       get(paste0("a.eq", k))[1:5], get(paste0("info.eq", k))[1:5])$root))
  }
}

## use the initial decision thresholds to obtain the tuned decision thresholds
tuneFinalDT <- function(samp, prev.gam, t1E){
  gam.low <- 0.5
  gam.up <- 1
  while (gam.up - gam.low > 0.0001){
    gam.next <- 0.5*(gam.low + gam.up)
    treshs <- c(prev.gam, gam.next)
    samp.bool <- NULL
    for (i in 1:ncol(samp)){
      samp.bool <- cbind(samp.bool, samp[,i] >= treshs[i])
    }
    t1E.temp <- mean(as.numeric(rowSums(samp.bool)) > 0)
    
    ## keep increasing final decision threshold until type I error
    ## rate is small enough
    if (t1E.temp > t1E){
      gam.low <- gam.next
    } else {
      gam.up <- gam.next
    }
  }
  return(gam.up)
}

## load in sampling distribution estimate from Algorithm 1 under H0
sampH0.back <- read.csv("comp_probs_H0_T5_back_scen1_200.csv")
sampH0.eq <- read.csv("comp_probs_H0_T5_equal_scen1_200.csv")

## complementary probabilities were computed before
sampH0.back <- 1 - sampH0.back
sampH0.eq <- 1 - sampH0.eq 

## get the final decision thresholds for the backloaded settings
for (k in 2:5){
  assign(paste0("gamma.back", k), head(get(paste0("gamma.back", k)), k-1))
  assign(paste0("gamma.back", k),
         c(get(paste0("gamma.back", k)), 
           tuneFinalDT(sampH0.back[, 1:k], head(get(paste0("gamma.back", k)), k-1), 0.025)))
  
}

## get the final decision thresholds for the equally spaced settings
for (k in 2:5){
  assign(paste0("gamma.eq", k), head(get(paste0("gamma.eq", k)), k-1))
  assign(paste0("gamma.eq", k),
         c(get(paste0("gamma.eq", k)), 
           tuneFinalDT(sampH0.eq[, 1:k], head(get(paste0("gamma.eq", k)), k-1), 0.025)))
  
}

## define criterion for power
target.pwr <- 0.8

## create linear approximations for Figure 1 (backloaded settings)
## load in sampling distribution estimates at na = 200 and nb = 400
probs.back200 <- read.csv(paste0("comp_probs_T5_back_scen1_", 200, ".csv"))
probs.back400 <- read.csv(paste0("comp_probs_T5_back_scen1_", 400, ".csv"))

## create matrix of cumulative stopping probabilities
for (k in 2:5){
  assign(paste0("lin.back", k),
         stop_mat(m1 = probs.back200[, 1:k], m2 = probs.back400[, 1:k], 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                  400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, get(paste0("gamma.back", k))))
  
}

## get bootstrap confidence intervals for the stopping probabilities
## define the number and size of the bootstrap samples
MM <- 2000
m <- 10000

## implement the process to construct bootstrap confidence intervals for T = 2 (backloaded)
## load in initial sampling distribution estimates
k <- 2
samp1 <- read.csv(paste0("comp_probs_T5_back_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_back_scen1_", 400, ".csv"))[, 1:k]

registerDoSNOW(cl)
pb <- txtProgressBar(max = MM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement bootstrap sampling
boot.back2 <- foreach(l=1:MM, .combine=rbind,
                      .options.snow=opts) %dopar% {
                        
                        ## resample from each of the initial sampling distribution estimates
                        samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                        samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                        
                        ## get stopping probabilities based on the resampled stopping probabilities
                        stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                                 400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.back2)
                      }

## get quantiles of the bootstrap stopping probabilities for each sample size
lb.back2 <- as.numeric(apply(boot.back2, 2, quantile, probs = 0.025))
ub.back2 <- as.numeric(apply(boot.back2, 2, quantile, probs = 0.975))

## save results to a .csv file
write.csv(lb.back2, "lb_back2.csv", row.names = FALSE)
write.csv(ub.back2, "ub_back2.csv", row.names = FALSE)

## repeat process for T = 3 with backloaded setting
k <- 3
samp1 <- read.csv(paste0("comp_probs_T5_back_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_back_scen1_", 400, ".csv"))[, 1:k]
boot.back3 <- foreach(l=1:MM, .combine=rbind,
                      .options.snow=opts) %dopar% {
                        
                        samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                        samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                        
                        stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                                 400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.back3)
                      }

lb.back3 <- as.numeric(apply(boot.back3, 2, quantile, probs = 0.025))
ub.back3 <- as.numeric(apply(boot.back3, 2, quantile, probs = 0.975))

write.csv(lb.back3, "lb_back3.csv", row.names = FALSE)
write.csv(ub.back3, "ub_back3.csv", row.names = FALSE)

## repeat process for T = 4 with backloaded setting
k <- 4
samp1 <- read.csv(paste0("comp_probs_T5_back_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_back_scen1_", 400, ".csv"))[, 1:k]
boot.back4 <- foreach(l=1:MM, .combine=rbind,
                      .options.snow=opts) %dopar% {
                        
                        samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                        samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                        
                        stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                                 400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.back4)
                      }

lb.back4 <- as.numeric(apply(boot.back4, 2, quantile, probs = 0.025))
ub.back4 <- as.numeric(apply(boot.back4, 2, quantile, probs = 0.975))

write.csv(lb.back4, "lb_back4.csv", row.names = FALSE)
write.csv(ub.back4, "ub_back4.csv", row.names = FALSE)

## repeat process for T = 5 with backloaded setting
k <- 5
samp1 <- read.csv(paste0("comp_probs_T5_back_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_back_scen1_", 400, ".csv"))[, 1:k]
boot.back5 <- foreach(l=1:MM, .combine=rbind,
                      .options.snow=opts) %dopar% {
                        
                        samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                        samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                        
                        stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                                 400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.back5)
                      }

lb.back5 <- as.numeric(apply(boot.back5, 2, quantile, probs = 0.025))
ub.back5 <- as.numeric(apply(boot.back5, 2, quantile, probs = 0.975))

write.csv(lb.back5, "lb_back5.csv", row.names = FALSE)
write.csv(ub.back5, "ub_back5.csv", row.names = FALSE)

## repeat with equally spaced settings
probs.eq200 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 200, ".csv"))
probs.eq400 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 400, ".csv"))

## create matrix of cumulative stopping probabilities
for (k in 2:5){
  assign(paste0("lin.eq", k),
         stop_mat(m1 = probs.eq200[, 1:k], m2 = probs.eq400[, 1:k], 200*c(1, 2, 3, 4, 5)[1:k], 
                  400*c(1, 2, 3, 4, 5)[1:k], 100, 500, by = 1, get(paste0("gamma.eq", k))))
  
}

## repeat process for T = 2 with equally spaced setting
k <- 2
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 400, ".csv"))[, 1:k]

boot.eq2 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                               400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.eq2)
                    }

lb.eq2 <- as.numeric(apply(boot.eq2, 2, quantile, probs = 0.025))
ub.eq2 <- as.numeric(apply(boot.eq2, 2, quantile, probs = 0.975))

write.csv(lb.eq2, "lb_eq2.csv", row.names = FALSE)
write.csv(ub.eq2, "ub_eq2.csv", row.names = FALSE)

## repeat process for T = 3 with equally spaced setting
k <- 3
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 400, ".csv"))[, 1:k]

boot.eq3 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                               400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.eq3)
                    }

lb.eq3 <- as.numeric(apply(boot.eq3, 2, quantile, probs = 0.025))
ub.eq3 <- as.numeric(apply(boot.eq3, 2, quantile, probs = 0.975))

write.csv(lb.eq3, "lb_eq3.csv", row.names = FALSE)
write.csv(ub.eq3, "ub_eq3.csv", row.names = FALSE)

## repeat process for T = 4 with equally spaced setting
k <- 4
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 400, ".csv"))[, 1:k]

boot.eq4 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                               400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.eq4)
                    }

lb.eq4 <- as.numeric(apply(boot.eq4, 2, quantile, probs = 0.025))
ub.eq4 <- as.numeric(apply(boot.eq4, 2, quantile, probs = 0.975))

write.csv(lb.eq4, "lb_eq4.csv", row.names = FALSE)
write.csv(ub.eq4, "ub_eq4.csv", row.names = FALSE)

## repeat process for T = 5 with equally spaced setting
k <- 5
samp1 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 200, ".csv"))[, 1:k]
samp2 <- read.csv(paste0("comp_probs_T5_equal_scen1_", 400, ".csv"))[, 1:k]

boot.eq5 <- foreach(l=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      stop_mat(m1 = samp1.temp, m2 = samp2.temp, 200*c(1, 1.5, 2, 2.5, 3)[1:k], 
                               400*c(1, 1.5, 2, 2.5, 3)[1:k], 100, 500, by = 1, gamma.eq5)
                    }

lb.eq5 <- as.numeric(apply(boot.eq5, 2, quantile, probs = 0.025))
ub.eq5 <- as.numeric(apply(boot.eq5, 2, quantile, probs = 0.975))

write.csv(lb.eq5, "lb_eq5.csv", row.names = FALSE)
write.csv(ub.eq5, "ub_eq5.csv", row.names = FALSE)

## create data frames for the first analysis (t = 1)
n.cols <- 401
j <- 1

## we create a data frame for each T value that will have a curve for the 
## estimated stopping probabilities, the simulated stopping probabilities, and the lower and upper
## bootstrap CIs

## we first implement this for the backloaded and equally spaced settings then combine the data frames
for (k in 2:5){
  assign(paste0("df1.back", k), data.frame(n = c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                           prob = c(get(paste0("lin.back", k))[1:n.cols], get(paste0("stop.back", k))[,1], 
                                                    read.csv(paste0("lb_back", k, ".csv"))$x[1:n.cols], read.csv(paste0("ub_back", k, ".csv"))$x[1:n.cols]),
                                           type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))
  
  assign(paste0("df1.eq", k), data.frame(n = c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                         prob = c(get(paste0("lin.eq", k))[1:n.cols], get(paste0("stop.eq", k))[,1], 
                                                  read.csv(paste0("lb_eq", k, ".csv"))$x[1:n.cols], read.csv(paste0("ub_eq", k, ".csv"))$x[1:n.cols]),
                                         type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
  assign(paste0("df1.t", k), rbind(get(paste0("df1.back", k)), get(paste0("df1.eq", k))))
}

## get the subpanels for T = 2 to 4
for (k in 2:4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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

## the legends are different for T = 5
for (k in 5){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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

## repeat this process to get the subplots for the second analysis (t = 2)
j <- 2
for (k in 2:5){
  assign(paste0("df", j, ".back", k), data.frame(n = 1.5*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                                 prob = c(get(paste0("lin.back", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.back", k))[,j], 
                                                          read.csv(paste0("lb_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                                 type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))
  
  assign(paste0("df", j, ".eq", k), data.frame(n = 2*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                               prob = c(get(paste0("lin.eq", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.eq", k))[,j], 
                                                        read.csv(paste0("lb_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                               type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
  assign(paste0("df", j, ".t", k), rbind(get(paste0("df", j, ".back", k)), get(paste0("df", j, ".eq", k))))
}

for (k in 2:4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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

## repeat this process to get the subplots for the third analysis (t = 3)
j <- 3
for (k in 3:5){
  assign(paste0("df", j, ".back", k), data.frame(n = 2*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                                 prob = c(get(paste0("lin.back", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.back", k))[,j], 
                                                          read.csv(paste0("lb_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                                 type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))
  
  assign(paste0("df", j, ".eq", k), data.frame(n = 3*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                               prob = c(get(paste0("lin.eq", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.eq", k))[,j], 
                                                        read.csv(paste0("lb_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                               type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
  assign(paste0("df", j, ".t", k), rbind(get(paste0("df", j, ".back", k)), get(paste0("df", j, ".eq", k))))
}

## these plots are only relevant for T >= 3
for (k in 3:4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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

## repeat this process to get the subplots for the fourth analysis (t = 4)
j <- 4
for (k in 4:5){
  assign(paste0("df", j, ".back", k), data.frame(n = 2.5*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                                 prob = c(get(paste0("lin.back", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.back", k))[,j], 
                                                          read.csv(paste0("lb_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                                 type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))
  
  assign(paste0("df", j, ".eq", k), data.frame(n = 4*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                               prob = c(get(paste0("lin.eq", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.eq", k))[,j], 
                                                        read.csv(paste0("lb_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                               type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
  assign(paste0("df", j, ".t", k), rbind(get(paste0("df", j, ".back", k)), get(paste0("df", j, ".eq", k))))
}

for (k in 4){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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

## repeat this process to get the subplots for the fifth analysis (t = 5)
j <- 5
for (k in 5){
  assign(paste0("df", j, ".back", k), data.frame(n = 3*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                                 prob = c(get(paste0("lin.back", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.back", k))[,j], 
                                                          read.csv(paste0("lb_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_back", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                                 type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))
  
  assign(paste0("df", j, ".eq", k), data.frame(n = 5*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                                               prob = c(get(paste0("lin.eq", k))[((j-1)*n.cols + 1):(j*n.cols)], get(paste0("stop.eq", k))[,j], 
                                                        read.csv(paste0("lb_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)], read.csv(paste0("ub_eq", k, ".csv"))$x[((j-1)*n.cols + 1):(j*n.cols)]),
                                               type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))
  
  assign(paste0("df", j, ".t", k), rbind(get(paste0("df", j, ".back", k)), get(paste0("df", j, ".eq", k))))
}

for (k in 5){
  assign(paste0("plot", j, ".t", k), 
         ggplot(get(paste0("df",j,".t", k)), 
                aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
           geom_line() +
           scale_color_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                    "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                    "Estimated (Equal)", "Simulated (Equal)", 
                                                    "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                              values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
           scale_linetype_manual(name = "", labels = c("Estimated (Backload)", "Simulated (Backload)", 
                                                       "95% Bootstrap CI (Backload)", "95% Bootstrap CI (Backload)",
                                                       "Estimated (Equal)", "Simulated (Equal)", 
                                                       "95% Bootstrap CI (Equal)", "95% Bootstrap CI (Equal)"),
                                 values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
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

## create matrix to plot grid of subplots
m <- matrix(NA, 4, 5)
m[1, 1:2] <- 1:2
m[2, 1:3] <- 3:5
m[3, 1:4] <- 6:9
m[4, 1:5] <- 10:14

## get simplified legend for plotting
df5.t5.leg <- subset(df5.t5, !(df5.t5$type %in% c(4, 8)))

## obtain plot with simplified legend
j <- 5; k <- 5
assign(paste0("plot", j, ".t", k, ".leg"), 
       ggplot(get(paste0("df",j,".t", k, ".leg")), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Backload)  ", "Simulated (Backload)  ", 
                                                  "95% Bootstrap CI (Backload)  ",
                                                  "Estimated (Equal)  ", "Simulated (Equal)  ", 
                                                  "95% Bootstrap CI (Equal)  "),
                            values = cbb[c(6, 7, 6, 4, 1, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Backload)  ", "Simulated (Backload)  ", 
                                                     "95% Bootstrap CI (Backload)  ",
                                                     "Estimated (Equal)  ", "Simulated (Equal)  ", 
                                                     "95% Bootstrap CI (Equal)  "),
                               values = rep(c("solid", "longdash", "dashed"), 2)) +
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
         guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
         ylim(0,1)) 


## arrange plot grid
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
                         layout_matrix = m, widths = c(1.2475, rep(1,4)), heights = c(rep(1, 3), 1.15))

## add legend to the bottom
fig_final <- plot_grid(fig.row1, get_legend(get(paste0("plot5.t5.leg"))), ncol = 1, rel_heights = c(1, .1))

# output as .pdf file for the article
pdf(file = paste0("Fig1Seq.pdf"),   # The directory you want to save the file in
    width = 1.25*12, # The width of the plot in inches
    height = 1.25*10) # The height of the plot in inches

fig_final

dev.off()

## get the sample sizes for the confirmatory sampling distribution estimates

## backloaded with T = 2
mat.back2 <- cbind(seq(100, 500, 1), lin.back2[1:401], lin.back2[402:802])
ind.back2 <- min(which(mat.back2[,3] >= 0.8))
n1.back2 <- mat.back2[ind.back2,1]

## Expected total sample size
n1.back2*mat.back2[ind.back2, 2] + 1.5*n1.back2*(1 - mat.back2[ind.back2, 2])

## backloaded with T = 3
mat.back3 <- cbind(seq(100, 500, 1), lin.back3[1:401], lin.back3[402:802], lin.back3[803:1203])
ind.back3 <- min(which(mat.back3[,4] >= 0.8))
n1.back3 <- mat.back3[ind.back3,1]

## Expected total sample size
n1.back3*mat.back3[ind.back3, 2] + 1.5*n1.back3*(mat.back3[ind.back3, 3] - mat.back3[ind.back3, 2]) +
  2*n1.back3*(1- mat.back3[ind.back3, 3])

## backloaded with T = 4
mat.back4 <- cbind(seq(100, 500, 1), lin.back4[1:401], lin.back4[402:802], lin.back4[803:1203], lin.back4[1204:1604])
ind.back4 <- min(which(mat.back4[,5] >= 0.8))
n1.back4 <- mat.back4[ind.back4,1]

## Expected total sample size
n1.back4*mat.back4[ind.back4, 2] + 1.5*n1.back4*(mat.back4[ind.back4, 3] - mat.back4[ind.back4, 2]) + 
  2*n1.back4*(mat.back4[ind.back4, 4] - mat.back4[ind.back4, 3]) +
  2.5*n1.back4*(1- mat.back4[ind.back4, 4])

## backloaded with T = 5
mat.back5 <- cbind(seq(100, 500, 1), lin.back5[1:401], lin.back5[402:802], lin.back5[803:1203], lin.back5[1204:1604],
                   lin.back5[1605:2005])
ind.back5 <- min(which(mat.back5[,6] >= 0.8))
n1.back5 <- mat.back5[ind.back5,1]

## Expected total sample size
n1.back5*mat.back5[ind.back5, 2] + 1.5*n1.back5*(mat.back5[ind.back5, 3] - mat.back5[ind.back5, 2]) + 
  2*n1.back5*(mat.back5[ind.back5, 4] - mat.back5[ind.back5, 3]) +
  2.5*n1.back5*(mat.back5[ind.back5, 5] - mat.back5[ind.back5, 4]) +
  3*n1.back5*(1- mat.back5[ind.back5, 5])

## equally spaced with T = 2
mat.eq2 <- cbind(seq(100, 500, 1), lin.eq2[1:401], lin.eq2[402:802])
ind.eq2 <- min(which(mat.eq2[,3] >= 0.8))
n1.eq2 <- mat.eq2[ind.eq2,1]

## Expected total sample size
n1.eq2*mat.eq2[ind.eq2, 2] + 2*n1.eq2*(1 - mat.eq2[ind.eq2, 2])

## equally spaced with T = 3
mat.eq3 <- cbind(seq(100, 500, 1), lin.eq3[1:401], lin.eq3[402:802], lin.eq3[803:1203])
ind.eq3 <- min(which(mat.eq3[,4] >= 0.8))
n1.eq3 <- mat.eq3[ind.eq3,1]

## Expected total sample size
n1.eq3*mat.eq3[ind.eq3, 2] + 2*n1.eq3*(mat.eq3[ind.eq3, 3] - mat.eq3[ind.eq3, 2]) +
  3*n1.eq3*(1- mat.eq3[ind.eq3, 3])

## equally spaced with T = 4
mat.eq4 <- cbind(seq(100, 500, 1), lin.eq4[1:401], lin.eq4[402:802], lin.eq4[803:1203], lin.eq4[1204:1604])
ind.eq4 <- min(which(mat.eq4[,5] >= 0.8))
n1.eq4 <- mat.eq4[ind.eq4,1]

## Expected total sample size
n1.eq4*mat.eq4[ind.eq4, 2] + 2*n1.eq4*(mat.eq4[ind.eq4, 3] - mat.eq4[ind.eq4, 2]) + 
  3*n1.eq4*(mat.eq4[ind.eq4, 4] - mat.eq4[ind.eq4, 3]) +
  4*n1.eq4*(1- mat.eq4[ind.eq4, 4])

## equally spaced with T = 5
mat.eq5 <- cbind(seq(100, 500, 1), lin.eq5[1:401], lin.eq5[402:802], lin.eq5[803:1203], lin.eq5[1204:1604],
                 lin.eq5[1605:2005])
ind.eq5 <- min(which(mat.eq5[,6] >= 0.8))
n1.eq5 <- mat.eq5[ind.eq5,1]

## Expected total sample size
n1.eq5*mat.eq5[ind.eq5, 2] + 2*n1.eq5*(mat.eq5[ind.eq5, 3] - mat.eq5[ind.eq5, 2]) + 
  3*n1.eq5*(mat.eq5[ind.eq5, 4] - mat.eq5[ind.eq5, 3]) +
  4*n1.eq5*(mat.eq5[ind.eq5, 5] - mat.eq5[ind.eq5, 4]) +
  5*n1.eq5*(1- mat.eq5[ind.eq5, 5])

## now get the confirmatory power and type I error rate estimate for the table
## in the main text

## implement for T = 2

## backloaded
ns <- 412
## confirm power
j <- 1
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 20)*100)
  write.csv(probs.temp, paste0("confirm_probs_T2_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## confirm the type I error rate
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 40)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T2_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## also implement for equally spaced
ns <- 319
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
  write.csv(probs.temp, paste0("confirm_probs_T2_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
}

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
  write.csv(probs.temp, paste0("confirm_probs_H0_T2_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## implement for T = 3

## backloaded
ns <- 306
## confirm power
j <- 1
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 100)*100)
  write.csv(probs.temp, paste0("confirm_probs_T3_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 120)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T3_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## also implement for equally spaced
ns <- 208
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
                                 seed = (i + 140)*100)
  write.csv(probs.temp, paste0("confirm_probs_T3_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
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
                                 seed = (i + 160)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T3_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## implement for T = 4
ns <- 234
## confirm power
j <- 1
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2, 2.5), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 180)*100)
  write.csv(probs.temp, paste0("confirm_probs_T4_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2, 2.5), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 200)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T4_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## also implement for equally spaced
ns <- 161
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
                                 seed = (i + 220)*100)
  write.csv(probs.temp, paste0("confirm_probs_T4_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
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
                                 seed = (i + 240)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T4_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## implement for T = 5

## backloaded
ns <- 197
## confirm power
j <- 1
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2, 2.5, 3), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 260)*100)
  write.csv(probs.temp, paste0("confirm_probs_T5_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2, 2.5, 3), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 1, ndraws = 10000,
                                 seed = (i + 280)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T5_back_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## also implement for equally spaced
ns <- 125
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
                                 seed = (i + 300)*100)
  write.csv(probs.temp, paste0("confirm_probs_T5_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
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
                                 seed = (i + 320)*100)
  write.csv(probs.temp, paste0("confirm_probs_H0_T5_equal_scen1_", ns[i], ".csv"), row.names = FALSE)
}



## this function is used to get the confirmatory estimates of the stopping probabilities
confirm_probs <- function(mat, gam){
  ## mat is a matrix of complementary probabilities (one column for each analysis)
  ## gam are the decision tresholds
  mat <- 1 - mat
  res <- NULL
  stop.temp <- mat[,1] >= gam[1]
  res[1] <- mean(stop.temp)
  for (j in 2:ncol(mat)){
    
    stop.temp <- ifelse(stop.temp, 1, ifelse(mat[,j] >= gam[j], 1, 0))
    
    res[j] <- mean(stop.temp)
  }
  return(res)
}

## implement for backloaded with T = 2
## type I error
confirm_probs(read.csv(paste0("confirm_probs_H0_T2_back_scen1_412.csv")), 
              gamma.back2)

## power
confirm_probs(read.csv(paste0("confirm_probs_T2_back_scen1_412.csv")), 
              gamma.back2)

## implement for backloaded with T = 3
confirm_probs(read.csv(paste0("confirm_probs_H0_T3_back_scen1_306.csv")), 
              gamma.back3)

confirm_probs(read.csv(paste0("confirm_probs_T3_back_scen1_306.csv")), 
              gamma.back3)

## implement for backloaded with T = 4
confirm_probs(read.csv(paste0("confirm_probs_H0_T4_back_scen1_234.csv")), 
              gamma.back4)

confirm_probs(read.csv(paste0("confirm_probs_T4_back_scen1_234.csv")), 
              gamma.back4)

## implement for backloaded with T = 5
confirm_probs(read.csv(paste0("confirm_probs_H0_T5_back_scen1_197.csv")), 
              gamma.back5)

confirm_probs(read.csv(paste0("confirm_probs_T5_back_scen1_197.csv")), 
              gamma.back5)

## implement for equally spaced with T = 2
confirm_probs(read.csv(paste0("confirm_probs_H0_T2_equal_scen1_319.csv")), 
              gamma.eq2)

confirm_probs(read.csv(paste0("confirm_probs_T2_equal_scen1_319.csv")), 
              gamma.eq2)

## implement for equally spaced with T = 3
confirm_probs(read.csv(paste0("confirm_probs_H0_T3_equal_scen1_208.csv")), 
              gamma.eq3)

confirm_probs(read.csv(paste0("confirm_probs_T3_equal_scen1_208.csv")), 
              gamma.eq3)

## implement for equally spaced with T = 4
confirm_probs(read.csv(paste0("confirm_probs_H0_T4_equal_scen1_161.csv")), 
              gamma.eq4)

confirm_probs(read.csv(paste0("confirm_probs_T4_equal_scen1_161.csv")), 
              gamma.eq4)

## implement for equally spaced with T = 5
confirm_probs(read.csv(paste0("confirm_probs_H0_T5_equal_scen1_125.csv")), 
              gamma.eq5)

confirm_probs(read.csv(paste0("confirm_probs_T5_equal_scen1_125.csv")), 
              gamma.eq5)