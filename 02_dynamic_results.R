## this code file contains the code to reproduce the numerical results and
## and plots in Section 6.1 of the main text. Please run "01_dynamic_functions.R" first
## to ensure that the necessary packages and functions are loaded

## load the packages
require(cowplot)
require(ggpubr)
require(mvtnorm)
require(ggplot2)
require(foreach)
require(doParallel)
require(doSNOW)
require(truncnorm)
source("01_dynamic_functions.R")

cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get parameter values for various binomial distributions;
## generate these parameter values for the sample sizes that will be explored later in
## this file
ns <- seq(100, 300, 20)

## binomial summaries for the conditional approach
seeds <- seq(1, length(ns), 1)
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.8), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen1_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 1000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.7), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen1_", ns[k], ".csv"), row.names = FALSE)
}

## scenario 2
seeds <- seq(1, length(ns), 1) + 2000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.75, 0.843), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen2_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 3000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.75, 0.75), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen2_", ns[k], ".csv"), row.names = FALSE)
}

## scenario 3
seeds <- seq(1, length(ns), 1) + 4000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.8), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen3_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 5000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.7), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen3_", ns[k], ".csv"), row.names = FALSE)
}

## scenario 4
seeds <- seq(1, length(ns), 1) + 6000
p1s <- sort(rtruncnorm(10000, 0.75, 0.85, 0.8, 0.03))
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummariesPred(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), cbind(0.7, p1s), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen4_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 7000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.7), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen4_", ns[k], ".csv"), row.names = FALSE)
}

## simulations for scenario 1 
ll <- "H1_scen1"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

ll <- "H0_scen1"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## simulations for scenario 2 
ll <- "H1_scen2"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

ll <- "H0_scen2"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## simulations for scenario 3
## switch posterior approximation to the commented out version of 
## getPostBin() in "01_dynamic_functions.R"
ll <- "H1_scen3"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

ll <- "H0_scen3"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## simulations for scenario 4
ll <- "H1_scen4"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

ll <- "H0_scen4"
c_vec <- c(1, 2, 3)
map.par <- c(134, 53)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("bin_summary_", ll, "_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:2){
    ni.temp <- ceiling(c_vec[i]*ns[k])
    nf.temp <- ceiling(tail(c_vec, 1)*ns[k])
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
      temp <- getCompPred(ni.temp, nf.temp, summaries = summary.tempj, deltaL = 0, mappar = map.par,
                          ratio = 0.5, ndraws = 5000, M = 1000)
      write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 3
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), seq(4*(i-1) + 1, 4*i, 1)]
    temp <- getCompPost(summaries = summary.tempj, deltaL = 0, mappar = map.par,
                        ratio = 0.5, ndraws = 5000, M = 1000)
    write.csv(temp, paste0(ll, "_n_",ns[k],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## get the stopping probabilities and plots
c_vec <- c(1, 2, 3)
a.eq3 <- c(0.0013, 0.0061, 0.025)

## construct sampling distribution of posterior probabilities

## combine the batches
ll <- "H1_scen1"
ns <- seq(100,300,20)
for (k in c(2, 10)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))[,1]
    print(c(j, 2))
    temp.res2 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))[,1]
    temp.res3 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 3, "_batch_", j, ".csv"))
    full <- rbind(full, cbind(temp.res1, temp.res2, temp.res3))
  }
  write.csv(full, paste0(ll, "_n_",ns[k],".csv"), row.names = FALSE)
}

ll <- "H0_scen1"
ns <- seq(100,300,20)
for (k in c(2, 10)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))[,1]
    print(c(j, 2))
    temp.res2 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))[,1]
    temp.res3 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 3, "_batch_", j, ".csv"))
    full <- rbind(full, cbind(temp.res1, temp.res2, temp.res3))
  }
  write.csv(full, paste0(ll, "_n_",ns[k],".csv"), row.names = FALSE)
}

ll <- "H1_scen1"
probs1 <- read.csv(paste0(ll, "_n_",120,".csv"))
probs2 <- read.csv(paste0(ll, "_n_",280,".csv"))

ll <- "H0_scen1"
probs10 <- read.csv(paste0(ll, "_n_",120,".csv"))
probs20 <- read.csv(paste0(ll, "_n_",280,".csv"))

stop.lin1 = stop_mat(m1 = probs1, m2 = probs2,
                     m10 = probs10, m20 = probs20, 120*c_vec,
                     280*c_vec, 100, 300, by = 1, a.eq3)

gam.sub1 <- stop.lin1[[3]]

## construct the sampling distributions of posterior summaries under H0
## and H1 at na and nb (for every sample size since the threshold changes);
## in particular, we extract the posterior probabilities resulting from the 
## posterior predictive distribution at analyses 1 and 2 used to compute the
## posterior predictive probabilities
ll <- "H1_scen1"
ns <- seq(100,300,20)
for (k in c(2, 10)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))[,2:1001]
    full <- rbind(full, temp.res1)
  }
  assign(paste0("H1_n_",ns[k],"_t_", 1), full)
}

for (k in c(2, 10)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))[,2:1001]
    full <- rbind(full, temp.res1)
  }
  assign(paste0("H1_n_",ns[k],"_t_", 2), full)
}

ll <- "H0_scen1"
for (k in c(2, 10)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))[,2:1001]
    full <- rbind(full, temp.res1)
  }
  assign(paste0("H0_n_",ns[k],"_t_", 1), full)
}

for (k in c(2, 10)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))[,2:1001]
    full <- rbind(full, temp.res1)
  }
  assign(paste0("H0_n_",ns[k],"_t_", 2), full)
}

## extract the two sampling distribution estimates used for constructing linear approximations
ll <- "H1_scen1"
H1_n_120 <- read.csv(paste0(ll, "_n_",120,".csv"))
H1_n_280 <- read.csv(paste0(ll, "_n_",280,".csv"))

ll <- "H0_scen1"
H0_n_120 <- read.csv(paste0(ll, "_n_",120,".csv"))
H0_n_280 <- read.csv(paste0(ll, "_n_",280,".csv"))

## get the stopping probabilities and thresholds for all n using linear approximations
ns_vec <- seq(100, 300, 1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = length(ns_vec), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

mat.temp <- foreach(j=1:length(ns_vec), .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      ## get the posterior predictive probabilities under H1 at analyses 1 and 2 using
                      ## the sample-size specific threshold
                      pp.temp1 <- apply(H1_n_120_t_1, 1, ppp.dens.logit, gam = gam.sub1[j, 3])
                      pp.temp1 <- cbind(pp.temp1, apply(H1_n_120_t_2, 1, ppp.dens.logit, gam = gam.sub1[j, 3]))
                      pp.temp2 <- apply(H1_n_280_t_1, 1, ppp.dens.logit, gam = gam.sub1[j, 3])
                      pp.temp2 <- cbind(pp.temp2, apply(H1_n_280_t_2, 1, ppp.dens.logit, gam = gam.sub1[j, 3]))
                      
                      ## get the posterior predictive probabilities under H0 at analyses 1 and 2 using
                      ## the sample-size specific threshold
                      pp.temp1.0 <- apply(H0_n_120_t_1, 1, ppp.dens.logit, gam = gam.sub1[j, 3])
                      pp.temp1.0 <- cbind(pp.temp1.0, apply(H0_n_120_t_2, 1, ppp.dens.logit, gam = gam.sub1[j, 3]))
                      pp.temp2.0 <- apply(H0_n_280_t_1, 1, ppp.dens.logit, gam = gam.sub1[j, 3])
                      pp.temp2.0 <- cbind(pp.temp2.0, apply(H0_n_280_t_2, 1, ppp.dens.logit, gam = gam.sub1[j, 3]))
                      
                      ## tune the failure thresholds
                      xi <- as.numeric(stop_mat(pp.temp1.0, pp.temp2.0, pp.temp1, pp.temp2, 
                                                 120*c_vec, 280*c_vec, ns_vec[j], ns_vec[j], 
                                                 by = 1, c(0.01, 0.05))[[3]])
                      
                      ## given thresholds, get linear approximations under H1
                      lines.H1 <- getLines(cbind(H1_n_120[,1], pp.temp1[,1], H1_n_120[,2], 
                                                 pp.temp1[,2], H1_n_120[,3]),
                                           cbind(H1_n_280[,1], pp.temp2[,1], H1_n_280[,2], 
                                                 pp.temp2[,2], H1_n_280[,3]),
                                           120*c(1, 1, 2, 2, 3), 280*c(1, 1, 2, 2, 3))
                      
                      ## given linear approximations under H1, get stopping probabilities
                      stop.H1 <- stopMatBoth(lines.H1, c(1, 1, 2, 2, 3), 
                                             ns_vec[j], ns_vec[j], 1, gam.sub1[j, ], xi)
                      
                      ## given thresholds, get linear approximations under H0
                      lines.H0 <- getLines(cbind(H0_n_120[,1], pp.temp1.0[,1], H0_n_120[,2], 
                                                 pp.temp1.0[,2], H0_n_120[,3]),
                                           cbind(H0_n_280[,1], pp.temp2.0[,1], H0_n_280[,2], 
                                                 pp.temp2.0[,2], H0_n_280[,3]),
                                           120*c(1, 1, 2, 2, 3), 280*c(1, 1, 2, 2, 3))
                      
                      ## given linear approximations under H0, get stopping probabilities
                      stop.H0 <- stopMatBoth(lines.H0, c(1, 1, 2, 2, 3), 
                                             ns_vec[j], ns_vec[j], 1, gam.sub1[j, ], xi)
                      
                      ## return stopping probabilities and recommended decision thresholds gamma and xi (mixed by analysis number)
                      c(stop.H0, stop.H1, c(gam.sub1[j, 1], xi[1], gam.sub1[j, 2], xi[2], gam.sub1[j, 3]))
                    }

## output results to .csv file
write.csv(mat.temp, "mat3_scen1_lin.csv", row.names = FALSE)

## now repeat the process using the sampling distribution estimates at the different sample sizes
registerDoSNOW(cl)
pb <- txtProgressBar(max = length(ns), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## extract the decision thresholds at the 11 values considered using simulation
thres.ns <- mat.temp[seq(1, 201, 20),11:15]

mat.temp.dat <- foreach(k=1:length(ns), .combine=rbind,
                        .options.snow=opts) %dopar% {
                          
                          ll <- "H1_scen1"
                          
                          ## get the two sampling distribution estimates under H1 (na = 120 and nb = 280)
                          full <- NULL
                          pp <- NULL
                          for (j in 1:10){
                            print(c(j, 1))
                            temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))
                            full <- rbind(full, temp.res1[, 2:1001])
                            pp <- c(pp, temp.res1[,1])
                          }
                          assign(paste0("H1_n_",ns[k],"_t_", 1), full)
                          assign(paste0("pp_H1_n_",ns[k],"_t_", 1), pp)
                          
                          
                          full <- NULL
                          pp <- NULL
                          for (j in 1:10){
                            print(c(j, 1))
                            temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))
                            full <- rbind(full, temp.res1[, 2:1001])
                            pp <- c(pp, temp.res1[,1])
                          }
                          assign(paste0("H1_n_",ns[k],"_t_", 2), full)
                          assign(paste0("pp_H1_n_",ns[k],"_t_", 2), pp)
                          
                          pp <- NULL
                          for (j in 1:10){
                            print(c(j, 1))
                            temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 3, "_batch_", j, ".csv"))
                            pp <- c(pp, temp.res1[,1])
                          }
                          assign(paste0("pp_H1_n_",ns[k],"_t_", 3), pp)
                          
                          ## get the two sampling distribution estimates under H0 (na = 120 and nb = 280)
                          ll <- "H0_scen1"
                          full <- NULL
                          pp <- NULL
                          for (j in 1:10){
                            print(c(j, 1))
                            temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))
                            full <- rbind(full, temp.res1[, 2:1001])
                            pp <- c(pp, temp.res1[,1])
                          }
                          assign(paste0("H0_n_",ns[k],"_t_", 1), full)
                          assign(paste0("pp_H0_n_",ns[k],"_t_", 1), pp)
                          
                          full <- NULL
                          pp <- NULL
                          for (j in 1:10){
                            print(c(j, 1))
                            temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))
                            full <- rbind(full, temp.res1[, 2:1001])
                            pp <- c(pp, temp.res1[,1])
                          }
                          assign(paste0("H0_n_",ns[k],"_t_", 2), full)
                          assign(paste0("pp_H0_n_",ns[k],"_t_", 2), pp)
                          
                          pp <- NULL
                          for (j in 1:10){
                            print(c(j, 1))
                            temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 3, "_batch_", j, ".csv"))
                            pp <- c(pp, temp.res1[,1])
                          }
                          assign(paste0("pp_H0_n_",ns[k],"_t_", 3), pp)
                          
                          ## use the recommended decision thresholds to compute posterior predictive probabilities under H1
                          pp.temp1 <- apply(get(paste0("H1_n_",ns[k],"_t_", 1)), 1, ppp.prob, gam = thres.ns[k, 5])
                          pp.temp1 <- cbind(pp.temp1, apply(get(paste0("H1_n_",ns[k],"_t_", 2)), 1, ppp.prob, gam = thres.ns[k, 5]))
                          
                          ## use the recommended decision thresholds to compute posterior predictive probabilities under H0
                          pp.temp1.0 <- apply(get(paste0("H0_n_",ns[k],"_t_", 1)), 1, ppp.prob, gam = thres.ns[k, 5])
                          pp.temp1.0 <- cbind(pp.temp1.0, apply(get(paste0("H0_n_",ns[k],"_t_", 2)), 1, ppp.prob, gam = thres.ns[k, 5]))
                          
                          ## use the sampling distribution estimates at the explored sample sizes to get the stopping probabilities under
                          ## H1 and H0
                          stop.H1 <- stopData(cbind(get(paste0("pp_H1_n_",ns[k],"_t_", 1)),
                                                    pp.temp1[,1], get(paste0("pp_H1_n_",ns[k],"_t_", 2)),
                                                    pp.temp1[,2], get(paste0("pp_H1_n_",ns[k],"_t_", 3))), 
                                              thres.ns[k, c(1, 3, 5)], thres.ns[k, c(2,4)])
                          
                          stop.H0 <- stopData(cbind(get(paste0("pp_H0_n_",ns[k],"_t_", 1)),
                                                    pp.temp1.0[,1], get(paste0("pp_H0_n_",ns[k],"_t_", 2)),
                                                    pp.temp1.0[,2], get(paste0("pp_H0_n_",ns[k],"_t_", 3))), 
                                              thres.ns[k, c(1, 3, 5)], thres.ns[k, c(2,4)])
                          
                          ## return stopping probabilities
                          c(stop.H0, stop.H1)
                        }

## return results
write.csv(mat.temp.dat, "mat3_scen1_dat.csv", row.names = FALSE)


## repeat the process in lines 352 to 618 for scenarios 2, 3, and 4
## it is only necessary to replace ll <- "H1_scen1" or ll <- "H1_scen1" with
## the correct scenario number and change the labels on the .csv files "mat3_scen1_lin.csv"
## and "mat3_scen1_dat.csv" to match the correct scenario number as well

## the following code is used to get the stopping probabilities for scenario 5, which
## is simpler since posterior predictive probabilities are not involved

## construct sampling distribution estimates under H1 and H0
## using the results from scenario 1 (same data-generating parameters and priors;
## just a different decision rule)
ll <- "H1_scen1"
ns <- seq(100,300,20)
for (k in c(1:11)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))[,1]
    print(c(j, 2))
    temp.res2 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))[,1]
    temp.res3 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 3, "_batch_", j, ".csv"))
    full <- rbind(full, cbind(temp.res1, temp.res2, temp.res3))
  }
  write.csv(full, paste0(ll, "_n_",ns[k],".csv"), row.names = FALSE)
}

ll <- "H0_scen1"
ns <- seq(100,300,20)
for (k in c(1:11)){
  full <- NULL
  for (j in 1:10){
    print(c(j, 1))
    temp.res1 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 1, "_batch_", j, ".csv"))[,1]
    print(c(j, 2))
    temp.res2 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 2, "_batch_", j, ".csv"))[,1]
    temp.res3 <- read.csv(paste0(ll, "_n_",ns[k],"_t_", 3, "_batch_", j, ".csv"))
    full <- rbind(full, cbind(temp.res1, temp.res2, temp.res3))
  }
  write.csv(full, paste0(ll, "_n_",ns[k],".csv"), row.names = FALSE)
}

## extract the results for the two sample sizes used to construct linear approximations
ll <- "H1_scen1"
probs1 <- read.csv(paste0(ll, "_n_",120,".csv"))
probs2 <- read.csv(paste0(ll, "_n_",280,".csv"))

ll <- "H0_scen1"
probs10 <- read.csv(paste0(ll, "_n_",120,".csv"))
probs20 <- read.csv(paste0(ll, "_n_",280,".csv"))

stop.lin1 = stop_mat(m1 = probs1, m2 = probs2,
                     m10 = probs10, m20 = probs20, 120*c_vec,
                     280*c_vec, 100, 300, by = 1, a.eq3)

## get results to compare to other scenarios (add 0 for the failure thresholds)
stop.H1.temp <- cbind(stop.lin1[[1]], 0, 0)
stop.H0.temp <- cbind(stop.lin1[[2]], 0, 0)
thres.temp <- cbind(stop.lin1[[3]], 0, 0)

stop.H1.temp <- stop.H1.temp[, c(1, 4, 2, 5, 3)]
stop.H0.temp <- stop.H0.temp[, c(1, 4, 2, 5, 3)]
thres.temp <- thres.temp[, c(1, 4, 2, 5, 3)]

## output results for the linear approximations
write.csv(cbind(stop.H0.temp, stop.H1.temp, thres.temp), 
          "mat3_scen5_lin.csv", row.names = FALSE)

## extract recommended success thresholds for sample sizes 
## considered via simulation
thres.dat <- thres.temp[seq(1, 201, 20), c(1, 3, 5)]

## get power estimates (under H1) from data
n_vals <- seq(100, 300, 20)
sim.pwr1 <- NULL
for (kk in 1:length(n_vals)){
  probs.temp <- read.csv(paste0("H1_scen1_n_", n_vals[kk], ".csv"))
  sim.pwr1 <- rbind(sim.pwr1,
                    stop_dat(probs.temp, thres.dat[kk,]))
}

## add failure stopping probabilities of 0 for first two analyses
sim.pwr1 <- cbind(sim.pwr1, 0, 0)
sim.pwr1 <- sim.pwr1[, c(1, 4, 2, 5, 3)]

## now consider the type I error rate (under H0) from data
n_vals <- seq(100, 300, 20)
sim.t1e1 <- NULL
for (kk in 1:length(n_vals)){
  probs.temp <- read.csv(paste0("H0_scen1_n_", n_vals[kk], ".csv"))
  sim.t1e1 <- rbind(sim.t1e1,
                    stop_dat(probs.temp, thres.dat[kk,]))
}

## add failure stopping probabilities of 0 for first two analyses
sim.t1e1 <- cbind(sim.t1e1, 0, 0)
sim.t1e1 <- sim.t1e1[, c(1, 4, 2, 5, 3)]

## save results to a .csv file
write.csv(cbind(sim.t1e1, sim.pwr1), 
          "mat3_scen5_dat.csv", row.names = FALSE)


## construct Figures 1 and 2 in Section 6.1
mat.temp <- read.csv("mat3_scen1_lin.csv")
mat.temp.dat <- read.csv("mat3_scen1_dat.csv")

mat.temp5 <- read.csv("mat3_scen5_lin.csv")
mat.temp.dat5 <- read.csv("mat3_scen5_dat.csv")

mat.temp2 <- read.csv("mat3_scen2_lin.csv")
mat.temp.dat2 <- read.csv("mat3_scen2_dat.csv")

mat.temp3 <- read.csv("mat3_scen3_lin.csv")
mat.temp.dat3 <- read.csv("mat3_scen3_dat.csv")

mat.temp4 <- read.csv("mat3_scen4_lin_prec.csv")
mat.temp.dat4 <- read.csv("mat3_scen4_dat_prec.csv")

## hard code c vector, sample sizes,
## upper and lower bounds for each subplot and the y-axis plot labels
cc_vec = rep(c(1, 1, 2, 2, 3), 3)
ns <- seq(100, 300, 20)
lbs <- c(0, 0, 0, 0, 0, 
         0, 0, 0, 0, 0,
         0.95, 0, 0.95, 0, 0.95)
ubs <- c(0.05, 1, 0.05, 1, 0.05,
         1, 0.1, 1, 0.1, 1,
         1, 0.35, 1, 0.35, 1)

ylabs <- c(expression("Success Probability ("*Psi[0]*")"),
  expression("Failure Probability ("*Psi[0]*")"),
  expression("Success Probability ("*Psi[0]*")"),
  expression("Failure Probability ("*Psi[0]*")"),
  expression("Success Probability ("*Psi[0]*")"),
  
  expression("Success Probability ("*Psi[1]*")"),
  expression("Failure Probability ("*Psi[1]*")"),
  expression("Stopping Probability ("*Psi[1]*")"),
  expression("Failure Probability ("*Psi[1]*")"),
  expression("Success Probability ("*Psi[1]*")"),
  
  expression("Success Threshold ("*gamma['t']*")"),
  expression("Failure Threshold ("*rho['t']*")"),
  expression(gamma[2]),
  expression(rho[2]),
  expression(gamma[3]))

## get the correct order of the plots for Figure 1 (success only)
for (kk in c(6, 8, 10, 1, 3, 5, 11, 13, 15)){
  
  ## get n-values for that subplot
  n_full <- cc_vec[kk]*seq(100, 300, 1)
  n_dat  <- cc_vec[kk]*ns
  
  ## for the stopping probabilities, we have "simulated" and "estimated" results
  if (kk < 11){
    assign(paste0("df", kk, "1"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp[,kk], mat.temp.dat[,kk]),
                      design = c(rep("D1", length(n_full)), rep("ZD1", length(n_dat)))))
    
    assign(paste0("df", kk, "2"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp2[,kk], mat.temp.dat2[,kk]),
                      design = c(rep("D2", length(n_full)), rep("ZD2", length(n_dat)))))
    
    assign(paste0("df", kk, "3"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp3[,kk], mat.temp.dat3[,kk]),
                      design = c(rep("D3", length(n_full)), rep("ZD3", length(n_dat)))))
    
    assign(paste0("df", kk, "4"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp4[,kk], mat.temp.dat4[,kk]),
                      design = c(rep("D4", length(n_full)), rep("ZD4", length(n_dat)))))
    
    assign(paste0("df", kk, "5"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp5[,kk], mat.temp.dat5[,kk]),
                      design = c(rep("D5", length(n_full)), rep("ZD5", length(n_dat)))))
    
    assign(paste0("df", kk),
           rbind(get(paste0("df", kk, "1")), get(paste0("df", kk, "2")),
                 get(paste0("df", kk, "3")), get(paste0("df", kk, "4")),
                 get(paste0("df", kk, "5"))))
    
    assign(paste0("plot", kk),
           ggplot(get(paste0("df", kk)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = y, color=as.factor(design), linetype = as.factor(design), alpha = as.factor(design)), 
                       size = 1) +
             labs(title='') +
             labs(x= bquote(italic(n)[.(cc_vec[kk])]), y= ylabs[kk]) +
             theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c("firebrick1", "steelblue1", "seagreen1", "darkorchid1", "orange1",
                                           "firebrick4", "steelblue4", "seagreen4", "darkorchid4", "orange4")) +
             scale_linetype_manual(name = " ", 
                                   values = c(rep(1, 5), rep(2, 5))) +
             scale_alpha_manual(name = " ", 
                                values = c(rep(0.7, 5), rep(0.9, 5))) +
             theme(legend.text=element_text(size=18)) +
             theme(legend.key.size = unit(1, "cm")) +
             ylim(lbs[kk],ubs[kk]) +
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))
  }
  
  ## we only have the "estimated" results for the recommended decision thresholds
  if (kk > 10){
    assign(paste0("df", kk, "1"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp[,kk]),
                      design = c(rep("D1", length(n_full)))))
    
    assign(paste0("df", kk, "2"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp2[,kk]),
                      design = c(rep("D2", length(n_full)))))
    
    assign(paste0("df", kk, "3"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp3[,kk]),
                      design = c(rep("D3", length(n_full)))))
    
    assign(paste0("df", kk, "4"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp4[,kk]),
                      design = c(rep("D4", length(n_full)))))
    
    assign(paste0("df", kk, "5"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp5[,kk]),
                      design = c(rep("D5", length(n_full)))))
    
    assign(paste0("df", kk),
           rbind(get(paste0("df", kk, "1")), get(paste0("df", kk, "2")),
                 get(paste0("df", kk, "3")), get(paste0("df", kk, "4")),
                 get(paste0("df", kk, "5"))))
    
    assign(paste0("plot", kk),
           ggplot(get(paste0("df", kk)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = y, color=as.factor(design), linetype = as.factor(design)), 
                       alpha = 0.7, size = 1) +
             labs(title='') +
             labs(x= bquote(italic(n)[.(cc_vec[kk])]), y= ylabs[kk]) +
             theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c("firebrick1", "steelblue1", "seagreen1", "darkorchid1", "orange1")) +
             scale_linetype_manual(name = " ", 
                                   values = c(rep(1, 5))) +
             theme(legend.text=element_text(size=18)) +
             theme(legend.key.size = unit(1, "cm")) +
             ylim(lbs[kk],ubs[kk]) +
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))
  }
  
  ## remove x-axis for plots not the bottom row
  if (kk %in% c(6, 8, 10, 1, 3, 5)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + xlab('') + theme(axis.text.x = element_blank(),
                                                      axis.ticks.x = element_blank()))
  }
  
  ## remove y-axis for plots not the left column
  if (kk %in% c(8, 10, 3, 5, 13, 15)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + ylab('') + theme(axis.text.y = element_blank(),
                                                      axis.ticks.y = element_blank()))
  }
  
  ## add titles to the top row
  if (kk %in% c(6, 8, 10)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + labs(title=bquote('Analysis '*.(cc_vec[kk]))) +
             theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                             margin = margin(t = 0, 0, 5, 0))))
  }
  
  ## force x-axis scaling for the right column
  if (kk %in% c(10, 5, 15)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + scale_x_continuous(breaks = seq(300, 900, 150)))
  }
  
}

## get the first legend (colours) to be common for all plots
kk <- 6
assign(paste0("plot", kk, ".leg1"),
       ggplot(subset(get(paste0("df", kk)),
                     get(paste0("df", kk))$design %in% c("D1", "D2", "D3", "D4", "D5")), aes(x=n)) + theme_bw() +
         geom_line(aes(y = y, color=as.factor(design), linetype = as.factor(design), alpha = as.factor(design)), 
                   size = 1) +
         labs(title='') +
         labs(x= bquote(italic(n)[.(cc_vec[kk])]), y= ylabs[kk]) +
         theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="bottom") +
         scale_color_manual(name = " ", 
                            values = c("firebrick1", "steelblue1", "seagreen1", "darkorchid1", "orange1",
                                       "firebrick4", "steelblue4", "seagreen4", "darkorchid4", "orange4")) +
         scale_linetype_manual(name = " ", 
                               values = c(rep(1, 5), rep(2, 5))) +
         scale_alpha_manual(name = " ", 
                            values = c(rep(0.7, 5), rep(0.9, 5))) +
         theme(legend.text=element_text(size=18)) +
         theme(legend.key.size = unit(1, "cm")) +
         ylim(lbs[kk],ubs[kk]) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))

plot6.leg1 

## get the second legend (linetype) to be common for all plots
assign(paste0("plot", kk, ".leg2"),
       ggplot(subset(get(paste0("df", kk)),
                     get(paste0("df", kk))$design %in% c("D1", "D2")), aes(x=n)) + theme_bw() +
         geom_line(aes(y = y, linetype = as.factor(design)), 
                   size = 1, col = "black") +
         labs(title='') +
         labs(x= bquote(italic(n)[.(cc_vec[kk])]), y= ylabs[kk]) +
         theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="bottom") +
         scale_linetype_manual(name = " ", 
                               labels = c("Estimated  ", "Simulated  "),
                               values = c(rep(1, 1), rep(2, 1))) +
         theme(legend.text=element_text(size=18)) +
         theme(legend.key.size = unit(1.25, "cm")) +
         ylim(lbs[kk],ubs[kk]) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))

plot6.leg2 

## construct plot into a 3X3 grid
figp.row1 <- plot_grid(plot6 + theme(plot.margin=unit(c(0.1, 0.15, -0.6, 0.1),"cm")), 
                       plot8 + theme(plot.margin=unit(c(0.1, 0.15, -0.6, -0.6),"cm")),
                       plot10 + theme(plot.margin=unit(c(0.1, 0.15, -0.6, -0.6),"cm")),
                       rel_widths = c(1.255, 1, 1), nrow = 1)
figp.row2 <- plot_grid(plot1 + theme(plot.margin=unit(c(-0.6, 0.15, -0.6, 0.1),"cm")), 
                       plot3 + theme(plot.margin=unit(c(-0.6, 0.15, -0.6, -0.6),"cm")),
                       plot5 + theme(plot.margin=unit(c(-0.6, 0.15, -0.6, -0.6),"cm")),
                       rel_widths = c(1.255, 1, 1), nrow = 1)
figp.row3 <- plot_grid(plot11 + theme(plot.margin=unit(c(-0.6, 0.15, 0.1, 0.1),"cm")), 
                       plot13 + theme(plot.margin=unit(c(-0.6, 0.15, 0.1, -0.6),"cm")),
                       plot15 + theme(plot.margin=unit(c(-0.6, 0.15, 0.1, -0.6),"cm")),
                       rel_widths = c(1.255, 1, 1), nrow = 1)

## add legends and output figure to PDF file
figp <- plot_grid(figp.row1, 
                  figp.row2, 
                  figp.row3, 
                  ggpubr::get_legend(plot6.leg1),
                  ggpubr::get_legend(plot6.leg2),
                  ncol = 1, rel_heights = c(1.105, 1, 1.2, 0.15, 0.15))

# output as .pdf file for the article
pdf(file = "FigEff.pdf",   # The directory you want to save the file in
    width =10.5, # The width of the plot in inches
    height = 10.5) # The height of the plot in inches

figp

dev.off()


## now repeat for Figure 2 (failure only)
for (kk in c(7, 9, 2, 4, 12, 14)){
  
  ## get n-values
  n_full <- cc_vec[kk]*seq(100, 300, 1)
  n_dat  <- cc_vec[kk]*ns
  
  if (kk < 11){
    ## only designs 1 to 4 in this plot
    assign(paste0("df", kk, "1"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp[,kk], mat.temp.dat[,kk]),
                      design = c(rep("D1", length(n_full)), rep("ZD1", length(n_dat)))))
    
    assign(paste0("df", kk, "2"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp2[,kk], mat.temp.dat2[,kk]),
                      design = c(rep("D2", length(n_full)), rep("ZD2", length(n_dat)))))
    
    assign(paste0("df", kk, "3"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp3[,kk], mat.temp.dat3[,kk]),
                      design = c(rep("D3", length(n_full)), rep("ZD3", length(n_dat)))))
    
    assign(paste0("df", kk, "4"), 
           data.frame(n = c(n_full, n_dat),
                      y = c(mat.temp4[,kk], mat.temp.dat4[,kk]),
                      design = c(rep("D4", length(n_full)), rep("ZD4", length(n_dat)))))
    
    assign(paste0("df", kk),
           rbind(get(paste0("df", kk, "1")), get(paste0("df", kk, "2")),
                 get(paste0("df", kk, "3")), get(paste0("df", kk, "4"))))
    
    assign(paste0("plot", kk),
           ggplot(get(paste0("df", kk)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = y, color=as.factor(design), linetype = as.factor(design), alpha = as.factor(design)), 
                       size = 1) +
             labs(title='') +
             labs(x= bquote(italic(n)[.(cc_vec[kk])]), y= ylabs[kk]) +
             theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c("firebrick1", "steelblue1", "seagreen1", "darkorchid1",
                                           "firebrick4", "steelblue4", "seagreen4", "darkorchid4")) +
             scale_linetype_manual(name = " ", 
                                   values = c(rep(1, 4), rep(2, 4))) +
             scale_alpha_manual(name = " ", 
                                values = c(rep(0.7, 5), rep(0.9, 5))) +
             theme(legend.text=element_text(size=18)) +
             theme(legend.key.size = unit(1, "cm")) +
             ylim(lbs[kk],ubs[kk]) +
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))
  }
  
  ## no "simulated" results for decision thresholds
  if (kk > 10){
    assign(paste0("df", kk, "1"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp[,kk]),
                      design = c(rep("D1", length(n_full)))))
    
    assign(paste0("df", kk, "2"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp2[,kk]),
                      design = c(rep("D2", length(n_full)))))
    
    assign(paste0("df", kk, "3"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp3[,kk]),
                      design = c(rep("D3", length(n_full)))))
    
    assign(paste0("df", kk, "4"), 
           data.frame(n = c(n_full),
                      y = c(mat.temp4[,kk]),
                      design = c(rep("D4", length(n_full)))))
    
    assign(paste0("df", kk),
           rbind(get(paste0("df", kk, "1")), get(paste0("df", kk, "2")),
                 get(paste0("df", kk, "3")), get(paste0("df", kk, "4"))))
    
    assign(paste0("plot", kk),
           ggplot(get(paste0("df", kk)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = y, color=as.factor(design), linetype = as.factor(design)), 
                       alpha = 0.7, size = 1) +
             labs(title='') +
             labs(x= bquote(italic(n)[.(cc_vec[kk])]), y= ylabs[kk]) +
             theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c("firebrick1", "steelblue1", "seagreen1", "darkorchid1")) +
             scale_linetype_manual(name = " ", 
                                   values = c(rep(1, 4))) +
             theme(legend.text=element_text(size=18)) +
             theme(legend.key.size = unit(1, "cm")) +
             ylim(lbs[kk],ubs[kk]) +
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))
  }
  
  ## tweak the plot axes and labels for different rows and columns
  if (kk %in% c(7, 9, 2, 4)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + xlab('') + theme(axis.text.x = element_blank(),
                                                      axis.ticks.x = element_blank()))
  }
  
  if (kk %in% c(9, 4, 14)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + ylab('') + theme(axis.text.y = element_blank(),
                                                      axis.ticks.y = element_blank()))
  }
  
  if (kk %in% c(7, 9)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + labs(title=bquote('Analysis '*.(cc_vec[kk]))) +
             theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                             margin = margin(t = 0, 0, 5, 0))))
  }
  
  if (kk %in% c(9, 4, 14)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + scale_x_continuous(breaks = seq(200, 600, 100)))
  }
  
  if (kk %in% c(7, 9)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + scale_y_continuous(breaks = seq(0, 0.1, 0.02),
                                                        limits = c(0, 0.1)))
  }
  
  if (kk %in% c(12, 14)){
    assign(paste0("plot", kk),
           get(paste0("plot", kk)) + scale_y_continuous(breaks = seq(0, 0.3, 0.1),
                                                        labels = c("0.00", "0.10", "0.20", "0.30"),
                                                        limits = c(0, 0.35)))
  }
  
}

## create a 2x3 grid (right column is left blank)
figp.row1 <- plot_grid(plot7 + theme(plot.margin=unit(c(0.1, 0.15, -0.6, 0.1),"cm")), 
                       plot9 + theme(plot.margin=unit(c(0.1, 0.15, -0.6, -0.6),"cm")),
                       NULL,
                       rel_widths = c(1.255, 1, 1), nrow = 1)
figp.row2 <- plot_grid(plot2 + theme(plot.margin=unit(c(-0.6, 0.15, -0.6, 0.1),"cm")), 
                       plot4 + theme(plot.margin=unit(c(-0.6, 0.15, -0.6, -0.6),"cm")),
                       NULL,
                       rel_widths = c(1.255, 1, 1), nrow = 1)
figp.row3 <- plot_grid(plot12 + theme(plot.margin=unit(c(-0.6, 0.15, 0.1, 0.1),"cm")), 
                       plot14 + theme(plot.margin=unit(c(-0.6, 0.15, 0.1, -0.6),"cm")),
                       NULL,
                       rel_widths = c(1.255, 1, 1), nrow = 1)

## need a third legend (for colour) that does not include design 5 (no stopping for failure)
kk <- 6
assign(paste0("plot", kk, ".leg3"),
       ggplot(subset(get(paste0("df", kk)),
                     get(paste0("df", kk))$design %in% c("D1", "D2", "D3", "D4")), aes(x=n)) + theme_bw() +
         geom_line(aes(y = y, color=as.factor(design), linetype = as.factor(design), alpha = as.factor(design)), 
                   size = 1) +
         labs(title='') +
         labs(x= bquote(italic(n)[.(cc_vec[kk])]), y= ylabs[kk]) +
         theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="bottom") +
         scale_color_manual(name = " ", 
                            values = c("firebrick1", "steelblue1", "seagreen1", "darkorchid1", "orange1",
                                       "firebrick4", "steelblue4", "seagreen4", "darkorchid4", "orange4")) +
         scale_linetype_manual(name = " ", 
                               values = c(rep(1, 5), rep(2, 5))) +
         scale_alpha_manual(name = " ", 
                            values = c(rep(0.7, 5), rep(0.9, 5))) +
         theme(legend.text=element_text(size=18)) +
         theme(legend.key.size = unit(1, "cm")) +
         ylim(lbs[kk],ubs[kk]) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))

plot6.leg3

## add legends and output figure to PDF file
figp2 <- plot_grid(figp.row1, 
                   figp.row2, 
                   figp.row3, 
                   ggpubr::get_legend(plot6.leg3),
                   ggpubr::get_legend(plot6.leg2),
                   ncol = 1, rel_heights = c(1.105, 1, 1.2, 0.15, 0.15))

# output as .pdf file for the article
pdf(file = "FigFut.pdf",   # The directory you want to save the file in
    width =10.5, # The width of the plot in inches
    height = 10.5) # The height of the plot in inches

figp2

dev.off()

## the results in Table 1 are obtained using the numbers from 
## "mat3_scen1_dat.csv", ... , "mat3_scen4_dat.csv"


## confirmatory simulations used to construct Table 2 are overviewed below
## generate the binomial summaries at the recommended sample sizes for each design

## scenario 5
## binomial summaries for the conditional approach
ns <- 194
seeds <- seq(1, length(ns), 1) + 8000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.8), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen1_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 9000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.7), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen1_", ns[k], ".csv"), row.names = FALSE)
}

## scenario 2
ns <- 221
seeds <- seq(1, length(ns), 1) + 10000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.75, 0.843), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen2_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 11000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.75, 0.75), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen2_", ns[k], ".csv"), row.names = FALSE)
}

## scenario 3
ns <- 239
seeds <- seq(1, length(ns), 1) + 12000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.8), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen3_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 13000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.7), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen3_", ns[k], ".csv"), row.names = FALSE)
}

## scenario 4
ns <- 252
require(truncnorm)
seeds <- seq(1, length(ns), 1) + 14000
p1s <- sort(rtruncnorm(10000, 0.75, 0.85, 0.8, 0.03))
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummariesPred(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), cbind(0.7, p1s), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H1_scen4_", ns[k], ".csv"), row.names = FALSE)
}

## get the binomial shape parameters and summaries for the model Psi0
seeds <- seq(1, length(ns), 1) + 15000
for (k in 1:length(ns)){
  assign(paste0("summary.cond", ns[k]), 
         getBinSummaries(ns[k], c(1, 2, 3), c(2/3, 1/2, 1/2), c(0.7, 0.7), R = 10000, seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("bin_summary_H0_scen4_", ns[k], ".csv"), row.names = FALSE)
}

## the code in lines 102 to 348 can be rerun with these new binomial summaries to get 
## the confirmatory estimates for the sampling distribution of posterior summaries;
## given these sampling distribution estimates, stopping probabilities can be computed
## using the process in Lines 515 to 620 using the decision thresholds from Table 1