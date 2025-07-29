## this code file contains the code to reproduce the numerical results and
## and plots in Section 5.2 of the text. Please run "04_ex2_functions.R" first
## to ensure that the necessary packages and functions are loaded

## load the more packages
require(cowplot)
require(ggpubr)
require(mvtnorm)

## define expit and logit functions
expit <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(x){
  return(log(x) - log(1-x))
}

## get parameter values for various beta distributions under the predictive approach;
## generate these parameter values for the sample sizes that will be explored later in
## this file
ns <- c(18, 34, 22, 26, 30, 36, 38)
for (k in 1:7){
  assign(paste0("shapes.pred", ns[k]), getShapes(10000, seed = k))
  write.csv(get(paste0("shapes.pred", ns[k])), 
            paste0("shapes_pred_", ns[k], ".csv"), row.names = FALSE)
}

## use the functions in the previous file to get shape paramters under the 
## conditional approach where the 0.99-quantile is 0.025
alpha.cond <- sapply(0.025, uuu, lower = 0.01, upper = 100, total = 400, q = 0.99)
shapes.cond <- cbind(rep(0.025, 10000), rep(alpha.cond, 10000), rep(400-alpha.cond, 10000))

## beta summaries for predictive approach
seeds <- 10000*c(0, 1, 4, 5, 6, 10, 14) + 1
for (k in 1:7){
  assign(paste0("summary.pred", ns[k]), 
         getBetaSummaries(ns[k], c(1, 1.5, 2, 2.5), get(paste0("shapes.pred", ns[k])), seed = seeds[k]))
  write.csv(get(paste0("summary.pred", ns[k])), 
            paste0("beta_summary_pred_", ns[k], ".csv"), row.names = FALSE)
}

## beta summaries for conditional approach
ns <- c(18, 34, 22, 26, 30, 32, 38)
seeds <- 10000*c(2, 3, 7, 8, 9, 11, 15) + 1
for (k in 1:7){
  assign(paste0("summary.cond", ns[k]), 
         getBetaSummaries(ns[k], c(1, 1.5, 2, 2.5), get(paste0("shapes.cond")), seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("beta_summary_cond_", ns[k], ".csv"), row.names = FALSE)
}

## get the beta shape parameters and summaries for the model Psi0
alpha.condH0 <- sapply(0.03, uuu, lower = 0.01, upper = 100, total = 400, q = 0.99)
shapes.condH0 <- cbind(rep(0.03, 10000), rep(alpha.condH0, 10000), rep(400-alpha.condH0, 10000))

seeds <- 18001
ns <- 37
for (k in 1){
  assign(paste0("summary.cond.H0", ns[k]), 
         getBetaSummariesData(ns[k], c(1, 1.5, 2, 2.5), shapes.condH0, seed = seeds[k]))
  write.csv(get(paste0("summary.cond.H0", ns[k])), 
            paste0("beta_summary_cond_H0_", ns[k], ".csv"), row.names = FALSE)
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## estimate sampling distributions under predictive approach
ns <- c(18, 34, 22, 26, 30, 38)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_pred_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "pred"
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:3){
    ni.temp <- summary.temp[1, 3*i + 1]
    nf.temp <- summary.temp[1, ncol(summary.temp)]
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
      temp <- getCompPred(ni.temp, nf.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                          burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
      write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 4
  n.temp <- summary.temp[1, ncol(summary.temp)]
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
    temp <- getCompPost(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                        burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
    write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}


## repeat process for the conditional approach
ns <- c(18, 34, 22, 26, 30, 38)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_cond_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "cond"
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:3){
    ni.temp <- summary.temp[1, 3*i + 1]
    nf.temp <- summary.temp[1, ncol(summary.temp)]
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
      temp <- getCompPred(ni.temp, nf.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                          burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
      write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 4
  n.temp <- summary.temp[1, ncol(summary.temp)]
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
    temp <- getCompPost(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                        burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
    write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## repeat for H0
ns <- 38
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_cond_H0_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "cond_H0"
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:3){
    ni.temp <- summary.temp[1, 3*i + 1]
    nf.temp <- summary.temp[1, ncol(summary.temp)]
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
      temp <- getCompPredMean(ni.temp, nf.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                              burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
      write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 4
  n.temp <- summary.temp[1, ncol(summary.temp)]
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
    temp <- getCompPostMean(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                            burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
    write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## use simulation results to get joint sampling distribution estimate under Psi0 to tune decision thresholds
ns <- 38
for (k in 1:length(ns)){
  assign(paste0("jointH0", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "cond_H0", gam = gam.final[4]))
  write.csv(get(paste0("jointH0", ns[k])), paste0("joint_cond_H0_", ns[k], ".csv"), row.names = FALSE)
}

## this function tunes the success thresholds by proportionally
## adjusting gamma.start toward 1 until the type I error estimate
## given fixed failure thresholds (fails) is lower than t1E
tuneDTsPP <- function(gamma.start, fails, t1E, factor.low = 0.5, 
                      factor.up = 1, nn = 38, ll = "cond_H0"){
  
  ## we need to construct joint sampling distibution estimates during this process
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
  
  
  alpha.start <- 1 - gamma.start
  ii <- 0
  while (factor.up - factor.low > 0.005 & ii < 15){
    factor.next <- 0.5*(factor.low + factor.up)
    alpha.next <- factor.next*alpha.start
    gam.next <- 1 - alpha.next
    
    samp.temp <- jointSampDat(nn, c(1, 1.5, 2, 2.5), ll, gam = gam.next[4])
    
    stop.temp <- stopData(samp.temp, gam.next, c(0.1, 0.2, 0.3))
    
    t1E.temp <- tail(stop.temp, 1)
    
    ## decrease factor until t1E is sufficiently low
    if (t1E.temp > t1E){
      factor.up <- factor.next
    } else {
      factor.low <- factor.next
    }
    ii <- ii + 1
  }
  return(1 - factor.low*alpha.start)
}

## obtain decision thresholds for this example
gam.final <- tuneDTsPP(c(0.999, 0.995, 0.99, 0.93), 
                       fails = c(0.1, 0.2, 0.3), t1E = 0.1)
write.csv(gam.final, "gam_final.csv", row.names = FALSE)

## given the selected final success threshold, estimate the joint sampling 
## distributions of posterior probabilities under the predictive approach
ns <- c(18, 22, 26, 30, 34, 38)
for (k in 1:length(ns)){
  assign(paste0("jointPred", ns[k]),
         jointSamp(ns[k], c(1, 1.5, 2, 2.5), ll = "pred", gam = gam.final[4]))
  write.csv(get(paste0("jointPred", ns[k])), paste0("joint_pred_", ns[k], ".csv"), row.names = FALSE)
}

## repeat estimation process under the conditional approach
for (k in 1:length(ns)){
  assign(paste0("jointCond", ns[k]),
         jointSamp(ns[k], c(1, 1.5, 2, 2.5), ll = "cond", gam = gam.final[4]))
  write.csv(get(paste0("jointCond", ns[k])), paste0("joint_cond_", ns[k], ".csv"), row.names = FALSE)
}

## get linear approximations under the predictive and conditional approaches
lines.pred <- getLinesPred(jointPred18, jointPred38, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                       38*c(1, 1, 1.5, 1.5, 2, 2, 2.5), J = 10)

lines.cond <- getLines(jointCond18, jointCond38, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                       38*c(1, 1, 1.5, 1.5, 2, 2, 2.5))

## get the matrices of stopping probabilities based on the linear approximations
xi.test <- c(0.1, 0.2, 0.3)
stop.pred <- stopMat(lines.pred, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                     18, 38, 1, gam.final, xi.test)
write.csv(stop.pred, "stop_mat_lin_pred_sec5.csv", row.names = FALSE)

stop.cond <- stopMat(lines.cond, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                     18, 38, 1, gam.final, xi.test)
write.csv(stop.cond, "stop_mat_lin_cond_sec5.csv", row.names = FALSE)

## now get the stopping probabilities based on independent estimates of the
## sampling distribution
stop.pred.dat <- NULL
for (k in 1:length(ns)){
  assign(paste0("stopVec", ns[k]),
         stopData(get(paste0("jointPred", ns[k])), gam.final, xi.test))
  stop.pred.dat <- rbind(stop.pred.dat, get(paste0("stopVec", ns[k])))
}
write.csv(stop.pred.dat, "stop_mat_dat_pred_sec5.csv", row.names = FALSE)

stop.cond.dat <- NULL
for (k in 1:length(ns)){
  assign(paste0("stopVec", ns[k]),
         stopData(get(paste0("jointCond", ns[k])), gam.final, xi.test))
  stop.cond.dat <- rbind(stop.cond.dat, get(paste0("stopVec", ns[k])))
}
write.csv(stop.cond.dat, "stop_mat_dat_cond_sec5.csv", row.names = FALSE)


## get the bootstrap confidence intervals for both approaches
## define the number and size of the bootstrap samples
MM <- 2000
m <- 10000

registerDoSNOW(cl)
pb <- txtProgressBar(max = MM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

samp1 <- jointPred18
samp2 <- jointPred38

## implement bootstrap sampling for the predictive approach
boot.pred <- foreach(k=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sort(sample(1:m, m, replace = TRUE)),]
                      samp2.temp <- samp2[sort(sample(1:m, m, replace = TRUE)),]
                      
                      lines.temp <- getLinesPred(samp1.temp, samp2.temp, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                                                 38*c(1, 1, 1.5, 1.5, 2, 2, 2.5), J = 10)
                      
                      stopMat(lines.temp, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                              18, 38, 1, gam.final, xi.test, boot = TRUE)
                    }

n.cols <- ncol(boot.pred)/7
for (k in 1:7){
  write.csv(boot.pred[, ((k-1)*n.cols + 1):(k*n.cols)], 
            paste0("CIs_pred_", k, ".csv"), row.names = FALSE)
}

## repeat bootstrap sampling for conditional approach
samp1 <- jointCond18
samp2 <- jointCond38

boot.cond <- foreach(k=1:MM, .combine=rbind,
                     .options.snow=opts) %dopar% {
                       
                       samp1.temp <- samp1[sort(sample(1:m, m, replace = TRUE)),]
                       samp2.temp <- samp2[sort(sample(1:m, m, replace = TRUE)),]
                       
                       lines.temp <- getLines(samp1.temp, samp2.temp, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                                              38*c(1, 1, 1.5, 1.5, 2, 2, 2.5))
                       
                       stopMat(lines.temp, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                               18, 38, 1, gam.final, xi.test, boot = TRUE)
                     }

n.cols <- ncol(boot.cond)/7
for (k in 1:7){
  write.csv(boot.cond[, ((k-1)*n.cols + 1):(k*n.cols)], 
            paste0("CIs_cond_", k, ".csv"), row.names = FALSE)
}

## get the bootstrap CIs for the stopping probabilities (predictive)
pred.u95 <- apply(boot.pred, 2, quantile, probs = 0.975)
pred.l95 <- apply(boot.pred, 2, quantile, probs = 0.025)

## get the bootstrap CIs for the stopping probabilities (conditional)
cond.u95 <- apply(boot.cond, 2, quantile, probs = 0.975)
cond.l95 <- apply(boot.cond, 2, quantile, probs = 0.025)

## now construct the Figure for Section 5.2

## get colour palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

c_vec <- c(1, 1, 1.5, 1.5, 2, 2, 2.5)
j <- 1
## start with the components of the data frames based on our linear approximations
assign(paste0("df", j, ".pred"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.pred[, j], stop.pred.dat[,j], pred.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     pred.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))

assign(paste0("df", j, ".cond"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.cond[, j], stop.cond.dat[,j], cond.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     cond.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))

assign(paste0("df", j),  rbind(get(paste0("df", j, ".pred")), get(paste0("df", j, ".cond"))))

## get first plot
j <- 1
assign(paste0("plot", j), 
       ggplot(get(paste0("df", j)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                  "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                  "Estimated (Conditional)", "Simulated (Conditional)", 
                                                  "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                            values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                     "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                     "Estimated (Conditional)", "Simulated (Conditional)", 
                                                     "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                               values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(n)[1]), y= bquote('Stopping Probability (Success)')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="none") +
         scale_y_continuous(limits = c(0.1, 0.9), breaks=seq(0.1, 0.9, by = 0.2)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j)))
plot1

## get next plot
j <- 3

assign(paste0("df", j, ".pred"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.pred[, j], stop.pred.dat[,j], pred.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     pred.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))

assign(paste0("df", j, ".cond"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.cond[, j], stop.cond.dat[,j], cond.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     cond.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))

assign(paste0("df", j),  rbind(get(paste0("df", j, ".pred")), get(paste0("df", j, ".cond"))))
df3$n <- c_vec[j]*df3$n
assign(paste0("plot", j), 
       ggplot(get(paste0("df", j)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                  "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                  "Estimated (Conditional)", "Simulated (Conditional)", 
                                                  "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                            values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                     "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                     "Estimated (Conditional)", "Simulated (Conditional)", 
                                                     "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                               values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(n)[2]), y= bquote('')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="none") +
         scale_y_continuous(limits = c(0.1, 0.9), breaks=seq(0.1, 0.9, by = 0.2)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2 + 0.5)))
plot3

## get next plot
j <- 5

assign(paste0("df", j, ".pred"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.pred[, j], stop.pred.dat[,j], pred.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     pred.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))

assign(paste0("df", j, ".cond"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.cond[, j], stop.cond.dat[,j], cond.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     cond.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))

assign(paste0("df", j),  rbind(get(paste0("df", j, ".pred")), get(paste0("df", j, ".cond"))))
df5$n <- c_vec[j]*df5$n
assign(paste0("plot", j), 
       ggplot(get(paste0("df", j)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                  "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                  "Estimated (Conditional)", "Simulated (Conditional)", 
                                                  "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                            values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                     "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                     "Estimated (Conditional)", "Simulated (Conditional)", 
                                                     "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                               values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(n)[3]), y= bquote('')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="none") +
         scale_y_continuous(limits = c(0.1, 0.9), breaks=seq(0.1, 0.9, by = 0.2)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2 + 0.5)))
plot5

# get next plot
j <- 7

assign(paste0("df", j, ".pred"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.pred[, j], stop.pred.dat[,j], pred.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     pred.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))

assign(paste0("df", j, ".cond"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.cond[, j], stop.cond.dat[,j], cond.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     cond.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))

assign(paste0("df", j),  rbind(get(paste0("df", j, ".pred")), get(paste0("df", j, ".cond"))))
df7$n <- c_vec[j]*df7$n
assign(paste0("plot", j), 
       ggplot(get(paste0("df", j)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                  "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                  "Estimated (Conditional)", "Simulated (Conditional)", 
                                                  "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                            values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Predictive)", "Simulated (Predictive)", 
                                                     "95% Bootstrap CI (Predictive)", "95% Bootstrap CI (Predictive)",
                                                     "Estimated (Conditional)", "Simulated (Conditional)", 
                                                     "95% Bootstrap CI (Conditional)", "95% Bootstrap CI (Conditional)"),
                               values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(n)[4]), y= bquote('')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="none") +
         scale_y_continuous(limits = c(0.1, 0.9), breaks=seq(0.1, 0.9, by = 0.2)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2 + 0.5)))
plot7

## get simplified legend for plotting
df11 <- subset(df1, !(df1$type %in% c(4, 8)))

j <- 1
assign(paste0("plot", j, "1"), 
       ggplot(get(paste0("df", j, 1)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Predictive)  ", "Simulated (Predictive)  ", 
                                                  "95% Bootstrap CI (Predictive)",
                                                  "Estimated (Conditional)  ", "Simulated (Conditional)  ", 
                                                  "95% Bootstrap CI (Conditional)"),
                            values = cbb[c(6, 7, 6, 4, 1, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Predictive)  ", "Simulated (Predictive)  ", 
                                                     "95% Bootstrap CI (Predictive)",
                                                     "Estimated (Conditional)  ", "Simulated (Conditional)  ", 
                                                     "95% Bootstrap CI (Conditional)"),
                               values = rep(c("solid", "longdash", "dashed"), 2)) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(n)[1]), y= bquote('Stopping Probability')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="bottom",
               legend.key.size = unit(0.925, "cm")) +
         ylim(0,1) + 
         guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j)))
plot11

fig.row1 <- plot_grid(plot1 + theme(plot.margin=unit(c(0.15,0.35,0.15,0.15),"cm")) +
                        theme(legend.position="none"), 
                      plot3 + theme(plot.margin=unit(c(0.15,0.35,0.15,0),"cm")) +
                        theme(legend.position="none"),
                      plot5 + theme(plot.margin=unit(c(0.15,0.35,0.15,0),"cm")) +
                        theme(legend.position="none"),
                      plot7 + theme(plot.margin=unit(c(0.15,0.35,0.15,0),"cm")) +
                        theme(legend.position="none"),
                      rel_widths = c(0.8445, 0.8, 0.8, 0.8), nrow = 1)

fig_final <- plot_grid(fig.row1, get_legend(get(paste0("plot11"))), ncol = 1, rel_heights = c(1, .2))

# output as .pdf file for the article
pdf(file = paste0("Fig2Seq.pdf"),   # The directory you want to save the file in
    width = 1.25*12.5, # The width of the plot in inches
    height = 1.25*4.5) # The height of the plot in inches

fig_final

dev.off()

## now get the statistics for Table 2

## let's get the beta shape parameters corresponding to H0
alpha.condH0 <- sapply(0.03, uuu, lower = 0.01, upper = 100, total = 400, q = 0.99)
shapes.condH0 <- cbind(rep(0.03, 10000), rep(alpha.condH0, 10000), rep(400-alpha.condH0, 10000))

## beta summaries for predictive approach
seeds <- 120001
ns <- 36
for (k in 1){
  assign(paste0("summary.pred.H0", ns[k]), 
         getBetaSummariesData(ns[k], c(1, 1.5, 2, 2.5), shapes.condH0, seed = seeds[k]))
  write.csv(get(paste0("summary.pred.H0", ns[k])), 
            paste0("beta_summary_pred_H0_", ns[k], ".csv"), row.names = FALSE)
}

## beta summaries for conditional approach
seeds <- 130001
ns <- 32
for (k in 1){
  assign(paste0("summary.cond.H0", ns[k]), 
         getBetaSummariesData(ns[k], c(1, 1.5, 2, 2.5), shapes.condH0, seed = seeds[k]))
  write.csv(get(paste0("summary.cond.H0", ns[k])), 
            paste0("beta_summary_cond_H0_", ns[k], ".csv"), row.names = FALSE)
}

## use these beta summaries to estimate sampling distributions
registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## estimate sampling distributions under predictive approach
ns <- c(36)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_pred_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "pred"
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:3){
    ni.temp <- summary.temp[1, 3*i + 1]
    nf.temp <- summary.temp[1, ncol(summary.temp)]
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
      temp <- getCompPredMean(ni.temp, nf.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                          burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
      write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 4
  n.temp <- summary.temp[1, ncol(summary.temp)]
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
    temp <- getCompPostMean(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                        burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
    write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## repeat for H0
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_pred_H0_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "pred_H0"
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:3){
    ni.temp <- summary.temp[1, 3*i + 1]
    nf.temp <- summary.temp[1, ncol(summary.temp)]
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
      temp <- getCompPredMean(ni.temp, nf.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                              burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
      write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 4
  n.temp <- summary.temp[1, ncol(summary.temp)]
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
    temp <- getCompPostMean(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                            burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
    write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## repeat process for the conditional approach
ns <- 32
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_cond_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "cond"
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:3){
    ni.temp <- summary.temp[1, 3*i + 1]
    nf.temp <- summary.temp[1, ncol(summary.temp)]
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
      temp <- getCompPredMean(ni.temp, nf.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                              burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
      write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 4
  n.temp <- summary.temp[1, ncol(summary.temp)]
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
    temp <- getCompPostMean(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                            burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
    write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}

## repeat for H0
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_cond_H0_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "cond_H0"
  
  ## get posterior and posterior predictive probabilities
  for (i in 1:3){
    ni.temp <- summary.temp[1, 3*i + 1]
    nf.temp <- summary.temp[1, ncol(summary.temp)]
    for (j in 1:10){
      summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
      temp <- getCompPredMean(ni.temp, nf.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                              burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
      write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
    }
  }
  
  ## get posterior probabilities for the last analysis
  i <- 4
  n.temp <- summary.temp[1, ncol(summary.temp)]
  for (j in 1:10){
    summary.tempj <- summary.temp[seq((j-1)*inc + 1, j*inc, by = 1), c(3*i - 1, 3*i)]
    temp <- getCompPostMean(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                            burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
    write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, "_batch_", j, ".csv"), row.names = FALSE)
  }
}


## now get the joint sampling distribution estimates and estimated stopping probabilities

## predictive approach
ns <- 36
k <- 1
assign(paste0("jointPred", ns[k]),
       jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "pred", gam = gam.final[4]))
write.csv(get(paste0("jointPred", ns[k])), paste0("joint_pred_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopPred", ns[k]),
       stopData(get(paste0("jointPred", ns[k])), gam.final, xi.test))

## repeat for H0
assign(paste0("jointPredH0", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "pred_H0", gam = gam.final[4]))
write.csv(get(paste0("jointPredH0", ns[k])), paste0("joint_pred_H0_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopPredH0", ns[k]),
         stopData(get(paste0("jointPredH0", ns[k])), gam.final, xi.test))
  
## repeat for conditional approach
ns <- 32
assign(paste0("jointCond", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "cond", gam = gam.final[4]))
write.csv(get(paste0("jointCond", ns[k])), paste0("joint_cond_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopCond", ns[k]),
         stopData(get(paste0("jointCond", ns[k])), gam.final, xi.test))
  
assign(paste0("jointCondH0", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "cond_H0", gam = gam.final[4]))
write.csv(get(paste0("jointCondH0", ns[k])), paste0("joint_cond_H0_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopCondH0", ns[k]),
         stopData(get(paste0("jointCondH0", ns[k])), gam.final, xi.test))

tab2 <- rbind(stopPred36, stopPredH036, stopCond32, stopCondH032)
write.csv(tab2, "tab2.csv", row.names = FALSE)

## get the plots for the design prior in Appendix E
## get the analytical form of the design prior for alpha

## this is the inverse transformation from alpha back to delta
finv <- function(x){
  return(qbeta(0.99, x, 400 - x))
}

## get the pdf for the uniform design prior
funif <- function(xxx){
  200
}

## create function to compute gradient
gg <- function(x){
  grad(finv, x)
}

## get the smallest and largest possible values of alpha
r1 <- uniroot(find_alpha, lower = 0.01, upper = 10, total = 400, q = 0.99, p =0.0225)$root
r2 <- uniroot(find_alpha, lower = 0.01, upper = 10, total = 400, q = 0.99, p =0.0275)$root

## get the density function values for various alpha values
alphas.x <- seq(r1, r2, length.out = 1000)
der.vec <- sapply(alphas.x, gg)
alphas.y <- funif(finv(alphas.x))*der.vec

x_dens <- NULL
x_s <- seq(0.0225,0.0275,0.00005)
for (i in 1:length(x_s)){
  x_dens[i] <- 1/0.005
}
x_df <- data.frame(X = c(0.0225,x_s,0.0275), Density = c(0,x_dens,0))

y_df <- data.frame(Y = c(r1,alphas.x,r2), Density = c(0,alphas.y,0))

z_df <- data.frame(Z = c(400-r1,400-alphas.x,400-r2), Density = c(0,alphas.y,0))

## generate plots
plot1x <- ggplot(data=x_df, aes(x=X)) + 
  geom_polygon(aes(y=Density), col=cbb[6], fill = NA, linewidth=1, alpha=0.35) +
  labs(title="") +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30)) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0.0225, 0.0275, 0.0025)) +
  labs(x= bquote(delta['r']), y= bquote('Density')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot1y <- ggplot(data=y_df, aes(x=Y)) + 
  geom_polygon(aes(y=Density), col=cbb[4], 
               fill= NA, size=1, alpha=0.35) +
  labs(title="") +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30)) +
  theme_bw() +
  labs(x= bquote(alpha), y= bquote('Density')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot1z <- ggplot(data=z_df, aes(x=Z)) + 
  geom_polygon(aes(y=Density), col=cbb[7], 
               fill=NA, size=1, alpha=0.35) +
  labs(title="") +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30)) +
  theme_bw() +
  labs(x= bquote(beta), y= bquote('Density')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

fig1 <- plot_grid(plot1x + theme(plot.margin=unit(c(0,0.5,0,0.5),"cm")), 
                  plot1y + theme(plot.margin=unit(c(0,0.5,0,0.5),"cm")),
                  plot1z + theme(plot.margin=unit(c(0,0.5,0,0.5),"cm")), nrow = 1)

fig1

# output as .pdf file for the article
pdf(file = paste0("FigDesign.pdf"),   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches

fig1

dev.off()