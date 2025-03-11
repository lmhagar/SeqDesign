## this code file contains the code to reproduce the numerical results and
## and plots in Appendix D. Please run "03_ex2_functions.R" and "04_ex2_results_sec5.R"
## first to ensure that the necessary packages and functions are loaded

## load the more packages
require(cowplot)
require(ggpubr)
require(mvtnorm)
require(numDeriv)

## define expit and logit functions
expit <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(x){
  return(log(x) - log(1-x))
}

## estimate new sampling distributions for Appendix D
ns <- c(38)
for (k in 1){
  assign(paste0("shapes.pred", ns[k]), getShapes(10000, seed = k))
  write.csv(get(paste0("shapes.pred", ns[k])), 
            paste0("shapes_pred_v2_", ns[k], ".csv"), row.names = FALSE)
}

## use the functions in the previous file to get shape paramters under the 
## conditional approach where the 0.99-quantile is 0.025
shapes.cond <- cbind(rep(0.025, 10000), rep(alpha.cond, 10000), rep(400-alpha.cond, 10000))

## beta summaries for predictive approach
seeds <- 160000 + 1
for (k in 1){
  assign(paste0("summary.pred", ns[k]), 
         getBetaSummaries(ns[k], c(1, 1.5, 2, 2.5), get(paste0("shapes.pred", ns[k])), seed = seeds[k]))
  write.csv(get(paste0("summary.pred", ns[k])), 
            paste0("beta_summary_pred_v2_", ns[k], ".csv"), row.names = FALSE)
}

## beta summaries for conditional approach
ns <- c(34)
seeds <- 170000 + 1
for (k in 1:7){
  assign(paste0("summary.cond", ns[k]), 
         getBetaSummaries(ns[k], c(1, 1.5, 2, 2.5), get(paste0("shapes.cond")), seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("beta_summary_cond_v2_", ns[k], ".csv"), row.names = FALSE)
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## estimate sampling distributions under predictive approach
ns <- c(38)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_pred_v2_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "pred_v2"
  
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
ns <- c(34)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_cond_v2_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "cond_v2"
  
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

## consider different success thresholds
gamma.new <- c(0.999, 0.995, 0.99, 0.94)


## given the selected final success threshold, estimate the joint sampling 
## distributions of posterior probabilities under the predictive approach
ns <- c(18, 22, 26, 30, 34, 38)
for (k in 1:length(ns)){
  assign(paste0("jointPredv2", ns[k]),
         jointSamp(ns[k], c(1, 1.5, 2, 2.5), ll = "pred", gam = gamma.new[4]))
  write.csv(get(paste0("jointPredv2", ns[k])), paste0("joint_pred_v2_", ns[k], ".csv"), row.names = FALSE)
}

## repeat estimation process under the conditional approach
for (k in 1:length(ns)){
  assign(paste0("jointCondv2", ns[k]),
         jointSamp(ns[k], c(1, 1.5, 2, 2.5), ll = "cond", gam = gamma.new[4]))
  write.csv(get(paste0("jointCondv2", ns[k])), paste0("joint_cond_v2_", ns[k], ".csv"), row.names = FALSE)
}

## get linear approximations under the predictive and conditional approaches
lines.pred2 <- getLinesPred(jointPredv218, jointPredv238, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                       38*c(1, 1, 1.5, 1.5, 2, 2, 2.5), J = 10)

lines.cond2 <- getLines(jointCondv218, jointCondv238, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                       38*c(1, 1, 1.5, 1.5, 2, 2, 2.5))

## get the matrices of stopping probabilities based on the linear approximations
xi.test <- c(0.1, 0.2, 0.3)
stop.pred2 <- stopMat(lines.pred2, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                     18, 38, 1, gamma.new, xi.test)
write.csv(stop.pred2, "stop_mat_lin_pred_appd.csv", row.names = FALSE)

stop.cond2 <- stopMat(lines.cond2, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                     18, 38, 1, gamma.new, xi.test)
write.csv(stop.cond2, "stop_mat_lin_cond_appd.csv", row.names = FALSE)

## now get the stopping probabilities based on independent estimates of the
## sampling distribution
stop.pred.dat <- NULL
for (k in 1:length(ns)){
  assign(paste0("stopVecv2", ns[k]),
         stopData(get(paste0("jointPredv2", ns[k])), gamma.new, xi.test))
  stop.pred.dat <- rbind(stop.pred.dat, get(paste0("stopVecv2", ns[k])))
}
write.csv(stop.pred.dat, "stop_mat_dat_pred_appd.csv", row.names = FALSE)

stop.cond.dat <- NULL
for (k in 1:length(ns)){
  assign(paste0("stopVecv2", ns[k]),
         stopData(get(paste0("jointCondv2", ns[k])), gamma.new, xi.test))
  stop.cond.dat <- rbind(stop.cond.dat, get(paste0("stopVecv2", ns[k])))
}
write.csv(stop.cond.dat, "stop_mat_dat_cond_appd.csv", row.names = FALSE)


## get the bootstrap confidence intervals for both approaches
## define the number and size of the bootstrap samples
MM <- 2000
m <- 10000

registerDoSNOW(cl)
pb <- txtProgressBar(max = MM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

samp1 <- jointPredv218
samp2 <- jointPredv238

## implement bootstrap sampling for the predictive approach
boot.pred <- foreach(k=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sort(sample(1:m, m, replace = TRUE)),]
                      samp2.temp <- samp2[sort(sample(1:m, m, replace = TRUE)),]
                      
                      lines.temp <- getLinesPred(samp1.temp, samp2.temp, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                                                 38*c(1, 1, 1.5, 1.5, 2, 2, 2.5), J = 10)
                      
                      stopMat(lines.temp, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                              18, 38, 1, gamma.new, xi.test, boot = TRUE)
                    }

n.cols <- ncol(boot.pred)/7
for (k in 1:7){
  write.csv(boot.pred[, ((k-1)*n.cols + 1):(k*n.cols)], 
            paste0("CIs_pred_v2_", k, ".csv"), row.names = FALSE)
}

## repeat bootstrap sampling for conditional approach
samp1 <- jointCondv218
samp2 <- jointCondv238

boot.cond <- foreach(k=1:MM, .combine=rbind,
                     .options.snow=opts) %dopar% {
                       
                       samp1.temp <- samp1[sort(sample(1:m, m, replace = TRUE)),]
                       samp2.temp <- samp2[sort(sample(1:m, m, replace = TRUE)),]
                       
                       lines.temp <- getLines(samp1.temp, samp2.temp, 18*c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                                              38*c(1, 1, 1.5, 1.5, 2, 2, 2.5))
                       
                       stopMat(lines.temp, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                               18, 38, 1, gamma.new, xi.test, boot = TRUE)
                     }

n.cols <- ncol(boot.cond)/7
for (k in 1:7){
  write.csv(boot.cond[, ((k-1)*n.cols + 1):(k*n.cols)], 
            paste0("CIs_cond_v2_", k, ".csv"), row.names = FALSE)
}

## get the bootstrap CIs for the stopping probabilities (predictive)
pred.u95 <- apply(boot.pred, 2, quantile, probs = 0.975)
pred.l95 <- apply(boot.pred, 2, quantile, probs = 0.025)

## get the bootstrap CIs for the stopping probabilities (conditional)
cond.u95 <- apply(boot.cond, 2, quantile, probs = 0.975)
cond.l95 <- apply(boot.cond, 2, quantile, probs = 0.025)


## get the CIs for the sample size recommendations

## start with predictive approach
target.pwr <- 0.8
boot_npred = apply(boot.pred[,(6*n.cols + 1):(7*n.cols)], 1, function(x, Gam){seq(18,38,1)[which.max(as.numeric(x) >= Gam)]}, 
                  Gam = target.pwr)
print(quantile(boot_npred, c(0.025, 0.975)))
## [37, 39]

## repeat for conditional approach
boot_ncond = apply(boot.cond[,(6*n.cols + 1):(7*n.cols)], 1, function(x, Gam){seq(18,38,1)[which.max(as.numeric(x) >= Gam)]}, 
                   Gam = target.pwr)
print(quantile(boot_ncond, c(0.025, 0.975)))
## [33, 35]

## now construct the Figure for Section 5.2

stop.pred <- stop.pred2
stop.cond <- stop.cond2
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
         scale_y_continuous(limits = c(0, 0.9), breaks=seq(0, 0.9, by = 0.3)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j)))
plot1

## get next plot
j <- 2

assign(paste0("df", j, ".pred"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.pred[, j], stop.pred.dat[,j], pred.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     pred.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))

assign(paste0("df", j, ".cond"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.cond[, j], stop.cond.dat[,j], cond.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     cond.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))

assign(paste0("df", j),  rbind(get(paste0("df", j, ".pred")), get(paste0("df", j, ".cond"))))

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
         labs(x= bquote(italic(n)[1]), y= bquote('Stopping Probability (Failure)')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="none") +
         ylim(0,0.3) + 
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2)))
plot2

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
         scale_y_continuous(limits = c(0, 0.9), breaks=seq(0, 0.9, by = 0.3)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2 + 0.5)))
plot3

## get next plot
j <- 4

assign(paste0("df", j, ".pred"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.pred[, j], stop.pred.dat[,j], pred.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     pred.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))

assign(paste0("df", j, ".cond"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.cond[, j], stop.cond.dat[,j], cond.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     cond.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))

assign(paste0("df", j),  rbind(get(paste0("df", j, ".pred")), get(paste0("df", j, ".cond"))))
df4$n <- c_vec[j]*df4$n

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
         ylim(0,0.3) + 
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2)))
plot4

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
         scale_y_continuous(limits = c(0, 0.9), breaks=seq(0, 0.9, by = 0.3)) +
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2 + 0.5)))
plot5

## get next plot
j <- 6

assign(paste0("df", j, ".pred"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.pred[, j], stop.pred.dat[,j], pred.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     pred.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns))))))

assign(paste0("df", j, ".cond"), data.frame(n = c(seq(18, 38, 1), ns, rep(seq(18, 38, 1), 2)),
                                            prob = c(stop.cond[, j], stop.cond.dat[,j], cond.l95[((j-1)*n.cols + 1):(j*n.cols)], 
                                                     cond.u95[((j-1)*n.cols + 1):(j*n.cols)]),
                                            type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns))))))

assign(paste0("df", j),  rbind(get(paste0("df", j, ".pred")), get(paste0("df", j, ".cond"))))
df6$n <- c_vec[j]*df6$n

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
         ylim(0,0.3) + 
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j/2)))
plot6

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
         scale_y_continuous(limits = c(0, 0.9), breaks=seq(0, 0.9, by = 0.3)) +
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

fig.row2 <- plot_grid(plot2 + theme(plot.margin=unit(c(0.15,0.35,0.15,0.15),"cm")) +
                        theme(legend.position="none"), 
                      plot4 + theme(plot.margin=unit(c(0.15,0.35,0.15,0),"cm")) +
                        theme(legend.position="none"),
                      plot6 + theme(plot.margin=unit(c(0.15,0.35,0.15,0),"cm")) +
                        theme(legend.position="none"),
                      NULL,
                      rel_widths = c(0.8445, 0.8, 0.8, 0.8), nrow = 1)

fig_final <- plot_grid(fig.row1, fig.row2, get_legend(get(paste0("plot11"))), ncol = 1, rel_heights = c(1, 1, .2))

# output as .pdf file for the article
pdf(file = paste0("FigDSeq.pdf"),   # The directory you want to save the file in
    width = 1.25*12.5, # The width of the plot in inches
    height = 1.25*8) # The height of the plot in inches

fig_final

dev.off()

## now get the statistics for Table in Appendix D

## let's get the beta shape parameters corresponding to H0
alpha.condH0 <- sapply(0.03, uuu, lower = 0.01, upper = 100, total = 400, q = 0.99)
shapes.condH0 <- cbind(rep(0.03, 10000), rep(alpha.condH0, 10000), rep(400-alpha.condH0, 10000))

## beta summaries for predictive approach
seeds <- 180001
ns <- 38
for (k in 1){
  assign(paste0("summary.pred.H0", ns[k]), 
         getBetaSummariesData(ns[k], c(1, 1.5, 2, 2.5), shapes.condH0, seed = seeds[k]))
  write.csv(get(paste0("summary.pred.H0", ns[k])), 
            paste0("beta_summary_pred_H0_", ns[k], ".csv"), row.names = FALSE)
}

## beta summaries for conditional approach
seeds <- 190001
ns <- 34
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
ns <- c(38)
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_pred_v2_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "pred_v2"
  
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
ns <- 34
for (k in 1:length(ns)){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_cond_v2_", ns[k], ".csv"))
  inc <- nrow(summary.temp)/10
  ll <- "cond_v2"
  
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
ns <- 38
k <- 1
assign(paste0("jointPredv2", ns[k]),
       jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "pred_v2", gam = gamma.new[4]))
write.csv(get(paste0("jointPredv2", ns[k])), paste0("joint_pred_v2_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopPredv2", ns[k]),
       stopData(get(paste0("jointPredv2", ns[k])), gamma.new, xi.test))

## repeat for H0
assign(paste0("jointPredH0", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "pred_H0", gam = gamma.new[4]))
write.csv(get(paste0("jointPredH0", ns[k])), paste0("joint_pred_H0_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopPredH0", ns[k]),
         stopData(get(paste0("jointPredH0", ns[k])), gamma.new, xi.test))
  
## repeat for conditional approach
ns <- 34
assign(paste0("jointCondv2", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "cond_v2", gam = gamma.new[4]))
write.csv(get(paste0("jointCondv2", ns[k])), paste0("joint_cond_v2_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopCondv2", ns[k]),
         stopData(get(paste0("jointCondv2", ns[k])), gamma.new, xi.test))
  
assign(paste0("jointCondH0", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "cond_H0", gam = gamma.new[4]))
write.csv(get(paste0("jointCondH0", ns[k])), paste0("joint_cond_H0_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopCondH0", ns[k]),
         stopData(get(paste0("jointCondH0", ns[k])), gamma.new, xi.test))

tabd <- rbind(stop.pred[21,], stopPredv238, stopPredH038, stop.cond[17,], stopCondv234, stopCondH034)
write.csv(tabd, "tabd.csv", row.names = FALSE)

## get the plots for the design prior
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
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20)) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0.0225, 0.0275, 0.0025)) +
  labs(x= bquote(delta['r']), y= bquote('Density')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot1y <- ggplot(data=y_df, aes(x=Y)) + 
  geom_polygon(aes(y=Density), col=cbb[4], 
               fill= NA, size=1, alpha=0.35) +
  labs(title="") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20)) +
  theme_bw() +
  labs(x= bquote(alpha), y= bquote('Density')) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot1z <- ggplot(data=z_df, aes(x=Z)) + 
  geom_polygon(aes(y=Density), col=cbb[7], 
               fill=NA, size=1, alpha=0.35) +
  labs(title="") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20)) +
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

## now we get the sampling distribution of posterior probabilities for various sample sizes
## under the hypothesis H0
seeds <- 10000*c(20,21,22,23,24) + 1
ns <- c(20, 100, 500, 200, 1000)
for (k in 1:3){
  assign(paste0("summary.cond.H0", ns[k]), 
         getBetaSummariesData(ns[k], c(1, 1.5, 2, 2.5), shapes.condH0, seed = seeds[k]))
  write.csv(get(paste0("summary.cond.H0", ns[k])), 
            paste0("beta_summary_cond_H0_", ns[k], ".csv"), row.names = FALSE)
}

## get the sampling distribution of posterior probabilities
registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

for (k in 1:5){
  ## extract the data summary for this sample size
  summary.temp <- read.csv(paste0("beta_summary_cond_H0_", ns[k], ".csv"))
  ll <- "cond_H0"
  n.temp <- summary.temp[1, 4]
  summary.tempj <- summary.temp[, c(3*i - 1, 3*i)]
  temp <- getCompPostMean(n.temp, summary.tempj, 0.03, "jags_beta.txt", hyper = 10000, 
                          burnin = 1000, nchains = 1, nthin = 1, ndraws = 5000)
  write.csv(temp, paste0(ll, "_n_",summary.temp[1,4],"_t_", i, ".csv"), row.names = FALSE)
}

## obtain vectors for plotting
for (k in 1:5){
  assign(paste0("post", ns[k]), 1 - read.csv(paste0(ll, "_n_",ns[k],"_t_", i, ".csv"))$V1)
}

## create data frames and histograms
df20 <- data.frame(x = post20)
df200 <- data.frame(x = post200)
df1000 <- data.frame(x = post1000)

## n = 20
plot1a <-
  ggplot(data=df20, aes(x = x)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))), breaks = seq(0,1, by = 0.05),
                 col="white", 
                 fill="grey40", 
                 alpha = 0.6, size = 1) + 
  labs(x=bquote(italic('Pr')*'('*italic('H')[1]*' | '*italic('data')*')'), 
       y="Proportion") + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) + 
  theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
  labs(title=bquote(italic(n)[1] == 20)) +
  theme(plot.title = element_text(size=20,face="bold", hjust = 0.5,
                                  margin = margin(t = 0, 0, 4, 0)))
plot1a

## n = 200
plot1b <-
  ggplot(data=df200, aes(x = x)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))), breaks = seq(0,1, by = 0.05),
                 col="white", 
                 fill="grey40", 
                 alpha = 0.6, size = 1) + 
  labs(x=bquote(italic('Pr')*'('*italic('H')[1]*' | '*italic('data')*')'), 
       y="Proportion") + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) + 
  theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
  labs(title=bquote(italic(n)[1] == 200)) +
  theme(plot.title = element_text(size=20,face="bold", hjust = 0.5,
                                  margin = margin(t = 0, 0, 4, 0)))
plot1b

# n = 1000
plot1c <-
  ggplot(data=df1000, aes(x = x)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))), breaks = seq(0,1, by = 0.05),
                 col="white", 
                 fill="grey40", 
                 alpha = 0.6, size = 1) + 
  labs(x=bquote(italic('Pr')*'('*italic('H')[1]*' | '*italic('data')*')'), 
       y="Proportion") + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) + 
  theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
  labs(title=bquote(italic(n)[1] == 1000)) +
  theme(plot.title = element_text(size=20,face="bold", hjust = 0.5,
                                  margin = margin(t = 0, 0, 4, 0)))
plot1c

fig.hist <- plot_grid(plot1a + theme(plot.margin=unit(c(0.05,0.5,0.05,0.5),"cm")), 
                  plot1b + theme(plot.margin=unit(c(0.05,0.5,0.05,0.5),"cm")),
                  plot1c + theme(plot.margin=unit(c(0.05,0.5,0.05,0.5),"cm")), nrow = 1)

fig.hist

# output as .pdf file for the article
pdf(file = paste0("FigHist.pdf"),   # The directory you want to save the file in
    width = 12.5, # The width of the plot in inches
    height = 4) # The height of the plot in inches

fig.hist

dev.off()

## now we verify the appropriateness of our linearity result for smaller sample sizes
l1 <- function(x){log(1-x)-log(x)}

ns <- seq(18, 38, 4)
qs <- NULL
for (i in 1:length(ns)){
  ## extract quantiles from independent estimates of the sampling distribution
  dat <- read.csv(paste0("joint_cond_", ns[i], ".csv"))  
  qs <- cbind(qs, quantile(l1(dat[,2]), seq(0.1, 0.9, 0.1)))
}

df.q <- data.frame(l = as.numeric(qs), n = rep(seq(18, 38, 4), each = 9),
                   q = rep(seq(0.1, 0.9, 0.1), 6))

## plot quantiles from independent estimates of the sampling distribution as 
## a function of the sample size
plot.final <- ggplot(df.q, 
              aes(x=n, y=l, color = as.factor(q))) + theme_bw() +
         geom_line(linewidth =1) +
         scale_color_manual(name = bquote(italic(q)), labels = seq(0.1, 0.9, 0.1),
                            values = c("black", 2:8, "brown")) +
         labs(color = bquote(italic(q))) +
         labs(x= bquote(italic(n)[1]), y= bquote('Logit of Posterior Predictive Probability')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(legend.title=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="bottom") +
         guides(color = guide_legend(nrow = 1, byrow = TRUE)) 
plot.final

# output as .pdf file for the article
pdf(file = paste0("FigLogit.pdf"),   # The directory you want to save the file in
    width =9, # The width of the plot in inches
    height = 6) # The height of the plot in inches

plot.final

dev.off()
