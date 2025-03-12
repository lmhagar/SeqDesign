## this code file contains the code to reproduce the numerical results and
## and plots in Section 5.2 of the text. Please run "03_ex2_functions.R" first
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
ns <- c(18, 34, 22, 26, 30, 29, 38)
for (k in 1:7){
  assign(paste0("shapes.pred", ns[k]), getShapes(10000, seed = k))
  write.csv(get(paste0("shapes.pred", ns[k])), 
            paste0("shapes_pred_", ns[k], ".csv"), row.names = FALSE)
}

## use the functions in the previous file to get shape paramters under the 
## conditional approach where the 0.99-quantile is 0.025
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
ns <- c(18, 34, 22, 26, 30, 25, 38)
seeds <- 10000*c(2, 3, 7, 8, 9, 11, 15) + 1
for (k in 1:7){
  assign(paste0("summary.cond", ns[k]), 
         getBetaSummaries(ns[k], c(1, 1.5, 2, 2.5), get(paste0("shapes.cond")), seed = seeds[k]))
  write.csv(get(paste0("summary.cond", ns[k])), 
            paste0("beta_summary_cond_", ns[k], ".csv"), row.names = FALSE)
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

## now get the decision thresholds based on the OBF alpha-spending function
a.obf <- 2*(1-pnorm(qnorm(1-0.1/2)/sqrt(2:5/5)))

## Define information fractions
info_times <- 2:5/5

## Define correlation matrix
corr_matrix <- outer(info_times, info_times, FUN = function(x, y) sqrt(pmin(x, y) / pmax(x, y)))

## function to get get the critical value for second analysis
getGamma2 <- function(gamma.prev, alpha.cum, cor.mat){
  error_fn <- function(gamma.prev, cor.mat, gamma.new, target.err){
    pmvnorm(lower = c(-Inf, qnorm(gamma.new)), upper = c(qnorm(gamma.prev[1]), Inf), mean = rep(0, 2), corr = cor.mat) - target.err
  }
  return(uniroot(error_fn, lower = 0.5, upper = 0.999999, gamma.prev = gamma.prev, cor.mat = cor.mat, target.err = alpha.cum[2] - alpha.cum[1]))
}

## function to get the critical value for third analysis
getGamma3 <- function(gamma.prev, alpha.cum, cor.mat){
  error_fn <- function(gamma.prev, cor.mat, gamma.new, target.err){
    pmvnorm(lower = c(-Inf, -Inf, qnorm(gamma.new)), upper = c(qnorm(gamma.prev[1]), qnorm(gamma.prev[2]), Inf), mean = rep(0, 3), corr = cor.mat) - target.err
  }
  return(uniroot(error_fn, lower = 0.5, upper = 0.999999, gamma.prev = gamma.prev, cor.mat = cor.mat, target.err = alpha.cum[3] - alpha.cum[2]))
}

## function to get the critical value for fourth analysis
getGamma4 <- function(gamma.prev, alpha.cum, cor.mat){
  error_fn <- function(gamma.prev, cor.mat, gamma.new, target.err){
    pmvnorm(lower = c(-Inf, -Inf, -Inf, qnorm(gamma.new)), upper = c(qnorm(gamma.prev[1]), qnorm(gamma.prev[2]), qnorm(gamma.prev[3]), Inf), mean = rep(0, 4), corr = cor.mat) - target.err
  }
  return(uniroot(error_fn, lower = 0.5, upper = 0.999999, gamma.prev = gamma.prev, cor.mat = cor.mat, target.err = alpha.cum[4] - alpha.cum[3]))
}

## get the bounds for the Pocock-type function
gamma.1 <- 1 - a.obf[1]
gamma.2 <- getGamma2(gamma.1, a.obf[1:2], corr_matrix[1:2, 1:2])$root
gamma.3 <- getGamma3(c(gamma.1, gamma.2), a.obf, corr_matrix[1:3, 1:3])$root
gamma.4 <- getGamma4(c(gamma.1, gamma.2, gamma.3), a.obf, corr_matrix)$root

gamma.obf <- round(c(gamma.1, gamma.2, gamma.3, gamma.4),4)
write.csv(gamma.obf, "gamma_obf.csv", row.names = FALSE)

## given the selected final success threshold, estimate the joint sampling 
## distributions of posterior probabilities under the predictive approach
ns <- c(18, 22, 26, 30, 34, 38)
for (k in 1:length(ns)){
  assign(paste0("jointPred", ns[k]),
         jointSamp(ns[k], c(1, 1.5, 2, 2.5), ll = "pred", gam = gamma.obf[4]))
  write.csv(get(paste0("jointPred", ns[k])), paste0("joint_pred_", ns[k], ".csv"), row.names = FALSE)
}

## repeat estimation process under the conditional approach
for (k in 1:length(ns)){
  assign(paste0("jointCond", ns[k]),
         jointSamp(ns[k], c(1, 1.5, 2, 2.5), ll = "cond", gam = gamma.obf[4]))
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
                     18, 38, 1, gamma.obf, xi.test)
write.csv(stop.pred, "stop_mat_lin_pred_sec5.csv", row.names = FALSE)

stop.cond <- stopMat(lines.cond, c(1, 1, 1.5, 1.5, 2, 2, 2.5), 
                     18, 38, 1, gamma.obf, xi.test)
write.csv(stop.cond, "stop_mat_lin_cond_sec5.csv", row.names = FALSE)

## now get the stopping probabilities based on independent estimates of the
## sampling distribution
stop.pred.dat <- NULL
for (k in 1:length(ns)){
  assign(paste0("stopVec", ns[k]),
         stopData(get(paste0("jointPred", ns[k])), gamma.obf, xi.test))
  stop.pred.dat <- rbind(stop.pred.dat, get(paste0("stopVec", ns[k])))
}
write.csv(stop.pred.dat, "stop_mat_dat_pred_sec5.csv", row.names = FALSE)

stop.cond.dat <- NULL
for (k in 1:length(ns)){
  assign(paste0("stopVec", ns[k]),
         stopData(get(paste0("jointCond", ns[k])), gamma.obf, xi.test))
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
                              18, 38, 1, gamma.obf, xi.test, boot = TRUE)
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
                               18, 38, 1, gamma.obf, xi.test, boot = TRUE)
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


## get the CIs for the sample size recommendations

## start with predictive approach
target.pwr <- 0.8
boot_npred = apply(boot.pred[,(6*n.cols + 1):(7*n.cols)], 1, function(x, Gam){seq(18,38,1)[which.max(as.numeric(x) >= Gam)]}, 
                  Gam = target.pwr)
print(quantile(boot_npred, c(0.025, 0.975)))
## [28, 30]

## repeat for conditional approach
boot_ncond = apply(boot.cond[,(6*n.cols + 1):(7*n.cols)], 1, function(x, Gam){seq(18,38,1)[which.max(as.numeric(x) >= Gam)]}, 
                   Gam = target.pwr)
print(quantile(boot_ncond, c(0.025, 0.975)))
## [25, 26]

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
ns <- 29
for (k in 1){
  assign(paste0("summary.pred.H0", ns[k]), 
         getBetaSummariesData(ns[k], c(1, 1.5, 2, 2.5), shapes.condH0, seed = seeds[k]))
  write.csv(get(paste0("summary.pred.H0", ns[k])), 
            paste0("beta_summary_pred_H0_", ns[k], ".csv"), row.names = FALSE)
}

## beta summaries for conditional approach
seeds <- 130001
ns <- 25
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
ns <- c(29)
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
ns <- 25
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
ns <- 29
k <- 1
assign(paste0("jointPred", ns[k]),
       jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "pred", gam = gamma.obf[4]))
write.csv(get(paste0("jointPred", ns[k])), paste0("joint_pred_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopPred", ns[k]),
       stopData(get(paste0("jointPred", ns[k])), gamma.obf, xi.test))

## repeat for H0
assign(paste0("jointPredH0", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "pred_H0", gam = gamma.obf[4]))
write.csv(get(paste0("jointPredH0", ns[k])), paste0("joint_pred_H0_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopPredH0", ns[k]),
         stopData(get(paste0("jointPredH0", ns[k])), gamma.obf, xi.test))
  
## repeat for conditional approach
ns <- 25
assign(paste0("jointCond", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "cond", gam = gamma.obf[4]))
write.csv(get(paste0("jointCond", ns[k])), paste0("joint_cond_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopCond", ns[k]),
         stopData(get(paste0("jointCond", ns[k])), gamma.obf, xi.test))
  
assign(paste0("jointCondH0", ns[k]),
         jointSampDat(ns[k], c(1, 1.5, 2, 2.5), ll = "cond_H0", gam = gamma.obf[4]))
write.csv(get(paste0("jointCondH0", ns[k])), paste0("joint_cond_H0_", ns[k], ".csv"), row.names = FALSE)
  
assign(paste0("stopCondH0", ns[k]),
         stopData(get(paste0("jointCondH0", ns[k])), gamma.obf, xi.test))

tab2 <- rbind(stop.pred[12,], stopPred29, stopPredH029, stop.cond[8,], stopCond25, stopCondH025)
write.csv(tab2, "tab2.csv", row.names = FALSE)