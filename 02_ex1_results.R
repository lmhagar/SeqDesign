## this code file contains the code to reproduce the numerical results and
## and plots in Section 5.1 of the text. Please run "01_ex_functions.R" first
## to ensure that the necessary packages and functions are loaded

## load the more packages
require(cowplot)
require(ggpubr)
require(mvtnorm)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get matrix of beta coefficients were each row corresponds to a scenario
betas.mat <- rbind(c(log(1.25), 0, log(1.4)-log(1.25)),
                   c(log(1.3), 0, log(1.6)-log(1.3)), 
                   c(log(1.35), 0, log(1.8)-log(1.35)))
## get probabilities for censoring vector
censors.v1 <- c(0.0250, 0.0378, 0.0768, 0.1097)

## estimate sampling distributions for scenario 1 (rate ratio = 1.3)
## we estimate many sampling distributions here for plotting
j <- 1
ns <- seq(100, 500, 50)
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- getCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                             hazards = c(0.055, 0.095, 0.04, 0.0200),
                             beta = as.numeric(betas.mat[j,]),
                             probs = c(0.5, 1/3), censors = censors.v1,
                             deltaL = 1, R = 10000, jags_model = "jags_surv.txt", 
                             prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                             burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000, 
                             seed = i*100)
  write.csv(probs.temp, paste0("comp_probs_scen", j,"_", ns[i], ".csv"), row.names = FALSE)
}

## repeat estimation process for scenario 2 (rate ratio = 1.4)
## we only need to estimate two sampling distributions for our method
j <- 2
ns <- c(200, 300)
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- getCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                             hazards = c(0.055, 0.095, 0.04, 0.0200),
                             beta = as.numeric(betas.mat[j,]),
                             probs = c(0.5, 1/3), censors = censors.v1,
                             deltaL = 1, R = 10000, jags_model = "jags_surv.txt", 
                             prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                             burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000, 
                             seed = 100 + (j-1)*10000)
  write.csv(probs.temp, paste0("comp_probs_scen", j,"_", ns[i], ".csv"), row.names = FALSE)
}

## repeat estimation process for scenario 3 (rate ratio = 1.5)
## we only need to estimate two sampling distributions for our method
j <- 3
ns <- c(100, 200)
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- getCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                             hazards = c(0.055, 0.095, 0.04, 0.0200),
                             beta = as.numeric(betas.mat[j,]),
                             probs = c(0.5, 1/3), censors = censors.v1,
                             deltaL = 1, R = 10000, jags_model = "jags_surv.txt", 
                             prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                             burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000, 
                             seed = 100 + (j-1)*10000)
  write.csv(probs.temp, paste0("comp_probs_scen", j,"_", ns[i], ".csv"), row.names = FALSE)
}

## now let's get the optimal sample size recommendations

## let's obtain the decision thresholds from the alpha-spending functions

## these are the cumulative alpha spent for each spending function (Pocock and OBF)
a.poc <- 0.025*log(1 + (exp(1)-1)*2:4/4)
a.obf <- 2*(1-pnorm(qnorm(1-0.025/2)/sqrt(2:4/4)))

## Define information fractions
info_times <- c(1/2, 3/4, 1)

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

## get the bounds for the Pocock-type function
gamma.1 <- 1 - a.poc[1]
gamma.2 <- getGamma2(1 - a.poc[1], a.poc[1:2], corr_matrix[1:2, 1:2])$root
gamma.3 <- getGamma3(c(1 - a.poc[1], gamma.2), a.poc, corr_matrix)$root

gamma.poc <- round(c(gamma.1, gamma.2, gamma.3),4)
## gamma.poc <- c(0.9845, 0.9896, 0.9900); these are the values

## get bounds for OBF-like function
gamma.1 <- 1 - a.obf[1]
gamma.2 <- getGamma2(gamma.1, a.obf[1:2], corr_matrix[1:2, 1:2])$root
gamma.3 <- getGamma3(c(gamma.1, gamma.2), a.obf, corr_matrix)$root

gamma.obf <- round(c(gamma.1, gamma.2, gamma.3),4)
## gamma.obf <- c(0.9985, 0.9908, 0.9780); these are the values

## define criterion for power
target.pwr <- 0.8

## construct the linear approximations and matrix of stopping probabilities
## for scenario 1
probs.200 <- read.csv(paste0("comp_probs_scen1_", 200, ".csv"))
probs.400 <- read.csv(paste0("comp_probs_scen1_", 400, ".csv"))
lin.poc1 <- stop_mat(m1 = probs.200, m2 = probs.400, c(200, 300, 400), 
                         c(400, 600, 800), 100, 500, by = 1, gamma.poc)

lin.obf1 <- stop_mat(m1 = probs.200, m2 = probs.400, c(200, 300, 400), 
                         c(400, 600, 800), 100, 500, by = 1, gamma.obf)

## get index of optimal sample size in final analysis
ind.poc1 <- which.min(split(lin.poc1, rep(1:3, each = length(lin.poc1)/3))[[3]]< target.pwr)
n.poc1 <- seq(100, 500, 1)[ind.poc1]

ind.obf1 <- which.min(split(lin.obf1, rep(1:3, each = length(lin.obf1)/3))[[3]]< target.pwr)
n.obf1 <- seq(100, 500, 1)[ind.obf1]

## construct the linear approximations and matrix of stopping probabilities
## for scenario 2
probs.200 <- read.csv(paste0("comp_probs_scen2_", 200, ".csv"))
probs.300 <- read.csv(paste0("comp_probs_scen2_", 300, ".csv"))
lin.poc2 <- stop_mat(m1 = probs.200, m2 = probs.300, c(200, 300, 400), 
                     c(300, 450, 600), 100, 500, by = 1, gamma.poc)

lin.obf2 <- stop_mat(m1 = probs.200, m2 = probs.300, c(200, 300, 400), 
                     c(300, 450, 600), 100, 500, by = 1, gamma.obf)

## get index of optimal sample size in final analysis
ind.poc2 <- which.min(split(lin.poc2, rep(1:3, each = length(lin.poc2)/3))[[3]]< target.pwr)
n.poc2 <- seq(100, 500, 1)[ind.poc2]

ind.obf2 <- which.min(split(lin.obf2, rep(1:3, each = length(lin.obf2)/3))[[3]]< target.pwr)
n.obf2 <- seq(100, 500, 1)[ind.obf2]

## construct the linear approximations and matrix of stopping probabilities
## for scenario 3
probs.200 <- read.csv(paste0("comp_probs_scen3_", 200, ".csv"))
probs.100 <- read.csv(paste0("comp_probs_scen3_", 100, ".csv"))
lin.poc3 <- stop_mat(m1 = probs.100, m2 = probs.200, c(100, 150, 200), 
                     c(200, 300, 400), 100, 500, by = 1, gamma.poc)

lin.obf3 <- stop_mat(m1 = probs.100, m2 = probs.200, c(100, 150, 200), 
                     c(200, 300, 400), 100, 500, by = 1, gamma.obf)

## get index of optimal sample size in final analysis
ind.poc3 <- which.min(split(lin.poc3, rep(1:3, each = length(lin.poc3)/3))[[3]]< target.pwr)
n.poc3 <- seq(100, 500, 1)[ind.poc3]

ind.obf3 <- which.min(split(lin.obf3, rep(1:3, each = length(lin.obf3)/3))[[3]]< target.pwr)
n.obf3 <- seq(100, 500, 1)[ind.obf3]

## now get the confirmatory power and type I error rate estimate for the table
## in the main text

## implement for scenario 1
j <- 1
ns <- c(n.obf1, n.poc1)
## confirm power
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000,
                                 seed = (i + 20)*100)
  write.csv(probs.temp, paste0("comp_probs_H1_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## confirm the type I error rate
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000,
                                 seed = (i + 22)*100)
  write.csv(probs.temp, paste0("comp_probs_H0_scen1_", ns[i], ".csv"), row.names = FALSE)
}

## repeat this process for scenario 2
j <- 2
ns <- c(n.obf2, n.poc2)
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000,
                                 seed = (i + 24)*100)
  write.csv(probs.temp, paste0("comp_probs_H1_scen2_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000,
                                 seed = (i + 26)*100)
  write.csv(probs.temp, paste0("comp_probs_H0_scen2_", ns[i], ".csv"), row.names = FALSE)
}

## repeat this process for scenario 3
j <- 3
ns <- c(n.obf3, n.poc3)
for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = as.numeric(betas.mat[j,]),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000,
                                 seed = (i + 28)*100)
  write.csv(probs.temp, paste0("comp_probs_H1_scen3_", ns[i], ".csv"), row.names = FALSE)
}

for (i in 1:length(ns)){
  print(ns[i])
  probs.temp <- confirmCompProbs(n = ns[i], c_vec = c(1, 1.5, 2), cutpoints = c(7, 14, 21, 28),
                                 hazards = c(0.055, 0.095, 0.04, 0.0200),
                                 beta = rep(0,3),
                                 probs = c(0.5, 1/3), censors = censors.v1,
                                 deltaL = 1, R = 10000, jags_model = "jags_surv.txt",
                                 prec_beta = rep(0.01, 3), prec_logh = rep(0.01, 4) ,
                                 burnin = 1000, nchains = 1, nthin = 4, ndraws = 8000,
                                 seed = (i + 30)*100)
  write.csv(probs.temp, paste0("comp_probs_H0_scen3_", ns[i], ".csv"), row.names = FALSE)
}

## implement the process to construct bootstrap confidence intervals for scenario 1
samp1 <- read.csv(paste0("comp_probs_scen1_", 200, ".csv"))
samp2 <- read.csv(paste0("comp_probs_scen1_", 400, ".csv"))

## define the number and size of the bootstrap samples
MM <- 2000
m <- 10000

registerDoSNOW(cl)
pb <- txtProgressBar(max = MM, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement bootstrap sampling for the Pocock bounds
boot.poc <- foreach(k=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      stop_mat(samp1.temp, samp2.temp, c(200, 300, 400), 
                               c(400, 600, 800), 100, 500, by = 1,
                               gamma.poc)
                    }

n.cols <- ncol(boot.poc)/3
write.csv(boot.poc[, 1:n.cols], paste0("CIs_poc_1.csv"), row.names = FALSE)
write.csv(boot.poc[, (n.cols + 1):(2*n.cols)], paste0("CIs_poc_2.csv"), row.names = FALSE)
write.csv(boot.poc[, (2*n.cols + 1):(3*n.cols)], paste0("CIs_poc_3.csv"), row.names = FALSE)

## repeat process for the OBF bounds
boot.obf <- foreach(k=1:MM, .combine=rbind,
                    .options.snow=opts) %dopar% {
                      
                      samp1.temp <- samp1[sample(1:m, m, replace = TRUE),]
                      samp2.temp <- samp2[sample(1:m, m, replace = TRUE),]
                      
                      stop_mat(samp1.temp, samp2.temp, c(200, 300, 400), 
                               c(400, 600, 800), 100, 500, by = 1,
                               gamma.obf)
                    }

n.cols <- ncol(boot.obf)/3
write.csv(boot.obf[, 1:n.cols], paste0("CIs_obf_1.csv"), row.names = FALSE)
write.csv(boot.obf[, (n.cols + 1):(2*n.cols)], paste0("CIs_obf_2.csv"), row.names = FALSE)
write.csv(boot.obf[, (2*n.cols + 1):(3*n.cols)], paste0("CIs_obf_3.csv"), row.names = FALSE)

## get the bootstrap CIs for the stopping probabilities (Pocock)
poc.u95 <- apply(boot.poc, 2, quantile, probs = 0.975)
poc.l95 <- apply(boot.poc, 2, quantile, probs = 0.025)

## get the bootstrap CIs for the stopping probabilities (OBF)
obf.u95 <- apply(boot.obf, 2, quantile, probs = 0.975)
obf.l95 <- apply(boot.obf, 2, quantile, probs = 0.025)

## get the CIs for the sample size recommendations

## start with the Pocock bounds
boot_npoc = apply(boot.poc[,(2*n.cols + 1):(3*n.cols)], 1, function(x, Gam){seq(100,500,1)[which.max(as.numeric(x) >= Gam)]}, 
                  Gam = target.pwr)
print(quantile(boot_npoc, c(0.025, 0.975)))
## [380, 393]
print(quantile(boot_npoc, c(0.025, 0.975)))

## repeat for the OBF bounds
boot_nobf = apply(boot.obf[,(2*n.cols + 1):(3*n.cols)], 1, function(x, Gam){seq(100,500,1)[which.max(as.numeric(x) >= Gam)]}, 
                  Gam = target.pwr)
print(quantile(boot_nobf, c(0.025, 0.975)))
## [330, 340]

## get the expected sample sizes for scenario 1

## split the linear approximations based on analysis
split.poc1 = split(lin.poc1, rep(1:3, each = length(lin.poc3)/3))
n.poc1*split.poc1[[1]][ind.poc1] + 1.5*n.poc1*(split.poc1[[2]][ind.poc1] - split.poc1[[1]][ind.poc1]) +  
  2*n.poc1*(1-split.poc1[[2]][ind.poc1]) ## 550.51

## repeat process for the OBF bounds
split.obf1 = split(lin.obf1, rep(1:3, each = length(lin.obf3)/3))
n.obf1*split.obf1[[1]][ind.obf1] + 1.5*n.obf1*(split.obf1[[2]][ind.obf1] - split.obf1[[1]][ind.obf1]) +  
  2*n.obf1*(1-split.obf1[[2]][ind.obf1]) ## 556.03

## create composite plots for Figure 1

## get colour palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## use this function to quickly get confirmation probabilities for all three analyses
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

## get confirmation probabilities for all three analyses
## begin with the Pocock bounds
ns <- seq(100, 500, 50)
stop.poc1 <- NULL
for (i in 1:length(ns)){
  print(ns[i])
  stop.temp <- confirm_probs(read.csv(paste0("comp_probs_scen1_", ns[i], ".csv")), gamma.poc)
  
  stop.poc1 <- rbind(stop.poc1, stop.temp)
}

## repeat with OBF bounds
ns <- seq(100, 500, 50)
stop.obf1 <- NULL
for (i in 1:length(ns)){
  print(ns[i])
  stop.temp <- confirm_probs(read.csv(paste0("comp_probs_scen1_", ns[i], ".csv")), gamma.obf)
  
  stop.obf1 <- rbind(stop.obf1, stop.temp)
}

## create data frames for the first analysis
## start with the components of the data frames based on our linear approximations
df1.poc <- data.frame(n = c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                      prob = c(lin.poc1[1:n.cols], stop.poc1[,1], poc.l95[1:n.cols], poc.u95[1:n.cols]),
                      type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns)))))

df1.obf <- data.frame(n = c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                      prob = c(lin.obf1[1:n.cols], stop.obf1[,1], obf.l95[1:n.cols], obf.u95[1:n.cols]),
                      type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns)))))

df1 <- rbind(df1.poc, df1.obf)

## get first plot
j <- 1
assign(paste0("plot", j), 
       ggplot(get(paste0("df", j)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Pocock)", "Simulated (Pocock)", 
                                                  "95% Bootstrap CI (Pocock)", "95% Bootstrap CI (Pocock)",
                                                  "Estimated (OBF)", "Simulated (OBF)", 
                                                  "95% Bootstrap CI (OBF)", "95% Bootstrap CI (OBF)"),
                            values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Pocock)", "Simulated (Pocock)", 
                                                     "95% Bootstrap CI (Pocock)", "95% Bootstrap CI (Pocock)",
                                                     "Estimated (OBF)", "Simulated (OBF)", 
                                                     "95% Bootstrap CI (OBF)", "95% Bootstrap CI (OBF)"),
                               values = rep(c("solid", "longdash", "dashed", "dashed"), 2)) +
         labs(color  = "", linetype = "") +
         labs(x= bquote(italic(n)[1]), y= bquote('Stopping Probability')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.text=element_text(size=18)) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
         theme(legend.position="none") +
         ylim(0,1) + 
         theme(plot.title = element_text(hjust = 0.5)) +
         labs(title=paste0("Analysis ", j)))
plot1

## repeat process for the second plot
df2.poc <- data.frame(n = 1.5*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                      prob = c(lin.poc1[(n.cols + 1):(2*n.cols)], stop.poc1[,2], poc.l95[(n.cols + 1):(2*n.cols)], poc.u95[(n.cols + 1):(2*n.cols)]),
                      type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns)))))

df2.obf <- data.frame(n = 1.5*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                      prob = c(lin.obf1[(n.cols + 1):(2*n.cols)], stop.obf1[,2], obf.l95[(n.cols + 1):(2*n.cols)], obf.u95[(n.cols + 1):(2*n.cols)]),
                      type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns)))))

df2 <- rbind(df2.poc, df2.obf)

j <- 2
assign(paste0("plot", j), 
       ggplot(get(paste0("df", j)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Pocock)", "Simulated (Pocock)", 
                                                  "95% Bootstrap CI (Pocock)", "95% Bootstrap CI (Pocock)",
                                                  "Estimated (OBF)", "Simulated (OBF)", 
                                                  "95% Bootstrap CI (OBF)", "95% Bootstrap CI (OBF)"),
                            values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Pocock)", "Simulated (Pocock)", 
                                                     "95% Bootstrap CI (Pocock)", "95% Bootstrap CI (Pocock)",
                                                     "Estimated (OBF)", "Simulated (OBF)", 
                                                     "95% Bootstrap CI (OBF)", "95% Bootstrap CI (OBF)"),
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
         ylim(0,1) + 
         theme(plot.title = element_text(hjust = 0.5)) +
         theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
         scale_x_continuous(breaks=1.5*seq(100, 500, 100)) +
         labs(title=paste0("Analysis ", j)))
plot2

## repeat process for the third plot
df3.poc <- data.frame(n = 2*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                      prob = c(lin.poc1[(2*n.cols + 1):(3*n.cols)], stop.poc1[,3], poc.l95[(2*n.cols + 1):(3*n.cols)], poc.u95[(2*n.cols + 1):(3*n.cols)]),
                      type = sort(c(rep(c(1, 3, 4), each = n.cols), rep(2, length(ns)))))

df3.obf <- data.frame(n = 2*c(seq(100, 500, 1), ns, rep(seq(100, 500, 1), 2)),
                      prob = c(lin.obf1[(2*n.cols + 1):(3*n.cols)], stop.obf1[,3], obf.l95[(2*n.cols + 1):(3*n.cols)], obf.u95[(2*n.cols + 1):(3*n.cols)]),
                      type = sort(c(rep(c(5, 7, 8), each = n.cols), rep(6, length(ns)))))

df3 <- rbind(df3.poc, df3.obf)

j <- 3
assign(paste0("plot", j), 
       ggplot(get(paste0("df", j)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Pocock)", "Simulated (Pocock)", 
                                                  "95% Bootstrap CI (Pocock)", "95% Bootstrap CI (Pocock)",
                                                  "Estimated (OBF)", "Simulated (OBF)", 
                                                  "95% Bootstrap CI (OBF)", "95% Bootstrap CI (OBF)"),
                            values = cbb[c(6, 7, 6, 6, 4, 1, 4, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Pocock)", "Simulated (Pocock)", 
                                                     "95% Bootstrap CI (Pocock)", "95% Bootstrap CI (Pocock)",
                                                     "Estimated (OBF)", "Simulated (OBF)", 
                                                     "95% Bootstrap CI (OBF)", "95% Bootstrap CI (OBF)"),
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
         ylim(0,1) + 
         theme(plot.title = element_text(hjust = 0.5)) +
         theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
         scale_x_continuous(breaks=2*seq(100, 500, 100)) +
         labs(title=paste0("Analysis ", j)))
plot3

## get simplified legend for plotting
df11 <- subset(df1, !(df1$type %in% c(4, 8)))

j <- 1
assign(paste0("plot", j, "1"), 
       ggplot(get(paste0("df", j, 1)), 
              aes(x=n, y=prob, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
         geom_line() +
         scale_color_manual(name = "", labels = c("Estimated (Pocock)  ", "Simulated (Pocock)  ", 
                                                  "95% Bootstrap CI (Pocock)",
                                                  "Estimated (OBF)  ", "Simulated (OBF)  ", 
                                                  "95% Bootstrap CI (OBF)"),
                            values = cbb[c(6, 7, 6, 4, 1, 4)]) +
         scale_linetype_manual(name = "", labels = c("Estimated (Pocock)  ", "Simulated (Pocock)  ", 
                                                     "95% Bootstrap CI (Pocock)",
                                                     "Estimated (OBF)  ", "Simulated (OBF)  ", 
                                                     "95% Bootstrap CI (OBF)"),
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
                       plot2 + theme(plot.margin=unit(c(0.15,0.35,0.15,0),"cm")) +
                         theme(legend.position="none"),
                       plot3 + theme(plot.margin=unit(c(0.15,0.35,0.15,0),"cm")) +
                         theme(legend.position="none"),
                       rel_widths = c(0.929, 0.8, 0.8), nrow = 1)

fig_final <- plot_grid(fig.row1, get_legend(get(paste0("plot11"))), ncol = 1, rel_heights = c(1, .2))

# output as .pdf file for the article
pdf(file = paste0("Fig1Seq.pdf"),   # The directory you want to save the file in
    width = 1.25*10.5, # The width of the plot in inches
    height = 1.25*5) # The height of the plot in inches

fig_final

dev.off()

## now obtain the numbers for Table 1

## output sample sizes
print(c(n.poc1, n.poc2, n.poc3))
print(c(n.obf1, n.obf2, n.obf3))

## get the stopping probabilities for the linear approximations
for (j in 1:3){
  assign(paste0("split.poc", j),
         split(get(paste0("lin.poc", j)), rep(1:3, each = length(get(paste0("lin.poc", j)))/3)))
  print(round(c(get(paste0("split.poc", j))[[1]][get(paste0("ind.poc", j))],
                get(paste0("split.poc", j))[[2]][get(paste0("ind.poc", j))],
                get(paste0("split.poc", j))[[3]][get(paste0("ind.poc", j))]), 4))
}

for (j in 1:3){
  assign(paste0("split.obf", j),
         split(get(paste0("lin.obf", j)), rep(1:3, each = length(get(paste0("lin.obf", j)))/3)))
  print(round(c(get(paste0("split.obf", j))[[1]][get(paste0("ind.obf", j))],
                get(paste0("split.obf", j))[[2]][get(paste0("ind.obf", j))],
                get(paste0("split.obf", j))[[3]][get(paste0("ind.obf", j))]), 4))
}

## now get stopping probabilities for confirmation simulations
for (j in 1:3){
  print(round(confirm_probs(read.csv(paste0("comp_probs_H1_scen",j, "_", 
                                            get(paste0("n.poc", j)), ".csv")), gamma.poc), 4))
  
  print(round(confirm_probs(read.csv(paste0("comp_probs_H0_scen",j, "_", 
                                            get(paste0("n.poc", j)), ".csv")), gamma.poc), 4))
  
}

for (j in 1:3){
  print(round(confirm_probs(read.csv(paste0("comp_probs_H1_scen",j, "_", 
                                            get(paste0("n.obf", j)), ".csv")), gamma.obf), 4))
  
  print(round(confirm_probs(read.csv(paste0("comp_probs_H0_scen",j, "_", 
                                            get(paste0("n.obf", j)), ".csv")), gamma.obf), 4))
  
}