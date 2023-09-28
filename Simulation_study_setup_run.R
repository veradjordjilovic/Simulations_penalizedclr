require(tidyverse)
require(penalizedclr)
data_generation <- function(p1 = 100, p2 = 100, p1_r = 10, p2_r = 10,
                            beta1 = 0.5, beta2 = 0.5, cor.coef = 0, n = 50){
  # corr.coef not implemented yet
  # number of covariates on two blocks
  p <- c(p1, p2)
  # number of covariates with signal in each block
  m_nz <- c(p1_r, p2_r)
  # number of different strata (case-control pairs)
  K <- n
  
  # number of cases and controls in each stratum (not necessarily 1:1 matching,
  # other designs are also allowed)
  n_cases <- 1
  n_ctrl <- 1
  
  # generating covariates
  X = cbind(matrix(rnorm(p[1] * K * (n_cases + n_ctrl), 0, 1), ncol = p[1]),
            matrix(rnorm(p[2] * K * (n_cases + n_ctrl), 0, 1), ncol = p[2]))
  
  # coefficients
  beta1_v <- rep(beta1, m_nz[1]) * sample(c(1, -1), m_nz[1], replace = T)
  beta2_v <- rep(beta2, m_nz[2]) * sample(c(1, -1), m_nz[2], replace = T)
  beta <- as.matrix(c(beta1_v, rep(0, p[1] - m_nz[1]),
                      beta2_v, rep(0, p[2] - m_nz[2])), ncol = 1)
  
  
  
 # beta_stratum <- rep(rnorm(K, 0, 2), each = n_cases+n_ctrl)
  
 
  # stratum membership
  stratum <- rep(1:K, each= n_cases+n_ctrl)
  
  # linear predictor
  lin_pred <- X %*% beta
  
  prob_case <- exp(lin_pred) / (1 + exp(lin_pred))
  
  
  # generate the response
  
  Y <- rep(0, length(stratum))
  
  data_sim <- as_tibble(data.frame(stratum = stratum,
                                   probability = prob_case,
                                   obs_id = 1 : length(stratum)))
  data_sim_cases <- data_sim %>%
    group_by(stratum)%>%
    sample_n(n_cases, weight = probability)
  
  Y[data_sim_cases$obs_id] <- 1
  
  return(list(Y = Y, X = X, stratum = stratum, beta = beta))
}

simulation_study <- function(p1 = 100, p2 = 100, p1_r = 10, p2_r = 10,
                             beta1 = 0.5, beta2 = 0.5, cor.coef = 0,
                             n = 50,
                             B = 50){
  # Input: 
  # p1... size of the first block of predictors
  # p2... size of the second block of predictors
  # p1_r... number of signal variables in the first block
  # p2_r... number of signal variables in the second block
  # beta1... coefficient of the signal variables in the first block
  # beta2... coefficient of the signal variables in the second block
  # corr.coef... correlation coefficient between correlated variables
  # reps... number of Monte Carlo runs
  # n.. number of strata
  # B... number of subsamples for stability selection
  # Output:
  # a list with two matrices: the first one giving selection status
  # in each run for each variable for IPF-lasso and the second one
  # giving selection probabilities from stable.clr.g for each run/variable

data_gen <- data_generation(p1 = p1, p2 = p2, p1_r = p1_r,
                            p2_r = p2_r, beta1 = beta1, beta2 = beta2,
                            cor.coef = cor.coef, n = n)

response <- data_gen$Y
penalized <- data_gen$X
stratum <- data_gen$stratum

# find data driven vector of penalty factors
pf  <- default.pf(response = response, 
                  stratum = stratum,
                  penalized = penalized, 
                  p = c(p1, p2),
                  type.step1 = "comb")

if(length(pf$pf) == 2) {pf_v <- c(1, pf$pf[2]/pf$pf[1])} else {pf_v <- c(1,1)
pf_v[pf$exc] = 4}

lambda <- find.default.lambda(response = response, 
                              stratum = stratum,
                              penalized = penalized, 
                              p = c(p1, p2),
                              pf.list = list(pf_v))

lambda.vector <- unlist(lambda) * pf_v

# run stability selection with penalized conditional logistic regression
stable_res <- stable.clr.g(response = response,
                           stratum = stratum, 
                           p = c(p1, p2),  
                           standardize = TRUE,
                           penalized = penalized,  
                           lambda.list = list(lambda.vector))
return(list(Pistab = stable_res$Pistab, lambda.vector = lambda.vector))
}

settings <- list(setting1 = list(p1 = 50, p2 = 50, p1_r = 10, p2_r = 10,
                              beta1 = 4, beta2 = 4, n = 200),
                 setting2 = list(p1 = 50, p2 = 50, p1_r = 3, p2_r = 17,
                                 beta1 = 4, beta2 = 4),
                 setting3 = list(p1 = 50, p2 = 50, p1_r = 20, p2_r = 0,
                                 beta1 = 4, beta2 = 0),
                 setting4 = list(p1 = 20, p2 = 80, p1_r = 10, p2_r = 10,
                                 beta1 = 4, beta2 = 1),
                 setting5 = list(p1 = 20, p2 = 80, p1_r = 15, p2_r = 5,
                                 beta1 = 4, beta2 = 4),
                 setting6 = list(p1 = 20, p2 = 80, p1_r = 5, p2_r = 15,
                                 beta1 = 4, beta2 = 4))

true_signal <- lapply(settings, function(x) c(rep(1, x$p1_r,),
                                              rep(0, x$p1 - x$p1_r),
                                              rep(1, x$p2_r,),
                                              rep(0, x$p2 - x$p2_r)))
# run simulations
set.seed(123)
nruns <- 100
temp <- lapply(settings, function(x) 
  replicate(nruns, simulation_study(p1 = x$p1, p2 = x$p2,
                                    p1_r = x$p1_r,
                                    p2_r = x$p2_r, beta1 = x$beta1,
                                    beta2 = x$beta2, n = 200)$Pistab))



threshold <- seq(from = 0.55, to = 0.9, by = 0.05)
final_results<- list()
res <- list()
for (j in 1:length(threshold)){
for (i in 1:length(settings)){
  true_status <- true_signal[[i]]
  est_status <- as.matrix(temp[[i]])
  est_status <- est_status >= threshold[j]
  res[[i]] <- t(apply(est_status, 2, function(x) c(TP = sum((x == 1) & true_status == 1), 
                                    TN = sum((x == 0) & true_status == 0),
                                    FP = sum((x == 1) & true_status == 0),
                                    FN = sum((x == 0) & true_status == 1),
                                    power = sum((x == 1) & true_status == 1)/max(1,
                                            (sum((x == 1) & true_status == 1) + 
                                             sum((x == 0) & true_status == 1))),
                                    FDR =  sum((x == 1) & true_status == 0)/
                                           max((sum((x == 1) & true_status == 0)+
                                           sum((x == 1) & true_status == 1)), 1))))
}

final_results[[j]] <- lapply(res, colMeans)}

power <- matrix(0, ncol = length(settings), nrow = length(threshold))
for (i in 1:length(threshold)){
  temporary_list = final_results[[i]]
  power[i, ] <-  sapply(temporary_list, function(x) x["power"])}

FDR <- matrix(0, ncol = length(settings), nrow = length(threshold))
for (i in 1:length(threshold)){
  temporary_list = final_results[[i]]
  FDR[i, ] <-  sapply(temporary_list, function(x) x["FDR"])}


### Create figure

setwd("~/Library/CloudStorage/Dropbox/WorkingProjects/MultiOmics/Simulations")
pdf("Power_fdr_sim_study.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(x = threshold, y = power[, 1], type = "b", col = "lightsalmon",
     xlab = "Threshold", ylab = "Power", 
     xlim = c(0.55, 0.9), 
     ylim = c(0,1), pch = 15)

points(x = threshold, y = power[, 2], type = "b", col = "darkolivegreen4", pch = 16)

points(x = threshold, y = power[, 3], type = "b", col = "darkorange1", pch = 17)
points(x = threshold, y = power[, 4], type = "b", col = "goldenrod2", pch = 18)

points(x = threshold, y = power[, 5], type = "b", col = "indianred3", pch = 19)
points(x = threshold, y = power[, 6], type = "b", col = "seagreen2", pch = 20)
legend("bottomleft",
       col = c("lightsalmon","darkolivegreen4","darkorange1",
                "goldenrod2", "indianred3", "seagreen2"),
       legend = c(1,2,3,4,5,6), 
     
       pch = c(15:20),
       bty = "n")

plot(x = threshold, y = FDR[, 1], type = "b", col = "lightsalmon",
     xlab = "Threshold", ylab = "FDR", 
     xlim = c(0.55, 0.9), 
     ylim = c(0, 0.29), pch = 15)

points(x = threshold, y = FDR[, 2], type = "b", col = "darkolivegreen4", pch = 16)

points(x = threshold, y = FDR[, 3], type = "b", col = "darkorange1", pch = 17)
points(x = threshold, y = FDR[, 4], type = "b", col = "goldenrod2", pch = 18)

points(x = threshold, y = FDR[, 5], type = "b", col = "indianred3", pch = 19)
points(x = threshold, y = FDR[, 6], type = "b", col = "seagreen2", pch = 20)
legend("topright",
       col = c("lightsalmon","darkolivegreen4","darkorange1",
               "goldenrod2", "indianred3", "seagreen2"),
       legend = c(1,2,3,4,5,6), 
       pch = c(15:20),
       bty = "n")
dev.off()

# Results for the Table
round(sapply(final_results[[1]], function(x)x[["power"]]),2)
round(sapply(final_results[[1]], function(x)x[["FDR"]]),2)

