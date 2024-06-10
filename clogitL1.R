require(clogitL1)
require(tidyverse)
# function to generate binary matched data 
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
  data_sim <- add_column(data_sim, Y = Y)
  data_sim_arranged<- data_sim %>% 
    group_by(stratum) %>% 
    arrange(desc(Y))
  
  return(list(Y = data_sim_arranged$Y, 
              X = X[data_sim_arranged$obs_id,], 
              stratum = data_sim_arranged$stratum,
              beta = beta))
}

# simulation settings
settings <- list(setting1 = list(p1 = 50, p2 = 50, p1_r = 10, p2_r = 10,
                                 beta1 = 4, beta2 = 4, n = 200),
                 setting4 = list(p1 = 20, p2 = 80, p1_r = 10, p2_r = 10,
                                 beta1 = 4, beta2 = 1))

# true status of covariates (0 for null, 1 for signal)
true_signal <- lapply(settings, function(x) c(rep(1, x$p1_r,),
                                              rep(0, x$p1 - x$p1_r),
                                              rep(1, x$p2_r,),
                                              rep(0, x$p2 - x$p2_r)))



# simulation_study_clogitL1 for 1 run
simulation_study_clogitL1 <- function(p1 = 100, p2 = 100, p1_r = 10, p2_r = 10,
                                         beta1 = 0.5, beta2 = 0.5, cor.coef = 0,
                                         n = 50){
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
  # Output:
  # a list with a vector of estimated coefficients
  
  data_gen <- data_generation(p1 = p1, p2 = p2, p1_r = p1_r,
                              p2_r = p2_r, beta1 = beta1, beta2 = beta2,
                              cor.coef = cor.coef, n = n)
  
  
  fit_clogitL1 <- clogitL1(x = data_gen$X,
                       y = data_gen$Y,
                       strata = data_gen$stratum)
  cv_fit <- cv.clogitL1(fit_clogitL1)
  lambda_clogitL1 <- exp(cv_fit$minCV1se_lambda)
  index <- which(abs(fit_clogitL1$lambda-lambda_clogitL1) == min(abs(fit_clogitL1$lambda-lambda_clogitL1)))
  beta_clogitL1 <- fit_clogitL1$beta[index,]
  return(list(coef = beta_clogitL1))
}


# run nruns of simulations
set.seed(123)
#nruns <- 2
nruns <-  100

start.time <- Sys.time()
temp <- lapply(settings, 
               function(x) 
                 replicate(nruns,
                           simulation_study_clogitL1(p1 = x$p1,
                                                        p2 = x$p2,
                                                        p1_r = x$p1_r,
                                                        p2_r = x$p2_r, 
                                                        beta1 = x$beta1,
                                                        beta2 = x$beta2, 
                                                        n = 200)$coef))


end.time <- Sys.time()
run.time <- end.time - start.time
# collect and summarize results
res <- list()
for(i in 1:length(settings)){
  true_status <- true_signal[[i]]
  est_status <- as.matrix(temp[[i]])
  est_status <- est_status != 0
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

# Power and FDR for the six settings
final_results <- lapply(res, colMeans)
final_results

