readRDS("simulation_results.rds")
attach(simulation_results)




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

threshold <- seq(from = 0.55, to = 0.9, by = 0.05)


final_results<- list()
res <- list()
for (j in 1:length(threshold)){
  for (i in 1:length(settings)){
    true_status <- true_signal[[i]]
    est_status <- as.matrix(simulation_results[[i]])
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


pdf("Power_fdr_sim_study_updated.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(x = threshold, y = power[, 1], type = "b", col = "lightsalmon",
     xlab = "Stability selection threshold", ylab = "Power", 
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
       bty = "o",
       title= "Settings",
       bg = "gray97",
       horiz=T
       )

plot(x = threshold, y = FDR[, 1], type = "b", col = "lightsalmon",
     xlab = "Stability selection threshold", ylab = "FDR", 
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
       bty = "o",
       title="Settings",
       bg = "gray97")
dev.off()

# Results for the Table
round(sapply(final_results[[1]], function(x)x[["power"]]),2)
round(sapply(final_results[[1]], function(x)x[["FDR"]]),2)
