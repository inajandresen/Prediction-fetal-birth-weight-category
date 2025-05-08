rm(list=ls())
library(pROC)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
#library(ggpattern)
library(ggpubr)
library(caret)

load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("AGA", "LGA")]

MoM_exprset_STORK_LGA_SGA_T1 <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$TimePoint %in% "1"]
samp_T1 <- pData(MoM_exprset_STORK_LGA_SGA_T1)
samp_T1 <- samp_T1 %>% dplyr::select(Bwkat)
MoM_exprset_STORK_LGA_SGA_T2 <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$TimePoint %in% "2"]
samp_T2 <- pData(MoM_exprset_STORK_LGA_SGA_T2)
samp_T2 <- samp_T2 %>% dplyr::select(Bwkat) 
MoM_exprset_STORK_LGA_SGA_T3 <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$TimePoint %in% "3"]
samp_T3 <- pData(MoM_exprset_STORK_LGA_SGA_T3)
samp_T3 <- samp_T3 %>% dplyr::select(Bwkat)

# T1 EN all proteins----------------------------------------------------
#All proteins EN
load("../Predict_LGA_new/All_proteins/Predict_LGA_T1_EN_all_proteins.RData")
T1_EN_freq <- freq
T1_EN_pile <- pile
T1_EN_f1_values <- f1_values
T1_EN_balanced_acc_values <- balanced_acc_values

T1_EN_a=NULL
for (i in 1:length(T1_EN_pile)){
  T1_EN_a=rbind(T1_EN_a,T1_EN_pile[[i]]$outmat)
}
T1_EN_RC=roc(response=T1_EN_a[,1],predictor=T1_EN_a[,2],direction="<")
T1_EN_auc=round(ci.auc(T1_EN_RC),3)
T1_EN_auc <- round(T1_EN_RC$auc, digits = 2)
T1_EN_ci <- ci(T1_EN_RC, of = "auc")
T1_EN_ci_txt <- "95% CI: 0.34-1.00"


#Visualize F1 score for each iteration
T1_EN_f1_df <- data.frame(Iteration = 1:length(T1_EN_f1_values), F1_Score = T1_EN_f1_values)

ggplot(T1_EN_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T1_EN_a_df <- as.data.frame(T1_EN_a)
T1_EN_a_df$pred_binary = ifelse(T1_EN_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T1_EN_tp = sum(T1_EN_a_df$true_labels == 1 & T1_EN_a_df$pred_binary == 1)
T1_EN_fp = sum(T1_EN_a_df$true_labels == 0 & T1_EN_a_df$pred_binary == 1)
T1_EN_fn = sum(T1_EN_a_df$true_labels == 1 & T1_EN_a_df$pred_binary == 0)

T1_EN_precision = ifelse((T1_EN_tp + T1_EN_fp) == 0, 0, T1_EN_tp / (T1_EN_tp + T1_EN_fp))
T1_EN_recall = ifelse((T1_EN_tp + T1_EN_fn) == 0, 0, T1_EN_tp / (T1_EN_tp + T1_EN_fn))

#Calculate F1 score
if (T1_EN_precision + T1_EN_recall == 0) {
  T1_EN_f1_score <- 0
} else {
  T1_EN_f1_score <- 2 * (T1_EN_precision * T1_EN_recall) / (T1_EN_precision + T1_EN_recall)
}
T1_EN_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T1_EN_ba_df <- data.frame(Iteration = 1:length(T1_EN_balanced_acc_values), Balanced_Accuracy = T1_EN_balanced_acc_values)

ggplot(T1_EN_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T1_EN_sensitivity_all = ifelse(sum(T1_EN_a_df$outcs == 1) == 0, 0, sum(T1_EN_a_df$outcs == 1 & T1_EN_a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
T1_EN_specificity_all = ifelse(sum(T1_EN_a_df$outcs == 0) == 0, 0, sum(T1_EN_a_df$outcs == 0 & T1_EN_a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))

T1_EN_balanced_accuracy = (T1_EN_sensitivity_all + T1_EN_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T1_EN_PR <- pr.curve(scores.class0 = T1_EN_a[,2], weights.class0 = T1_EN_a[,1], curve = TRUE)

# Extract AUC-PR value
T1_EN_auc_pr <- round(T1_EN_PR$auc.integral, digits = 2)
#0.18

plot(T1_EN_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")


T1_EN_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T1_EN_a[,2]>=pr,1,0))), factor((T1_EN_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T1_EN_a[,2]>pr,1,0))), factor(as.character(T1_EN_a[,1])),positive="1")
  T1_EN_tabout=rbind(T1_EN_tabout,c(pr,aa$byClass[1:4]))
}
#            Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.1           0   0.9310345              0      0.9152542
#[2,] 0.1           0   0.9310345              0      0.9152542
> 

# T1 RF All proteins ------------------------------------------------------

load("../Predict_LGA_new/All_proteins/Predict_LGA_T1_RF_all_proteins.RData")
T1_RF_freq <- freq
T1_RF_pile <- pile
T1_RF_f1_values <- f1_values
T1_RF_balanced_acc_values <- balanced_acc_values

T1_RF_a=NULL
for (i in 1:length(T1_RF_pile)){
  T1_RF_a=rbind(T1_RF_a,T1_RF_pile[[i]]$outmat)
}
T1_RF_RC=roc(response=T1_RF_a[,1],predictor=T1_RF_a[,2],direction="<")
T1_RF_auc <- round(T1_RF_RC$auc, digits = 2)
T1_RF_ci <- ci(T1_RF_RC, of = "auc")
T1_RF_ci_txt <- "95% CI: 0.40-0.90"

#Visualize F1 score for each iteration
T1_RF_f1_df <- data.frame(Iteration = 1:length(T1_RF_f1_values), F1_Score = T1_RF_f1_values)

ggplot(T1_RF_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T1_RF_a_df <- as.data.frame(T1_RF_a)
T1_RF_a_df$pred_binary = ifelse(T1_RF_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T1_RF_tp = sum(T1_RF_a_df$true_labels == 1 & T1_RF_a_df$pred_binary == 1)
T1_RF_fp = sum(T1_RF_a_df$true_labels == 0 & T1_RF_a_df$pred_binary == 1)
T1_RF_fn = sum(T1_RF_a_df$true_labels == 1 & T1_RF_a_df$pred_binary == 0)

T1_RF_precision = ifelse((T1_RF_tp + T1_RF_fp) == 0, 0, T1_RF_tp / (T1_RF_tp + T1_RF_fp))
T1_RF_recall = ifelse((T1_RF_tp + T1_RF_fn) == 0, 0, T1_RF_tp / (T1_RF_tp + T1_RF_fn))

#Calculate F1 score
if (T1_RF_precision + T1_RF_recall == 0) {
  T1_RF_f1_score <- 0
} else {
  T1_RF_f1_score <- 2 * (T1_RF_precision * T1_RF_recall) / (T1_RF_precision + T1_RF_recall)
}
T1_RF_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T1_RF_ba_df <- data.frame(Iteration = 1:length(T1_RF_balanced_acc_values), Balanced_Accuracy = T1_RF_balanced_acc_values)

ggplot(T1_RF_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T1_RF_sensitivity_all = ifelse(sum(T1_RF_a_df$outcs == 1) == 0, 0, sum(T1_RF_a_df$outcs == 1 & T1_RF_a_df$pred_binary == 1) / sum(T1_RF_a_df$outcs == 1))
T1_RF_specificity_all = ifelse(sum(T1_RF_a_df$outcs == 0) == 0, 0, sum(T1_RF_a_df$outcs == 0 & T1_RF_a_df$pred_binary == 0) / sum(T1_RF_a_df$outcs == 0))

T1_RF_balanced_accuracy = (T1_RF_sensitivity_all + T1_RF_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T1_RF_PR <- pr.curve(scores.class0 = T1_RF_a[,2], weights.class0 = T1_RF_a[,1], curve = TRUE)

# Extract AUC-PR value
T1_RF_auc_pr <- round(T1_RF_PR$auc.integral, digits = 2)
#0.15

plot(T1_RF_PR, col = "blue", main = "Precision-Recall Curve (Random Forest)")


T1_RF_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T1_RF_a[,2]>=pr,1,0))), factor((T1_RF_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T1_RF_a[,2]>pr,1,0))), factor(as.character(T1_RF_a[,1])),positive="1")
  T1_RF_tabout=rbind(T1_RF_tabout,c(pr,aa$byClass[1:4]))
}
#           Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.1           0   0.9310345              0      0.9152542
#[2,] 0.1           0   0.9310345              0      0.9152542

# T1 EN 20 proteins ----------------------------------------------------------

#20 proteins EN
load("../Predict_LGA_new/20_proteins/Predict_LGA_T1_EN_20_proteins.RData")
T1_EN_20_freq <- freq
T1_EN_20_pile <- pile
T1_EN_20_f1_values <- f1_values
T1_EN_20_balanced_acc_values <- balanced_acc_values

T1_EN_20_a=NULL
for (i in 1:length(T1_EN_20_pile)){
  T1_EN_20_a=rbind(T1_EN_20_a,T1_EN_20_pile[[i]]$outmat)
}
T1_EN_20_RC=roc(response=T1_EN_20_a[,1],predictor=T1_EN_20_a[,2],direction="<")
T1_EN_20_auc <- round(T1_EN_20_RC$auc, digits = 2)
T1_EN_20_ci <- ci(T1_EN_20_RC, of = "auc")
T1_EN_20_ci_txt <- "95% CI: 0.48-1.00"

#Visualize F1 score for each iteration
T1_EN_20_f1_df <- data.frame(Iteration = 1:length(T1_EN_20_f1_values), F1_Score = T1_EN_20_f1_values)

ggplot(T1_EN_20_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T1_EN_20_a_df <- as.data.frame(T1_EN_20_a)
T1_EN_20_a_df$pred_binary = ifelse(T1_EN_20_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T1_EN_20_tp = sum(T1_EN_20_a_df$true_labels == 1 & T1_EN_20_a_df$pred_binary == 1)
T1_EN_20_fp = sum(T1_EN_20_a_df$true_labels == 0 & T1_EN_20_a_df$pred_binary == 1)
T1_EN_20_fn = sum(T1_EN_20_a_df$true_labels == 1 & T1_EN_20_a_df$pred_binary == 0)

T1_EN_20_precision = ifelse((T1_EN_20_tp + T1_EN_20_fp) == 0, 0, T1_EN_20_tp / (T1_EN_20_tp + T1_EN_20_fp))
T1_EN_20_recall = ifelse((T1_EN_20_tp + T1_EN_20_fn) == 0, 0, T1_EN_20_tp / (T1_EN_20_tp + T1_EN_20_fn))

#Calculate F1 score
if (T1_EN_20_precision + T1_EN_20_recall == 0) {
  T1_EN_20_f1_score <- 0
} else {
  T1_EN_20_f1_score <- 2 * (T1_EN_20_precision * T1_EN_20_recall) / (T1_EN_20_precision + T1_EN_20_recall)
}
T1_EN_20_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T1_EN_20_ba_df <- data.frame(Iteration = 1:length(T1_EN_20_balanced_acc_values), Balanced_Accuracy = T1_EN_20_balanced_acc_values)

ggplot(T1_EN_20_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T1_EN_20_sensitivity_all = ifelse(sum(T1_EN_20_a_df$outcs == 1) == 0, 0, sum(T1_EN_20_a_df$outcs == 1 & T1_EN_20_a_df$pred_binary == 1) / sum(T1_EN_20_a_df$outcs == 1))
T1_EN_20_specificity_all = ifelse(sum(T1_EN_20_a_df$outcs == 0) == 0, 0, sum(T1_EN_20_a_df$outcs == 0 & T1_EN_20_a_df$pred_binary == 0) / sum(T1_EN_20_a_df$outcs == 0))

T1_EN_20_balanced_accuracy = (T1_EN_20_sensitivity_all + T1_EN_20_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T1_EN_20_PR <- pr.curve(scores.class0 = T1_EN_20_a[,2], weights.class0 = T1_EN_20_a[,1], curve = TRUE)

# Extract AUC-PR value
T1_EN_20_auc_pr <- round(T1_EN_20_PR$auc.integral, digits = 2)
#0.18

plot(T1_EN_20_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T1_EN_20_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T1_EN_20_a[,2]>=pr,1,0))), factor((T1_EN_20_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T1_EN_20_a[,2]>pr,1,0))), factor(as.character(T1_EN_20_a[,1])),positive="1")
  T1_EN_20_tabout=rbind(T1_EN_20_tabout,c(pr,aa$byClass[1:4]))
}
#            Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.1         0.2           1              1      0.9354839
#[2,] 0.1         0.2           1              1      0.9354839


# T1 RF 20 proteins -------------------------------------------------------


#20 proteins RF
load("../Predict_LGA_new/20_proteins/Predict_LGA_T1_RF_20_proteins.RData")
T1_RF_20_freq <- freq
T1_RF_20_pile <- pile
T1_RF_20_f1_values <- f1_values
T1_RF_20_balanced_acc_values <- balanced_acc_values

T1_RF_20_a=NULL
for (i in 1:length(T1_RF_20_pile)){
  T1_RF_20_a=rbind(T1_RF_20_a,T1_RF_20_pile[[i]]$outmat)
}
T1_RF_20_RC=roc(response=T1_RF_20_a[,1],predictor=T1_RF_20_a[,2],direction="<")
T1_RF_20_auc <- round(T1_RF_20_RC$auc, digits = 2)
T1_RF_20_auc <- "0.85"
T1_RF_20_ci <- ci(T1_RF_20_RC, of = "auc")
T1_RF_20_ci_txt <- "95% CI: 0.72-0.98"

#Visualize F1 score for each iteration
T1_RF_20_f1_df <- data.frame(Iteration = 1:length(T1_RF_20_f1_values), F1_Score = T1_RF_20_f1_values)

ggplot(T1_RF_20_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T1_RF_20_a_df <- as.data.frame(T1_RF_20_a)
T1_RF_20_a_df$pred_binary = ifelse(T1_RF_20_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T1_RF_20_tp = sum(T1_RF_20_a_df$true_labels == 1 & T1_RF_20_a_df$pred_binary == 1)
T1_RF_20_fp = sum(T1_RF_20_a_df$true_labels == 0 & T1_RF_20_a_df$pred_binary == 1)
T1_RF_20_fn = sum(T1_RF_20_a_df$true_labels == 1 & T1_RF_20_a_df$pred_binary == 0)

T1_RF_20_precision = ifelse((T1_RF_20_tp + T1_RF_20_fp) == 0, 0, T1_RF_20_tp / (T1_RF_20_tp + T1_RF_20_fp))
T1_RF_20_recall = ifelse((T1_RF_20_tp + T1_RF_20_fn) == 0, 0, T1_RF_20_tp / (T1_RF_20_tp + T1_RF_20_fn))

#Calculate F1 score
if (T1_RF_20_precision + T1_RF_20_recall == 0) {
  T1_RF_20_f1_score <- 0
} else {
  T1_RF_20_f1_score <- 2 * (T1_RF_20_precision * T1_RF_20_recall) / (T1_RF_20_precision + T1_RF_20_recall)
}
T1_RF_20_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T1_RF_20_ba_df <- data.frame(Iteration = 1:length(T1_RF_20_balanced_acc_values), Balanced_Accuracy = T1_RF_20_balanced_acc_values)

ggplot(T1_RF_20_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T1_RF_20_sensitivity_all = ifelse(sum(T1_RF_20_a_df$outcs == 1) == 0, 0, sum(T1_RF_20_a_df$outcs == 1 & T1_RF_20_a_df$pred_binary == 1) / sum(T1_RF_20_a_df$outcs == 1))
T1_RF_20_specificity_all = ifelse(sum(T1_RF_20_a_df$outcs == 0) == 0, 0, sum(T1_RF_20_a_df$outcs == 0 & T1_RF_20_a_df$pred_binary == 0) / sum(T1_RF_20_a_df$outcs == 0))

T1_RF_20_balanced_accuracy = (T1_RF_20_sensitivity_all + T1_RF_20_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T1_RF_20_PR <- pr.curve(scores.class0 = T1_RF_20_a[,2], weights.class0 = T1_RF_20_a[,1], curve = TRUE)

# Extract AUC-PR value
T1_RF_20_auc_pr <- round(T1_RF_20_PR$auc.integral, digits = 2)
#0.18

plot(T1_RF_20_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T1_RF_20_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T1_RF_20_a[,2]>=pr,1,0))), factor((T1_RF_20_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T1_RF_20_a[,2]>pr,1,0))), factor(as.character(T1_RF_20_a[,1])),positive="1")
  T1_RF_20_tabout=rbind(T1_RF_20_tabout,c(pr,aa$byClass[1:4]))
}
#           Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.12         0.6   0.8965517      0.3333333       0.962963
#[2,] 0.10         0.6   0.8275862      0.2307692       0.960000

# Plot T1 -----------------------------------------------------------------

#Plot
ROC_plot_T1 <- ggroc(list("Random Forest (all proteins)" = T1_RF_RC,
                          "Random Forest (20 proteins)" = T1_RF_20_RC,
                          "Elastic net (all proteins)" = T1_EN_RC,
                          "Elastic net (20 proteins)" = T1_EN_20_RC),
                     aes = c("colour", "size"),
                     legacy.axes = TRUE) +
  #scale_color_manual(values = c("mediumblue", "dodgerblue", "maroon4", "mediumpurple2")) +
  scale_color_manual(values = c("navy", "dodgerblue", "maroon4", "maroon2")) +
  #scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Visit 1\nWeek 12-19") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.text=element_text(size=11),
        axis.title = element_text(size=14)) +
  annotate("text", 
           label= c(paste0('AUC = ', T1_RF_auc, ', ', T1_RF_ci_txt), 
                    paste0('AUC = ', T1_RF_20_auc, ', ', T1_RF_20_ci_txt), 
                    paste0('AUC = ', T1_EN_auc, ', ', T1_EN_ci_txt), 
                    paste0('AUC = ', T1_EN_20_auc, ', ', T1_EN_20_ci_txt)),
           #color = c("mediumblue", "dodgerblue", "maroon4", "mediumpurple2"), 
           color = c("navy", "dodgerblue", "maroon4", "maroon2"),
           size = 4, 
           x = 0.65, 
           y = c(0.2, 0.15, 0.10, 0.05))


#Boxplot
T1_EN_a_mod <- as.data.frame(T1_EN_a)
T1_EN_a_mod$Model <- "Elastic net (All proteins)"

T1_EN_20_a_mod <- as.data.frame(T1_EN_20_a)
T1_EN_20_a_mod$Model <- "Elastic net (20 proteins)"

T1_RF_a_mod <- as.data.frame(T1_RF_a)
T1_RF_a_mod$Model <- "Random Forest (All proteins)"

T1_RF_20_a_mod <- as.data.frame(T1_RF_20_a)
T1_RF_20_a_mod$Model <- "Random Forest (20 proteins)"


T1_EN_RF_res <- rbind(T1_EN_a_mod, T1_EN_20_a_mod, T1_RF_a_mod, T1_RF_20_a_mod)
T1_EN_RF_res$Bwkat <- ifelse(T1_EN_RF_res$true_labels == 1, "LGA", "AGA")
T1_EN_RF_res$Model <- factor(T1_EN_RF_res$Model, levels = c("Random Forest (All proteins)",
                                                            "Random Forest (20 proteins)",
                                                            "Elastic net (All proteins)",
                                                            "Elastic net (20 proteins)"))

T1_box <- ggplot(T1_EN_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "navy",
                              "Random Forest (20 proteins)" = "dodgerblue",
                              "Elastic net (All proteins)" = "maroon4",
                              "Elastic net (20 proteins)" = "maroon2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Visit 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))


# T2 EN All proteins ------------------------------------------------------
#All proteins EN
load("../Predict_LGA_new/All_proteins/Predict_LGA_T2_EN_all_proteins.RData")
T2_EN_freq <- freq
T2_EN_pile <- pile

T2_EN_a=NULL
for (i in 1:length(T2_EN_pile)){
  T2_EN_a=rbind(T2_EN_a,T2_EN_pile[[i]]$outmat)
}
T2_EN_RC=roc(response=T2_EN_a[,1],predictor=T2_EN_a[,2],direction="<")
T2_EN_auc=round(ci.auc(T2_EN_RC),3)

T2_EN_auc <- round(T2_EN_RC$auc, digits = 2)
T2_EN_auc <- "0.70"
T2_EN_ci <- ci(T2_EN_RC, of = "auc")
T2_EN_ci_txt <- "95% CI: 0.42-0.97"

T2_EN_a=NULL
for (i in 1:length(T2_EN_pile)){
  T2_EN_a=rbind(T2_EN_a,T2_EN_pile[[i]]$outmat)
}
T2_EN_RC=roc(response=T2_EN_a[,1],predictor=T2_EN_a[,2],direction="<")
T2_EN_auc=round(ci.auc(T2_EN_RC),3)
T2_EN_auc <- round(T2_EN_RC$auc, digits = 2)
T2_EN_auc <- "0.70"
T2_EN_ci <- ci(T2_EN_RC, of = "auc")
T2_EN_ci_txt <- "95% CI: 0.42-0.97"


#Visualize F1 score for each iteration
T2_EN_f1_df <- data.frame(Iteration = 1:length(T2_EN_f1_values), F1_Score = T2_EN_f1_values)

ggplot(T2_EN_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T2_EN_a_df <- as.data.frame(T2_EN_a)
T2_EN_a_df$pred_binary = ifelse(T2_EN_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T2_EN_tp = sum(T2_EN_a_df$true_labels == 1 & T2_EN_a_df$pred_binary == 1)
T2_EN_fp = sum(T2_EN_a_df$true_labels == 0 & T2_EN_a_df$pred_binary == 1)
T2_EN_fn = sum(T2_EN_a_df$true_labels == 1 & T2_EN_a_df$pred_binary == 0)

T2_EN_precision = ifelse((T2_EN_tp + T2_EN_fp) == 0, 0, T2_EN_tp / (T2_EN_tp + T2_EN_fp))
T2_EN_recall = ifelse((T2_EN_tp + T2_EN_fn) == 0, 0, T2_EN_tp / (T2_EN_tp + T2_EN_fn))

#Calculate F1 score
if (T2_EN_precision + T2_EN_recall == 0) {
  T2_EN_f1_score <- 0
} else {
  T2_EN_f1_score <- 2 * (T2_EN_precision * T2_EN_recall) / (T2_EN_precision + T2_EN_recall)
}
T2_EN_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T2_EN_ba_df <- data.frame(Iteration = 1:length(T2_EN_balanced_acc_values), Balanced_Accuracy = T2_EN_balanced_acc_values)

ggplot(T2_EN_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T2_EN_sensitivity_all = ifelse(sum(T2_EN_a_df$outcs == 1) == 0, 0, sum(T2_EN_a_df$outcs == 1 & T2_EN_a_df$pred_binary == 1) / sum(T2_EN_a_df$outcs == 1))
T2_EN_specificity_all = ifelse(sum(T2_EN_a_df$outcs == 0) == 0, 0, sum(T2_EN_a_df$outcs == 0 & T2_EN_a_df$pred_binary == 0) / sum(T2_EN_a_df$outcs == 0))

T2_EN_balanced_accuracy = (T2_EN_sensitivity_all + T2_EN_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T2_EN_PR <- pr.curve(scores.class0 = T2_EN_a[,2], weights.class0 = T2_EN_a[,1], curve = TRUE)

# Extract AUC-PR value
T2_EN_auc_pr <- round(T2_EN_PR$auc.integral, digits = 2)
#0.18

plot(T2_EN_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T2_EN_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_EN_a[,2]>=pr,1,0))), factor((T2_EN_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_EN_a[,2]>pr,1,0))), factor(as.character(T2_EN_a[,1])),positive="1")
  T2_EN_tabout=rbind(T2_EN_tabout,c(pr,aa$byClass[1:4]))
}
#            Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.1           0   0.9310345              0      0.9152542
#[2,] 0.1           0   0.9310345              0      0.9152542

# T2 RF All proteins ------------------------------------------------------

#All proteins RF
load("../Predict_LGA_new/All_proteins/Predict_LGA_T2_RF_all_proteins.RData")
T2_RF_freq <- freq
T2_RF_pile <- pile

T2_RF_a=NULL
for (i in 1:length(T2_RF_pile)){
  T2_RF_a=rbind(T2_RF_a,T2_RF_pile[[i]]$outmat)
}
T2_RF_RC=roc(response=T2_RF_a[,1],predictor=T2_RF_a[,2],direction="<")
T2_RF_auc <- round(T2_RF_RC$auc, digits = 2)
T2_RF_ci <- ci(T2_RF_RC, of = "auc")
T2_RF_ci_txt <- "95% CI: 0.62-0.93"

T2_RF_a=NULL
for (i in 1:length(T2_RF_pile)){
  T2_RF_a=rbind(T2_RF_a,T2_RF_pile[[i]]$outmat)
}
T2_RF_RC=roc(response=T2_RF_a[,1],predictor=T2_RF_a[,2],direction="<")
T2_RF_auc=round(ci.auc(T2_RF_RC),3)
T2_RF_auc <- round(T2_RF_RC$auc, digits = 2)
T2_RF_ci <- ci(T2_RF_RC, of = "auc")
T2_RF_ci_txt <- "95% CI: 0.62-0.94"


#Visualize F1 score for each iteration
T2_RF_f1_df <- data.frame(Iteration = 1:length(T2_RF_f1_values), F1_Score = T2_RF_f1_values)

ggplot(T2_RF_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T2_RF_a_df <- as.data.frame(T2_RF_a)
T2_RF_a_df$pred_binary = ifelse(T2_RF_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T2_RF_tp = sum(T2_RF_a_df$true_labels == 1 & T2_RF_a_df$pred_binary == 1)
T2_RF_fp = sum(T2_RF_a_df$true_labels == 0 & T2_RF_a_df$pred_binary == 1)
T2_RF_fn = sum(T2_RF_a_df$true_labels == 1 & T2_RF_a_df$pred_binary == 0)

T2_RF_precision = ifelse((T2_RF_tp + T2_RF_fp) == 0, 0, T2_RF_tp / (T2_RF_tp + T2_RF_fp))
T2_RF_recall = ifelse((T2_RF_tp + T2_RF_fn) == 0, 0, T2_RF_tp / (T2_RF_tp + T2_RF_fn))

#Calculate F1 score
if (T2_RF_precision + T2_RF_recall == 0) {
  T2_RF_f1_score <- 0
} else {
  T2_RF_f1_score <- 2 * (T2_RF_precision * T2_RF_recall) / (T2_RF_precision + T2_RF_recall)
}
T2_RF_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T2_RF_ba_df <- data.frame(Iteration = 1:length(T2_RF_balanced_acc_values), Balanced_Accuracy = T2_RF_balanced_acc_values)

ggplot(T2_RF_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T2_RF_sensitivity_all = ifelse(sum(T2_RF_a_df$outcs == 1) == 0, 0, sum(T2_RF_a_df$outcs == 1 & T2_RF_a_df$pred_binary == 1) / sum(T2_RF_a_df$outcs == 1))
T2_RF_specificity_all = ifelse(sum(T2_RF_a_df$outcs == 0) == 0, 0, sum(T2_RF_a_df$outcs == 0 & T2_RF_a_df$pred_binary == 0) / sum(T2_RF_a_df$outcs == 0))

T2_RF_balanced_accuracy = (T2_RF_sensitivity_all + T2_RF_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T2_RF_PR <- pr.curve(scores.class0 = T2_RF_a[,2], weights.class0 = T2_RF_a[,1], curve = TRUE)

# Extract AUC-PR value
T2_RF_auc_pr <- round(T2_RF_PR$auc.integral, digits = 2)
#0.18

plot(T2_RF_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T2_RF_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_RF_a[,2]>=pr,1,0))), factor((T2_RF_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_RF_a[,2]>pr,1,0))), factor(as.character(T2_RF_a[,1])),positive="1")
  T2_RF_tabout=rbind(T2_RF_tabout,c(pr,aa$byClass[1:4]))
}
#             Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.175         0.2   0.8965517      0.1428571      0.9285714
#[2,] 0.115         0.6   0.8103448      0.2142857      0.9591837

# T2 EN 20 proteins ------------------------------------------------------

#20 proteins EN
load("../Predict_LGA_new/20_proteins/Predict_LGA_T2_EN_20_proteins.RData")
T2_EN_20_freq <- freq
T2_EN_20_pile <- pile

T2_EN_20_a=NULL
for (i in 1:length(T2_EN_20_pile)){
  T2_EN_20_a=rbind(T2_EN_20_a,T2_EN_20_pile[[i]]$outmat)
}
T2_EN_20_RC=roc(response=T2_EN_20_a[,1],predictor=T2_EN_20_a[,2],direction="<")
T2_EN_20_auc <- round(T2_EN_20_RC$auc, digits = 2)
T2_EN_20_ci <- ci(T2_EN_20_RC, of = "auc")
T2_EN_20_ci_txt <- "95% CI: 0.70-1.00"

T2_EN_20_a=NULL
for (i in 1:length(T2_EN_20_pile)){
  T2_EN_20_a=rbind(T2_EN_20_a,T2_EN_20_pile[[i]]$outmat)
}
T2_EN_20_RC=roc(response=T2_EN_20_a[,1],predictor=T2_EN_20_a[,2],direction="<")
T2_EN_20_auc=round(ci.auc(T2_EN_20_RC),3)
T2_EN_20_auc <- round(T2_EN_20_RC$auc, digits = 2)
T2_EN_20_ci <- ci(T2_EN_20_RC, of = "auc")
T2_EN_20_ci_txt <- "95% CI: 0.70-1.00"


#Visualize F1 score for each iteration
T2_EN_20_f1_df <- data.frame(Iteration = 1:length(T2_EN_20_f1_values), F1_Score = T2_EN_20_f1_values)

ggplot(T2_EN_20_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T2_EN_20_a_df <- as.data.frame(T2_EN_20_a)
T2_EN_20_a_df$pred_binary = ifelse(T2_EN_20_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T2_EN_20_tp = sum(T2_EN_20_a_df$true_labels == 1 & T2_EN_20_a_df$pred_binary == 1)
T2_EN_20_fp = sum(T2_EN_20_a_df$true_labels == 0 & T2_EN_20_a_df$pred_binary == 1)
T2_EN_20_fn = sum(T2_EN_20_a_df$true_labels == 1 & T2_EN_20_a_df$pred_binary == 0)

T2_EN_20_precision = ifelse((T2_EN_20_tp + T2_EN_20_fp) == 0, 0, T2_EN_20_tp / (T2_EN_20_tp + T2_EN_20_fp))
T2_EN_20_recall = ifelse((T2_EN_20_tp + T2_EN_20_fn) == 0, 0, T2_EN_20_tp / (T2_EN_20_tp + T2_EN_20_fn))

#Calculate F1 score
if (T2_EN_20_precision + T2_EN_20_recall == 0) {
  T2_EN_20_f1_score <- 0
} else {
  T2_EN_20_f1_score <- 2 * (T2_EN_20_precision * T2_EN_20_recall) / (T2_EN_20_precision + T2_EN_20_recall)
}
T2_EN_20_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T2_EN_20_ba_df <- data.frame(Iteration = 1:length(T2_EN_20_balanced_acc_values), Balanced_Accuracy = T2_EN_20_balanced_acc_values)

ggplot(T2_EN_20_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T2_EN_20_sensitivity_all = ifelse(sum(T2_EN_20_a_df$outcs == 1) == 0, 0, sum(T2_EN_20_a_df$outcs == 1 & T2_EN_20_a_df$pred_binary == 1) / sum(T2_EN_20_a_df$outcs == 1))
T2_EN_20_specificity_all = ifelse(sum(T2_EN_20_a_df$outcs == 0) == 0, 0, sum(T2_EN_20_a_df$outcs == 0 & T2_EN_20_a_df$pred_binary == 0) / sum(T2_EN_20_a_df$outcs == 0))

T2_EN_20_balanced_accuracy = (T2_EN_20_sensitivity_all + T2_EN_20_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T2_EN_20_PR <- pr.curve(scores.class0 = T2_EN_20_a[,2], weights.class0 = T2_EN_20_a[,1], curve = TRUE)

# Extract AUC-PR value
T2_EN_20_auc_pr <- round(T2_EN_20_PR$auc.integral, digits = 2)
#0.18

plot(T2_EN_20_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T2_EN_20_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_EN_20_a[,2]>=pr,1,0))), factor((T2_EN_20_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_EN_20_a[,2]>pr,1,0))), factor(as.character(T2_EN_20_a[,1])),positive="1")
  T2_EN_20_tabout=rbind(T2_EN_20_tabout,c(pr,aa$byClass[1:4]))
}
#          Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.1         0.2   0.9655172      0.3333333      0.9333333
#[2,] 0.1         0.2   0.9655172      0.3333333      0.9333333

# T2 RF 20 proteins ------------------------------------------------------

#20 proteins RF
load("../Predict_LGA_new/20_proteins/Predict_LGA_T2_RF_20_proteins.RData")
T2_RF_20_freq <- freq
T2_RF_20_pile <- pile

T2_RF_20_a=NULL
for (i in 1:length(T2_RF_20_pile)){
  T2_RF_20_a=rbind(T2_RF_20_a,T2_RF_20_pile[[i]]$outmat)
}
T2_RF_20_RC=roc(response=T2_RF_20_a[,1],predictor=T2_RF_20_a[,2],direction="<")
T2_RF_20_auc=round(ci.auc(T2_RF_20_RC),3)
T2_RF_20_auc <- round(T2_RF_20_RC$auc, digits = 2)
T2_RF_20_ci <- ci(T2_RF_20_RC, of = "auc")
T2_RF_20_ci_txt <- "95% CI: 0.73-1.00"


#Visualize F1 score for each iteration
T2_RF_20_f1_df <- data.frame(Iteration = 1:length(T2_RF_20_f1_values), F1_Score = T2_RF_20_f1_values)

ggplot(T2_RF_20_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T2_RF_20_a_df <- as.data.frame(T2_RF_20_a)
T2_RF_20_a_df$pred_binary = ifelse(T2_RF_20_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T2_RF_20_tp = sum(T2_RF_20_a_df$true_labels == 1 & T2_RF_20_a_df$pred_binary == 1)
T2_RF_20_fp = sum(T2_RF_20_a_df$true_labels == 0 & T2_RF_20_a_df$pred_binary == 1)
T2_RF_20_fn = sum(T2_RF_20_a_df$true_labels == 1 & T2_RF_20_a_df$pred_binary == 0)

T2_RF_20_precision = ifelse((T2_RF_20_tp + T2_RF_20_fp) == 0, 0, T2_RF_20_tp / (T2_RF_20_tp + T2_RF_20_fp))
T2_RF_20_recall = ifelse((T2_RF_20_tp + T2_RF_20_fn) == 0, 0, T2_RF_20_tp / (T2_RF_20_tp + T2_RF_20_fn))

#Calculate F1 score
if (T2_RF_20_precision + T2_RF_20_recall == 0) {
  T2_RF_20_f1_score <- 0
} else {
  T2_RF_20_f1_score <- 2 * (T2_RF_20_precision * T2_RF_20_recall) / (T2_RF_20_precision + T2_RF_20_recall)
}
T2_RF_20_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T2_RF_20_ba_df <- data.frame(Iteration = 1:length(T2_RF_20_balanced_acc_values), Balanced_Accuracy = T2_RF_20_balanced_acc_values)

ggplot(T2_RF_20_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T2_RF_20_sensitivity_all = ifelse(sum(T2_RF_20_a_df$outcs == 1) == 0, 0, sum(T2_RF_20_a_df$outcs == 1 & T2_RF_20_a_df$pred_binary == 1) / sum(T2_RF_20_a_df$outcs == 1))
T2_RF_20_specificity_all = ifelse(sum(T2_RF_20_a_df$outcs == 0) == 0, 0, sum(T2_RF_20_a_df$outcs == 0 & T2_RF_20_a_df$pred_binary == 0) / sum(T2_RF_20_a_df$outcs == 0))

T2_RF_20_balanced_accuracy = (T2_RF_20_sensitivity_all + T2_RF_20_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T2_RF_20_PR <- pr.curve(scores.class0 = T2_RF_20_a[,2], weights.class0 = T2_RF_20_a[,1], curve = TRUE)

# Extract AUC-PR value
T2_RF_20_auc_pr <- round(T2_RF_20_PR$auc.integral, digits = 2)
#0.18

plot(T2_RF_20_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T2_RF_20_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_RF_20_a[,2]>=pr,1,0))), factor((T2_RF_20_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_RF_20_a[,2]>pr,1,0))), factor(as.character(T2_RF_20_a[,1])),positive="1")
  T2_RF_20_tabout=rbind(T2_RF_20_tabout,c(pr,aa$byClass[1:4]))
}
#           Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.14         0.6   0.9137931          0.375      0.9636364
#[2,] 0.10         0.6   0.7931034          0.200      0.9583333

# T2 plot -----------------------------------------------------------------

#Plot
ROC_plot_T2 <- ggroc(list("Random Forest (all proteins)" = T2_RF_RC,
                          "Random Forest (20 proteins)" = T2_RF_20_RC,
                          "Elastic net (all proteins)" = T2_EN_RC,
                          "Elastic net (20 proteins)" = T2_EN_20_RC),
                     aes = c("colour", "size"),
                     legacy.axes = TRUE) +
  #scale_color_manual(values = c("mediumblue", "dodgerblue", "maroon4", "mediumpurple2")) +
  scale_color_manual(values = c("navy", "dodgerblue", "maroon4", "maroon2")) +
  #scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Visit 2\nWeek 21-27") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.text=element_text(size=11),
        axis.title = element_text(size=14)) +
  annotate("text", 
           label= c(paste0('AUC = ', T2_RF_auc, ', ', T2_RF_ci_txt), 
                    paste0('AUC = ', T2_RF_20_auc, ', ', T2_RF_20_ci_txt), 
                    paste0('AUC = ', T2_EN_auc, ', ', T2_EN_ci_txt), 
                    paste0('AUC = ', T2_EN_20_auc, ', ', T2_EN_20_ci_txt)),
           #color = c("mediumblue", "dodgerblue", "maroon4", "mediumpurple2"),
           color = c("navy", "dodgerblue", "maroon4", "maroon2"),
           size = 4, 
           x = 0.65, 
           y = c(0.2, 0.15, 0.10, 0.05))


#Boxplot
T2_EN_a_mod <- as.data.frame(T2_EN_a)
T2_EN_a_mod$Model <- "Elastic net (All proteins)"

T2_EN_20_a_mod <- as.data.frame(T2_EN_20_a)
T2_EN_20_a_mod$Model <- "Elastic net (20 proteins)"

T2_RF_a_mod <- as.data.frame(T2_RF_a)
T2_RF_a_mod$Model <- "Random Forest (All proteins)"

T2_RF_20_a_mod <- as.data.frame(T2_RF_20_a)
T2_RF_20_a_mod$Model <- "Random Forest (20 proteins)"

T2_EN_RF_res <- rbind(T2_EN_a_mod, T2_EN_20_a_mod, T2_RF_a_mod, T2_RF_20_a_mod)
T2_EN_RF_res$Bwkat <- ifelse(T2_EN_RF_res$true_labels == 1, "LGA", "AGA")
T2_EN_RF_res$Model <- factor(T2_EN_RF_res$Model, levels = c("Random Forest (All proteins)",
                                                            "Random Forest (20 proteins)",
                                                            "Elastic net (All proteins)",
                                                            "Elastic net (20 proteins)"))

T2_box <- ggplot(T2_EN_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "navy",
                              "Random Forest (20 proteins)" = "dodgerblue",
                              "Elastic net (All proteins)" = "maroon4",
                              "Elastic net (20 proteins)" = "maroon2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Visit 2\nWeek 21-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

# T3 EN All proteins  ----------------------------------------------------
#All proteins EN
load("../Predict_LGA_new/All_proteins/Predict_LGA_T3_EN_all_proteins.RData")
T3_EN_freq <- freq
T3_EN_pile <- pile

T3_EN_a=NULL
for (i in 1:length(T3_EN_pile)){
  T3_EN_a=rbind(T3_EN_a,T3_EN_pile[[i]]$outmat)
}
T3_EN_RC=roc(response=T3_EN_a[,1],predictor=T3_EN_a[,2],direction="<")
T3_EN_auc=round(ci.auc(T3_EN_RC),3)
T3_EN_auc <- round(T3_EN_RC$auc, digits = 2)
T3_EN_ci <- ci(T3_EN_RC, of = "auc")
T3_EN_ci_txt <- "95% CI: 0.55-0.87"


#Visualize F1 score for each iteration
T3_EN_f1_df <- data.frame(Iteration = 1:length(T3_EN_f1_values), F1_Score = T3_EN_f1_values)

ggplot(T3_EN_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T3_EN_a_df <- as.data.frame(T3_EN_a)
T3_EN_a_df$pred_binary = ifelse(T3_EN_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T3_EN_tp = sum(T3_EN_a_df$true_labels == 1 & T3_EN_a_df$pred_binary == 1)
T3_EN_fp = sum(T3_EN_a_df$true_labels == 0 & T3_EN_a_df$pred_binary == 1)
T3_EN_fn = sum(T3_EN_a_df$true_labels == 1 & T3_EN_a_df$pred_binary == 0)

T3_EN_precision = ifelse((T3_EN_tp + T3_EN_fp) == 0, 0, T3_EN_tp / (T3_EN_tp + T3_EN_fp))
T3_EN_recall = ifelse((T3_EN_tp + T3_EN_fn) == 0, 0, T3_EN_tp / (T3_EN_tp + T3_EN_fn))

#Calculate F1 score
if (T3_EN_precision + T3_EN_recall == 0) {
  T3_EN_f1_score <- 0
} else {
  T3_EN_f1_score <- 2 * (T3_EN_precision * T3_EN_recall) / (T3_EN_precision + T3_EN_recall)
}
T3_EN_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T3_EN_ba_df <- data.frame(Iteration = 1:length(T3_EN_balanced_acc_values), Balanced_Accuracy = T3_EN_balanced_acc_values)

ggplot(T3_EN_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T3_EN_sensitivity_all = ifelse(sum(T3_EN_a_df$outcs == 1) == 0, 0, sum(T3_EN_a_df$outcs == 1 & T3_EN_a_df$pred_binary == 1) / sum(T3_EN_a_df$outcs == 1))
T3_EN_specificity_all = ifelse(sum(T3_EN_a_df$outcs == 0) == 0, 0, sum(T3_EN_a_df$outcs == 0 & T3_EN_a_df$pred_binary == 0) / sum(T3_EN_a_df$outcs == 0))

T3_EN_balanced_accuracy = (T3_EN_sensitivity_all + T3_EN_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T3_EN_PR <- pr.curve(scores.class0 = T3_EN_a[,2], weights.class0 = T3_EN_a[,1], curve = TRUE)

# Extract AUC-PR value
T3_EN_auc_pr <- round(T3_EN_PR$auc.integral, digits = 2)
#0.18

plot(T3_EN_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T3_EN_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_EN_a[,2]>=pr,1,0))), factor((T3_EN_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_EN_a[,2]>pr,1,0))), factor(as.character(T3_EN_a[,1])),positive="1")
  T3_EN_tabout=rbind(T3_EN_tabout,c(pr,aa$byClass[1:4]))
}
#            Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.1           0   0.9137931              0      0.9137931
#[2,] 0.1           0   0.9137931              0      0.9137931

# T3 RF All proteins  ----------------------------------------------------

#All proteins RF
load("../Predict_LGA_new/All_proteins/Predict_LGA_T3_RF_all_proteins.RData")
T3_RF_freq <- freq
T3_RF_pile <- pile

T3_RF_a=NULL
for (i in 1:length(T3_RF_pile)){
  T3_RF_a=rbind(T3_RF_a,T3_RF_pile[[i]]$outmat)
}
T3_RF_RC=roc(response=T3_RF_a[,1],predictor=T3_RF_a[,2],direction="<")
T3_RF_auc <- round(T3_RF_RC$auc, digits = 2)
T3_RF_ci <- ci(T3_RF_RC, of = "auc")
T3_RF_ci_txt <- "95% CI: 0.58-0.88"

#Visualize F1 score for each iteration
T3_RF_f1_df <- data.frame(Iteration = 1:length(T3_RF_f1_values), F1_Score = T3_RF_f1_values)

ggplot(T3_RF_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T3_RF_a_df <- as.data.frame(T3_RF_a)
T3_RF_a_df$pred_binary = ifelse(T3_RF_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T3_RF_tp = sum(T3_RF_a_df$true_labels == 1 & T3_RF_a_df$pred_binary == 1)
T3_RF_fp = sum(T3_RF_a_df$true_labels == 0 & T3_RF_a_df$pred_binary == 1)
T3_RF_fn = sum(T3_RF_a_df$true_labels == 1 & T3_RF_a_df$pred_binary == 0)

T3_RF_precision = ifelse((T3_RF_tp + T3_RF_fp) == 0, 0, T3_RF_tp / (T3_RF_tp + T3_RF_fp))
T3_RF_recall = ifelse((T3_RF_tp + T3_RF_fn) == 0, 0, T3_RF_tp / (T3_RF_tp + T3_RF_fn))

#Calculate F1 score
if (T3_RF_precision + T3_RF_recall == 0) {
  T3_RF_f1_score <- 0
} else {
  T3_RF_f1_score <- 2 * (T3_RF_precision * T3_RF_recall) / (T3_RF_precision + T3_RF_recall)
}
T3_RF_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T3_RF_ba_df <- data.frame(Iteration = 1:length(T3_RF_balanced_acc_values), Balanced_Accuracy = T3_RF_balanced_acc_values)

ggplot(T3_RF_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T3_RF_sensitivity_all = ifelse(sum(T3_RF_a_df$outcs == 1) == 0, 0, sum(T3_RF_a_df$outcs == 1 & T3_RF_a_df$pred_binary == 1) / sum(T3_RF_a_df$outcs == 1))
T3_RF_specificity_all = ifelse(sum(T3_RF_a_df$outcs == 0) == 0, 0, sum(T3_RF_a_df$outcs == 0 & T3_RF_a_df$pred_binary == 0) / sum(T3_RF_a_df$outcs == 0))

T3_RF_balanced_accuracy = (T3_RF_sensitivity_all + T3_RF_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T3_RF_PR <- pr.curve(scores.class0 = T3_RF_a[,2], weights.class0 = T3_RF_a[,1], curve = TRUE)

# Extract AUC-PR value
T3_RF_auc_pr <- round(T3_RF_PR$auc.integral, digits = 2)
#0.18

plot(T3_RF_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T3_RF_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_RF_a[,2]>=pr,1,0))), factor((T3_RF_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_RF_a[,2]>pr,1,0))), factor(as.character(T3_RF_a[,1])),positive="1")
  T3_RF_tabout=rbind(T3_RF_tabout,c(pr,aa$byClass[1:4]))
}
#             Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.160         0.2   0.8965517     0.14285714      0.9285714
#[2,] 0.105         0.2   0.8103448     0.08333333      0.9215686

# T3 EN 20 proteins  ----------------------------------------------------

#20 proteins EN
load("../Predict_LGA_new/20_proteins/Predict_LGA_T3_EN_20_proteins.RData")
T3_EN_20_freq <- freq
T3_EN_20_pile <- pile

T3_EN_20_a=NULL
for (i in 1:length(T3_EN_20_pile)){
  T3_EN_20_a=rbind(T3_EN_20_a,T3_EN_20_pile[[i]]$outmat)
}
T3_EN_20_RC=roc(response=T3_EN_20_a[,1],predictor=T3_EN_20_a[,2],direction="<")
T3_EN_20_auc=round(ci.auc(T3_EN_20_RC),3)
T3_EN_20_auc <- round(T3_EN_20_RC$auc, digits = 2)
T3_EN_20_ci <- ci(T3_EN_20_RC, of = "auc")
T3_EN_20_ci_txt <- "95% CI: 0.64-0.92"


#Visualize F1 score for each iteration
T3_EN_20_f1_df <- data.frame(Iteration = 1:length(T3_EN_20_f1_values), F1_Score = T3_EN_20_f1_values)

ggplot(T3_EN_20_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T3_EN_20_a_df <- as.data.frame(T3_EN_20_a)
T3_EN_20_a_df$pred_binary = ifelse(T3_EN_20_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T3_EN_20_tp = sum(T3_EN_20_a_df$true_labels == 1 & T3_EN_20_a_df$pred_binary == 1)
T3_EN_20_fp = sum(T3_EN_20_a_df$true_labels == 0 & T3_EN_20_a_df$pred_binary == 1)
T3_EN_20_fn = sum(T3_EN_20_a_df$true_labels == 1 & T3_EN_20_a_df$pred_binary == 0)

T3_EN_20_precision = ifelse((T3_EN_20_tp + T3_EN_20_fp) == 0, 0, T3_EN_20_tp / (T3_EN_20_tp + T3_EN_20_fp))
T3_EN_20_recall = ifelse((T3_EN_20_tp + T3_EN_20_fn) == 0, 0, T3_EN_20_tp / (T3_EN_20_tp + T3_EN_20_fn))

#Calculate F1 score
if (T3_EN_20_precision + T3_EN_20_recall == 0) {
  T3_EN_20_f1_score <- 0
} else {
  T3_EN_20_f1_score <- 2 * (T3_EN_20_precision * T3_EN_20_recall) / (T3_EN_20_precision + T3_EN_20_recall)
}
T3_EN_20_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T3_EN_20_ba_df <- data.frame(Iteration = 1:length(T3_EN_20_balanced_acc_values), Balanced_Accuracy = T3_EN_20_balanced_acc_values)

ggplot(T3_EN_20_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T3_EN_20_sensitivity_all = ifelse(sum(T3_EN_20_a_df$outcs == 1) == 0, 0, sum(T3_EN_20_a_df$outcs == 1 & T3_EN_20_a_df$pred_binary == 1) / sum(T3_EN_20_a_df$outcs == 1))
T3_EN_20_specificity_all = ifelse(sum(T3_EN_20_a_df$outcs == 0) == 0, 0, sum(T3_EN_20_a_df$outcs == 0 & T3_EN_20_a_df$pred_binary == 0) / sum(T3_EN_20_a_df$outcs == 0))

T3_EN_20_balanced_accuracy = (T3_EN_20_sensitivity_all + T3_EN_20_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T3_EN_20_PR <- pr.curve(scores.class0 = T3_EN_20_a[,2], weights.class0 = T3_EN_20_a[,1], curve = TRUE)

# Extract AUC-PR value
T3_EN_20_auc_pr <- round(T3_EN_20_PR$auc.integral, digits = 2)
#0.18

plot(T3_EN_20_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T3_EN_20_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_EN_20_a[,2]>=pr,1,0))), factor((T3_EN_20_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_EN_20_a[,2]>pr,1,0))), factor(as.character(T3_EN_20_a[,1])),positive="1")
  T3_EN_20_tabout=rbind(T3_EN_20_tabout,c(pr,aa$byClass[1:4]))
}
#Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.1         0.2   0.9482759           0.25      0.9322034
#[2,] 0.1         0.2   0.9482759           0.25      0.9322034

# T3 RF 20 proteins  ----------------------------------------------------

#20 proteins RF
load("../Predict_LGA_new/20_proteins/Predict_LGA_T3_RF_20_proteins.RData")
T3_RF_20_freq <- freq
T3_RF_20_pile <- pile

T3_RF_20_a=NULL
for (i in 1:length(T3_RF_20_pile)){
  T3_RF_20_a=rbind(T3_RF_20_a,T3_RF_20_pile[[i]]$outmat)
}
T3_RF_20_RC=roc(response=T3_RF_20_a[,1],predictor=T3_RF_20_a[,2],direction="<")
T3_RF_20_auc=round(ci.auc(T3_RF_20_RC),3)
T3_RF_20_auc <- round(T3_RF_20_RC$auc, digits = 2)
T3_RF_20_auc <- "0.80"
T3_RF_20_ci <- ci(T3_RF_20_RC, of = "auc")
T3_RF_20_ci_txt <- "95% CI: 0.59-1.00"


#Visualize F1 score for each iteration
T3_RF_20_f1_df <- data.frame(Iteration = 1:length(T3_RF_20_f1_values), F1_Score = T3_RF_20_f1_values)

ggplot(T3_RF_20_f1_df, aes(x = Iteration, y = F1_Score)) +
  geom_line(color = "green") +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "F1 Score Across LOOCV Iterations",
       x = "Iteration",
       y = "F1 Score")

# Compute overall F1 score across all iterations
T3_RF_20_a_df <- as.data.frame(T3_RF_20_a)
T3_RF_20_a_df$pred_binary = ifelse(T3_RF_20_a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)

T3_RF_20_tp = sum(T3_RF_20_a_df$true_labels == 1 & T3_RF_20_a_df$pred_binary == 1)
T3_RF_20_fp = sum(T3_RF_20_a_df$true_labels == 0 & T3_RF_20_a_df$pred_binary == 1)
T3_RF_20_fn = sum(T3_RF_20_a_df$true_labels == 1 & T3_RF_20_a_df$pred_binary == 0)

T3_RF_20_precision = ifelse((T3_RF_20_tp + T3_RF_20_fp) == 0, 0, T3_RF_20_tp / (T3_RF_20_tp + T3_RF_20_fp))
T3_RF_20_recall = ifelse((T3_RF_20_tp + T3_RF_20_fn) == 0, 0, T3_RF_20_tp / (T3_RF_20_tp + T3_RF_20_fn))

#Calculate F1 score
if (T3_RF_20_precision + T3_RF_20_recall == 0) {
  T3_RF_20_f1_score <- 0
} else {
  T3_RF_20_f1_score <- 2 * (T3_RF_20_precision * T3_RF_20_recall) / (T3_RF_20_precision + T3_RF_20_recall)
}
T3_RF_20_f1_score
#0

#Balanced accuracy 
# Convert stored balanced accuracy values into a data frame
T3_RF_20_ba_df <- data.frame(Iteration = 1:length(T3_RF_20_balanced_acc_values), Balanced_Accuracy = T3_RF_20_balanced_acc_values)

ggplot(T3_RF_20_ba_df, aes(x = Iteration, y = Balanced_Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Balanced Accuracy Over LOOCV Iterations",
       x = "Iteration",
       y = "Balanced Accuracy")

# Compute overall Balanced Accuracy across all iterations
T3_RF_20_sensitivity_all = ifelse(sum(T3_RF_20_a_df$outcs == 1) == 0, 0, sum(T3_RF_20_a_df$outcs == 1 & T3_RF_20_a_df$pred_binary == 1) / sum(T3_RF_20_a_df$outcs == 1))
T3_RF_20_specificity_all = ifelse(sum(T3_RF_20_a_df$outcs == 0) == 0, 0, sum(T3_RF_20_a_df$outcs == 0 & T3_RF_20_a_df$pred_binary == 0) / sum(T3_RF_20_a_df$outcs == 0))

T3_RF_20_balanced_accuracy = (T3_RF_20_sensitivity_all + T3_RF_20_specificity_all) / 2
#0

# Compute Precision-Recall AUC
T3_RF_20_PR <- pr.curve(scores.class0 = T3_RF_20_a[,2], weights.class0 = T3_RF_20_a[,1], curve = TRUE)

# Extract AUC-PR value
T3_RF_20_auc_pr <- round(T3_RF_20_PR$auc.integral, digits = 2)
#0.18

plot(T3_RF_20_PR, col = "blue", main = "Precision-Recall Curve (Elastic Net)")

T3_RF_20_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_RF_20_a[,2]>=pr,1,0))), factor((T3_RF_20_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_RF_20_a[,2]>pr,1,0))), factor(as.character(T3_RF_20_a[,1])),positive="1")
  T3_RF_20_tabout=rbind(T3_RF_20_tabout,c(pr,aa$byClass[1:4]))
}
#              Sensitivity Specificity Pos Pred Value Neg Pred Value
#[1,] 0.15         0.4   0.8965517           0.25      0.9454545
#[2,] 0.10         0.6   0.8448276           0.25      0.9607843

# T3 Plot  ----------------------------------------------------

#Plot
ROC_plot_T3 <- ggroc(list("Random Forest (all proteins)" = T3_RF_RC,
                          "Random Forest (20 proteins)" = T3_RF_20_RC,
                          "Elastic net (all proteins)" = T3_EN_RC,
                          "Elastic net (20 proteins)" = T3_EN_20_RC),
                     aes = c("colour", "size"),
                     legacy.axes = TRUE) +
  #scale_color_manual(values = c("mediumblue", "dodgerblue", "maroon4", "mediumpurple2")) +
  scale_color_manual(values = c("navy", "dodgerblue", "maroon4", "maroon2")) +
  #scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Visit 3\nWeek 28-34") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.text=element_text(size=11),
        axis.title = element_text(size=14)) +
  annotate("text", 
           label= c(paste0('AUC = ', T3_RF_auc, ', ', T3_RF_ci_txt), 
                    paste0('AUC = ', T3_RF_20_auc, ', ', T3_RF_20_ci_txt), 
                    paste0('AUC = ', T3_EN_auc, ', ', T3_EN_ci_txt), 
                    paste0('AUC = ', T3_EN_20_auc, ', ', T3_EN_20_ci_txt)),
           #color = c("mediumblue", "dodgerblue", "maroon4", "mediumpurple2"), 
           color = c("navy", "dodgerblue", "maroon4", "maroon2"),
           size = 4, 
           x = 0.65, 
           y = c(0.2, 0.15, 0.10, 0.05))


#Boxplot
T3_EN_a_mod <- as.data.frame(T3_EN_a)
T3_EN_a_mod$Model <- "Elastic net (All proteins)"

T3_EN_20_a_mod <- as.data.frame(T3_EN_20_a)
T3_EN_20_a_mod$Model <- "Elastic net (20 proteins)"

T3_RF_a_mod <- as.data.frame(T3_RF_a)
T3_RF_a_mod$Model <- "Random Forest (All proteins)"

T3_RF_20_a_mod <- as.data.frame(T3_RF_20_a)
T3_RF_20_a_mod$Model <- "Random Forest (20 proteins)"

T3_EN_RF_res <- rbind(T3_EN_a_mod, T3_EN_20_a_mod, T3_RF_a_mod, T3_RF_20_a_mod)
T3_EN_RF_res$Bwkat <- ifelse(T3_EN_RF_res$true_labels == 1, "LGA", "AGA")
T3_EN_RF_res$Model <- factor(T3_EN_RF_res$Model, levels = c("Random Forest (All proteins)",
                                                            "Random Forest (20 proteins)",
                                                            "Elastic net (All proteins)",
                                                            "Elastic net (20 proteins)"))


T3_box <- ggplot(T3_EN_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "navy",
                              "Random Forest (20 proteins)" = "dodgerblue",
                              "Elastic net (All proteins)" = "maroon4",
                              "Elastic net (20 proteins)" = "maroon2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Visit 3\nWeek 28-34") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

# Arrange and save plots --------------------------------------------------
get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

for_legend_roc <- ggroc(list("Random Forest (all proteins)" = T3_RF_RC,
                             "Random Forest (20 proteins)" = T3_RF_20_RC,
                             "Elastic net (all proteins)" = T3_EN_RC,
                             "Elastic net (20 proteins)" = T3_EN_20_RC),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  #scale_color_manual(values = c("mediumblue", "dodgerblue", "maroon4", "mediumpurple2")) +
  scale_color_manual(values = c("navy", "dodgerblue", "maroon4", "maroon2")) +
  #scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=13)) 


legend_roc <- get_legend(for_legend_roc)


for_legend_box <- ggplot(T3_EN_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "navy",
                              "Random Forest (20 proteins)" = "dodgerblue",
                              "Elastic net (All proteins)" = "maroon4",
                              "Elastic net (20 proteins)" = "maroon2")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),)

legend_box <- get_legend(for_legend_box)


ROCplots <- grid.arrange(ROC_plot_T1, ROC_plot_T2, ROC_plot_T3, legend_roc,
                         ncol = 4, 
                         widths = c(0.8, 0.8, 0.8, 0.6))

BOXplots <- grid.arrange(T1_box, T2_box, T3_box, legend_box,
                         ncol = 4, 
                         widths = c(0.8, 0.8, 0.8, 0.6))


all_plots <- ggarrange(ROCplots, BOXplots,
                       nrow = 2, 
                       #heights = c(1.3, 1),
                       heights = c(1, 1),
                       labels = c("A)", "B)"),
                       font.label = list(size = 20))


ggsave(filename= "ROC_Box_LGA_all_and_20_proteins.png",
       plot = all_plots,
       device = "png",
       path = "../Predict_LGA_new/",
       width = 35,
       height = 20,
       units = "cm",
       bg = "white")

ROCplots_LGA <- ROCplots
save(ROCplots_LGA, file = "../Predict_LGA_new/ROCplots_LGA.RData")

