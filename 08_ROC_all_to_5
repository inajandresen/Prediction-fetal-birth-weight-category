rm(list=ls())
library(pROC)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
#library(ggpattern)
library(ggpubr)
library(PRROC)

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

# T1 Elastic net ----------------------------------------------------
#All proteins EN
load("../All_proteins/Predict_LGA_T1_EN_all_proteins.RData")
T1_EN_freq <- freq
T1_EN_pile <- pile

T1_EN_a=NULL
for (i in 1:length(T1_EN_pile)){
  T1_EN_a=rbind(T1_EN_a,T1_EN_pile[[i]]$outmat)
}
T1_EN_RC=roc(response=T1_EN_a[,1],predictor=T1_EN_a[,2],direction="<")
T1_EN_auc <- round(T1_EN_RC$auc, digits = 2)
T1_EN_ci <- ci(T1_EN_RC, of = "auc")
T1_EN_ci_txt <- "95% CI: 0.34-1.00"



#50 proteins EN
load("../50_proteins/Predict_LGA_T1_EN_50_proteins.RData")

T1_EN_50_freq <- freq
T1_EN_50_pile <- pile

T1_EN_50_a=NULL
for (i in 1:length(T1_EN_50_pile)){
  T1_EN_50_a=rbind(T1_EN_50_a,T1_EN_50_pile[[i]]$outmat)
}
T1_EN_50_RC=roc(response=T1_EN_50_a[,1],predictor=T1_EN_50_a[,2],direction="<")

T1_EN_50_auc <- round(T1_EN_50_RC$auc, digits = 2)
T1_EN_50_ci <- ci(T1_EN_50_RC, of = "auc")
T1_EN_50_ci_txt <- "95% CI: 0.43-1.00"


#40 proteins
load("../40_proteins/Predict_LGA_T1_EN_40_proteins.RData")

T1_EN_40_freq <- freq
T1_EN_40_pile <- pile

T1_EN_40_a=NULL
for (i in 1:length(T1_EN_40_pile)){
  T1_EN_40_a=rbind(T1_EN_40_a,T1_EN_40_pile[[i]]$outmat)
}
T1_EN_40_RC=roc(response=T1_EN_40_a[,1],predictor=T1_EN_40_a[,2],direction="<")
T1_EN_40_auc <- round(T1_EN_40_RC$auc, digits = 2)
T1_EN_40_ci <- ci(T1_EN_40_RC, of = "auc")
T1_EN_40_ci_txt <- "95% CI: 0.52-1.00"

#30 proteins EN
load("../30_proteins/Predict_LGA_T1_EN_30_proteins.RData")

T1_EN_30_freq <- freq
T1_EN_30_pile <- pile

T1_EN_30_a=NULL
for (i in 1:length(T1_EN_30_pile)){
  T1_EN_30_a=rbind(T1_EN_30_a,T1_EN_30_pile[[i]]$outmat)
}
T1_EN_30_RC=roc(response=T1_EN_30_a[,1],predictor=T1_EN_30_a[,2],direction="<")
T1_EN_30_auc <- round(T1_EN_30_RC$auc, digits = 2)
T1_EN_30_ci <- ci(T1_EN_30_RC, of = "auc")
#T1_EN_30_ci_txt <- "95% CI: 0.55-0.85"
T1_EN_30_ci_txt <- "95% CI: 0.61-1.00"


#20 proteins EN
load("../20_proteins/Predict_LGA_T1_EN_20_proteins.RData")

T1_EN_20_freq <- freq
T1_EN_20_pile <- pile

T1_EN_20_a=NULL
for (i in 1:length(T1_EN_20_pile)){
  T1_EN_20_a=rbind(T1_EN_20_a,T1_EN_20_pile[[i]]$outmat)
}
T1_EN_20_RC=roc(response=T1_EN_20_a[,1],predictor=T1_EN_20_a[,2],direction="<")
T1_EN_20_auc <- round(T1_EN_20_RC$auc, digits = 2)
T1_EN_20_ci <- ci(T1_EN_20_RC, of = "auc")
T1_EN_20_ci_txt <- "95% CI: 0.48-1.00"


#10 proteins EN
load("../10_proteins/Predict_LGA_T1_EN_10_proteins.RData")

T1_EN_10_freq <- freq
T1_EN_10_pile <- pile

T1_EN_10_a=NULL
for (i in 1:length(T1_EN_10_pile)){
  T1_EN_10_a=rbind(T1_EN_10_a,T1_EN_10_pile[[i]]$outmat)
}
T1_EN_10_RC=roc(response=T1_EN_10_a[,1],predictor=T1_EN_10_a[,2],direction="<")
T1_EN_10_auc <- round(T1_EN_10_RC$auc, digits = 2)
T1_EN_10_auc <- "0.90"
T1_EN_10_ci <- ci(T1_EN_10_RC, of = "auc")
#T1_EN_10_ci_txt <- "95% CI: 0.19-0.64"
T1_EN_10_ci_txt <- "95% CI: 0.73-1.00"

#5 proteins EN
load("../5_proteins/Predict_LGA_T1_EN_5_proteins.RData")

T1_EN_5_freq <- freq
T1_EN_5_pile <- pile

T1_EN_5_a=NULL
for (i in 1:length(T1_EN_5_pile)){
  T1_EN_5_a=rbind(T1_EN_5_a,T1_EN_5_pile[[i]]$outmat)
}
T1_EN_5_RC=roc(response=T1_EN_5_a[,1],predictor=T1_EN_5_a[,2],direction="<")
T1_EN_5_auc <- round(T1_EN_5_RC$auc, digits = 2)
T1_EN_5_ci <- ci(T1_EN_5_RC, of = "auc")
T1_EN_5_ci_txt <- "95% CI: 0.45-1.00"


#Plot
ROC_plot_EN_T1 <- ggroc(list("All proteins" = T1_EN_RC,
                             "50 proteins" = T1_EN_50_RC,
                             "40 proteins" = T1_EN_40_RC,
                             "30 proteins" = T1_EN_30_RC,
                             "20 proteins" = T1_EN_20_RC,
                             "10 proteins" = T1_EN_10_RC,
                             "5 proteins" = T1_EN_5_RC),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown3")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Elastic net\nVisit 1 (Week 12-19)") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T1_EN_auc, ', ', T1_EN_ci_txt),
                    paste0('AUC = ', T1_EN_50_auc, ', ', T1_EN_50_ci_txt), 
                    paste0('AUC = ', T1_EN_40_auc, ', ', T1_EN_40_ci_txt), 
                    paste0('AUC = ', T1_EN_30_auc, ', ', T1_EN_30_ci_txt), 
                    paste0('AUC = ', T1_EN_20_auc, ', ', T1_EN_20_ci_txt),
                    paste0('AUC = ', T1_EN_10_auc, ', ', T1_EN_10_ci_txt),
                    paste0('AUC = ', T1_EN_5_auc, ', ', T1_EN_5_ci_txt)),
           color = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2"), 
           size = 3, 
           x = 0.75, 
           y = c(0.30,0.25, 0.2, 0.15, 0.10, 0.05, 0.0))

#Boxplot
T1_EN_a_mod <- as.data.frame(T1_EN_a)
T1_EN_a_mod$Model <- "All proteins"

T1_EN_50_a_mod <- as.data.frame(T1_EN_50_a)
T1_EN_50_a_mod$Model <- "50 proteins"

T1_EN_40_a_mod <- as.data.frame(T1_EN_40_a)
T1_EN_40_a_mod$Model <- "40 proteins"

T1_EN_30_a_mod <- as.data.frame(T1_EN_30_a)
T1_EN_30_a_mod$Model <- "30 proteins"

T1_EN_20_a_mod <- as.data.frame(T1_EN_20_a)
T1_EN_20_a_mod$Model <- "20 proteins"

T1_EN_10_a_mod <- as.data.frame(T1_EN_10_a)
T1_EN_10_a_mod$Model <- "10 proteins"

T1_EN_5_a_mod <- as.data.frame(T1_EN_5_a)
T1_EN_5_a_mod$Model <- "5 proteins"


T1_EN_res <- rbind(T1_EN_a_mod, T1_EN_50_a_mod, T1_EN_40_a_mod, T1_EN_30_a_mod, T1_EN_20_a_mod, T1_EN_10_a_mod, T1_EN_5_a_mod)
T1_EN_res$Bwkat <- ifelse(T1_EN_res$true_labels == 1, "LGA", "AGA")
T1_EN_res$Model <- factor(T1_EN_res$Model, levels = c("All proteins",
                                                            "50 proteins",
                                                            "40 proteins",
                                                            "30 proteins",
                                                            "20 proteins",
                                                            "10 proteins",
                                                            "5 proteins"))

T1_EN_box <- ggplot(T1_EN_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "maroon4",
                              "50 proteins" = "mediumpurple4",
                              "40 proteins" = "mediumorchid2",
                              "30 proteins" = "mediumpurple2",
                              "20 proteins" = "maroon2",
                              "10 proteins" = "pink3",
                              "5 proteins" = "brown2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  #ggtitle("Visit 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

# Compute Precision-Recall AUC
T1_EN_PR <- pr.curve(scores.class0 = T1_EN_a[,2], weights.class0 = T1_EN_a[,1], curve = TRUE)
T1_EN_PR_auc <- round(T1_EN_PR$auc.integral, digits = 2)
T1_EN_PR_df <- data.frame(Recall = T1_EN_PR$curve[, 1], Precision = T1_EN_PR$curve[, 2], Model = "All proteins")

T1_EN_50_PR <- pr.curve(scores.class0 = T1_EN_50_a[,2], weights.class0 = T1_EN_50_a[,1], curve = TRUE)
T1_EN_50_PR_auc <- round(T1_EN_50_PR$auc.integral, digits = 2)
T1_EN_50_PR_df <- data.frame(Recall = T1_EN_50_PR$curve[, 1], Precision = T1_EN_50_PR$curve[, 2], Model = "50 proteins")

T1_EN_40_PR <- pr.curve(scores.class0 = T1_EN_40_a[,2], weights.class0 = T1_EN_40_a[,1], curve = TRUE)
T1_EN_40_PR_auc <- round(T1_EN_40_PR$auc.integral, digits = 2)
T1_EN_40_PR_df <- data.frame(Recall = T1_EN_40_PR$curve[, 1], Precision = T1_EN_40_PR$curve[, 2], Model = "40 proteins")

T1_EN_30_PR <- pr.curve(scores.class0 = T1_EN_30_a[,2], weights.class0 = T1_EN_30_a[,1], curve = TRUE)
T1_EN_30_PR_auc <- round(T1_EN_30_PR$auc.integral, digits = 2)
T1_EN_30_PR_df <- data.frame(Recall = T1_EN_30_PR$curve[, 1], Precision = T1_EN_30_PR$curve[, 2], Model = "30 proteins")

T1_EN_20_PR <- pr.curve(scores.class0 = T1_EN_20_a[,2], weights.class0 = T1_EN_20_a[,1], curve = TRUE)
T1_EN_20_PR_auc <- round(T1_EN_20_PR$auc.integral, digits = 2)
T1_EN_20_PR_df <- data.frame(Recall = T1_EN_20_PR$curve[, 1], Precision = T1_EN_20_PR$curve[, 2], Model = "20 proteins")

T1_EN_10_PR <- pr.curve(scores.class0 = T1_EN_10_a[,2], weights.class0 = T1_EN_10_a[,1], curve = TRUE)
T1_EN_10_PR_auc <- round(T1_EN_10_PR$auc.integral, digits = 2)
T1_EN_10_PR_df <- data.frame(Recall = T1_EN_10_PR$curve[, 1], Precision = T1_EN_10_PR$curve[, 2], Model = "10 proteins")

T1_EN_5_PR <- pr.curve(scores.class0 = T1_EN_5_a[,2], weights.class0 = T1_EN_5_a[,1], curve = TRUE)
T1_EN_5_PR_auc <- round(T1_EN_5_PR$auc.integral, digits = 2)
T1_EN_5_PR_df <- data.frame(Recall = T1_EN_5_PR$curve[, 1], Precision = T1_EN_5_PR$curve[, 2], Model = "5 proteins")

T1_EN_pr_data <- bind_rows(T1_EN_PR_df, T1_EN_50_PR_df, T1_EN_40_PR_df, T1_EN_30_PR_df, T1_EN_20_PR_df, T1_EN_10_PR_df, T1_EN_5_PR_df)

T1_EN_pr_data$Model <- factor(T1_EN_pr_data$Model, levels = c("All proteins",
                                                              "50 proteins",
                                                              "40 proteins",
                                                              "30 proteins",
                                                              "20 proteins",
                                                              "10 proteins",
                                                              "5 proteins"))


T1_EN_PR <- ggplot(T1_EN_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "maroon4",
                              "50 proteins" = "mediumpurple4",
                              "40 proteins" = "mediumorchid2",
                              "30 proteins" = "mediumpurple2",
                              "20 proteins" = "maroon2",
                              "10 proteins" = "pink3",
                              "5 proteins" = "brown2")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('PR-AUC = ', T1_EN_PR_auc),
                    paste0('PR-AUC = ', T1_EN_50_PR_auc), 
                    paste0('PR-AUC = ', T1_EN_40_PR_auc), 
                    paste0('PR-AUC = ', T1_EN_30_PR_auc), 
                    paste0('PR-AUC = ', T1_EN_20_PR_auc),
                    paste0('PR-AUC = ', T1_EN_10_PR_auc),
                    paste0('PR-AUC = ', T1_EN_5_PR_auc)),
           color = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2"), 
           size = 3, 
           x = 0.9, 
           y = c(0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4)) 
  


# T2 Elastic net ----------------------------------------------------
#All proteins EN
load("../All_proteins/Predict_LGA_T2_EN_all_proteins.RData")

T2_EN_freq <- freq
T2_EN_pile <- pile

T2_EN_a=NULL
for (i in 1:length(T2_EN_pile)){
  T2_EN_a=rbind(T2_EN_a,T2_EN_pile[[i]]$outmat)
}
T2_EN_RC=roc(response=T2_EN_a[,1],predictor=T2_EN_a[,2],direction="<")
T2_EN_auc <- round(T2_EN_RC$auc, digits = 2)
T2_EN_auc <- "0.70"
T2_EN_ci <- ci(T2_EN_RC, of = "auc")
T2_EN_ci_txt <- "95% CI: 0.42-0.97"


#50 proteins EN
load("../50_proteins/Predict_LGA_T2_EN_50_proteins.RData")
T2_EN_50_freq <- freq
T2_EN_50_pile <- pile

T2_EN_50_a=NULL
for (i in 1:length(T2_EN_50_pile)){
  T2_EN_50_a=rbind(T2_EN_50_a,T2_EN_50_pile[[i]]$outmat)
}
T2_EN_50_RC=roc(response=T2_EN_50_a[,1],predictor=T2_EN_50_a[,2],direction="<")

T2_EN_50_auc <- round(T2_EN_50_RC$auc, digits = 2)
T2_EN_50_ci <- ci(T2_EN_50_RC, of = "auc")
T2_EN_50_ci_txt <- "95% CI: 0.64-0.98"


#40 proteins
load("../40_proteins/Predict_LGA_T2_EN_40_proteins.RData")
T2_EN_40_freq <- freq
T2_EN_40_pile <- pile

T2_EN_40_a=NULL
for (i in 1:length(T2_EN_40_pile)){
  T2_EN_40_a=rbind(T2_EN_40_a,T2_EN_40_pile[[i]]$outmat)
}
T2_EN_40_RC=roc(response=T2_EN_40_a[,1],predictor=T2_EN_40_a[,2],direction="<")

T2_EN_40_auc <- round(T2_EN_40_RC$auc, digits = 2)
T2_EN_40_ci <- ci(T2_EN_40_RC, of = "auc")
T2_EN_40_ci_txt <- "95% CI: 0.64-0.98"

#30 proteins EN
load("../30_proteins/Predict_LGA_T2_EN_30_proteins.RData")
T2_EN_30_freq <- freq
T2_EN_30_pile <- pile

T2_EN_30_a=NULL
for (i in 1:length(T2_EN_30_pile)){
  T2_EN_30_a=rbind(T2_EN_30_a,T2_EN_30_pile[[i]]$outmat)
}
T2_EN_30_RC=roc(response=T2_EN_30_a[,1],predictor=T2_EN_30_a[,2],direction="<")

T2_EN_30_auc <- round(T2_EN_30_RC$auc, digits = 2)
T2_EN_30_auc <- "0.80"
T2_EN_30_ci <- ci(T2_EN_30_RC, of = "auc")
T2_EN_30_ci_txt <- "95% CI: 0.64-0.96"


#20 proteins EN
load("../20_proteins/Predict_LGA_T2_EN_20_proteins.RData")
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

#10 proteins EN
load("../10_proteins/Predict_LGA_T2_EN_10_proteins.RData")
T2_EN_10_freq <- freq
T2_EN_10_pile <- pile

T2_EN_10_a=NULL
for (i in 1:length(T2_EN_10_pile)){
  T2_EN_10_a=rbind(T2_EN_10_a,T2_EN_10_pile[[i]]$outmat)
}
T2_EN_10_RC=roc(response=T2_EN_10_a[,1],predictor=T2_EN_10_a[,2],direction="<")
T2_EN_10_auc <- round(T2_EN_10_RC$auc, digits = 2)
T2_EN_10_auc <- "0.80"
T2_EN_10_ci <- ci(T2_EN_10_RC, of = "auc")
T2_EN_10_ci_txt <- "95% CI: 0.65-0.95"

#5 proteins EN
load("../5_proteins/Predict_LGA_T2_EN_5_proteins.RData")
T2_EN_5_freq <- freq
T2_EN_5_pile <- pile

T2_EN_5_a=NULL
for (i in 1:length(T2_EN_5_pile)){
  T2_EN_5_a=rbind(T2_EN_5_a,T2_EN_5_pile[[i]]$outmat)
}
T2_EN_5_RC=roc(response=T2_EN_5_a[,1],predictor=T2_EN_5_a[,2],direction="<")
T2_EN_5_auc <- round(T2_EN_5_RC$auc, digits = 2)
T2_EN_5_ci <- ci(T2_EN_5_RC, of = "auc")
T2_EN_5_ci_txt <- "95% CI: 0.18-0.79"


#Plot
ROC_plot_EN_T2 <- ggroc(list("All proteins" = T2_EN_RC,
                             "50 proteins" = T2_EN_50_RC,
                             "40 proteins" = T2_EN_40_RC,
                             "30 proteins" = T2_EN_30_RC,
                             "20 proteins" = T2_EN_20_RC,
                             "10 proteins" = T2_EN_10_RC,
                             "5 proteins" = T2_EN_5_RC),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Elastic net\nVisit 2 (Week 21-27)") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T2_EN_auc, ', ', T2_EN_ci_txt),
                    paste0('AUC = ', T2_EN_50_auc, ', ', T2_EN_50_ci_txt), 
                    paste0('AUC = ', T2_EN_40_auc, ', ', T2_EN_40_ci_txt), 
                    paste0('AUC = ', T2_EN_30_auc, ', ', T2_EN_30_ci_txt), 
                    paste0('AUC = ', T2_EN_20_auc, ', ', T2_EN_20_ci_txt),
                    paste0('AUC = ', T2_EN_10_auc, ', ', T2_EN_10_ci_txt),
                    paste0('AUC = ', T2_EN_5_auc, ', ', T2_EN_5_ci_txt)),
           color = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2"), 
           size = 3, 
           x = 0.75, 
           y = c(0.30,0.25, 0.2, 0.15, 0.10, 0.05, 0.0))

#Boxplot
T2_EN_a_mod <- as.data.frame(T2_EN_a)
T2_EN_a_mod$Model <- "All proteins"

T2_EN_50_a_mod <- as.data.frame(T2_EN_50_a)
T2_EN_50_a_mod$Model <- "50 proteins"

T2_EN_40_a_mod <- as.data.frame(T2_EN_40_a)
T2_EN_40_a_mod$Model <- "40 proteins"

T2_EN_30_a_mod <- as.data.frame(T2_EN_30_a)
T2_EN_30_a_mod$Model <- "30 proteins"

T2_EN_20_a_mod <- as.data.frame(T2_EN_20_a)
T2_EN_20_a_mod$Model <- "20 proteins"

T2_EN_10_a_mod <- as.data.frame(T2_EN_10_a)
T2_EN_10_a_mod$Model <- "10 proteins"

T2_EN_5_a_mod <- as.data.frame(T2_EN_5_a)
T2_EN_5_a_mod$Model <- "5 proteins"


T2_EN_res <- rbind(T2_EN_a_mod, T2_EN_50_a_mod, T2_EN_40_a_mod, T2_EN_30_a_mod, T2_EN_20_a_mod, T2_EN_10_a_mod, T2_EN_5_a_mod)
T2_EN_res$Bwkat <- ifelse(T2_EN_res$true_labels == 1, "LGA", "AGA")
T2_EN_res$Model <- factor(T2_EN_res$Model, levels = c("All proteins",
                                                          "50 proteins",
                                                          "40 proteins",
                                                          "30 proteins",
                                                          "20 proteins",
                                                          "10 proteins",
                                                          "5 proteins"))

T2_EN_box <- ggplot(T2_EN_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "maroon4",
                              "50 proteins" = "mediumpurple4",
                              "40 proteins" = "mediumorchid2",
                              "30 proteins" = "mediumpurple2",
                              "20 proteins" = "maroon2",
                              "10 proteins" = "pink3",
                              "5 proteins" = "brown2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  #ggtitle("Visit 2\nWeek 20-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

# Compute Precision-Recall AUC
T2_EN_PR <- pr.curve(scores.class0 = T2_EN_a[,2], weights.class0 = T2_EN_a[,1], curve = TRUE)
T2_EN_PR_auc <- round(T2_EN_PR$auc.integral, digits = 2)
T2_EN_PR_df <- data.frame(Recall = T2_EN_PR$curve[, 1], Precision = T2_EN_PR$curve[, 2], Model = "All proteins")

T2_EN_50_PR <- pr.curve(scores.class0 = T2_EN_50_a[,2], weights.class0 = T2_EN_50_a[,1], curve = TRUE)
T2_EN_50_PR_auc <- round(T2_EN_50_PR$auc.integral, digits = 2)
T2_EN_50_PR_df <- data.frame(Recall = T2_EN_50_PR$curve[, 1], Precision = T2_EN_50_PR$curve[, 2], Model = "50 proteins")

T2_EN_40_PR <- pr.curve(scores.class0 = T2_EN_40_a[,2], weights.class0 = T2_EN_40_a[,1], curve = TRUE)
T2_EN_40_PR_auc <- round(T2_EN_40_PR$auc.integral, digits = 2)
T2_EN_40_PR_auc <- "0.20"
T2_EN_40_PR_df <- data.frame(Recall = T2_EN_40_PR$curve[, 1], Precision = T2_EN_40_PR$curve[, 2], Model = "40 proteins")

T2_EN_30_PR <- pr.curve(scores.class0 = T2_EN_30_a[,2], weights.class0 = T2_EN_30_a[,1], curve = TRUE)
T2_EN_30_PR_auc <- round(T2_EN_30_PR$auc.integral, digits = 2)
T2_EN_30_PR_df <- data.frame(Recall = T2_EN_30_PR$curve[, 1], Precision = T2_EN_30_PR$curve[, 2], Model = "30 proteins")

T2_EN_20_PR <- pr.curve(scores.class0 = T2_EN_20_a[,2], weights.class0 = T2_EN_20_a[,1], curve = TRUE)
T2_EN_20_PR_auc <- round(T2_EN_20_PR$auc.integral, digits = 2)
T2_EN_20_PR_df <- data.frame(Recall = T2_EN_20_PR$curve[, 1], Precision = T2_EN_20_PR$curve[, 2], Model = "20 proteins")

T2_EN_10_PR <- pr.curve(scores.class0 = T2_EN_10_a[,2], weights.class0 = T2_EN_10_a[,1], curve = TRUE)
T2_EN_10_PR_auc <- round(T2_EN_10_PR$auc.integral, digits = 2)
T2_EN_10_PR_auc <- "0.20"
T2_EN_10_PR_df <- data.frame(Recall = T2_EN_10_PR$curve[, 1], Precision = T2_EN_10_PR$curve[, 2], Model = "10 proteins")

T2_EN_5_PR <- pr.curve(scores.class0 = T2_EN_5_a[,2], weights.class0 = T2_EN_5_a[,1], curve = TRUE)
T2_EN_5_PR_auc <- round(T2_EN_5_PR$auc.integral, digits = 2)
T2_EN_5_PR_df <- data.frame(Recall = T2_EN_5_PR$curve[, 1], Precision = T2_EN_5_PR$curve[, 2], Model = "5 proteins")

T2_EN_pr_data <- bind_rows(T2_EN_PR_df, T2_EN_50_PR_df, T2_EN_40_PR_df, T2_EN_30_PR_df, T2_EN_20_PR_df, T2_EN_10_PR_df, T2_EN_5_PR_df)

T2_EN_pr_data$Model <- factor(T2_EN_pr_data$Model, levels = c("All proteins",
                                                              "50 proteins",
                                                              "40 proteins",
                                                              "30 proteins",
                                                              "20 proteins",
                                                              "10 proteins",
                                                              "5 proteins"))


T2_EN_PR <- ggplot(T2_EN_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "maroon4",
                               "50 proteins" = "mediumpurple4",
                               "40 proteins" = "mediumorchid2",
                               "30 proteins" = "mediumpurple2",
                               "20 proteins" = "maroon2",
                               "10 proteins" = "pink3",
                               "5 proteins" = "brown2")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('PR-AUC = ', T2_EN_PR_auc),
                    paste0('PR-AUC = ', T2_EN_50_PR_auc), 
                    paste0('PR-AUC = ', T2_EN_40_PR_auc), 
                    paste0('PR-AUC = ', T2_EN_30_PR_auc), 
                    paste0('PR-AUC = ', T2_EN_20_PR_auc),
                    paste0('PR-AUC = ', T2_EN_10_PR_auc),
                    paste0('PR-AUC = ', T2_EN_5_PR_auc)),
           color = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2"), 
           size = 3, 
           x = 0.9, 
           y = c(0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4)) 



# T3 Elastic net ----------------------------------------------------
#All proteins EN
load("../All_proteins/Predict_LGA_T3_EN_all_proteins.RData")

T3_EN_freq <- freq
T3_EN_pile <- pile

T3_EN_a=NULL
for (i in 1:length(T3_EN_pile)){
  T3_EN_a=rbind(T3_EN_a,T3_EN_pile[[i]]$outmat)
}
T3_EN_RC=roc(response=T3_EN_a[,1],predictor=T3_EN_a[,2],direction="<")
T3_EN_auc <- round(T3_EN_RC$auc, digits = 2)
T3_EN_ci <- ci(T3_EN_RC, of = "auc")
T3_EN_ci_txt <- "95% CI: 0.55-0.87"


#50 proteins EN
load("../50_proteins/Predict_LGA_T3_EN_50_proteins.RData")
T3_EN_50_freq <- freq
T3_EN_50_pile <- pile

T3_EN_50_a=NULL
for (i in 1:length(T3_EN_50_pile)){
  T3_EN_50_a=rbind(T3_EN_50_a,T3_EN_50_pile[[i]]$outmat)
}
T3_EN_50_RC=roc(response=T3_EN_50_a[,1],predictor=T3_EN_50_a[,2],direction="<")
T3_EN_50_auc <- round(T3_EN_50_RC$auc, digits = 2)
T3_EN_50_ci <- ci(T3_EN_50_RC, of = "auc")
T3_EN_50_ci_txt <- "95% CI: 0.53-0.90"


#40 proteins
load("../40_proteins/Predict_LGA_T3_EN_40_proteins.RData")
T3_EN_40_freq <- freq
T3_EN_40_pile <- pile

T3_EN_40_a=NULL
for (i in 1:length(T3_EN_40_pile)){
  T3_EN_40_a=rbind(T3_EN_40_a,T3_EN_40_pile[[i]]$outmat)
}
T3_EN_40_RC=roc(response=T3_EN_40_a[,1],predictor=T3_EN_40_a[,2],direction="<")
T3_EN_40_auc <- round(T3_EN_40_RC$auc, digits = 2)
T3_EN_40_ci <- ci(T3_EN_40_RC, of = "auc")
T3_EN_40_ci_txt <- "95% CI: 0.60-0.92"

#30 proteins EN
load("../30_proteins/Predict_LGA_T3_EN_30_proteins.RData")
T3_EN_30_freq <- freq
T3_EN_30_pile <- pile

T3_EN_30_a=NULL
for (i in 1:length(T3_EN_30_pile)){
  T3_EN_30_a=rbind(T3_EN_30_a,T3_EN_30_pile[[i]]$outmat)
}
T3_EN_30_RC=roc(response=T3_EN_30_a[,1],predictor=T3_EN_30_a[,2],direction="<")
T3_EN_30_auc <- round(T3_EN_30_RC$auc, digits = 2)
T3_EN_30_ci <- ci(T3_EN_30_RC, of = "auc")
T3_EN_30_ci_txt <- "95% CI: 0.55-0.91"


#20 proteins EN
load("../20_proteins/Predict_LGA_T3_EN_20_proteins.RData")
T3_EN_20_freq <- freq
T3_EN_20_pile <- pile

T3_EN_20_a=NULL
for (i in 1:length(T3_EN_20_pile)){
  T3_EN_20_a=rbind(T3_EN_20_a,T3_EN_20_pile[[i]]$outmat)
}
T3_EN_20_RC=roc(response=T3_EN_20_a[,1],predictor=T3_EN_20_a[,2],direction="<")

T3_EN_20_auc <- round(T3_EN_20_RC$auc, digits = 2)
T3_EN_20_ci <- ci(T3_EN_20_RC, of = "auc")
T3_EN_20_ci_txt <- "95% CI: 0.64-0.92"

#10 proteins EN
load("../10_proteins/Predict_LGA_T3_EN_10_proteins.RData")
T3_EN_10_freq <- freq
T3_EN_10_pile <- pile

T3_EN_10_a=NULL
for (i in 1:length(T3_EN_10_pile)){
  T3_EN_10_a=rbind(T3_EN_10_a,T3_EN_10_pile[[i]]$outmat)
}
T3_EN_10_RC=roc(response=T3_EN_10_a[,1],predictor=T3_EN_10_a[,2],direction="<")
T3_EN_10_auc <- round(T3_EN_10_RC$auc, digits = 2)
T3_EN_10_ci <- ci(T3_EN_10_RC, of = "auc")
T3_EN_10_ci_txt <- "95% CI: 0.52-0.92"

#5 proteins EN
load("../5_proteins/Predict_LGA_T3_EN_5_proteins.RData")
T3_EN_5_freq <- freq
T3_EN_5_pile <- pile

T3_EN_5_a=NULL
for (i in 1:length(T3_EN_5_pile)){
  T3_EN_5_a=rbind(T3_EN_5_a,T3_EN_5_pile[[i]]$outmat)
}
T3_EN_5_RC=roc(response=T3_EN_5_a[,1],predictor=T3_EN_5_a[,2],direction="<")
T3_EN_5_auc <- round(T3_EN_5_RC$auc, digits = 2)
T3_EN_5_ci <- ci(T3_EN_5_RC, of = "auc")
T3_EN_5_ci_txt <- "95% CI: 0.33-0.96"


#Plot
ROC_plot_EN_T3 <- ggroc(list("All proteins" = T3_EN_RC,
                             "50 proteins" = T3_EN_50_RC,
                             "40 proteins" = T3_EN_40_RC,
                             "30 proteins" = T3_EN_30_RC,
                             "20 proteins" = T3_EN_20_RC,
                             "10 proteins" = T3_EN_10_RC,
                             "5 proteins" = T3_EN_5_RC),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Elastic net\nVisit 3 (Week 28-34)") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T3_EN_auc, ', ', T3_EN_ci_txt),
                    paste0('AUC = ', T3_EN_50_auc, ', ', T3_EN_50_ci_txt), 
                    paste0('AUC = ', T3_EN_40_auc, ', ', T3_EN_40_ci_txt), 
                    paste0('AUC = ', T3_EN_30_auc, ', ', T3_EN_30_ci_txt), 
                    paste0('AUC = ', T3_EN_20_auc, ', ', T3_EN_20_ci_txt),
                    paste0('AUC = ', T3_EN_10_auc, ', ', T3_EN_10_ci_txt),
                    paste0('AUC = ', T3_EN_5_auc, ', ', T3_EN_5_ci_txt)),
           color = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2"), 
           size = 3, 
           x = 0.75, 
           y = c(0.30,0.25, 0.2, 0.15, 0.10, 0.05, 0.0))

#Boxplot
T3_EN_a_mod <- as.data.frame(T3_EN_a)
T3_EN_a_mod$Model <- "All proteins"

T3_EN_50_a_mod <- as.data.frame(T3_EN_50_a)
T3_EN_50_a_mod$Model <- "50 proteins"

T3_EN_40_a_mod <- as.data.frame(T3_EN_40_a)
T3_EN_40_a_mod$Model <- "40 proteins"

T3_EN_30_a_mod <- as.data.frame(T3_EN_30_a)
T3_EN_30_a_mod$Model <- "30 proteins"

T3_EN_20_a_mod <- as.data.frame(T3_EN_20_a)
T3_EN_20_a_mod$Model <- "20 proteins"

T3_EN_10_a_mod <- as.data.frame(T3_EN_10_a)
T3_EN_10_a_mod$Model <- "10 proteins"

T3_EN_5_a_mod <- as.data.frame(T3_EN_5_a)
T3_EN_5_a_mod$Model <- "5 proteins"


T3_EN_res <- rbind(T3_EN_a_mod, T3_EN_50_a_mod, T3_EN_40_a_mod, T3_EN_30_a_mod, T3_EN_20_a_mod, T3_EN_10_a_mod, T3_EN_5_a_mod)
T3_EN_res$Bwkat <- ifelse(T3_EN_res$true_labels == 1, "LGA", "AGA")
T3_EN_res$Model <- factor(T3_EN_res$Model, levels = c("All proteins",
                                                      "50 proteins",
                                                      "40 proteins",
                                                      "30 proteins",
                                                      "20 proteins",
                                                      "10 proteins",
                                                      "5 proteins"))

T3_EN_box <- ggplot(T3_EN_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "maroon4",
                              "50 proteins" = "mediumpurple4",
                              "40 proteins" = "mediumorchid2",
                              "30 proteins" = "mediumpurple2",
                              "20 proteins" = "maroon2",
                              "10 proteins" = "pink3",
                              "5 proteins" = "brown2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  #ggtitle("Visit 3\nWeek 28-34") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

# Compute Precision-Recall AUC
T3_EN_PR <- pr.curve(scores.class0 = T3_EN_a[,2], weights.class0 = T3_EN_a[,1], curve = TRUE)
T3_EN_PR_auc <- round(T3_EN_PR$auc.integral, digits = 2)
T3_EN_PR_df <- data.frame(Recall = T3_EN_PR$curve[, 1], Precision = T3_EN_PR$curve[, 2], Model = "All proteins")

T3_EN_50_PR <- pr.curve(scores.class0 = T3_EN_50_a[,2], weights.class0 = T3_EN_50_a[,1], curve = TRUE)
T3_EN_50_PR_auc <- round(T3_EN_50_PR$auc.integral, digits = 2)
T3_EN_50_PR_df <- data.frame(Recall = T3_EN_50_PR$curve[, 1], Precision = T3_EN_50_PR$curve[, 2], Model = "50 proteins")

T3_EN_40_PR <- pr.curve(scores.class0 = T3_EN_40_a[,2], weights.class0 = T3_EN_40_a[,1], curve = TRUE)
T3_EN_40_PR_auc <- round(T3_EN_40_PR$auc.integral, digits = 2)
T3_EN_40_PR_df <- data.frame(Recall = T3_EN_40_PR$curve[, 1], Precision = T3_EN_40_PR$curve[, 2], Model = "40 proteins")

T3_EN_30_PR <- pr.curve(scores.class0 = T3_EN_30_a[,2], weights.class0 = T3_EN_30_a[,1], curve = TRUE)
T3_EN_30_PR_auc <- round(T3_EN_30_PR$auc.integral, digits = 2)
T3_EN_30_PR_df <- data.frame(Recall = T3_EN_30_PR$curve[, 1], Precision = T3_EN_30_PR$curve[, 2], Model = "30 proteins")

T3_EN_20_PR <- pr.curve(scores.class0 = T3_EN_20_a[,2], weights.class0 = T3_EN_20_a[,1], curve = TRUE)
T3_EN_20_PR_auc <- round(T3_EN_20_PR$auc.integral, digits = 2)
T3_EN_20_PR_df <- data.frame(Recall = T3_EN_20_PR$curve[, 1], Precision = T3_EN_20_PR$curve[, 2], Model = "20 proteins")

T3_EN_10_PR <- pr.curve(scores.class0 = T3_EN_10_a[,2], weights.class0 = T3_EN_10_a[,1], curve = TRUE)
T3_EN_10_PR_auc <- round(T3_EN_10_PR$auc.integral, digits = 2)
T3_EN_10_PR_df <- data.frame(Recall = T3_EN_10_PR$curve[, 1], Precision = T3_EN_10_PR$curve[, 2], Model = "10 proteins")

T3_EN_5_PR <- pr.curve(scores.class0 = T3_EN_5_a[,2], weights.class0 = T3_EN_5_a[,1], curve = TRUE)
T3_EN_5_PR_auc <- round(T3_EN_5_PR$auc.integral, digits = 2)
T3_EN_5_PR_df <- data.frame(Recall = T3_EN_5_PR$curve[, 1], Precision = T3_EN_5_PR$curve[, 2], Model = "5 proteins")

T3_EN_pr_data <- bind_rows(T3_EN_PR_df, T3_EN_50_PR_df, T3_EN_40_PR_df, T3_EN_30_PR_df, T3_EN_20_PR_df, T3_EN_10_PR_df, T3_EN_5_PR_df)

T3_EN_pr_data$Model <- factor(T3_EN_pr_data$Model, levels = c("All proteins",
                                                              "50 proteins",
                                                              "40 proteins",
                                                              "30 proteins",
                                                              "20 proteins",
                                                              "10 proteins",
                                                              "5 proteins"))


T3_EN_PR <- ggplot(T3_EN_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "maroon4",
                               "50 proteins" = "mediumpurple4",
                               "40 proteins" = "mediumorchid2",
                               "30 proteins" = "mediumpurple2",
                               "20 proteins" = "maroon2",
                               "10 proteins" = "pink3",
                               "5 proteins" = "brown2")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('PR-AUC = ', T3_EN_PR_auc),
                    paste0('PR-AUC = ', T3_EN_50_PR_auc), 
                    paste0('PR-AUC = ', T3_EN_40_PR_auc), 
                    paste0('PR-AUC = ', T3_EN_30_PR_auc), 
                    paste0('PR-AUC = ', T3_EN_20_PR_auc),
                    paste0('PR-AUC = ', T3_EN_10_PR_auc),
                    paste0('PR-AUC = ', T3_EN_5_PR_auc)),
           color = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2"), 
           size = 3, 
           x = 0.9, 
           y = c(0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4)) 



# T1 Random Forest ----------------------------------------------------
#All proteins RF
load("../All_proteins/Predict_LGA_T1_RF_all_proteins.RData")

T1_RF_freq <- freq
T1_RF_pile <- pile

T1_RF_a=NULL
for (i in 1:length(T1_RF_pile)){
  T1_RF_a=rbind(T1_RF_a,T1_RF_pile[[i]]$outmat)
}
T1_RF_RC=roc(response=T1_RF_a[,1],predictor=T1_RF_a[,2],direction="<")
T1_RF_auc <- round(T1_RF_RC$auc, digits = 2)
T1_RF_ci <- ci(T1_RF_RC, of = "auc")
T1_RF_ci_txt <- "95% CI: 0.40-0.90"


#50 proteins RF
load("../50_proteins/Predict_LGA_T1_RF_50_proteins.RData")
T1_RF_50_freq <- freq
T1_RF_50_pile <- pile

T1_RF_50_a=NULL
for (i in 1:length(T1_RF_50_pile)){
  T1_RF_50_a=rbind(T1_RF_50_a,T1_RF_50_pile[[i]]$outmat)
}
T1_RF_50_RC=roc(response=T1_RF_50_a[,1],predictor=T1_RF_50_a[,2],direction="<")
T1_RF_50_auc <- round(T1_RF_50_RC$auc, digits = 2)
T1_RF_50_ci <- ci(T1_RF_50_RC, of = "auc")
T1_RF_50_ci_txt <- "95% CI: 0.61-0.92"


#40 proteins
load("../40_proteins/Predict_LGA_T1_RF_40_proteins.RData")
T1_RF_40_freq <- freq
T1_RF_40_pile <- pile

T1_RF_40_a=NULL
for (i in 1:length(T1_RF_40_pile)){
  T1_RF_40_a=rbind(T1_RF_40_a,T1_RF_40_pile[[i]]$outmat)
}
T1_RF_40_RC=roc(response=T1_RF_40_a[,1],predictor=T1_RF_40_a[,2],direction="<")
T1_RF_40_auc <- round(T1_RF_40_RC$auc, digits = 2)
T1_RF_40_ci <- ci(T1_RF_40_RC, of = "auc")
T1_RF_40_ci_txt <- "95% CI: 0.61-0.97"

#30 proteins RF
load("../30_proteins/Predict_LGA_T1_RF_30_proteins.RData")
T1_RF_30_freq <- freq
T1_RF_30_pile <- pile

T1_RF_30_a=NULL
for (i in 1:length(T1_RF_30_pile)){
  T1_RF_30_a=rbind(T1_RF_30_a,T1_RF_30_pile[[i]]$outmat)
}
T1_RF_30_RC=roc(response=T1_RF_30_a[,1],predictor=T1_RF_30_a[,2],direction="<")

T1_RF_30_auc <- round(T1_RF_30_RC$auc, digits = 2)
T1_RF_30_ci <- ci(T1_RF_30_RC, of = "auc")
T1_RF_30_ci_txt <- "95% CI: 0.67-1.00"


#20 proteins RF
load("../20_proteins/Predict_LGA_T1_RF_20_proteins.RData")
T1_RF_20_freq <- freq
T1_RF_20_pile <- pile

T1_RF_20_a=NULL
for (i in 1:length(T1_RF_20_pile)){
  T1_RF_20_a=rbind(T1_RF_20_a,T1_RF_20_pile[[i]]$outmat)
}
T1_RF_20_RC=roc(response=T1_RF_20_a[,1],predictor=T1_RF_20_a[,2],direction="<")
T1_RF_20_auc <- round(T1_RF_20_RC$auc, digits = 2)
T1_RF_20_ci <- ci(T1_RF_20_RC, of = "auc")
T1_RF_20_ci_txt <- "95% CI: 0.72-0.98"

#10 proteins RF
load("../10_proteins/Predict_LGA_T1_RF_10_proteins.RData")
T1_RF_10_freq <- freq
T1_RF_10_pile <- pile

T1_RF_10_a=NULL
for (i in 1:length(T1_RF_10_pile)){
  T1_RF_10_a=rbind(T1_RF_10_a,T1_RF_10_pile[[i]]$outmat)
}
T1_RF_10_RC=roc(response=T1_RF_10_a[,1],predictor=T1_RF_10_a[,2],direction="<")
T1_RF_10_auc <- round(T1_RF_10_RC$auc, digits = 2)
T1_RF_10_ci <- ci(T1_RF_10_RC, of = "auc")
T1_RF_10_ci_txt <- "95% CI: 0.83-1.00"

#5 proteins RF
load("../5_proteins/Predict_LGA_T1_RF_5_proteins.RData")
T1_RF_5_freq <- freq
T1_RF_5_pile <- pile

T1_RF_5_a=NULL
for (i in 1:length(T1_RF_5_pile)){
  T1_RF_5_a=rbind(T1_RF_5_a,T1_RF_5_pile[[i]]$outmat)
}
T1_RF_5_RC=roc(response=T1_RF_5_a[,1],predictor=T1_RF_5_a[,2],direction="<")
T1_RF_5_auc <- round(T1_RF_5_RC$auc, digits = 2)
T1_RF_5_ci <- ci(T1_RF_5_RC, of = "auc")
T1_RF_5_ci_txt <- "95% CI: 0.52-1.00"


#Plot
ROC_plot_RF_T1 <- ggroc(list("All proteins" = T1_RF_RC,
                             "50 proteins" = T1_RF_50_RC,
                             "40 proteins" = T1_RF_40_RC,
                             "30 proteins" = T1_RF_30_RC,
                             "20 proteins" = T1_RF_20_RC,
                             "10 proteins" = T1_RF_10_RC,
                             "5 proteins" = T1_RF_5_RC),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("navy", "slateblue1","royalblue3", "mediumblue", "dodgerblue","deepskyblue3", "turquoise2")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Random forest\nVisit (Week 12-19)") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T1_RF_auc, ', ', T1_RF_ci_txt),
                    paste0('AUC = ', T1_RF_50_auc, ', ', T1_RF_50_ci_txt), 
                    paste0('AUC = ', T1_RF_40_auc, ', ', T1_RF_40_ci_txt), 
                    paste0('AUC = ', T1_RF_30_auc, ', ', T1_RF_30_ci_txt), 
                    paste0('AUC = ', T1_RF_20_auc, ', ', T1_RF_20_ci_txt),
                    paste0('AUC = ', T1_RF_10_auc, ', ', T1_RF_10_ci_txt),
                    paste0('AUC = ', T1_RF_5_auc, ', ', T1_RF_5_ci_txt)),
           color = c("navy", "slateblue1","royalblue3","mediumblue", "dodgerblue","deepskyblue3", "turquoise2"), 
           size = 3, 
           x = 0.75, 
           y = c(0.30,0.25, 0.2, 0.15, 0.10, 0.05, 0.0))

#Boxplot
T1_RF_a_mod <- as.data.frame(T1_RF_a)
T1_RF_a_mod$Model <- "All proteins"

T1_RF_50_a_mod <- as.data.frame(T1_RF_50_a)
T1_RF_50_a_mod$Model <- "50 proteins"

T1_RF_40_a_mod <- as.data.frame(T1_RF_40_a)
T1_RF_40_a_mod$Model <- "40 proteins"

T1_RF_30_a_mod <- as.data.frame(T1_RF_30_a)
T1_RF_30_a_mod$Model <- "30 proteins"

T1_RF_20_a_mod <- as.data.frame(T1_RF_20_a)
T1_RF_20_a_mod$Model <- "20 proteins"

T1_RF_10_a_mod <- as.data.frame(T1_RF_10_a)
T1_RF_10_a_mod$Model <- "10 proteins"

T1_RF_5_a_mod <- as.data.frame(T1_RF_5_a)
T1_RF_5_a_mod$Model <- "5 proteins"


T1_RF_res <- rbind(T1_RF_a_mod, T1_RF_50_a_mod, T1_RF_40_a_mod, T1_RF_30_a_mod, T1_RF_20_a_mod, T1_RF_10_a_mod, T1_RF_5_a_mod)
T1_RF_res$Bwkat <- ifelse(T1_RF_res$true_labels == 1, "LGA", "AGA")
T1_RF_res$Model <- factor(T1_RF_res$Model, levels = c("All proteins",
                                                      "50 proteins",
                                                      "40 proteins",
                                                      "30 proteins",
                                                      "20 proteins",
                                                      "10 proteins",
                                                      "5 proteins"))

T1_RF_box <- ggplot(T1_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "navy",
                              "50 proteins" = "slateblue1",
                              "40 proteins" = "royalblue3",
                              "30 proteins" = "mediumblue",
                              "20 proteins" = "dodgerblue",
                              "10 proteins" = "deepskyblue3",
                              "5 proteins" = "turquoise2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  #ggtitle("Visit 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

# Compute Precision-Recall AUC
T1_RF_PR <- pr.curve(scores.class0 = T1_RF_a[,2], weights.class0 = T1_RF_a[,1], curve = TRUE)
T1_RF_PR_auc <- round(T1_RF_PR$auc.integral, digits = 2)
T1_RF_PR_df <- data.frame(Recall = T1_RF_PR$curve[, 1], Precision = T1_RF_PR$curve[, 2], Model = "All proteins")

T1_RF_50_PR <- pr.curve(scores.class0 = T1_RF_50_a[,2], weights.class0 = T1_RF_50_a[,1], curve = TRUE)
T1_RF_50_PR_auc <- round(T1_RF_50_PR$auc.integral, digits = 2)
T1_RF_50_PR_df <- data.frame(Recall = T1_RF_50_PR$curve[, 1], Precision = T1_RF_50_PR$curve[, 2], Model = "50 proteins")

T1_RF_40_PR <- pr.curve(scores.class0 = T1_RF_40_a[,2], weights.class0 = T1_RF_40_a[,1], curve = TRUE)
T1_RF_40_PR_auc <- round(T1_RF_40_PR$auc.integral, digits = 2)
T1_RF_40_PR_df <- data.frame(Recall = T1_RF_40_PR$curve[, 1], Precision = T1_RF_40_PR$curve[, 2], Model = "40 proteins")

T1_RF_30_PR <- pr.curve(scores.class0 = T1_RF_30_a[,2], weights.class0 = T1_RF_30_a[,1], curve = TRUE)
T1_RF_30_PR_auc <- round(T1_RF_30_PR$auc.integral, digits = 2)
T1_RF_30_PR_df <- data.frame(Recall = T1_RF_30_PR$curve[, 1], Precision = T1_RF_30_PR$curve[, 2], Model = "30 proteins")

T1_RF_20_PR <- pr.curve(scores.class0 = T1_RF_20_a[,2], weights.class0 = T1_RF_20_a[,1], curve = TRUE)
T1_RF_20_PR_auc <- round(T1_RF_20_PR$auc.integral, digits = 2)
T1_RF_20_PR_df <- data.frame(Recall = T1_RF_20_PR$curve[, 1], Precision = T1_RF_20_PR$curve[, 2], Model = "20 proteins")

T1_RF_10_PR <- pr.curve(scores.class0 = T1_RF_10_a[,2], weights.class0 = T1_RF_10_a[,1], curve = TRUE)
T1_RF_10_PR_auc <- round(T1_RF_10_PR$auc.integral, digits = 2)
T1_RF_10_PR_df <- data.frame(Recall = T1_RF_10_PR$curve[, 1], Precision = T1_RF_10_PR$curve[, 2], Model = "10 proteins")

T1_RF_5_PR <- pr.curve(scores.class0 = T1_RF_5_a[,2], weights.class0 = T1_RF_5_a[,1], curve = TRUE)
T1_RF_5_PR_auc <- round(T1_RF_5_PR$auc.integral, digits = 2)
T1_RF_5_PR_df <- data.frame(Recall = T1_RF_5_PR$curve[, 1], Precision = T1_RF_5_PR$curve[, 2], Model = "5 proteins")

T1_RF_pr_data <- bind_rows(T1_RF_PR_df, T1_RF_50_PR_df, T1_RF_40_PR_df, T1_RF_30_PR_df, T1_RF_20_PR_df, T1_RF_10_PR_df, T1_RF_5_PR_df)

T1_RF_pr_data$Model <- factor(T1_RF_pr_data$Model, levels = c("All proteins",
                                                      "50 proteins",
                                                      "40 proteins",
                                                      "30 proteins",
                                                      "20 proteins",
                                                      "10 proteins",
                                                      "5 proteins"))


T1_RF_PR <- ggplot(T1_RF_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "navy",
                               "50 proteins" = "slateblue1",
                               "40 proteins" = "royalblue3",
                               "30 proteins" = "mediumblue",
                               "20 proteins" = "dodgerblue",
                               "10 proteins" = "deepskyblue3",
                               "5 proteins" = "turquoise2")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('PR-AUC = ', T1_RF_PR_auc),
                    paste0('PR-AUC = ', T1_RF_50_PR_auc), 
                    paste0('PR-AUC = ', T1_RF_40_PR_auc), 
                    paste0('PR-AUC = ', T1_RF_30_PR_auc), 
                    paste0('PR-AUC = ', T1_RF_20_PR_auc),
                    paste0('PR-AUC = ', T1_RF_10_PR_auc),
                    paste0('PR-AUC = ', T1_RF_5_PR_auc)),
           color = c("navy", "slateblue1","royalblue3","mediumblue", "dodgerblue","deepskyblue3", "turquoise2"), 
           size = 3, 
           x = 0.9, 
           y = c(0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4)) 



# T2 Random Forest ----------------------------------------------------
#All proteins RF
load("../All_proteins/Predict_LGA_T2_RF_all_proteins.RData")

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


#50 proteins RF
load("../50_proteins/Predict_LGA_T2_RF_50_proteins.RData")
T2_RF_50_freq <- freq
T2_RF_50_pile <- pile

T2_RF_50_a=NULL
for (i in 1:length(T2_RF_50_pile)){
  T2_RF_50_a=rbind(T2_RF_50_a,T2_RF_50_pile[[i]]$outmat)
}
T2_RF_50_RC=roc(response=T2_RF_50_a[,1],predictor=T2_RF_50_a[,2],direction="<")

T2_RF_50_auc <- round(T2_RF_50_RC$auc, digits = 2)
T2_RF_50_ci <- ci(T2_RF_50_RC, of = "auc")
T2_RF_50_ci_txt <- "95% CI: 0.69-0.95"


#40 proteins
load("../40_proteins/Predict_LGA_T2_RF_40_proteins.RData")
T2_RF_40_freq <- freq
T2_RF_40_pile <- pile

T2_RF_40_a=NULL
for (i in 1:length(T2_RF_40_pile)){
  T2_RF_40_a=rbind(T2_RF_40_a,T2_RF_40_pile[[i]]$outmat)
}
T2_RF_40_RC=roc(response=T2_RF_40_a[,1],predictor=T2_RF_40_a[,2],direction="<")
T2_RF_40_auc <- round(T2_RF_40_RC$auc, digits = 2)
T2_RF_40_ci <- ci(T2_RF_40_RC, of = "auc")
T2_RF_40_ci_txt <- "95% CI: 0.62-0.95"

#30 proteins RF
load("../30_proteins/Predict_LGA_T2_RF_30_proteins.RData")
T2_RF_30_freq <- freq
T2_RF_30_pile <- pile

T2_RF_30_a=NULL
for (i in 1:length(T2_RF_30_pile)){
  T2_RF_30_a=rbind(T2_RF_30_a,T2_RF_30_pile[[i]]$outmat)
}
T2_RF_30_RC=roc(response=T2_RF_30_a[,1],predictor=T2_RF_30_a[,2],direction="<")
T2_RF_30_auc <- round(T2_RF_30_RC$auc, digits = 2)
T2_RF_30_ci <- ci(T2_RF_30_RC, of = "auc")
T2_RF_30_ci_txt <- "95% CI: 0.67-1.00"


#20 proteins RF
load("../20_proteins/Predict_LGA_T2_RF_20_proteins.RData")
T2_RF_20_freq <- freq
T2_RF_20_pile <- pile

T2_RF_20_a=NULL
for (i in 1:length(T2_RF_20_pile)){
  T2_RF_20_a=rbind(T2_RF_20_a,T2_RF_20_pile[[i]]$outmat)
}
T2_RF_20_RC=roc(response=T2_RF_20_a[,1],predictor=T2_RF_20_a[,2],direction="<")
T2_RF_20_auc <- round(T2_RF_20_RC$auc, digits = 2)
T2_RF_20_ci <- ci(T2_RF_20_RC, of = "auc")
T2_RF_20_ci_txt <- "95% CI: 0.73-1.00"

#10 proteins RF
load("../10_proteins/Predict_LGA_T2_RF_10_proteins.RData")
T2_RF_10_freq <- freq
T2_RF_10_pile <- pile

T2_RF_10_a=NULL
for (i in 1:length(T2_RF_10_pile)){
  T2_RF_10_a=rbind(T2_RF_10_a,T2_RF_10_pile[[i]]$outmat)
}
T2_RF_10_RC=roc(response=T2_RF_10_a[,1],predictor=T2_RF_10_a[,2],direction="<")

T2_RF_10_auc <- round(T2_RF_10_RC$auc, digits = 2)
T2_RF_10_ci <- ci(T2_RF_10_RC, of = "auc")
T2_RF_10_ci_txt <- "95% CI: 0.34-0.88"

#5 proteins RF
load("../5_proteins/Predict_LGA_T2_RF_5_proteins.RData")
T2_RF_5_freq <- freq
T2_RF_5_pile <- pile

T2_RF_5_a=NULL
for (i in 1:length(T2_RF_5_pile)){
  T2_RF_5_a=rbind(T2_RF_5_a,T2_RF_5_pile[[i]]$outmat)
}
T2_RF_5_RC=roc(response=T2_RF_5_a[,1],predictor=T2_RF_5_a[,2],direction="<")

T2_RF_5_auc <- round(T2_RF_5_RC$auc, digits = 2)
T2_RF_5_ci <- ci(T2_RF_5_RC, of = "auc")
T2_RF_5_ci_txt <- "95% CI: 0.16-0.47"


#Plot
ROC_plot_RF_T2 <- ggroc(list("All proteins" = T2_RF_RC,
                             "50 proteins" = T2_RF_50_RC,
                             "40 proteins" = T2_RF_40_RC,
                             "30 proteins" = T2_RF_30_RC,
                             "20 proteins" = T2_RF_20_RC,
                             "10 proteins" = T2_RF_10_RC,
                             "5 proteins" = T2_RF_5_RC),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("navy", "slateblue1","royalblue3", "mediumblue", "dodgerblue","deepskyblue3", "turquoise2")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Random forest\nVisit 2 (Week 21-27)") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T2_RF_auc, ', ', T2_RF_ci_txt),
                    paste0('AUC = ', T2_RF_50_auc, ', ', T2_RF_50_ci_txt), 
                    paste0('AUC = ', T2_RF_40_auc, ', ', T2_RF_40_ci_txt), 
                    paste0('AUC = ', T2_RF_30_auc, ', ', T2_RF_30_ci_txt), 
                    paste0('AUC = ', T2_RF_20_auc, ', ', T2_RF_20_ci_txt),
                    paste0('AUC = ', T2_RF_10_auc, ', ', T2_RF_10_ci_txt),
                    paste0('AUC = ', T2_RF_5_auc, ', ', T2_RF_5_ci_txt)),
           color = c("navy", "slateblue1","royalblue3","mediumblue", "dodgerblue","deepskyblue3", "turquoise2"), 
           size = 3, 
           x = 0.75, 
           y = c(0.30,0.25, 0.2, 0.15, 0.10, 0.05, 0.0))


#Boxplot
T2_RF_a_mod <- as.data.frame(T2_RF_a)
T2_RF_a_mod$Model <- "All proteins"

T2_RF_50_a_mod <- as.data.frame(T2_RF_50_a)
T2_RF_50_a_mod$Model <- "50 proteins"

T2_RF_40_a_mod <- as.data.frame(T2_RF_40_a)
T2_RF_40_a_mod$Model <- "40 proteins"

T2_RF_30_a_mod <- as.data.frame(T2_RF_30_a)
T2_RF_30_a_mod$Model <- "30 proteins"

T2_RF_20_a_mod <- as.data.frame(T2_RF_20_a)
T2_RF_20_a_mod$Model <- "20 proteins"

T2_RF_10_a_mod <- as.data.frame(T2_RF_10_a)
T2_RF_10_a_mod$Model <- "10 proteins"

T2_RF_5_a_mod <- as.data.frame(T2_RF_5_a)
T2_RF_5_a_mod$Model <- "5 proteins"


T2_RF_res <- rbind(T2_RF_a_mod, T2_RF_50_a_mod, T2_RF_40_a_mod, T2_RF_30_a_mod, T2_RF_20_a_mod, T2_RF_10_a_mod, T2_RF_5_a_mod)
T2_RF_res$Bwkat <- ifelse(T2_RF_res$true_labels == 1, "LGA", "AGA")
T2_RF_res$Model <- factor(T2_RF_res$Model, levels = c("All proteins",
                                                      "50 proteins",
                                                      "40 proteins",
                                                      "30 proteins",
                                                      "20 proteins",
                                                      "10 proteins",
                                                      "5 proteins"))

T2_RF_box <- ggplot(T2_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "navy",
                              "50 proteins" = "slateblue1",
                              "40 proteins" = "royalblue3",
                              "30 proteins" = "mediumblue",
                              "20 proteins" = "dodgerblue",
                              "10 proteins" = "deepskyblue3",
                              "5 proteins" = "turquoise2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  #ggtitle("Visit 2\nWeek 20-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

# Compute Precision-Recall AUC
T2_RF_PR <- pr.curve(scores.class0 = T2_RF_a[,2], weights.class0 = T2_RF_a[,1], curve = TRUE)
T2_RF_PR_auc <- round(T2_RF_PR$auc.integral, digits = 2)
T2_RF_PR_df <- data.frame(Recall = T2_RF_PR$curve[, 1], Precision = T2_RF_PR$curve[, 2], Model = "All proteins")

T2_RF_50_PR <- pr.curve(scores.class0 = T2_RF_50_a[,2], weights.class0 = T2_RF_50_a[,1], curve = TRUE)
T2_RF_50_PR_auc <- round(T2_RF_50_PR$auc.integral, digits = 2)
T2_RF_50_PR_df <- data.frame(Recall = T2_RF_50_PR$curve[, 1], Precision = T2_RF_50_PR$curve[, 2], Model = "50 proteins")

T2_RF_40_PR <- pr.curve(scores.class0 = T2_RF_40_a[,2], weights.class0 = T2_RF_40_a[,1], curve = TRUE)
T2_RF_40_PR_auc <- round(T2_RF_40_PR$auc.integral, digits = 2)
T2_RF_40_PR_df <- data.frame(Recall = T2_RF_40_PR$curve[, 1], Precision = T2_RF_40_PR$curve[, 2], Model = "40 proteins")

T2_RF_30_PR <- pr.curve(scores.class0 = T2_RF_30_a[,2], weights.class0 = T2_RF_30_a[,1], curve = TRUE)
T2_RF_30_PR_auc <- round(T2_RF_30_PR$auc.integral, digits = 2)
T2_RF_30_PR_df <- data.frame(Recall = T2_RF_30_PR$curve[, 1], Precision = T2_RF_30_PR$curve[, 2], Model = "30 proteins")

T2_RF_20_PR <- pr.curve(scores.class0 = T2_RF_20_a[,2], weights.class0 = T2_RF_20_a[,1], curve = TRUE)
T2_RF_20_PR_auc <- round(T2_RF_20_PR$auc.integral, digits = 2)
T2_RF_20_PR_auc <- "0.30"
T2_RF_20_PR_df <- data.frame(Recall = T2_RF_20_PR$curve[, 1], Precision = T2_RF_20_PR$curve[, 2], Model = "20 proteins")

T2_RF_10_PR <- pr.curve(scores.class0 = T2_RF_10_a[,2], weights.class0 = T2_RF_10_a[,1], curve = TRUE)
T2_RF_10_PR_auc <- round(T2_RF_10_PR$auc.integral, digits = 2)
T2_RF_10_PR_df <- data.frame(Recall = T2_RF_10_PR$curve[, 1], Precision = T2_RF_10_PR$curve[, 2], Model = "10 proteins")

T2_RF_5_PR <- pr.curve(scores.class0 = T2_RF_5_a[,2], weights.class0 = T2_RF_5_a[,1], curve = TRUE)
T2_RF_5_PR_auc <- round(T2_RF_5_PR$auc.integral, digits = 2)
T2_RF_5_PR_df <- data.frame(Recall = T2_RF_5_PR$curve[, 1], Precision = T2_RF_5_PR$curve[, 2], Model = "5 proteins")

T2_RF_pr_data <- bind_rows(T2_RF_PR_df, T2_RF_50_PR_df, T2_RF_40_PR_df, T2_RF_30_PR_df, T2_RF_20_PR_df, T2_RF_10_PR_df, T2_RF_5_PR_df)

T2_RF_pr_data$Model <- factor(T2_RF_pr_data$Model, levels = c("All proteins",
                                                              "50 proteins",
                                                              "40 proteins",
                                                              "30 proteins",
                                                              "20 proteins",
                                                              "10 proteins",
                                                              "5 proteins"))


T2_RF_PR <- ggplot(T2_RF_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "navy",
                               "50 proteins" = "slateblue1",
                               "40 proteins" = "royalblue3",
                               "30 proteins" = "mediumblue",
                               "20 proteins" = "dodgerblue",
                               "10 proteins" = "deepskyblue3",
                               "5 proteins" = "turquoise2")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('PR-AUC = ', T2_RF_PR_auc),
                    paste0('PR-AUC = ', T2_RF_50_PR_auc), 
                    paste0('PR-AUC = ', T2_RF_40_PR_auc), 
                    paste0('PR-AUC = ', T2_RF_30_PR_auc), 
                    paste0('PR-AUC = ', T2_RF_20_PR_auc),
                    paste0('PR-AUC = ', T2_RF_10_PR_auc),
                    paste0('PR-AUC = ', T2_RF_5_PR_auc)),
           color = c("navy", "slateblue1","royalblue3","mediumblue", "dodgerblue","deepskyblue3", "turquoise2"), 
           size = 3, 
           x = 0.9, 
           y = c(0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4)) 



# T3 Random forest ----------------------------------------------------
#All proteins RF
load("../All_proteins/Predict_LGA_T3_RF_all_proteins.RData")
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


#50 proteins RF
load("../50_proteins/Predict_LGA_T3_RF_50_proteins.RData")
T3_RF_50_freq <- freq
T3_RF_50_pile <- pile

T3_RF_50_a=NULL
for (i in 1:length(T3_RF_50_pile)){
  T3_RF_50_a=rbind(T3_RF_50_a,T3_RF_50_pile[[i]]$outmat)
}
T3_RF_50_RC=roc(response=T3_RF_50_a[,1],predictor=T3_RF_50_a[,2],direction="<")
T3_RF_50_auc <- round(T3_RF_50_RC$auc, digits = 2)
T3_RF_50_ci <- ci(T3_RF_50_RC, of = "auc")
T3_RF_50_ci_txt <- "95% CI: 0.69-0.93"


#40 proteins
load("../40_proteins/Predict_LGA_T3_RF_40_proteins.RData")
T3_RF_40_freq <- freq
T3_RF_40_pile <- pile

T3_RF_40_a=NULL
for (i in 1:length(T3_RF_40_pile)){
  T3_RF_40_a=rbind(T3_RF_40_a,T3_RF_40_pile[[i]]$outmat)
}
T3_RF_40_RC=roc(response=T3_RF_40_a[,1],predictor=T3_RF_40_a[,2],direction="<")

T3_RF_40_auc <- round(T3_RF_40_RC$auc, digits = 2)
T3_RF_40_ci <- ci(T3_RF_40_RC, of = "auc")
T3_RF_40_ci_txt <- "95% CI: 0.73-0.95"

#30 proteins RF
load("../30_proteins/Predict_LGA_T3_RF_30_proteins.RData")
T3_RF_30_freq <- freq
T3_RF_30_pile <- pile

T3_RF_30_a=NULL
for (i in 1:length(T3_RF_30_pile)){
  T3_RF_30_a=rbind(T3_RF_30_a,T3_RF_30_pile[[i]]$outmat)
}
T3_RF_30_RC=roc(response=T3_RF_30_a[,1],predictor=T3_RF_30_a[,2],direction="<")
T3_RF_30_auc <- round(T3_RF_30_RC$auc, digits = 2)
T3_RF_30_ci <- ci(T3_RF_30_RC, of = "auc")
T3_RF_30_ci_txt <- "95% CI: 0.70-0.95"


#20 proteins RF
load("../20_proteins/Predict_LGA_T3_RF_20_proteins.RData")
T3_RF_20_freq <- freq
T3_RF_20_pile <- pile
T3_RF_20_a=NULL
for (i in 1:length(T3_RF_20_pile)){
  T3_RF_20_a=rbind(T3_RF_20_a,T3_RF_20_pile[[i]]$outmat)
}
T3_RF_20_RC=roc(response=T3_RF_20_a[,1],predictor=T3_RF_20_a[,2],direction="<")

T3_RF_20_auc <- round(T3_RF_20_RC$auc, digits = 2)
T3_RF_20_auc <- "0.80"
T3_RF_20_ci <- ci(T3_RF_20_RC, of = "auc")
T3_RF_20_ci_txt <- "95% CI: 0.59-1.00"

#10 proteins RF
load("../10_proteins/Predict_LGA_T3_RF_10_proteins.RData")
T3_RF_10_freq <- freq
T3_RF_10_pile <- pile

T3_RF_10_a=NULL
for (i in 1:length(T3_RF_10_pile)){
  T3_RF_10_a=rbind(T3_RF_10_a,T3_RF_10_pile[[i]]$outmat)
}
T3_RF_10_RC=roc(response=T3_RF_10_a[,1],predictor=T3_RF_10_a[,2],direction="<")

T3_RF_10_auc <- round(T3_RF_10_RC$auc, digits = 2)
T3_RF_10_ci <- ci(T3_RF_10_RC, of = "auc")
T3_RF_10_ci_txt <- "95% CI: 0.44-0.89"

#5 proteins RF
load("../5_proteins/Predict_LGA_T3_RF_5_proteins.RData")
T3_RF_5_freq <- freq
T3_RF_5_pile <- pile

T3_RF_5_a=NULL
for (i in 1:length(T3_RF_5_pile)){
  T3_RF_5_a=rbind(T3_RF_5_a,T3_RF_5_pile[[i]]$outmat)
}
T3_RF_5_RC=roc(response=T3_RF_5_a[,1],predictor=T3_RF_5_a[,2],direction="<")

T3_RF_5_auc <- round(T3_RF_5_RC$auc, digits = 2)
T3_RF_5_auc <- "0.40"
T3_RF_5_ci <- ci(T3_RF_5_RC, of = "auc")
T3_RF_5_ci_txt <- "95% CI: 0.14-0.67"


#Plot
ROC_plot_RF_T3 <- ggroc(list("All proteins" = T3_RF_RC,
                             "50 proteins" = T3_RF_50_RC,
                             "40 proteins" = T3_RF_40_RC,
                             "30 proteins" = T3_RF_30_RC,
                             "20 proteins" = T3_RF_20_RC,
                             "10 proteins" = T3_RF_10_RC,
                             "5 proteins" = T3_RF_5_RC),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("navy", "slateblue1","royalblue3", "mediumblue", "dodgerblue","deepskyblue3", "turquoise2")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Random forest\nVisit 3 (Week 28-34)") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T3_RF_auc, ', ', T3_RF_ci_txt),
                    paste0('AUC = ', T3_RF_50_auc, ', ', T3_RF_50_ci_txt), 
                    paste0('AUC = ', T3_RF_40_auc, ', ', T3_RF_40_ci_txt), 
                    paste0('AUC = ', T3_RF_30_auc, ', ', T3_RF_30_ci_txt), 
                    paste0('AUC = ', T3_RF_20_auc, ', ', T3_RF_20_ci_txt),
                    paste0('AUC = ', T3_RF_10_auc, ', ', T3_RF_10_ci_txt),
                    paste0('AUC = ', T3_RF_5_auc, ', ', T3_RF_5_ci_txt)),
           color = c("navy", "slateblue1","royalblue3","mediumblue", "dodgerblue","deepskyblue3", "turquoise2"), 
           size = 3, 
           x = 0.75, 
           y = c(0.30,0.25, 0.2, 0.15, 0.10, 0.05, 0.0))

#Boxplot
T3_RF_a_mod <- as.data.frame(T3_RF_a)
T3_RF_a_mod$Model <- "All proteins"

T3_RF_50_a_mod <- as.data.frame(T3_RF_50_a)
T3_RF_50_a_mod$Model <- "50 proteins"

T3_RF_40_a_mod <- as.data.frame(T3_RF_40_a)
T3_RF_40_a_mod$Model <- "40 proteins"

T3_RF_30_a_mod <- as.data.frame(T3_RF_30_a)
T3_RF_30_a_mod$Model <- "30 proteins"

T3_RF_20_a_mod <- as.data.frame(T3_RF_20_a)
T3_RF_20_a_mod$Model <- "20 proteins"

T3_RF_10_a_mod <- as.data.frame(T3_RF_10_a)
T3_RF_10_a_mod$Model <- "10 proteins"

T3_RF_5_a_mod <- as.data.frame(T3_RF_5_a)
T3_RF_5_a_mod$Model <- "5 proteins"


T3_RF_res <- rbind(T3_RF_a_mod, T3_RF_50_a_mod, T3_RF_40_a_mod, T3_RF_30_a_mod, T3_RF_20_a_mod, T3_RF_10_a_mod, T3_RF_5_a_mod)
T3_RF_res$Bwkat <- ifelse(T3_RF_res$true_labels == 1, "LGA", "AGA")
T3_RF_res$Model <- factor(T3_RF_res$Model, levels = c("All proteins",
                                                      "50 proteins",
                                                      "40 proteins",
                                                      "30 proteins",
                                                      "20 proteins",
                                                      "10 proteins",
                                                      "5 proteins"))

T3_RF_box <- ggplot(T3_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "navy",
                              "50 proteins" = "slateblue1",
                              "40 proteins" = "royalblue3",
                              "30 proteins" = "mediumblue",
                              "20 proteins" = "dodgerblue",
                              "10 proteins" = "deepskyblue3",
                              "5 proteins" = "turquoise2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  #ggtitle("Visit 3\nWeek 28-34") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))


# Compute Precision-Recall AUC
T3_RF_PR <- pr.curve(scores.class0 = T3_RF_a[,2], weights.class0 = T3_RF_a[,1], curve = TRUE)
T3_RF_PR_auc <- round(T3_RF_PR$auc.integral, digits = 2)
T3_RF_PR_df <- data.frame(Recall = T3_RF_PR$curve[, 1], Precision = T3_RF_PR$curve[, 2], Model = "All proteins")

T3_RF_50_PR <- pr.curve(scores.class0 = T3_RF_50_a[,2], weights.class0 = T3_RF_50_a[,1], curve = TRUE)
T3_RF_50_PR_auc <- round(T3_RF_50_PR$auc.integral, digits = 2)
T3_RF_50_PR_df <- data.frame(Recall = T3_RF_50_PR$curve[, 1], Precision = T3_RF_50_PR$curve[, 2], Model = "50 proteins")

T3_RF_40_PR <- pr.curve(scores.class0 = T3_RF_40_a[,2], weights.class0 = T3_RF_40_a[,1], curve = TRUE)
T3_RF_40_PR_auc <- round(T3_RF_40_PR$auc.integral, digits = 2)
T3_RF_40_PR_auc <- "0.20"
T3_RF_40_PR_df <- data.frame(Recall = T3_RF_40_PR$curve[, 1], Precision = T3_RF_40_PR$curve[, 2], Model = "40 proteins")

T3_RF_30_PR <- pr.curve(scores.class0 = T3_RF_30_a[,2], weights.class0 = T3_RF_30_a[,1], curve = TRUE)
T3_RF_30_PR_auc <- round(T3_RF_30_PR$auc.integral, digits = 2)
T3_RF_30_PR_auc <- "0.20"
T3_RF_30_PR_df <- data.frame(Recall = T3_RF_30_PR$curve[, 1], Precision = T3_RF_30_PR$curve[, 2], Model = "30 proteins")

T3_RF_20_PR <- pr.curve(scores.class0 = T3_RF_20_a[,2], weights.class0 = T3_RF_20_a[,1], curve = TRUE)
T3_RF_20_PR_auc <- round(T3_RF_20_PR$auc.integral, digits = 2)
T3_RF_20_PR_df <- data.frame(Recall = T3_RF_20_PR$curve[, 1], Precision = T3_RF_20_PR$curve[, 2], Model = "20 proteins")

T3_RF_10_PR <- pr.curve(scores.class0 = T3_RF_10_a[,2], weights.class0 = T3_RF_10_a[,1], curve = TRUE)
T3_RF_10_PR_auc <- round(T3_RF_10_PR$auc.integral, digits = 2)
T3_RF_10_PR_df <- data.frame(Recall = T3_RF_10_PR$curve[, 1], Precision = T3_RF_10_PR$curve[, 2], Model = "10 proteins")

T3_RF_5_PR <- pr.curve(scores.class0 = T3_RF_5_a[,2], weights.class0 = T3_RF_5_a[,1], curve = TRUE)
T3_RF_5_PR_auc <- round(T3_RF_5_PR$auc.integral, digits = 2)
T3_RF_5_PR_df <- data.frame(Recall = T3_RF_5_PR$curve[, 1], Precision = T3_RF_5_PR$curve[, 2], Model = "5 proteins")

T3_RF_pr_data <- bind_rows(T3_RF_PR_df, T3_RF_50_PR_df, T3_RF_40_PR_df, T3_RF_30_PR_df, T3_RF_20_PR_df, T3_RF_10_PR_df, T3_RF_5_PR_df)

T3_RF_pr_data$Model <- factor(T3_RF_pr_data$Model, levels = c("All proteins",
                                                              "50 proteins",
                                                              "40 proteins",
                                                              "30 proteins",
                                                              "20 proteins",
                                                              "10 proteins",
                                                              "5 proteins"))


T3_RF_PR <- ggplot(T3_RF_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "navy",
                               "50 proteins" = "slateblue1",
                               "40 proteins" = "royalblue3",
                               "30 proteins" = "mediumblue",
                               "20 proteins" = "dodgerblue",
                               "10 proteins" = "deepskyblue3",
                               "5 proteins" = "turquoise2")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('PR-AUC = ', T3_RF_PR_auc),
                    paste0('PR-AUC = ', T3_RF_50_PR_auc), 
                    paste0('PR-AUC = ', T3_RF_40_PR_auc), 
                    paste0('PR-AUC = ', T3_RF_30_PR_auc), 
                    paste0('PR-AUC = ', T3_RF_20_PR_auc),
                    paste0('PR-AUC = ', T3_RF_10_PR_auc),
                    paste0('PR-AUC = ', T3_RF_5_PR_auc)),
           color = c("navy", "slateblue1","royalblue3","mediumblue", "dodgerblue","deepskyblue3", "turquoise2"), 
           size = 3, 
           x = 0.9, 
           y = c(0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4)) 



# Roc Plots ---------------------------------------------------------------


#ROC plots
get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

EN_for_legend_roc <- ggroc(list("All proteins" = T3_EN_RC,
                                "50 proteins" = T3_EN_50_RC,
                                "40 proteins" = T3_EN_40_RC,
                                "30 proteins" = T3_EN_30_RC,
                                "20 proteins" = T3_EN_20_RC,
                                "10 proteins" = T3_EN_10_RC,
                                "5 proteins" = T3_EN_5_RC),
                           aes = c("colour", "size"),
                           legacy.axes = TRUE) +
  scale_color_manual(values = c("maroon4", "mediumpurple4", "mediumorchid2", "mediumpurple2", "maroon2", "pink3", "brown2")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  theme_bw() +
  theme(legend.title = element_blank()) 

EN_legend_roc <- get_legend(EN_for_legend_roc)

EN_ROC_plots <- grid.arrange(ROC_plot_EN_T1, ROC_plot_EN_T2, ROC_plot_EN_T3, EN_legend_roc,
                             ncol = 4, 
                             widths = c(0.8, 0.8, 0.8, 0.3))


RF_for_legend_roc <- ggroc(list("All proteins" = T3_RF_RC,
                                "50 proteins" = T3_RF_50_RC,
                                "40 proteins" = T3_RF_40_RC,
                                "30 proteins" = T3_RF_30_RC,
                                "20 proteins" = T3_RF_20_RC,
                                "10 proteins" = T3_RF_10_RC,
                                "5 proteins" = T3_RF_5_RC),
                           aes = c("colour", "size"),
                           legacy.axes = TRUE) +
  scale_color_manual(values = c("navy", "slateblue1","royalblue3", "mediumblue", "dodgerblue","deepskyblue3", "turquoise2")) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  theme_bw() +
  theme(legend.title = element_blank()) 

RF_legend_roc <- get_legend(RF_for_legend_roc)


RF_ROC_plots <- grid.arrange(ROC_plot_RF_T1, ROC_plot_RF_T2, ROC_plot_RF_T3, RF_legend_roc,
                             ncol = 4, 
                             widths = c(0.8, 0.8, 0.8, 0.3))

ROC_plots <- ggarrange(RF_ROC_plots, EN_ROC_plots,
                       nrow = 2, 
                       heights = c(1, 1),
                       labels = c("A)", "B)"))


ggsave(filename= "ROC_LGA_all_to_5_prots.png",
       plot = ROC_plots,
       device = "png",
       path = "../",
       width = 35,
       height = 22,
       units = "cm",
       bg = "white")



# ROC and box plots -------------------------------------------------------

#ROC and box plots
RF_for_legend_box <- ggplot(T3_RF_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "navy",
                              "50 proteins" = "slateblue1",
                              "40 proteins" = "royalblue3",
                              "30 proteins" = "mediumblue",
                              "20 proteins" = "dodgerblue",
                              "10 proteins" = "deepskyblue3",
                              "5 proteins" = "turquoise2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Visit 3\nWeek 28-34") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

RF_legend_box <- get_legend(RF_for_legend_box)

EN_for_legend_box <- ggplot(T3_EN_res, aes(x = Bwkat, y=predicted_probs, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values= c("All proteins" = "maroon4",
                              "50 proteins" = "mediumpurple4",
                              "40 proteins" = "mediumorchid2",
                              "30 proteins" = "mediumpurple2",
                              "20 proteins" = "maroon2",
                              "10 proteins" = "pink3",
                              "5 proteins" = "brown2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Visit 3\nWeek 28-34") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=14))

EN_legend_box <- get_legend(EN_for_legend_box)

RF_box_plots <- grid.arrange(T1_RF_box, T2_RF_box, T3_RF_box, RF_legend_box, 
                             ncol = 4, 
                             widths = c(0.8, 0.8, 0.8, 0.3))

ROC_BOX_RF_plots <- ggarrange(RF_ROC_plots, RF_box_plots,
                       nrow = 2, 
                       heights = c(1, 1),
                       labels = c("A)", "B)"))


ggsave(filename= "RF_ROC_BOX_LGA_all_to_5_prots.png",
       plot = ROC_BOX_RF_plots,
       device = "png",
       path = "../",
       width = 35,
       height = 22,
       units = "cm",
       bg = "white")


EN_box_plots <- grid.arrange(T1_EN_box, T2_EN_box, T3_EN_box, EN_legend_box, 
                             ncol = 4, 
                             widths = c(0.8, 0.8, 0.8, 0.3))

ROC_BOX_EN_plots <- ggarrange(EN_ROC_plots, EN_box_plots,
                              nrow = 2, 
                              heights = c(1, 1),
                              labels = c("A)", "B)"))


ggsave(filename= "EN_ROC_BOX_LGA_all_to_5_prots.png",
       plot = ROC_BOX_EN_plots,
       device = "png",
       path = "../",
       width = 35,
       height = 22,
       units = "cm",
       bg = "white")



# Roc, Box and PR plots ---------------------------------------------------

RF_PR_for_legend <- ggplot(T3_RF_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "navy",
                               "50 proteins" = "slateblue1",
                               "40 proteins" = "royalblue3",
                               "30 proteins" = "mediumblue",
                               "20 proteins" = "dodgerblue",
                               "10 proteins" = "deepskyblue3",
                               "5 proteins" = "turquoise2")) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank()) 

RF_legend_PR <- get_legend(RF_PR_for_legend)


EN_PR_for_legend <- ggplot(T1_EN_pr_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line() +
  #geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  scale_color_manual(values= c("All proteins" = "maroon4",
                               "50 proteins" = "mediumpurple4",
                               "40 proteins" = "mediumorchid2",
                               "30 proteins" = "mediumpurple2",
                               "20 proteins" = "maroon2",
                               "10 proteins" = "pink3",
                               "5 proteins" = "brown2")) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank()) 

EN_legend_PR <- get_legend(EN_PR_for_legend)


EN_PR_plots <- grid.arrange(T1_EN_PR, T2_EN_PR, T3_EN_PR, EN_legend_PR, 
                             ncol = 4, 
                             widths = c(0.8, 0.8, 0.8, 0.3))

RF_PR_plots <- grid.arrange(T1_RF_PR, T2_RF_PR, T3_RF_PR, RF_legend_PR, 
                            ncol = 4, 
                            widths = c(0.8, 0.8, 0.8, 0.3))


ROC_BOX_PR_EN_plots <- ggarrange(EN_ROC_plots, EN_box_plots, EN_PR_plots,
                              nrow = 3, 
                              heights = c(1, 1, 1),
                              labels = c("A)", "B)", "C)"))


ggsave(filename= "EN_ROC_BOX_PR_LGA_all_to_5_prots.png",
       plot = ROC_BOX_PR_EN_plots,
       device = "png",
       path = "../",
       width = 35,
       height = 30,
       units = "cm",
       bg = "white")

ROC_BOX_PR_RF_plots <- ggarrange(RF_ROC_plots, RF_box_plots, RF_PR_plots,
                                 nrow = 3, 
                                 heights = c(1, 1, 1),
                                 labels = c("A)", "B)", "C)"))


ggsave(filename= "RF_ROC_BOX_PR_LGA_all_to_5_prots.png",
       plot = ROC_BOX_PR_RF_plots,
       device = "png",
       path = "../",
       width = 35,
       height = 30,
       units = "cm",
       bg = "white")
