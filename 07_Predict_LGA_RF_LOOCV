rm(list=ls())
#
library(pROC)
#require(ROCR)
require(caret)
require(MASS)
library(tidymodels) 
library(recipes)
library(limma)
library(randomForest)
library(Biobase)
#library(tune)
library(PRROC)


# All proteins ------------------------------------------------------------


rm(list=ls())

modtype="RF" #RF #EN
tpoints=c("1","2","3")
set.seed(1)
#trim = 1
for(trim in tpoints){
  
  load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
  MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("LGA", "AGA")]
  
  esetbig=exprs(MoM_exprset_STORK_LGA_SGA)
  anobig=pData(MoM_exprset_STORK_LGA_SGA)
  anobig <- anobig %>%
    dplyr::rename("bw_zscore" = `Birthweight z-score`)
  
  anosafe=anobig[anobig$TimePoint==trim,]
  esetsafe=esetbig[,rownames(anosafe)]
  
  all(rownames(anosafe)==colnames(esetsafe))
  anosafe$G=factor(ifelse(anosafe$Bwkat=="LGA","D","C"))
  #anosafe have G as "C", "D"
  
  pile=list()
  all_true_labels <- c()
  all_predicted_probs <- c()
  f1_values <- c() # Store F1 score for each LOOCV iteration
  balanced_acc_values <- c() # Store Balanced Accuracy for each LOOCV iteration
  
  # Leave-One-Out Cross-Validation
  for(ite in 1:nrow(anosafe)){
    cat(ite);cat("\n")
    #cat(1);cat("\n")
    ano=anosafe
    esetn=esetsafe
    
    ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[ite],"validation","discovery")
    #ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[1],"validation","discovery")
    anothis=ano
    
    ano_tr=ano[ano$subset_Bwkat=="discovery",]
    esetn_tr=esetn[,rownames(ano_tr)]
    
    #####prepare test data 
    ano_t=anothis[anothis$subset_Bwkat=="validation",]
    esetn_t=matrix(esetsafe[,rownames(ano_t)], ncol = 1)
    colnames(esetn_t) <- rownames(ano_t)
    rownames(esetn_t) <-rownames(esetn_tr)
    
    #select top proteins by limma test
    #design <- model.matrix(~0+G,ano_tr) 
    #colnames(design)<-gsub("G","",colnames(design))
    
    #pe_v_control_cont = makeContrasts(
    #  D - C,
    #  levels = design)
    #fit_noBayes = lmFit(esetn_tr, design)
    #fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
    #fit_pe_v_control = eBayes(fit_pe_v_control)
    #deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
    #sel=rownames(deT1)[1:20]
    
    #Choose proteins for reduced models
    #train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
    #test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
    #rownames(test) <- colnames(esetn_t)
    
    #Models with all proteins
    train=data.frame(G=factor(ano_tr$G),t(esetn_tr))
    test=data.frame(G=factor(ano_t$G),t(esetn_t))
    
    
    if(modtype=="EN"){
      ####select genes by lasso
      #EN
      p_recipe <- recipe(G ~ ., data = train) %>% 
        step_zv(all_numeric(), -all_outcomes()) 
      
      lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      lasso_wf <- workflow() %>%
        add_recipe(p_recipe)
      
      lasso_fit <- lasso_wf %>%
        add_model(lasso) %>%
        fit(data = train)
      
      lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      # Define the search grid for hyperparameter tuning 
      lasso_grid <- grid_regular(penalty(), mixture(),levels = 10) 
      
      # Tune the hyperparameters using split sample
      lasso_res <- tune_grid(lasso_wf %>% 
                               add_model(lasso_tun), 
                             resamples=apparent(train),
                             grid = lasso_grid,
                             control=control_grid(parallel_over="everything", 
                                                  verbose = TRUE))
      
      #Tune hyperparameters by 5-fold CV
      #folds <- vfold_cv(train, v = 5)
      #lasso_res <- tune_grid(
      #  lasso_wf %>% add_model(lasso_tun), 
      #  resamples = folds,
      #  grid = lasso_grid,
      #  control = control_grid(parallel_over = "everything", verbose = TRUE)
      #)
      
      lowest_auc <- lasso_res %>%
        select_best("roc_auc")
      
      lasso_model_tun <- finalize_workflow(
        lasso_wf %>% add_model(lasso_tun),
        lowest_auc
      )%>%  fit(data = train)
      
      final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
        tidy()
      #
      ###Predictions
      lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
      #final$term[final$estimate!=0][-1]
      
    }else{
      modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
      lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
    }
    
    # Collect the true labels and predicted probabilities
    true_labels = as.numeric(test$G == "D")
    predicted_probs = lasso_pred_tun
    all_true_labels = c(all_true_labels, true_labels)
    all_predicted_probs = c(all_predicted_probs, predicted_probs)
    pred_binary = ifelse(predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
    
    # Confusion matrix components
    tp = sum(true_labels == 1 & pred_binary == 1)
    fp = sum(true_labels == 0 & pred_binary == 1)
    fn = sum(true_labels == 1 & pred_binary == 0)
    tn = sum(true_labels == 0 & pred_binary == 0)
    
    # Compute precision and recall
    precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    # Compute F1 score for this fold
    f1_score = ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    # Store F1 score for visualization
    f1_values[ite] <- f1_score
    
    
    # Compute sensitivity (recall for positive class)
    sensitivity = ifelse(sum(true_labels == 1) == 0, 0, sum(true_labels == 1 & pred_binary == 1) / sum(true_labels == 1))
    
    # Compute specificity (recall for negative class)
    specificity = ifelse(sum(true_labels == 0) == 0, 0, sum(true_labels == 0 & pred_binary == 0) / sum(true_labels == 0))
    
    # Compute balanced accuracy for this fold
    balanced_acc = (sensitivity + specificity) / 2
    
    # Store the balanced accuracy for visualization later
    balanced_acc_values[ite] <- balanced_acc
    
    #outcs=as.numeric(test$G=="D")
    #out=lasso_pred_tun
    #outcs = true_labels
    #out = predicted_probs
    outmat=cbind(true_labels,predicted_probs) #true outcomes for current test set side by side with predictited risk scores
    rownames(outmat) = rownames(test)
    sel <- rownames(esetn_tr)
    pile[[ite]]<-list(outmat=outmat,feat=sel)
    #pile[[1]]<-list(outmat=outmat,feat=sel)
    
  }
  
  #create the full matrix of predicted risk scores and true outcomes
  a=NULL
  for (i in 1:length(pile)){
    a=rbind(a,pile[[i]]$outmat)
  }
  
  # Compute overall F1 score across all iterations
  a_df <- as.data.frame(a)
  a_df$pred_binary = ifelse(a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
  
  tp = sum(a_df$true_labels == 1 & a_df$pred_binary == 1)
  fp = sum(a_df$true_labels == 0 & a_df$pred_binary == 1)
  fn = sum(a_df$true_labels == 1 & a_df$pred_binary == 0)
  
  precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  
  #Calculate F1 score
  if (precision + recall == 0) {
    f1_score <- 0
  } else {
    f1_score <- 2 * (precision * recall) / (precision + recall)
  }

  # Compute overall Balanced Accuracy across all iterations
  sensitivity_all = ifelse(sum(a_df$outcs == 1) == 0, 0, sum(a_df$outcs == 1 & a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
  specificity_all = ifelse(sum(a_df$outcs == 0) == 0, 0, sum(a_df$outcs == 0 & a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))
  
  balanced_accuracy = (sensitivity_all + specificity_all) / 2
  
  
  # Compute Precision-Recall AUC
  PR <- pr.curve(scores.class0 = a[,2], weights.class0 = a[,1], curve = TRUE)
  
  # Extract AUC-PR value
  auc_pr <- round(PR$auc.integral, digits = 2)
  
  
  
  # Cross validation stats
  stat=table(unlist(lapply(pile,function(x){(x$feat)})))
  stat=sort(stat,decreasing=TRUE)
  
  freq=cbind(data.frame(stat))
  save(pile,freq,f1_values,f1_score,balanced_acc_values,balanced_accuracy,a_df,auc_pr,file=paste("../Predict_LGA_new/All_proteins/Predict_LGA_T",trim,"_",modtype,"_all_proteins.RData",sep=""))
}


# 50 proteins -------------------------------------------------------------

rm(list=ls())

modtype="RF" #RF #EN
tpoints=c("1","2","3")
set.seed(1)
#trim = 1
for(trim in tpoints){
  
  load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
  MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("LGA", "AGA")]
  
  esetbig=exprs(MoM_exprset_STORK_LGA_SGA)
  anobig=pData(MoM_exprset_STORK_LGA_SGA)
  anobig <- anobig %>%
    dplyr::rename("bw_zscore" = `Birthweight z-score`)
  
  anosafe=anobig[anobig$TimePoint==trim,]
  esetsafe=esetbig[,rownames(anosafe)]
  
  all(rownames(anosafe)==colnames(esetsafe))
  anosafe$G=factor(ifelse(anosafe$Bwkat=="LGA","D","C"))
  #anosafe have G as "C", "D"
  
  pile=list()
  all_true_labels <- c()
  all_predicted_probs <- c()
  f1_values <- c() # Store F1 score for each LOOCV iteration
  balanced_acc_values <- c() # Store Balanced Accuracy for each LOOCV iteration
  
  # Leave-One-Out Cross-Validation
  for(ite in 1:nrow(anosafe)){
    cat(ite);cat("\n")
    #cat(1);cat("\n")
    ano=anosafe
    esetn=esetsafe
    
    ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[ite],"validation","discovery")
    #ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[1],"validation","discovery")
    anothis=ano
    
    ano_tr=ano[ano$subset_Bwkat=="discovery",]
    esetn_tr=esetn[,rownames(ano_tr)]
    
    #####prepare test data 
    ano_t=anothis[anothis$subset_Bwkat=="validation",]
    esetn_t=matrix(esetsafe[,rownames(ano_t)], ncol = 1)
    colnames(esetn_t) <- rownames(ano_t)
    rownames(esetn_t) <-rownames(esetn_tr)
    
    #select top proteins by limma test
    design <- model.matrix(~0+G,ano_tr) 
    colnames(design)<-gsub("G","",colnames(design))
    
    pe_v_control_cont = makeContrasts(
      D - C,
      levels = design)
    fit_noBayes = lmFit(esetn_tr, design)
    fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
    fit_pe_v_control = eBayes(fit_pe_v_control)
    deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
    sel=rownames(deT1)[1:50]
    
    #Choose proteins for reduced models
    train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
    test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
    rownames(test) <- colnames(esetn_t)
    
    #Models with all proteins
    #train=data.frame(G=factor(ano_tr$G),t(esetn_tr))
    #test=data.frame(G=factor(ano_t$G),t(esetn_t))
    
    
    if(modtype=="EN"){
      ####select genes by lasso
      #EN
      p_recipe <- recipe(G ~ ., data = train) %>% 
        step_zv(all_numeric(), -all_outcomes()) 
      
      lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      lasso_wf <- workflow() %>%
        add_recipe(p_recipe)
      
      lasso_fit <- lasso_wf %>%
        add_model(lasso) %>%
        fit(data = train)
      
      lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      # Define the search grid for hyperparameter tuning 
      lasso_grid <- grid_regular(penalty(), mixture(),levels = 10) 
      
      # Tune the hyperparameters using split sample
      lasso_res <- tune_grid(lasso_wf %>% 
                               add_model(lasso_tun), 
                             resamples=apparent(train),
                             grid = lasso_grid,
                             control=control_grid(parallel_over="everything", 
                                                  verbose = TRUE))
      
      #Tune hyperparameters by 5-fold CV
      #folds <- vfold_cv(train, v = 5)
      #lasso_res <- tune_grid(
      #  lasso_wf %>% add_model(lasso_tun), 
      #  resamples = folds,
      #  grid = lasso_grid,
      #  control = control_grid(parallel_over = "everything", verbose = TRUE)
      #)
      
      lowest_auc <- lasso_res %>%
        select_best("roc_auc")
      
      lasso_model_tun <- finalize_workflow(
        lasso_wf %>% add_model(lasso_tun),
        lowest_auc
      )%>%  fit(data = train)
      
      final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
        tidy()
      #
      ###Predictions
      lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
      #final$term[final$estimate!=0][-1]
      
    }else{
      modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
      lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
    }
    
    # Collect the true labels and predicted probabilities
    true_labels = as.numeric(test$G == "D")
    predicted_probs = lasso_pred_tun
    all_true_labels = c(all_true_labels, true_labels)
    all_predicted_probs = c(all_predicted_probs, predicted_probs)
    pred_binary = ifelse(predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
    
    # Confusion matrix components
    tp = sum(true_labels == 1 & pred_binary == 1)
    fp = sum(true_labels == 0 & pred_binary == 1)
    fn = sum(true_labels == 1 & pred_binary == 0)
    tn = sum(true_labels == 0 & pred_binary == 0)
    
    # Compute precision and recall
    precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    # Compute F1 score for this fold
    f1_score = ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    # Store F1 score for visualization
    f1_values[ite] <- f1_score
    
    
    # Compute sensitivity (recall for positive class)
    sensitivity = ifelse(sum(true_labels == 1) == 0, 0, sum(true_labels == 1 & pred_binary == 1) / sum(true_labels == 1))
    
    # Compute specificity (recall for negative class)
    specificity = ifelse(sum(true_labels == 0) == 0, 0, sum(true_labels == 0 & pred_binary == 0) / sum(true_labels == 0))
    
    # Compute balanced accuracy for this fold
    balanced_acc = (sensitivity + specificity) / 2
    
    # Store the balanced accuracy for visualization later
    balanced_acc_values[ite] <- balanced_acc
    
    #outcs=as.numeric(test$G=="D")
    #out=lasso_pred_tun
    #outcs = true_labels
    #out = predicted_probs
    outmat=cbind(true_labels,predicted_probs) #true outcomes for current test set side by side with predictited risk scores
    rownames(outmat) = rownames(test)
    sel <- rownames(esetn_tr)
    pile[[ite]]<-list(outmat=outmat,feat=sel)
    #pile[[1]]<-list(outmat=outmat,feat=sel)
    
  }
  
  #create the full matrix of predicted risk scores and true outcomes
  a=NULL
  for (i in 1:length(pile)){
    a=rbind(a,pile[[i]]$outmat)
  }
  
  # Compute overall F1 score across all iterations
  a_df <- as.data.frame(a)
  a_df$pred_binary = ifelse(a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
  
  tp = sum(a_df$true_labels == 1 & a_df$pred_binary == 1)
  fp = sum(a_df$true_labels == 0 & a_df$pred_binary == 1)
  fn = sum(a_df$true_labels == 1 & a_df$pred_binary == 0)
  
  precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  
  #Calculate F1 score
  if (precision + recall == 0) {
    f1_score <- 0
  } else {
    f1_score <- 2 * (precision * recall) / (precision + recall)
  }
  
  # Compute overall Balanced Accuracy across all iterations
  sensitivity_all = ifelse(sum(a_df$outcs == 1) == 0, 0, sum(a_df$outcs == 1 & a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
  specificity_all = ifelse(sum(a_df$outcs == 0) == 0, 0, sum(a_df$outcs == 0 & a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))
  
  balanced_accuracy = (sensitivity_all + specificity_all) / 2
  
  
  # Compute Precision-Recall AUC
  PR <- pr.curve(scores.class0 = a[,2], weights.class0 = a[,1], curve = TRUE)
  
  # Extract AUC-PR value
  auc_pr <- round(PR$auc.integral, digits = 2)
  
  
  # Cross validation stats
  stat=table(unlist(lapply(pile,function(x){(x$feat)})))
  stat=sort(stat,decreasing=TRUE)
  
  freq=cbind(data.frame(stat))
  save(pile,freq,f1_values,f1_score,balanced_acc_values,balanced_accuracy,a_df,auc_pr,file=paste("../Predict_LGA_new/50_proteins/Predict_LGA_T",trim,"_",modtype,"_50_proteins.RData",sep=""))
}


# 40 proteins -------------------------------------------------------------

rm(list=ls())

modtype="RF" #RF #EN
tpoints=c("1","2","3")
set.seed(1)
#trim = 1
for(trim in tpoints){
  
  load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
  MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("LGA", "AGA")]
  
  esetbig=exprs(MoM_exprset_STORK_LGA_SGA)
  anobig=pData(MoM_exprset_STORK_LGA_SGA)
  anobig <- anobig %>%
    dplyr::rename("bw_zscore" = `Birthweight z-score`)
  
  anosafe=anobig[anobig$TimePoint==trim,]
  esetsafe=esetbig[,rownames(anosafe)]
  
  all(rownames(anosafe)==colnames(esetsafe))
  anosafe$G=factor(ifelse(anosafe$Bwkat=="LGA","D","C"))
  #anosafe have G as "C", "D"
  
  pile=list()
  all_true_labels <- c()
  all_predicted_probs <- c()
  f1_values <- c() # Store F1 score for each LOOCV iteration
  balanced_acc_values <- c() # Store Balanced Accuracy for each LOOCV iteration
  
  # Leave-One-Out Cross-Validation
  for(ite in 1:nrow(anosafe)){
    cat(ite);cat("\n")
    #cat(1);cat("\n")
    ano=anosafe
    esetn=esetsafe
    
    ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[ite],"validation","discovery")
    #ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[1],"validation","discovery")
    anothis=ano
    
    ano_tr=ano[ano$subset_Bwkat=="discovery",]
    esetn_tr=esetn[,rownames(ano_tr)]
    
    #####prepare test data 
    ano_t=anothis[anothis$subset_Bwkat=="validation",]
    esetn_t=matrix(esetsafe[,rownames(ano_t)], ncol = 1)
    colnames(esetn_t) <- rownames(ano_t)
    rownames(esetn_t) <-rownames(esetn_tr)
    
    #select top proteins by limma test
    design <- model.matrix(~0+G,ano_tr) 
    colnames(design)<-gsub("G","",colnames(design))
    
    pe_v_control_cont = makeContrasts(
      D - C,
      levels = design)
    fit_noBayes = lmFit(esetn_tr, design)
    fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
    fit_pe_v_control = eBayes(fit_pe_v_control)
    deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
    sel=rownames(deT1)[1:40]
    
    #Choose proteins for reduced models
    train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
    test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
    rownames(test) <- colnames(esetn_t)
    
    #Models with all proteins
    #train=data.frame(G=factor(ano_tr$G),t(esetn_tr))
    #test=data.frame(G=factor(ano_t$G),t(esetn_t))
    
    
    if(modtype=="EN"){
      ####select genes by lasso
      #EN
      p_recipe <- recipe(G ~ ., data = train) %>% 
        step_zv(all_numeric(), -all_outcomes()) 
      
      lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      lasso_wf <- workflow() %>%
        add_recipe(p_recipe)
      
      lasso_fit <- lasso_wf %>%
        add_model(lasso) %>%
        fit(data = train)
      
      lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      # Define the search grid for hyperparameter tuning 
      lasso_grid <- grid_regular(penalty(), mixture(),levels = 10) 
      
      # Tune the hyperparameters using split sample
      lasso_res <- tune_grid(lasso_wf %>% 
                               add_model(lasso_tun), 
                             resamples=apparent(train),
                             grid = lasso_grid,
                             control=control_grid(parallel_over="everything", 
                                                  verbose = TRUE))
      
      #Tune hyperparameters by 5-fold CV
      #folds <- vfold_cv(train, v = 5)
      #lasso_res <- tune_grid(
      #  lasso_wf %>% add_model(lasso_tun), 
      #  resamples = folds,
      #  grid = lasso_grid,
      #  control = control_grid(parallel_over = "everything", verbose = TRUE)
      #)
      
      lowest_auc <- lasso_res %>%
        select_best("roc_auc")
      
      lasso_model_tun <- finalize_workflow(
        lasso_wf %>% add_model(lasso_tun),
        lowest_auc
      )%>%  fit(data = train)
      
      final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
        tidy()
      #
      ###Predictions
      lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
      #final$term[final$estimate!=0][-1]
      
    }else{
      modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
      lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
    }
    
    # Collect the true labels and predicted probabilities
    true_labels = as.numeric(test$G == "D")
    predicted_probs = lasso_pred_tun
    all_true_labels = c(all_true_labels, true_labels)
    all_predicted_probs = c(all_predicted_probs, predicted_probs)
    pred_binary = ifelse(predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
    
    # Confusion matrix components
    tp = sum(true_labels == 1 & pred_binary == 1)
    fp = sum(true_labels == 0 & pred_binary == 1)
    fn = sum(true_labels == 1 & pred_binary == 0)
    tn = sum(true_labels == 0 & pred_binary == 0)
    
    # Compute precision and recall
    precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    # Compute F1 score for this fold
    f1_score = ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    # Store F1 score for visualization
    f1_values[ite] <- f1_score
    
    
    # Compute sensitivity (recall for positive class)
    sensitivity = ifelse(sum(true_labels == 1) == 0, 0, sum(true_labels == 1 & pred_binary == 1) / sum(true_labels == 1))
    
    # Compute specificity (recall for negative class)
    specificity = ifelse(sum(true_labels == 0) == 0, 0, sum(true_labels == 0 & pred_binary == 0) / sum(true_labels == 0))
    
    # Compute balanced accuracy for this fold
    balanced_acc = (sensitivity + specificity) / 2
    
    # Store the balanced accuracy for visualization later
    balanced_acc_values[ite] <- balanced_acc
    
    #outcs=as.numeric(test$G=="D")
    #out=lasso_pred_tun
    #outcs = true_labels
    #out = predicted_probs
    outmat=cbind(true_labels,predicted_probs) #true outcomes for current test set side by side with predictited risk scores
    rownames(outmat) = rownames(test)
    sel <- rownames(esetn_tr)
    pile[[ite]]<-list(outmat=outmat,feat=sel)
    #pile[[1]]<-list(outmat=outmat,feat=sel)
    
  }
  
  #create the full matrix of predicted risk scores and true outcomes
  a=NULL
  for (i in 1:length(pile)){
    a=rbind(a,pile[[i]]$outmat)
  }
  
  # Compute overall F1 score across all iterations
  a_df <- as.data.frame(a)
  a_df$pred_binary = ifelse(a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
  
  tp = sum(a_df$true_labels == 1 & a_df$pred_binary == 1)
  fp = sum(a_df$true_labels == 0 & a_df$pred_binary == 1)
  fn = sum(a_df$true_labels == 1 & a_df$pred_binary == 0)
  
  precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  
  #Calculate F1 score
  if (precision + recall == 0) {
    f1_score <- 0
  } else {
    f1_score <- 2 * (precision * recall) / (precision + recall)
  }
  
  # Compute overall Balanced Accuracy across all iterations
  sensitivity_all = ifelse(sum(a_df$outcs == 1) == 0, 0, sum(a_df$outcs == 1 & a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
  specificity_all = ifelse(sum(a_df$outcs == 0) == 0, 0, sum(a_df$outcs == 0 & a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))
  
  balanced_accuracy = (sensitivity_all + specificity_all) / 2
  
  
  # Compute Precision-Recall AUC
  PR <- pr.curve(scores.class0 = a[,2], weights.class0 = a[,1], curve = TRUE)
  
  # Extract AUC-PR value
  auc_pr <- round(PR$auc.integral, digits = 2)
  
  
  # Cross validation stats
  stat=table(unlist(lapply(pile,function(x){(x$feat)})))
  stat=sort(stat,decreasing=TRUE)
  
  freq=cbind(data.frame(stat))
  save(pile,freq,f1_values,f1_score,balanced_acc_values,balanced_accuracy,a_df,auc_pr,file=paste("../Predict_LGA_new/40_proteins/Predict_LGA_T",trim,"_",modtype,"_40_proteins.RData",sep=""))
}


# 30 proteins -------------------------------------------------------------

rm(list=ls())

modtype="RF" #RF #EN
tpoints=c("1","2","3")
set.seed(1)
#trim = 1
for(trim in tpoints){
  
  load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
  MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("LGA", "AGA")]
  
  esetbig=exprs(MoM_exprset_STORK_LGA_SGA)
  anobig=pData(MoM_exprset_STORK_LGA_SGA)
  anobig <- anobig %>%
    dplyr::rename("bw_zscore" = `Birthweight z-score`)
  
  anosafe=anobig[anobig$TimePoint==trim,]
  esetsafe=esetbig[,rownames(anosafe)]
  
  all(rownames(anosafe)==colnames(esetsafe))
  anosafe$G=factor(ifelse(anosafe$Bwkat=="LGA","D","C"))
  #anosafe have G as "C", "D"
  
  pile=list()
  all_true_labels <- c()
  all_predicted_probs <- c()
  f1_values <- c() # Store F1 score for each LOOCV iteration
  balanced_acc_values <- c() # Store Balanced Accuracy for each LOOCV iteration
  
  # Leave-One-Out Cross-Validation
  for(ite in 1:nrow(anosafe)){
    cat(ite);cat("\n")
    #cat(1);cat("\n")
    ano=anosafe
    esetn=esetsafe
    
    ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[ite],"validation","discovery")
    #ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[1],"validation","discovery")
    anothis=ano
    
    ano_tr=ano[ano$subset_Bwkat=="discovery",]
    esetn_tr=esetn[,rownames(ano_tr)]
    
    #####prepare test data 
    ano_t=anothis[anothis$subset_Bwkat=="validation",]
    esetn_t=matrix(esetsafe[,rownames(ano_t)], ncol = 1)
    colnames(esetn_t) <- rownames(ano_t)
    rownames(esetn_t) <-rownames(esetn_tr)
    
    #select top proteins by limma test
    design <- model.matrix(~0+G,ano_tr) 
    colnames(design)<-gsub("G","",colnames(design))
    
    pe_v_control_cont = makeContrasts(
      D - C,
      levels = design)
    fit_noBayes = lmFit(esetn_tr, design)
    fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
    fit_pe_v_control = eBayes(fit_pe_v_control)
    deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
    sel=rownames(deT1)[1:30]
    
    #Choose proteins for reduced models
    train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
    test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
    rownames(test) <- colnames(esetn_t)
    
    #Models with all proteins
    #train=data.frame(G=factor(ano_tr$G),t(esetn_tr))
    #test=data.frame(G=factor(ano_t$G),t(esetn_t))
    
    
    if(modtype=="EN"){
      ####select genes by lasso
      #EN
      p_recipe <- recipe(G ~ ., data = train) %>% 
        step_zv(all_numeric(), -all_outcomes()) 
      
      lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      lasso_wf <- workflow() %>%
        add_recipe(p_recipe)
      
      lasso_fit <- lasso_wf %>%
        add_model(lasso) %>%
        fit(data = train)
      
      lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      # Define the search grid for hyperparameter tuning 
      lasso_grid <- grid_regular(penalty(), mixture(),levels = 10) 
      
      # Tune the hyperparameters using split sample
      lasso_res <- tune_grid(lasso_wf %>% 
                               add_model(lasso_tun), 
                             resamples=apparent(train),
                             grid = lasso_grid,
                             control=control_grid(parallel_over="everything", 
                                                  verbose = TRUE))
      
      #Tune hyperparameters by 5-fold CV
      #folds <- vfold_cv(train, v = 5)
      #lasso_res <- tune_grid(
      #  lasso_wf %>% add_model(lasso_tun), 
      #  resamples = folds,
      #  grid = lasso_grid,
      #  control = control_grid(parallel_over = "everything", verbose = TRUE)
      #)
      
      lowest_auc <- lasso_res %>%
        select_best("roc_auc")
      
      lasso_model_tun <- finalize_workflow(
        lasso_wf %>% add_model(lasso_tun),
        lowest_auc
      )%>%  fit(data = train)
      
      final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
        tidy()
      #
      ###Predictions
      lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
      #final$term[final$estimate!=0][-1]
      
    }else{
      modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
      lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
    }
    
    # Collect the true labels and predicted probabilities
    true_labels = as.numeric(test$G == "D")
    predicted_probs = lasso_pred_tun
    all_true_labels = c(all_true_labels, true_labels)
    all_predicted_probs = c(all_predicted_probs, predicted_probs)
    pred_binary = ifelse(predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
    
    # Confusion matrix components
    tp = sum(true_labels == 1 & pred_binary == 1)
    fp = sum(true_labels == 0 & pred_binary == 1)
    fn = sum(true_labels == 1 & pred_binary == 0)
    tn = sum(true_labels == 0 & pred_binary == 0)
    
    # Compute precision and recall
    precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    # Compute F1 score for this fold
    f1_score = ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    # Store F1 score for visualization
    f1_values[ite] <- f1_score
    
    
    # Compute sensitivity (recall for positive class)
    sensitivity = ifelse(sum(true_labels == 1) == 0, 0, sum(true_labels == 1 & pred_binary == 1) / sum(true_labels == 1))
    
    # Compute specificity (recall for negative class)
    specificity = ifelse(sum(true_labels == 0) == 0, 0, sum(true_labels == 0 & pred_binary == 0) / sum(true_labels == 0))
    
    # Compute balanced accuracy for this fold
    balanced_acc = (sensitivity + specificity) / 2
    
    # Store the balanced accuracy for visualization later
    balanced_acc_values[ite] <- balanced_acc
    
    #outcs=as.numeric(test$G=="D")
    #out=lasso_pred_tun
    #outcs = true_labels
    #out = predicted_probs
    outmat=cbind(true_labels,predicted_probs) #true outcomes for current test set side by side with predictited risk scores
    rownames(outmat) = rownames(test)
    sel <- rownames(esetn_tr)
    pile[[ite]]<-list(outmat=outmat,feat=sel)
    #pile[[1]]<-list(outmat=outmat,feat=sel)
    
  }
  
  #create the full matrix of predicted risk scores and true outcomes
  a=NULL
  for (i in 1:length(pile)){
    a=rbind(a,pile[[i]]$outmat)
  }
  
  # Compute overall F1 score across all iterations
  a_df <- as.data.frame(a)
  a_df$pred_binary = ifelse(a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
  
  tp = sum(a_df$true_labels == 1 & a_df$pred_binary == 1)
  fp = sum(a_df$true_labels == 0 & a_df$pred_binary == 1)
  fn = sum(a_df$true_labels == 1 & a_df$pred_binary == 0)
  
  precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  
  #Calculate F1 score
  if (precision + recall == 0) {
    f1_score <- 0
  } else {
    f1_score <- 2 * (precision * recall) / (precision + recall)
  }
  
  # Compute overall Balanced Accuracy across all iterations
  sensitivity_all = ifelse(sum(a_df$outcs == 1) == 0, 0, sum(a_df$outcs == 1 & a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
  specificity_all = ifelse(sum(a_df$outcs == 0) == 0, 0, sum(a_df$outcs == 0 & a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))
  
  balanced_accuracy = (sensitivity_all + specificity_all) / 2
  
  
  # Compute Precision-Recall AUC
  PR <- pr.curve(scores.class0 = a[,2], weights.class0 = a[,1], curve = TRUE)
  
  # Extract AUC-PR value
  auc_pr <- round(PR$auc.integral, digits = 2)
  
  
  # Cross validation stats
  stat=table(unlist(lapply(pile,function(x){(x$feat)})))
  stat=sort(stat,decreasing=TRUE)
  
  freq=cbind(data.frame(stat))
  save(pile,freq,f1_values,f1_score,balanced_acc_values,balanced_accuracy,a_df,auc_pr,file=paste("../Predict_LGA_new/30_proteins/Predict_LGA_T",trim,"_",modtype,"_30_proteins.RData",sep=""))
}



# 20 proteins -------------------------------------------------------------


rm(list=ls())

modtype="RF" #RF #EN
tpoints=c("1","2","3")
set.seed(1)
#trim = 1
for(trim in tpoints){
  
  load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
  MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("LGA", "AGA")]
  
  esetbig=exprs(MoM_exprset_STORK_LGA_SGA)
  anobig=pData(MoM_exprset_STORK_LGA_SGA)
  anobig <- anobig %>%
    dplyr::rename("bw_zscore" = `Birthweight z-score`)
  
  anosafe=anobig[anobig$TimePoint==trim,]
  esetsafe=esetbig[,rownames(anosafe)]
  
  all(rownames(anosafe)==colnames(esetsafe))
  anosafe$G=factor(ifelse(anosafe$Bwkat=="LGA","D","C"))
  #anosafe have G as "C", "D"
  
  pile=list()
  all_true_labels <- c()
  all_predicted_probs <- c()
  f1_values <- c() # Store F1 score for each LOOCV iteration
  balanced_acc_values <- c() # Store Balanced Accuracy for each LOOCV iteration
  
  # Leave-One-Out Cross-Validation
  for(ite in 1:nrow(anosafe)){
    cat(ite);cat("\n")
    #cat(1);cat("\n")
    ano=anosafe
    esetn=esetsafe
    
    ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[ite],"validation","discovery")
    #ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[1],"validation","discovery")
    anothis=ano
    
    ano_tr=ano[ano$subset_Bwkat=="discovery",]
    esetn_tr=esetn[,rownames(ano_tr)]
    
    #####prepare test data 
    ano_t=anothis[anothis$subset_Bwkat=="validation",]
    esetn_t=matrix(esetsafe[,rownames(ano_t)], ncol = 1)
    colnames(esetn_t) <- rownames(ano_t)
    rownames(esetn_t) <-rownames(esetn_tr)
    
    #select top proteins by limma test
    design <- model.matrix(~0+G,ano_tr) 
    colnames(design)<-gsub("G","",colnames(design))
    
    pe_v_control_cont = makeContrasts(
      D - C,
      levels = design)
    fit_noBayes = lmFit(esetn_tr, design)
    fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
    fit_pe_v_control = eBayes(fit_pe_v_control)
    deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
    sel=rownames(deT1)[1:20]
    
    #Choose proteins for reduced models
    train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
    test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
    rownames(test) <- colnames(esetn_t)
    
    #Models with all proteins
    #train=data.frame(G=factor(ano_tr$G),t(esetn_tr))
    #test=data.frame(G=factor(ano_t$G),t(esetn_t))
    
    
    if(modtype=="EN"){
      ####select genes by lasso
      #EN
      p_recipe <- recipe(G ~ ., data = train) %>% 
        step_zv(all_numeric(), -all_outcomes()) 
      
      lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      lasso_wf <- workflow() %>%
        add_recipe(p_recipe)
      
      lasso_fit <- lasso_wf %>%
        add_model(lasso) %>%
        fit(data = train)
      
      lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      # Define the search grid for hyperparameter tuning 
      lasso_grid <- grid_regular(penalty(), mixture(),levels = 10) 
      
      # Tune the hyperparameters using split sample
      lasso_res <- tune_grid(lasso_wf %>% 
                               add_model(lasso_tun), 
                             resamples=apparent(train),
                             grid = lasso_grid,
                             control=control_grid(parallel_over="everything", 
                                                  verbose = TRUE))
      
      #Tune hyperparameters by 5-fold CV
      #folds <- vfold_cv(train, v = 5)
      #lasso_res <- tune_grid(
      #  lasso_wf %>% add_model(lasso_tun), 
      #  resamples = folds,
      #  grid = lasso_grid,
      #  control = control_grid(parallel_over = "everything", verbose = TRUE)
      #)
      
      lowest_auc <- lasso_res %>%
        select_best("roc_auc")
      
      lasso_model_tun <- finalize_workflow(
        lasso_wf %>% add_model(lasso_tun),
        lowest_auc
      )%>%  fit(data = train)
      
      final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
        tidy()
      #
      ###Predictions
      lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
      #final$term[final$estimate!=0][-1]
      
    }else{
      modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
      lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
    }
    
    # Collect the true labels and predicted probabilities
    true_labels = as.numeric(test$G == "D")
    predicted_probs = lasso_pred_tun
    all_true_labels = c(all_true_labels, true_labels)
    all_predicted_probs = c(all_predicted_probs, predicted_probs)
    pred_binary = ifelse(predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
    
    # Confusion matrix components
    tp = sum(true_labels == 1 & pred_binary == 1)
    fp = sum(true_labels == 0 & pred_binary == 1)
    fn = sum(true_labels == 1 & pred_binary == 0)
    tn = sum(true_labels == 0 & pred_binary == 0)
    
    # Compute precision and recall
    precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    # Compute F1 score for this fold
    f1_score = ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    # Store F1 score for visualization
    f1_values[ite] <- f1_score
    
    
    # Compute sensitivity (recall for positive class)
    sensitivity = ifelse(sum(true_labels == 1) == 0, 0, sum(true_labels == 1 & pred_binary == 1) / sum(true_labels == 1))
    
    # Compute specificity (recall for negative class)
    specificity = ifelse(sum(true_labels == 0) == 0, 0, sum(true_labels == 0 & pred_binary == 0) / sum(true_labels == 0))
    
    # Compute balanced accuracy for this fold
    balanced_acc = (sensitivity + specificity) / 2
    
    # Store the balanced accuracy for visualization later
    balanced_acc_values[ite] <- balanced_acc
    
    #outcs=as.numeric(test$G=="D")
    #out=lasso_pred_tun
    #outcs = true_labels
    #out = predicted_probs
    outmat=cbind(true_labels,predicted_probs) #true outcomes for current test set side by side with predictited risk scores
    rownames(outmat) = rownames(test)
    sel <- rownames(esetn_tr)
    pile[[ite]]<-list(outmat=outmat,feat=sel)
    #pile[[1]]<-list(outmat=outmat,feat=sel)
    
  }
  
  #create the full matrix of predicted risk scores and true outcomes
  a=NULL
  for (i in 1:length(pile)){
    a=rbind(a,pile[[i]]$outmat)
  }
  
  # Compute overall F1 score across all iterations
  a_df <- as.data.frame(a)
  a_df$pred_binary = ifelse(a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
  
  tp = sum(a_df$true_labels == 1 & a_df$pred_binary == 1)
  fp = sum(a_df$true_labels == 0 & a_df$pred_binary == 1)
  fn = sum(a_df$true_labels == 1 & a_df$pred_binary == 0)
  
  precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  
  #Calculate F1 score
  if (precision + recall == 0) {
    f1_score <- 0
  } else {
    f1_score <- 2 * (precision * recall) / (precision + recall)
  }
  
  # Compute overall Balanced Accuracy across all iterations
  sensitivity_all = ifelse(sum(a_df$outcs == 1) == 0, 0, sum(a_df$outcs == 1 & a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
  specificity_all = ifelse(sum(a_df$outcs == 0) == 0, 0, sum(a_df$outcs == 0 & a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))
  
  balanced_accuracy = (sensitivity_all + specificity_all) / 2
  
  
  # Compute Precision-Recall AUC
  PR <- pr.curve(scores.class0 = a[,2], weights.class0 = a[,1], curve = TRUE)
  
  # Extract AUC-PR value
  auc_pr <- round(PR$auc.integral, digits = 2)
  
  
  # Cross validation stats
  stat=table(unlist(lapply(pile,function(x){(x$feat)})))
  stat=sort(stat,decreasing=TRUE)
  
  freq=cbind(data.frame(stat))
  save(pile,freq,f1_values,f1_score,balanced_acc_values,balanced_accuracy,a_df,auc_pr,file=paste("../Predict_LGA_new/20_proteins/Predict_LGA_T",trim,"_",modtype,"_20_proteins.RData",sep=""))
}

# 10 proteins -------------------------------------------------------------


rm(list=ls())

modtype="RF" #RF #EN
tpoints=c("1","2","3")
set.seed(1)
#trim = 1
for(trim in tpoints){
  
  load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
  MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("LGA", "AGA")]
  
  esetbig=exprs(MoM_exprset_STORK_LGA_SGA)
  anobig=pData(MoM_exprset_STORK_LGA_SGA)
  anobig <- anobig %>%
    dplyr::rename("bw_zscore" = `Birthweight z-score`)
  
  anosafe=anobig[anobig$TimePoint==trim,]
  esetsafe=esetbig[,rownames(anosafe)]
  
  all(rownames(anosafe)==colnames(esetsafe))
  anosafe$G=factor(ifelse(anosafe$Bwkat=="LGA","D","C"))
  #anosafe have G as "C", "D"
  
  pile=list()
  all_true_labels <- c()
  all_predicted_probs <- c()
  f1_values <- c() # Store F1 score for each LOOCV iteration
  balanced_acc_values <- c() # Store Balanced Accuracy for each LOOCV iteration
  
  # Leave-One-Out Cross-Validation
  for(ite in 1:nrow(anosafe)){
    cat(ite);cat("\n")
    #cat(1);cat("\n")
    ano=anosafe
    esetn=esetsafe
    
    ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[ite],"validation","discovery")
    #ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[1],"validation","discovery")
    anothis=ano
    
    ano_tr=ano[ano$subset_Bwkat=="discovery",]
    esetn_tr=esetn[,rownames(ano_tr)]
    
    #####prepare test data 
    ano_t=anothis[anothis$subset_Bwkat=="validation",]
    esetn_t=matrix(esetsafe[,rownames(ano_t)], ncol = 1)
    colnames(esetn_t) <- rownames(ano_t)
    rownames(esetn_t) <-rownames(esetn_tr)
    
    #select top proteins by limma test
    design <- model.matrix(~0+G,ano_tr) 
    colnames(design)<-gsub("G","",colnames(design))
    
    pe_v_control_cont = makeContrasts(
      D - C,
      levels = design)
    fit_noBayes = lmFit(esetn_tr, design)
    fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
    fit_pe_v_control = eBayes(fit_pe_v_control)
    deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
    sel=rownames(deT1)[1:10]
    
    #Choose proteins for reduced models
    train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
    test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
    rownames(test) <- colnames(esetn_t)
    
    #Models with all proteins
    #train=data.frame(G=factor(ano_tr$G),t(esetn_tr))
    #test=data.frame(G=factor(ano_t$G),t(esetn_t))
    
    
    if(modtype=="EN"){
      ####select genes by lasso
      #EN
      p_recipe <- recipe(G ~ ., data = train) %>% 
        step_zv(all_numeric(), -all_outcomes()) 
      
      lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      lasso_wf <- workflow() %>%
        add_recipe(p_recipe)
      
      lasso_fit <- lasso_wf %>%
        add_model(lasso) %>%
        fit(data = train)
      
      lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      # Define the search grid for hyperparameter tuning 
      lasso_grid <- grid_regular(penalty(), mixture(),levels = 10) 
      
      # Tune the hyperparameters using split sample
      lasso_res <- tune_grid(lasso_wf %>% 
                               add_model(lasso_tun), 
                             resamples=apparent(train),
                             grid = lasso_grid,
                             control=control_grid(parallel_over="everything", 
                                                  verbose = TRUE))
      
      #Tune hyperparameters by 5-fold CV
      #folds <- vfold_cv(train, v = 5)
      #lasso_res <- tune_grid(
      #  lasso_wf %>% add_model(lasso_tun), 
      #  resamples = folds,
      #  grid = lasso_grid,
      #  control = control_grid(parallel_over = "everything", verbose = TRUE)
      #)
      
      lowest_auc <- lasso_res %>%
        select_best("roc_auc")
      
      lasso_model_tun <- finalize_workflow(
        lasso_wf %>% add_model(lasso_tun),
        lowest_auc
      )%>%  fit(data = train)
      
      final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
        tidy()
      #
      ###Predictions
      lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
      #final$term[final$estimate!=0][-1]
      
    }else{
      modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
      lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
    }
    
    # Collect the true labels and predicted probabilities
    true_labels = as.numeric(test$G == "D")
    predicted_probs = lasso_pred_tun
    all_true_labels = c(all_true_labels, true_labels)
    all_predicted_probs = c(all_predicted_probs, predicted_probs)
    pred_binary = ifelse(predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
    
    # Confusion matrix components
    tp = sum(true_labels == 1 & pred_binary == 1)
    fp = sum(true_labels == 0 & pred_binary == 1)
    fn = sum(true_labels == 1 & pred_binary == 0)
    tn = sum(true_labels == 0 & pred_binary == 0)
    
    # Compute precision and recall
    precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    # Compute F1 score for this fold
    f1_score = ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    # Store F1 score for visualization
    f1_values[ite] <- f1_score
    
    
    # Compute sensitivity (recall for positive class)
    sensitivity = ifelse(sum(true_labels == 1) == 0, 0, sum(true_labels == 1 & pred_binary == 1) / sum(true_labels == 1))
    
    # Compute specificity (recall for negative class)
    specificity = ifelse(sum(true_labels == 0) == 0, 0, sum(true_labels == 0 & pred_binary == 0) / sum(true_labels == 0))
    
    # Compute balanced accuracy for this fold
    balanced_acc = (sensitivity + specificity) / 2
    
    # Store the balanced accuracy for visualization later
    balanced_acc_values[ite] <- balanced_acc
    
    #outcs=as.numeric(test$G=="D")
    #out=lasso_pred_tun
    #outcs = true_labels
    #out = predicted_probs
    outmat=cbind(true_labels,predicted_probs) #true outcomes for current test set side by side with predictited risk scores
    rownames(outmat) = rownames(test)
    sel <- rownames(esetn_tr)
    pile[[ite]]<-list(outmat=outmat,feat=sel)
    #pile[[1]]<-list(outmat=outmat,feat=sel)
    
  }
  
  #create the full matrix of predicted risk scores and true outcomes
  a=NULL
  for (i in 1:length(pile)){
    a=rbind(a,pile[[i]]$outmat)
  }
  
  # Compute overall F1 score across all iterations
  a_df <- as.data.frame(a)
  a_df$pred_binary = ifelse(a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
  
  tp = sum(a_df$true_labels == 1 & a_df$pred_binary == 1)
  fp = sum(a_df$true_labels == 0 & a_df$pred_binary == 1)
  fn = sum(a_df$true_labels == 1 & a_df$pred_binary == 0)
  
  precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  
  #Calculate F1 score
  if (precision + recall == 0) {
    f1_score <- 0
  } else {
    f1_score <- 2 * (precision * recall) / (precision + recall)
  }
  
  # Compute overall Balanced Accuracy across all iterations
  sensitivity_all = ifelse(sum(a_df$outcs == 1) == 0, 0, sum(a_df$outcs == 1 & a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
  specificity_all = ifelse(sum(a_df$outcs == 0) == 0, 0, sum(a_df$outcs == 0 & a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))
  
  balanced_accuracy = (sensitivity_all + specificity_all) / 2
  
  
  # Compute Precision-Recall AUC
  PR <- pr.curve(scores.class0 = a[,2], weights.class0 = a[,1], curve = TRUE)
  
  # Extract AUC-PR value
  auc_pr <- round(PR$auc.integral, digits = 2)
  
  
  # Cross validation stats
  stat=table(unlist(lapply(pile,function(x){(x$feat)})))
  stat=sort(stat,decreasing=TRUE)
  
  freq=cbind(data.frame(stat))
  save(pile,freq,f1_values,f1_score,balanced_acc_values,balanced_accuracy,a_df,auc_pr,file=paste("../Predict_LGA_new/10_proteins/Predict_LGA_T",trim,"_",modtype,"_10_proteins.RData",sep=""))
}


# 5 proteins --------------------------------------------------------------


rm(list=ls())

modtype="RF" #RF #EN
tpoints=c("1","2","3")
set.seed(1)
#trim = 1
for(trim in tpoints){
  
  load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
  MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$Bwkat %in% c("LGA", "AGA")]
  
  esetbig=exprs(MoM_exprset_STORK_LGA_SGA)
  anobig=pData(MoM_exprset_STORK_LGA_SGA)
  anobig <- anobig %>%
    dplyr::rename("bw_zscore" = `Birthweight z-score`)
  
  anosafe=anobig[anobig$TimePoint==trim,]
  esetsafe=esetbig[,rownames(anosafe)]
  
  all(rownames(anosafe)==colnames(esetsafe))
  anosafe$G=factor(ifelse(anosafe$Bwkat=="LGA","D","C"))
  #anosafe have G as "C", "D"
  
  pile=list()
  all_true_labels <- c()
  all_predicted_probs <- c()
  f1_values <- c() # Store F1 score for each LOOCV iteration
  balanced_acc_values <- c() # Store Balanced Accuracy for each LOOCV iteration
  
  # Leave-One-Out Cross-Validation
  for(ite in 1:nrow(anosafe)){
    cat(ite);cat("\n")
    #cat(1);cat("\n")
    ano=anosafe
    esetn=esetsafe
    
    ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[ite],"validation","discovery")
    #ano$subset_Bwkat=ifelse(rownames(ano)==rownames(ano)[1],"validation","discovery")
    anothis=ano
    
    ano_tr=ano[ano$subset_Bwkat=="discovery",]
    esetn_tr=esetn[,rownames(ano_tr)]
    
    #####prepare test data 
    ano_t=anothis[anothis$subset_Bwkat=="validation",]
    esetn_t=matrix(esetsafe[,rownames(ano_t)], ncol = 1)
    colnames(esetn_t) <- rownames(ano_t)
    rownames(esetn_t) <-rownames(esetn_tr)
    
    #select top proteins by limma test
    design <- model.matrix(~0+G,ano_tr) 
    colnames(design)<-gsub("G","",colnames(design))
    
    pe_v_control_cont = makeContrasts(
      D - C,
      levels = design)
    fit_noBayes = lmFit(esetn_tr, design)
    fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
    fit_pe_v_control = eBayes(fit_pe_v_control)
    deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
    sel=rownames(deT1)[1:5]
    
    #Choose proteins for reduced models
    train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
    test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
    rownames(test) <- colnames(esetn_t)
    
    #Models with all proteins
    #train=data.frame(G=factor(ano_tr$G),t(esetn_tr))
    #test=data.frame(G=factor(ano_t$G),t(esetn_t))
    
    
    if(modtype=="EN"){
      ####select genes by lasso
      #EN
      p_recipe <- recipe(G ~ ., data = train) %>% 
        step_zv(all_numeric(), -all_outcomes()) 
      
      lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      lasso_wf <- workflow() %>%
        add_recipe(p_recipe)
      
      lasso_fit <- lasso_wf %>%
        add_model(lasso) %>%
        fit(data = train)
      
      lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
        set_engine("glmnet") %>% 
        set_mode("classification") 
      
      # Define the search grid for hyperparameter tuning 
      lasso_grid <- grid_regular(penalty(), mixture(),levels = 10) 
      
      # Tune the hyperparameters using split sample
      lasso_res <- tune_grid(lasso_wf %>% 
                               add_model(lasso_tun), 
                             resamples=apparent(train),
                             grid = lasso_grid,
                             control=control_grid(parallel_over="everything", 
                                                  verbose = TRUE))
      
      #Tune hyperparameters by 5-fold CV
      #folds <- vfold_cv(train, v = 5)
      #lasso_res <- tune_grid(
      #  lasso_wf %>% add_model(lasso_tun), 
      #  resamples = folds,
      #  grid = lasso_grid,
      #  control = control_grid(parallel_over = "everything", verbose = TRUE)
      #)
      
      lowest_auc <- lasso_res %>%
        select_best("roc_auc")
      
      lasso_model_tun <- finalize_workflow(
        lasso_wf %>% add_model(lasso_tun),
        lowest_auc
      )%>%  fit(data = train)
      
      final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
        tidy()
      #
      ###Predictions
      lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
      #final$term[final$estimate!=0][-1]
      
    }else{
      modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
      lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
    }
    
    # Collect the true labels and predicted probabilities
    true_labels = as.numeric(test$G == "D")
    predicted_probs = lasso_pred_tun
    all_true_labels = c(all_true_labels, true_labels)
    all_predicted_probs = c(all_predicted_probs, predicted_probs)
    pred_binary = ifelse(predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
    
    # Confusion matrix components
    tp = sum(true_labels == 1 & pred_binary == 1)
    fp = sum(true_labels == 0 & pred_binary == 1)
    fn = sum(true_labels == 1 & pred_binary == 0)
    tn = sum(true_labels == 0 & pred_binary == 0)
    
    # Compute precision and recall
    precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    # Compute F1 score for this fold
    f1_score = ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    # Store F1 score for visualization
    f1_values[ite] <- f1_score
    
    
    # Compute sensitivity (recall for positive class)
    sensitivity = ifelse(sum(true_labels == 1) == 0, 0, sum(true_labels == 1 & pred_binary == 1) / sum(true_labels == 1))
    
    # Compute specificity (recall for negative class)
    specificity = ifelse(sum(true_labels == 0) == 0, 0, sum(true_labels == 0 & pred_binary == 0) / sum(true_labels == 0))
    
    # Compute balanced accuracy for this fold
    balanced_acc = (sensitivity + specificity) / 2
    
    # Store the balanced accuracy for visualization later
    balanced_acc_values[ite] <- balanced_acc
    
    #outcs=as.numeric(test$G=="D")
    #out=lasso_pred_tun
    #outcs = true_labels
    #out = predicted_probs
    outmat=cbind(true_labels,predicted_probs) #true outcomes for current test set side by side with predictited risk scores
    rownames(outmat) = rownames(test)
    sel <- rownames(esetn_tr)
    pile[[ite]]<-list(outmat=outmat,feat=sel)
    #pile[[1]]<-list(outmat=outmat,feat=sel)
    
  }
  
  #create the full matrix of predicted risk scores and true outcomes
  a=NULL
  for (i in 1:length(pile)){
    a=rbind(a,pile[[i]]$outmat)
  }
  
  # Compute overall F1 score across all iterations
  a_df <- as.data.frame(a)
  a_df$pred_binary = ifelse(a_df$predicted_probs > 0.5, 1, 0)  # Convert probabilities to binary predictions (threshold at 0.5)
  
  tp = sum(a_df$true_labels == 1 & a_df$pred_binary == 1)
  fp = sum(a_df$true_labels == 0 & a_df$pred_binary == 1)
  fn = sum(a_df$true_labels == 1 & a_df$pred_binary == 0)
  
  precision = ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall = ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  
  #Calculate F1 score
  if (precision + recall == 0) {
    f1_score <- 0
  } else {
    f1_score <- 2 * (precision * recall) / (precision + recall)
  }
  
  # Compute overall Balanced Accuracy across all iterations
  sensitivity_all = ifelse(sum(a_df$outcs == 1) == 0, 0, sum(a_df$outcs == 1 & a_df$pred_binary == 1) / sum(T1_EN_a_df$outcs == 1))
  specificity_all = ifelse(sum(a_df$outcs == 0) == 0, 0, sum(a_df$outcs == 0 & a_df$pred_binary == 0) / sum(T1_EN_a_df$outcs == 0))
  
  balanced_accuracy = (sensitivity_all + specificity_all) / 2
  
  
  # Compute Precision-Recall AUC
  PR <- pr.curve(scores.class0 = a[,2], weights.class0 = a[,1], curve = TRUE)
  
  # Extract AUC-PR value
  auc_pr <- round(PR$auc.integral, digits = 2)
  
  
  # Cross validation stats
  stat=table(unlist(lapply(pile,function(x){(x$feat)})))
  stat=sort(stat,decreasing=TRUE)
  
  freq=cbind(data.frame(stat))
  save(pile,freq,f1_values,f1_score,balanced_acc_values,balanced_accuracy,a_df,auc_pr,file=paste("../Predict_LGA_new/5_proteins/Predict_LGA_T",trim,"_",modtype,"_5_proteins.RData",sep=""))
}

