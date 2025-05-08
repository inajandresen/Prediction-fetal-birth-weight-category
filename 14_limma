library(tidyverse)
library(writexl)
library(readxl)
library(Biobase)
library(limma)
library(ggVennDiagram)
library(pcaMethods)
library(scales)
library(gridExtra)
library(ggforce)
library(ggrepel)
library(ggvenn)
library(ggpubr)


# Limma t-test, adjust nulliparity and BMI ------------------------------------------------------------

load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
samp <- pData(MoM_exprset_STORK_LGA_SGA)

#How many patients of each group?
samp_unique_ids <- samp %>%  #Extract one row per patient
  distinct(ID, .keep_all = TRUE)
table(unlist(samp_unique_ids$Bwkat))

rm(samp_unique_ids)


#Change the 1, 2, 3 to T1, T2, T3 in the time column
MoM_exprset_STORK_LGA_SGA$TimePoint[MoM_exprset_STORK_LGA_SGA$TimePoint == "1"] <- "T1"
MoM_exprset_STORK_LGA_SGA$TimePoint[MoM_exprset_STORK_LGA_SGA$TimePoint == "2"] <- "T2"
MoM_exprset_STORK_LGA_SGA$TimePoint[MoM_exprset_STORK_LGA_SGA$TimePoint == "3"] <- "T3"

#Merge health (preeclamptic or normal) with Time
Time_health <- as.factor(paste(MoM_exprset_STORK_LGA_SGA$Bwkat, MoM_exprset_STORK_LGA_SGA$TimePoint, sep = "."))

#Extract the sample information
samp <- pData(MoM_exprset_STORK_LGA_SGA)

#Make the design matrix and adjust for age, parity, smoking and BMI
design <- model.matrix(~0+Time_health+as.numeric(BMI)+as.factor(Nulliparity),samp)
colnames(design) <- c("AGA_T1", "AGA_T2", "AGA_T3", "LGA_T1", "LGA_T2", "LGA_T3", "SGA_T1", "SGA_T2", "SGA_T3", "BMI", "Nulliparity")


#See this about blocking: https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit <- duplicateCorrelation(MoM_exprset_STORK_LGA_SGA, design, block = MoM_exprset_STORK_LGA_SGA$ID)
corfit$consensus 
#0.5679379

#Do the limma test, and block for correlations effects from patients
fit <- lmFit(MoM_exprset_STORK_LGA_SGA, design, block = MoM_exprset_STORK_LGA_SGA$ID, correlation = corfit$consensus)

#Make the contrast: Which proteins are differentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast <- makeContrasts(T1_LGAvsAGA = LGA_T1 - AGA_T1, 
                          T2_LGAvsAGA = LGA_T2 - AGA_T2, 
                          T3_LGAvsAGA = LGA_T3 - AGA_T3,
                          T1_SGAvsAGA = SGA_T1 - AGA_T1, 
                          T2_SGAvsAGA = SGA_T2 - AGA_T2, 
                          T3_SGAvsAGA = SGA_T3 - AGA_T3,
                          T1_LGAvsSGA = LGA_T1 - SGA_T1, 
                          T2_LGAvsSGA = LGA_T2 - SGA_T2, 
                          T3_LGAvsSGA = LGA_T3 - SGA_T3,
                          levels = design)


#Compute the moderated t-tests
fit <- contrasts.fit(fit, contrast)
efit <- eBayes(fit)

save(fit, design, efit, file = "../Limma_groups/Adjust_BMI_nulli/efit_LGA_AGA_SGA_adj_nulli_bmi.RData")
load("../Limma_groups/Adjust_BMI_nulli/efit_LGA_AGA_SGA_adj_nulli_bmi.RData")

#Negative fold change means higher in AGA
#Positive fold change means higher in LGA

#No cut off values
T1_LGAvsAGA <- topTable(efit, coef = "T1_LGAvsAGA", number = "all", adjust.method="BH")
T2_LGAvsAGA <- topTable(efit, coef = "T2_LGAvsAGA", number = "all", adjust.method="BH")
T3_LGAvsAGA <- topTable(efit, coef = "T3_LGAvsAGA", number = "all", adjust.method="BH")

T1_SGAvsAGA <- topTable(efit, coef = "T1_SGAvsAGA", number = "all", adjust.method="BH")
T2_SGAvsAGA <- topTable(efit, coef = "T2_SGAvsAGA", number = "all", adjust.method="BH")
T3_SGAvsAGA <- topTable(efit, coef = "T3_SGAvsAGA", number = "all", adjust.method="BH")

T1_LGAvsSGA <- topTable(efit, coef = "T1_LGAvsSGA", number = "all", adjust.method="BH")
T2_LGAvsSGA <- topTable(efit, coef = "T2_LGAvsSGA", number = "all", adjust.method="BH")
T3_LGAvsSGA <- topTable(efit, coef = "T3_LGAvsSGA", number = "all", adjust.method="BH")


#Save LGA vs AGA
T1_LGAvsAGA <- T1_LGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T2_LGAvsAGA <- T2_LGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T3_LGAvsAGA <- T3_LGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T1_LGAvsAGA, path = "../Limma_groups/Adjust_BMI_nulli/T1_LGA_AGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)
write_xlsx(T2_LGAvsAGA, path = "..//Limma_groups/Adjust_BMI_nulli/T2_LGA_AGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)
write_xlsx(T3_LGAvsAGA, path = "../Limma_groups/Adjust_BMI_nulli/T3_LGA_AGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)

T1_LGAvsAGA_sig_prot <- T1_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T1_LGAvsAGA_sig_prot)

T2_LGAvsAGA_sig_prot <- T2_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T2_LGAvsAGA_sig_prot)

T3_LGAvsAGA_sig_prot <- T3_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T3_LGAvsAGA_sig_prot)

#Save SGA vs AGA
T1_SGAvsAGA <- T1_SGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T2_SGAvsAGA <- T2_SGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T3_SGAvsAGA <- T3_SGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T1_SGAvsAGA, path = "../Limma_groups/Adjust_BMI_nulli/T1_SGA_AGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)
write_xlsx(T2_SGAvsAGA, path = "../Limma_groups/Adjust_BMI_nulli/T2_SGA_AGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)
write_xlsx(T3_SGAvsAGA, path = "../Limma_groups/Adjust_BMI_nulli/T3_SGA_AGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)

T1_SGAvsAGA_sig_prot <- T1_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T1_SGAvsAGA_sig_prot)

T2_SGAvsAGA_sig_prot <- T2_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T2_SGAvsAGA_sig_prot)

T3_SGAvsAGA_sig_prot <- T3_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T3_SGAvsAGA_sig_prot)

#Save LGA vs SGA
T1_LGAvsSGA <- T1_LGAvsSGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T2_LGAvsSGA <- T2_LGAvsSGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T3_LGAvsSGA <- T3_LGAvsSGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T1_LGAvsSGA, path = "../Limma_groups/Adjust_BMI_nulli/T1_LGA_SGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)
write_xlsx(T2_LGAvsSGA, path = "../Limma_groups/Adjust_BMI_nulli/T2_LGA_SGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)
write_xlsx(T3_LGAvsSGA, path = "../Limma_groups/Adjust_BMI_nulli/T3_LGA_SGA_adj_nulli_bmi.xlsx",
           col_names = TRUE)


T1_LGAvsSGA_sig_prot <- T1_LGAvsSGA %>%
  filter(adj.P.Val < 0.05)
dim(T1_LGAvsSGA_sig_prot)

T2_LGAvsSGA_sig_prot <- T2_LGAvsSGA %>%
  filter(adj.P.Val < 0.05)
dim(T2_LGAvsSGA_sig_prot)

T3_LGAvsSGA_sig_prot <- T3_LGAvsSGA %>%
  filter(adj.P.Val < 0.05)
dim(T3_LGAvsSGA_sig_prot)


# Limma t-test, no adjustment ------------------------------------------------------------

load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
samp <- pData(MoM_exprset_STORK_LGA_SGA)

#How many patients of each group?
samp_unique_ids <- samp %>%  #Extract one row per patient
  distinct(ID, .keep_all = TRUE)
table(unlist(samp_unique_ids$Bwkat))

rm(samp_unique_ids)


#Change the 1, 2, 3 to T1, T2, T3 in the time column
MoM_exprset_STORK_LGA_SGA$TimePoint[MoM_exprset_STORK_LGA_SGA$TimePoint == "1"] <- "T1"
MoM_exprset_STORK_LGA_SGA$TimePoint[MoM_exprset_STORK_LGA_SGA$TimePoint == "2"] <- "T2"
MoM_exprset_STORK_LGA_SGA$TimePoint[MoM_exprset_STORK_LGA_SGA$TimePoint == "3"] <- "T3"

#Merge health (preeclamptic or normal) with Time
Time_health <- as.factor(paste(MoM_exprset_STORK_LGA_SGA$Bwkat, MoM_exprset_STORK_LGA_SGA$TimePoint, sep = "."))

#Extract the sample information
samp <- pData(MoM_exprset_STORK_LGA_SGA)

#Make the design matrix and adjust for age, parity, smoking and BMI
design <- model.matrix(~0+Time_health,samp)
colnames(design) <- c("AGA_T1", "AGA_T2", "AGA_T3", "LGA_T1", "LGA_T2", "LGA_T3", "SGA_T1", "SGA_T2", "SGA_T3")


#See this about blocking: https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit <- duplicateCorrelation(MoM_exprset_STORK_LGA_SGA, design, block = MoM_exprset_STORK_LGA_SGA$ID)
corfit$consensus 

#Do the limma test, and block for correlations effects from patients
fit <- lmFit(MoM_exprset_STORK_LGA_SGA, design, block = MoM_exprset_STORK_LGA_SGA$ID, correlation = corfit$consensus)

#Make the contrast: Which proteins are deferentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast <- makeContrasts(T1_LGAvsAGA = LGA_T1 - AGA_T1, 
                          T2_LGAvsAGA = LGA_T2 - AGA_T2, 
                          T3_LGAvsAGA = LGA_T3 - AGA_T3,
                          T1_SGAvsAGA = SGA_T1 - AGA_T1, 
                          T2_SGAvsAGA = SGA_T2 - AGA_T2, 
                          T3_SGAvsAGA = SGA_T3 - AGA_T3,
                          T1_LGAvsSGA = LGA_T1 - SGA_T1, 
                          T2_LGAvsSGA = LGA_T2 - SGA_T2, 
                          T3_LGAvsSGA = LGA_T3 - SGA_T3,
                          levels = design)


#Compute the moderated t-tests
fit <- contrasts.fit(fit, contrast)
efit <- eBayes(fit)

save(fit, design, efit, file = "../Limma_groups/No_adjustments/efit_LGA_AGA_SGA_noadj.RData")

#Negative fold change means higher in AGA
#Positive fold change means higher in LGA

#No cut off values
T1_LGAvsAGA <- topTable(efit, coef = "T1_LGAvsAGA", number = "all", adjust.method="BH")
T2_LGAvsAGA <- topTable(efit, coef = "T2_LGAvsAGA", number = "all", adjust.method="BH")
T3_LGAvsAGA <- topTable(efit, coef = "T3_LGAvsAGA", number = "all", adjust.method="BH")

T1_SGAvsAGA <- topTable(efit, coef = "T1_SGAvsAGA", number = "all", adjust.method="BH")
T2_SGAvsAGA <- topTable(efit, coef = "T2_SGAvsAGA", number = "all", adjust.method="BH")
T3_SGAvsAGA <- topTable(efit, coef = "T3_SGAvsAGA", number = "all", adjust.method="BH")

T1_LGAvsSGA <- topTable(efit, coef = "T1_LGAvsSGA", number = "all", adjust.method="BH")
T2_LGAvsSGA <- topTable(efit, coef = "T2_LGAvsSGA", number = "all", adjust.method="BH")
T3_LGAvsSGA <- topTable(efit, coef = "T3_LGAvsSGA", number = "all", adjust.method="BH")


#Save LGA vs AGA
T1_LGAvsAGA <- T1_LGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T2_LGAvsAGA <- T2_LGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T3_LGAvsAGA <- T3_LGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T1_LGAvsAGA, path = "../Limma_groups/No_adjustments/T1_LGA_AGA_noadj.xlsx",
           col_names = TRUE)
write_xlsx(T2_LGAvsAGA, path = "../Limma_groups/No_adjustments/T2_LGA_AGA_noadj.xlsx",
           col_names = TRUE)
write_xlsx(T3_LGAvsAGA, path = "../Limma_groups/No_adjustments/T3_LGA_AGA_noadj.xlsx",
           col_names = TRUE)

T1_LGAvsAGA_sig_prot <- T1_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T1_LGAvsAGA_sig_prot)

T2_LGAvsAGA_sig_prot <- T2_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T2_LGAvsAGA_sig_prot)

T3_LGAvsAGA_sig_prot <- T3_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T3_LGAvsAGA_sig_prot)

#Save SGA vs AGA 
T1_SGAvsAGA <- T1_SGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T2_SGAvsAGA <- T2_SGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T3_SGAvsAGA <- T3_SGAvsAGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T1_SGAvsAGA, path = "../Limma_groups/No_adjustments/T1_SGA_AGA_noadj.xlsx",
           col_names = TRUE)
write_xlsx(T2_SGAvsAGA, path = "../Limma_groups/No_adjustments/T2_SGA_AGA_noadj.xlsx",
           col_names = TRUE)
write_xlsx(T3_SGAvsAGA, path = "../Limma_groups/No_adjustments/T3_SGA_AGA_noadj.xlsx",
           col_names = TRUE)

T1_SGAvsAGA_sig_prot <- T1_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T1_SGAvsAGA_sig_prot)

T2_SGAvsAGA_sig_prot <- T2_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T2_SGAvsAGA_sig_prot)
 
T3_SGAvsAGA_sig_prot <- T3_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
dim(T3_SGAvsAGA_sig_prot)


#Save LGA vs SGA
T1_LGAvsSGA <- T1_LGAvsSGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T2_LGAvsSGA <- T2_LGAvsSGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
T3_LGAvsSGA <- T3_LGAvsSGA %>%
  dplyr::select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

write_xlsx(T1_LGAvsSGA, path = "../Limma_groups/No_adjustments/T1_LGA_SGA_noadj.xlsx",
           col_names = TRUE)
write_xlsx(T2_LGAvsSGA, path = "../Limma_groups/No_adjustments/T2_LGA_SGA_noadj.xlsx",
           col_names = TRUE)
write_xlsx(T3_LGAvsSGA, path = "../Limma_groups/No_adjustments/T3_LGA_SGA_noadj.xlsx",
           col_names = TRUE)

T1_LGAvsSGA_sig_prot <- T1_LGAvsSGA %>%
  filter(adj.P.Val < 0.05)
dim(T1_LGAvsSGA_sig_prot)

T2_LGAvsSGA_sig_prot <- T2_LGAvsSGA %>%
  filter(adj.P.Val < 0.05)
dim(T2_LGAvsSGA_sig_prot)

T3_LGAvsSGA_sig_prot <- T3_LGAvsSGA %>%
  filter(adj.P.Val < 0.05)
dim(T3_LGAvsSGA_sig_prot)

