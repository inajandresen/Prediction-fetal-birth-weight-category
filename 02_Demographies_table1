library(tidyverse)
library(Biobase)
library(reshape)
library(haven)




load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
MOM_exprset_STORK_LGA_SGA <- MOM_exprset_STORK_LGA_SGA[, MOM_exprset_STORK_LGA_SGA$PE %in% "0"]
samp <- pData(MOM_exprset_STORK_LGA_SGA)


#Add clinical factors 
#Find out which patient has missing BMI (has been imputed and needs to be imputet in another way)
STORK.cl <- read_sav("file.sav")
colnames(STORK.cl)[1] <- "ID"
colnames(STORK.cl)[6] <- "Weightgain_1_2"
colnames(STORK.cl)[7] <- "Weightgain_1_3"
colnames(STORK.cl)[8] <- "Weightgain_1_4"

STORK.cl2 <- STORK.cl %>%
  dplyr::select("ID", "alder", "bmi1", "Weightgain_1_2", "Weightgain_1_3", "Weightgain_1_4", "sivilst", "utd", "yrkesakt", 
                "bydel", "royk", "daglig_royk", "fvektmor", "fvektfar", "ivf_II", "hoyde", "pregravid_vekt_II", "vekt1", 
                "vekt2", "vekt3", "vekt4", "uke_dag1", "uke_dag4", "lengdeb", "gest_dia_I")

samp2 <- merge(samp, STORK.cl2, by = "ID")

#Remove dupicated patients 
samp_unique <- samp2 %>%  #Extract one row per patient
  distinct(ID, .keep_all = TRUE)
samp_unique_SGA <- samp_unique[samp_unique$Bwkat == "SGA", ]
samp_unique_AGA <- samp_unique[samp_unique$Bwkat == "AGA", ]
samp_unique_LGA <- samp_unique[samp_unique$Bwkat == "LGA", ]

#How many of each group 
table(unlist(samp_unique$Bwkat))
#AGA 58
#LGA 5
#SGA 7


#Maternal age -------------------------------------------------------------------------

mean(samp_unique_SGA$Age)
sd(samp_unique_SGA$Age)
mean(samp_unique_AGA$Age)
sd(samp_unique_AGA$Age)
mean(samp_unique_LGA$Age)
sd(samp_unique_LGA$Age)

t.test(x = samp_unique_SGA$Age, y = samp_unique_AGA$Age)
t.test(x = samp_unique_LGA$Age, y = samp_unique_AGA$Age)
t.test(x = samp_unique_LGA$Age, y = samp_unique_SGA$Age)

# Maternal heigth ---------------------------------------------------------

mean(samp_unique_SGA$hoyde)
sd(samp_unique_SGA$hoyde)
mean(samp_unique_AGA$hoyde)
sd(samp_unique_AGA$hoyde)
mean(samp_unique_LGA$hoyde)
sd(samp_unique_LGA$hoyde)

t.test(x = samp_unique_SGA$hoyde, y = samp_unique_AGA$hoyde)
t.test(x = samp_unique_LGA$hoyde, y = samp_unique_AGA$hoyde)
t.test(x = samp_unique_LGA$hoyde, y = samp_unique_SGA$hoyde)


# Maternal weight ---------------------------------------------------------

mean(samp_unique_SGA$vekt1)
sd(samp_unique_SGA$vekt1)
mean(samp_unique_AGA$vekt1, na.rm = TRUE)
sd(samp_unique_AGA$vekt1, na.rm = TRUE)
mean(samp_unique_LGA$vekt1)
sd(samp_unique_LGA$vekt1)

t.test(x = samp_unique_SGA$vekt1, y = samp_unique_AGA$vekt1)
t.test(x = samp_unique_LGA$vekt1, y = samp_unique_AGA$vekt1)
t.test(x = samp_unique_LGA$vekt1, y = samp_unique_SGA$vekt1)



# Maternal BMI ------------------------------------------------------------
#Based on weight at visit 1 

mean(samp_unique_SGA$BMI)
sd(samp_unique_SGA$BMI)
mean(samp_unique_AGA$BMI)
sd(samp_unique_AGA$BMI)
mean(samp_unique_LGA$BMI)
sd(samp_unique_LGA$BMI)

t.test(x = samp_unique_SGA$BMI, y = samp_unique_AGA$BMI)
t.test(x = samp_unique_LGA$BMI, y = samp_unique_AGA$BMI)
t.test(x = samp_unique_LGA$BMI, y = samp_unique_SGA$BMI)


# Maternal weight gain ----------------------------------------------------
test_weightgain <- samp_unique %>%
  dplyr::select(Weightgain_1_4, vekt1, vekt4, uke_dag1, uke_dag4)

mean(samp_unique_SGA$Weightgain_1_4)
sd(samp_unique_SGA$Weightgain_1_4)
mean(samp_unique_AGA$Weightgain_1_4, na.rm = TRUE)
sd(samp_unique_AGA$Weightgain_1_4, na.rm = TRUE)
mean(samp_unique_LGA$Weightgain_1_4, na.rm = TRUE)
sd(samp_unique_LGA$Weightgain_1_4, na.rm = TRUE)

t.test(x = samp_unique_SGA$Weightgain_1_4, y = samp_unique_AGA$Weightgain_1_4)
t.test(x = samp_unique_LGA$Weightgain_1_4, y = samp_unique_AGA$Weightgain_1_4)
t.test(x = samp_unique_LGA$Weightgain_1_4, y = samp_unique_SGA$Weightgain_1_4)

# Placental weight --------------------------------------------------------

mean(samp_unique_SGA$PlWeigth, na.rm = TRUE)
sd(samp_unique_SGA$PlWeigth, na.rm = TRUE)
mean(samp_unique_AGA$PlWeigth, na.rm = TRUE)
sd(samp_unique_AGA$PlWeigth, na.rm = TRUE)
mean(samp_unique_LGA$PlWeigth, na.rm = TRUE)
sd(samp_unique_LGA$PlWeigth, na.rm = TRUE)

t.test(x = samp_unique_SGA$PlWeigth, y = samp_unique_AGA$PlWeigth)
t.test(x = samp_unique_LGA$PlWeigth, y = samp_unique_AGA$PlWeigth)
t.test(x = samp_unique_LGA$PlWeigth, y = samp_unique_SGA$PlWeigth)


# Nulliparity -------------------------------------------------------------


Nulliparity_table <- table(samp_unique$Bwkat, samp_unique$Nulliparity)
Nulliparity_table <- as.data.frame(Nulliparity_table)
colnames(Nulliparity_table) <- c("Bwkat", "Nulliparity", "Freq")
Nulliparity_table <- cast(Nulliparity_table, Nulliparity~Bwkat)
Nulliparity_table <- Nulliparity_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Nulliparity")
#   AGA LGA SGA
#0  29   4   2
#1  29   1   5

#SGA percentage
(5/7)*100
(2/7)*100

#AGA percentage
(29/58)*100
(29/58)*100

#LGA percetage
(1/5)*100
(4/5)*100

#SGA vs AGA
Nulliparity_table_SGA_AGA <- Nulliparity_table %>%
  dplyr::select(SGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(Nulliparity_table_SGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_nulliparity_SGA_AGA <- fisher.test(Nulliparity_table_SGA_AGA)
fisher_nulliparity_SGA_AGA

#LGA vs AGA
Nulliparity_table_LGA_AGA <- Nulliparity_table %>%
  dplyr::select(LGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(Nulliparity_table_LGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_nulliparity_LGA_AGA <- fisher.test(Nulliparity_table_LGA_AGA)
fisher_nulliparity_LGA_AGA

#LGA vs SGA
Nulliparity_table_LGA_SGA <- Nulliparity_table %>%
  dplyr::select(LGA, SGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(Nulliparity_table_LGA_SGA)$expected
#No
#Do the fisher's exact test 
fisher_nulliparity_LGA_SGA <- fisher.test(Nulliparity_table_LGA_SGA)
fisher_nulliparity_LGA_SGA



# Smoking yes vs no -----------------------------------------------------------------
#Royk
#0 = r?yker ikke 
#1 = R?yker >1=<10
#2 = sluttet i svangerskapet
#3 = r?yer >10

#Make "r?yker ikke" and "sluttet i svangerskapet) = 1 and "R?yker >1=<10" and "r?yer >10" = 0
samp_unique$Smoking <- as.character(samp_unique$royk)
samp_unique <- samp_unique %>%
  dplyr::mutate(Smoking = ifelse(Smoking == "1", "1", 0))
test2 <- samp_unique %>%
  dplyr::select(Smoking, royk)

smoking_table <- table(samp_unique$Bwkat, samp_unique$Smoking)
smoking_table <- as.data.frame(smoking_table)
colnames(smoking_table) <- c("Bwkat", "Smoking", "Freq")
smoking_table <- cast(smoking_table, Smoking~Bwkat)
smoking_table <- smoking_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Smoking")
#   AGA LGA SGA
#0  57   5   7
#1   1   0   0

#SGA percentage
(0/7)*100
(7/7)*100

#AGA percentage
(1/58)*100
(57/58)*100

#LGA percetage
(0/5)*100
(5/5)*100

#SGA vs AGA
smoking_table_SGA_AGA <- smoking_table %>%
  dplyr::select(SGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(smoking_table_SGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_smoking_SGA_AGA <- fisher.test(smoking_table_SGA_AGA)
fisher_smoking_SGA_AGA

#LGA vs AGA
smoking_table_LGA_AGA <- smoking_table %>%
  dplyr::select(LGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(smoking_table_LGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_smoking_LGA_AGA <- fisher.test(smoking_table_LGA_AGA)
fisher_smoking_LGA_AGA

#LGA vs SGA
smoking_table_LGA_SGA <- smoking_table %>%
  dplyr::select(LGA, SGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(smoking_table_LGA_SGA)$expected
#No
#Do the fisher's exact test 
fisher_smoking_LGA_SGA <- fisher.test(smoking_table_LGA_SGA)
fisher_smoking_LGA_SGA

# Smoking yes, quit, no -----------------------------------------------------------------
#Royk
#0 = r?yker ikke 
#1 = R?yker >1=<10
#2 = sluttet i svangerskapet
#3 = r?yer >10

#Make "r?yker ikke" and "sluttet i svangerskapet) = 1 and "R?yker >1=<10" and "r?yer >10" = 0
samp_unique$Smoking2 <- as.character(samp_unique$royk)

smoking_table2 <- table(samp_unique$Bwkat, samp_unique$Smoking2)
smoking_table2 <- as.data.frame(smoking_table2)
colnames(smoking_table2) <- c("Bwkat", "Smoking", "Freq")
smoking_table2 <- cast(smoking_table2, Smoking~Bwkat)
smoking_table2 <- smoking_table2 %>%
  remove_rownames() %>%
  column_to_rownames(var = "Smoking")
#   AGA LGA SGA
#0  45   5   4
#1   1   0   0
#2   12  0   3

#SGA percentage
(4/7)*100 #No
(0/7)*100 #Yes
(3/7)*100 #Quit

#AGA percentage
(45/58)*100 #No
(1/58)*100 #Yes
(12/58)*100 #Quit

#LGA percetage
(5/5)*100 #No
(0/5)*100 #Yes
(0/5)*100 #Quit

#SGA vs AGA
smoking_table_SGA_AGA <- smoking_table2 %>%
  dplyr::select(SGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(smoking_table_SGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_smoking_SGA_AGA <- fisher.test(smoking_table_SGA_AGA)
fisher_smoking_SGA_AGA

#LGA vs AGA
smoking_table_LGA_AGA <- smoking_table2 %>%
  dplyr::select(LGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(smoking_table_LGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_smoking_LGA_AGA <- fisher.test(smoking_table_LGA_AGA)
fisher_smoking_LGA_AGA

#LGA vs SGA
smoking_table_LGA_SGA <- smoking_table2 %>%
  dplyr::select(LGA, SGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(smoking_table_LGA_SGA)$expected
#No
#Do the fisher's exact test 
fisher_smoking_LGA_SGA <- fisher.test(smoking_table_LGA_SGA)
fisher_smoking_LGA_SGA



# Education ---------------------------------------------------------------
#Utd
#1 = grunnskole
#2 = videreg?ende
#3 = H?yere utdannelse

#Make "h?yere utdannelse" = 1 and "grunnskole" and "videreg?ende" = 0
samp_unique$utd <- as.character(samp_unique$utd)
samp_unique <- samp_unique %>%
  dplyr::mutate(edu = ifelse(utd == "3", "1", 0))

edu_table <- table(samp_unique$Bwkat, samp_unique$edu)
edu_table <- as.data.frame(edu_table)
colnames(edu_table) <- c("Bwkat", "Education", "Freq")
edu_table <- cast(edu_table, Education~Bwkat)
edu_table <- edu_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Education")
#   AGA LGA SGA
#0  8   1   3
#1  50  4   4

#SGA percentage
(4/7)*100
(3/7)*100

#AGA percentage
(50/58)*100
(8/58)*100

#LGA percetage
(4/5)*100
(1/5)*100

#SGA vs AGA
edu_table_SGA_AGA <- edu_table %>%
  dplyr::select(SGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(edu_table_SGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_edu_SGA_AGA <- fisher.test(edu_table_SGA_AGA)
fisher_edu_SGA_AGA

#LGA vs AGA
edu_table_LGA_AGA <- edu_table %>%
  dplyr::select(LGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(edu_table_LGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_edu_LGA_AGA <- fisher.test(edu_table_LGA_AGA)
fisher_edu_LGA_AGA

#LGA vs SGA
edu_table_LGA_SGA <- edu_table %>%
  dplyr::select(LGA, SGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(edu_table_LGA_SGA)$expected
#No
#Do the fisher's exact test 
fisher_edu_LGA_SGA <- fisher.test(edu_table_LGA_SGA)
fisher_edu_LGA_SGA


# Martial status ----------------------------------------------------------
#Sivil status 
#1 = gift
#2 = samboer
#3 = ugift/enslig

samp_unique$M_stat <- as.character(samp_unique$sivilst)
samp_unique <- samp_unique %>%
  dplyr::mutate(M_stat = ifelse(M_stat == "3", "0", 1))
test_mstatus <- samp_unique %>%
  dplyr::select(M_stat, sivilst)

M_stat_table <- table(samp_unique$Bwkat, samp_unique$M_stat)
M_stat_table <- as.data.frame(M_stat_table)
colnames(M_stat_table) <- c("Bwkat", "M_stat", "Freq")
M_stat_table <- cast(M_stat_table, M_stat~Bwkat)
M_stat_table <- M_stat_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "M_stat")
#   AGA LGA SGA
#0   0   0   0
#1  58   5   7


#SGA vs AGA
M_stat_table_SGA_AGA <- M_stat_table %>%
  dplyr::select(SGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(M_stat_table_SGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_M_stat_SGA_AGA <- fisher.test(M_stat_table_SGA_AGA)
fisher_M_stat_SGA_AGA

#LGA vs AGA
M_stat_table_LGA_AGA <- M_stat_table %>%
  dplyr::select(LGA, AGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(M_stat_table_LGA_AGA)$expected
#No
#Do the fisher's exact test 
fisher_M_stat_LGA_AGA <- fisher.test(M_stat_table_LGA_AGA)
fisher_M_stat_LGA_AGA

#LGA vs SGA
M_stat_table_LGA_SGA <- M_stat_table %>%
  dplyr::select(LGA, SGA)
#Get the expected frequencies (can we do a chi squared test?)
chisq.test(M_stat_table_LGA_SGA)$expected
#No
#Do the fisher's exact test 
fisher_M_stat_LGA_SGA <- fisher.test(M_stat_table_LGA_SGA)
fisher_M_stat_LGA_SGA



# Gestational age at visit 1 ----------------------------------------------

T1_samp <- samp2 %>%
  filter(TimePoint == 1)

T1_samp_SGA <- T1_samp %>%
  filter(Bwkat == "SGA")
T1_samp_AGA <- T1_samp %>%
  filter(Bwkat == "AGA")
T1_samp_LGA <- T1_samp %>%
  filter(Bwkat == "LGA")

mean(T1_samp_SGA$GAWeeks, na.rm = TRUE)
sd(T1_samp_SGA$GAWeeks, na.rm = TRUE)
mean(T1_samp_AGA$GAWeeks, na.rm = TRUE)
sd(T1_samp_AGA$GAWeeks, na.rm = TRUE)
mean(T1_samp_LGA$GAWeeks, na.rm = TRUE)
sd(T1_samp_LGA$GAWeeks, na.rm = TRUE)

t.test(x = T1_samp_SGA$GAWeeks, y = T1_samp_AGA$GAWeeks)
t.test(x = T1_samp_LGA$GAWeeks, y = T1_samp_AGA$GAWeeks)
t.test(x = T1_samp_LGA$GAWeeks, y = T1_samp_SGA$GAWeeks)

# Gestational age at visit 2 ----------------------------------------------

T2_samp <- samp2 %>%
  filter(TimePoint == 2)

T2_samp_SGA <- T2_samp %>%
  filter(Bwkat == "SGA")
T2_samp_AGA <- T2_samp %>%
  filter(Bwkat == "AGA")
T2_samp_LGA <- T2_samp %>%
  filter(Bwkat == "LGA")

mean(T2_samp_SGA$GAWeeks, na.rm = TRUE)
sd(T2_samp_SGA$GAWeeks, na.rm = TRUE)
mean(T2_samp_AGA$GAWeeks, na.rm = TRUE)
sd(T2_samp_AGA$GAWeeks, na.rm = TRUE)
mean(T2_samp_LGA$GAWeeks, na.rm = TRUE)
sd(T2_samp_LGA$GAWeeks, na.rm = TRUE)

t.test(x = T2_samp_SGA$GAWeeks, y = T2_samp_AGA$GAWeeks)
t.test(x = T2_samp_LGA$GAWeeks, y = T2_samp_AGA$GAWeeks)
t.test(x = T2_samp_LGA$GAWeeks, y = T2_samp_SGA$GAWeeks)

# Gestational age at visit 3 ----------------------------------------------

T3_samp <- samp2 %>%
  filter(TimePoint == 3)

T3_samp_SGA <- T3_samp %>%
  filter(Bwkat == "SGA")
T3_samp_AGA <- T3_samp %>%
  filter(Bwkat == "AGA")
T3_samp_LGA <- T3_samp %>%
  filter(Bwkat == "LGA")

mean(T3_samp_SGA$GAWeeks, na.rm = TRUE)
sd(T3_samp_SGA$GAWeeks, na.rm = TRUE)
mean(T3_samp_AGA$GAWeeks, na.rm = TRUE)
sd(T3_samp_AGA$GAWeeks, na.rm = TRUE)
mean(T3_samp_LGA$GAWeeks, na.rm = TRUE)
sd(T3_samp_LGA$GAWeeks, na.rm = TRUE)

t.test(x = T3_samp_SGA$GAWeeks, y = T3_samp_AGA$GAWeeks)
t.test(x = T3_samp_LGA$GAWeeks, y = T3_samp_AGA$GAWeeks)
t.test(x = T3_samp_LGA$GAWeeks, y = T3_samp_SGA$GAWeeks)

#Fetal birth weight -----------------------------------------------------

mean(samp_unique_SGA$`Birth weight`, na.rm = TRUE)
sd(samp_unique_SGA$`Birth weight`, na.rm = TRUE)
mean(samp_unique_AGA$`Birth weight`, na.rm = TRUE)
sd(samp_unique_AGA$`Birth weight`, na.rm = TRUE)
mean(samp_unique_LGA$`Birth weight`, na.rm = TRUE)
sd(samp_unique_LGA$`Birth weight`, na.rm = TRUE)

t.test(x = samp_unique_SGA$`Birth weight`, y = samp_unique_AGA$`Birth weight`)
t.test(x = samp_unique_LGA$`Birth weight`, y = samp_unique_AGA$`Birth weight`)
t.test(x = samp_unique_LGA$`Birth weight`, y = samp_unique_SGA$`Birth weight`)


#Fetal birth weight (z-score) -----------------------------------------------------

mean(samp_unique_SGA$`Birthweight z-score`, na.rm = TRUE)
sd(samp_unique_SGA$`Birthweight z-score`, na.rm = TRUE)
mean(samp_unique_AGA$`Birthweight z-score`, na.rm = TRUE)
sd(samp_unique_AGA$`Birthweight z-score`, na.rm = TRUE)
mean(samp_unique_LGA$`Birthweight z-score`, na.rm = TRUE)
sd(samp_unique_LGA$`Birthweight z-score`, na.rm = TRUE)

t.test(x = samp_unique_SGA$`Birthweight z-score`, y = samp_unique_AGA$`Birthweight z-score`)
t.test(x = samp_unique_LGA$`Birthweight z-score`, y = samp_unique_AGA$`Birthweight z-score`)
t.test(x = samp_unique_LGA$`Birthweight z-score`, y = samp_unique_SGA$`Birthweight z-score`)

# Fetal lenght ------------------------------------------------------------

mean(samp_unique_SGA$lengdeb, na.rm = TRUE)
sd(samp_unique_SGA$lengdeb, na.rm = TRUE)
mean(samp_unique_AGA$lengdeb, na.rm = TRUE)
sd(samp_unique_AGA$lengdeb, na.rm = TRUE)
mean(samp_unique_LGA$lengdeb, na.rm = TRUE)
sd(samp_unique_LGA$lengdeb, na.rm = TRUE)

t.test(x = samp_unique_SGA$lengdeb, y = samp_unique_AGA$lengdeb)
t.test(x = samp_unique_LGA$lengdeb, y = samp_unique_AGA$lengdeb)
t.test(x = samp_unique_LGA$lengdeb, y = samp_unique_SGA$lengdeb)


# Fetal weight/Placental weight ratio (placenta efficiency) ---------------

samp_unique_SGA$fw_pw_ratio <- samp_unique_SGA$`Birth weight`/samp_unique_SGA$PlWeigth
mean(samp_unique_SGA$fw_pw_ratio, na.rm = TRUE)
sd(samp_unique_SGA$fw_pw_ratio, na.rm = TRUE)
samp_unique_AGA$fw_pw_ratio <- samp_unique_AGA$`Birth weight`/samp_unique_AGA$PlWeigth
mean(samp_unique_AGA$fw_pw_ratio, na.rm = TRUE)
sd(samp_unique_AGA$fw_pw_ratio, na.rm = TRUE)
samp_unique_LGA$fw_pw_ratio <- samp_unique_LGA$`Birth weight`/samp_unique_LGA$PlWeigth
mean(samp_unique_LGA$fw_pw_ratio, na.rm = TRUE)
sd(samp_unique_LGA$fw_pw_ratio, na.rm = TRUE)

t.test(x = samp_unique_SGA$fw_pw_ratio, y = samp_unique_AGA$fw_pw_ratio)
t.test(x = samp_unique_LGA$fw_pw_ratio, y = samp_unique_AGA$fw_pw_ratio)
t.test(x = samp_unique_LGA$fw_pw_ratio, y = samp_unique_SGA$fw_pw_ratio)


# Fetal gender  -----------------------------------------------------------
#kjonn
#1 = jente
#2 = gutt

#1 = female, 0 = male
Gender_table <- table(samp_unique$Bwkat, samp_unique$Gender)
Gender_table <- as.data.frame(Gender_table)
colnames(Gender_table) <- c("Bwkat", "Gender", "Freq")
Gender_table <- cast(Gender_table, Gender~Bwkat)
Gender_table <- Gender_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Gender")

#   AGA LGA SGA
#0  33   4   3
#1  25   1   4

#SGA percentage
(4/7)*100
(3/7)*100

#AGA percentage
(25/58)*100
(33/58)*100

#LGA percetage
(1/5)*100
(4/5)*100

#SGA vs AGA
Gender_table_SGA_AGA <- Gender_table %>%
  dplyr::select(SGA, AGA)
#Get the exGendercted frequencies (can we do a chi squared test?)
chisq.test(Gender_table_SGA_AGA)$exGendercted
#No
#Do the fisher's exact test 
fisher_Gender_SGA_AGA <- fisher.test(Gender_table_SGA_AGA)
fisher_Gender_SGA_AGA


#LGA vs AGA
Gender_table_LGA_AGA <- Gender_table %>%
  dplyr::select(LGA, AGA)
#Get the exGendercted frequencies (can we do a chi squared test?)
chisq.test(Gender_table_LGA_AGA)$exGendercted
#No
#Do the fisher's exact test 
fisher_Gender_LGA_AGA <- fisher.test(Gender_table_LGA_AGA)
fisher_Gender_LGA_AGA

#LGA vs SGA
Gender_table_LGA_SGA <- Gender_table %>%
  dplyr::select(LGA, SGA)
#Get the exGendercted frequencies (can we do a chi squared test?)
chisq.test(Gender_table_LGA_SGA)$exGendercted
#No
#Do the fisher's exact test 
fisher_Gender_LGA_SGA <- fisher.test(Gender_table_LGA_SGA)
fisher_Gender_LGA_SGA


