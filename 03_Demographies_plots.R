library(tidyverse)
library(Biobase)
library(writexl)
library(gridExtra)
library(psych)
library(pcaMethods)
library(scales)
library(grid)
library(ggpubr)

#Load the data and format so you have the birth weight z scores as one column, and the proteins abundances as separate columns
load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$PE %in% "0"]

samp <- pData(MoM_exprset_STORK_LGA_SGA)
samp_unique <- samp %>%  #Extract one row per patient
  distinct(ID, .keep_all = TRUE)
table(samp_unique$Bwkat)
#AGA: 58
#LGA: 5
#SGA: 7



# Placenta weight against birth weight (absolute values) EXCLUDED FROM PLOT ---------------------------------------------------------

#Plot Birth weight against placenta weight
bweight_plweight_plot <- ggplot(samp_unique, aes(PlWeigth, `Birth weight`, color = Bwkat)) +
  geom_point() +
  geom_smooth(method = "lm",
              mapping = aes(PlWeigth, `Birth weight`), 
              inherit.aes = FALSE, color = "black") +
  #annotate("text", x = 1050, y = 4000, label = "R = 0.55\np-value < 0.001", size = 3.5) +
  xlab("Placenta weigth (g)") +
  ylab("Birthweight (g)") +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_cor(method = "spearman", mapping = aes(PlWeigth, `Birth weight`), 
        inherit.aes = FALSE, label.y = 5000)


# Placenta weight against birth weight (z-score) EXCLUDED FROM PLOT -----------------------------------------------

zbweight_zplweight_plot <- ggplot(samp_unique, aes(`PWZ`, `Birthweight z-score`, color = Bwkat)) +
  geom_point() +
  geom_smooth(method = "lm",
              mapping = aes(`PWZ`, `Birthweight z-score`), 
              inherit.aes = FALSE, color = "black") +
  #annotate("text", x = 2.5, y = 1.9, label = "R = 0.58\np-value < 0.001", size = 3.5) +
  labs(tag = "A)") +
  xlab("Placenta weight z-score") +
  ylab("Birthweight z-score") +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
   stat_cor(method = "spearman", mapping = aes(`PWZ`, `Birthweight z-score`), 
        inherit.aes = FALSE) 


# BMI ---------------------------------------------------------------------

#Plot BMI against birth weight z-score
bmi_plot <- ggplot(samp_unique, aes(BMI, `Birthweight z-score`, color = Bwkat)) +
  geom_point() +
  geom_smooth(method = "lm",
              mapping = aes(BMI, `Birthweight z-score`), 
              inherit.aes = FALSE, color = "black") +
  #annotate("text", x = 36.5, y = 1.4, label = "R = 0.32\np-value < 0.01", size = 3.5) +
  labs(tag = "A)") +
  xlab("BMI (kg/m2)") +
  ylab("Birthweight z-score") +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  #scale_x_continuous(breaks = seq(15,35,5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  stat_cor(method = "spearman", mapping = aes(BMI, `Birthweight z-score`), 
           inherit.aes = FALSE) 
  #stat_cor(method = "pearson", mapping = aes(BMI, `Birthweight z-score`), 
  #      inherit.aes = FALSE)


# Age ---------------------------------------------------------------------

#Plot age against birth weight z-score
age_plot <- ggplot(samp_unique, aes(Age, `Birthweight z-score`, color = Bwkat)) +
  geom_point() +
  geom_smooth(method = "lm",
              mapping = aes(Age, `Birthweight z-score`), 
              inherit.aes = FALSE, color = "black") +
  #annotate("text", x = 39, y = 0.5, label = "R = 0.081\np-value > 0.4", size = 3.5) +
  labs(tag = "B)") +
  xlab("Age (years)") +
  ylab("Birthweight z-score") +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  stat_cor(method = "spearman", mapping = aes(Age, `Birthweight z-score`), 
           inherit.aes = FALSE)
  #stat_cor(method = "pearson", mapping = aes(Age, `Birthweight z-score`), 
  #         inherit.aes = FALSE)




# Nulliparity -------------------------------------------------------------

samp_unique$Nulliparity <- as.factor(samp_unique$Nulliparity)
samp_unique <- samp_unique %>%
  mutate("Nulliparity2" = Nulliparity)
samp_unique$Nulliparity2 <- sub("1", "Nulliparous", samp_unique$Nulliparity2)
samp_unique$Nulliparity2 <- sub("0", "Multiparous", samp_unique$Nulliparity2)

nulli_ttest <- list(c("Nulliparous", "Multiparous"))

#Plot nulliparity and birth weigth z-score
nulli_plot <- ggplot(samp_unique, aes(Nulliparity2, `Birthweight z-score`)) +
  geom_boxplot() +
  labs(tag = "C)") +
  xlab("") +
  ylab("Birthweight z-score") +
  #ggtitle("t statistics = 3.0, p-value < 0.01") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  stat_compare_means(method = "t.test", comparisons = nulli_ttest) 





# Gender ------------------------------------------------------------------

samp_unique$Gender <- as.factor(samp_unique$Gender)
samp_unique <- samp_unique %>%
  mutate("Gender2" = Gender)
samp_unique$Gender2 <- sub("1", "Female", samp_unique$Gender2)
samp_unique$Gender2 <- sub("0", "Male", samp_unique$Gender2)

gender_ttest <- list(c("Female", "Male"))


#Plot gender and birth weigth z-score
Gender_plot <- ggplot(samp_unique, aes(Gender2, `Birthweight z-score`)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  labs(tag = "D)") +
  xlab("") +
  ylab("Birthweigth z-score") +
  #ggtitle("t statistics = 0.82, p-value > 0.4") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  stat_compare_means(method = "t.test", comparisons = gender_ttest) 



# Combine plots -----------------------------------------------------------

#Extract legend 
#Extract the legend
for_lenged <- ggplot(samp_unique, aes(Age, `Birthweight z-score`, color = Bwkat)) +
  geom_point() +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10))

get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(for_lenged)


text <- textGrob("")

plots <- grid.arrange(bmi_plot, age_plot, legend,
                      nulli_plot, Gender_plot, text,
                      ncol = 3, widths = c(1.0, 1.0, 0.2))


ggsave(filename= "Clinical_variables_plots_without_PE_spearman.png",
       plot = plots,
       device = "png",
       path = "../Demographics",
       #width = 40,
       #height = 25,
       width = 20,
       height = 19,
       units = "cm")




# Interaction between nulliparity and BMI ---------------------------------

#Plot nulliparity and birth weigth z-score
nulli_plot <- ggplot(samp_unique, aes(Nulliparity2, BMI)) +
  geom_boxplot() +
  #labs(tag = "D)") +
  xlab("") +
  ylab("BMI") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  stat_compare_means(method = "t.test", comparisons = nulli_ttest) 

#Although multiparous woman have a higher mean BMI than nulliparous women, the means are not statistically 
#significantly different


# BMI, seperate groups  ---------------------------------------------------------------------

#Plot BMI against birth weight z-score
bmi_sep_plot <- ggplot(samp_unique, aes(BMI, `Birthweight z-score`, color = Bwkat, fill = Bwkat)) +
  geom_point() +
  labs(tag = "B)") +
  xlab("BMI (kg/m2)") +
  ylab("Birthweight z-score") +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  scale_fill_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  geom_smooth(aes(linetype = Bwkat),
              method = "lm", formula = y ~ x, se = FALSE) +
  #scale_x_continuous(breaks = seq(15,35,5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  stat_cor(method = "spearman", r.digits = 2, p.digits = 2, group = 3) 
