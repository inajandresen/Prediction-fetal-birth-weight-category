library(tidyverse)
library(Biobase)
library(pcaMethods)
library(scales)
library(gridExtra)
library(ggpubr)

load("../Data/MOM_exprset_STORK_SGA_LGA.RData")
T1_MOM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$TimePoint %in% "1"]
T1_MOM <- exprs(T1_MOM_exprset_STORK_LGA_SGA)

pc1 <- pca(T1_MOM_exprset_STORK_LGA_SGA)
df1 <- merge(scores(pc1), pData(T1_MOM_exprset_STORK_LGA_SGA), by = 0)
t1_pc1_r2 <- 100*pc1@R2[1]
t1_pc1_r2 <- t1_pc1_r2
pca1 <- ggplot(df1, aes(PC1, PC2, color = Bwkat)) +
  geom_point(size = 2, alpha = 0.8) +
  #stat_ellipse(type = "t") +
  labs(title = "Visit 1", tag = "A)", subtitle = "Week 12-19") + 
  xlab(paste0("PC1 (R2 = ",round(100*pc1@R2[1], digits = 2),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc1@R2[2], digits = 2),"%)")) +
  xlim(c(-55, 80)) +
  ylim(c(-30, 40)) +
  #scale_y_continuous(labels = comma) +
  #scale_x_continuous(labels = comma) +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.position = "none",
        #legend.title = element_blank(),
        #legend.text = element_text(size = 12),
        title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


T2_MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$TimePoint %in% "2"]

pc2 <- pca(T2_MoM_exprset_STORK_LGA_SGA)
df2 <- merge(scores(pc2), pData(T2_MoM_exprset_STORK_LGA_SGA), by = 0)
pca2 <- ggplot(df2, aes(PC1, PC2, color = Bwkat)) +
  geom_point(size = 2, alpha = 0.8) +
  #stat_ellipse(type = "t") +
  labs(title = "Visit 2", tag = "B)", subtitle = "Week 21-27") + 
  xlab(paste0("PC1 (R2 = ",round(100*pc2@R2[1], digits = 2),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc2@R2[2], digits = 2),"%)")) +
  xlim(c(-55, 80)) +
  ylim(c(-30, 40)) +
  #scale_y_continuous(labels = comma) +
  #scale_x_continuous(labels = comma) +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.position = "none",
        #legend.title = element_blank(),
        #legend.text = element_text(size = 12),
        title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

T3_MoM_exprset_STORK_LGA_SGA <- MoM_exprset_STORK_LGA_SGA[, MoM_exprset_STORK_LGA_SGA$TimePoint %in% "3"]

pc3 <- pca(T3_MoM_exprset_STORK_LGA_SGA)
df3 <- merge(scores(pc3), pData(T3_MoM_exprset_STORK_LGA_SGA), by = 0)
pca3 <- ggplot(df3, aes(PC1, PC2, color = Bwkat)) +
  geom_point(size = 2, alpha = 0.8) +
  #stat_ellipse(type = "t") +
  labs(title = "Visit 3", tag = "C)", subtitle = "Week 28-34") + 
  xlab(paste0("PC1 (R2 = ",round(100*pc3@R2[1], digits = 2),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc3@R2[2], digits = 2),"%)")) +
  xlim(c(-55, 80)) +
  ylim(c(-30, 40)) +
  #scale_y_continuous(labels = comma) +
  #scale_x_continuous(labels = comma) +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.position = "none",
        #legend.title = element_blank(),
        #legend.text = element_text(size = 12),
        title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


#Extract the legend
for_lenged1 <- ggplot(df3, aes(PC1, PC2, color = Bwkat)) +
  geom_point(size = 2, alpha = 0.8) +
  #stat_ellipse(type = "t") +
  scale_color_manual(values = c("Grey", "dodgerblue", "orangered2")) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        #legend.text = element_text(size = 12),
        title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size=12))

get_legend1 <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend1 <- get_legend(for_lenged1)

pca <- grid.arrange(pca1, pca2, pca3, legend1,
                    widths = c(1.0, 1.0, 1.0, 0.2))


ggsave(filename= "PCA_samples_same_axis.png",
       plot = pca,
       device = "png",
       path = "../Sample_clustring/",
       width = 35,
       height = 12,
       units = "cm")

