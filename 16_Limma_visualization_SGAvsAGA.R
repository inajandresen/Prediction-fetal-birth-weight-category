library(tidyverse)
library(writexl)
library(readxl)
library(Biobase)
library(limma)
#library(ggVennDiagram)
library(pcaMethods)
library(scales)
library(gridExtra)
library(ggforce)
library(ggrepel)
library(ggvenn)
library(ggpubr)
library(pheatmap)
library(ggplotify)


# P-value histograms -----------------------------------------
T1_SGAvsAGA <- read_excel("../T1_SGA_AGA_adj_nulli_bmi.xlsx")
T2_SGAvsAGA <- read_excel("../T2_SGA_AGA_adj_nulli_bmi.xlsx")
T3_SGAvsAGA <- read_excel("../T3_SGA_AGA_adj_nulli_bmi.xlsx")

T1_SGAvsAGA_adj_nulli_bmi_hist <- ggplot(T1_SGAvsAGA, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  #geom_histogram(binwidth = 0.5, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Visit 1") +
  theme_classic()

T2_SGAvsAGA_adj_nulli_bmi_hist <- ggplot(T2_SGAvsAGA, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  #geom_histogram(binwidth = 0.5, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Visit 2") +
  theme_classic()

T3_SGAvsAGA_adj_nulli_bmi_hist <- ggplot(T3_SGAvsAGA, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  #geom_histogram(binwidth = 0.5, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Visit 3") +
  theme_classic()

histogram_adj_nulli_bmi_SGAvsAGA <- grid.arrange(T1_SGAvsAGA_adj_nulli_bmi_hist, T2_SGAvsAGA_adj_nulli_bmi_hist, T3_SGAvsAGA_adj_nulli_bmi_hist, nrow = 1, top = "SGA vs. AGA")

ggsave(filename= "Histogram_p-values_SGA_AGA_adj_nulli_bmi.png",
       plot = histogram_adj_nulli_bmi_SGAvsAGA,
       device = "png",
       path = "../",
       width = 9,
       height = 3,
       units = "in")




# Venn diagram SGA vs AGA  ------------------------------------------------
T1_SGAvsAGA <- read_excel("../T1_SGA_AGA_adj_nulli_bmi.xlsx")
T2_SGAvsAGA <- read_excel("../T2_SGA_AGA_adj_nulli_bmi.xlsx")
T3_SGAvsAGA <- read_excel("../T3_SGA_AGA_adj_nulli_bmi.xlsx")

#Venn diagram of significant proteins 
T1_SGAvsAGA_sig <- T1_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)

T2_SGAvsAGA_sig <- T2_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)

T3_SGAvsAGA_sig <- T3_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)

#Plot venn diagram, LGA vs AGA
T1T2T3_SGAvsAGA <- list("Visit 1" = T1_SGAvsAGA_sig$AptName, "Visit 2" = T2_SGAvsAGA_sig$AptName, "Visit 3" = T3_SGAvsAGA_sig$AptName)

vennplot_T1T2T3_SGAvsAGA <- ggvenn(T1T2T3_SGAvsAGA, 
                                   fill_color = c("lightblue", "lightblue", "lightblue"),
                                   stroke_size = 0.5, 
                                   set_name_size = 5,
                                   show_percentage = FALSE,
                                   text_size = 5)



ggsave(filename= "Venn_SGA_AGA_adj_nulli_bmi.png",
       plot = vennplot_T1T2T3_SGAvsAGA,
       device = "png",
       path = "../",
       width = 10,
       height = 10,
       units = "cm",
       bg = "white")



# Volcano plots ------------------------------------------------
load("../efit_LGA_AGA_SGA_adj_nulli_bmi.RData")

#LGA vs AGA
T1_SGAvsAGA_adj_nulli_bmi <- topTable(efit, coef = "T1_SGAvsAGA", number = "all", adjust.method="BH")
sig_T1_SGAvsAGA_adj_nulli_bmi <- T1_SGAvsAGA_adj_nulli_bmi %>%
  filter(adj.P.Val < 0.05)
dim(sig_T1_SGAvsAGA_adj_nulli_bmi)
T2_SGAvsAGA_adj_nulli_bmi <- topTable(efit, coef = "T2_SGAvsAGA", number = "all", adjust.method="BH")
sig_T2_SGAvsAGA_adj_nulli_bmi <- T2_SGAvsAGA_adj_nulli_bmi %>%
  filter(adj.P.Val < 0.05)
dim(sig_T2_SGAvsAGA_adj_nulli_bmi)
T3_SGAvsAGA_adj_nulli_bmi <- topTable(efit, coef = "T3_SGAvsAGA", number = "all", adjust.method="BH")
sig_T3_SGAvsAGA_adj_nulli_bmi <- T3_SGAvsAGA_adj_nulli_bmi %>%
  filter(adj.P.Val < 0.05)
dim(sig_T3_SGAvsAGA_adj_nulli_bmi)

#Time interval 1
#Add a column displaying the differential expressed proteins to the data frame
T1_SGAvsAGA_adj_nulli_bmi <- T1_SGAvsAGA_adj_nulli_bmi %>%
  remove_rownames()
T1_SGAvsAGA_adj_nulli_bmi$DE_0.05 <- "NO"
T1_SGAvsAGA_adj_nulli_bmi$DE_0.05[T1_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_SGAvsAGA_adj_nulli_bmi$logFC > 0] <- "UP"
T1_SGAvsAGA_adj_nulli_bmi$DE_0.05[T1_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_SGAvsAGA_adj_nulli_bmi$logFC < 0] <- "DOWN"
T1_SGAvsAGA_adj_nulli_bmi$label <- NA
T1_SGAvsAGA_adj_nulli_bmi$label[T1_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_SGAvsAGA_adj_nulli_bmi$logFC > 0] <- T1_SGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T1_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_SGAvsAGA_adj_nulli_bmi$logFC > 0]
T1_SGAvsAGA_adj_nulli_bmi$label[T1_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_SGAvsAGA_adj_nulli_bmi$logFC < 0] <- T1_SGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T1_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_SGAvsAGA_adj_nulli_bmi$logFC < 0]

#Plot volcanoplot
Volcano_T1_SGAvsAGA_adj_nulli_bmi <- ggplot(data = T1_SGAvsAGA_adj_nulli_bmi, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 4.5,
                  color = "black") +
  theme_classic() +
  #labs(tag = "A)") +
  xlab("logFC") +
  ylab("-log10 q-value") +
  scale_y_continuous(limits = c(0, 2.5), 
                     breaks = seq(0, 2.5, by = 0.5)) +
  scale_x_continuous(limits = c(-3.5,2.1), 
                     breaks = seq(-4, 2, by = 1)) +
  scale_color_manual(values = c("grey", "brown2")) +
  theme(legend.position = "none") +
  ggtitle("Visit 1", subtitle = "Week 12-19") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))


#Time interval 2
#Add a column displaying the differential expressed proteins to the data frame
T2_SGAvsAGA_adj_nulli_bmi <- T2_SGAvsAGA_adj_nulli_bmi %>%
  remove_rownames()
T2_SGAvsAGA_adj_nulli_bmi$DE_0.05 <- "NO"
T2_SGAvsAGA_adj_nulli_bmi$DE_0.05[T2_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_SGAvsAGA_adj_nulli_bmi$logFC > 0] <- "UP"
T2_SGAvsAGA_adj_nulli_bmi$DE_0.05[T2_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_SGAvsAGA_adj_nulli_bmi$logFC < 0] <- "DOWN"
T2_SGAvsAGA_adj_nulli_bmi$label <- NA
T2_SGAvsAGA_adj_nulli_bmi$label[T2_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_SGAvsAGA_adj_nulli_bmi$logFC > 0] <- T2_SGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T2_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_SGAvsAGA_adj_nulli_bmi$logFC > 0]
T2_SGAvsAGA_adj_nulli_bmi$label[T2_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_SGAvsAGA_adj_nulli_bmi$logFC < 0] <- T2_SGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T2_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_SGAvsAGA_adj_nulli_bmi$logFC < 0]

#Plot volcanoplot
Volcano_T2_SGAvsAGA_adj_nulli_bmi <- ggplot(data = T2_SGAvsAGA_adj_nulli_bmi, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 15,
                  segment.color = 'grey50',
                  size = 4.5,
                  color = "black") +
  theme_classic() +
  #labs(tag = "B)") +
  xlab("logFC") +
  ylab("-log10 q-value") +
  scale_y_continuous(limits = c(0, 2.5), 
                     breaks = seq(0, 2.5, by = 0.5)) +
  scale_x_continuous(limits = c(-3.5,2.1), 
                     breaks = seq(-4, 2, by = 1)) +
  scale_color_manual(values = c("grey", "brown2")) +
  theme(legend.position = "none") +
  ggtitle("Visit 2", subtitle = "Week 21-27") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))




#Time interval 3
#Add a column displaying the differential expressed proteins to the data frame
T3_SGAvsAGA_adj_nulli_bmi <- T3_SGAvsAGA_adj_nulli_bmi %>%
  remove_rownames()
T3_SGAvsAGA_adj_nulli_bmi$DE_0.05 <- "NO"
T3_SGAvsAGA_adj_nulli_bmi$DE_0.05[T3_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_SGAvsAGA_adj_nulli_bmi$logFC > 0] <- "UP"
T3_SGAvsAGA_adj_nulli_bmi$DE_0.05[T3_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_SGAvsAGA_adj_nulli_bmi$logFC < 0] <- "DOWN"
T3_SGAvsAGA_adj_nulli_bmi$label <- NA
T3_SGAvsAGA_adj_nulli_bmi$label[T3_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_SGAvsAGA_adj_nulli_bmi$logFC > 0] <- T3_SGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T3_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_SGAvsAGA_adj_nulli_bmi$logFC > 0]
T3_SGAvsAGA_adj_nulli_bmi$label[T3_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_SGAvsAGA_adj_nulli_bmi$logFC < 0] <- T3_SGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T3_SGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_SGAvsAGA_adj_nulli_bmi$logFC < 0]

#Plot volcanoplot
Volcano_T3_SGAvsAGA_adj_nulli_bmi <- ggplot(data = T3_SGAvsAGA_adj_nulli_bmi, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  #geom_text_repel(aes(label = label),
  #                #label.padding = 1, 
  #                point.padding = 0.5,
  #                max.overlaps = 15,
  #               segment.color = 'grey50',
  #                size = 4.5,
  #                color = "black") +
  theme_classic() +
  #labs(tag = "C)") +
  xlab("logFC") +
  ylab("-log10 q-value") +
  scale_y_continuous(limits = c(0, 2.5), 
                     breaks = seq(0, 2.5, by = 0.5)) +
  scale_x_continuous(limits = c(-3.5,2.1), 
                     breaks = seq(-4, 2, by = 1)) +
  scale_color_manual(values = c("cornflowerblue", "grey", "brown2")) +
  theme(legend.position = "none") +
  ggtitle("Visit 3", subtitle = "Week 28-34") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#Combine volcanoplots only
SGAvsAGA_Volcanoplots_adj_nulli_bmi <- grid.arrange(Volcano_T1_SGAvsAGA_adj_nulli_bmi, Volcano_T2_SGAvsAGA_adj_nulli_bmi, Volcano_T3_SGAvsAGA_adj_nulli_bmi, 
                                                    nrow = 1, widths = c(1.0, 1.0, 1.0), 
                                                    top = textGrob("SGA vs. AGA", gp = gpar(fontsize = 15)))

ggsave(filename= "Volcano_SGA_AGA_adj_nulli_bmi.png",
       plot = SGAvsAGA_Volcanoplots_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 35,
       height = 13,
       units = "cm")



# Heat map, SGA vs AGA ----------------------------------------------------
#Get a list of the significant proteins, LGA vs AGA
T1_SGAvsAGA <- read_excel("../T1_SGA_AGA_adj_nulli_bmi.xlsx")
T1_SGAvsAGA_sig <- T1_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
T1_SGAvsAGA_sig <- T1_SGAvsAGA_sig$AptName

T2_SGAvsAGA <- read_excel("../T2_SGA_AGA_adj_nulli_bmi.xlsx")
T2_SGAvsAGA_sig <- T2_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
T2_SGAvsAGA_sig <- T2_SGAvsAGA_sig$AptName

T3_SGAvsAGA <- read_excel("../T3_SGA_AGA_adj_nulli_bmi.xlsx")
T3_SGAvsAGA_sig <- T3_SGAvsAGA %>%
  filter(adj.P.Val < 0.05)
T3_SGAvsAGA_sig <- T3_SGAvsAGA_sig$AptName


SGAvsAGA_SigProteins <- c(T1_SGAvsAGA_sig, T2_SGAvsAGA_sig, T3_SGAvsAGA_sig)
SGAvsAGA_SigProteins <- unique(SGAvsAGA_SigProteins)

#Format a data frame with logFC
T1_SGAvsAGA <- T1_SGAvsAGA %>%
  dplyr::select(c(logFC, AptName, EntrezGeneSymbol, TargetFullName)) %>%
  dplyr::rename("Visit 1" = logFC)

T2_SGAvsAGA <- T2_SGAvsAGA %>%
  dplyr::select(c(logFC, AptName, EntrezGeneSymbol, TargetFullName)) %>%
  dplyr::rename("Visit 2" = logFC)

T3_SGAvsAGA <- T3_SGAvsAGA %>%
  dplyr::select(c(logFC, AptName, EntrezGeneSymbol, TargetFullName)) %>%
  dplyr::rename("Visit 3" = logFC)

LogFC_SGAvsAGA <- merge(T1_SGAvsAGA, T2_SGAvsAGA, by = "AptName", all = "TRUE") %>%
  merge(T3_SGAvsAGA, by = "AptName", all = "TRUE") 

#Filter on the significant proteins 
LogFC_SGAvsAGA <- filter(LogFC_SGAvsAGA, grepl(paste(SGAvsAGA_SigProteins, collapse = "|"), AptName)) #Filter on common SomaIds
LogFC_SGAvsAGA <- LogFC_SGAvsAGA %>%
  dplyr::mutate(EntrezGeneSymbol2 = paste0("(", EntrezGeneSymbol,  ")"))
LogFC_SGAvsAGA$TargetFullName_EntrezGeneSymbol <- paste(LogFC_SGAvsAGA$TargetFullName, LogFC_SGAvsAGA$EntrezGeneSymbol2)

#Check duplicates, Targetfullname
SGAvsAGA_TargetFullName_EntrezGeneSymbol <- LogFC_SGAvsAGA$TargetFullName_EntrezGeneSymbol #Extract to check for duplicate names
SGAvsAGA_dup <- SGAvsAGA_TargetFullName_EntrezGeneSymbol[duplicated(SGAvsAGA_TargetFullName_EntrezGeneSymbol)] #Extract the duplicates

LogFC_SGAvsAGA <- LogFC_SGAvsAGA %>% 
  remove_rownames() %>% 
  column_to_rownames(var="TargetFullName_EntrezGeneSymbol") %>%
  #column_to_rownames(var="EntrezGeneSymbol2") %>%
  select(`Visit 1`, `Visit 2`, `Visit 3`)

#Set value zero to be white color
Breaks <- c(seq(min(LogFC_SGAvsAGA), 0, length.out=ceiling(100/2) + 1), 
            seq(max(LogFC_SGAvsAGA)/100, max(LogFC_SGAvsAGA), length.out=floor(100/2)))


heatmap_SGAvsAGA <- LogFC_SGAvsAGA %>%
  pheatmap(border_color = NA,
           breaks = Breaks,
           clustering_method = "complete",
           cluster_cols = FALSE,
           show_rownames = TRUE,
           color = colorRampPalette(c("cornflowerblue", "white", "brown2"))(100),
           fontsize_row = 12,
           fontsize_col = 14,
           treeheight_row = 25,
           cellwidth = 20,
           cellheight = 18,
           main = "")

ggsave(filename= "Heatmap_SGAvsAGA_DE_proteins.png",
       plot = heatmap_SGAvsAGA,
       device = "png",
       path = "../",
       width = 18,
       height = 6,
       units = "cm")



# Combine plots SGA vs AGA ------------------------------------------------


#Combine volcano and venn
SGAvsAGA_heatmap_gr <- as.grob(heatmap_SGAvsAGA)

SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- ggarrange(Volcano_T1_SGAvsAGA_adj_nulli_bmi, Volcano_T2_SGAvsAGA_adj_nulli_bmi, 
                                                              Volcano_T3_SGAvsAGA_adj_nulli_bmi, vennplot_T1T2T3_SGAvsAGA, SGAvsAGA_heatmap_gr,
                                                              ncol = 2, nrow = 2, widths = c(1, 1),
                                                              labels = c("A)", "B)", "C)", "D)"))

first_row <- ggarrange(Volcano_T1_SGAvsAGA_adj_nulli_bmi, Volcano_T2_SGAvsAGA_adj_nulli_bmi, Volcano_T3_SGAvsAGA_adj_nulli_bmi,
                       ncol = 3, widths = c(1, 1, 1),
                       labels = c("A)", "B)", "C)"))

second_row <- ggarrange(vennplot_T1T2T3_SGAvsAGA, SGAvsAGA_heatmap_gr,
                        ncol = 2, widths = c(0.8, 1),
                        labels = c("D)", "E)"))

SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- ggarrange(first_row, second_row, ncol = 1, nrow = 2,
                                                              heights = c(1, 0.8))


ggsave(filename= "Volcano_venn_heatmap_SGA_AGA_adj_nulli_bmi.png",
       plot = SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 35,
       height = 20,
       units = "cm",
       bg = "white")


#Combine volcano, venn and heatmap
SGAvsAGA_heatmap_gr <- as.grob(heatmap_SGAvsAGA)

Column_1 <- ggarrange(Volcano_T1_SGAvsAGA_adj_nulli_bmi, Volcano_T2_SGAvsAGA_adj_nulli_bmi, 
                      Volcano_T3_SGAvsAGA_adj_nulli_bmi, vennplot_T1T2T3_SGAvsAGA, 
                      ncol = 1, labels = c("A)", "B)", "C)", "D)"))

Column_2 <- ggarrange(SGAvsAGA_heatmap_gr,
                      ncol = 1, labels = c("E)"))

SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- ggarrange(Column_1, Column_2, ncol = 2, 
                                                              nrow = 1, widths = c(0.7, 1))


#SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- annotate_figure(SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi, 
#                                                                    top = textGrob("LGA vs. AGA", gp = gpar(fontsize = 25)))

ggsave(filename= "Volcano_venn_heatmap_LGA_AGA_adj_nulli_bmi.png",
       plot = SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 40,
       height = 50,
       units = "cm",
       bg = "white")


# Combine plots SGA vs AGA ------------------------------------------------


#Combine volcano and venn
SGAvsAGA_heatmap_gr <- as.grob(heatmap_SGAvsAGA)

SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- ggarrange(Volcano_T1_SGAvsAGA_adj_nulli_bmi, Volcano_T2_SGAvsAGA_adj_nulli_bmi, 
                                                              Volcano_T3_SGAvsAGA_adj_nulli_bmi, SGAvsAGA_heatmap_gr,
                                                              ncol = 2, nrow = 2, widths = c(1, 1))



ggsave(filename= "Volcano_venn_heatmap_SGA_AGA_adj_nulli_bmi_2.png",
       plot = SGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 25,
       height = 20,
       units = "cm",
       bg = "white")


