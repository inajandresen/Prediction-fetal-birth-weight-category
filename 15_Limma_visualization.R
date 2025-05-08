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


# P-value histograms, LGA vs AGA  -----------------------------------------
T1_LGAvsAGA <- read_excel("../T1_LGA_AGA_adj_nulli_bmi.xlsx")
T2_LGAvsAGA <- read_excel("../T2_LGA_AGA_adj_nulli_bmi.xlsx")
T3_LGAvsAGA <- read_excel("../T3_LGA_AGA_adj_nulli_bmi.xlsx")

T1_LGAvsAGA_adj_nulli_bmi_hist <- ggplot(T1_LGAvsAGA, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  #geom_histogram(binwidth = 0.5, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Visit 1") +
  theme_classic()

T2_LGAvsAGA_adj_nulli_bmi_hist <- ggplot(T2_LGAvsAGA, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  #geom_histogram(binwidth = 0.5, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Visit 2") +
  theme_classic()

T3_LGAvsAGA_adj_nulli_bmi_hist <- ggplot(T3_LGAvsAGA, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  #geom_histogram(binwidth = 0.5, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Visit 3") +
  theme_classic()

histogram_adj_nulli_bmi_LGAvsAGA <- grid.arrange(T1_LGAvsAGA_adj_nulli_bmi_hist, T2_LGAvsAGA_adj_nulli_bmi_hist, T3_LGAvsAGA_adj_nulli_bmi_hist, nrow = 1, top = "LGA vs. AGA")

ggsave(filename= "Histogram_p-values_LGA_AGA_adj_nulli_bmi.png",
       plot = histogram_adj_nulli_bmi_LGAvsAGA,
       device = "png",
       path = "../",
       width = 9,
       height = 3,
       units = "in")


# Venn diagram LGA vs AGA  ------------------------------------------------
T1_LGAvsAGA <- read_excel("../T1_LGA_AGA_adj_nulli_bmi.xlsx")
T2_LGAvsAGA <- read_excel("../T2_LGA_AGA_adj_nulli_bmi.xlsx")
T3_LGAvsAGA <- read_excel("../T3_LGA_AGA_adj_nulli_bmi.xlsx")

#Venn diagram of significant proteins 
T1_LGAvsAGA_sig <- T1_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)

T2_LGAvsAGA_sig <- T2_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)

T3_LGAvsAGA_sig <- T3_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)

#Plot venn diagram, LGA vs AGA
T1T2T3_LGAvsAGA <- list("Visit 1" = T1_LGAvsAGA_sig$AptName, "Visit 2" = T2_LGAvsAGA_sig$AptName, "Visit 3" = T3_LGAvsAGA_sig$AptName)

vennplot_T1T2T3_LGAvsAGA <- ggvenn(T1T2T3_LGAvsAGA, 
                                   fill_color = c("lightblue", "lightblue", "lightblue"),
                                   stroke_size = 0.5, 
                                   set_name_size = 7,
                                   show_percentage = FALSE,
                                   text_size = 6.5)



ggsave(filename= "Venn_LGA_AGA_adj_nulli_bmi.png",
       plot = vennplot_T1T2T3_LGAvsAGA,
       device = "png",
       path = "../",
       width = 10,
       height = 10,
       units = "cm",
       bg = "white")

Common_60 <- intersect(T1_LGAvsAGA_sig$AptName, T2_LGAvsAGA_sig$AptName)
Common_60 <- intersect(Common_60, T3_LGAvsAGA_sig$AptName)


# Volcano plots LGA vs AGA ------------------------------------------------
load("../efit_LGA_AGA_SGA_adj_nulli_bmi.RData")

#LGA vs AGA
T1_LGAvsAGA_adj_nulli_bmi <- topTable(efit, coef = "T1_LGAvsAGA", number = "all", adjust.method="BH")
sig_T1_LGAvsAGA_adj_nulli_bmi <- T1_LGAvsAGA_adj_nulli_bmi %>%
  filter(adj.P.Val < 0.05)
dim(sig_T1_LGAvsAGA_adj_nulli_bmi)
T2_LGAvsAGA_adj_nulli_bmi <- topTable(efit, coef = "T2_LGAvsAGA", number = "all", adjust.method="BH")
sig_T2_LGAvsAGA_adj_nulli_bmi <- T2_LGAvsAGA_adj_nulli_bmi %>%
  filter(adj.P.Val < 0.05)
dim(sig_T2_LGAvsAGA_adj_nulli_bmi)
T3_LGAvsAGA_adj_nulli_bmi <- topTable(efit, coef = "T3_LGAvsAGA", number = "all", adjust.method="BH")
sig_T3_LGAvsAGA_adj_nulli_bmi <- T3_LGAvsAGA_adj_nulli_bmi %>%
  filter(adj.P.Val < 0.05)
dim(sig_T3_LGAvsAGA_adj_nulli_bmi)

#Time interval 1
#Add a column displaying the differential expressed proteins to the data frame
T1_LGAvsAGA_adj_nulli_bmi <- T1_LGAvsAGA_adj_nulli_bmi %>%
  remove_rownames()
T1_LGAvsAGA_adj_nulli_bmi$DE_0.05 <- "NO"
T1_LGAvsAGA_adj_nulli_bmi$DE_0.05[T1_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_LGAvsAGA_adj_nulli_bmi$logFC > 0] <- "UP"
T1_LGAvsAGA_adj_nulli_bmi$DE_0.05[T1_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_LGAvsAGA_adj_nulli_bmi$logFC < 0] <- "DOWN"
T1_LGAvsAGA_adj_nulli_bmi$label <- NA
T1_LGAvsAGA_adj_nulli_bmi$label[T1_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_LGAvsAGA_adj_nulli_bmi$logFC > 0] <- T1_LGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T1_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_LGAvsAGA_adj_nulli_bmi$logFC > 0]
T1_LGAvsAGA_adj_nulli_bmi$label[T1_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_LGAvsAGA_adj_nulli_bmi$logFC < 0] <- T1_LGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T1_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T1_LGAvsAGA_adj_nulli_bmi$logFC < 0]

#Plot volcanoplot
Volcano_T1_LGAvsAGA_adj_nulli_bmi <- ggplot(data = T1_LGAvsAGA_adj_nulli_bmi, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 5,
                  color = "black") +
  theme_classic() +
  #labs(tag = "A)") +
  xlab("logFC") +
  ylab("-log10 q-value") +
  scale_y_continuous(limits = c(0, 11), 
                     breaks = seq(0, 12, by = 2)) +
  scale_x_continuous(limits = c(-3,2.1), 
                     breaks = seq(-3, 2, by = 1)) +
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
T2_LGAvsAGA_adj_nulli_bmi <- T2_LGAvsAGA_adj_nulli_bmi %>%
  remove_rownames()
T2_LGAvsAGA_adj_nulli_bmi$DE_0.05 <- "NO"
T2_LGAvsAGA_adj_nulli_bmi$DE_0.05[T2_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_LGAvsAGA_adj_nulli_bmi$logFC > 0] <- "UP"
T2_LGAvsAGA_adj_nulli_bmi$DE_0.05[T2_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_LGAvsAGA_adj_nulli_bmi$logFC < 0] <- "DOWN"
T2_LGAvsAGA_adj_nulli_bmi$label <- NA
T2_LGAvsAGA_adj_nulli_bmi$label[T2_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_LGAvsAGA_adj_nulli_bmi$logFC > 0] <- T2_LGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T2_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_LGAvsAGA_adj_nulli_bmi$logFC > 0]
T2_LGAvsAGA_adj_nulli_bmi$label[T2_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_LGAvsAGA_adj_nulli_bmi$logFC < 0] <- T2_LGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T2_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T2_LGAvsAGA_adj_nulli_bmi$logFC < 0]

#Plot volcanoplot
Volcano_T2_LGAvsAGA_adj_nulli_bmi <- ggplot(data = T2_LGAvsAGA_adj_nulli_bmi, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 15,
                  segment.color = 'grey50',
                  size = 5,
                  color = "black") +
  theme_classic() +
  #labs(tag = "B)") +
  xlab("logFC") +
  ylab("-log10 q-value") +
  scale_y_continuous(limits = c(0, 11), 
                     breaks = seq(0, 12, by = 2)) +
  scale_x_continuous(limits = c(-3,2.1), 
                     breaks = seq(-3, 2, by = 1)) +
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
T3_LGAvsAGA_adj_nulli_bmi <- T3_LGAvsAGA_adj_nulli_bmi %>%
  remove_rownames()
T3_LGAvsAGA_adj_nulli_bmi$DE_0.05 <- "NO"
T3_LGAvsAGA_adj_nulli_bmi$DE_0.05[T3_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_LGAvsAGA_adj_nulli_bmi$logFC > 0] <- "UP"
T3_LGAvsAGA_adj_nulli_bmi$DE_0.05[T3_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_LGAvsAGA_adj_nulli_bmi$logFC < 0] <- "DOWN"
T3_LGAvsAGA_adj_nulli_bmi$label <- NA
T3_LGAvsAGA_adj_nulli_bmi$label[T3_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_LGAvsAGA_adj_nulli_bmi$logFC > 0] <- T3_LGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T3_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_LGAvsAGA_adj_nulli_bmi$logFC > 0]
T3_LGAvsAGA_adj_nulli_bmi$label[T3_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_LGAvsAGA_adj_nulli_bmi$logFC < 0] <- T3_LGAvsAGA_adj_nulli_bmi$EntrezGeneSymbol[T3_LGAvsAGA_adj_nulli_bmi$adj.P.Val < 0.05 & T3_LGAvsAGA_adj_nulli_bmi$logFC < 0]

#Plot volcanoplot
Volcano_T3_LGAvsAGA_adj_nulli_bmi <- ggplot(data = T3_LGAvsAGA_adj_nulli_bmi, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 15,
                  segment.color = 'grey50',
                  size = 5,
                  color = "black") +
  theme_classic() +
  #labs(tag = "C)") +
  xlab("logFC") +
  ylab("-log10 q-value") +
  scale_y_continuous(limits = c(0, 11), 
                     breaks = seq(0, 12, by = 2)) +
  scale_x_continuous(limits = c(-3,2.1), 
                     breaks = seq(-3, 2, by = 1)) +
  scale_color_manual(values = c("grey", "brown2")) +
  theme(legend.position = "none") +
  ggtitle("Visit 3", subtitle = "Week 28-34") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

#Combine volcanoplots only
LGAvsAGA_Volcanoplots_adj_nulli_bmi <- grid.arrange(Volcano_T1_LGAvsAGA_adj_nulli_bmi, Volcano_T2_LGAvsAGA_adj_nulli_bmi, Volcano_T3_LGAvsAGA_adj_nulli_bmi, 
                                                    nrow = 1, widths = c(1.0, 1.0, 1.0), 
                                                    top = textGrob("LGA vs. AGA", gp = gpar(fontsize = 15)))

ggsave(filename= "Volcano_LGA_AGA_adj_nulli_bmi.png",
       plot = LGAvsAGA_Volcanoplots_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 38,
       height = 13,
       units = "cm")



# Heat map, LGA vs AGA ----------------------------------------------------
#Get a list of the significant proteins, LGA vs AGA
T1_LGAvsAGA <- read_excel("../T1_LGA_AGA_adj_nulli_bmi.xlsx")
T1_LGAvsAGA_sig <- T1_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
T1_LGAvsAGA_sig <- T1_LGAvsAGA_sig$AptName

T2_LGAvsAGA <- read_excel("../T2_LGA_AGA_adj_nulli_bmi.xlsx")
T2_LGAvsAGA_sig <- T2_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
T2_LGAvsAGA_sig <- T2_LGAvsAGA_sig$AptName

T3_LGAvsAGA <- read_excel("../T3_LGA_AGA_adj_nulli_bmi.xlsx")
T3_LGAvsAGA_sig <- T3_LGAvsAGA %>%
  filter(adj.P.Val < 0.05)
T3_LGAvsAGA_sig <- T3_LGAvsAGA_sig$AptName


LGAvsAGA_SigProteins <- c(T1_LGAvsAGA_sig, T2_LGAvsAGA_sig, T3_LGAvsAGA_sig)
LGAvsAGA_SigProteins <- unique(LGAvsAGA_SigProteins)

#Format a data frame with logFC
T1_LGAvsAGA <- T1_LGAvsAGA %>%
  dplyr::select(c(logFC, AptName, EntrezGeneSymbol, TargetFullName)) %>%
  dplyr::rename("Visit 1" = logFC)

T2_LGAvsAGA <- T2_LGAvsAGA %>%
  dplyr::select(c(logFC, AptName, EntrezGeneSymbol, TargetFullName)) %>%
  dplyr::rename("Visit 2" = logFC)

T3_LGAvsAGA <- T3_LGAvsAGA %>%
  dplyr::select(c(logFC, AptName, EntrezGeneSymbol, TargetFullName)) %>%
  dplyr::rename("Visit 3" = logFC)

LogFC_LGAvsAGA <- merge(T1_LGAvsAGA, T2_LGAvsAGA, by = "AptName", all = "TRUE") %>%
  merge(T3_LGAvsAGA, by = "AptName", all = "TRUE") 

#Filter on the significant proteins 
LogFC_LGAvsAGA <- filter(LogFC_LGAvsAGA, grepl(paste(LGAvsAGA_SigProteins, collapse = "|"), AptName)) #Filter on common SomaIds
LogFC_LGAvsAGA <- LogFC_LGAvsAGA %>%
  dplyr::mutate(EntrezGeneSymbol2 = paste0("(", EntrezGeneSymbol,  ")"))
LogFC_LGAvsAGA$TargetFullName_EntrezGeneSymbol <- paste(LogFC_LGAvsAGA$TargetFullName, LogFC_LGAvsAGA$EntrezGeneSymbol2)

#Check duplicates, Targetfullname
LGAvsAGA_TargetFullName_EntrezGeneSymbol <- LogFC_LGAvsAGA$TargetFullName_EntrezGeneSymbol #Extract to check for duplicate names
LGAvsAGA_dup <- LGAvsAGA_TargetFullName_EntrezGeneSymbol[duplicated(LGAvsAGA_TargetFullName_EntrezGeneSymbol)] #Extract the duplicates

LogFC_LGAvsAGA <- LogFC_LGAvsAGA %>% 
  remove_rownames() %>% 
  column_to_rownames(var="TargetFullName_EntrezGeneSymbol") %>%
  #column_to_rownames(var="EntrezGeneSymbol2") %>%
  select(`Visit 1`, `Visit 2`, `Visit 3`)

#Set value zero to be white color
Breaks <- c(seq(min(LogFC_LGAvsAGA), 0, length.out=ceiling(100/2) + 1), 
            seq(max(LogFC_LGAvsAGA)/100, max(LogFC_LGAvsAGA), length.out=floor(100/2)))


heatmap_LGAvsAGA <- LogFC_LGAvsAGA %>%
  pheatmap(border_color = NA,
           breaks = Breaks,
           clustering_method = "complete",
           cluster_cols = FALSE,
           show_rownames = TRUE,
           color = colorRampPalette(c("cornflowerblue", "white", "brown2"))(100),
           fontsize_row = 14,
           fontsize_col = 17,
           treeheight_row = 0,
           cellwidth = 20,
           cellheight = 18,
           main = "")

heatmap_LGAvsAGA <- LogFC_LGAvsAGA %>%
  pheatmap(border_color = NA,
           breaks = Breaks,
           clustering_method = "complete",
           cluster_cols = FALSE,
           show_rownames = FALSE,
           color = colorRampPalette(c("cornflowerblue", "white", "brown2"))(100),
           fontsize_row = 14,
           fontsize_col = 17,
           treeheight_row = 50,
           cellwidth = 20,
           cellheight = 5,
           main = "")


ggsave(filename= "Heatmap_LGAvsAGA_DE_proteins.png",
       plot = heatmap_LGAvsAGA,
       device = "png",
       path = "../",
       width = 10,
       height = 30,
       units = "cm")



# Combine plots LGA vs AGA ------------------------------------------------


#Combine volcano and venn
LGAvsAGA_Volcanoplots_venn_adj_nulli_bmi <- ggarrange(Volcano_T1_LGAvsAGA_adj_nulli_bmi, Volcano_T2_LGAvsAGA_adj_nulli_bmi, 
                                                      Volcano_T3_LGAvsAGA_adj_nulli_bmi, vennplot_T1T2T3_LGAvsAGA,
                                                      ncol = 2, nrow = 2, widths = c(1.0, 1.0),
                                                      labels = c("A)", "B)", "C)", "D)"))

LGAvsAGA_Volcanoplots_venn_adj_nulli_bmi <- annotate_figure(LGAvsAGA_Volcanoplots_venn_adj_nulli_bmi, 
                                                            top = textGrob("LGA vs. AGA", gp = gpar(fontsize = 15)))

ggsave(filename= "Volcano_venn_LGA_AGA_adj_nulli_bmi.png",
       plot = LGAvsAGA_Volcanoplots_venn_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 23,
       height = 20,
       #width = 45,
       #height = 18,
       units = "cm",
       bg = "white")


#Combine volcano, venn and heatmap
LGAvsAGA_heatmap_gr <- as.grob(heatmap_LGAvsAGA)

Column_1 <- ggarrange(Volcano_T1_LGAvsAGA_adj_nulli_bmi, 
                      Volcano_T3_LGAvsAGA_adj_nulli_bmi, 
                      ncol = 1, labels = c("A)", "C)"), font.label = list(size = 20))

Column_2 <- ggarrange(Volcano_T2_LGAvsAGA_adj_nulli_bmi, 
                      vennplot_T1T2T3_LGAvsAGA, 
                      ncol = 1, labels = c("B)", "D)"), font.label = list(size = 20))


Column_3 <- ggarrange(LGAvsAGA_heatmap_gr,
                      ncol = 1, labels = c("E)"), font.label = list(size = 20))

LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- ggarrange(Column_1, Column_2, Column_3, ncol = 3, 
                                                              nrow = 1, widths = c(1, 1, 0.7))


#LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- annotate_figure(LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi, 
#                                                                    top = textGrob("LGA vs. AGA", gp = gpar(fontsize = 25)))

ggsave(filename= "Volcano_venn_heatmap_LGA_AGA_adj_nulli_bmi.png",
       plot = LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 40,
       height = 30,
       units = "cm",
       bg = "white")

# Combine plots LGA vs AGA , Figure for biorender------------------------------------------------


#Combine volcano, venn and heatmap
LGAvsAGA_heatmap_gr <- as.grob(heatmap_LGAvsAGA)

Column_1 <- ggarrange(Volcano_T1_LGAvsAGA_adj_nulli_bmi, Volcano_T2_LGAvsAGA_adj_nulli_bmi, 
                      Volcano_T3_LGAvsAGA_adj_nulli_bmi, vennplot_T1T2T3_LGAvsAGA, 
                      ncol = 1)

Column_2 <- ggarrange(LGAvsAGA_heatmap_gr,
                      ncol = 1)

LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- ggarrange(Column_1, Column_2, ncol = 2, 
                                                              nrow = 1, widths = c(0.7, 1))


#LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi <- annotate_figure(LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi, 
#                                                                    top = textGrob("LGA vs. AGA", gp = gpar(fontsize = 25)))

ggsave(filename= "Volcano_venn_heatmap_LGA_AGA_adj_nulli_bmi_without_PE.png",
       plot = LGAvsAGA_Volcanoplots_venn_heatmap_adj_nulli_bmi,
       device = "png",
       path = "../",
       width = 40,
       height = 50,
       units = "cm",
       bg = "white")




