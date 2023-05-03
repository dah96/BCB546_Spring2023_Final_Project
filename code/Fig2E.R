#Load data and packages
library(tidyverse)
library(data.table)
Ribodiff.ZIKV <- fread("data/GSE207347_A1B1_vs_A2B2_ZIKV_ribodiff_name.txt.gz")
Ribodff.DENV <- fread("data/GSE207347_A1B1_vs_A3B3_DENV_ribodiff_name.txt.gz")
DESeq2.ZIKV <- fread("data/GSE207347_ZIKV_DESeq2_result_name.txt.gz")
DESeq2.DENV <- fread("data/GSE207347_DENV_DESeq2_result_name.txt.gz")

# Figure 2 E


#Grabbing tRNA Synthetases
ZIKV.trna <- subset(DESeq2.ZIKV, na %in% c("CARS", "YARS", "WARS", "AARS", "MARS", "SARS", "GARS", "VARS", "TARS", "IARS"))
DENV.trna <- subset(DESeq2.DENV, na %in% c("CARS", "YARS", "WARS", "AARS", "MARS", "SARS", "GARS", "VARS", "TARS", "IARS"))

#Combining ZIKV and DENV data
ZIKV.trna$virus <- "ZIKV"
DENV.trna$virus <- "DENV"
complete.trna <- rbind(ZIKV.trna, DENV.trna)

ggplot(complete.trna, aes(x=na, y=log2FoldChange, label=log2FoldChange)) +
  geom_bar(stat = "identity", aes(fill=virus), width=0.75, show.legend = FALSE) +
  ggtitle("Aminoacyl tRNA Sythetases") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab("Gene Name") + ylab("RNA log2 Fold Change") +
  facet_wrap(~virus)
#So weirdly, the paper shows TARS decreasing for DENV, but examining the data shows otherwise? Maybe they used a different variant of TARS i.e. TARSL2?
#Also can't seem to adjust color of the bars.
