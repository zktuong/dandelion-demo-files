##### volcano plot to look at overlap of CD4 vs CD8 genes in panfetal & panimmune
library(dplyr)
library(ggplot2)
library(ggrepel)
setwd('/lustre/scratch117/cellgen/team205/cs42/VDJ_collab_manuscript/')

panfetal = read.csv('csv/CD4_CD8_LR.csv')
panimmune = read.csv('csv/panimmune_CD4_CD8_LR.csv',row.names = 1)
panfetal$panimmune_estimate = panimmune[panfetal$gene, 'logOddsRatio']
panfetal$panimmune_adjp = panimmune[panfetal$gene, 'adj_pval']

# add a column of sig
panfetal$sig <- "NO"
# if significant in both list
tmp = panfetal$adj_pval < 0.05 & panfetal$panimmune_adjp < 0.05
panfetal$sig[panfetal$logOddsRatio > 0 & panfetal$panimmune_estimate > 0 & tmp] <- 'INCREASE'
panfetal$sig[panfetal$logOddsRatio < 0 & panfetal$panimmune_estimate < 0 & tmp] <- 'DECREASE'

panfetal$sig <- factor(panfetal$sig , levels = c("DECREASE", "NO", "INCREASE"))

# add a column for labelling - only label those that are diffed
panfetal$delabel <- NA
panfetal$delabel[panfetal$sig != "NO"] <- panfetal$gene[panfetal$sig != "NO"]

## plot
pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"LR_CD4_CD8_overlap_volcano_plot.pdf"),width=8, height=4)
ggplot(data=panfetal, aes(x=logOddsRatio, y=-log10(adj_pval), col=sig, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#1F77B4", "black", "#D62728")) +
  xlim(-2.5,2.5) +
  xlab('log(Odds Ratio)')
dev.off()