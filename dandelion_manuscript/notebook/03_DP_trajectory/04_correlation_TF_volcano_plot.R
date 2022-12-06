library(dplyr)
library(ggplot2)
library(ggrepel)
setwd('/lustre/scratch117/cellgen/team205/cs42/VDJ_collab_manuscript/')

results = read.csv('csv/abtentry_cor_result_tf.csv',row.names = 1)

# add a column of NAs
results$significantly_correlated <- "NO"
# if cor_tcr > 0.1 and adjp_tcr < 0.05, set as "UP" 
results$significantly_correlated[results$cor_tcr > 0.1 & results$adjp_tcr < 0.05] <- "UP"
# if cor_tcr <-0.1 and adjp_tcr < 0.05, set as "DOWN"
results$significantly_correlated[results$cor_tcr < -0.1 & results$adjp_tcr < 0.05] <- "DOWN"

results$delabel <- NA
results$delabel[results$significantly_correlated != "NO"] <- rownames(results)[results$significantly_correlated != "NO"]

pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"TF_volcano_plot.pdf"),width=8, height=4)
ggplot(data=results, aes(x=cor_tcr, y=-log10(adjp_tcr), col=significantly_correlated, label=delabel)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("#1F77B4", "black", "#D62728")) +
  xlab('correlation coefficient')+
  ylab('-log10(adjusted p value)')
dev.off()
