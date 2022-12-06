library(dplyr)
library(ggplot2)
library(ggrepel)
setwd('/lustre/scratch117/cellgen/team205/cs42/VDJ_collab_manuscript/')

results = read.csv('csv/LR_results.csv',row.names = 1)
results_pbmc = read.csv('csv/LR_results_combined.csv',row.names = 1)
sig_pbmc_list = rownames(results_pbmc[results_pbmc$adj_pval<0.05,])

# add a column of j_multimapper
results$j_multimapper <- "NO"
# if significant in both list (results and results_pbmc), set as "DIFF" 
results$j_multimapper[results$adj_pval < 0.05 & rownames(results) %in% sig_pbmc_list] <- "DIFF"

results$j_multimapper[results$j_multimapper=='DIFF' & results$Estimate > 0] <- 'INCREASE'
results$j_multimapper[results$j_multimapper=='DIFF' & results$Estimate < 0] <- 'DECREASE'

results$j_multimapper <- factor(results$j_multimapper, levels = c("DECREASE", "NO", "INCREASE"))

# add a column for labelling - only label those that are diffed
results$delabel <- NA
results$delabel[results$j_multimapper != "NO"] <- rownames(results)[results$j_multimapper != "NO"]

# update labels to remove j_call_leftmost
for (i in 1:nrow(results)){
  label = results$delabel[i]
  if (!is.na(label) & startsWith(label, 'j_call_leftmost')){
    results$delabel[i] = substr(label, 16,nchar(label)+1)
  }
}
# change 'v_gene_presentTrue' to 'v_gene_present'
results$delabel[results$delabel == 'v_gene_presentTrue'] = 'v_gene_present'

## plot
pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"LR_volcano_plot.pdf"),width=8, height=4)
ggplot(data=results, aes(x=Estimate, y=-log10(adj_pval), col=j_multimapper, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#1F77B4", "black", "#D62728")) +
  xlab('log(Odds Ratio)')+
  ylim(0,250)
dev.off()

##### write the list for Kris # all the J genes that were significant in the panfetal list 
diff_j_list = results[results$adj_pval<0.05 & startsWith(rownames(results), 'j_call_leftmost'),]
diff_j_list$label = rownames(diff_j_list)
for (i in 1:nrow(diff_j_list)){
  label = diff_j_list$label[i]
  diff_j_list$label[i] = substr(label, 16,nchar(label)+1)
}
up_list = diff_j_list$label[diff_j_list$Estimate>0]
down_list = diff_j_list$label[diff_j_list$Estimate<0]
# remove *01 isotype informaiton
for (i in 1:length(up_list)){
  up_list[i] = strsplit(up_list[i], split='*',fixed=TRUE)[[1]][1]
}
for (i in 1:length(down_list)){
  down_list[i] = strsplit(down_list[i], split='*',fixed=TRUE)[[1]][1]
}
# remove entries that were both in increase or decrease list - only IGHJ6 is 
up_list_update = unique(up_list[!up_list %in% down_list])
down_list_update = unique(down_list[!down_list %in% up_list])
cat(up_list_update, file='csv/j_multimapper_increase.txt', sep='\n')
cat(down_list_update, file='csv/j_multimapper_decrease.txt', sep='\n')

# 26 out of 27 significant Js from PBMC data are also significant in panfetal data
diff_j_list_pbmc = results_pbmc[results_pbmc$adj_pval<0.05 & startsWith(rownames(results_pbmc), 'j_call_leftmost'),]
rownames(diff_j_list_pbmc) %in% rownames(diff_j_list)

# and all are of the same direction i.e. increase or decrease J multimappers
rownames(diff_j_list_pbmc)[diff_j_list_pbmc$Estimate<0] %in% rownames(diff_j_list)[diff_j_list$Estimate<0]
rownames(diff_j_list_pbmc)[diff_j_list_pbmc$Estimate>0] %in% rownames(diff_j_list)[diff_j_list$Estimate>0]
