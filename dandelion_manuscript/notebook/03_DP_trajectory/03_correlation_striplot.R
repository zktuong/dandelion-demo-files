library(ggplot2)

setwd('/lustre/scratch117/cellgen/team205/cs42/VDJ_collab_manuscript/')

results = read.csv('~/mount/gdrive/VDJ_collab/manuscript/supplement/abtentry_cor_result.csv',row.names = 1)

resu_melt = data.frame()
for (con in c('tcr','gex')){
  res = results[,c(paste0('cor_',con),paste0('adjp_',con))]
  colnames(res) = c('cor','adjp')
  res$pseudotime = con
  res$genes = rownames(res)
  
  resu_melt = rbind(resu_melt, res)
}
resu_melt$pseudotime = factor(resu_melt$pseudotime, levels = c('tcr','gex'))


genes = c('CD8A','CD8B','CD40LG','ITM2A','CD4','RUNX3','ZBTB7B','TOX','GATA3')

resu_melt$label = NA
resu_melt$label[resu_melt$gene %in% genes] = resu_melt$gene[resu_melt$gene %in% genes]
resu_melt$significantly_correlated <- 'others'
resu_melt$significantly_correlated[resu_melt$cor > 0 & resu_melt$adjp < 0.05 & resu_melt$gene %in% genes] <- "POSITIVE"
resu_melt$significantly_correlated[resu_melt$cor < 0 & resu_melt$adjp < 0.05 & resu_melt$gene %in% genes] <- "NEGATIVE"
resu_melt$significantly_correlated[resu_melt$adjp >= 0.05 & resu_melt$gene %in% genes] <- 'NO'
resu_melt$significantly_correlated = factor(resu_melt$significantly_correlated, levels = c('others','NO','NEGATIVE','POSITIVE'))

ggplot(data=resu_melt, aes(x=pseudotime, y=cor,col=significantly_correlated), label = label) + 
  geom_jitter(position=position_jitter(0.1), cex=0.8)+
  scale_color_manual(values=c( 'black',"#D89000","#0072B2","#D62728")) +
  theme_minimal() +
  geom_text_repel(label=resu_melt$label)+
  xlab('')+
  ylab('correlation coefficient')+ 
  scale_x_discrete(labels=c("tcr" = "VDJ pseudotime", "gex" = "GEX pseudotime"))

# only plot background dots 
resu_melt1 = resu_melt
resu_melt1$significantly_correlated[resu_melt1$significantly_correlated %in% c('NO','NEGATIVE','POSITIVE')] = NaN 
#resu_melt1$significantly_correlated = factor(resu_melt1$significantly_correlated, levels = c('NO','NEGATIVE','POSITIVE'))

pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"cor_volcano_plot_background.pdf"),width=6, height=6)
ggplot(data=resu_melt1, aes(x=pseudotime, y=cor,col=significantly_correlated), label = label) + 
  geom_jitter(position=position_jitter(0.1), cex=0.8)+
  scale_color_manual(values=c('grey')) +
  theme_minimal() +
  geom_text_repel(label=resu_melt1$label)+
  xlab('')+
  ylab('correlation coefficient')+ 
  scale_x_discrete(labels=c("tcr" = "VDJ pseudotime", "gex" = "GEX pseudotime"))
dev.off()

# put others as NaN and plot the colored dots
resu_melt2 = resu_melt
resu_melt2$significantly_correlated[resu_melt2$significantly_correlated=='others'] = NaN 
resu_melt2$significantly_correlated = factor(resu_melt2$significantly_correlated, levels = c('NO','NEGATIVE','POSITIVE'))

pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"cor_volcano_plot.pdf"),width=6, height=6)
ggplot(data=resu_melt2, aes(x=pseudotime, y=cor,col=significantly_correlated), label = label) + 
  geom_jitter(position=position_jitter(0.1), cex=0.8)+
  scale_color_manual(values=c("#D89000","#0072B2","#D62728")) +
  theme_minimal() +
  geom_text_repel(label=resu_melt2$label)+
  xlab('')+
  ylab('correlation coefficient')+ 
  scale_x_discrete(labels=c("tcr" = "VDJ pseudotime", "gex" = "GEX pseudotime"))
dev.off()