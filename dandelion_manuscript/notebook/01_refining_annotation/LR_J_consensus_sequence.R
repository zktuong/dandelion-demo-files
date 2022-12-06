setwd('/lustre/scratch117/cellgen/team205/cs42/VDJ_collab_manuscript/')

##### read in sequences from Krzystof
path = '/lustre/scratch117/cellgen/team205/kp9/jupyter/cs42/22-09-30-j-tails/'

j1 = read.csv(paste0(path, 'j+seq.fa'), header = FALSE)
j2 = read.csv(paste0(path, 'j-seq-revcomp.fa'), header = FALSE)
j_combine = rbind(j1, j2)
j_sequence = data.frame (genes  = j_combine$V1[seq(1, nrow(j_combine), by=2)],
                       sequence = j_combine$V1[seq(2, nrow(j_combine), by=2)]
)

# remove > at the front of each gene name
j_sequence$genes = substr(j_sequence$genes, 2,100)
# remove ' reverse complement' at the end of gene name
tmp = j_sequence$genes[51:nrow(j_sequence)]
j_sequence$genes[51:nrow(j_sequence)] = substr(tmp, 1, nchar(tmp)-19)

##### read in list of Js that increase/decrease J multimapper
j_increase = read.csv('csv/j_multimapper_increase.txt', header = FALSE)
#tmp <- strsplit(x=j_increase$V1,split = '\\*')
j_increase_list = j_increase$V1

j_decrease = read.csv('csv/j_multimapper_decrease.txt', header = FALSE)
#tmp <- strsplit(x=j_decrease$V1,split = '\\*')
j_decrease_list = j_decrease$V1

# everything is in j_sequence, with an added column of effect
j_sequence$effect = NA
j_sequence$effect[j_sequence$genes %in% j_increase_list] = 'increase'
j_sequence$effect[j_sequence$genes %in% j_decrease_list] = 'decrease'
# remove IGHJ6
j_sequence = j_sequence[j_sequence$genes != 'IGHJ6',]
write.csv(j_sequence, file='csv/j_sequence_affect_j_multimapper.csv')
write.csv(j_sequence, file='/home/jovyan/mount/gdrive/VDJ_collab/manuscript/supplement/j_sequence_affect_j_multimapper.csv')

### perform logistic regression to look at which sequences are different btw increase vs decrease list
# look at T in position 6 & T in position 17
lr_result = data.frame(n_start = c(6,17,12), n_end = c(6,17,17), nt = c('T','T','GTAAGT'), logOddsRatio = rep(NA,3), p_val = rep(NA, 3))
for (i in 1:nrow(lr_result)){
  position = substr(j_sequence$sequence,lr_result$n_start[i],lr_result$n_end[i])
  is_nt = position==lr_result$nt[i]
  
  fit <- glm(is_nt ~ effect, data=j_sequence, family='binomial')
  
  lr_result$logOddsRatio[i] = summary(fit)$coefficients['effectincrease', 1]
  lr_result$p_val[i] = summary(fit)$coefficients['effectincrease', 4]
}
# pval is 0.019 for position 6, 0.052 for position 17

j_sequence$motif = substr(j_sequence$sequence, 12, 17)
table(j_sequence$motif[j_sequence$effect=='increase'])
table(j_sequence$motif[j_sequence$effect=='decrease'])
