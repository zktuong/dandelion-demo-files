### Plotting nhood graph with proper legend

# load libraries
install.packages('ggraph')
library(ggraph)
library(igraph)

# use seaborn colorblind pallette for plotting
sb_colorblind = c('#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9')

# set directory
setwd('/home/jovyan/VDJ_collab_manuscript/')

### plot the initial neighbourhood sampled from GEX embedding

# load csv saved from DP_VDJ_feature_space.ipynb
nhood_conn = read.csv('csv/nhood_conn_dn.csv', row.names = 1)
layout = data.matrix(read.csv('csv/layout_dn.csv', row.names = 1))
obs = read.csv('csv/nhood_obs_dn.csv', row.names = 1)

# change nhood_conn to matrix
adjmatrix = data.matrix(nhood_conn)

# make adjacency matrix graph
nh_graph = graph_from_adjacency_matrix(
  adjmatrix,
  mode = 'directed',
  weighted = TRUE,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)

# set parameters for plotting
node_stroke= 0.3
size_range=c(0.5,3)
size = obs$Nhood_size
colour_by = as.factor(obs$nhood_annotation)
colour_by <- factor(colour_by, levels = c('DN(early)_T','DN(P)_T', 'DN(Q)_T', 'DP(P)_T', 'DP(Q)_T',
                                          'ILC2', 'ILC3','CYCLING_ILC','NK', 'CYCLING_NK'))

# plotting
pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"nhood_umap_dn.pdf"),width=9, height=6)
ggraph(nh_graph, layout = layout) +
  geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
  geom_node_point(aes(fill = colour_by, size = size), shape=21, stroke=node_stroke) +
  scale_size(range = size_range, name="Nhood size") +
  scale_edge_width(range = c(0.2,3), name="overlap size") +
  theme_classic(base_size=14) +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())+
  scale_fill_manual(values=sb_colorblind)
dev.off()

### plot the new neighbourhood VDJ feature space
# PCA
layout_vdj_pca = data.matrix(read.csv('csv/layout_vdj_pca_dn.csv', row.names = 1))

pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"nhood_vdj_pca_dn.pdf"),width=9, height=6)
ggraph(nh_graph, layout = layout_vdj_pca) +
  geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
  geom_node_point(aes(fill = colour_by, size = size), shape=21, stroke=node_stroke) +
  scale_size(range = size_range, name="Nhood size") +
  scale_edge_width(range = c(0.2,3), name="overlap size") +
  theme_classic(base_size=14) +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())+
  scale_fill_manual(values=sb_colorblind)
dev.off()

# UMAP
layout_vdj_umap = data.matrix(read.csv('csv/layout_vdj_umap_dn.csv', row.names = 1))

pdf(paste0('/home/jovyan/mount/gdrive/VDJ_collab/plots_output/chenqu_jhub/',"nhood_vdj_umap_dn.pdf"),width=9, height=6)
ggraph(nh_graph, layout = layout_vdj_umap) +
  geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
  geom_node_point(aes(fill = colour_by, size = size), shape=21, stroke=node_stroke) +
  scale_size(range = size_range, name="Nhood size") +
  scale_edge_width(range = c(0.2,3), name="overlap size") +
  theme_classic(base_size=14) +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())+
  scale_fill_manual(values=sb_colorblind)
dev.off()
