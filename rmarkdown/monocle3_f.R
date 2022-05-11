#set path
setwd("/fs/ess/PCON0022/guoqi/Spinal cord")
#install packages
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("satijalab/seurat-wrappers")
#load packages
library(SeuratWrappers)
library(monocle3)


#load data
#monocle
genemeta<-data.frame(gene_short_name=rownames(merge[["RNA"]]@meta.features))
rownames(genemeta)<-rownames(merge@assays$RNA@meta.features)
a<-new_cell_data_set(merge@assays$RNA@counts,
                     cell_metadata = merge@meta.data,
                     gene_metadata = genemeta)
#PCA OR LCA
a <- preprocess_cds(a, num_dim = 10)
a <- align_cds(a, alignment_group = "orig.ident")
a<-reduce_dimension(a,umap.metric="manhattan",umap.n_neighbors=20,
                    umap.min_dist = 0.1)#like findclusters, partition is called seurat clusters in seurat
a <- cluster_cells(a)
a <- learn_graph(a)

selfroot_cds <- order_cells(a, root_pr_nodes=get_earliest_principal_node(a))
setwd("/fs/ess/PCON0022/guoqi/Spinal cord/result/Monocle_prommatical_root")
saveRDS(selfroot_cds,"monocleobject_monocleumap_specificparameter_prommaticalroot.rds")
s1<-plot_cells(selfroot_cds,color_cells_by = "pseudotime",show_trajectory_graph = T,
               trajectory_graph_segment_size = 0.5,
               label_leaves=F,
               label_branch_points =F,
               label_cell_groups = F,
               graph_label_size=2)

ggsave(
  plot=s1,
  filename ="./specific_proroot_pseudotime.tiff",
  device = "tiff",
  dpi = 150,
  width = 15,
  height = 10,
  units = "in"
)
s2<-plot_cells(selfroot_cds,color_cells_by = "orig.ident",show_trajectory_graph = T,
               trajectory_graph_segment_size = 0.4,
               label_leaves=F,
               label_branch_points =F,
               label_cell_groups = F,
               graph_label_size=1.5)
ggsave(
  plot=s2,
  filename ="./specific_proroot_sample.tiff",
  device = "tiff",
  dpi = 150,
  width = 15,
  height = 10,
  units = "in"
)

s3<-plot_cells(selfroot_cds,color_cells_by = "Celltype",show_trajectory_graph = T,
               trajectory_graph_segment_size = 0.4,
               label_leaves=F,
               label_branch_points =F,
               label_cell_groups = F,
               graph_label_size=1.5)
ggsave(
  plot=s3,
  filename ="./specific_proroot_celltype.tiff",
  device = "tiff",
  dpi = 150,
  width = 15,
  height = 10,
  units = "in"
)

#specific gene
s4<-plot_cells(selfroot_cds, genes=c("P2ry12", "Siglech", "Gpnmb", "Spp1"),
               trajectory_graph_segment_size = 0.4,
               label_leaves=F,
               label_branch_points =F,
               label_cell_groups = F,
               graph_label_size=1.5)
ggsave(
  plot=s4,
  filename ="./specific_proroot_specificgene.tiff",
  device = "tiff",
  dpi = 150,
  width = 15,
  height = 10,
  units = "in"
)
#density
selfroot_cds_pseudo<-selfroot_cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
selfroot_cds_pseudo<-data.frame()
saveRDS(selfroot_cds_pseudo,"selfroot.rds")
library(ggpubr)
#add cell information
self_data<-data.frame(row.names = names(selfroot_cds_pseudo),pseud=selfroot_cds_pseudo,group=rep(0,length(selfroot_cds_pseudo)))
mdm<-colnames(merge)[merge$Celltype=="MDM"]
microglia<-colnames(merge)[merge$Celltype=="microglia"]
self_data[mdm,2]<-"mdm"
self_data[microglia,2]<-"microglia"
s<-ggdensity(self_data, x = "pseud", rug = TRUE,
             color = "group", palette = c("#00AFBB", "#E7B800"),fill = "group")+xlim(0,35)
ggsave(
  plot=s,
  filename ="psedutime_selfroot.tiff",
  device = "tiff",
  dpi = 150,
  width = 15,
  height = 10,
  units = "in"
)

#choose intermediate pseudotime which has lowest density
mdm<-self_data[1:6791,]
m<-plot(density(mdm$psedu))
mdm_d<-density(mdm$psedu)
mdm_d<-data.frame(mdm_d$x,mdm_d$y)
10.582074  21.94708
microglia<-self_data[6792:nrow(self_data),]
mi<-plot(density(microglia$psedu))
microglia_d<-density(microglia$psedu)
microglia_d<-data.frame(microglia_d$x,microglia_d$y)
20.65509




