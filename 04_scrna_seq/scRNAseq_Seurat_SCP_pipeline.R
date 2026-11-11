.libPaths("/home/data/t020406/R/x86_64-pc-linux-gnu-library/4.2/")
library(data.table)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(plyr)
library(limma)
library(AUCell)
library(SeuratData)
library(patchwork)
library(monocle)
library(SeuratObject)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
options(future.globals.maxSize = 1e9)
library(harmony)
library(cowplot)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(SCENIC)
library(ggplot2)
#library(monocle3)
library(dplyr)
library(ggrepel)
library(future)
library(scCustomize)
library(circlize)
library(ComplexHeatmap)
library(VennDiagram)
library(scales)
library(ggnewscale)
library(gtable)
library(BiocGenerics)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpointdensity)
library(DoubletFinder)
setwd("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/")


source("/home/data/t020406/R/WXY_R_function.R")

#if (!require("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("zhanghao-njmu/SCP")

#options(reticulate.conda_binary = "/home/data/t020406/miniconda3/condabin/conda")
#SCP::PrepareEnv()


##K5 activate TF####
k5_ac_ETP_tf=read.csv("iK5_Day2_activate_TFs.csv")
pass_gene=k5_ac_ETP_tf$X

k5_ac_ETP_tf_invivo=read.csv("iK5_Day2_activate_TFs_invivo.csv")
pass_gene=k5_ac_ETP_tf_invivo$X

k5_bind_ac=read.table("K5_L_bind_active_anno.bed")
k5_bind_ac=k5_bind_ac[k5_bind_ac$V5=="Bind","V4"]


GO_fun




mm_tf=read.table("mm_mgi_tfs.txt")

k5_bind_ac_tf=intersect(k5_bind_ac,mm_tf$V1)

test$is_TF="No"
test[test$gene%in%mm_tf$V1,"is_TF"]="Yes"

test$is_K5_bind_and_activate="No"
test[test$gene%in%k5_bind_ac_tf,"is_K5_bind_and_activate"]="Yes"

test$is_K5_activate="No"
test[test$gene%in%k5_bind_ac$V4,"is_K5_activate"]="Yes"

#TE ICM PGC regulon####

ICM_TE_regulon=read.csv("ICM_TE2_sorted_regulon.csv",row.names = 1)
PGC_regulon=read.csv("PGC_sorted_regulon.csv",row.names = 1)


ICM_TE_PGC_regulon=na.omit(unique(as.charater(unlist(ICM_TE_regulon[,c(1,2,4,6,8)])),PGC_regulon[,2]))

k5_bind_ac_tf_regulon=intersect(k5_bind_ac_tf,ICM_TE_PGC_regulon)

write.csv(k5_bind_ac_tf_regulon,"k5_bind_ac_tf_regulon_fin.csv")

##SCP####
library(SCP)
library(BiocParallel)
register(MulticoreParam(workers = 30, timeout = 1000000000,progressbar = TRUE))

CellDimPlot(
  srt = seurat, group.by = c("sample"),
  reduction = "UMAP", theme_use = "theme_blank"
)

FeatureDimPlot(
  srt = seurat, features = c("Klf5"),
  reduction = "UMAP", theme_use = "theme_blank"
)



raw_seurat=seurat

table(seurat$sample)
Epi=subset(seurat,sample=="EpiSCs")
cutoff=quantile(Epi@assays$RNA@data["Klf5",],0.95)
seurat=seurat[,seurat@assays$RNA@data["Klf5",]>cutoff]
seurat=subset(seurat,sample!="EpiSCs")
seurat=merge(seurat,Epi)
table(seurat$sample)
Epi=NULL
gc()


Idents(raw_seurat) <- "sample"
# Downsample the number of cells per identity class
seurat_merged<- subset(x = raw_seurat, downsample = 11000)
#seurat_merged=seurat_merged[,seurat_merged$sample!="EpiSCs"]
table(seurat_merged$sample)

pancreas_sub1=pancreas_sub

pancreas_sub <- Standard_SCP(srt = seurat_merged)
pancreas_sub=monocle
gc()

pancreas_sub <- NormalizeData(pancreas_sub)
all.genes <- rownames(pancreas_sub)
pancreas_sub <- ScaleData(pancreas_sub, features = all.genes)
pancreas_sub <- FindVariableFeatures(pancreas_sub, selection.method = "vst", nfeatures = 2000)

pancreas_sub <- RunPCA(pancreas_sub, features = VariableFeatures(object = pancreas_sub))
 
# Determine percent of variation associated with each PC
pct <- pancreas_sub[["pca"]]@stdev / sum(pancreas_sub[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

pancreas_sub <- FindNeighbors(pancreas_sub, dims = 1:pcs)


###Cluster####
for(seed in 50:100){
  for(re in 1.8){
    filename1=paste("K5_D0-5_cl_UMAP_sample_seed",seed,"re",re,".png",sep="_")
pancreas_sub <- FindClusters(pancreas_sub, resolution =re, random.seed = seed)
CellDimPlot(
  srt = pancreas_sub, group.by = c("seurat_clusters","sample"),
  reduction = "umap", label = T,theme_use = "theme_blank"
)
ggsave(
  filename = filename1, # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 20,             # 宽
  height = 11,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)
system.time(2)
  }
}

pancreas_sub <- FindClusters(pancreas_sub, resolution =2, random.seed =5)
CellDimPlot( pancreas_sub, group.by = c("seurat_clusters"),
  reduction = "umap", label = T)

color = c(colorRampPalette(colors = rev(brewer.pal(n = 10, name = "RdBu"))[1:5])((50)),
          colorRampPalette(colors = rev(brewer.pal(n = 10, name = "RdBu"))[6:10])((50)))
FeaturePlot(pancreas_sub, features = c("ESC_scale","TSC_scale","PGCLC_scale","Primed_scale"), 
               reduction = "umap",min.cutoff = -2,max.cutoff = 2)&
  scale_color_gradientn(colors = color,  limits = c(-2, 2))


pancreas_sub <- FindClusters(pancreas_sub, resolution =2, random.seed = 5)

###UMAP####
pancreas_sub@meta.data=meta
CellDimPlot(monocle, reduction = "umap",label = T,group.by = c("seurat_clusters","lineage_cl"),
        label.size = 5)


write.csv(pancreas_sub@meta.data,"pancreas_sub_ver5_meta.data.csv")


for(s in 1){
  for(n in 40){
    for(seed in seq(1,100,3)){
    #seurat <- FindClusters(seurat, resolution =re, random.seed = i)
    filename=paste("K5_D0-5_filterK5_UMAP_dist0.14_seed",seed,"n.neighbors",n,"spread",s,".png",sep="_")
    
pancreas_sub=RunUMAP(pancreas_sub,
                     dims = 1:pcs,
                     seed.use = seed,
                     n.neighbors = n,#Larger values will result in more global structure 5 to 50.
                     n.epochs = 200,#200 for large datasets
                     min.dist = 0.14,#Larger values ensure embedded points are moreevenly distributed 0.001 to 0.5
                     spread=s#lager will be mixed
)
CellDimPlot(
  srt = pancreas_sub, group.by = c("sample"),
  reduction = "umap", label = T,theme_use = "theme_blank"
)
ggsave(
  filename = filename, # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 15,             # 宽
  height = 10,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)
Sys.sleep(2)
    }
  }
}


all.genes <- rownames(monocle)
monocle=FindVariableFeatures(monocle,selection.method = "vst",nfeatures = 2000)#Pgc1100
monocle <- ScaleData(monocle, features = all.genes)

# pancreas_sub=RunUMAP(pancreas_sub, 
#                      dims = 1:pcs,
#                      seed.use = 1,
#                      verbose = F,
#                      n.neighbors = 50, #Larger values will result in more global structure 5 to 50.
#                      n.epochs = 200,#200 for large datasets
#                      min.dist = 0.14, #Larger values ensure embedded points are moreevenly distributed 0.001 to 0.5
#                      spread=1#lager will be mixed
# )
pancreas_sub=RunUMAP(pancreas_sub,
                     dims = 1:pcs,
                     seed.use = 3,#37
                     verbose = F,
                     n.neighbors =50,#15 #Larger values will result in more global structure 5 to 50.
                     n.epochs = 200,#200 for large datasets
                     min.dist = 0.2,#0.11 #Larger values ensure embedded points are moreevenly distributed 0.001 to 0.5
                     spread=1#lager will be mixed
)
CellDimPlot(
  srt = pancreas_sub, group.by = c("lineage_cl"),
  reduction = "umap", label = T,theme_use = "theme_blank"
)



FeatureDimPlot(seurat, features = c("ESC_AUC","TSC_AUC","PGCLC_AUC","Primed_AUC","CHIR_ac_1.5_AUC","CHIR_ac_2_AUC"),
                  reduction = "umap", theme_use = "theme_blank")
FeatureDimPlot(monocle, features = c("Ets2"), 
               reduction = "umap", theme_use = "theme_blank")


pancreas_sub=readRDS("pancreas_sub_ver6.rds")


# length(intersect(colnames(pancreas_sub),colnames(random_seurat)))
# random_seurat@meta.data=pancreas_sub@meta.data
# random_seurat@tools=pancreas_sub@tools
# random_seurat@reductions=pancreas_sub@reductions
# 
# pancreas_sub=random_seurat

saveRDS(pancreas_sub,"pancreas_sub_ver6.rds")


AUC_scale=data.frame(ESC_scale=pancreas_sub$ESC_AUC,
           TSC_scale=pancreas_sub$TSC_AUC,
           PGCLC_scale=pancreas_sub$PGCLC_AUC,
           Primed_scale=pancreas_sub$Primed_AUC)
AUC_scale=scale(AUC_scale)
pancreas_sub@meta.data=cbind(pancreas_sub@meta.data,AUC_scale)

pancreas_sub1=subset(pancreas_sub,sample!="EpiSCs")

###Linage####
for(i in 0:38){
end=c(33,34,16,32,22,30,21,27,29)#5
#end=c(12,28,11,31,18,16,21)#5
#end=c(17,22,30,20,26,11,i)#4
#end=c(22,28,12,29,15,26,25)#9
end=end[!duplicated(end)]
filename=paste("K5_D0-5_lineage_ver2_end",5,21,27,".png",sep="_")
pancreas_sub <- RunSlingshot(srt = pancreas_sub, group.by = "seurat_clusters", 
                             reduction = "umap",seed = 2,show_plot = F,
                             start = c(5),end = end )
print(i)
ggsave(
  filename = filename, # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 20,             # 宽
  height = 11,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)
Sys.sleep(2)
}
CellDimPlot(pancreas_sub, group.by = "seurat_clusters", reduction = "umap",
            lineages = paste0("Lineage", 1:17), lineages_span =0.3)
FeatureDimPlot(
  srt = pancreas_sub,
  features = #c("ESC_scale","TSC_scale","PGCLC_scale","Primed_scale"),
    c("Cdx2","Cdcp1"),
  compare_features = TRUE, label = TRUE, label_insitu = TRUE,
  reduction = "umap", theme_use = "theme_blank"
)



FeatureDimPlot(pancreas_sub, features = paste0("Lineage", c( 1:17)), 
               reduction = "umap", theme_use = "theme_blank")


want_lineage= c(2,4,5,9,10,12,14,16)

CellDimPlot(monocle, group.by = "lineage_cl", reduction = "umap",
            lineages_palcolor = c( "#6699CC","#661100",
                                  "#CC6677","#332288" ,"#117733" ,
                                  "#999933","#888888","#AA4499"),
            lineages = paste0("Lineage", want_lineage), lineages_span =0.3)
ggsave(
  filename = "D0-5_K5_lineage_umap_cl.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 8,            # 高
  units = "in",          # 单位
  dpi = 1000              # 分辨率DPI
)




DynamicPlot(
  srt = pancreas_sub, lineages = paste0("Lineage", c(12,10)), group.by = "seurat_clusters",
  features = c("Klf4","Prdm1"),
  compare_lineages = TRUE, compare_features = FALSE
)

FeaturePlot(pancreas_sub,reduction = "umap",
            features =c("Gata3", "Nanog", "Tfap2c", "Prdm1"),
            ncol=2,order=T)


for(i in 7){
  print(want_linage[i:c(i+1)])
  pancreas_sub <- RunDynamicFeatures(srt = pancreas_sub,minfreq = 1,features=VariableFeatures(lineage_sub),
                                     lineages =paste0("Lineage",c("10_12_15_8","10_12_8","10_12")))#features = pass_gene
}


pancreas_sub$Lineage10_12_15_8=na.omit(c(pancreas_sub$Lineage10,pancreas_sub$Lineage12,
                                         pancreas_sub$Lineage15,pancreas_sub$Lineage8))

length(pancreas_sub$Lineage10_12_15_8)
head(pancreas_sub$Lineage10_12_15_8)
max(pancreas_sub$Lineage10_12_15_8)

RunDynamic=
  function(pancreas_sub,want_lineage,study_gene){
    
    # features=c()
    # for(i in want_lineage){
    # cat("\n"," Finding VariableFeatures in all lineage...","\n")
    # ln=na.omit(as.data.frame(pancreas_sub[[paste0("Lineage",c(i))]]))
    # lineage_sub=pancreas_sub[,colnames(pancreas_sub)%in%rownames(ln)]
    # cat(" Cells ",length(lineage_sub[[paste0("Lineage",c(i))]][[1]]),"\n")
    # lineage_sub=FindVariableFeatures(lineage_sub,selection.method = "vst",nfeatures = 2000)
    # features=c(features,VariableFeatures(lineage_sub))
    # }
    # features=unique(features)
    # cat("VariableFeatures in all lineage ",length(features),"\n")
    # cat(" Without ",setdiff(study_gene,features),"\n")
    
for(i in want_lineage){

  
  cat("\n"," Doing...",paste0("Lineage",c(i)),"\n",
        "Output in",paste0("DynamicFeatures_Lineage",c(i),"\n"))
  ln=na.omit(as.data.frame(pancreas_sub[[paste0("Lineage",c(i))]]))
  lineage_sub=pancreas_sub[,colnames(pancreas_sub)%in%rownames(ln)]
  
  cat(" Cells ",length(lineage_sub[[paste0("Lineage",c(i))]][[1]]),"\n")
  lineage_sub=FindVariableFeatures(lineage_sub,selection.method = "vst",nfeatures = 2000)
  cat(" Without ",setdiff(study_gene,VariableFeatures(lineage_sub)),"\n")
  
  
  pancreas_sub <- RunDynamicFeatures(srt = pancreas_sub,minfreq = 5,features=VariableFeatures(lineage_sub),
                                     lineages =paste0("Lineage",c(i)))
  
  cat(" Slimming...","\n","Before:",
      names(pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]]),"\n")
  #print(head(pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]][[2]][,1:6]))
  
  pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]][[2]]=pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]][[2]][,1:2]
  pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]][2:8]=NULL
  
  cat("\n","After:",
      names(pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]]),"\n")
  #print(head(pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]][[2]]))
  
}
    return(pancreas_sub)
}


meta=readRDS("monocle_meta_maked_lineage.rds")



monocle@meta.data=meta
monocle=RunDynamic(monocle,c(#"_9v10_9","_9v10_10","_9v12_9","_9v12_12",
                     #"_9v14_9","_9v14_14","_9v16_9","_9v16_16",
                     "_16v10_16","_16v10_10","_16v12_16","_16v12_12","_16v14_16","_16v14_14",
                     "_14v10_14","_14v10_10","_14v12_14","_14v12_12",
                     "_10v12_10","_10v12_12"
),study_gene=study_gene)


monocle=RunDynamic(monocle,want_lineage,study_gene)




saveRDS(monocle_annoed,"monocle2_new.rds")

monocle=readRDS("monocle2_new.rds")


normalize_pseudotimeta.datanormalize_pseudotime=
  function(pancreas_sub,want_lineage){
    meta=pancreas_sub@meta.data
    
    for(i in want_lineage){
    nor_pseudo=pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]][[2]][,1]
    nor_pseudo=c(nor_pseudo-min(nor_pseudo))/c(max(nor_pseudo)-min(nor_pseudo))
    nor_pseudo=as.data.frame(nor_pseudo)
    colnames(nor_pseudo)=paste0("Lineage",i)
    
    
    meta[[paste0("Lineage",i)]]=NULL
    meta=Btable_insert_into_Atable(meta,nor_pseudo)
    meta[[paste0("Lineage",i)]]=as.numeric(meta[[paste0("Lineage",i)]])
    
    }
    pancreas_sub@meta.data=meta
    
    return(pancreas_sub)
  }



pancreas_sub=normalize_pseudotime(pancreas_sub,want_lineage)


select_feature_in_lineage=
  function(pancreas_sub,want_lineage,cutoff=c(0.1,0.1,0.01)){

    pass_feature=c()
    for(i in want_lineage){
      feature_mat=pancreas_sub@tools[[paste0("DynamicFeatures_Lineage",c(i))]][[1]]
      feature_mat=feature_mat[which(feature_mat$r.sq>cutoff[1] &
                                      feature_mat$dev.expl>cutoff[2] &
                                      feature_mat$padjust<cutoff[3] &
                                      feature_mat$exp_ncells >20),]
      pass_feature=c(pass_feature,feature_mat$features)
      
    }
    pass_feature=pass_feature[!duplicated(pass_feature)]
    return(pass_feature)
  }


pass_feature=select_feature_in_lineage(pancreas_sub,c(10,12))




monocle$lineage_9v10_lineage_cl=NA
meta_monocle=monocle@meta.data
meta_monocle[meta_monocle$lineage_cl%in%c(31),"lineage_9v10_lineage_cl"]="Lineage_9_cellfate"
meta_monocle[meta_monocle$lineage_cl%in%c(2,18,26,30),"lineage_9v10_lineage_cl"]="Lineage_10_cellfate"
meta_monocle[meta_monocle$lineage_cl%in%c(1),"lineage_9v10_lineage_cl"]="Lineage_prebranch"
monocle@meta.data=meta_monocle

 
c(#"_9v10_9","_9v10_10","_9v12_9","_9v12_12",
  #"_9v14_9","_9v14_14","_9v16_9","_9v16_16",
  "_16v10_16","_16v10_10","_16v12_16","_16v12_12","_16v14_16","_16v14_14",
  "_14v10_14","_14v10_10","_14v12_14","_14v12_12",
  "_10v12_10","_10v12_12")

ht1 <- DynamicHeatmap_modify(
  srt = monocle, slot = "data",study_gene_set=co16,
  lineages = paste0("Lineage", c("_16v10_16","_16v10_10")),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  reverse_ht="Lineage_16v10_16",
  show_row_names=T,
  n_split = 5,
  split_method = "mfuzz",#mfuzz,kmeans-peaktime
  order_by = "peaktime",
  limits=c(-2,2),
  heatmap_palette = "viridis", 
  cell_annotation = c("sample","Lineage_16v10_branch"),
  cell_annotation_palette=c("Paired","Set2"),
  width=2,height=5,raster_by_magick=T
)
ht2 <- DynamicHeatmap_modify(
  srt = monocle, slot = "data",study_gene_set=co16,
  lineages = paste0("Lineage", c("_16v12_16","_16v12_12")),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  reverse_ht="Lineage_16v12_16",
  show_row_names=T,
  n_split = 5,
  split_method = "mfuzz",#mfuzz,kmeans-peaktime
  order_by = "peaktime",
  limits=c(-2,2),
  heatmap_palette = "viridis", 
  cell_annotation = c("sample","Lineage_16v12_branch"),
  cell_annotation_palette=c("Paired","Set2"),
  width=2,height=5,raster_by_magick=T
)
ht3 <- DynamicHeatmap_modify(
  srt = monocle, slot = "data",study_gene_set=co16,
  lineages = paste0("Lineage", c("_16v14_16","_16v14_14")),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  reverse_ht="Lineage_16v14_16",
  show_row_names=T,
  n_split = 5,
  split_method = "mfuzz",#mfuzz,kmeans-peaktime
  order_by = "peaktime",
  limits=c(-2,2),
  heatmap_palette = "viridis", 
  cell_annotation = c("sample","Lineage_16v14_branch"),
  cell_annotation_palette=c("Paired","Set2"),
  width=2,height=5,raster_by_magick=T
)

co16=intersect(intersect(ht1$feature_metadata$features,ht2$feature_metadata$features),ht3$feature_metadata$features)

htco16 <- DynamicHeatmap_modify(
  srt = monocle, slot = "data",study_gene_set=co16,
  lineages = paste0("Lineage", c("_16v10_16","_16v10_10","_16v12_12","_16v14_14")),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  reverse_ht="Lineage_16v10_16",
  show_row_names=T,
  n_split = 5,
  split_method = "mfuzz",#mfuzz,kmeans-peaktime
  order_by = "peaktime",
  limits=c(-2,2),
  heatmap_palette = "viridis", 
  cell_annotation = c("sample","Lineage_16v14_branch"),
  cell_annotation_palette=c("Paired","Set2"),
  width=2,height=5,raster_by_magick=T
)

ht4 <- DynamicHeatmap_modify(
  srt = monocle, slot = "data",study_gene_set=co14,
  lineages = paste0("Lineage", c("_14v10_14","_14v10_10")),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  reverse_ht="Lineage_14v10_14",
  show_row_names=T,
  n_split = 5,
  split_method = "mfuzz",#mfuzz,kmeans-peaktime
  order_by = "peaktime",
  limits=c(-2,2),
  heatmap_palette = "viridis", 
  cell_annotation = c("sample","Lineage_14v10_branch"),
  cell_annotation_palette=c("Paired","Set2"),
  width=2,height=5,raster_by_magick=T
)
ht5 <- DynamicHeatmap_modify(
  srt = monocle, slot = "data",study_gene_set=co14,
  lineages = paste0("Lineage", c("_14v12_14","_14v12_12")),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  reverse_ht="Lineage_14v12_14",
  show_row_names=T,
  n_split = 5,
  split_method = "mfuzz",#mfuzz,kmeans-peaktime
  order_by = "peaktime",
  limits=c(-2,2),
  heatmap_palette = "viridis", 
  cell_annotation = c("sample","Lineage_14v12_branch"),
  cell_annotation_palette=c("Paired","Set2"),
  width=2,height=5,raster_by_magick=T
)


co14=intersect(ht4$feature_metadata$features,ht5$feature_metadata$features)


ht6 <- DynamicHeatmap_modify(
  srt = monocle, slot = "data",study_gene_set=k5_bind_ac_tf,
  lineages = paste0("Lineage", c("_10v12_10","_10v12_12")),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  reverse_ht="Lineage_10v12_10",
  show_row_names=T,
  n_split = 5,
  split_method = "mfuzz",#mfuzz,kmeans-peaktime
  order_by = "peaktime",
  limits=c(-2,2),
  heatmap_palette = "viridis", 
  cell_annotation = c("sample","Lineage_10v12_branch"),
  cell_annotation_palette=c("Paired","Set2"),
  width=2,height=5,raster_by_magick=T
)
#7x9
ht1$plot
ht2$plot
ht3$plot
ht4$plot
ht5$plot
ht6$plot



RunDynamicHeatmap=
  function(monocle,want_lineage,split_method="mfuzz",palcolor=palette_list$Paired){
Lineage_DyHeatmap=list()
for(i in 1:length(want_lineage)){
  print(paste("Runing",paste0("Lineage", want_lineage[i]),"DyHeatmap"))
  monocle@tools[[paste0("DyHeatmap_Lineage", want_lineage[i])]]
  ht <- DynamicHeatmap(
    srt = monocle,heatmap_palette = "RdBu", slot = "data",
    lineages = paste0("Lineage", want_lineage[i]),
    r.sq = 0.05,dev.expl=0.05, padjust =0.01,cell_bins = 4000,
    n_split = 4,min_expcells = round(length(na.omit(monocle@meta.data[[paste0("Lineage",want_lineage[i])]]))*0.05),
    split_method = split_method,order_by = "peaktime",cell_annotation_palcolor = list(palcolor),
    #cluster_rows = T,
    cell_annotation = "sample",limits=c(-2,2),use_raster = F,
    width = 2,height=4
  )
  Lineage_DyHeatmap[[paste0("DyHeatmap_Lineage", want_lineage[i])]]=ht
  #print(ht$plot)
}

return(Lineage_DyHeatmap)
}


Lineage_DyHeatmap1=RunDynamicHeatmap(monocle,c(10,12,16,14,9),palcolor=sample_palette_col[2:6])
Lineage_DyHeatmap2=RunDynamicHeatmap(monocle,c(2,4,5),palcolor=sample_palette_col)


for(i in c(2,4,5)){
panel_fix(Lineage_DyHeatmap2[[paste0("DyHeatmap_Lineage",i)]]$plot, height = 6, width = 4, raster = TRUE, dpi = 1000)
ggsave(
  filename = paste0("DyHeatmap_lineage",i,".pdf"),
  width = 6,             # 宽
  height = 6,            # 高
  units = "in",          # 单位
  dpi = 1000              # 分辨率DPI
)

}

saveRDS(Lineage_DyHeatmap2,"Lineage_DyHeatmap2.rds")


#change to my computer to do cluster GO



ht1 <- DynamicHeatmap_modify(
  srt = pancreas_sub, slot = "data",
  lineages = paste0("Lineage", c(9,12)),
  #lineages = paste0("Lineage", c(10)),
  r.sq = 0.1,dev.expl=0.1, padjust =1e-10,
  #use_fitted = TRUE,
  reverse_ht="Lineage9",
  show_row_names=F,
  split_method = "kmeans",
  limits=c(-3,3),
  n_split=5,
  #feature_split=data.frame(row.names=dynamic_new_anno$gene,
  #                         feature_split=factor(dynamic_new_anno$cluster)),
  cluster_rows = TRUE,
  #reverse_ht = "Lineage15",
  #species = "Mus_musculus", db = "GO_BP", #anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "viridis", cell_annotation = "sample",#cell_annotation_palcolor = list(palette_list$Paired[2:6]),
  #separate_annotation = list("seurat_clusters", c("Nanog", "Gata3")), 
  #separate_annotation_palette = c("Paired", "Set1"),
  #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  #pseudotime_label = 0.7, pseudotime_label_color = "red",
  width=2,height=7,raster_by_magick=T
)


head(sub_seurat_marker)

sub_seurat_marker_dynamic=sub_seurat_marker[sub_seurat_marker$gene%in%ht$feature_metadata$features,]
dim(sub_seurat_marker_dynamic)


sub_seurat_marker_dynamic_new=sub_seurat_marker_dynamic[1,]
for( i in 1:length(levels(factor(sub_seurat_marker_dynamic$gene)))){
  gene=levels(factor(sub_seurat_marker_dynamic$gene))[i]
  sub_mat=sub_seurat_marker_dynamic[sub_seurat_marker_dynamic$gene%in%gene,]
  sub_mat=sub_mat[which.max(sub_mat$avg_log2FC),]
  sub_seurat_marker_dynamic_new=rbind(sub_seurat_marker_dynamic_new,sub_mat)
}
sub_seurat_marker_dynamic_new=sub_seurat_marker_dynamic_new[-1,]


nofind=data.frame(gene=setdiff(ht$feature_metadata$features,sub_seurat_marker$gene),
                  cluster=0)
dynamic_new_anno=rbind(sub_seurat_marker_dynamic_new[,c(7,6)],nofind)


feature_split=data.frame(row.names=dynamic_new_anno$gene,
                         feature_split=dynamic_new_anno$cluster,
                         gene=dynamic_new_anno$gene)
feature_split=feature_split[order(feature_split$feature_split),]
head(feature_split)

htmat=as.data.frame(ht$matrix)


htlist=ht$cell_order
htlist=lapply(htlist, function(x) x=strsplit2(x,"Lineage")[,1] )
venn=get.venn.partitions(htlist)
paste(names(venn$..values..),venn$..count..,venn$..set..)



meta=pancreas_sub@meta.data
htlist=list(lineage9=rownames(meta[!is.na(meta$Lineage9),]),
            lineage10=rownames(meta[!is.na(meta$Lineage10),]),
            lineage12=rownames(meta[!is.na(meta$Lineage12),]),
            lineage14=rownames(meta[!is.na(meta$Lineage14),]),
            lineage16=rownames(meta[!is.na(meta$Lineage16),]))
venn=get.venn.partitions(htlist)
paste(names(venn$..values..),venn$..count..,venn$..set..)


#lineage cluster
meta$lineage_cl=0
for(i in 1:length(venn$..values..)){
  if(venn$..count..[[i]]>10){
  meta[rownames(meta)%in%venn$..values..[[i]],"lineage_cl"]=i
  }
}

#normalize pseudotime
for(i in want_lineage){
  print(paste0("Lineage",i))
  nor_pseudo=na.omit(data.frame(row.names =rownames(meta),
                        pstime=meta[[paste0("Lineage",i)]]))
  nor_pseudo$pstime=c(nor_pseudo$pstime-min(nor_pseudo$pstime))/c(max(nor_pseudo$pstime)-min(nor_pseudo$pstime))
  colnames(nor_pseudo)=paste0("Lineage",i)
  
  meta[[paste0("Lineage",i)]]=NULL
  meta=Btable_insert_into_Atable(meta,nor_pseudo)
  meta[[paste0("Lineage",i)]]=as.numeric(meta[[paste0("Lineage",i)]])
}

table(meta$lineage_cl)

b1=meta[meta$lineage_cl==1,c(36:40)]
0:max(b1)

b2=meta[meta$lineage_cl==2,c(37:40)]
min(b2):max(b2)

b3=meta[meta$lineage_cl==18,c(37:39)]
min(b3):max(b3)

b4=meta[meta$lineage_cl==18,c(37:38)]
min(b4):max(b4)

meta$pseudotime_cl=


sub_seurat=subset(pancreas_sub,lineage_cl!=0)
Idents(sub_seurat)="lineage_cl"
sub_seurat_marker <- FindAllMarkers(sub_seurat, only.pos = TRUE, logfc.threshold = 0.25)

pancreas_sub@meta.data=meta

row_meta=meta
meta$lineage_cl=gsub(23,31,meta$lineage_cl)
meta$lineage_cl=gsub(15,1,meta$lineage_cl)
meta$lineage_cl=gsub(10,2,meta$lineage_cl)
meta$lineage_cl=gsub(22,18,meta$lineage_cl)


ht=readRDS("DynamicHeatmap_lineage12.rds")
#saveRDS(ht,"DynamicHeatmap_lineage12.rds")





htmat$gene=rownames(htmat)
plot_mat=Btable_insert_into_Atable(feature_split,htmat)
plot_mat$gene=NULL
plot_mat$feature_split=NULL
head(plot_mat[,1:10])


  
       
binary_plot_mat=
  function(plot_mat,bins=100){
    plot_mat_bin=plot_mat[,1]
    for(i in c(10,12,15,7,8)){
      print(i)
      sub_htmat=htmat[,which(strsplit2(colnames(htmat),split = "Lineage")[,2]==i)]
      step=round(dim(sub_htmat)[2]/bins)
      
      #sub_htmat_binmat <- sapply(split(1:dim(sub_htmat)[2],
      #                     c(rep(1:c(bins-1),step),rep(bins,length(c(length(rep(1:c(bins-1),step))+1):dim(sub_htmat)[2])))),
      #               function(cells) rowMeans(sub_htmat[,cells]))
      down=sample(1:dim(sub_htmat)[2],1000)
      down=down[order(down)]
      sub_htmat_binmat=sub_htmat[,down]
      plot_mat_bin=cbind(plot_mat_bin,sub_htmat_binmat)
    }
    plot_mat_bin=plot_mat_bin[,-1]
    return(plot_mat_bin)
  }

plot_mat_bin=binary_plot_mat(plot_mat)
plot_mat_bin=Btable_insert_into_Atable(feature_split,plot_mat_bin)
plot_mat_bin$gene=NULL
plot_mat_bin$feature_split=NULL
head(plot_mat_bin[,1:10])
plot_mat_bin[,1:5000]=apply(plot_mat_bin,2,as.numeric)

anno=data.frame(row.names=feature_split$gene,group=feature_split$feature_split)

bk=seq(-3,3,0.1)
pheatmap::pheatmap(plot_mat_bin,
  #ht$matrix[rownames(ht$matrix)%in%formative,],
                   scale="row",
                   breaks = bk,annotation_row =anno,gaps_col=seq(1000,5000,1000),
                   cluster_rows=T,clustering_method="ward.D2",
                   cluster_cols=F,
                   show_rownames = F,
                   show_colnames =F,
                   #color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"))((length(bk))*0.5),
                  #           colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5))
                    color  = c(colorRampPalette(colors = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF" ))((length(bk))*0.5),
                               colorRampPalette(colors = c("#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF"))((length(bk))*0.5))
)



Complexheatmap=function(){}

colum_split_fun=c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000))

bk=seq(-3,3,0.1)
col_fun = colorRamp2(bk,
                     c(colorRampPalette(colors =c("#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF" ))(length(bk)*0.5-1),
                       colorRampPalette(colors =  c("#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF"))(length(bk)*0.5)),space = "RGB")
left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill=brewer.pal(n = length(levels(factor(anno$group))), name = 'Set1')),
                                                 labels = as.character(rownames(anno)), 
                                                 labels_gp = gpar(col = "black", #fontface="bold",
                                                                  fontsize = 8)))

Heatmap(#cluster_mat_scale[,c(1:2,27:34,3:26,35:67)], 
  plot_mat_bin,
  border = T,#是否开启聚类的边框
  border_gp = gpar(lwd=1.5),
  #split = gene_cluster_factor,
  name = "z-score", 
  #heatmap_legend_param = list(at = bk_legend),
  column_title = paste("Total",length(plot_mat_bin[,1]),"genes"),
  column_title_gp = gpar(fontsize=12),
  column_names_gp = gpar(fontsize=10),
  column_split = colum_split_fun,
  row_gap = unit(1.5, "mm"),
  #row_title_gp = gpar(fontsize = anno_font_size),
  #row_names_gp = gpar(col = gene_name_color_mark),
  show_row_names = F,
  cluster_rows = F,
  cluster_columns = F,
  use_raster = T,#raster_device = c("CairoPNG"),
  #raster_quality = 10,raster_resize_mat = TRUE,
  raster_by_magick = TRUE,raster_device_param =  raster_magick_filter = "Gaussian",
  left_annotation=left_annotation,
  col=col_fun)


ht=NULL
gc()


ht_feature=ht$feature_metadata

meta=pancreas_sub@meta.data

data=pancreas_sub@assays$RNA$data
lineage12=pancreas_sub@tools$DynamicFeatures_Lineage12

sub=data[,colnames(data)%in%order2]

order1=rownames(pancreas_sub@tools$DynamicFeatures_Lineage10$raw_matrix)
sub1=sub[,which(colnames(sub)%in%order1)]

study_gene=c(formative,"Prdm1","Dppa3","Nanos3",
             "Sox15","Nanog","Tfap2c","Rhox6","Rhox9")
study_gene=ht_feature$features

epi_meta=meta[meta$sample=="EpiSCs",]
epi=data[,which(colnames(data)%in%rownames(epi_meta))]
epi=epi[rownames(epi)%in%study_gene,]
epi=as.data.frame(epi)
epi$gene=factor(rownames(epi),levels = study_gene)
epi=epi[order(epi$gene),]
epi$gene=NULL


order2=rownames(lineage12$raw_matrix)
sub2=data[,which(colnames(data)%in%order2)]
sub2=sub2[rownames(sub2)%in%study_gene,]
sub2=as.data.frame(sub2)
sub2$gene=factor(rownames(sub2),levels = study_gene)
sub2=sub2[order(sub2$gene),]
sub2$gene=NULL
sub2=sub2[,1:c(dim(sub2)[2]/2)]


study_feature=ht_feature[ht_feature$features%in%study_gene,c(1,8)]


anno=data.frame(features=setdiff(study_gene,study_feature$features),
                feature_split=rep("Outlier",length(setdiff(study_gene,study_feature$features))))
anno=rbind(study_feature,anno)
anno=data.frame(row.names = anno$features,cl=anno$feature_split)

anno=data.frame(row.names = ht_feature$features,cl=ht_feature$feature_split)


plotmat=cbind(epi,sub2)
plotmat$gene=factor(rownames(plotmat),levels = rownames(anno))
plotmat=plotmat[order(plotmat$gene),]
plotmat$gene=NULL
bk=seq(-2,2,0.1)
pheatmap::pheatmap(plotmat[,c(dim(epi)[2]-1000):dim(plotmat)[2]],scale="row",show_colnames =F,breaks = bk,gaps_col=1000,
                   cluster_cols=F,cluster_rows=F,show_rownames = F,annotation_row = anno,
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5))
                   )



ht <- GroupHeatmap(
  srt = pancreas_sub,
  features =pass_feature,
  group.by = c("seurat_clusters", "sample"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("ESC_scale","TSC_scale","PGCLC_scale","Primed_scale"),
  #cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE,cluster_rows = F,show_column_names = T
)
print(ht$plot)

###Linage DEG####



Find_lineage_DEG=
  function(monocle,b1_vec,b2_vec){
    
    
    deg_fin=data.frame( p_val=1,
                         avg_log2FC =0,
                         pct.1 =0,
                         pct.2=0,
                         p_val_adj=1,
                         cluster="A",
                         gene="A",
                         lineage_branch="A")
    for(i in 1:length(b1_vec)){
    b1=b1_vec[i]
    b2=b2_vec[i]
    lineage_branch=paste0("Lineage_",b1,"v",b2,"_branch")
    print(lineage_branch)
    Idents(monocle)=lineage_branch
    
    cl=unique(monocle@meta.data[[lineage_branch]])
    cl=cl[grep("CellFate",cl)]
    
    deg_test1=FindMarkers(monocle,only.pos = T,
                            ident.1=cl[1],
                            ident.2=cl[2],
                            logfc.threshold = 0.5,
                            min.pct = 0.01)
    deg_test2=FindMarkers(monocle,only.pos = T,
                          ident.1=cl[2],
                          ident.2=cl[1],
                          logfc.threshold = 0.5,
                          min.pct = 0.01)
    deg_test1$cluster=cl[1]
    deg_test2$cluster=cl[2]
    
    deg_test1$gene=rownames(deg_test1)
    deg_test2$gene=rownames(deg_test2)
    
    deg_test1$lineage_branch=lineage_branch
    deg_test2$lineage_branch=lineage_branch
    
    deg_test=rbind(deg_test1,deg_test2)
    
    print(head(deg_test,3))
    deg_fin=rbind(deg_fin,deg_test)
    }
    return(deg_fin[-1,])
}



vec=strsplit2(colnames(monocle@meta.data)[grepl("branch",colnames(monocle@meta.data))],"_")[,2]
b1_vec=strsplit2(vec,"v")[,1]
b2_vec=strsplit2(vec,"v")[,2]

deg_fin_mat=Find_lineage_DEG(monocle,b1_vec,b2_vec)



saveRDS(deg_fin_mat,"deg_fin_mat.rds")

deg_fin_mat=readRDS("deg_fin_mat.rds")

table(deg_fin_mat$lineage_branch)

deg_fin_mat_filter=deg_fin_mat[deg_fin_mat$p_val_adj<0.01,]
deg_fin_mat_filter=deg_fin_mat_filter[deg_fin_mat_filter$lineage_branch%in%c("Lineage_16v12_branch"),]
table(deg_fin_mat_filter$cluster)

deg_list1=deg_fin_mat_filter[deg_fin_mat_filter$cluster==unique(deg_fin_mat_filter$cluster)[1],"gene"]
deg_list2=deg_fin_mat_filter[deg_fin_mat_filter$cluster==unique(deg_fin_mat_filter$cluster)[2],"gene"]



sub_monocle=monocle[,grepl("CellFate",monocle$Lineage_10v12_branch)]
table(sub_monocle$Lineage_10v12_branch)

deg_mat1=sub_monocle[rownames(sub_monocle)%in%deg_list1,]
deg_mat2=sub_monocle[rownames(sub_monocle)%in%deg_list2,]

deg_mat=as.data.frame(rbind(deg_mat1@assays$RNA@data,deg_mat2@assays$RNA@data))

bk=seq(-2,2,0.1)
col_fun = colorRamp2(bk,
                     c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"))(length(bk)/2-1),
                       colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))(length(bk)/2)),space = "RGB")
left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill=brewer.pal(n = 2, name = 'Set1')),
                                                 labels = as.character(unique(deg_fin_mat_filter$cluster)), 
                                                 labels_gp = gpar(col = "black", #fontface="bold",
                                                                  fontsize = 8)))

top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill=brewer.pal(n = 2, name = 'Set2')),
                                                 labels = as.character(unique(sub_monocle$Lineage_10v12_branch)), 
                                                 labels_gp = gpar(col = "black", #fontface="bold",
                                                                  fontsize = 8)))

ha = rowAnnotation(link = anno_mark(at = which(rownames(cluster_mat_scale) %in% study_genes),
                                    labels =  rownames(cluster_mat_scale)[which(rownames(cluster_mat_scale) %in% study_genes)], 
                                    labels_gp = gpar(fontsize = 8,fontface = "italic")))


deg_mat_scale=t(scale(t(deg_mat)))
p=Heatmap( #cluster_mat_scale[,c(1:2,27:34,3:26,35:67)], 
  deg_mat_scale,
  border = T,#是否开启聚类的边框
  border_gp = gpar(lwd=1.5),
  split = c(rep(1,sum(deg_fin_mat_filter$cluster==unique(deg_fin_mat_filter$cluster)[1])),
            rep(2,sum(deg_fin_mat_filter$cluster==unique(deg_fin_mat_filter$cluster)[2]))),
  name = "z-score", 
  #heatmap_legend_param = list(at = bk_legend),
  #column_title = paste("Total",length(cluster_mat_scale[,1]),"genes"),
  column_split = sub_monocle$Lineage_10v12_branch,
  column_title_gp = gpar(fontsize=10),
  column_names_gp = gpar(fontsize=10),
  #row_title_gp = gpar(fontsize = anno_font_size),
  row_gap = unit(1.5, "mm"),
  #right_annotation = ha,
  show_row_names = F,
  show_column_names = F,
  cluster_rows = F,use_raster = F,
  cluster_columns = F,height = 10,width = 10,
  top_annotation =top_annotation,
  #left_annotation=left_annotation,
  col=col_fun)
p





DEG_in_monocle=function(){}


k5_bind_ac_tf_regulon=read.csv("k5_bind_ac_tf_regulon_fin.csv")

deg_fin_mat_filter=deg_fin_mat[deg_fin_mat$p_val_adj<0.01,]
deg_fin_mat_filter=deg_fin_mat_filter[deg_fin_mat_filter$gene%in%k5_bind_ac_tf_regulon$x,]
deg_fin_mat_filter=deg_fin_mat_filter[deg_fin_mat_filter$lineage_branch%in%c("Lineage_14v12_branch"),]
table(deg_fin_mat_filter$cluster)


TSC_unique=unique(deg_fin_mat_filter[deg_fin_mat_filter$cluster=="TSC_CellFate","gene"])
TSC2_unique=unique(deg_fin_mat_filter[deg_fin_mat_filter$cluster=="TSC2_CellFate","gene"])
ESC_unique=unique(deg_fin_mat_filter[deg_fin_mat_filter$cluster=="ESC_CellFate","gene"])
PGC_unique=unique(deg_fin_mat_filter[deg_fin_mat_filter$cluster=="PGC_CellFate","gene"])
length(TSC_unique)
length(TSC2_unique)
length(ESC_unique)
length(PGC_unique)

cds_10_12=readRDS("cds_10_12_.rds")
cds_16_12=readRDS("cds_16_12_.rds")
cds_16_10=readRDS("cds_16_10_.rds")

cds_14_10=readRDS("cds_14_10_.rds")
cds_14_12=readRDS("cds_14_12_.rds")
cds_14_16=readRDS("cds_16_14_.rds")


cds=cds_14_12


pdata=as.data.frame(t(cds@reducedDimS))
# pdata$Prdm1=cds$Prdm1
# 
# color = c(colorRampPalette(colors = rev(brewer.pal(name = "RdYlBu"))[1:5])((50)),
#           colorRampPalette(colors = rev(brewer.pal(n = 10, name = "RdYlBu"))[6:10])((50)))
# ggplot(pdata,aes(x=V1,y=V2,alpha=Prdm1,color=Prdm1))+
#   geom_point()+
#   scale_color_gradientn(colors=palette_list$RdYlBu,limits=c(0,1.5))+
#   scale_alpha_continuous(limits=c(0,1.5))+
#   theme_classic()+
#   ylab("Component2")+
#   xlab("Component1")

monocle$Pseudotime_14v12=as.numeric(monocle$Pseudotime_14v12)
monocle[["Lineage14v12"]]=CreateDimReducObject(embeddings = as.matrix(pdata[,1:2]),
                                        key = "Component_",assay = "RNA")

FeatureDimPlot(monocle,reduction = "Lineage14v12",features = c("Pseudotime_14v12",TSC2_unique)
               ,pt.size = .5,ncol = 4)

  scale_color_gradientn(colors=palette_list$RdYlBu)


 
###Celltype annotation####

data("ref_scMCA")

ETP_ref=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/ETP_log2_TPM.rds")
ETP_ref[ETP_ref<0]=0
ETP_ref=as(as.matrix(ETP_ref),"dgCMatrix")
# 
# Bulk_in_vivo_ref=readRDS("Bulk_in_vivo.rds")
# Bulk_in_vivo_ref[Bulk_in_vivo_ref<0]=0
# Bulk_in_vivo_ref=as(as.matrix(Bulk_in_vivo_ref),"dgCMatrix")
# 
# icm_te=readRDS("in_vivo_data/E3.5_6.5_ICM_TE.Rds")
# icm_te=NormalizeData(icm_te)
# pgc=readRDS("in_vivo_data/PGC_E6.5-14.5.Rds")
# min.cells=dim(srt_query@assays$RNA$counts)[2]*0.01
# srt_query <- RunKNNPredict(srt_query = srt_query,
#                              #srt_ref = icm_te,ref_group = "Celltype",
#                              bulk_ref = Bulk_in_vivo_ref[,c(1,7:9)],
#                              filter_lowfreq = min.cells)
# CellDimPlot(srt_query,reduction = "tsne", group.by =c("KNNPredict_classification"),raster = F,
#             label = F,label.fg = "white",label_insitu = T)



min.cells=dim(monocle@assays$RNA$counts)[2]*0.01

monocle <- RunKNNPredict(srt_query = monocle, bulk_ref = ref_scMCA,filter_lowfreq = min.cells)

predictmeta=monocle$KNNPredict_classification

monocle <- RunKNNPredict(srt_query = monocle, bulk_ref = ETP_ref,filter_lowfreq = min.cells)

predictmeta2=monocle$KNNPredict_classification


predictmeta_modify=data.frame(predictmeta=strsplit2(strsplit2(predictmeta,"\\(")[,1],"\\_")[,1],predictmeta2=predictmeta2,predictmeta_modify="Others")
table(predictmeta_modify$predictmeta)
predictmeta_modify$predictmeta=gsub("ES","Others",predictmeta_modify[,1])
predictmeta_modify$predictmeta=gsub("TS","Others",predictmeta_modify[,1])
table(predictmeta_modify$predictmeta)
table(predictmeta_modify$predictmeta2)


for(i in 1:dim(predictmeta_modify)[1]){
  if(predictmeta_modify[i,1]=="Others"){
    predictmeta_modify[i,1]=predictmeta_modify[i,2]
  }
  if(predictmeta_modify[i,1]=="unreliable" & predictmeta_modify[i,2]!="EpiSCs"){
    predictmeta_modify[i,1]=predictmeta_modify[i,2]
  }
}
table(predictmeta_modify[,1])


monocle$annotated_celltype=predictmeta_modify[,1]
monocle$annotated_celltype=gsub("unreliable","Others",monocle$annotated_celltype)
monocle$annotated_celltype=factor(monocle$annotated_celltype,levels = c("ESC","TSC","PGCLC","EpiSCs",
                                                                        setdiff(unique(monocle$annotated_celltype),c("ESC","TSC","PGCLC","EpiSCs","Others")),
                                                                        "Others"))

DimPlot(monocle,reduction = "umap", group.by = "annotated_celltype", label = TRUE,pt.size = .1,
        cols = colorRampPalette(colors =as.character(SCP::palette_list$Set2#,as.character(SCP::palette_list$Set2),
                                          ))(12))
CellDimPlot(monocle,reduction = "umap", group.by =c("annotated_celltype","sample"),
            label = F,pt.size = .1,palcolor =c("#66C2A5" ,"#FB6A4A","#87CEFA", "#CCCCCC" ,
            "#DF92B6", "#8DA0CB" ,"#FCCDE5" ,"#F3CF5B" ,"#BEBADA", "#FFFFB3", "#DA9870","#666666"))


meta=read.csv("meta.csv",header = T,row.names = 1)


monocle$invivo_annoted_NR=factor(monocle$invivo_annoted_NR,levels = c("Astrocytes","Neuron","Radial glia cell",
                                                                      "Glomerular epithelial cell","Vascular smooth muscle cell","Stromal cell","Others","unreliable","None"))
DimPlot(monocle,reduction = "umap", group.by = "invivo_annoted_NR", label = TRUE,pt.size = .1,
        cols = c("#8DD3C7","#BC80BD","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#666666","#D9D9D9"))



ggsave(
  filename = "K5_D0-5_umap_ref_scMCA_annoted.pdf",
  width = 16,             # 宽
  height = 8,            # 高
  units = "in",          # 单位
  dpi = 1000              # 分辨率DPI
)




seurat_D0_5

monocle=readRDS("monocle2.rds")




seurat_D0_5$anno_celltype="Others"

meta=seurat_D0_5@meta.data

meta_anno=anno_seurat@meta.data
meta_anno$id=rownames(meta_anno)
meta_anno=meta_anno[,c("id","anno_celltype")]
meta_anno$anno_celltype=as.character(meta_anno$anno_celltype)
head(meta_anno)


meta$anno_celltype="Others"
for(i in 1:dim(meta_anno)[1]){
  meta[rownames(meta)==meta_anno[i,"id"],"anno_celltype"]=meta_anno[i,"anno_celltype"]
}

table(meta$anno_celltype)
monocle$anno_celltype=meta$anno_celltype
monocle$anno_celltype=factor(monocle$anno_celltype,levels=c("ESC","PGCLC","TSC","Others"))

CellDimPlot(seurat_D0_5,reduction = "tsne", group.by =c("anno_celltype2"),palette="Set2",
            label = F,label_insitu = T, label.size = 2,pt.size = .1#,palcolor =  c("#66C2A5" ,"#6699CC","#CC6677","grey")
            )

table(anno_seurat$seurat_clusters)

k5_0_5_etp_cl=anno_seurat

sub_k5_0_5_etp_cl=k5_0_5_etp_cl[,!k5_0_5_etp_cl$seurat_clusters%in%c(23,4,5,18,21,6,13,20)]
table(sub_k5_0_5_etp_cl$seurat_clusters)
sub_k5_0_5_etp_cl$Celltype=paste(sub_k5_0_5_etp_cl$Celltype,sub_k5_0_5_etp_cl$seurat_clusters,sep="_")
table(sub_k5_0_5_etp_cl$Celltype)
sub_k5_0_5_etp_cl=sub_k5_0_5_etp_cl[,!sub_k5_0_5_etp_cl$Celltype%in%c("ESC_0","ESC_1","ESC_17")]
table(sub_k5_0_5_etp_cl$Celltype)


new_sub=sub_k5_0_5_etp_cl[,grepl("TSC_11",sub_k5_0_5_etp_cl$Celltype)]
table(srt_query$KNNPredict_classification)
meta_anno=srt_query@meta.data
meta_anno$id=rownames(meta_anno)
meta_anno=meta_anno[,c("id","KNNPredict_classification")]
meta_anno$KNNPredict_classification=as.character(meta_anno$KNNPredict_classification)
table(meta_anno$KNNPredict_classification)

seurat_D0_5$anno_celltype2="Others"
for(i in 1:dim(meta_anno)[1]){
  seurat_D0_5@meta.data[colnames(seurat_D0_5)==meta_anno[i,"id"],"anno_celltype2"]=meta_anno[i,"KNNPredict_classification"]
}
table(seurat_D0_5$anno_celltype2)


#PAGA####
pancreas_sub <- RunPAGA(
  srt = pancreas_sub, group_by = "seurat_clusters",
  linear_reduction = "PCA", nonlinear_reduction = "UMAP"
)
PAGAPlot(srt = pancreas_sub, reduction = "UMAP", label = TRUE,
         label_insitu = TRUE, label_repel = TRUE)

####Cell fate branch analysis#####


PAGA_raw=readRDS("K5_PAGA_leiden_add_EpiSC_EpiLC.rds")
table(PAGA$sample)
table(PAGA$leiden)
PAGA=PAGA_raw[,!PAGA_raw$sample%in%c("EpiLC","EpiSC")]

episc=seurat_D0_5[,seurat_D0_5$seurat_cluster%in%c(13,4)]
episc=episc[,episc$sample=="EpiSCs"]
table(episc$sample)
episc$leiden=episc$seurat_cluster
table(episc$leiden)
episc$leiden=gsub(4,42,episc$leiden)
episc$leiden=gsub(13,43,episc$leiden)


PAGA=merge(PAGA,episc)
PAGA=JoinLayers(PAGA)
PAGA=NormalizeData(PAGA)


monocle_repgrm=PAGA[,PAGA$leiden%in%c(0,40,5,15,9,16,24,41,42,43,21,28)]

expr=as(as.matrix(monocle_repgrm@assays$RNA$counts),"sparseMatrix")
cellinfo=monocle_repgrm@meta.data
f_data=data.frame(gene_short_name=rownames(monocle_repgrm),row.names = rownames(monocle_repgrm))

cellinfo=new("AnnotatedDataFrame",data=cellinfo)
f_data=new("AnnotatedDataFrame",data=f_data)
cds=newCellDataSet(expr,
                   phenoData = cellinfo,
                   featureData=f_data,
                   lowerDetectionLimit=0.5,
                   expressionFamily=negbinomial.size()#for raw counts
)


cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)


cds=detectGenes(cds,min_expr=0.1)
#print(head(fData(cds)))
expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))
express_genes=VariableFeatures(FindVariableFeatures(monocle_repgrm,nfeatures = 2000))
express_genes=intersect(expressed_genes,express_genes)


cds=setOrderingFilter(cds,express_genes)
print(plot_ordering_genes(cds))

cds = reduceDimension(cds,
                      max_components=5,
                      method = "DDRTree")

cds=orderCells(cds)


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$leiden)[,"42"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))

cds@phenoData@data$Nanog=monocle_repgrm@assays$RNA$data["Nanog",]


cds$leiden=as.character(cds$leiden)
plot_cell_trajectory(cds,color_by="leiden",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Nanog",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="sample",size=1,show_backbone=T)










####D1 analysis#####



submonocle=readRDS("D1_PAGA.rds")
submonocle=monocle[,colnames(monocle)%in%c(rownames(monocle@meta.data[monocle@meta.data$sample=="EpiSCs",]),
                                           rownames(seurat_D0_5@meta.data[seurat_D0_5@meta.data$leiden%in%c(0,40,5,15,24,41,16),]))]


submonocle@meta.data=Btable_insert_into_Atable(submonocle@meta.data,seurat_D0_5@meta.data[,c(4,118)])
submonocle@meta.data=Btable_insert_into_Atable(submonocle@meta.data,meta[,c(4,27)])
submonocle$sample_cluster=paste(submonocle$sample,submonocle$seurat_cluster,sep="_")


submonocle$leiden[submonocle$sample=="EpiSCs"]="EpiSCs"

ssubmonocle=submonocle[,submonocle$leiden%in%c(0,15,16,24,40,41,5)]
ssubmonocle$leiden=factor(ssubmonocle$leiden,levels=c(40,0,5,15,24,16,41))



Stacked_VlnPlot(ssubmonocle, c("Pou5f1","Sall4","Nanog","Sox2",
                               #"Klf2","Zfp42",
                               "Prdm1","Ifitm1",
                               "Tead4","Eomes","Tfap2c","Gata3","Klf5"),#y.max = 2,
                x_lab_rotate = T,
                colors_use=c("#90EE90","#66C2A5","#8DD3C7", "#CCEBC5",
                             "#FDB462","#FFD92F", "#FFED6F"),
                group.by = "leiden")
"#AAAAAA"


sssubmonocle=submonocle[,submonocle$leiden%in%c(40,24,"EpiSCs")]
Idents(submonocle)="leiden"
deg_res_sss=FindMarkers(submonocle,ident.1 = c(15),ident.2 = "EpiSCs",
                      logfc.threshold = 1,only.pos = T,
                      min.pct = 0.25)

intersect(rownames(deg_res_sss),esc_relative)


####DEG from scanpy leiden clusters####
monocle=readRDS("monocle_for_PAGA.rds")
monocle=monocle[,monocle$sample!="EpiSCs"]
table(monocle$sample)

monocle@tools[1:length(monocle@tools)]=NULL

monocle$Klf5.=monocle@assays$RNA$scale.data[rownames(monocle@assays$RNA$scale.data)=="Klf5",]
monocle$Klf5.[monocle$Klf5.>1.5]=1.5


write.csv(monocle@meta.data,"meta.csv")

table(monocle$sample)
table(monocle$leiden)





monocle=readRDS("monocle_for_PAGA.rds")
monocle=CreateSeuratObject(monocle@assays$RNA$counts,meta.data=monocle@meta.data)

EpiLC=as.data.frame(fread("/home/data/t020406/8CLC_DDi/scRNAseq/GSE155088_UMI_counts_Bleckwehl_2020.csv"))
colnames(EpiLC)=c("gene",colnames(EpiLC)[2:dim(EpiLC)[2]])

colna=gsub("\"\"","",colnames(EpiLC)[1])
colna=as.character(strsplit2(colna,","))


anno=bitr(EpiLC$gene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Mm.eg.db")
EpiLC_anno=Btable_insert_into_Atable(anno,EpiLC,index_A =1,index_B=1 )

EpiLC_anno=readRDS("PGCLC_diff_scRNAseq.rds")
EpiLC_anno=EpiLC_anno[!duplicated(EpiLC_anno$SYMBOL),]
rownames(EpiLC_anno)=EpiLC_anno$SYMBOL
EpiLC_anno$SYMBOL=NULL
EpiLC_anno$gene=NULL
EpiLC_anno$ENSEMBL=NULL
EpiLC_anno=na.omit(EpiLC_anno)

barcode=read.csv("t-SNE-Projection.csv")
barcode=barcode$Barcode

colnames(EpiLC_anno)=barcode


EpiLC=CreateSeuratObject(EpiLC_anno)
EpiLC$sample="EpiLC"
EpiLC$leiden="EpiLC"

epilc_d2=read.csv("LibraryID.csv")

EpiLC_d2=EpiLC[,colnames(EpiLC)%in%epilc_d2$Barcode]

episc=read.csv("LibraryID_EpiSC.csv")
episc=EpiLC[,colnames(EpiLC)%in%episc$Barcode]
episc$sample="EpiSC"
episc$leiden="EpiSC"

obj=EpiLC
Seurat_stander_process{
  
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
  obj=subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 10000)
  
  
  
  obj <- NormalizeData(obj)
  obj <- ScaleData(obj)
  
  obj=FindVariableFeatures(obj,nfeatures = 2000)
  
  
  #obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  
  
  
  obj <- RunPCA(obj)
  
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
  obj <- FindClusters(obj, resolution = 1,random.seed =1, cluster.name = c("seurat_clusters"))
  #obj <- FindClusters(obj, resolution = 0.65,random.seed =1, cluster.name = c("seurat_clusters"))
  
  
  obj <- RunUMAP(obj, reduction = "pca",
                 dims = 1:30, seed.use =2,
                 n.neighbors = 20, #20更多粘连和小块,50偏整体http://biotrainee.cn:9013/graphics/plot_zoom_png?width=912&height=990
                 #n.epochs = 200,#200大数据，500小数据
                 min.dist = 1,#距离
                 #spread=1,#揉
                 reduction.name = "umap",verbose = FALSE)
  
  
  p1=DimPlot(obj,reduction = "umap",group.by = c("seurat_clusters"))
  obj$seurat_clusters=paste(obj$seurat_clusters,obj$sample)
  
  p2=Stacked_VlnPlot(obj,zto4_marker2,x_lab_rotate = T,#y.max = 2,
                     # colors_use=c("#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA","#CCCCCC"),
                     #colors_use=c("#FFD700","#DE2D26","#FF8C00","#90EE90","#87CEFA"),
                     group.by = "seurat_clusters")
  p1+p2
  
  DoHeatmap(obj,features = zto4_gene,group.by = "seurat_clusters",disp.min = -1,disp.max = 1)
  table(obj$seurat_clusters)
}






epi=readRDS("/home/data/t020406/8CLC_DDi/scRNAseq/EpiSCs_QC.rds")
epi$leiden="EpiSCs"

monocle=merge(monocle,c(EpiLC_d2,episc))
monocle=JoinLayers(monocle)
monocle=NormalizeData(monocle)
monocle=ScaleData(monocle)


monocle=AddModuleScore(monocle,features = formative,
                          name = "formative")

VlnPlot(monocle, features = c("formative1"),group.by = "leiden",pt.size = 0)
VlnPlot(monocle, features = c("Zyx"),group.by = "leiden",pt.size = 0)


monocle_without=monocle[,monocle$sample!="EpiSCs"]
write.csv(monocle_without@meta.data,"PAGA_leiden_meta_ver4.csv")


#####Non-Reprogramming part####

meta=read.csv("PAGA_leiden_meta.csv")
table(meta$leiden)


monocle$leiden=meta$leiden
monocle$annotated_celltype=meta$annotated_celltype

meta=monocle@meta.data
meta$PAGA_big_group="Others"
meta[meta$leiden%in%c(35,23,32,34,37),"PAGA_big_group"]="Reprogramming_branch"
meta[meta$leiden%in%c(36,30,26,27,22,17,33,1,2,3,
                      18,10,8,14,12,31,29),"PAGA_big_group"]="Non-Reprogramming_branch"
meta[meta$leiden%in%c(20,25,39),"PAGA_big_group"]="Trophoblaste_branch"
meta[meta$leiden%in%c(0,40,5,15,9,16,24,41),"PAGA_big_group"]="Precursor"
monocle@meta.data=meta

min.cells=dim(monocle@assays$RNA$counts)[2]*0.01

monocle <- RunKNNPredict(srt_query = monocle, bulk_ref = ref_scMCA,filter_lowfreq = min.cells)


meta$KNNPredict_classification=strsplit2(strsplit2(meta$KNNPredict_classification,"\\(")[,1],"\\_")[,1]

meta$invivo_annoted_NR="None"
meta[meta$PAGA_big_group=="Non-Reprogramming_branch","invivo_annoted_NR"]=meta[meta$PAGA_big_group=="Non-Reprogramming_branch","KNNPredict_classification"]

meta[meta$invivo_annoted_NR%in%c("Epithelial cell","ES","TS"),"invivo_annoted_NR"]="Others"

meta[meta$invivo_annoted_NR=="MEF","invivo_annoted_NR"]="Vascular smooth muscle cell"

table(meta$invivo_annoted_NR)

monocle@meta.data=meta

write.csv(meta,"meta.csv")

use_monocle=monocle[,monocle$PAGA_big_group!="Others"]
table(use_monocle$PAGA_big_group)

Idents(use_monocle)="PAGA_big_group"

use_monocle=NormalizeData(use_monocle)
all.genes <- rownames(use_monocle)
use_monocle <- ScaleData(use_monocle, features = all.genes)


deg_res=FindAllMarkers(use_monocle,only.pos = T,
               logfc.threshold = 1,
               min.pct = 0.1)


saveRDS(deg_res,"PAGA_big_group_DEG_resluts_minpct0.1.rds")


deg_res_filter=deg_res[deg_res$p_val_adj<0.01,]
deg_res_filter=deg_res_filter[deg_res_filter$avg_log2FC>2,]

table(deg_res_filter$cluster)


#####leiden part####
exprMat <- t(as.matrix(use_monocle@assays$RNA@data))
write.csv(use_monocle@meta.data, "leiden_c40,21,24,28_cellInfo.csv", row.names = T)
gc()
write.csv(exprMat, "leiden_c40,21,24,28_exprMat.csv")

exprMat <- t(as.matrix(use_monocle2@assays$RNA@data))
write.csv(use_monocle2@meta.data, "leiden_c6,11,4,7,19_cellInfo.csv", row.names = T)
gc()
write.csv(exprMat, "leiden_c6,11,4,7,19_exprMat.csv")


use_monocle=monocle[,monocle$leiden%in%c(24,28,40,21)]
table(use_monocle$leiden)


use_monocle2=monocle[,monocle$leiden%in%c(4,7)]
table(use_monocle2$leiden)


use_monocle3=monocle[,monocle$leiden%in%c(35,37)]
table(use_monocle3$leiden)


Idents(use_monocle)="leiden"
Idents(use_monocle2)="leiden"
Idents(use_monocle3)="leiden"


deg_res24=FindMarkers(use_monocle,ident.1 = 24,ident.2 = 40,
                       logfc.threshold = 1,only.pos = T,
                       min.pct = 0.25)
deg_res40=FindMarkers(use_monocle,ident.1 = 40,ident.2 = 24,
                     logfc.threshold = 1,only.pos = T,
                     min.pct = 0.25)

deg_res28=FindMarkers(use_monocle,ident.1 = 28,ident.2 = 21,
                      logfc.threshold = 1,only.pos = T,
                      min.pct = 0.25)
deg_res21=FindMarkers(use_monocle,ident.1 = 21,ident.2 = 28,
                      logfc.threshold = 1,only.pos = T,
                      min.pct = 0.25)


deg_res=FindAllMarkers(use_monocle2,only.pos = T,
                       logfc.threshold = 1,
                       min.pct = 0.25)

deg_res4=FindMarkers(use_monocle2,ident.1 = 4,ident.2 = 7,
                      logfc.threshold = 1,only.pos = T,
                      min.pct = 0.25)
deg_res7=FindMarkers(use_monocle2,ident.1 = 7,ident.2 = 4,
                      logfc.threshold = 1,only.pos = T,
                      min.pct = 0.25)


deg_res35=FindMarkers(use_monocle3,ident.1 = 35,ident.2 = 37,
                     logfc.threshold = 1,only.pos = T,
                     min.pct = 0.25)
deg_res37=FindMarkers(use_monocle3,ident.1 = 37,ident.2 = 35,
                     logfc.threshold = 1,only.pos = T,
                     min.pct = 0.25)



test=c("Sox15","Dnmt3l","Tex19.1","Eomes","Apela","St14","Esrrb","Nrk","Ascl2","Llgl2","Cdh1","Gcm1","Lif","Tfeb","
  Epas1","Plcd3","Etnk2","Ovol2","Plcd1","Gjb3","Gjb5","Tppp3","Htra1","Wnt7b","Spint2","Fzd5","Grhl2","
  Foxd3","Ncoa3","Nodal","Peg10","Phlda2","Ncoa1","Birc2")

saveRDS(deg_res,"leiden_24,28,40,21_DEG_resluts_minpct0.1.rds")


filter_gene=
  function(deg_res){
    deg_res_filter=deg_res[deg_res$p_val_adj<0.01,]
    deg_res_filter=deg_res_filter[deg_res_filter$avg_log2FC>1,]
    deg_res_filter=deg_res_filter[order(deg_res_filter$avg_log2FC,decreasing = T),]
    return(deg_res_filter)
  }


deg_res_filter24=filter_gene(deg_res24)
deg_res_filter40=filter_gene(deg_res40)
deg_res_filter28=filter_gene(deg_res28)
deg_res_filter21=filter_gene(deg_res21)


deg_res_filter=filter_gene(deg_res)
deg_res_filter4=deg_res_filter[deg_res_filter$cluster==4,"gene"]
deg_res_filter7=deg_res_filter[deg_res_filter$cluster==7,"gene"]
deg_res_filter19=deg_res_filter[deg_res_filter$cluster==19,"gene"]

deg_res_filter4=filter_gene(deg_res4)
deg_res_filter7=filter_gene(deg_res7)


deg_res_filter35=filter_gene(deg_res35)
deg_res_filter37=filter_gene(deg_res37)


deg_res_list1=list(C24=rownames(deg_res_filter24),
                  C40=rownames(deg_res_filter40),
                  C28=rownames(deg_res_filter28),
                  C21=rownames(deg_res_filter21))

deg_res_list2=list(C4=rownames(deg_res_filter4),
                   C7=rownames(deg_res_filter7)#,
                   #C19=deg_res_filter19
                   )

deg_res_list3=list(C37=rownames(deg_res_filter37),
                   C35=rownames(deg_res_filter35))


saveRDS(deg_res_list1,"deg_res_list_40,21,24,28.rds")
saveRDS(deg_res_list2,"deg_res_list_4,7.rds")
saveRDS(deg_res_list3,"deg_res_list_35,37.rds")


VennPlot(deg_res_list1)


write.table(rownames(deg_res_filter40),"deg_res_filter40.txt")
write.csv(deg_res_filter21,"deg_res_filter21.csv")
write.csv(deg_res_filter24,"deg_res_filter24.csv")
write.csv(deg_res_filter28,"deg_res_filter28.csv")



deg_res_top50_list=list(deg_res_filter24=rownames(deg_res_filter24)[1:50],
                  deg_res_filter40=rownames(deg_res_filter40)[1:50],
                  deg_res_filter28=rownames(deg_res_filter28)[1:50],
                  deg_res_filter21=rownames(deg_res_filter21)[1:50])


deg_res_tf_list=list(deg_res_filter24=intersect(rownames(deg_res_filter24),mm_tf$V1),
                  deg_res_filter40=intersect(rownames(deg_res_filter40),mm_tf$V1),
                  deg_res_filter28=intersect(rownames(deg_res_filter28),mm_tf$V1),
                  deg_res_filter21=intersect(rownames(deg_res_filter21),mm_tf$V1))


#####GO####
Compare_Cluster_GO_results(deg_res_list1,top_go_term = 10,simplify_ver2 = T)
Compare_Cluster_GO_results(deg_res_list2,top_go_term = 10,simplify_ver2 = T)
Compare_Cluster_GO_results(deg_res_list3,top_go_term = 10,simplify_ver2 = T)



term=rev(c("DNA replication","histone modification","chromatin remodeling","stem cell population maintenance",
           "positive regulation of cell cycle",
           "regulation of cell cycle phase transition",
           "actin filament organization","regulation of neurogenesis","cell-substrate adhesion","muscle cell differentiation",
           "gliogenesis",
           "ossification","renal system development"))

Compare_Cluster_GO_results(deg_res_list1[c(2,4,1,3)],top_go_term = 5,pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = list.files(path = "PAGA_DEG_LFC_1_0.3_GO_results/",pattern = "BP_Simplify_go.csv",full.names = T)[c(6,2,3,4)],
                           custom_design_go_term = term)


term=rev(c("cellular response to leukemia inhibitory factor","placenta development","gland morphogenesis","blastocyst formation",
           'stem cell differentiation',"actin filament organization","response to oxidative stress",
           "pattern specification process","embryonic organ morphogenesis","regulation of neurogenesis",
           "glial cell differentiation","ossification","urogenital system development","muscle cell differentiation"))

term=rev(c("cellular response to leukemia inhibitory factor","placenta development","gland morphogenesis","blastocyst formation",
           'stem cell differentiation',"actin filament organization","response to oxidative stress",
           "pattern specification process","embryonic organ morphogenesis","regulation of neurogenesis",
           "glial cell differentiation","ossification","urogenital system development","muscle cell differentiation"))
#[c(7,5,1,8)]


term=rev(c("actin filament organization","cellular response to leukemia inhibitory factor",
           "gland morphogenesis","response to oxidative stress","epithelial cell development",
          "DNA methylation or demethylation","stem cell population maintenance",
          "pattern specification process","heart morphogenesis","urogenital system development",
          "regulation of neuron differentiation","muscle tissue development",
          "Wnt signaling pathway","response to BMP"))

Compare_Cluster_GO_results(deg_res_list2,top_go_term = 7,pvalue_max_limit = 25,
                           input_go_cluster_resluts_pathlist = list.files(path = "PAGA_DEG_LFC_1_0.25_GO_results/",pattern = "BP_Simplify_go.csv",full.names = T)[c(6,8)],
                           custom_design_go_term = term)




Compare_Cluster_GO_results(inter_deg_list,top_go_term = 10,simplify_ver2 = T)


term=rev(c("cell junction assembly","regulation of histone modification",
           "insulin receptor signaling pathway","positive regulation of cell cycle",
           "ribonucleoprotein complex biogenesis","DNA replication",
           "regulation of chromosome organization","telomere maintenance","stem cell population maintenance",
           "mitotic cell cycle phase transition","meiotic cell cycle",
           "axonogenesis","cell-substrate adhesion","wound healing","skin development",
           "muscle cell differentiation","regulation of neurogenesis","renal system development","ossification","heart process",
           "stem cell differentiation","regulation of Wnt signaling pathway","kidney epithelium development",
           "neuron fate commitment","regulation of apoptotic signaling pathway"))


Compare_Cluster_GO_results(inter_deg_list,top_go_term = 10,pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = list.files(pattern = "BP_Simplify_go.csv")[c(5,6,1,
                                                                                                                 3,4,2)],
                           custom_design_go_term = term)


xx=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C40 (449)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(4,19,23),],display_top = 3,save = F)

xx=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C24 (1667)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(1,4,8),],display_top = 3,save = F)

xx21=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C21 (1688)_BP_Simplify_go.csv")
Plot_GO_result(xx21[c(2,12,39),],display_top = 3,save = F)
xx28=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C28 (2219)_BP_Simplify_go.csv")
Plot_GO_result(xx28[c(2,10,15),],display_top = 3,save = F)


Compare_Cluster_GO_results(list(C21=rownames(deg_res_filter21),
                                C28=rownames(deg_res_filter28)),
                           top_go_term = 5,#pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = c("PAGA_DEG_LFC_1_0.3_GO_results/C21 (1688)_BP_Simplify_go.csv",
                                                                 "PAGA_DEG_LFC_1_0.3_GO_results/C28 (2219)_BP_Simplify_go.csv"),
                           custom_design_go_term = rev(c("DNA replication","regulation of cell cycle phase transition",
                                                         "chromatin remodeling","histone modification","stem cell population maintenance",
                                                     "actin filament organization","regulation of neurogenesis",
                                                     "cell-substrate adhesion","muscle cell differentiation","ossification")))
#3.3x4.7

xx4=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C4 (584)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(1,2,5),],display_top = 3,save = F)
xx7=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C7 (460)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(1,9,17,49),],display_top = 3,save = F)
xx19=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C19 (402)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(1,5,18),],display_top = 3,save = F)

Compare_Cluster_GO_results(deg_res_list2,
                           top_go_term = 5,#pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = c("PAGA_DEG_LFC_1_0.3_GO_results/C4 (584)_BP_Simplify_go.csv",
                                                                 "PAGA_DEG_LFC_1_0.3_GO_results/C7 (460)_BP_Simplify_go.csv",
                                                                 "PAGA_DEG_LFC_1_0.3_GO_results/C19 (402)_BP_Simplify_go.csv"),
                           custom_design_go_term = rev(c("cellular response to leukemia inhibitory factor",
                                                         "stem cell population maintenance","maintenance of cell number",
                                                         "stem cell differentiation","epigenetic regulation of gene expression",
                                                         "pattern specification process","mesenchyme development",
                                                         "kidney development","regulation of neurogenesis","heart morphogenesis",
                                                         "embryonic placenta development","epidermis development",
                                                         "cell-substrate adhesion","keratinocyte differentiation",
                                                         "ossification"
                                                         )))
#3.6x5.2

xx37=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C37 (805)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(1,15,18),],display_top = 3,save = F)
xx35=read.csv("PAGA_DEG_LFC_1_0.3_GO_results/C35 (565)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(1,6,97),],display_top = 3,save = F)


Compare_Cluster_GO_results(deg_res_list3,
                           top_go_term = 5,#pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = c("PAGA_DEG_LFC_1_0.3_GO_results/C37 (805)_BP_Simplify_go.csv",
                                                                 "PAGA_DEG_LFC_1_0.3_GO_results/C35 (565)_BP_Simplify_go.csv"),
                           custom_design_go_term = rev(c("response to leukemia inhibitory factor",
                                                         "cell-substrate adhesion","stem cell population maintenance",
                                                         "maintenance of cell number","response to oxidative stress",
                                                         "actin filament organization","mitochondrial transmembrane transport",
                                                         "cellular lipid catabolic process","histone modification",
                                                         "germ cell development")))
#3.2X4.7





#####DEG heatmap####

deg_res_list1=readRDS("deg_res_list_40,21,24,28.rds")
readRDS("deg_res_list_4,7.rds")
readRDS("deg_res_list_35,37.rds")


venn=get.venn.partitions(deg_res_list1)
names(venn$..values..)=venn$..count..
as.data.frame(names(venn$..values..))
names(venn$..values..[c(14,6,8,15,11,12)])

names(venn$..values..)=1:15

inter_deg_list=venn$..values..[c(14,6,8,15,11,12)]
names(inter_deg_list)=c("40","40&21","21","24","28&24","28")



sub_use_monocle=use_monocle3#[,use_monocle3$leiden%in%c(28,21)]
table(sub_use_monocle$leiden)


deg_mat=sub_use_monocle[rownames(sub_use_monocle)%in%setdiff(deg_res_list1$C28,deg_res_list1$C24),]
deg_mat$leiden=factor(deg_mat$leiden,levels=c("40","21","24","28"))

deg_mat=sub_use_monocle[rownames(sub_use_monocle)%in%deg_res_list2$C7,]
deg_mat$leiden=factor(deg_mat$leiden,levels=c("4","7"))

deg_mat=sub_use_monocle[rownames(sub_use_monocle)%in%deg_res_list3$C35,]
deg_mat$leiden=factor(deg_mat$leiden,levels=c("37","35"))



venn=get.venn.partitions(deg_res_list2)
names(venn$..values..)=venn$..count..
as.data.frame(names(venn$..values..))
names(venn$..values..[c(15,11,12,14,6,8)])

names(venn$..values..)=1:15

inter_deg_list=venn$..values..[c(15,11,12,14,6,8)]
names(inter_deg_list)=c("6","6&4","4","11","11&7","7")

deg_mat=use_monocle2[rownames(use_monocle2)%in%inter_deg_list$`7`,]
deg_mat=deg_mat[,!deg_mat$leiden==19]
deg_mat$leiden=factor(deg_mat$leiden,levels=c("6","4","11","7"))


deg_mat=use_monocle2[rownames(use_monocle2)%in%deg_res_list$C19,]
deg_mat$leiden=factor(deg_mat$leiden,levels=c("4","7","19"))



dim(deg_mat)
DoHeatmap(deg_mat,group.by = "leiden",features = rownames(deg_mat),
          disp.max = 2,disp.min = -2,lines.width =round(dim(deg_mat@assays$RNA$data)[2]*0.01))+
  scale_fill_gradientn(colors = c("#4B0082","black","yellow"))

c(palette_list$RdYlBu)



pluripotency=c("Pou5f1","Nanog","Zfp42","Esrrb","Sox2","Tfcp2l1",
               "Nr0b1","Tbx3","Fgf4","Sall4","Klf5","Nr5a2","Prdm14",
               "Zic2","Zic3","Zfp281","Lifr","Bmp4","Dppa3",
               "Pou5f1","Zfp42","Nanog","Dppa5a","Sall4","Tfcp2l1","Sox2","Dppa4","Dppa3",
               "Nr5a2","Tbx3","Nanog","Klf2","Klf4",
               "Esrrb","Nr0b1","Prdm14","Prdm1","Lin28a",
               "Myc","Utf1","Td1","Fgf4","Lin28a"#oct4,L80A文章
)
tsc=c("Cdx2","Tfap2c","Elf5","Gata3","Hand1","Eomes","Tead4","Fgfr2","Krt7","Krt8",
      "Ascl2","Plet1","Rin3","Tinagl1","Lasp1","Map3k8","Sh2d3c")
tsc_stem=c("Cdx2","Tfap2c","Elf5","Gata3","Eomes","Tead4","Fgfr2","Id1","Id2","Sox2",
           "Plet1","Rin3","Lasp1","Map3k8","Sh2d3c","Krt7",
           "Ly6a","Ddah1","Gsto1","Epcam","Smad6","Smad7","E2f4","Smad3" #epi inducer文章
)
tsc_diff=c("Prl3d1","Junb","Prl2c2","Dsp","Klf6","Cebpb","Krt19","Krt18","Jund","Cdkn1c","Tinagl1",
           "Cgn","Gcm1","Synb","Syna","Dlx3","Prl3b1","Prl8a9","Tpbpa","Ascl2","Hand1","Krt8")
test=c(study_gene,"Zic3","Sox2","Hdac2","Cebpz","Tead4","Id1","Id2","Id3","Id4",
       "Bmp1","Bmp5","Bmp4","Bmp6","Bmp7","Smad1","Smad3",'Smad5',"Smad6","Smad7","Ifitm1","Ifitm3","Ifitm2",
       signal[signal$signaling=="sig_BMP","gene"],tsc,tsc_stem,tsc_diff,
       pluripotency)
test=test[!duplicated(test)]
Study_Gene_Distribution_in_gene_list(test,
                                     deg_res_list2)


#####Klf5 expression and K5_ac_2_bind_gene relasionship####

DotPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,41)],features = c("k5_ac_2_bind_AUC"),group.by = "leiden")
DotPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,41)],features = c("Klf5"),group.by = "leiden")

VlnPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,41)],features = c("k5_ac_2_bind_AUC"),group.by = "leiden")
VlnPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,41)],features = c("Klf5"),group.by = "leiden")

saveRDS(seurat_D0_5,"monocle_for_PAGA.rds")

#DefaultAssay(sub1)="RNA"
DefaultAssay(seurat_D0_5)="RNA"
#DefaultAssay(immune.combined)="integrated"
expr=AverageExpression(seurat_D0_5,group.by = "leiden",assays = "RNA")
#vf=esc_relative
expr$RNA[rownames(expr$RNA)%in%"Klf5",]



#vf=c(feature1,feature2,feature3,feature4)
options(repr.plot.width=12,repr.plot.height=12)
cor_mat=cor(data.frame(ac=seurat_D0_5$k5_ac_2_bind_AUC,
            k5=seurat_D0_5@assays$RNA@data[rownames(seurat_D0_5@assays$RNA@data)=="Klf5",]),
            method="spearman")
pheatmap::pheatmap(cor_mat,display_numbers = T)



Caculate_expr_and_acgene_correlation=
  function(seurat_D0_5,cluster){
    
    cor_mat=data.frame(cluster=cluster,cor=rep(0,length(cluster)))
    
    for(n in 1:length(cluster)){
      c=cluster[n]
      subseurat=seurat_D0_5[,seurat_D0_5$leiden==c]
      subcor=cor(data.frame(ac=subseurat$k5_ac_tsc_2_bind_AUC,
                            k5=subseurat@assays$RNA@data[rownames(subseurat@assays$RNA@data)=="Klf5",]))
      cor_mat[n,2]=subcor[2,1]
    }
    return(cor_mat)
  }


cor_mat=Caculate_expr_and_acgene_correlation(seurat_D0_5,c(0,5,9,15,16,24,40,41))

ggplot(cor_mat,aes(x=factor(cluster,levels=cluster),y=cor))+
  geom_bar(stat = "identity")


c=cluster[8]
subseurat=seurat_D0_5[,seurat_D0_5$leiden==c]
pmat=data.frame(ac=subseurat$k5_ac_2_bind_AUC,
           k5=subseurat@assays$RNA@data[rownames(subseurat@assays$RNA@data)=="Klf5",])
ggplot(pmat,aes(x=ac,y=k5))+
  geom_pointdensity()+
  scale_color_gradientn(colours=palette_list$Spectral)+
  ggtitle(c)+
  scale_x_continuous(limits = c(0,0.1))+
scale_y_continuous(limits = c(0,4))+
  theme_bw()



####DEG heatmap####

deg_res=readRDS("PAGA_big_group_DEG_resluts_minpct0.05.rds")
deg_res_filter=deg_res[deg_res$p_val_adj<0.01,]
deg_res_filter=deg_res_filter[deg_res_filter$avg_log2FC>1,]


use_monocle=monocle[,monocle$PAGA_big_group!="Others"]
table(use_monocle$PAGA_big_group)

use_monocle=use_monocle[rownames(use_monocle)%in%deg_res_filter$gene,]

use_monocle=NormalizeData(use_monocle)
all.genes <- rownames(use_monocle)
use_monocle <- ScaleData(use_monocle, features = all.genes)


deg_mat=as.data.frame(use_monocle@assays$RNA@counts)
deg_mat=as.data.frame(t(deg_mat))
deg_mat$PAGA_big_group=use_monocle$PAGA_big_group
deg_mat$PAGA_big_group=factor(deg_mat$PAGA_big_group,levels = c("Reprogramming_branch",
                                                  "Non-Reprogramming_branch",
                                                  "Trophoblaste_branch"))
deg_mat=deg_mat[order(deg_mat$PAGA_big_group),]
deg_mat$PAGA_big_group=NULL


deg_mat=as.data.frame(t(deg_mat))


fin_deg_mat=deg_mat[1,]
for(cl in c("Reprogramming_branch","Non-Reprogramming_branch","Trophoblaste_branch")){
  sub_gene=deg_res_filter[deg_res_filter$cluster==cl,"gene"]
  sbu_mat=deg_mat[rownames(deg_mat)%in%sub_gene,]
  fin_deg_mat=rbind(fin_deg_mat,sbu_mat)
}
fin_deg_mat=fin_deg_mat[-1,]



bk=seq(-1,2,0.1)
col_fun = colorRamp2(bk,
                     c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"))(length(bk)/2-1),
                       colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))(length(bk)/2)),space = "RGB")

col_fun = colorRamp2(bk,c(colorRampPalette(colors = c(palette_list$RdGy))(length(bk))),space = "RGB")




top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill=brewer.pal(n = 3, name = 'Set2')),
                                                   labels = c("Reprogramming_branch","Non-Reprogramming_branch","Trophoblaste_branch"), 
                                                   labels_gp = gpar(col = "black", #fontface="bold",
                                                                    fontsize = 8)))


fin_deg_seurat=CreateSeuratObject(fin_deg_mat)
fin_deg_seurat$PAGA_big_group=factor(c(rep("Reprogramming_branch",3122),
                                rep("Non-Reprogramming_branch",21807),
                                rep("Trophoblaste_branch",2317)),levels = c("Reprogramming_branch","Non-Reprogramming_branch","Trophoblaste_branch"))

saveRDS(fin_deg_seurat,"fin_deg_seurat.rds")



ht_opt$raster_temp_image_max_width=10
ht_opt$raster_temp_image_max_height=10

deg_mat_scale=t(scale(t(deg_mat)))
p=Heatmap( #cluster_mat_scale[,c(1:2,27:34,3:26,35:67)], 
  fin_deg_mat,
  border = T,#是否开启聚类的边框
  border_gp = gpar(lwd=1.5),
  split = c(rep(1,2664),rep(2,1911),rep(3,3119)),
  name = "z-score",
  #heatmap_legend_param = list(at = bk_legend),
  #column_title = paste("Total",length(cluster_mat_scale[,1]),"genes"),
  column_split = c(rep(1,3122),rep(2,21807),rep(3,2317)),
  column_title_gp = gpar(fontsize=10),
  column_names_gp = gpar(fontsize=10),
  #row_title_gp = gpar(fontsize = anno_font_size),
  row_gap = unit(1.5, "mm"),
  #right_annotation = ha,
  show_row_names = F,
  show_column_names = F,
  cluster_rows = F,use_raster = T,raster_by_magick = T,
  cluster_columns = F,height = 10,width = 10,
  top_annotation =top_annotation,
  #left_annotation=left_annotation,
  col=col_fun)
p








####GO####

deg_list=Convert_Long_dataframe_to_list(deg_res_filter,group_index = 6,variable_index = 7)

Compare_Cluster_GO_results(deg_list,simplify_ver2 = T,top_go_term = 10)

term=rev(c("cellular response to leukemia inhibitory factor",
           "germ cell development",
           "DNA methylation or demethylation",
           "stem cell population maintenance",
           "histone modification",
           "regulation of neuron differentiation","muscle tissue development",
           "artery development","canonical Wnt signaling pathway","cell junction assembly",
           "regulation of epithelial cell proliferation","response to BMP",
           "epidermis development",
           "placenta development",
           "regulation of vasculature development",
           "heart process",
           "cell-cell junction organization",
           "muscle cell differentiation"))



Compare_Cluster_GO_results(deg_list[c(3,1,2)],custom_design_go_term = term,pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = list.files(pattern = "BP_Simplify_ver2_go.csv")[c(2,1,3)],
                           simplify_ver2 = T,top_go_term = 5)



####Marker of each cell type####

"Astrocytes","Neuron","Radial glia cell",
"Glomerular epithelial cell","MEF","Stromal cell"



deg_res_filter_NR=deg_res_filter[deg_res_filter$cluster=="Non-Reprogramming_branch",]


intersect(Stromal_cell,deg_res_filter_NR$gene)

intersect(Glomerular_epithelial_cell,rownames(monocle))



use_monocle=monocle[,monocle$PAGA_big_group=="Non-Reprogramming_branch"]
table(use_monocle$PAGA_big_group)


use_monocle=NormalizeData(use_monocle)
all.genes <- rownames(use_monocle)
use_monocle <- ScaleData(use_monocle, features = all.genes)


Idents(use_monocle)="invivo_annoted_NR"
deg_res_annoted=FindAllMarkers(use_monocle,only.pos = T,
                                       logfc.threshold = 1,
                                       min.pct = 0.01)

top10 <- deg_res_annoted %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(use_monocle, features = top10$gene) + NoLegend()



Astrocytes=c("Fabp7","Gfap","Sox9")
Neuron=c("Map2","Dcx","Snap25","Tubb3")
Radial_glial_cell=c("Nes","Notch2","Hes5","Fzd10")

Glomerular_epithelial_cell=c("Foxd1","Grrp1")
# "Col1a1","Col3a1","Dcn","Des","Fgf7",
# "Lef1","Serpine2","Gdnf",
# "Fabp2","Reg1",
VSMC=c("Acta2","Myl9","Tagln","Cnn1")#Vascular smooth muscle cell(VSMC)


Stromal_cell=c("Cd248","Fstl1","Col5a2","S100a4","Col1a1","Col1a2","Col3a1","Fap","Cxcl12",
               "Pdgfra","Pdgfrb","Prrx1")


all=c(Astrocytes,Neuron,Radial_glial_cell,Glomerular_epithelial_cell,VSMC,Stromal_cell)




deg_res_annoted_list=Convert_Long_dataframe_to_list(deg_res_annoted,group_index = 6,variable_index = 7)


Compare_Cluster_GO_results(deg_res_annoted_list,
                           simplify_ver2 = T,top_go_term = 10)



DotPlot(monocle,group.by = "invivo_annoted_NR",features =Glomerular_epithelial_cell)




fetal_marker=read.csv("Fetal kidney.csv")


intersect(deg_res_filter_NR$gene,
deg_res_annoted[deg_res_annoted$cluster=="Glomerular epithelial cell","gene"])




sub_monocle=monocle[,monocle$PAGA_big_group!="Others"]
table(sub_monocle$PAGA_big_group)


meta=sub_monocle@meta.data
meta[meta$invivo_annoted_NR=="Astrocytes","PAGA_big_group"]="Astrocytes"
meta[meta$invivo_annoted_NR=="Neuron","PAGA_big_group"]="Neuron"
meta[meta$invivo_annoted_NR=="Radial glia cell","PAGA_big_group"]="Radial glia cell"
meta[meta$invivo_annoted_NR=="Glomerular epithelial cell","PAGA_big_group"]="Glomerular epithelial cell"
meta[meta$invivo_annoted_NR=="Stromal cell","PAGA_big_group"]="Stromal cell"
sub_monocle@meta.data=meta

ssub_monocle=sub_monocle[,sub_monocle$PAGA_big_group!="Non-Reprogramming_branch"]
table(ssub_monocle$PAGA_big_group)
ssub_monocle$PAGA_big_group=factor(ssub_monocle$PAGA_big_group,
                                  levels =c("Precursor","Reprogramming_branch","Trophoblaste_branch",
                                             "Neuron","Astrocytes","Radial glia cell",
                                             "Glomerular epithelial cell","Stromal cell"))

use_ftr=c("Id2","Smad5","Bmpr2","Bmpr1a",
          "Map2","Tubb3","Fabp7","Sox9",
          "Hes5","Fzd10","Foxd1","Grrp1",
          "Col5a2","Pdgfrb","Fap")
# Stacked_VlnPlot(sub_monocle, use_ftr,x_lab_rotate = T,#y.max = 2,
#                 colors_use=c("#2178AE","#6EC5A3","#F74015","#D8BFD1"),
#                 group.by = "PAGA_big_group")

p=DotPlot(ssub_monocle,features =  use_ftr,
        group.by = "PAGA_big_group")
p+ggplot2::coord_flip()
p_data=p$data
p_anno=data.frame(id=levels(p_data$id),
                  label=factor(c(rep("D1",1),rep("Reprogramming",1),
                                 rep("TSCLC",1),rep("Non-Reprogramming",5)),
                               levels=c("D1","Reprogramming","TSCLC","Non-Reprogramming")))

p_data=plyr::join(p_data,p_anno)
# p_data$features.plot=factor(p_data$features.plot,levels =
#                               rev(c("Klf1","Klf15","Klf12","Klf14",
#                                     "Klf6","Klf17","Klf5","Klf4","Klf2","Klf7",
#                                     "Klf9","Klf3","Klf13","Klf8","Klf10","Klf11",
#                                     "Klf16")))
p_data$features.plot=factor(p_data$features.plot,levels =
                              rev(use_ftr))
p_data
scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size, #该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
limit_max=2
limit_max=2.5

name="Figure.1.A_Stringency_Single_cell_data_new_ver1"
pdf(paste(name,".pdf",sep=""),width=5.2,height =5)

color = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"   ))((50)),
          colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026" ))((50)))
#p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
#该函数可以完全替代scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  scale.func(range = c(0,6.5))+#点大小
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,#limits=c(0,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))




###RunMonocle3####


pancreas_sub <- RunMonocle3(srt = pancreas_sub, reduction = "umap",
                            root_cells=c(colnames(pancreas_sub[,pancreas_sub$seurat_clusters==5])),
                            clusters="seurat_clusters")
names(pancreas_sub@tools$Monocle3)
trajectory <- pancreas_sub@tools$Monocle3$trajectory
milestones <- pancreas_sub@tools$Monocle3$milestones

#CellDimPlot(pancreas_sub, group.by = "Monocle3_partitions", reduction = "umap", label = TRUE, theme_use = "theme_blank") + trajectory + milestones
CellDimPlot(pancreas_sub, group.by = "Monocle3_clusters", reduction = "umap", label = TRUE, theme_use = "theme_blank") + trajectory
FeatureDimPlot(pancreas_sub, features = "Monocle3_Pseudotime", reduction = "umap", theme_use = "theme_blank") + trajectory



###Visulization####
Lineage_UMAP=function(){}----------------------

want_lineage= c(2,4,5,10,12,16,14,9)
CellDimPlot(monocle, group.by = "sample", reduction = "umap",palcolor =sample_palette_col,
            lineages_palcolor = c("#8DA0CB", "#888888","#999933",
                                  "#E78AC3" ,"#66C2A5" ,"#6699CC",
                                  "#E5C494","#CC6677"),
            lineages = paste0("Lineage", want_lineage), lineages_span =0.3)

CellDimPlot(monocle, group.by = "seurat_clusters", reduction = "umap",
            lineages_palcolor =c("#8DA0CB", "#888888","#999933",
                                 "#E78AC3" ,"#66C2A5" ,"#6699CC",
                                 "#E5C494","#CC6677"),
            lineages = paste0("Lineage", want_lineage), lineages_span =0.3)

CellDimPlot(monocle, group.by = "invivo_celltype", reduction = "umap",palcolor =  c("#66C2A5" ,"#6699CC","#CC6677","grey"),
            lineages_palcolor =c("#8DA0CB", "#888888","#999933",
                                 "#E78AC3" ,"#66C2A5" ,"#6699CC",
                                 "#E5C494","#CC6677"),
            lineages = paste0("Lineage", want_lineage), lineages_span =0.3)
#8x8

FeatureDimPlot(monocle, features =c("CHIR_ac_2_AUC"),lower_cutoff = 0.04,upper_cutoff = 0.08,
               reduction = "umap", theme_use = "theme_blank")

FeatureDimPlot(monocle,raster =F, pt.size = .1,
               features = c("ESC_AUC"),lower_cutoff = 0.1,upper_cutoff = 0.18,
               reduction = "umap", theme_use = "theme_blank")

FeatureDimPlot(monocle,raster = F, pt.size = .1,
               features = c("TSC_AUC"),lower_cutoff = 0.13,upper_cutoff = 0.19,
               reduction = "umap", theme_use = "theme_blank")

FeatureDimPlot(monocle,raster = F, pt.size = .1,
               features = c("PGCLC_AUC"),lower_cutoff = 0.195,upper_cutoff = 0.26, 
               reduction = "umap", theme_use = "theme_blank")

FeatureDimPlot(monocle,raster = F, pt.size = .1,
               features = c("Primed_AUC"),lower_cutoff = 0.165,upper_cutoff = 0.225,  
               reduction = "umap", theme_use = "theme_blank")


FeatureDimPlot(monocle,raster = F, pt.size = .1,
               features = c("Bmp_sig_AUC"),lower_cutoff = 0.05,upper_cutoff = 0.28,  
               reduction = "umap", theme_use = "theme_blank")
FeatureDimPlot(monocle,raster = F, pt.size = .1,
               features = c("Wnt_sig_AUC"),lower_cutoff = 0.1,upper_cutoff = 0.25,  
               reduction = "umap", theme_use = "theme_blank")
FeatureDimPlot(monocle,raster = F, pt.size = .1,
               features = c("Hippo_sig_AUC"),lower_cutoff = 0.08,upper_cutoff = 0.28,  
               reduction = "umap", theme_use = "theme_blank")

#8x8

ggsave(
  #filename = "K5_D0-5_umap_lineage_sample.pdf",
  filename = "K5_D0-5_umap_lineage_split.png",
  width = 20,             # 宽
  height = 20,            # 高
  units = "in",          # 单位
  dpi = 1000              # 分辨率DPI
)



FeatureDimPlot(monocle,raster = F, pt.size = .1,lower_cutoff = 0,upper_cutoff =3,  
               features = signal$gene,keep_scale = NULL,
               reduction = "umap", theme_use = "theme_blank")



Dynamic_AUC_and_marker=function(){}


monocle@meta.data$Klf5.=monocle@assays$RNA$data[rownames(monocle@assays$RNA$data)=="Klf5",]

pancreas_sub@meta.data$Klf4.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Klf4",]
pancreas_sub@meta.data$Tbx3.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Tbx3",]
pancreas_sub@meta.data$Prdm1.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Prdm1",]
pancreas_sub@meta.data$Nanos3.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Nanos3",]
pancreas_sub@meta.data$Gata3.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Gata3",]
pancreas_sub@meta.data$Elf5.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Elf5",]
pancreas_sub@meta.data$Cdx2.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Cdx2",]
pancreas_sub@meta.data$Cdcp1.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Cdcp1",]
pancreas_sub@meta.data$Wnt3a.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Wnt3a",]
pancreas_sub@meta.data$Axin2.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Axin2",]
pancreas_sub@meta.data$Fgf5.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Fgf5",]
pancreas_sub@meta.data$Dnmt3b.=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Dnmt3b",]


ht2 <- DynamicHeatmap(
  srt = monocle_annoed,limits =c(-2,2),exp_method="raw",
  lineages = paste0("Lineage", want_lineage),heatmap_palette  = "RdBu",
  features = c("Bmp_sig_scale","Fgf_sig_scale","Hippo_sig_scale","Nodal_sig_scale",
               "Notch_sig_scale","Wnt_sig_scale","Hedgehog_sig_scale"),
  cell_annotation = "sample",width=2,height=3,raster_by_magick=F,use_raster = F,
)
ht3 <- DynamicHeatmap(
  srt = monocle_annoed,limits =c(-4,4),exp_method="raw",
  lineages = paste0("Lineage", want_lineage),heatmap_palette  = "RdBu",
  features = c("ESC_scale","PGCLC_scale","TSC_scale","Primed_scale","CHIR_ac_2_scale"),
  cell_annotation = "sample",width=2,height=2,raster_by_magick=F,use_raster = F
)
ht2 <- DynamicHeatmap(
  srt = monocle,limits =c(0,2),
  exp_method="raw",
  lineages = paste0("Lineage", want_lineage),heatmap_palette  = "RdBu",
  features = c("Klf5."),
  cell_annotation = "sample",width=2,height=1#raster_by_magick=T,use_raster = T
)
panel_fix(ht2$plot, height = 4, width = 30, raster = TRUE, dpi = 1000)

panel_fix(ht3$plot, height = 4, width = 30, raster = TRUE, dpi = 1000)





ht3 <- DynamicHeatmap(
  srt = pancreas_sub,limits =c(-4,4),
  lineages = paste0("Lineage", want_lineage),
  # features = c("ESC_AUC","Klf4.","Tbx3.",
  #              "PGCLC_AUC","Prdm1.","Nanos3.",
  #              "TSC_AUC","Gata3.","Elf5.","Cdcp1.","Cdx2.",
  #              "CHIR_ac_2_AUC","Wnt3a.","Axin2.",
  #              "Primed_AUC","Fgf5.","Dnmt3b."
  #              ),#n_split = 5,
  features=c("ESC_AUC","PGCLC_AUC","TSC_AUC"),
  show_row_names=T,heatmap_palette  = "RdBu",
  #split_method = "kmeans",
  cell_annotation = "sample",width=2,height=3,raster_by_magick=T,use_raster = T
)
print(ht2$plot)
saveRDS(ht2,"DynamicHeatmap_step1.rds")
ggsave(
  filename = "DynamicHeatmap_PGCLC_AUC_marker.pdf",
  width = 15,             # 宽
  height = 3,            # 高
  units = "in",          # 单位
  dpi = 1000              # 分辨率DPI
)



FeatureDimPlot(monocle,
               #features = c("ESC_scale","TSC_scale","PGCLC_scale","Primed_scale"),
               features = c("Bmp_sig_AUC","Fgf_sig_AUC","Hippo_sig_AUC","Nodal_sig_AUC",
                          "Notch_sig_AUC","Wnt_sig_AUC","Hedgehog_sig_AUC"),
               #features = c("Bmp_sig_scale","Fgf_sig_scale","Hippo_sig_scale","Nodal_sig_scale",
              #              "Notch_sig_scale","Wnt_sig_scale","Hedgehog_sig_scale"),
               pt.size = .3, reduction="umap",lower_cutoff = 0,upper_quantile = 1,
               palette = "RdBu",#upper_cutoff = 4,lower_cutoff = -1,
               #color_blend_mode = "blend",  label_insitu = TRUE,title = "blend",compare_features = TRUE,pt.alpha=0.75,
               #palcolor = c("#51FF51","red","blue","purple"),
               label = F)+
  
FeaturePlot(monocle,features = signal$gene,pt.size = .3, reduction="umap")



pancreas_sub@meta.data$Cdx2_data=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Cdx2",]
pancreas_sub@meta.data$Cdcp1_data=pancreas_sub@assays$RNA$data[rownames(pancreas_sub@assays$RNA$data)=="Cdcp1",]
FeatureDimPlot(pancreas_sub,
               features = c("Cdx2_data","Cdcp1_data"),
               pt.size = .1, reduction="umap",#slot = "scale.data",
               compare_features = TRUE,pt.alpha=0.55,
               color_blend_mode = "blend", title = "blend",
               palcolor = c("green","red"),
               label = F, label_insitu = TRUE)
FeaturePlot(pancreas_sub, features = c("Cdx2","Cdcp1"),order=T,
            #cols = c("lightgrey", "blue", "blue"),
            blend=T,blend.threshold=0)
  #scale_color_gradientn(limites=c(-4,4))
ggsave(
  filename = "K5_D0-5_umap_signal_pathway_markergene.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 50,             # 宽
  height = 50,            # 高
  units = "in",          # 单位
  dpi = 1000  ,limitsize = F           # 分辨率DPI
)
use=c("Klf5","Sox15","Zfp42","Esrrb",
      "Tead4","Nanog","Pou5f1","Klf2","Nr5a2","Dppa5a","Tfap2c","Rhox6",
"Rhox9","Msx1",
"Klf4","Tfcp2l1","Tbx3",
"Gata3","Cdx2","Plet1","Elf5","Gata2","Eomes",
"Prdm1","Dppa3","Nanos3",
"Fgf5","Otx2","Dnmt3a","Dnmt3b","T","Lefty1",
"Gata6","Foxa2","Cfc1","Sox7","Gata5",
"T","Msx2","Wnt3","Cebpa","Bmp6","Cdx4",
"Tcf15","Pax6","Dsp","Nes"
)
FeatureDimPlot(pancreas_sub,c("Zyx","Hoxb1","Sox11","Serpinh1","Fgf17",
                              "Scd2","Mdk","Hspa5","Epha2",formative) ,
               reduction = "umap",pt.size = .1)
ggsave(
  filename = "K5_D0-5_umap_formative.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 50,             # 宽
  height = 50,            # 高
  units = "in",          # 单位
  dpi = 500  ,limitsize=F            # 分辨率DPI
)

#monocle2####

# pancreas_sub_meta=pancreas_sub@meta.data
# pancreas_sub_counts=pancreas_sub@assays$RNA@counts
# pancreas_sub_data=pancreas_sub@assays$RNA@data
# pancreas_sub=NULL
# 
# monocle=CreateSeuratObject(counts = pancreas_sub_counts)
# monocle@meta.data=pancreas_sub_meta
# monocle@assays$RNA$data=pancreas_sub_data
#monocle=NormalizeData(monocle)
#pancreas_sub_meta=monocle@meta.data
#monocle[["RNA"]] <- as(object = monocle[["RNA"]], Class = "Assay")


monocle=readRDS("monocle_for_PAGA.rds")

monocle=CreateSeuratObject(monocle@assays$RNA@counts,meta.data=monocle@meta.data)

monocle_repgrm=monocle[,monocle$leiden%in%c(40,21,6,4,32,34,35,37,23)]

expr=as(as.matrix(monocle_repgrm@assays$RNA$counts),"sparseMatrix")
cellinfo=monocle_repgrm@meta.data
f_data=data.frame(gene_short_name=rownames(monocle_repgrm),row.names = rownames(monocle_repgrm))

cellinfo=new("AnnotatedDataFrame",data=cellinfo)
f_data=new("AnnotatedDataFrame",data=f_data)
cds=newCellDataSet(expr,
                   phenoData = cellinfo,
                   featureData=f_data,
                   lowerDetectionLimit=0.5,
                   expressionFamily=negbinomial.size()#for raw counts
)


cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)

cds=detectGenes(cds,min_expr=0.1)
#print(head(fData(cds)))
expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))
express_genes=VariableFeatures(FindVariableFeatures(monocle_repgrm,nfeatures = 2000))
express_genes=intersect(expressed_genes,express_genes)

cds=setOrderingFilter(cds,express_genes)
print(plot_ordering_genes(cds))

cds = reduceDimension(cds,
                      max_components=3,
                      method = "DDRTree")

cds=orderCells(cds)


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$leiden)[,"40"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))



cds$leiden=as.character(cds$leiden)
plot_cell_trajectory(cds,color_by="leiden",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="invivo_celltype",size=1,show_backbone=T)





Run_compare_2_branch_by_monocle2=
  function(monocle,branch1_vec,branch2_vec){
    
    monocle_meta=monocle@meta.data
    for(i in 1:length(branch1_vec)){
      branch1=branch1_vec[i]
      branch2=branch2_vec[i]
      
      message("Comparing ",paste0("Lineage",c(branch1))," and ",paste0("Lineage",c(branch2)))
      ln1=na.omit(as.data.frame(monocle[[paste0("Lineage",c(branch1))]]))
      ln2=na.omit(as.data.frame(monocle[[paste0("Lineage",c(branch2))]]))
      lineage_sub=monocle[,colnames(monocle)%in%c(rownames(ln1),rownames(ln2))]
      #lineage_sub
      expr=as(as.matrix(lineage_sub@assays$RNA@counts),"sparseMatrix")
      cellinfo=lineage_sub@meta.data
      f_data=data.frame(gene_short_name=rownames(lineage_sub),row.names = rownames(lineage_sub))
      
      cellinfo=new("AnnotatedDataFrame",data=cellinfo)
      f_data=new("AnnotatedDataFrame",data=f_data)
      cds=newCellDataSet(expr,
                         phenoData = cellinfo,
                         featureData=f_data,
                         lowerDetectionLimit=0.5,
                         expressionFamily=negbinomial.size()#for raw counts
      )
      
      
      cds=estimateSizeFactors(cds)
      cds=estimateDispersions(cds)
      lineage_cl=cds$lineage_cl
      lineage_cl[lineage_cl==0]=1
      cds$lineage_cl=lineage_cl
      
      print("check lineage_cl")
      print(table(cds$lineage_cl))
      
      
      cds=detectGenes(cds,min_expr=0.1)
      #print(head(fData(cds)))
      expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))
      express_genes=VariableFeatures(FindVariableFeatures(lineage_sub,nfeatures = 2000))
      express_genes=intersect(expressed_genes,express_genes)
      message("Use ",length(express_genes)," Variable Features")
      cds=setOrderingFilter(cds,express_genes)
      print(plot_ordering_genes(cds))
      
      #因为根据两个lineage之间的差异基因所以自然可以分成两个分支
      
      cds = reduceDimension(cds,
                            max_components=3,
                            method = "DDRTree")
      
      cds=orderCells(cds)
      
      
      GM_state <- function(cds){
        if (length(unique(pData(cds)$State)) > 1){
          T0_counts <- table(pData(cds)$State, pData(cds)$lineage_cl)[,"1"]
          return(as.numeric(names(T0_counts)[which
                                             (T0_counts == max(T0_counts))]))
        } else {
          return (1)
        }
      }
      cds <- orderCells(cds, root_state = GM_state(cds))
      
      
      message("Saving cds to ",paste("cds",branch1,branch2,".rds",sep = "_"))
      saveRDS(cds,paste("cds",branch1,branch2,".rds",sep = "_"))
      message("Saving lineage_sub to ",paste("lineage_sub",branch1,branch2,".rds",sep = "_"))
      saveRDS(lineage_sub,paste("lineage_sub",branch1,branch2,".rds",sep = "_"))
      
      plot_cell_trajectory(cds,color_by="lineage_cl",size=1,show_backbone=T)
      ggsave(
        filename = paste("monocle2_cds",branch1,branch2,"lineage_cl.png",sep = "_"), 
        width = 10,   height = 10,   units = "in",  dpi = 500)
      plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=T)
      ggsave(
        filename = paste("monocle2_cds",branch1,branch2,"Pseudotime.png",sep = "_"),  
        width = 10,   height = 10,   units = "in",  dpi = 500)
      plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=T)
      ggsave(
        filename = paste("monocle2_cds",branch1,branch2,"State.png",sep = "_"),  
        width = 10,   height = 10,   units = "in",  dpi = 500)
      
      meta=pData(cds)
      meta=meta[,c("State","Pseudotime")]
      colnames(meta)=c(paste0("State_",branch1,"v",branch2),paste0("Pseudotime_",branch1,"v",branch2))
      
      monocle_meta=Btable_insert_into_Atable(monocle_meta,meta)
      print(colnames(monocle_meta))
    }
    return(monocle_meta)
  }


monocle_meta=Run_compare_2_branch_by_monocle2(monocle,
                                              branch1_vec=c(9, 9, 9, 9, 16,16,16,14,14,10),
                                              branch2_vec=c(10,12,14,16,10,12,14,10,12,12)
                                              )

cds=readRDS("cds_10_12_.rds")
branch1=10
branch2=12



cds=setOrderingFilter(cds,deg_fin_mat_filter$gene)
plot_ordering_genes(cds)
cds = reduceDimension(cds,
                      max_components=2,
                      method = "DDRTree")

cds=orderCells(cds)


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$lineage_cl)[,"1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))

color_fun=data.frame(seurat_clusters=0:38,color=palette_scp(n=39,palette = "Paired"))

plot_cell_trajectory(cds,color_by="seurat_clusters",size=1,show_backbone=T,show_branch_points=F)+
  scale_colour_manual(
    values =color_fun$color[which(table(cds$seurat_clusters)!=0)] )

ggsave(
  filename = paste("monocle2_cds",branch1,branch2,"seurat_clusters.png",sep = "_"), 
  width = 10,   height = 10,   units = "in",  dpi = 500)



plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=T)
ggsave(
  filename = paste("monocle2_cds",branch1,branch2,"Pseudotime.png",sep = "_"),  
  width = 10,   height = 10,   units = "in",  dpi = 500)
plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=T)
ggsave(
  filename = paste("monocle2_cds",branch1,branch2,"State.png",sep = "_"),  
  width = 10,   height = 10,   units = "in",  dpi = 500)



meta=pData(cds)
meta=meta[,c("State","Pseudotime")]
colnames(meta)=c(paste0("State_",branch1,"v",branch2),paste0("Pseudotime_",branch1,"v",branch2))

monocle_meta=Btable_insert_into_Atable(monocle_meta,meta)
print(colnames(monocle_meta))



make_lineage=
  function(meta,b1,b2,b1_state,b2_state,cf1,cf2){
  
branch=paste0("Lineage_",b1,"v",b2,"_branch")
l1=paste0("Lineage_",b1,"v",b2,"_",b1)
l2=paste0("Lineage_",b1,"v",b2,"_",b2)

state=paste0("State","_",b1,"v",b2)
pseudo=paste0("Pseudotime","_",b1,"v",b2)

meta[meta[,state]%in%b1_state,l1]=meta[meta[,state]%in%b1_state,pseudo]
meta[meta[,state]%in%b2_state,l2]=meta[meta[,state]%in%b2_state,pseudo]

meta[,branch]=NA
meta[meta[,state]%in%setdiff(b1_state,b2_state),branch]=cf1
meta[meta[,state]%in%setdiff(b2_state,b1_state),branch]=cf2
meta[meta[,state]%in%intersect(b1_state,b2_state),branch]="Pre-branch"

return(meta)
}

meta=make_lineage(9,10,c(5,2,4,3),c(5,1),"Day1_Diff_CellFate","ESC_CellFate")
meta=make_lineage(9,12,c(3,1),c(3,2),"Day1_Diff_CellFate","PGC_CellFate")
meta=make_lineage(9,14,c(3,2,7,1),c(3,4,6,5),"Day1_Diff_CellFate","TSC2_CellFate")
meta=make_lineage(9,16,c(3,2,7,1),c(3,2,6,5),"Day1_Diff_CellFate","TSC_CellFate")

meta=make_lineage(16,10,c(1,3),c(1,2),"TSC_CellFate","ESC_CellFate")
meta=make_lineage(16,12,c(1,3,4,5),c(1,2),"TSC_CellFate","PGC_CellFate")
meta=make_lineage(16,14,c(4,3),c(4,2,5,1),"TSC_CellFate","TSC2_CellFate")

meta=make_lineage(14,10,c(1,2),c(1,3),"TSC2_CellFate","ESC_CellFate")
meta=make_lineage(14,12,c(4,2,5,1),c(4,3),"TSC2_CellFate","PGC_CellFate")

meta=make_lineage(meta,10,12,c(2,1),c(2,3),"ESC_CellFate","PGC_CellFate")



meta=normalize_pseudotime2(meta,c(#"_9v10_9","_9v10_10","_9v12_9","_9v12_12",
                                  #"_9v14_9","_9v14_14","_9v16_9","_9v16_16",
                                  #"_16v10_16","_16v10_10","_16v12_16","_16v12_12","_16v14_16","_16v14_14",
                                  #"_14v10_14","_14v10_10","_14v12_14","_14v12_12",
                                  "_10v12_10","_10v12_12"
                                  ))
#saveRDS(meta,"monocle_meta_maked_lineage.rds")



deg=read.table("train.monocle.DEG.txt",header = T)

ordergene=deg[deg$qval==0,5]
#ordergene=deg[deg$use_for_ordering==T,5]
length(ordergene)
#ordergene=rownames(deg)[order(deg$qval)][1:2000]
cds=setOrderingFilter(cds,ordergene)


#select high variable feature from
plot_ordering_genes(cds)
cds=detectGenes(cds,min_expr=0.1)
print(head(fData(cds)))
expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))


#monocle deg
cds=detectGenes(cds,min_expr=0.1)
#print(head(fData(cds)))
expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))

diff=differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~lineage_cl",cores=30)
#head(diff)

deg=subset(diff,qval<0.01)
deg=deg[order(deg$qval,decreasing=F),]
head(deg)

write.table(deg,"train.monocle.DEG_9_10.txt",col.name=T,row.names=F,sep="\t",quote=F)
#deg=read.table("train.monocle.DEG.txt",header = T)


ordergene=deg[deg$qval==0,5]
#ordergene=deg[deg$use_for_ordering==T,5]
length(ordergene)
#ordergene=rownames(deg)[order(deg$qval)][1:2000]
cds=setOrderingFilter(cds,ordergene)

plot_ordering_genes(cds)



#select high variable feature from seurat
expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))
express_genes=VariableFeatures(FindVariableFeatures(lineage_sub,nfeatures = 2000))
#express_genes=VariableFeatures(lineage_sub)
express_genes=intersect(expressed_genes,express_genes)
length(express_genes)
cds=setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)



#select clusters DEG
Idents(lineage_sub)="lineage_cl"
deg.cluster=FindMarkers(lineage_sub,ident.1 = 1,ident.2 = c(31,),only.pos = T,logfc.threshold = 0.25)
express_genes=unique(subset(deg.cluster,p_val_adj<1e-10 & abs(avg_log2FC)>3)$gene)
saveRDS(deg.cluster,"deg.cluster.rds")
express_genes=unique(deg.cluster[order(deg.cluster$p_val_adj),]$gene)[1:2000]
length(express_genes)
cds=setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)


#select high variable feature by monocle #no commend
disp_table=dispersionTable(cds)
disp.genes=subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds=setOrderingFilter(cds,disp.genes)
plot_ordering_genes(cds)


disp.genes=ht1$feature_metadata$features
cds=setOrderingFilter(cds,disp.genes)
plot_ordering_genes(cds)


#reduce dimension
cds = reduceDimension(cds,
                      max_components=3,
                      method = "DDRTree")

cds=orderCells(cds)


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$lineage_cl)[,"1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))


cds=readRDS("cds_9_10.rds")
lineage_sub=readRDS("lineage_sub_9_10.rds")


#cds=orderCells(cds,root_state=5) set the start point 5


#plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="TSC_AUC",size=1,show_backbone=T)
#plot_cell_trajectory(cds,color_by="seurat_clusters",size=1,show_backbone=T)
#plot_cell_trajectory(cds,color_by="sample")+facet_wrap("~sample",nrow=1)

plot_cell_trajectory(cds,color_by="lineage_cl",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="seurat_clusters",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=T)

gc()



BEAM_res <- BEAM(cds, branch_point = 2, cores = 1, progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval <  5.170534e-164)),],
                            branch_point = 2,
                            num_clusters = 5,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

meta=pData(cds)


meta[meta$State%in%c(5,2,4,3),"Lineage_9v10_9"]=meta[meta$State%in%c(5,2,4,3),"Pseudotime"]
meta[meta$State%in%c(5,1),"Lineage_9v10_10"]=meta[meta$State%in%c(5,1),"Pseudotime"]


normalize_pseudotime2=
  function(meta,want_lineage){
    
    for(i in want_lineage){
      nor_pseudo=as.numeric(meta[[paste0("Lineage",i)]])
      meta[[paste0("Lineage",i)]]=c(nor_pseudo-min(na.omit(nor_pseudo)))/c(max(na.omit(nor_pseudo))-min(na.omit(nor_pseudo)))
      
    }
    
    return(meta)
  }

meta=normalize_pseudotime2(meta,c("_9v10_9","_9v10_10"))

meta[meta$State%in%c(2,4,3),"Lineage_9v10_9only"]=meta[meta$State%in%c(2,4,3),"Lineage_9v10_9"]
meta[meta$State%in%c(1),"Lineage_9v10_10only"]=meta[meta$State%in%c(1),"Lineage_9v10_10"]
meta[meta$State%in%c(5,2,4,3),"Lineage_9v10_9"]=meta[meta$State%in%c(5,2,4,3),c("Lineage_9v10_9")]
meta[meta$State%in%c(5,2,4,3),"Lineage_9v10_10"]=meta[meta$State%in%c(5,2,4,3),c("Lineage_9v10_10")]

meta$lineage_9v10_branch=1
meta[meta$State%in%c(2,4,3),"lineage_9v10_branch"]="Lineage_9_cellfate"
meta[meta$State%in%c(1),"lineage_9v10_branch"]="Lineage_10_cellfate"
meta[meta$State%in%c(5),"lineage_9v10_branch"]="Lineage_prebranch"


pancreas_sub_meta=Btable_insert_into_Atable(pancreas_sub_meta,meta[,c(42:45,48)])
colnames(pancreas_sub_meta)
pancreas_sub_meta[,39:42]=apply(pancreas_sub_meta[,39:42],2,as.numeric)

monocle@meta.data=pancreas_sub_meta



##load from cluster####
sub_cds=readRDS("sub_cds_5000_random.rds")
sub_seurat_merged=readRDS("sub_seurat_5000_random.rds")

sub_seurat_merged=sub_seurat_merged[,colnames(sub_seurat_merged)%in%colnames(sub_cds)]
sub_seurat_merged
dim(sub_cds)


#sub_seurat_merged <- FindVariableFeatures(sub_seurat_merged, selection.method = "vst", nfeatures = 2000)
#sub_seurat_merged <- RunPCA(sub_seurat_merged, features = VariableFeatures(object = seurat_merged))
sub_seurat_merged <- FindNeighbors(sub_seurat_merged, dims = 1:25)
sub_seurat_merged <- FindClusters(sub_seurat_merged, resolution = 0.5)
table(Idents(sub_seurat_merged))
gc()
sub_seurat_merged <- RunUMAP(sub_seurat_merged, dims = 1:30,n.neighbors = 20, 
                             n.epochs = 200,min.dist = 0.1,spread=1#lager will be mixed
)

color=brewer.pal(n = 9, name = "Set1")
options(repr.plot.width=15,repr.plot.height=15)
p1=DimPlot(sub_seurat_merged, reduction = "umap",label = T,#group.by = "sample",
        label.size = 10)
p2=DimPlot(sub_seurat_merged, reduction = "umap",label = T,group.by = "sample",
        label.size = 10)
p1+p2

options(repr.plot.width=15,repr.plot.height=15)
FeaturePlot(sub_seurat_merged,reduction = "umap",features = c("ESC_AUC", "TSC_AUC","PGCLC_AUC","Primed_AUC"
                                                              # "ET_AUC","EP_AUC","TP_AUC"
),ncol=2,
cols = rev(brewer.pal(n = 10, name = "RdBu")))
options(repr.plot.width=13,repr.plot.height=13)
FeaturePlot(sub_seurat_merged, reduction = "umap",pt.size = 0.001,features = c("Cdx2", "Gata3","Elf5",
                                                                               "Prdm1","Dppa3","Nanos3",
                                                                               "Tbx3","Tfcp2l1","Dppa5a",
                                                                               "Pou5f1","Nanog","Fgf5"),
            ncol=3,max.cutoff="q75")


##select clusters DEG1####
#select high variable feature from ?
sub_cds=detectGenes(sub_cds,min_expr=1)
print(head(fData(sub_cds)))
expressed_genes=rownames(subset(fData(sub_cds),num_cells_expressed >=10))

diff=differentialGeneTest(sub_cds[expressed_genes,],fullModelFormulaStr="~sample",cores=1)
head(diff)

deg=subset(diff,qval<0.01)
deg=deg[order(deg$qval,decreasing=F),]
head(deg)

write.table(deg,"sub_cds_train.monocle.DEG.txt",col.name=T,row.names=F,sep="\t",quote=F)
deg=read.table("sub_cds_train.monocle.DEG.txt",header = T)


ordergene=deg[deg$qval==0,5]
#ordergene=deg[deg$use_for_ordering==T,5]
length(ordergene)
#ordergene=rownames(deg)[order(deg$qval)][1:2000]
sub_cds=setOrderingFilter(sub_cds,ordergene)

plot_ordering_genes(sub_cds)

##select clusters DEG2####
deg_seurat=sub_seurat_merged[,sub_seurat_merged@meta.data$seurat_clusters%in%c("10","8","4","1","0")]
table(Idents(deg_seurat))
deg.cluster=FindAllMarkers(deg_seurat,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(deg.cluster,"sub_seurat_merged_random_5000_findallmarker_deg.cluster_c014810cells.txt")
deg.cluster=read.table("sub_seurat_merged_random_5000_findallmarker_deg.cluster_c014810cells.txt",header = T)
express_genes=subset(deg.cluster,deg.cluster$p_val<0.01)$gene
length(express_genes)
express_genes=express_genes[1:2000]
sub_cds=setOrderingFilter(sub_cds,express_genes)
length(express_genes)
plot_ordering_genes(sub_cds)


##DDRTree####
sub_cds = reduceDimension(sub_cds,max_components=5,verbose = F,method="ICA")
sub_cds=orderCells(sub_cds)
sub_cds <- orderCells(sub_cds,num_paths=10, root_state = GM_state(sub_cds))
gc()


options(repr.plot.width=15,repr.plot.height=15)
gc()
#plot_cell_trajectory(sub_cds,color_by="Pseudotime",size=1,show_backbone=T)
plot_cell_trajectory(sub_cds,color_by="TSC_AUC",size=1,show_backbone=T)
plot_cell_trajectory(sub_cds,color_by="sample",size=0.1,show_backbone=T)
plot_cell_trajectory(sub_cds,color_by="State")+facet_wrap("~State",nrow=4)

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$sample)[,"EpiSCs"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}




seurat_merged

seurat=readRDS("K5_D0-5.rds")
gc()

seurat_merged <- RunUMAP(seurat_merged, dims = 1:30,n.neighbors = 20, 
                         n.epochs = 200,min.dist = 0.1,spread=1#lager will be mixed
)

DimPlot(seurat_merged, reduction = "umap",label = T,group.by = "sample",
        label.size = 8)

FeaturePlot(seurat_merged,reduction = "umap",features = c("ESC_AUC", "TSC_AUC","PGCLC_AUC","Primed_AUC"
                                                          # "ET_AUC","EP_AUC","TP_AUC"
),ncol=2,
cols = rev(brewer.pal(n = 10, name = "RdBu")))

FeaturePlot(seurat_merged, reduction = "umap",features = c("Cdx2", #"Gata3","Elf5",
                                                                           "Prdm1","Dppa3",#"Nanos3",
                                                                           #"Tbx3","Tfcp2l1","Dppa5a",
                                                                           #"Pou5f1","Nanog",
                                                                           "Fgf5"),
            ncol=2,max.cutoff="q75")




seurat=readRDS("K5_D0-5.rds")


Idents(seurat) <- "sample"
# Downsample the number of cells per identity class
random_seurat<- subset(x = seurat, downsample = 10000)
random_seurat
table(random_seurat@meta.data$sample)


seurat_merged=random_seurat
table(seurat_merged@meta.data$sample)
seurat_merged@assays$RNA@counts[1:5,1:5]
random_seurat=NULL
gc()

expr=as(as.matrix(seurat_merged@assays$RNA@counts),"sparseMatrix")
cellinfo=seurat_merged@meta.data
f_data=data.frame(gene_short_name=rownames(seurat_merged),row.names = rownames(seurat_merged))
expr[1:5,1:5]

#cellinfo=new("AnnotatedDataFrame",data=cellinfo)
#f_data=new("AnnotatedDataFrame",data=f_data)
cds <- new_cell_data_set(expr,
                         cell_metadata = cellinfo,
                         gene_metadata = f_data)
expr=NULL
gc()


cds <- preprocess_cds(cds, num_dim = 150)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,reduction_method = "UMAP",
                        umap.n_neighbors = 5,
                        umap.min_dist = 0.1,
                        cores = 2)
plot_cells(cds,color_cells_by = "sample")


ciliated_genes <- c("Cdx2", "Gata3","Elf5",
                    "Prdm1","Dppa3","Nanos3",
                    "Tbx3","Tfcp2l1","Dppa5a",
                    "Pou5f1","Nanog","Fgf5")
plot_cells(cds,
           genes=c("Fgf5"),
           label_cell_groups=FALSE,
           show_trajectory_graph=T)

cds <- cluster_cells(cds,k = 4,reduction_method = "UMAP")
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "sample",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# Now that we have a sense of where the early cells fall,
# we can call order_cells(), which will calculate where 
# each cell falls in pseudotime. In order to do so 
# order_cells()needs you to specify the root nodes of 
# the trajectory graph. If you don't provide them as an
# argument, it will launch a graphical user interface 
# for selecting one or more root nodes.
cds <- order_cells(cds)
cds_sub <- choose_graph_segments(cds,clear_cds =T)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)



plot_cells(cds,
           color_cells_by = "TSC_AUC",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)







cds <- as.cell_data_set(random_seurat)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "sample",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds_sub)
cds_sub <- choose_graph_segments(cds,clear_cds =T)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)



plot_cells(cds,genes =c("Fgf5") ,
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)









#seurat visulization####
seurat=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/K5_D0-5.rds")



raw_seurat=seurat
seurat=raw_seurat[,raw_seurat@meta.data$sample%in%c("D1_K5","D2_K5","EpiSCs")]



seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
setdiff(study_gene,VariableFeatures(seurat))


##select dim####
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

# Determine percent of variation associated with each PC
pct <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

seurat <- FindNeighbors(seurat, dims = 1:pcs)

##select resolution and seed####

for(i in 4:100){
  for(re in seq(0.7,0.9,0.01)){
seurat <- FindClusters(seurat, resolution =re, random.seed = i)
filename=paste("D5_K245_findcl_seed",i,"re",re,".png",sep="_")

p4=DimPlot(seurat, reduction = "tsne",label = T,group.by = "seurat_clusters",
           label.size =4)
p5=DimPlot(seurat, reduction = "tsne",label = T,group.by = "sample",
           label.size =4)
p1=FeaturePlot(seurat,reduction = "tsne",features = c("ESC_AUC"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
p2=FeaturePlot(seurat,reduction = "tsne",features = c("TSC_AUC"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
p3=FeaturePlot(seurat,reduction = "tsne",features = c("PGCLC_AUC"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
P6=VlnPlot(seurat,pt.size = 0, features = c("ESC_AUC"))
P7=VlnPlot(seurat,pt.size = 0, features = c("TSC_AUC"))
P8=VlnPlot(seurat,pt.size = 0, features = c("PGCLC_AUC"))
ggpubr::ggarrange(p1,p2,p3,P6,P7,P8,p4,p5,ncol = 3,nrow =3)
ggsave(
  filename = filename, # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 15,             # 宽
  height = 11,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)
system.time(3)
  }
}

##check and decide cluster parameter####
st=seurat[,seurat$seurat_clusters=="11"]
table(st$sample)/table(seurat$sample)[1:3]

df=data.frame(table(st$sample))
colnames(df)=c("condition","num")
Plot_Pie(df,color =rev(c("#90EE90","#1E90FF","#FC8D62")))



p1=FeaturePlot(st,reduction = "tsne",features = c("Cdx2"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

p=DotPlot(st,group.by = "sample",features=c("TSC_AUC","Cdx2","Gata3","Elf5","Krt7"))
p+ggplot2::coord_flip()

seurat <- FindClusters(seurat, resolution =0.73, random.seed = 3)
table(Idents(seurat))>length(colnames(seurat))*0.005
gc()


##modif UMAP or tSNE####

seurat@meta.data$sample=factor(seurat@meta.data$sample,levels =c("D5_K5","D5_K2","D5_K4","EpiSCs"))
seurat <- RunUMAP(seurat, dims = 1:pcs,n.neighbors = 20, 
                  n.epochs = 200,min.dist = 0.2,spread=1#lager will be mixed
)
seurat <- RunTSNE(seurat, dims = 1:pcs,seed.use=5)

redu="umap"
redu="tsne"

p1=DimPlot(seurat, reduction = redu,label = T,group.by = "seurat_clusters",
           label.size =4)
p2=DimPlot(seurat, reduction = redu,label = T,group.by = "sample",cols = c("#FC8D62","#90EE90","#1E90FF","grey"),
           label.size = 0,)
p3=FeaturePlot(seurat,reduction = redu,features = c("ESC_AUC"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
p4=FeaturePlot(seurat,reduction = redu,features = c("TSC_AUC"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
p5=FeaturePlot(seurat,reduction = redu,features = c("PGCLC_AUC"),min.cutoff = "q25")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
p6=FeaturePlot(seurat,reduction = redu,features = c("Primed_AUC"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,ncol = 3,nrow =2)
ggsave(
  filename = "D5_K245_cl_ETC_AUC.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 14,             # 宽
  height = 8,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)



##check gene in UMAP or tSNE####
FeaturePlot(seurat, reduction = redu,features = study_gene,
            ncol=7,max.cutoff="q95")
ggsave(
  filename = "K5_1-2D_cl_all_marker_gene.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 35,             # 宽
  height = 20,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)




##AUC and Markergene dotplot####
p1=DotPlot(seurat,scale = F,features=c("PGCLC_AUC","ESC_AUC","TSC_AUC"),dot.scale = 12)+
  scale_colour_gradientn(limits=c(0,0.25),colours = c( rev(brewer.pal(n = 5, name = "Greys")),
                                                         brewer.pal(n = 5, name = "RdPu")))
p1=p1+ggplot2::coord_flip()
p2=DotPlot(seurat,features=study_gene)+
  scale_colour_gradientn(limits=c(-3,3),colours = rev(brewer.pal(n = 10, name = "RdBu")))
p2=p2+ggplot2::coord_flip()
p1+p2


table(dot_seurat$seurat_clusters)
p1=DotPlot(seurat,scale = F,features=c("PGCLC_AUC"),dot.scale = 12)+
  scale_colour_gradientn(limits=c(0.1,0.25),colours = c( rev(brewer.pal(n = 5, name = "Greys")),
                                                        brewer.pal(n = 5, name = "RdPu")))
p1=p1+ggplot2::coord_flip()
ggsave(
  filename = "D5_K245_cl_PGCLC_AUC.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 1,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)
p2=DotPlot(seurat,features=c("Prdm1","Dppa3","Nanos3","Prdm14","Sox15"))+
  scale_colour_gradientn(limits=c(-1,3),colours = rev(brewer.pal(n = 10, name = "RdBu")))
p2=p2+ggplot2::coord_flip()
ggsave(
  filename = "D5_K245_cl_PGCLC_Markergene.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 1.7,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)


p1=DotPlot(seurat,scale = F,features=c("ESC_AUC"),dot.scale = 12)+
  scale_colour_gradientn(limits=c(0.05,0.2),colours = c( rev(brewer.pal(n = 5, name = "Greys")),
                                                         brewer.pal(n = 5, name = "RdPu")))
p1=p1+ggplot2::coord_flip()
ggsave(
  filename = "D5_K245_cl_ESC_AUC.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7.8,             # 宽
  height = 1,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)
p2=DotPlot(seurat,features=c("Tbx3","Tfcp2l1","Esrrb","Dppa5a","Zfp42"))+
  scale_colour_gradientn(limits=c(-1,3),colours = rev(brewer.pal(n = 10, name = "RdBu")))
p2=p2+ggplot2::coord_flip()
ggsave(
  filename = "D5_K245_cl_ESC_Markergene.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 1.7,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)


p1=DotPlot(seurat,scale = F,features=c("TSC_AUC"),dot.scale = 12)+
  scale_colour_gradientn(limits=c(0.12,0.18),colours = c( rev(brewer.pal(n = 5, name = "Greys")),
                                                         brewer.pal(n = 5, name = "RdPu")))
p1=p1+ggplot2::coord_flip()
ggsave(
  filename = "D5_K245_cl_TSC_AUC.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7.8,             # 宽
  height = 1,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)
p2=DotPlot(seurat,features=c("Cdx2","Gata3","Elf5","Plet1","Hand1"))+
  scale_colour_gradientn(limits=c(-1,3),colours = rev(brewer.pal(n = 10, name = "RdBu")))
p2=p2+ggplot2::coord_flip()
ggsave(
  filename = "D5_K245_cl_TSC_Markergene.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 1.7,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)


##statistic ETP like cell####
ips=c(10,12)
tsclc=c(11)
pgclc=c(17)
ips_seurat=seurat[,seurat$seurat_clusters%in%ips]
tsclc_seurat=seurat[,seurat$seurat_clusters%in%tsclc]
pgclc_seurat=seurat[,seurat$seurat_clusters%in%pgclc]

ips_accu_srt=ips_seurat[rownames(ips_seurat@assays$RNA$data)%in%c("Zfp42","Tfcp2l1"),]
ips_accu_srt=ips_accu_srt[,apply(ips_accu_srt@assays$RNA$data, 2, function(x) sum(x!=0)==2)]
table(ips_accu_srt@meta.data$sample)


tsclc_accu_srt=tsclc_seurat[rownames(tsclc_seurat@assays$RNA$data)%in%c("Gata3","Elf5"),]
tsclc_accu_srt=tsclc_accu_srt[,apply(tsclc_accu_srt@assays$RNA$data, 2, function(x) sum(x!=0)==2)]
table(tsclc_accu_srt@meta.data$sample)


pgclc_accu_srt=pgclc_seurat[rownames(pgclc_seurat@assays$RNA$data)%in%c("Prdm1","Dppa3"),]
pgclc_accu_srt=pgclc_accu_srt[,apply(pgclc_accu_srt@assays$RNA$data, 2, function(x) sum(x!=0)==2)]
table(pgclc_accu_srt@meta.data$sample)


#D5 K245 ETP distribution
df_bar=data.frame(iPSCs=table(ips_accu_srt$sample)/table(seurat$sample),
           TSCLCs=table(tsclc_accu_srt$sample)/table(seurat$sample),
           PGCLCs=table(pgclc_accu_srt$sample)/table(seurat$sample))
rownames(df_bar)=df_bar$iPSCs.Var1
df_bar=df_bar[,c(1,2,4,6)]
df_bar=reshape::melt(df_bar)
colnames(df_bar)=c("Sample","celltype","Freq")


df_bar2=data.frame(iPSCs=table(ips_accu_srt$sample),
TSCLCs=table(tsclc_accu_srt$sample),
PGCLCs=table(pgclc_accu_srt$sample))
rownames(df_bar2)=df_bar2$iPSCs.Var1
df_bar2=df_bar2[,c(1,2,4,6)]
df_bar2=reshape::melt(df_bar2)
colnames(df_bar2)=c("Sample","celltype","Num")
df_bar2$Freq=df_bar$Freq
df_bar2[df_bar2$Num<10,"Num"]=""

ggplot(data = df_bar,aes(x=Sample,y=Freq,fill=celltype))+
  geom_bar(stat = "identity",
           position = "stack")+
  scale_fill_manual(values = c("iPSCs.Freq"="#32CD32",
                               "TSCLCs.Freq"="#FFA500",
                               "PGCLCs.Freq"="#1E90FF"))+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.1)),
                     labels = scales::percent_format())+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.line = element_line(),
        legend.position = "right",
        legend.direction = "vertical")+
  labs(x=NULL,y="Cells in Total Cells (%)")+
  guides(fill=guide_legend(title = NULL,nrow = 3,byrow = FALSE))+
  geom_text(data=df_bar2,position = "stack",stat = "identity",
            aes(x=Sample,y=Freq,label=Num),vjust=-0.1)
ggsave(
  filename = "D5_K245_ETP_distribution.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4,             # 宽
  height = 3.5,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)


#ETP cluster K245 distribution
df=data.frame(table(ips_seurat$sample)[1:3])
colnames(df)=c("condition","num")
p1=Plot_Pie(df,color =rev(c("#FC8D62","#90EE90","#1E90FF")),title = "iPSCs")

df=data.frame(table(tsclc_seurat$sample)[1:3])
colnames(df)=c("condition","num")
p2=Plot_Pie(df,color =rev(c("#FC8D62","#90EE90","#1E90FF")),title = "TSCLCs")

df=data.frame(table(pgclc_accu_srt$sample)[1:3])
colnames(df)=c("condition","num")
p3=Plot_Pie(df,color =rev(c("#FC8D62","#90EE90","#1E90FF")),title = "PGCLCs")

ggpubr::ggarrange(p1,p2,p3,ncol = 3,nrow =1)
ggsave(
  filename = "ETP_K245_distribution.pdf", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 15,             # 宽
  height = 5,            # 高
  units = "in",          # 单位
  dpi = 500              # 分辨率DPI
)












FeaturePlot(seurat,reduction = "umap",features = c("Cdx2"),max.cutoff = "q75")


pData(sub_cds)$Prdm1=log2(exprs(sub_cds)["Prdm1",]+1)
p1=plot_cell_trajectory(sub_cds,color_by = "Prdm1")+scale_color_viridis()
p2=plot_cell_trajectory(sub_cds,color_by="sample",size=0.1,show_backbone=T)
p1+p2



FeaturePlot(seurat,reduction = "umap",features = c("ESC_AUC", "TSC_AUC","PGCLC_AUC","Primed_AUC"
                                                          # "ET_AUC","EP_AUC","TP_AUC"
),ncol=2,
cols = rev(brewer.pal(n = 10, name = "RdBu")))

FeaturePlot(seurat, reduction = "tsne",features = c("Cdx2", "Gata3","Elf5",
                                                           "Prdm1","Dppa3","Nanos3",
                                                           "Tbx3","Tfcp2l1","Dppa5a",
                                                           "Pou5f1","Nanog","Fgf5"),
            ncol=4,max.cutoff="q90")
VlnPlot(seurat,pt.size = 0.1, features=c("Cdx2", "Gata3","Elf5",
                              "Prdm1","Dppa3","Nanos3",
                              "Tbx3","Tfcp2l1","Dppa5a",
                              "Pou5f1","Nanog","Fgf5"))
p=DotPlot(seurat,features=c("Cdx2", "Gata3","Elf5",
                          "Prdm1","Dppa3","Nanos3",
                          "Tbx3","Tfcp2l1","Dppa5a",
                          "Pou5f1","Nanog","Fgf5"))
p=DotPlot(seurat,features=c("PGCLC_AUC","Prdm1","Dppa3","Nanos3","Tfap2c"))
p+ggplot2::coord_flip()







#AUC####
###do AUC####

signal=read.csv("mouse.signaling.cgz.csv")


chir_ac_1.5=read.csv("RFP_CL_vs_L_Day2_up_gene_heatmap_fold1.5_tpm2_1178_.csv")
chir_ac_2=read.csv("RFP_CL_vs_L_Day2_up_gene_heatmap_fold2_tpm2_578_.csv")


d2_K5_cl_1.5=read.table("Day2_CL_iKlf5 vs Day2_L_iKlf5_degene_list_1.5.txt",header = T)
d2_K5_cl_up_1.5=intersect(d2_K5_cl_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(10,12)],1,function(x) max(x)>2),]))
d2_K5_cl_dn_1.5=intersect(d2_K5_cl_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(10,12)],1,function(x) max(x)>2),]))

d2_K5_cl_1.5=readRDS("d2_K5_cl_vs_l_deg_1.5.rds")

monocle=AddModuleScore(monocle,features = list(d2_K5_cl_up_1.5=d2_K5_cl_1.5$d2_K5_cl_dn_1.5),
               name = "d2_K5_cl_dn_1.5")


VlnPlot(monocle, features = c("bmp_go_gene1"),group.by = "leiden",pt.size = 0)






k5_ac_2_bind=read.table("K5_L_anno_fin.bed")
k5_ac_2_bind=k5_ac_2_bind[k5_ac_2_bind$V5=="Bind",4]



ETP=read.csv("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/ETP_gene_WGCNA_ver3_4431_.csv")
head(ETP)
dim(ETP)
table(ETP$Cluster)

primed_relative=ETP[grepl("Primed",ETP$Cluster),]
primed_relative=primed_relative[order(primed_relative[,3],decreasing = T),1]

esc_relative=ETP[grepl("ESC",ETP$Cluster),]
esc_relative=esc_relative[order(esc_relative[,3],decreasing = T),1]

tsc_relative=ETP[grepl("TSC",ETP$Cluster),]
tsc_relative=tsc_relative[order(tsc_relative[,4],decreasing = T),1]

pgclc_relative=ETP[grepl("PGCLC",ETP$Cluster),]
pgclc_relative=pgclc_relative[order(pgclc_relative[,5],decreasing = T),1]

pgclc_relative=ETP[ETP$Cluster=="PGCLC",]
pgclc_relative=pgclc_relative[order(pgclc_relative[,5],decreasing = T),1]



setdiff(study_gene,esc_relative[1:500])
intersect(study_gene,esc_relative[1:500])
setdiff(study_gene,tsc_relative[1:500])
intersect(study_gene,tsc_relative[1:100])
setdiff(study_gene,pgclc_relative[1:1000])
intersect(study_gene,pgclc_relative[1:1000])


gene_set=list(ESC=esc_relative,
              TSC=tsc_relative,
              PGCLC=pgclc_relative,
              #Primed=ETP[ETP$Cluster=="Primed",1],
              episc=episc$FeatureName,
              epilc_d2=epilc_d2$FeatureName)

gene_set=list(#CHIR_ac_1.5=chir_ac_1.5$X,
              CHIR_ac_2=chir_ac_2$X)



gene_set=list(BMP=signal[signal$signaling=="sig_BMP","gene"],
              FGF=signal[signal$signaling=="sig_Fgf","gene"],
              Hippo=signal[signal$signaling=="sig_Hippo","gene"],
              Nodal=signal[signal$signaling=="sig_Nodal","gene"],
              Notch=signal[signal$signaling=="sig_Notch","gene"],
              Wnt=signal[signal$signaling=="sig_Wnt","gene"],
              Hedgehog=signal[signal$signaling=="sig_Hedgehog","gene"]
              )

gene_set=list(k5_ac_2_bind=k5_ac_2_bind)

gene_set=list(k5_ac_esc=intersect(k5_ac_2_bind,esc_relative),
              k5_ac_tsc=intersect(k5_ac_2_bind,tsc_relative))

epilc_d2=read.csv("Graph-Based Features.csv")
epilc_d2=epilc_d2[epilc_d2$Cluster.4.Log2.Fold.Change>1,]
epilc_d2=epilc_d2[epilc_d2$Cluster.4.P.Value<0.01,]

episc=read.csv("Graph-Based Features_episc.csv")
episc=episc[episc$Cluster.2.Log2.Fold.Change>1,]
episc=episc[episc$Cluster.2.P.Value<0.01,]

gene_set=list(epilc_d2=epilc_d2$FeatureName)

seurat_D0_5=monocle
seurat_D0_5=seurat_D0_5[,!colnames(seurat_D0_5)%in%c("AGGCCGTAGGCCCGTT-2","ATGAGGGCAGCTCGAC-2","CCTACCACACAACGCC-2","CGTGTCTAGTAGTGCG-2","GCGGGTTCATCACAAC-2","GGTGAAGTCGGAGGTA-2",
                                                     "AGGGATGCATGGTTGT-5","CCTATTAAGGCAGGTT-5","CGGACACAGAATTGTG-5","CTTAGGAGTGTGCCTG-5","GTCGGGTCATATACCG-5","TGACAACAGCTAGTCT-5")]
#sub_seurat_merged=seurat_merged[,seurat_merged@meta.data$sample%in%c("D3_K5")]
dg <- as(as.matrix(seurat_D0_5@assays$RNA$data),"sparseMatrix")
cells_rankings <- AUCell_buildRankings(dg,plotStats = F,verbose = F)
class(dg)
dg=NULL
gc()
#saveRDS(cells_rankings,"pancreas_sub_cells_rankings.rds")

#cells_rankings=readRDS("pancreas_sub_cells_rankings.rds")
cells_AUC <- AUCell_calcAUC(gene_set, cells_rankings,aucMaxRank = 3000, nCores=10)


seurat_D0_5$EpiLC_signature  <- as.numeric(getAUC(cells_AUC)["epilc_d2", ])
seurat_D0_5$EpiSC_signature  <- as.numeric(getAUC(cells_AUC)["episc", ])

seurat_D0_5$ESC_AUC  <- as.numeric(getAUC(cells_AUC)["ESC", ])
seurat_D0_5$TSC_AUC  <- as.numeric(getAUC(cells_AUC)["TSC", ])
seurat_D0_5$PGCLC_AUC  <- as.numeric(getAUC(cells_AUC)["PGCLC", ])
seurat_D0_5$Primed_AUC  <- as.numeric(getAUC(cells_AUC)["Primed", ])

seurat$CHIR_ac_1.5_AUC  <- as.numeric(getAUC(cells_AUC)["CHIR_ac_1.5", ])

pancreas_sub$CHIR_ac_2_AUC  <- as.numeric(getAUC(cells_AUC)["CHIR_ac_2", ])


seurat_D0_5$Bmp_sig_AUC  <- as.numeric(getAUC(cells_AUC)["BMP", ])
seurat_D0_5$Fgf_sig_AUC  <- as.numeric(getAUC(cells_AUC)["FGF", ])
seurat_D0_5$Hippo_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Hippo", ])
seurat_D0_5$Nodal_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Nodal", ])
seurat_D0_5$Notch_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Notch", ])
seurat_D0_5$Wnt_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Wnt", ])
seurat_D0_5$Hedgehog_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Hedgehog", ])

seurat_D0_5$k5_ac_2_bind_AUC  <- as.numeric(getAUC(cells_AUC)["k5_ac_2_bind", ])
seurat_D0_5$k5_ac_esc_2_bind_AUC  <- as.numeric(getAUC(cells_AUC)["k5_ac_esc", ])
seurat_D0_5$k5_ac_tsc_2_bind_AUC  <- as.numeric(getAUC(cells_AUC)["k5_ac_tsc", ])

gc()

FeaturePlot(seurat,reduction = "tsne",features = c("Prdm1"),ncol=2,max.cutoff = "q90")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

DotPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,24,41)],features = c("k5_ac_2_bind_AUC"),group.by = "leiden")
DotPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,24,41)],features = c("Klf5"),group.by = "leiden")

VlnPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,24,41)],features = c("k5_ac_2_bind_AUC"),group.by = "leiden")
VlnPlot(seurat_D0_5[,seurat_D0_5$leiden%in%c(40,24,0,5,15,9,16,24,41)],features = c("Klf5"),group.by = "leiden")


sub=seurat_D0_5[,seurat_D0_5$leiden%in%c("EpiSC","EpiLC",6,4,23,35,32,34,37)]
sub$leiden=factor(sub$leiden,levels = c("EpiSC","EpiLC",6,4,32,23,35,34,37))
VlnPlot(sub,features = c("EpiLC_signature","Dnmt3b","Dnmt3a","Pou3f1"),ncol = 4,
        group.by = "leiden",pt.size = 0)
VlnPlot(sub,features = c("PGCLC_AUC","Prdm1","Dppa3","Nanos3"),ncol = 4,
        #features = c("PGCLC_AUC","Tmem54","Dnmt3b","Dnmt3a","Pou3f1"),
        group.by = "leiden",pt.size = 0)
VlnPlot(sub,features = c("ESC_AUC","Klf4","Tbx3","Fgf4"),ncol = 4,
        #features = c("PGCLC_AUC","Tmem54","Dnmt3b","Dnmt3a","Pou3f1"),
        group.by = "leiden",pt.size = 0)
VlnPlot(sub,features = c("EpiSC_signature","Fgf8","Fgf5","Lefty1"),ncol = 4,
        #features = c("PGCLC_AUC","Tmem54","Dnmt3b","Dnmt3a","Pou3f1"),
        group.by = "leiden",pt.size = 0)

Stacked_VlnPlot(sub, c("PGCLC_AUC","ESC_AUC","EpiLC_signature","EpiSC_signature"),x_lab_rotate = T,y.max = 0.3,
                # colors_use=c("#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA","#CCCCCC"),
                #colors_use=c("#FFD700","#DE2D26","#FF8C00","#90EE90","#87CEFA"),
                group.by = "leiden")


VlnPlot(seurat_D0_5,features = c("EpiLC_signature"),pt.size = 0,group.by = "leiden")

VlnPlot(seurat_D0_5,features = c(formative[1:12]),pt.size = 0,group.by = "leiden")





saveRDS(seurat,"K5_D0-5.rds")


###scale_AUC####


AUC_scale=data.frame(ESC_scale=monocle$ESC_AUC,
                     TSC_scale=monocle$TSC_AUC,
                     PGCLC_scale=monocle_annoed$PGCLC_AUC,
                     Primed_scale=monocle_annoed$Primed_AUC,
                     CHIR_ac_2_scale=monocle_annoed$CHIR_ac_2_AUC)
AUC_scale=scale(AUC_scale)
monocle_annoed@meta.data=cbind( monocle_annoed@meta.data,AUC_scale)


monocle$Bmp_sig_AUC  <- as.numeric(getAUC(cells_AUC)["BMP", ])
monocle$Fgf_sig_AUC  <- as.numeric(getAUC(cells_AUC)["FGF", ])
monocle$Hippo_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Hippo", ])
monocle$Nodal_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Nodal", ])
monocle$Notch_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Notch", ])
monocle$Wnt_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Wnt", ])
monocle$Hedgehog_sig_AUC  <- as.numeric(getAUC(cells_AUC)["Hedgehog", ])


AUC_scale=data.frame(Bmp_sig_scale=monocle$Bmp_sig_AUC,
                     Fgf_sig_scale=monocle$Fgf_sig_AUC,
                     Hippo_sig_scale=monocle$Hippo_sig_AUC,
                     Nodal_sig_scale=monocle$Nodal_sig_AUC,
                     Notch_sig_scale=monocle$Notch_sig_AUC,
                     Wnt_sig_scale=monocle$Wnt_sig_AUC,
                     Hedgehog_sig_scale=monocle$Hedgehog_sig_AUC)
AUC_scale=scale(AUC_scale)
monocle@meta.data=monocle@meta.data[,c(1:95,100)]
monocle@meta.data=cbind( monocle@meta.data,AUC_scale)


###do threshold####

cellauc=data.frame(ESC_AUC=seurat$ESC_AUC,
                   TSC_AUC=seurat$TSC_AUC,
                   PGCLC_AUC=seurat$PGCLC_AUC,
                   Primed_AUC=seurat$Primed_AUC)
write.table(cellauc,"D5_K245_cellauc.txt",sep = "\t",quote = F)
cellauc <- importAUCfromText(fileName = "D5_K245_cellauc.txt")
cellauc <- cellauc[onlyNonDuplicatedExtended(rownames(cellauc)),]

cellauc=AUCell_exploreThresholds(cellauc,plotHist = T,assignCells = T)

ESC_threshold=0.147
TSC_threshold=as.numeric(cellauc$TSC_AUC$aucThr$selected)
PGCLC_threshold=0.275

seurat@meta.data$ESC_AUC_define=ifelse(seurat@meta.data$ESC_AUC>=ESC_threshold,
                                       "ESC","Others")
seurat@meta.data$TSC_AUC_define=ifelse(seurat@meta.data$TSC_AUC>=TSC_threshold,
                                       "TSC","Others")
seurat@meta.data$PGCLC_AUC_define=ifelse(seurat@meta.data$PGCLC_AUC>=PGCLC_threshold,
                                       "PGCLC","Others")

###statistic celltype####

VlnPlot(seurat,pt.size = 0.1, features = c("TSC_AUC"))



meta.data=seurat@meta.data

ESC_meta.data=meta.data[meta.data$ESC_AUC_define=="ESC",]
table(ESC_meta.data$sample)/table(meta.data$sample)[1:3]*100
PGCLC_meta.data=meta.data[meta.data$PGCLC_AUC_define=="PGCLC",]
table(PGCLC_meta.data$sample)/table(meta.data$sample)[1:3]
ESC_meta.data=meta.data[meta.data$ESC_AUC_define=="ESC",]
table(ESC_meta.data$sample)/table(meta.data$sample)[1:3]

st=meta.data[meta.data$seurat_clusters==17,]
table(st$sample)/table(seurat$sample)[1:3]











#study_gene####

Klfs=c("Klf2","Klf5")
ESC_Marker=c("Nanog","Tfcp2l1","Zfp42","Klf4","Tbx3")
PGC_Marker=c("Dppa3","Prdm1","Nanos3","Sox15","Id3")
TSC_and_ESC=c("Esrrb","Socs3","Tcf7l1","Zbtb18")
ESC_and_PGC=c("Pou5f1","Dppa5a","Nr5a2","Prdm14")
PGC_and_TSC=c("Tfap2c","Rhox6","Rhox9","Msx2","Msx1")
TSC_Marker=c("Gata3","Cdx2","Krt7","Gata2","Rin3","Plet1","Hand1","Elf5")

ETS=c("Klf5","Zfp42","Esrrb")
TSC_and_ESC=c("Tead4","Sox21")
ESC_and_PGC=c("Nanog","Pou5f1","Klf2","Nr5a2","Dppa5a")
PGC_and_TSC=c("Tfap2c","Rhox6","Rhox9","Msx2","Msx1")
ESC_Marker=c("Klf4","Tfcp2l1","Tbx3")
TSC_Marker=c("Gata3","Cdx2","Ddah1","Plet1","Elf5")
PGC_Marker=c("Prdm1","Dppa3","Nanos3","Id3","Sox15")

study_list=list(ETS=ETS,
                TSC_and_ESC=TSC_and_ESC,
                ESC_and_PGC=ESC_and_PGC,
                PGC_and_TSC=PGC_and_TSC,
                ESC=ESC_Marker,
                TSC=TSC_Marker,
                PGC=PGC_Marker)

study_gene=as.character(unlist(study_list))


en=c("Gdf3","Nckap1","Nodal","Ext1","Arc","Hhex","Gsc","Hnf1b","Sox17","Gata4",
     "Gata6","Foxa2","Cfc1","Sox7","Gata5",
     "Claudin6","Dab2","Foxa2","Fxyd3","Tm4sf20","Ttf1","Amn","Cxcr4","Emp2","Sdc4" #K4,K5,KD NC文章
)
me=c("Srf","Bmpr2","Eomes","T","Msx2","Wnt3","Cebpa","Bmp6","Cdx4","Mixl1","Runx1","Tal1",
     "Flk1","Myt6","Myod1","Spef2","Opn4","Tbx20","Vegfr1","Bmp2","Hand1" #K4,K5,KD NC文章
)
ex=c("Tcf15","Pax6","Dsp","Nes","Otx2","Gja1","Krt18","Tgm3","Col5a1","Sox18",
     "Runx3","Fgf5","Gfap",
     "Msi2","Cdh2","Neurod1"#K4,K5,KD NC文章
)
diff=c(en,me,ex)

formative=c("Utf1","Pou3f1","Fgf5","Sox3","Lef1","Dnmt3a","Dnmt3b","Dnmt3l",
            "Tead2","Otx2","Sox12","Pim2","Sall2","Hes6","Zic2","Zic3","Zic5","Lin28a",
            "Lin28b","Zyx","Arid2","Foxp1","Cebpb","Mapk12","Fzd2","Rspo1","Smo")


#EpiSCs DEG analysis####


seurat_D0_5=readRDS("K5_D0-5.rds")

Idents(seurat_D0_5)="seurat_cluster"
deg_res=FindMarkers(seurat_D0_5,ident.1 = 13,ident.2 = 4,
                      logfc.threshold = 0,only.pos = F,
                      min.pct = 0.1)

Idents(monocle)="seurat_cluster"
table(monocle$seurat_cluster)
deg_res=FindMarkers(monocle,ident.1 = 13,ident.2 = 4,
                    logfc.threshold = 0,only.pos = F,
                    min.pct = 0.1)

dim(deg_res)
deg_res[rownames(deg_res)%in%c(study_gene),]


EpiSC=seurat_D0_5[,seurat_D0_5$sample=="EpiSCs"]
EpiSC$leiden=EpiSC$seurat_cluster
EpiSC=EpiSC[,EpiSC$leiden%in%c(4,13)]
EpiSC$leiden=gsub(4,"Epi-4",EpiSC$leiden)
EpiSC$leiden=gsub(13,"Epi-13",EpiSC$leiden)
table(EpiSC$leiden)

VlnPlot(seurat_D0_5, features = c("Klf5", "Nanog"),pt.size = 0)

Stacked_VlnPlot(EpiSC, c(study_gene,"Klf5"),x_lab_rotate = T,#y.max = 2,
                # colors_use=c("#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA","#CCCCCC"),
                #colors_use=c("#FFD700","#DE2D26","#FF8C00","#90EE90","#87CEFA"),
                group.by = "leiden")
p <- DotPlot(EpiSC, features=study_gene,
             scale = F,
             group.by='leiden',
             cols =  c("RdYlBu"))
p+ggplot2::coord_flip()

library(ggplot2)
library(ggrepel)
target_genes <- c("Gata3","Cdx2","Fgfr2","Hand1","EIf5","Tfap2c",
                  "Rhox6","Rhox9","Msx2","Prdm1","Dppa3","Sox15","Dmrt1",
                  "Nanos3","Dazl","Pou5f1","Nanog","Dppa5a","Nr5a2",
                  "Klf4","Tbx3","Zfp42",PGC_all,TSC_all)
PGC_all=c("Dppa3", "Prdm1", "Tfap2c", "Nanos3","Id3", "Rhox6", "Irf1", "Sox15", "Rhox9", "Zic3", "Kit", 
          "Rbm47", "Phf19", "Msx1",  "Msx2","Dmrt1", "Dnd1", "Prdm14", "Dazl", "Ddx4")

TSC_all=c("Cdx2","Gata3","Tfap2c","Rhox5","Fgfr2","Gata2","Krt18","Ddah1","Rin3","Plet1","Hand1","Elf5",
          "Tead4","Eomes","Ets2","Yap1","Erf","Ly6a")

ESC_all=c("Klf4","Tbx3","Zfp42","Klf2","Tfcp2l1","Dppa5a","Esrrb","Prdm14","Nanog","Nr5a2",
          "Gbx2","Pbx1","Tcf15","Tcf3","Tgif1")

# 确保deg_res行名包含这些基因
# 如果deg_res行名不匹配，可能需要调整大小写或格式

deg_res$gene <- rownames(deg_res)
deg_res$log_pval <- -log10(deg_res$p_val_adj)
deg_res$log_pval[is.infinite(deg_res$log_pval)]=300

# 标记目标基因（即谱系标记基因）
deg_res$is_target <- deg_res$gene %in% target_genes

# 创建更详细的分类
deg_res$category <- "Not significant"
deg_res$category[deg_res$p_val_adj < 0.00001 & deg_res$avg_log2FC > log2(1.5)] <- "Cluster 13 up"
deg_res$category[deg_res$p_val_adj < 0.00001 & deg_res$avg_log2FC < -log2(1.5)] <- "Cluster 4 up"
# 谱系标记基因覆盖其他分类
deg_res$category[deg_res$is_target] <- "Expressed 3 lineage markers"

# 确保分类因子顺序正确
deg_res$category <- factor(deg_res$category, 
                           levels = c("Not significant", "Cluster 13 up", "Cluster 4 up", "Expressed 3 lineage markers"))

# 统计基因数量（不包括谱系标记基因）
cluster13_up_count <- sum(deg_res$category == "Cluster 13 up", na.rm = TRUE)
cluster4_up_count <- sum(deg_res$category == "Cluster 4 up", na.rm = TRUE)
lineage_count <- sum(deg_res$category == "Expressed 3 lineage markers", na.rm = TRUE)

table(deg_res$category)
# 绘制点图
ggplot(deg_res, aes(x = avg_log2FC, y = log_pval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  # 特别突出谱系标记基因
  geom_point(data = subset(deg_res, is_target),
             size = 3, shape = 17, alpha = 0.9) +  # 菱形形状，更大更明显
  geom_text_repel(data = subset(deg_res, is_target),
                  aes(label = gene),
                  size = 3.5,
                  color = "red",
                  fontface = "bold",
                  box.padding = 0.8,    # 增加框内边距
                  point.padding = 0.5,  # 增加点边距
                  min.segment.length = 0.2,  # 减少最小线段长度
                  force = 2,            # 增加排斥力
                  force_pull = 1,       # 增加拉力
                  direction = "both",   # 允许双向移动
                  nudge_y = -5,        # 向下偏移，利用负Y轴空间
                  max.time = 10,        # 增加计算时间
                  max.overlaps = Inf) +
  # 添加阈值线
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.00001), linetype = "dashed", alpha = 0.5) +
  # 设置颜色，突出谱系标记
  scale_color_manual(values = c("Not significant" = "gray",
                                "Cluster 13 up" = "#ff7f0e",  # 橙色
                                "Cluster 4 up" = "#1f77b4",   # 蓝色
                                "Expressed 3 lineage markers" = "red"), # 红色突出
                     name = "Gene Category") +
  # 添加基因数量标注
  annotate("text", x = max(deg_res$avg_log2FC, na.rm = TRUE) * 0.8, 
           y = 250,
           label = paste("Cluster 13 up:", cluster13_up_count),
           color = "#ff7f0e", size = 4, fontface = "bold") +
  annotate("text", x = min(deg_res$avg_log2FC, na.rm = TRUE) * 0.8, 
           y = 250,
           label = paste("Cluster 4 up:", cluster4_up_count),
           color = "#1f77b4", size = 4, fontface = "bold") +
  # 设置Y轴范围，提供更多标签空间
  scale_y_continuous(limits = c(-50, NA)) +
  theme_bw() +
  labs(title = "EpiSCs clusters 13 vs 4 DEG",
       x = "Log2FC",
       y = "-Log10(P.adj)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.title = element_text(face = "bold"))


#Make gastrulation in vivo data####

setwd("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/in_vivo_data/gastrulation_GSE86586/")

seurat=readRDS()


p <- DotPlot(seurat_fin, features=c("Klf5","Cdx2"),
             scale.min = 0,
             scale.max = 50,
             scale = F,
             group.by='Celltype',
             cols =  c("RdYlBu"))
p+ggplot2::coord_flip()
p_data=p$data
p_anno=data.frame(id=levels(p_data$id),
                  label=factor(c(rep("Zygote-Morula",8),rep("1st Segregation",2),rep("Post-implant",1)),
                               levels=c("Zygote-Morula","1st Segregation","Post-implant"))
)
p_data=plyr::join(p_data,p_anno)
p_data$features.plot=factor(p_data$features.plot,levels =
                              rev(c("Klf1","Klf15","Klf12","Klf14",
                                    "Klf6","Klf17","Klf5","Klf4","Klf2","Klf7",
                                    "Klf9","Klf3","Klf13","Klf8","Klf10","Klf11",
                                    "Klf16")))
p_data
scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size, #该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
limit_max=2
limit_max=2.5

name="Figure.1.A_Stringency_Single_cell_data_new_ver1"
pdf(paste(name,".pdf",sep=""),width=5.2,height =5)

p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
#该函数可以完全替代scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  scale.func(range = c(0,6.5))+#点大小
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(0,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))
dev.off()
dev.new()



#K5 D0-5 and in vivo data####

###D0-5 tsne####
#sample_palette_col=c("grey","#EF3B2C","#4DAF4A","#FFB6C1","#377EB8","#FFA500")
#sample_palette_col=c("grey","#BC80BD","#4DAF4A","#FFA500","#377EB8","#FB8072")
sample_palette_col=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#EF3B2C")


seurat_D0_5=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/K5_D0-5.rds")
#saveRDS(seurat_D0_5,"K5_D0-5.rds")


seurat_cluster=readRDS("seurat_clusters.rds")
seurat_D0_5$seurat_cluster=seurat_cluster

CellDimPlot(seurat_D0_5,reduction = "tsne", group.by =c("sample"),raster = F,
            label = F,label.fg = "white",label_insitu = T,palcolor = sample_palette_col )


DimPlot(seurat_D0_5,reduction = "tsne", group.by =c("seurat_cluster"),raster = F,label.size = 6,
          label =T,cols=colorRampPalette(colors=palette_list$Set2)(32))
#10x10

FeatureDimPlot(seurat_D0_5,raster =F, pt.size = .1,
               features = c("ESC_AUC"),lower_cutoff = 0.1,upper_cutoff = 0.18,
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("TSC_AUC"),lower_cutoff = 0.13,upper_cutoff = 0.19,
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("PGCLC_AUC"),lower_cutoff =0.25,upper_cutoff = 0.31, 
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("Primed_AUC"),lower_cutoff = 0.15,upper_cutoff = 0.2,  
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("Bmp_sig_AUC"),lower_quantile = 0.01,upper_quantile = 0.99,
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("Id2"),lower_quantile = 0.01,upper_quantile = 0.99,  
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("Nanog"),#lower_cutoff = 0.15,upper_cutoff = 0.2,  
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("Klf5"),#lower_cutoff = 0.15,upper_cutoff = 0.2,  
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = c("Sox2"),#lower_cutoff = 0.15,upper_cutoff = 0.2,  
               reduction = "tsne", theme_use = "theme_blank")
#10x10


bmp_sig=signal[signal$signaling=="sig_BMP","gene"]

marker=c("Klf4","Tbx3","Nanog","Zfp42",
  "Cdx2","Gata3","Tead4","Elf5","Plet1","Ascl2",
  "Dppa3","Prdm1","Nanos3",
  intersect(k5_bind_ac,bmp_sig))

for(i in marker){
color = c(colorRampPalette(colors = "gray")((40)),
          colorRampPalette(colors =brewer.pal(n = 9, name = "Reds"))((60)))

FeatureDimPlot(seurat_D0_5,raster = F, pt.size = .1,
               features = i,palcolor = color,upper_quantile = 0.999,
               reduction = "tsne", theme_use = "theme_blank")
  name=paste0("K5_D0-5_tsne_",i,".pdf")
  ggsave(
    filename = name, # 保存的文件名称。通过后缀来决定生成什么格式的图片
    width = 10,             # 宽
    height = 10,            # 高
    units = "in",          # 单位
    dpi = 500)              # 分辨率DPI
}



###D0-5 clusters dotplot####

study_gene=rev(c("ESC_AUC","TSC_AUC","PGCLC_AUC","Primed_AUC",
             "Klf4","Tbx3","Fgf4","Zfp42",
             "Cdx2","Elf5","Gata3","Gata2","Krt7","Plet1",
             "Ascl2","Hand1","Gcm1","Krt18","Krt8","Klf6","Dsp","Cgn",
             "Prdm1","Nanos3","Dppa3","Ifitm1","Kit",
             "Pou5f1","Nanog","Klf2","Dppa5a","Nr5a2","Prdm14",
             "Esrrb","Tfap2c","Rhox6","Rhox9","Sox15","Tead4"))

study_gene=rev(c("ESC_AUC","PGCLC_AUC","TSC_AUC","Primed_AUC",
                 "Klf4","Tbx3","Fgf4",
                 "Prdm1","Nanos3","Kit","Dppa3",
                 "Pou5f1","Nanog","Klf2","Esrrb",
                 "Tfap2c","Rhox6","Rhox9",
                 "Tead4","Gata3","Cdx2","Elf5",
                 "Krt7","Plet1","Ascl2","Hand1","Krt18","Krt8","Klf6","Dsp","Cgn"))

color = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"   ))((50)),
          colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026" ))((50)))

p <- DotPlot(seurat_D0_5, features=study_gene,
             scale.min = 0,
             scale.max = 100,
             scale = T,
             group.by='seurat_cluster',
             cols =  c("RdYlBu"))
p+ggplot2::coord_flip()
p_data=p$data
# p_anno=data.frame(id=levels(p_data$id),
#                   label=factor(c(rep("Zygote-Morula",8),rep("1st Segregation",2),rep("Post-implant",1)),
#                                levels=c("Zygote-Morula","1st Segregation","Post-implant"))
# )
# p_data=plyr::join(p_data,p_anno)

p_data$features.plot=factor(p_data$features.plot,levels = study_gene)
p_data
scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size, #该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
limit_max=2
limit_max=3
limit_min=-1

p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
p_data[p_data$avg.exp.scaled<=limit_min,"avg.exp.scaled"]=as.numeric(limit_min)
#该函数可以完全替代scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  scale.func(range = c(0,6.5))+#点大小
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(limit_min,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  #facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))




###D0-5 ETP clusters dotplot####

select_etp=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/K5_0-5D_select_etp_cl.rds")
#sub_select_etp=select_etp[,select_etp$Celltype=="iPSC"]

anno_seurat=select_etp

min.cells=dim(anno_seurat@assays$RNA$counts)[2]*0.01

select_etp=FindNeighbors(select_etp,dims = 1:16)
select_etp=FindClusters(select_etp,resolution = 1,random.seed =2)
p1=DimPlot(select_etp,reduction="tsne",group.by = "Celltype",label = T,label.size =5,pt.size = 3,
           cols=colorRampPalette(colors = brewer.pal(n=7,"Set2"))(24))
p2=DimPlot(select_etp,reduction="tsne",group.by = "seurat_clusters",label = T,pt.size = 3,label.size =5,
           cols=colorRampPalette(colors = brewer.pal(n=8,"Set2"))(length(levels(factor(select_etp$seurat_clusters)))))

ggpubr::ggarrange(p1,p2,ncol = 2,nrow =1)


anno_seurat <- RunKNNPredict(srt_query = anno_seurat, bulk_ref = ETP_ref,filter_lowfreq = min.cells)
#anno_seurat <- RunKNNPredict(srt_query = anno_seurat, bulk_ref = ref_scMCA,filter_lowfreq = min.cells)

anno_seurat$Celltype=gsub("iPSC","ESC",anno_seurat$Celltype)
anno_seurat$Celltype=gsub("TSCLC","TSC",anno_seurat$Celltype)

anno_seurat$etp_cl=gsub("28","28-TSC",anno_seurat$etp_cl)
anno_seurat$etp_cl=gsub("26","26-TSC",anno_seurat$etp_cl)
anno_seurat$etp_cl=gsub("27","27-TSC",anno_seurat$etp_cl)
anno_seurat$etp_cl=gsub("20","20-PGCLC",anno_seurat$etp_cl)
anno_seurat$etp_cl=gsub("18","18-ESC",anno_seurat$etp_cl)

table(anno_seurat$etp_cl)

meta=anno_seurat@meta.data
meta$anno_celltype="Others"
meta[meta$Celltype==meta$KNNPredict_classification,"anno_celltype"]=meta[meta$Celltype==meta$KNNPredict_classification,"Celltype"]

table(meta$anno_celltype)
anno_seurat@meta.data=meta
anno_seurat$anno_celltype=factor(anno_seurat$anno_celltype,levels = c("ESC","PGCLC","TSC","Others"))

CellDimPlot(anno_seurat,reduction="tsne",group.by = c("anno_celltype"),
            pt.size = 1,label.fg = "white",palcolor =  c("#66C2A5" ,"#6699CC","#CC6677","grey"))

CellDimPlot(anno_seurat,reduction="tsne",group.by = c("etp_cl"),label = T,label_insitu = T,pt.size = 1,palette = "Set2")



p <- DotPlot(select_etp, features=study_gene,
             scale.min = 0,
             scale.max = 100,
             scale = T,
             group.by='seurat_clusters',
             cols =  c("RdYlBu"))
p+ggplot2::coord_flip()

p_data=p$data
p_data$id=factor(p_data$id,levels = c(10,8,3,7,
                                       6,13,20,
                                       17,1,0,2,
                                       16,11,22,12,9,15,
                                       14,19,23,4,5,18,21))
p_anno=data.frame(id=levels(p_data$id),
                   label=factor(c(rep("iPSC",4),
                                  rep("intermediate",3),
                                  rep("PGCLC",4),rep("TSCLC",8),
                                  rep("Others",5)),
                                levels=c("iPSC","intermediate","PGCLC","TSCLC","Others")))
p_data=plyr::join(p_data,p_anno)


p_data$features.plot=factor(p_data$features.plot,levels = study_gene)
p_anno=data.frame(features.plot=rev(levels(p_data$features.plot)),
                  marker=factor(c(rep("AUC",4),
                                 rep("ESC",3),
                                 rep("PGC",3),rep("ESC&PGC",4),
                                 rep("PGC&TSC",3),rep("TSC",13)),
                               levels=c("AUC","ESC","PGC","ESC&PGC","PGC&TSC","TSC")))

p_data=plyr::join(p_data,p_anno)
#p_data
scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size, #该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
limit_max=2
limit_max=2.5
limit_min=-2.5

p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
p_data[p_data$avg.exp.scaled<=limit_min,"avg.exp.scaled"]=as.numeric(limit_min)
#该函数可以完全替代scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  scale.func(range = c(0,6.5))+#点大小
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(limit_min,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  #facet_grid(rows = vars(marker),scales = "free",space = "free")+theme_classic()+
  facet_grid(cols=vars(label),rows = vars(marker),scales = "free",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))


saveRDS(anno_seurat,"K5_0-5D_select_etp_cl.rds")


###D0-5 ETP marker gene exp dynamic####

select_etp=readRDS("K5_0-5D_select_etp_cl.rds")

npsc=select_etp[,select_etp$Celltype=="ESC"]
pgc=select_etp[,select_etp$Celltype=="PGCLC"]
tsc26=select_etp[,select_etp$etp_cl=="26-TSC"]
tsc27=select_etp[,select_etp$etp_cl=="27-TSC"]
tsc28=select_etp[,select_etp$etp_cl=="28-TSC"]
tsc=select_etp[,select_etp$Celltype=="TSC"]

use_ftr=c("Klf5","Eomes","Tead4","Gata3","Id2",
          "Krt7","Plet1","Tfap2c","Ets2",
          "Krt18","Ascl2","Krt19",
          "Klf4","Tbx3","Nanog",
          "Dppa3","Prdm1","Nanos3")


use_ftr=c("Klf4","Tbx3","Nanog","Dppa5a")
use_ftr=rev(c("Prdm1","Dppa3","Kit","Prdm14"))

use_ftr=rev(c("Tead4","Ly6a","Eomes","Ddah1","Gata3","Ets2",
          "Ascl2","Krt19","Hand1"))


p=DotPlot(tsc,features =  use_ftr,
          group.by = "sample")
p+ggplot2::coord_flip()
p_data=p$data
# p_data$id=factor(p_data$id,levels=c("Tang_Oocyte","Tang_Zygote","Tang_2cell","Tang_4cell","Tang_8cell",
#                                     "Tang_Morula","Tang_LateBlastocyst","Tang_hESCpassage#0","Tang_hESCpassage#10"))
# p_anno=data.frame(id=levels(p_data$id),
#                   label=factor(c(rep("D1",1),rep("Reprogramming",1),
#                                  rep("TSCLC",1),rep("Non-Reprogramming",5)),
#                                levels=c("D1","Reprogramming","TSCLC","Non-Reprogramming")))
# 
# p_data=plyr::join(p_data,p_anno)

p_data$features.plot=factor(p_data$features.plot,levels =
                              c(use_ftr))
#p_data
scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size ,#该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

name="Figure.1.A_Stringency_Single_cell_data_new_ver1"
#pdf(paste(name,".pdf",sep=""),width=5.2,height =5)

color = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"   ))((50)),
          colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026" ))((50)))
#p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
#该函数可以完全替代scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  #scale.func(range = c(0,6.5))+#点大小
  scale_size(limits = c(0,100),range = c(0,6.5))+
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(-2,2),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  #facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))





###D0-5 ETP statistic####

seurat_D0_5$invivo_celltype="Others"
meta=seurat_D0_5@meta.data
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(22)]),"invivo_celltype"]="TSCLC-subTE3.5"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(19)]),"invivo_celltype"]="TSCLC-subTE4.5"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(20)]),"invivo_celltype"]="TSCLC-subExE"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(2,4)]),"invivo_celltype"]="PGCLC-sub"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(0,3)]),"invivo_celltype"]="iPSC-sub"
table(meta$invivo_celltype)
seurat_D0_5@meta.data=meta

seurat_D0_5=readRDS("K5_0-5D_select_etp_cl.rds")

Bar_Line_plot=
  function(anno_seurat,celltype,limits,title=celltype,
           sample_palette_col=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#EF3B2C")){
    finmat=data.frame(Var1=paste0("Day",seq(0,5)),Freq=0)
    
    xx=anno_seurat[,anno_seurat$etp_cl%in%celltype]
    xx=as.data.frame(table(xx$sample))
    xx$Var1=gsub("D","Day",xx$Var1)
    xx$Var1=strsplit2(xx$Var1,"_")[,1]
    
    finmat[which(finmat$Var1%in%xx$Var1),"Freq"]=xx[which(xx$Var1%in%finmat$Var1),"Freq"]
    
    finmat$Percentage=finmat$Freq/c(9470,11898,15769,20001,21329,38946)
    #finmat$Percentage=finmat$Freq/c(9470,11898,15769,20001,21329,38946)*3e4
    
    p=ggplot(finmat)+
      geom_bar(aes(x=Var1,y=Percentage,fill=cl),stat = "identity",position = "stack")+
      scale_fill_manual(values = sample_palette_col)+
      #geom_line(aes(x=Var1,y=Percentage,group = 1))+
      #geom_point(aes(x=Var1,y=Percentage,group = 1))+
      scale_y_continuous(limits = limits,labels=scales::percent)+
      labs(x=NULL,y="Percentage %",title=title)+
      #geom_text(aes(x=Var1,y=Percentage,label = paste0("n=",Freq)),size=2, position = position_dodge(0.9), vjust = -0.8)+
      #guides(fill = FALSE)+
      theme(panel.grid = element_blank())+
      theme_classic()+
      theme(axis.text.x = element_text(angle=90,hjust=1)
            #axis.text.y = element_text(size = 5)
      )
    return(p)
    
  }


Bar_Line_plot(seurat_D0_5,"18-ESC",c(0,0.06))

Bar_Line_plot(seurat_D0_5,"20-PGCLC",c(0,0.05))

Bar_Line_plot(seurat_D0_5,c("26-TSC","27-TSC","28-TSC"),c(0,0.1),"26+27+28-TSC",
              sample_palette_col=c("#DE2D26","#FF8C00","#FFD700"))


##Celltype annotation####

data("ref_scMCA")

ETP_ref=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/ETP_log2_TPM.rds")
ETP_ref[ETP_ref<0]=0
ETP_ref=as(as.matrix(ETP_ref),"dgCMatrix")

# 
min.cells=dim(obj@assays$RNA$counts)[2]*0.01


obj <- RunKNNPredict(srt_query = obj, bulk_ref = ref_scMCA,filter_lowfreq = min.cells)

table(obj$KNNPredict_classification)

saveRDS(obj$KNNPredict_classification,"K5_D0-5_celltype_anno.rds")

obj$KNNPredict_classification=readRDS("K5_D0-5_celltype_anno.rds")

meta2=obj@meta.data

meta2=Btable_insert_into_Atable(meta2,meta[,c("orig.ident","KNNPredict_classification")])


meta2[rownames(meta2)==rownames(meta),"KNNPredict_classification"]=meta[rownames(meta)==rownames(meta2),"KNNPredict_classification"]
dim(meta2)

table(meta2[,40])
meta2$celltype=meta2$KNNPredict_classification
for(i in unique(meta2[,39][!is.na(meta2[,39])])){
  meta2[rownames(meta2)%in%rownames(meta2[grepl(i,meta2[,39]),]),40]=i
}

for(i in 1:dim(meta2)[1]){
  if(meta2[i,"sample"]!="D5_K5"){
    meta2$celltype[i]="Others"
  }
}


seurat_D0_5$Celltype=readRDS("K5_D0-5_celltype_anno.rds")

seurat_D0_5$Celltype[seurat_D0_5$Celltype=="ES_Actb_high(Embryonic-Stem-Cell)"]="EpiSCs"
seurat_D0_5$Celltype[seurat_D0_5$Celltype=="ES_Nedd4_high(Embryonic-Stem-Cell)"]="Others"
seurat_D0_5$Celltype[seurat_D0_5$Celltype=="ES_Rps28_high(Embryonic-Stem-Cell)"]="EpiSCs"

seurat_D0_5$Celltype=strsplit2(seurat_D0_5$Celltype,"\\_")[,1]
seurat_D0_5$Celltype=strsplit2(seurat_D0_5$Celltype,"\\(")[,1]

seurat_D0_5$Celltype[seurat_D0_5$Celltype=="TS"]="TSC"
seurat_D0_5$Celltype=gsub("MEF","TSC",seurat_D0_5$Celltype)
table(seurat_D0_5$Celltype)

CellDimPlot(seurat_D0_5,reduction = "tsne", group.by =c("Celltype"),raster = F,
            label = F,label.fg = "white",label_insitu = T)

obj=seurat_D0_5
obj <- RunKNNPredict(srt_query = obj, bulk_ref = ETP_ref,filter_lowfreq = min.cells)


meta=obj@meta.data

meta2=seurat_D0_5@meta.data


meta2[meta2$Celltype%in%c("EpiSCs","unreliable","Others"),"Celltype"]=meta[meta2$Celltype%in%c("EpiSCs","Others"),"KNNPredict_classification"]

meta2 <- meta2 %>%
  mutate(Celltype = case_when(
    Celltype %in% c("EpiSCs", "unreliable", "Others") & 
      meta$KNNPredict_classification != "EpiSCs" & 
      !is.na(meta$KNNPredict_classification) ~ meta$KNNPredict_classification,
    TRUE ~ Celltype
  ))

for(i in 1:dim(meta2)[1]){
  if(meta2[i,"sample"]!="D5_K5"){
    meta2$Celltype[i]="Others"
  }
}



rows_to_replace <- rownames(meta2) %in% rownames(d5_meta)

# 进行替换
meta2[rows_to_replace, "Celltype"] <- d5_meta[rownames(meta2)[rows_to_replace], "celltype"]

seurat_D0_5$celltype=meta2$Celltype
seurat_D0_5$celltype=factor(seurat_D0_5$celltype,levels=c("Astrocytes","Dopaminergic amacrine neuron",
                                                          "Glomerular epithelial cell","Epithelial cell",
                                                          "GABAergic neuron","Ganglion cell",
                                                          "Immature neuron","Neuron","Radial glia cell",
                                                          "Stromal cell","ESC","EpiSCs","TSC","PGCLC",
                                                          "unreliable","Others"))
CellDimPlot(seurat_D0_5,reduction = "tsne", group.by =c("celltype"),raster = F,
            label = F,label.fg = "white",label_insitu = T,palcolor = c("#FDBF6F","#FFFF99","#FB9A99",
                                                                       "#B15928","pink","#191970","#E31A1C",
                                                                       "#B2DF8A","#CAB2D6","#6A3D9A",
                                                                       "#33A02C","#A6CEE3","#FF7F00",
                                                                       "#1F78B4","#666666","grey"))

DimPlot(seurat_D0_5,reduction = "tsne", group.by =c("celltype"),raster = F,label.size = 6,
        label =T,cols=c("#FDBF6F","#FFFF99","#FB9A99",
                        "#B15928","pink","#191970","#E31A1C",
                        "#B2DF8A","#CAB2D6","#6A3D9A",
                        "#33A02C","#A6CEE3","#FF7F00",
                        "#1F78B4","#666666","grey"))
#8x10


seurat_D5=seurat_D0_5[,seurat_D0_5$sample=="D5_K5"]
seurat_D5=seurat_D5[,seurat_D5$seurat_cluster%in%c(0,2,5,6,7,10,14,15,23,25)]

min.cells=dim(seurat_D5@assays$RNA$counts)[2]*0.005
seurat_D5=RunKNNPredict(srt_query = seurat_D5, bulk_ref = ref_scMCA,nn_method = "raw",
                        filter_lowfreq = min.cells)

CellDimPlot(seurat_D5,reduction = "tsne", group.by =c("KNNPredict_classification"),raster = F,
            label = F,label.fg = "white",label_insitu = T)


d5_meta=seurat_D5@meta.data

d5_meta <- d5_meta %>%
  mutate(celltype = case_when(
    celltype == "unreliable" ~ KNNPredict_classification,
    TRUE ~ celltype
  ))


d5_meta$celltype=strsplit2(d5_meta$celltype,"\\_")[,1]
d5_meta$celltype=strsplit2(d5_meta$celltype,"\\(")[,1]

d5_meta$celltype[d5_meta$celltype=="TS"]="TSC"
d5_meta$celltype[d5_meta$celltype=="Stromal Cell"]="Stromal cell"
d5_meta$celltype[d5_meta$celltype=="GABAergic neuron"]="Neuron"
# 获取cluster和cell type信息
cluster_info <- d5_meta$seurat_cluster
celltype_info <- d5_meta$celltype

# 创建数据框
plot_data <- data.frame(
  Cluster = cluster_info,
  CellType = celltype_info
)

# 计算每个cluster中cell type的比例
prop_data <- plot_data %>%
  group_by(Cluster, CellType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Cluster) %>%
  mutate(Proportion = Count / sum(Count))

my_colors <- c("#FDBF6F","#FFFF99","#FB9A99",
               "#B15928","pink","#191970","#E31A1C",
               "#B2DF8A","#CAB2D6","#6A3D9A",
               "#A6CEE3","#FF7F00","#666666")

prop_data$CellType=gsub("PGCLC","unreliable",prop_data$CellType)

prop_data$CellType=factor(prop_data$CellType,levels=c("Astrocytes","Dopaminergic amacrine neuron",
                                                      "Glomerular epithelial cell","Epithelial cell",
                                                      "GABAergic neuron","Ganglion cell",
                                                      "Immature neuron","Neuron","Radial glia cell",
                                                      "Stromal cell","ESC","EpiSCs","TSC","PGCLC",
                                                      "unreliable","Others"))
ggplot(prop_data, aes(x = Cluster, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Cell Type Proportion by Cluster",
       x = "Cluster",
       y = "Percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10))





###compare with EPSC####

AA=readRDS("monocle2_new.rds")
Idents(AA)="sample"
A=subset(AA,downsample=600)


A=CreateSeuratObject(counts =A@assays$RNA$counts,
                   meta.data =A@meta.data,project = "k5")
#A=NormalizeData(A,normalization.method = "RC")
#A@assays$RNA$counts=A@assays$RNA$data


EPSC=readRDS("GSE145609-LD.EPSCs-raw.Rds")
EPSC=as.data.frame(fread("GSE145609_SmartSeq2_RPKM_matrix.txt"))
EPSC$`GeneID/SampleID`=gsub("_","",EPSC$`GeneID/SampleID`)
EPSC$`GeneID/SampleID`=strsplit2(EPSC$`GeneID/SampleID`,split = "\\+")[,1]
EPSC=EPSC[!duplicated(EPSC$`GeneID/SampleID`),]
rownames(EPSC)=EPSC$`GeneID/SampleID`
dim(EPSC)

D.EPSC=EPSC[,grep("D.EPSC",colnames(EPSC))]
L.EPSC=EPSC[,grep("L.EPSC",colnames(EPSC))]
ESC=EPSC[,grep("_2i",colnames(EPSC))]

D.EPSC=CreateSeuratObject(counts =D.EPSC,project = "D.EPSC")
D.EPSC$sample="D.EPSC"
L.EPSC=CreateSeuratObject(counts =L.EPSC,project = "L.EPSC")
L.EPSC$sample="L.EPSC"
ESC=CreateSeuratObject(counts =ESC,project = "ESC")
ESC$sample="ESC"


D.EPSC$stim="EPSC"
L.EPSC$stim="EPSC"
ESC$stim="EPSC"
A$stim="K5"

merge=merge(A[,A$sample%in%c("EpiSCs","D1_K5","D2_K5","D3_K5")],c(D.EPSC,L.EPSC,ESC))
merge=merge(A,c(D.EPSC,L.EPSC,ESC))

merge=JoinLayers(merge)

obj=merge

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$stim)
table(obj$stim)

VlnPlot(merge, features = c("nFeature_RNA", "nCount_RNA"),group.by = "sample", ncol = 3)


obj=NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2,random.seed = 1, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, min.dist =1,n.neighbors = 20,seed.use =4,
               reduction = "pca", reduction.name = "umap.unintegrated",verbose = F)

CellDimPlot(obj, reduction = "umap.unintegrated", group.by = c("sample","unintegrated_clusters"),palette = "Set2")

FeaturePlot(obj, reduction = "umap.unintegrated",features=c("Pou5f1","Dppa5a","Tfcp2l1","Dppa3",
                                                            "Nr5a2","Klf2","Klf4","Esrrb","Nr0b1",
                                                            "Prdm1","Tfap2c","Tead4",
                                                            "Nanog","Tbx3","Zfp42","Fgf5"),order = T)
FeaturePlot(obj, reduction = "umap.unintegrated",features=c("Cdx2","Gata3","Krt7","Elf5"),order = T,ncol = 4)





obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE#,features = c(pre_gene2
  #,pre_gene[1:100],
  #,sub_dmm_gene[1:100],
  #,post_gene
  #)
)
obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")




DimPlot(obj, reduction = "umap.harmony", group.by = c("sample"))

FeaturePlot(obj, reduction = "umap.harmony",features=c("Nanog","Tead4","Tfap2c","Esrrb"),order = T)


sub_obj=obj[,obj$unintegrated_clusters%in%c(0,2,4,7,8,9)]
table(sub_obj$sample)
sub_obj$new_sample=sub_obj$sample
sub_obj$new_sample[grep("K5",sub_obj$new_sample)]="C4"

Stacked_VlnPlot(sub_obj, c("Pou5f1","Dppa5a","Tfcp2l1","Dppa3",
                       "Nr5a2","Klf2","Klf4","Esrrb","Nr0b1",
                       "Prdm1","Tfap2c","Tead4"),x_lab_rotate = T,#y.max = 2,
                # colors_use=c("#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA","#CCCCCC"),
                #colors_use=c("#FFD700","#DE2D26","#FF8C00","#90EE90","#87CEFA"),
                group.by = "new_sample")


sub_obj=FindVariableFeatures(sub_obj,nfeatures = 2000)

av.exp <- AverageExpression(sub_obj,slot = "data",group.by = "new_sample"
                            ,features = VariableFeatures(sub_obj)
                            )$RNA
cor.exp <- as.data.frame(cor(as.matrix(av.exp)))

bk=seq(0,1,0.01)
pheatmap::pheatmap(cor.exp,breaks = bk,display_numbers = T)



C4=obj[,obj$unintegrated_clusters%in%c(4)]
pt.mat=data.frame(sample=unique(C4$sample),cell_num=table(C4$sample)/length(C4$sample))


ggplot(pt.mat,aes(x=cell_num.Var1,y=cell_num.Freq))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  scale.func(range = c(0,6.5))+#点大小
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(0,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))


ggplot(pt.mat)+
  geom_bar(aes(x=cell_num.Var1,y=cell_num.Freq),stat = "identity")+
  coord_cartesian(ylim = c(0,0.5))+
  #scale_fill_manual(values = c(brewer.pal(n=5,"Greens")[2:5],rep("white",2)))+#ESC
  #scale_color_manual(values = c(brewer.pal(n=5,"Greens")[2:5],"#90EE90","#98FB98"))+#ESC
  
  #scale_color_manual(values = c(brewer.pal(n=5,"Oranges")[2:5],brewer.pal(n=4,"Reds")[2:4]))+#TSC
  #scale_fill_manual(values = c(brewer.pal(n=5,"Oranges")[2:5],rep("white",3)))+#TSC
  
  #scale_color_manual(values = c(color,"#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA"))+#PGCLC
  #scale_fill_manual(values = c(color,rep("white",5)))+#PGCLC
  
  labs(x=NULL,y="Percentage of time point in C4")+
  guides(fill = FALSE,color=F)+
  theme(panel.grid = element_blank())+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,hjust=1))



###After integrate in vivo data####

etp_in_vivo=readRDS("K5_ETP_new_integrate_invivo_add_EpiSC.rds")

etp_in_vivo

meta=etp_in_vivo@meta.data
meta[grepl("PGCLC",meta$Celltype),"Celltype"]=meta[grepl("PGCLC",meta$Celltype),"etp_cl"]
meta[grepl("TSC",meta$Celltype),"Celltype"]=meta[grepl("TSC",meta$Celltype),"etp_cl"]
meta[grepl("ESC",meta$Celltype),"Celltype"]=meta[grepl("ESC",meta$Celltype),"etp_cl"]
meta[grepl("ESC",meta$Celltype),"Celltype"]=meta[grepl("ESC",meta$Celltype),"etp_cl"]
meta$Celltype=gsub("18-ESC","iPSC(18)",meta$Celltype)
meta$Celltype=gsub("20-PGCLC","PGCLC(20)",meta$Celltype)
meta$Celltype=gsub("26-TSC","TSCLC(26)",meta$Celltype)
meta$Celltype=gsub("27-TSC","TSCLC(27)",meta$Celltype)
meta$Celltype=gsub("28-TSC","TSCLC(28)",meta$Celltype)


etp_in_vivo@meta.data$Celltype2=meta$Celltype
etp_in_vivo@meta.data$Celltype2=factor(etp_in_vivo@meta.data$Celltype2,
                                       levels = c("E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI",
                                                  "E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE",
                                                  "E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC",
                                                  "iPSC(18)","PGCLC(20)","TSCLC(26)","TSCLC(28)",
                                                  "TSCLC(27)","EpiSC"))
table(etp_in_vivo$Celltype)
table(etp_in_vivo$Celltype2)

#####UMAP####
etp_in_vivo=RunUMAP(etp_in_vivo,reduction="harmony",seed.use =6,
                   n.neighbors =50,
                   n.epochs =500, min.dist=0.96,
                   dims=1:30,verbose=F)


DimPlot(etp_in_vivo,reduction = "umap",group.by = c("Celltype2"),pt.size = 1.5,
        shape.by = "stim",label = T,label.size = 5,
        cols=c(brewer.pal(n=7,"Greens")[2:5],brewer.pal(n=6,"Oranges")[2:5],brewer.pal(n=7,"Blues")[2:5],
               #"#90EE90","#87CEFA",brewer.pal(n=4,"Reds")[2:4]
               colorRampPalette(brewer.pal(n=9,"Greys")[c(2:3,5,4)])(6)
               ))+
scale_shape_manual(values = c(17,1,18))+
  scale_x_continuous(limits = c(-20,20))
#7x12

seurat_clusters=readRDS("seurat_clusters_1.7_43.rds")
etp_in_vivo$seurat_clusters=seurat_clusters

DimPlot(etp_in_vivo,reduction = "umap",group.by = c("seurat_clusters"),pt.size = 1.5,label = T,label.size = 5,
        cols=colorRampPalette(colors=palette_list$Set2)(29))+
  scale_shape_manual(values = c(17,1,18))+
  scale_x_continuous(limits = c(-20,20))
#7x12

etp_in_vivo$sample=as.character(etp_in_vivo$sample)
etp_in_vivo$sample[is.na(etp_in_vivo$sample)]="in_vivo_sample"
etp_in_vivo$sample=factor(etp_in_vivo$sample,levels = unique(etp_in_vivo$sample)[c(7,5,1:4,6)])

DimPlot_scCustom(etp_in_vivo,reduction = "umap",group.by = c("sample"),pt.size = 1,
                 label = F,label.size = 5,shuffle = T,
        colors_use =c(sample_palette_col,"#DDDDDD"))+
  scale_shape_manual(values = c(17,1,18))+
  scale_x_continuous(limits = c(-20,20))
#7x12


immune.combined=readRDS("immune.combined_add_EpiSC.rds")
DefaultAssay(immune.combined)="RNA"


#####Identity Score####
p=DotPlot(cols = "RdYlBu",immune.combined,features=c("ESC_AUC","TSC_AUC","PGCLC_AUC","Primed_AUC","Klf5"),
          group.by = "Celltype",#scale =F,#
          assay = "RNA",scale=T)
p_data=p$data

our=c("iPSC-sub","TSCLC-subTE3.5", "TSCLC-subTE4.5","TSCLC-subExE","PGCLC-sub")

sub_p=p_data[p_data$features.plot=="PGCLC_AUC",]
sub_p=sub_p[sub_p$id%in%c("E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI",our),]
sub_p=sub_p[sub_p$id%in%c("E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE",our),]
sub_p=sub_p[sub_p$id%in%c("E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC",our),]

sub_p$id=factor(sub_p$id,levels =c("E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC",our))

color= brewer.pal(n=5,"Greens")[2:5]#ESC
color= brewer.pal(n=5,"Oranges")[2:5]#TSC
color= brewer.pal(n=5,"Blues")[2:5]#PGCLC
           
ggplot(sub_p)+
  geom_bar(aes(x=id,y=avg.exp,fill=id,color=id),stat = "identity")+
  coord_cartesian(ylim = c(0.25,0.35))+
  #scale_fill_manual(values = c(brewer.pal(n=5,"Greens")[2:5],rep("white",2)))+#ESC
  #scale_color_manual(values = c(brewer.pal(n=5,"Greens")[2:5],"#90EE90","#98FB98"))+#ESC
  
  #scale_color_manual(values = c(brewer.pal(n=5,"Oranges")[2:5],brewer.pal(n=4,"Reds")[2:4]))+#TSC
  #scale_fill_manual(values = c(brewer.pal(n=5,"Oranges")[2:5],rep("white",3)))+#TSC

  scale_color_manual(values = c(color,"#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA"))+#PGCLC
  scale_fill_manual(values = c(color,rep("white",5)))+#PGCLC
  
  labs(x=NULL,y="PGCLC Identity Score")+
  guides(fill = FALSE,color=F)+
  theme(panel.grid = element_blank())+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,hjust=1))

#####Marker dotplot####

invivo=c("E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI")
invivo=c("E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE")
invivo=c("E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC")
#our=c("iPSC-sub","iPSC","TSCLC-subTE3.5", "TSCLC-subTE4.5","TSCLC-subExE","TSCLC","PGCLC-sub","PGCLC")
our=c("iPSC-sub","TSCLC-subTE3.5", "TSCLC-subTE4.5","TSCLC-subExE","PGCLC-sub")

use_ftr=c("Klf4","Tbx3","Nanog","Zfp42")
use_ftr=c("Cdx2","Gata3","Tead4","Ly6a","Elf5","Hand1","Ets2","Ascl2")
use_ftr=c("Prdm1","Prdm14","Dppa3","Nanos3")

in_vivo_part=immune.combined[,immune.combined$Celltype%in%invivo]
our_part=immune.combined[,immune.combined$Celltype%in%our]

p=DotPlot(cols = "RdYlBu",in_vivo_part,features=use_ftr,group.by = "Celltype",#scale =F,
          assay = "RNA",scale=T)
p=DotPlot(cols = "RdYlBu",our_part,features=use_ftr,group.by = "Celltype",#scale =F,
          assay = "RNA",scale=T)
p_data=p$data

p_data$features.plot=factor(p_data$features.plot,levels =rev(use_ftr))
p_data$id=factor(p_data$id,levels =our)
p_data$id=factor(p_data$id,levels =invivo)

scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size, #该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

#name="Figure.2_Stringency_Single_cell_data_new_ver1"
#pdf(paste(name,".pdf",sep=""),width=5.2,height =5)
color = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"   ))((50)),
          colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026" ))((50)))

color = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((50)),
          colorRampPalette(colors = c("#FFFFFF","#D7B5D8","#DF65B0","#DD1C77","#980043" ))((50)))

#p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
#该函数可以完全替代scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size=p_data$pct.exp),color="black",shape=21)+
  scale.func(range = c(0,10),limits = c(0,100))+#点大小
  guides(size = guide_legend(title = "Percent Expressed"))+
  #guides(size = F)+
  scale_fill_gradientn(colors = color,limits=c(-2,2),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  #facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  #facet_grid(vars(label),scales = "free",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))


#####Marker Vlnplot####

immune.combined=readRDS("K5_0-5D_select_etp_cl.rds")

use_ftr=c("Klf5","Eomes","Tead4","Gata3","Id2",
          "Krt7","Plet1","Tfap2c","Ets2",
          "Krt18","Ascl2","Krt19",
          "Klf4","Tbx3","Nanog",
          "Dppa3","Prdm1","Nanos3")

our=c("28-TSC","26-TSC","27-TSC","18-ESC","20-PGCLC")

immune.combined$etp_cl=factor(our_part$etp_cl,levels = our)

Stacked_VlnPlot(immune.combined, use_ftr,x_lab_rotate = T,#y.max = 2,
                # colors_use=c("#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA","#CCCCCC"),
                colors_use=c("#FFD700","#DE2D26","#FF8C00","#90EE90","#87CEFA"),
                group.by = "etp_cl")
#10x3
immune.combined=readRDS("immune.combined_add_EpiSC.rds")
DefaultAssay(immune.combined)="RNA"

invivo=c("E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE",
         "E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI",
         "E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC","EpiSC")

#our=c("iPSC-sub","TSCLC-subTE3.5", "TSCLC-subTE4.5","TSCLC-subExE","PGCLC-sub","EpiSC")

in_vivo_part=immune.combined[,immune.combined$Celltype%in%invivo]
#our_part=immune.combined[,immune.combined$Celltype%in%our]
in_vivo_part$Celltype=factor(in_vivo_part$Celltype,levels = invivo)
#our_part$Celltype=factor(our_part$Celltype,levels = our)

color= c(brewer.pal(n=5,"Oranges")[2:5],#TSC
         brewer.pal(n=5,"Greens")[2:5],#ESC
brewer.pal(n=5,"Blues")[2:5],"grey")#PGCLC


# feature="Klf4"
# VlnPlot(in_vivo_part,feature,
#         cols=c(color),
#         group.by = "Celltype") + 
#   xlab("") + ylab(feature) + ggtitle("") + 
#   theme(legend.position = "none", 
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank(), 
#         axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5), 
#         axis.text.y = element_text(size = rel(1)) ) 

Stacked_VlnPlot(in_vivo_part, use_ftr,#y.max = 2,
               colors_use=c(color),
               x_lab_rotate = T,
               group.by = "Celltype")
#10x5



#####Sub Celltype in Lineage umap####
monocle$invivo_celltype="Others"

meta=monocle@meta.data
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(25)]),"invivo_celltype"]="TSCLC-subTE3.5"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(19)]),"invivo_celltype"]="TSCLC-subTE4.5"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(16)]),"invivo_celltype"]="TSCLC-subExE"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(3,2)]),"invivo_celltype"]="PGCLC-sub"
meta[rownames(meta)%in%colnames(etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(0,1)]),"invivo_celltype"]="iPSC-sub"
table(meta$invivo_celltype)
monocle@meta.data=meta


CellDimPlot(monocle, group.by = "invivo_celltype", reduction = "umap",pt.size = .3)



sub=monocle[,monocle$invivo_celltype%in%c("iPSC-sub1")]
#sub=testseurat[,testseurat$Celltype%in%"PGCLC_3"]
table(sub$sample)
table(sub$stim)



#####sub cell type percentage in clusters####

meta=read.csv("PAGA_annodata_meta.csv")



table(meta[meta$invivo_celltype=="iPSC-sub","leiden"])/table(meta$invivo_celltype)[1]
table(meta[meta$invivo_celltype=="PGCLC-sub","leiden"])/table(meta$invivo_celltype)[3]
table(meta[meta$invivo_celltype=="TSCLC-subExE","leiden"])/table(meta$invivo_celltype)[4]
table(meta[meta$invivo_celltype=="TSCLC-subTE3.5","leiden"])/table(meta$invivo_celltype)[5]
table(meta[meta$invivo_celltype=="TSCLC-subTE4.5","leiden"])/table(meta$invivo_celltype)[6]



#K2,4,5 D5 and in vivo data####
seurat=readRDS("K245_D5_ver2.rds")


table(seurat$sample)

seurat <- NormalizeData(seurat)
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

###Findclusters and Run tSNE####

# Determine percent of variation associated with each PC
pct <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

seurat <- RunTSNE(seurat, dims = 1:pcs,seed.use=111)
seurat <- FindNeighbors(seurat, dims = 1:pcs)
seurat <- FindClusters(seurat, resolution =1, random.seed =90)


for(i in 100:200){
  seurat <- RunTSNE(seurat, dims = 1:pcs,seed.use=i)
  file_name=paste("K2,4,5_D5_cl_tsneseed",i,".png",sep="_")
  p1=DimPlot(seurat,reduction = "tsne", group.by =c("seurat_clusters"),raster = F,label.size = 6,
             label =T,cols=colorRampPalette(colors=palette_list$Set2)(25))
  p2=FeatureDimPlot(seurat,raster =F, pt.size = .1,
                 features = c("ESC_AUC"),lower_cutoff = 0.1,upper_cutoff = 0.2,
                 reduction = "tsne", theme_use = "theme_blank")
  p3=FeatureDimPlot(seurat,raster = F, pt.size = .1,
                 features = c("TSC_AUC"),lower_cutoff = 0.1,upper_cutoff = 0.21,
                 reduction = "tsne", theme_use = "theme_blank")
  p4=FeatureDimPlot(seurat,raster = F, pt.size = .1,
                 features = c("PGCLC_AUC"),lower_cutoff =0.17,upper_cutoff = 0.27, 
                 reduction = "tsne", theme_use = "theme_blank")
  ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow =2)
  ggsave(filename = file_name,width = 10,height = 10,dpi = 200)
}


seurat <- RunTSNE(seurat, dims = 1:16,seed.use=7)
seurat <- FindNeighbors(seurat, dims = 1:16)

  for(n in 1:200){
    for(r in seq(0.6,1.2,0.1)){
      file_name=paste("K2,4,5_D5_cl_tsneseed",7,"clseed",n,"resolution",r,".png",sep="_")
      print(c(file_name))
      seurat <- FindClusters(seurat, resolution =r, random.seed =n)
      
      p1=DimPlot(seurat,reduction = "tsne", group.by =c("seurat_clusters"),raster = F,label.size = 6,
                 label =T,cols=colorRampPalette(colors=palette_list$Set2)(33))
      p2=FeatureDimPlot(seurat,raster =F, pt.size = .1,
                        features = c("ESC_AUC"),lower_cutoff = 0.11,upper_cutoff = 0.205,
                        reduction = "tsne", theme_use = "theme_blank")
      p3=FeatureDimPlot(seurat,raster = F, pt.size = .1,
                        features = c("TSC_AUC"),lower_cutoff = 0.155,upper_cutoff = 0.20,
                        reduction = "tsne", theme_use = "theme_blank")
      p4=FeatureDimPlot(seurat,raster = F, pt.size = .1,
                        features = c("PGCLC_AUC"),lower_cutoff =0.19,upper_cutoff = 0.265, 
                        reduction = "tsne", theme_use = "theme_blank")
      ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow =2)
      ggsave(filename = file_name,width = 15,height = 10,dpi = 200)
    }}


###tSNE plot####
seurat <- RunTSNE(seurat, dims = 1:pcs,seed.use=111)
seurat <- FindNeighbors(seurat, dims = 1:pcs)
seurat <- FindClusters(seurat, resolution =1, random.seed =90)


sample_palette_col=c("#DC143C","#9400D3","#FFB6C1","grey")



CellDimPlot(seurat,reduction = "tsne", group.by =c("sample"),raster = F,
            label = F,label.fg = "white",label_insitu = T,palcolor = sample_palette_col )

DimPlot(seurat,reduction = "tsne", group.by =c("seurat_clusters"),raster = F,label.size = 6,
        label =T,cols=colorRampPalette(colors=palette_list$Set2)(33))
#10x10


FeatureDimPlot(seurat,raster =F, pt.size = .1,
               features = c("ESC_AUC"),lower_cutoff = 0.12,upper_cutoff = 0.20,
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat,raster = F, pt.size = .1,
               features = c("TSC_AUC"),lower_cutoff = 0.15,upper_cutoff = 0.20,
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat,raster = F, pt.size = .1,
               features = c("PGCLC_AUC"),lower_cutoff =0.25,upper_cutoff = 0.31, 
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat,raster = F, pt.size = .1,
               features = c("Primed_AUC"),lower_cutoff = 0.15,upper_cutoff = 0.2,  
               reduction = "tsne", theme_use = "theme_blank")

FeatureDimPlot(seurat,raster = F, pt.size = .1,
               features = c("Bmp_sig_AUC"),#lower_cutoff = 0.15,upper_cutoff = 0.2,  
               reduction = "tsne", theme_use = "theme_blank")
#10x10


subseurat=seurat[,seurat$seurat_clusters%in%c(21,11,9,14)]
meta=subseurat@meta.data
meta[meta$seurat_clusters%in%c(11,9),"Celltype"]="iPSC"
meta[meta$seurat_clusters%in%c(14),"Celltype"]="PGCLC"
meta[meta$seurat_clusters%in%c(21),"Celltype"]="TSCLC"
meta$Celltype=paste0(meta$Celltype,"(",meta$seurat_clusters,")")
subseurat@meta.data=meta
table(subseurat$Celltype)
saveRDS(subseurat,"K245_sub_etp.rds")


p=DotPlot(seurat,group.by = "seurat_clusters",features=c("TSC_AUC","ESC_AUC","Zfp42","Tbx3","Dppa5a",
                                                         "Cdx2","Gata3","Tead4","Plet1",
                                                         "Krt7","Elf5","Ascl2","Ets2"))
p+ggplot2::coord_flip()


subseurat14=seurat[,seurat$seurat_clusters%in%c(9,11)]
meta=subseurat14@meta.data
meta[meta$seurat_clusters%in%c(9,11),"Celltype"]="iPSC"
subseurat14@meta.data=meta
table(subseurat14$Celltype)
saveRDS(subseurat14,"K245_sub_iPSC,9,11.rds")



###After integrate in vivo data####

pgc=readRDS("K245_PGC_in_vivo.rds")
tsc=readRDS("K245_TSCLC_in_vivo.rds")
ipsc=readRDS("K245_iPSC_in_vivo.rds")


#####UMAP PGC####
etp_in_vivo=pgc

meta=etp_in_vivo@meta.data
meta[grepl("D5_K2",meta$sample),"Celltype"]="iKlf2-PGCLC"
meta[grepl("D5_K4",meta$sample),"Celltype"]="iKlf4-PGCLC"
meta[grepl("D5_K5",meta$sample),"Celltype"]="iKlf5-PGCLC"
meta$Celltype=factor(meta$Celltype,levels=c("E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI",
                                            "E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE",
                                            "E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC",
                                            "iKlf5-PGCLC","iKlf2-PGCLC","iKlf4-PGCLC","EpiSC"))
etp_in_vivo@meta.data=meta

pct <- etp_in_vivo[["harmony"]]@stdev / sum(etp_in_vivo[["harmony"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs

# etp_in_vivo <- FindNeighbors(etp_in_vivo, dims = 1:pcs)
# etp_in_vivo <- FindClusters(etp_in_vivo, resolution = 2,random.seed = 43)


etp_in_vivo=RunUMAP(etp_in_vivo,reduction="harmony",seed.use =3,
                    n.neighbors =50,
                    n.epochs =500, min.dist=0.8,
                    dims=1:pcs,verbose=F)

#PGC
p1=DimPlot(etp_in_vivo,reduction = "umap",group.by = c("Celltype"),pt.size = 1.5,
        shape.by = "stim",label = T,label.size = 5,
        cols=c(#colorRampPalette(brewer.pal(n=9,"Greys")[3:7])(8),brewer.pal(n=7,"Blues")[2:5],
          brewer.pal(n=7,"Greens")[2:5],brewer.pal(n=6,"Oranges")[2:5],brewer.pal(n=7,"Blues")[2:5],
               #"#90EE90","#87CEFA",brewer.pal(n=4,"Reds")[2:4]
               c("#DC143C","#D2B48C","#008B8B","grey")
        ))+
  scale_shape_manual(values = c(17,20,18))
  #+scale_x_continuous(limits = c(-20,20))
#8x10

p2=DimPlot(etp_in_vivo,reduction = "umap",group.by = c("seurat_clusters"),pt.size = 1.5,label = T,label.size = 5,
        cols=colorRampPalette(colors=palette_list$Set2)(23))+
  scale_shape_manual(values = c(17,20,18))
#8x10
p1+p2


FeaturePlot(etp_in_vivo,reduction = "umap",features = c("Prdm1","Dppa3","Nanos3"),min.cutoff = "q1",
            pt.size = 1.5,label = F,order = T,ncol = 3)
#23x8


sub=etp_in_vivo[,etp_in_vivo$seurat_clusters%in%c(0,8,19)]
sub=sub[,sub$sample%in%c("D5_K2","D5_K4","D5_K5","EpiSCs")]
table(sub$sample)
p=DotPlot(sub,group.by = "sample",features=c("PGC_AUC","Prdm1","Dppa3","Nanos3","Tfap2c"))
p+ggplot2::coord_flip()




#####UMAP TSC####
etp_in_vivo=tsc

meta=etp_in_vivo@meta.data
meta[grepl("D5_K2",meta$sample),"Celltype"]="iKlf2-TSCLC"
meta[grepl("D5_K4",meta$sample),"Celltype"]="iKlf4-TSCLC"
meta[grepl("D5_K5",meta$sample),"Celltype"]="iKlf5-TSCLC"
meta$Celltype=factor(meta$Celltype,levels=c("E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI",
                                            "E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE",
                                            "E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC",
                                            "iKlf5-TSCLC","iKlf2-TSCLC","iKlf4-TSCLC","EpiSC"))
etp_in_vivo@meta.data=meta

pct <- etp_in_vivo[["harmony"]]@stdev / sum(etp_in_vivo[["harmony"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs

# etp_in_vivo <- FindNeighbors(etp_in_vivo, dims = 1:pcs)
# etp_in_vivo <- FindClusters(etp_in_vivo, resolution = 1.5,random.seed = 1)


etp_in_vivo=RunUMAP(etp_in_vivo,reduction="harmony",seed.use =6,
                    n.neighbors =20,
                    n.epochs =500, min.dist=0.8,
                    dims=1:pcs,verbose=F)

p1=DimPlot(etp_in_vivo,reduction = "umap",group.by = c("Celltype"),pt.size = 1.5,
           shape.by = "stim",label = T,label.size = 5,
           cols=c(#colorRampPalette(brewer.pal(n=9,"Greys")[3:7])(8),brewer.pal(n=6,"Oranges")[2:5],
             brewer.pal(n=7,"Greens")[2:5],brewer.pal(n=6,"Oranges")[2:5],brewer.pal(n=7,"Blues")[2:5],
             #"#90EE90","#87CEFA",brewer.pal(n=4,"Reds")[2:4]
                  c("#DC143C","#D2B48C","#008B8B","grey")
           ))+
  scale_shape_manual(values = c(17,20,18))
#+scale_x_continuous(limits = c(-20,20))
#8x10

p2=DimPlot(etp_in_vivo,reduction = "umap",group.by = c("seurat_clusters"),pt.size = 1.5,label = T,label.size = 5,
           cols=colorRampPalette(colors=palette_list$Set2)(23))+
  scale_shape_manual(values = c(17,20,18))
#8x10
p1+p2
#16x8


FeaturePlot(etp_in_vivo,reduction = "umap",features = c("Gata3","Elf5","Ascl2"),
            pt.size = 1.5,label = F,order = T,ncol = 3)
#23x8


#####UMAP iPSC####
etp_in_vivo=ipsc

meta=etp_in_vivo@meta.data
meta[grepl("D5_K2",meta$sample),"Celltype"]="iKlf2-iPSC"
meta[grepl("D5_K4",meta$sample),"Celltype"]="iKlf4-iPSC"
meta[grepl("D5_K5",meta$sample),"Celltype"]="iKlf5-iPSC"
meta$Celltype=factor(meta$Celltype,levels=c("E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI",
                                            "E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE",
                                            "E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC",
                                            "iKlf5-iPSC","iKlf2-iPSC","iKlf4-iPSC","EpiSC"))
 etp_in_vivo@meta.data=meta
 
pct <- etp_in_vivo[["harmony"]]@stdev / sum(etp_in_vivo[["harmony"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs


etp_in_vivo=RunUMAP(etp_in_vivo,reduction="harmony",seed.use =6,
                    n.neighbors =50,
                    n.epochs =500, min.dist=0.8,
                    dims=1:pcs,verbose=F)
 
p1=DimPlot(etp_in_vivo,reduction = "umap",group.by = c("Celltype"),pt.size = 1.5,
           shape.by = "stim",label = T,label.size = 5,
           cols=c(#colorRampPalette(brewer.pal(n=9,"Greys")[3:7])(8),brewer.pal(n=7,"Greens")[2:5],
                  brewer.pal(n=7,"Greens")[2:5],brewer.pal(n=6,"Oranges")[2:5],brewer.pal(n=7,"Blues")[2:5],
                  #"#90EE90","#87CEFA",brewer.pal(n=4,"Reds")[2:4]
                  c("#DC143C","#D2B48C","#008B8B","grey")
           ))+
  scale_shape_manual(values = c(17,20,18))
#+scale_x_continuous(limits = c(-20,20))
#8x10

p2=DimPlot(etp_in_vivo,reduction = "umap",group.by = c("seurat_clusters"),pt.size = 1.5,label = T,label.size = 5,
           cols=colorRampPalette(colors=palette_list$Set2)(23))+
  scale_shape_manual(values = c(17,20,18))
#8x10
p1+p2
#16x8


FeaturePlot(etp_in_vivo,reduction = "umap",features = c("Klf4","Tbx3","Zfp42",
                                                        "Prdm1","Dppa3","Nanos3"),min.cutoff = "q50",
            pt.size = 1.5,label = F,order = T,ncol = 3)
#23x8


#####K2,4,5 ETP percentage####


sub_pgc=pgc[,pgc$seurat_clusters%in%c(8,19,0)]

table(sub_pgc$sample)[3]/table(seurat$sample)[1]
table(sub_pgc$sample)[1]/table(seurat$sample)[2]
table(sub_pgc$sample)[2]/table(seurat$sample)[3]
table(sub_pgc$sample)[4]/table(seurat$sample)[4]


finmat=data.frame(Var1=factor(c("iKlf5","iKlf2","iKlf4","EpiSCs"),
                              levels = c("iKlf5","iKlf2","iKlf4","EpiSCs")),
                  Percentage=c(0.0235711,0.004659643,0.001390369,0.0002111932))
    
 ggplot(finmat)+
      geom_bar(aes(x=Var1,y=Percentage,fill=Var1),stat = "identity")+
      scale_fill_manual(values = c("#DC143C","#D2B48C","#008B8B","#FFB6C1"))+
      #geom_line(aes(x=Var1,y=Percentage,group = 1))+
      #geom_point(aes(x=Var1,y=Percentage,group = 1))+
      scale_y_continuous(limits = c(0,0.025),labels=scales::percent)+
      labs(x=NULL,y="Percentage %",title="PGCLC")+
      #geom_text(aes(x=Var1,y=Percentage,label = paste0("n=",Freq)),size=2, position = position_dodge(0.9), vjust = -0.8)+
      #guides(fill = FALSE)+
      theme(panel.grid = element_blank())+
      theme_classic()+
      theme(axis.text.x = element_text(angle=90,hjust=1)
            #axis.text.y = element_text(size = 5)
      )

#3x3

#Compare to BMP4i TSC-like####
 
###Select TSC group####
 setwd("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/BMP4/")
 list=list.files()

 
 for( f in list){
   
   name=strsplit2(f,"_")[,2]
 bmp=Read10X(paste0(f,"/",name,"Solo.out/"))
 bmp <- bmp[!grepl("_", rownames(bmp)),]
 bmp=CreateSeuratObject(counts = bmp)
 bmp[["percent.mt"]] <- PercentageFeatureSet(bmp, pattern = "^mt-")
 bmp=subset(bmp, subset = nFeature_RNA > 500 & percent.mt < 20 & nCount_RNA > 1000)
 bmp$sample=name
saveRDS(bmp,paste0(name,"_BMP4.rds"))
 }
D0=readRDS("PNTD0_BMP4.rds")
D1=readRDS("PNTD1_BMP4.rds")
D2=readRDS("PNTD2_BMP4.rds")
D3=readRDS("PNTD3_BMP4.rds")
D5=readRDS("PNTD5_BMP4.rds")
D8=readRDS("PNTD8_BMP4.rds")
bmp=merge(D0,c(D1,D2,D3,D5,D8))
saveRDS(bmp,"BMP4_D0-8.rds")
bmp=JoinLayers(bmp)


obj=readRDS("BMP4_D0-8.rds")
obj=obj[!rownames(obj)%in%grep(pattern = "\\.",rownames(obj),value=T),]


obj=merge(D5,D8)

obj=JoinLayers(obj)
obj=obj[rownames(obj)%in%c(rownames(seurat),rownames(episc)),]

VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)

  
  obj <- NormalizeData(obj)
  obj <- ScaleData(obj)
  
  obj=FindVariableFeatures(obj,nfeatures = 4000)
  
  
  #obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  
  
  
  obj <- RunPCA(obj)
  
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:50)
  obj <- FindClusters(obj, resolution = 1,random.seed =1, cluster.name = c("seurat_clusters"))
  #obj <- FindClusters(obj, resolution = 0.65,random.seed =1, cluster.name = c("seurat_clusters"))
  
  obj = RunTSNE(obj)
  obj <- RunUMAP(obj, reduction = "pca",
                 dims = 1:50, seed.use =1,
                 n.neighbors = 20, #20更多粘连和小块,50偏整体http://biotrainee.cn:9013/graphics/plot_zoom_png?width=912&height=990
                 n.epochs = 200,#200大数据，500小数据
                 min.dist = 1,#距离
                 #spread=1,#揉
                 reduction.name = "umap",verbose = FALSE)
  use_ftr=c("Klf4","Tbx3","Nanog","Zfp42")
  
  use_ftr=c("Cdx2","Gata3","Tead4","Ly6a",
            "Elf5","Ddah1","Tfap2c","Plac1",
            "Hand1","Ets2","Peg10","Phlda2",
            "Fos","Klf6","Ctsb","Ascl2")
  use_ftr=c("Cdx2", "Elf5", "Gata3", "Tead4", "Tfap2c")
  
  degres_vivo_filter[degres_vivo_filter$cluster=="E4.5_TE",]
  use_ftr=c("Prdm1","Prdm14","Dppa3","Nanos3","Kit","Klf5")
  
  DimPlot(obj,reduction = "umap",group.by = c("sample"),label = T,label.size = 8)
  
  FeaturePlot(obj,features = use_ftr,order = T,ncol = 3,min.cutoff = "q1",max.cutoff = "q99")
  FeaturePlot(obj,features = "Prdm1",order = T,ncol = 1,min.cutoff = "q1",max.cutoff = "q99")
  
  
  obj$seurat_clusters=paste(obj$seurat_clusters,obj$sample)
  

  VlnPlot(obj, "TSC_AUC",group.by = "seurat_clusters",pt.size = 0)
  Stacked_VlnPlot(obj, "TSC_AUC",x_lab_rotate = T,#y.max = 2,
                     # colors_use=c("#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA","#CCCCCC"),
                     #colors_use=c("#FFD700","#DE2D26","#FF8C00","#90EE90","#87CEFA"),
                     group.by = "seurat_clusters")
  sub=obj[,obj$seurat_clusters%in%c(9,14)]
  table(sub$sample)
  
  
  
  
  
  
  
  
  p1+p2
  
  DoHeatmap(obj,features = zto4_gene,group.by = "seurat_clusters",disp.min = -1,disp.max = 1)
  table(obj$seurat_clusters)
}

gene_set=list(TSC=tsc_relative)
              
dg <- as(as.matrix(obj@assays$RNA$data),"sparseMatrix")
cells_rankings <- AUCell_buildRankings(dg,plotStats = F,verbose = F)
class(dg)
dg=NULL
gc()
#saveRDS(cells_rankings,"pancreas_sub_cells_rankings.rds")

#cells_rankings=readRDS("pancreas_sub_cells_rankings.rds")
cells_AUC <- AUCell_calcAUC(gene_set, cells_rankings,aucMaxRank = 3000, nCores=10)
obj$TSC_AUC  <- as.numeric(getAUC(cells_AUC)["TSC", ])
sub=obj[,obj$seurat_clusters%in%c(9,14)]
sub$etp_cl="BMPi-TSC"
saveRDS(sub,"bmp_tsc_D5,8.rds")


###Comparing Marker gene####
select_etp=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/K5_0-5D_select_etp_cl.rds")
select_etp=select_etp[,!select_etp$etp_cl%in%c("20-PGCLC")]
bmp_tsc=readRDS("bmp_tsc_D5,8.rds")


obj=merge(select_etp,bmp_tsc)
obj=JoinLayers(obj)


VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)


obj <- NormalizeData(obj)
obj <- ScaleData(obj)

obj=FindVariableFeatures(obj,nfeatures = 4000)


color = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"   ))((50)),
          colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026" ))((50)))


use_ftr=rev(use_ftr)
p <- DotPlot(invivo_TE, features=use_ftr,
             scale.min = 0,
             scale.max = 100,
             scale = T,
             group.by='Celltype',
             cols =  c("RdYlBu"))
p <- DotPlot(obj, features=use_ftr,
             scale.min = 0,
             scale.max = 100,
             scale = T,
             group.by='etp_cl',
             cols =  c("RdYlBu"))
p+ggplot2::coord_flip()
p_data=p$data
p_anno=data.frame(features.plot=levels(p_data$features.plot),
                  label=factor(c(rep("E3.5 TE Markers",4),
                                 rep("E4.5 TE Markers",4),
                                 rep("E5.5 ExE Markers",4),
                                 rep("E6.5 ExE Markers",4))))

p_data=plyr::join(p_data,p_anno)

p_data$features.plot=factor(p_data$features.plot,levels = rev(use_ftr))
p_data$id=factor(p_data$id,levels = c("18-ESC","28-TSC","26-TSC","27-TSC","BMPi-TSC"))
p_data
scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size, #该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
limit_max=2
limit_min=-2

p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
p_data[p_data$avg.exp.scaled<=limit_min,"avg.exp.scaled"]=as.numeric(limit_min)
#该函数可以完全替代scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  #scale.func(range = c(0,6.5))+#点大小
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(limit_min,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#添加边框
  facet_grid(rows = vars(label), scales = "free_y")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))

Stacked_VlnPlot(obj, use_ftr,x_lab_rotate = T,#y.max = 2,
                # colors_use=c("#90EE90",brewer.pal(n=4,"Reds")[2:4],"#87CEFA","#CCCCCC"),
                #colors_use=c("#FFD700","#DE2D26","#FF8C00","#90EE90","#87CEFA"),
                group.by = "etp_cl")

###Comparing GSVA####
library(GSVA)

invivo=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/in_vivo_data/E3.5_6.5_ICM_TE.Rds")
invivo=CreateSeuratObject(invivo@assays$RNA$counts,meta.data = invivo@meta.data)
invivo=NormalizeData(invivo)

invivo=NormalizeData(invivo)
Idents(invivo)=invivo$etp_cl
degres_vivo=FindAllMarkers(invivo,min.pct = 0.5,
                           logfc.threshold = 0.25,only.pos = T)

filter_gene=
  function(deg_res){
    deg_res_filter=deg_res[deg_res$p_val_adj<0.01,]
    deg_res_filter=deg_res_filter[deg_res_filter$avg_log2FC>1,]
    deg_res_filter=deg_res_filter[order(deg_res_filter$avg_log2FC,decreasing = T),]
    return(deg_res_filter)
  }

degres_vivo_filter=filter_gene(degres_vivo)
table(degres_vivo_filter$cluster)

genelist = split(degres_vivo_filter$gene, degres_vivo_filter$cluster)
genelist_1000 <- lapply(genelist, function(x) head(x, 500))

av.exp <- AverageExpression(obj,group.by = "etp_cl")$RNA
gsvapar =gsvaParam(av.exp,genelist_1000)
gsva_es <- as.data.frame(gsva(gsvapar))

bk=seq(-0.4,0.5,0.01)
gsva_es$celltype=rownames(gsva_es)
gsva_es$celltype=factor(gsva_es$celltype, levels = c("E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE",
                                                     "E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI"))
gsva_es=gsva_es[order(gsva_es$celltype),]
gsva_es$celltype=NULL
pheatmap::pheatmap(gsva_es[1:4,c(1,4,2,3,5)],show_rownames = T,
                   cluster_cols = F,cluster_rows =F,cellwidth = 20,cellheight =20,
                   breaks = bk,main = "GSVA",border_color = "black",
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5))
)

STUD



###Comparing in vivodata####

invivo=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/in_vivo_data/E3.5_6.5_ICM_TE.Rds")
invivo=CreateSeuratObject(invivo@assays$RNA$counts,meta.data = invivo@meta.data)
invivo=NormalizeData(invivo)

av.exp <- AverageExpression(invivo,group.by = "Celltype")$RNA

min.cells=dim(obj@assays$RNA$counts)[2]*0.01

obj <- RunKNNPredict(srt_query = obj, bulk_ref = av.exp,filter_lowfreq = min.cells)

sub1=obj[,obj$etp_cl=="28-TSC"]
table(sub1$KNNPredict_classification)


sub_sddip2=subset(sddip2,downsample=2000)
sub_sddip2$sample="sddip2"


bmp_tsc=readRDS("bmp_tsc_D5,8.rds")
bmp_tsc$stim="iBMP"
select_etp$stim="iKlf5"
invivo$stim="invivo"
invivo$etp_cl=invivo$Celltype

invivo_TE=invivo[,invivo$Celltype%in%c("E3.5_TE","E4.5_TE" ,"E5.5_ExE","E6.5_ExE" )]

obj=merge(select_etp,c(bmp_tsc,invivo))
obj=JoinLayers(obj)



table(obj$stim)
table(obj$etp_cl)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$stim)


obj <- NormalizeData(obj)
obj <- ScaleData(obj)


VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "sample", ncol = 3)



obj=FindVariableFeatures(obj,nfeatures =c(20*100))

use_feature=VariableFeatures(obj)


obj <- RunPCA(obj
              #,features = in_vivo_gene[1:3000]
)



obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,#normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "integrated",
  #features =use_feature[1:2000],
  verbose = FALSE#,k.weight=48
)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


obj <- FindNeighbors(obj, reduction = "integrated", dims = 1:50)
obj <- FindClusters(obj, resolution = 1,random.seed = 1, cluster.name = "integrated_clusters")

obj <- RunUMAP(obj, reduction = "integrated",
               dims = 1:50, seed.use =2,
               n.neighbors = 20, #20更多粘连和小块,50偏整体http://biotrainee.cn:9013/graphics/plot_zoom_png?width=912&height=990
               n.epochs = 200,#200大数据，500小数据
               min.dist = 1,#距离
               #spread=1,#揉
               reduction.name = "umap.integrated",verbose = FALSE)
DimPlot(
  obj,
  reduction = "umap.integrated",label = T,
  group.by = c("etp_cl", "stim"),label.size = 5
)





 tail(rownames(bmpd5))
#Seurat to loomR####,
install.packages("hdf5r")
devtools::install_local("loomR.zip")
library(loomR)

sdata <- readRDS(file = "/Users/yuanzan/Desktop/tmp/seurat_project.rds")
# seurat对象转换为loop文件
sdata.loom <- as.loom(x = sdata, filename = "/Users/yuanzan/Desktop/tmp/sdata.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()


###monocle2####
###1.BMP4 part####
monocle=readRDS("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/BMP4/BMP4_D0-8.rds")


obj=readRDS("K5_PAGA_leiden_add_EpiSC_EpiLC.rds")

monocle=obj


expr=as(as.matrix(monocle@assays$RNA$counts),"sparseMatrix")
cellinfo=monocle@meta.data
f_data=data.frame(gene_short_name=rownames(monocle),row.names = rownames(monocle))

cellinfo=new("AnnotatedDataFrame",data=cellinfo)
f_data=new("AnnotatedDataFrame",data=f_data)
cds=newCellDataSet(expr,
                   phenoData = cellinfo,
                   featureData=f_data,
                   lowerDetectionLimit=0.5,
                   expressionFamily=negbinomial.size()#for raw counts
)


cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)

cds=detectGenes(cds,min_expr=0.1)
#print(head(fData(cds)))
expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))
express_genes=VariableFeatures(FindVariableFeatures(monocle,nfeatures = 1000))
express_genes=as.character(unlist(MarkerGenes_embryo_Genelist))
express_genes=intersect(expressed_genes,express_genes)

cds=setOrderingFilter(cds,express_genes)
print(plot_ordering_genes(cds))

cds = reduceDimension(cds,
                      #max_components=5,
                      method = "DDRTree")

cds=orderCells2(cds)


# GM_state <- function(cds){
#   if (length(unique(pData(cds)$State)) > 1){
#     T0_counts <- table(pData(cds)$State, pData(cds)$sample)[,"PNTD0"]
#     return(as.numeric(names(T0_counts)[which
#                                        (T0_counts == max(T0_counts))]))
#   } else {
#     return (1)
#   }
# }



cds <- orderCells(cds, root_state = GM_state(cds))

cds$leiden=as.character(cds$leiden)

cds@phenoData@data$TPRX1=monocle@assays$RNA$data["TPRX1",]

plot_cell_trajectory(cds,color_by="Klf5",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Kit",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Prdm1",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="invivo_celltype",size=1,show_backbone=T)


###2.K5 0-5D part####

obj=readRDS("K5_PAGA_leiden_add_EpiSC_EpiLC.rds")
obj=obj[,!obj$sample%in%c("EpiSC","EpiLC")]
table(obj$sample)

EpiSC=seurat_D0_5[,seurat_D0_5$sample=="EpiSCs"]
EpiSC$leiden=EpiSC$seurat_cluster
EpiSC=EpiSC[,EpiSC$leiden%in%c(4,13)]
EpiSC$leiden=gsub(4,"Epi-4",EpiSC$leiden)
EpiSC$leiden=gsub(13,"Epi-13",EpiSC$leiden)
table(EpiSC$leiden)

monocle=merge(obj,EpiSC)
monocle=JoinLayers(monocle)


expr=as(as.matrix(monocle@assays$RNA$counts),"sparseMatrix")
cellinfo=monocle@meta.data
f_data=data.frame(gene_short_name=rownames(monocle),row.names = rownames(monocle))

cellinfo=new("AnnotatedDataFrame",data=cellinfo)
f_data=new("AnnotatedDataFrame",data=f_data)
cds=newCellDataSet(expr,
                   phenoData = cellinfo,
                   featureData=f_data,
                   lowerDetectionLimit=0.5,
                   expressionFamily=negbinomial.size()#for raw counts
)

cds=detectGenes(cds,min_expr=0.1)
#print(head(fData(cds)))
expressed_genes=rownames(subset(fData(cds),num_cells_expressed >=10))
express_genes=VariableFeatures(FindVariableFeatures(monocle,nfeatures = 3000))
express_genes=as.character(unlist(MarkerGenes_embryo_Genelist))
express_genes=intersect(expressed_genes,express_genes)

cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)



cds=setOrderingFilter(cds,express_genes)
print(plot_ordering_genes(cds))

cds = reduceDimension(cds,
                      #max_components=5,
                      method = "DDRTree")

cds=orderCells2(cds)


# GM_state <- function(cds){
#   if (length(unique(pData(cds)$State)) > 1){
#     T0_counts <- table(pData(cds)$State, pData(cds)$sample)[,"EpiSCs"]
#     return(as.numeric(names(T0_counts)[which
#                                        (T0_counts == max(T0_counts))]))
#   } else {
#     return (1)
#   }
# }



cds <- orderCells(cds, root_state = GM_state(cds))

cds$leiden=as.character(cds$leiden)


cds@phenoData@data$TPRX1=monocle@assays$RNA$data["TPRX1",]

plot_cell_trajectory(cds,color_by="Klf5",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Kit",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="Prdm1",size=1,show_backbone=T)
plot_cell_trajectory(cds,color_by="invivo_celltype",size=1,show_backbone=T)



setwd("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/")
.libPaths("/home/data/t020406/R/x86_64-pc-linux-gnu-library/4.4/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

seurat_D0_5=readRDS("K5_D0-5.rds")

obj=readRDS("K5_PAGA_leiden_add_EpiSC_EpiLC.rds")
obj=obj[,!obj$sample%in%c("EpiSC","EpiLC")]
table(obj$sample)

EpiSC=seurat_D0_5[,seurat_D0_5$sample=="EpiSCs"]
EpiSC$leiden=EpiSC$seurat_cluster
EpiSC=EpiSC[,EpiSC$leiden%in%c(4,13)]
EpiSC$leiden=gsub(4,"Epi-4",EpiSC$leiden)
EpiSC$leiden=gsub(13,"Epi-13",EpiSC$leiden)
table(EpiSC$leiden)

seurat_D0_5=NULL

monocle=merge(obj,EpiSC)
monocle=JoinLayers(monocle)
monocle=NormalizeData(monocle)
monocle=ScaleData(monocle)


result_df=readRDS("result_df.rds")
monocle$pseudotime=result_df$pseudotime_final
  
# Standard_process=
#   function(obj,nfeatures=2000,resolution=1){
#     obj <- NormalizeData(obj)
#     obj <- ScaleData(obj)
#     obj=FindVariableFeatures(obj,nfeatures = nfeatures)
#     obj <- RunPCA(obj)
#     
#     
#     
#     obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
#     obj <- FindClusters(obj, resolution =resolution,random.seed =1, cluster.name = c("seurat_clusters"))
#     #obj <- FindClusters(obj, resolution = 0.65,random.seed =1, cluster.name = c("seurat_clusters"))
#     
#     
#     obj <- RunUMAP(obj, reduction = "pca",
#                    dims = 1:30, seed.use =1,
#                    n.neighbors = 20, #20更多粘连和小块,50偏整体http://biotrainee.cn:9013/graphics/plot_zoom_png?width=912&height=990
#                    n.epochs = 200,#200大数据，500小数据
#                    #min.dist = 1,#距离
#                    #spread=1,#揉
#                    #return.model = TRUE,
#                    reduction.name = "umap",verbose = FALSE)
# 
#     
#     return(obj)
#   }

monocle=Standard_process(monocle)


options(future.globals.maxSize = 100 * 1024^3) # 10GB
monocle_raw[["percent.mt"]] <- PercentageFeatureSet(monocle_raw, pattern = "^mt-")
monocle_raw <- SCTransform(monocle_raw, vars.to.regress = "percent.mt", verbose = FALSE
                           # ,variable.features.n = 2000
)



####1.Nanog use gene####
ETP=read.csv("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/ETP_gene_WGCNA_ver3_4431_.csv")
head(ETP)
dim(ETP)
table(ETP$Cluster)


esc_relative=ETP[ETP$Cluster=="ESC",1]
tsc_relative=ETP[ETP$Cluster=="TSC",1]
pgclc_relative=ETP[ETP$Cluster=="PGCLC",1]
esc_pgc_relative=ETP[ETP$Cluster%in%c("ESC&PGCLC","ESC","PGCLC"),1]

PGC_all=c("Dppa3", "Prdm1", "Tfap2c", "Nanos3","Id3", "Rhox6", "Irf1", "Sox15", "Rhox9", "Zic3", "Kit", 
          "Rbm47", "Phf19", "Msx1",  "Msx2","Dmrt1", "Dnd1", "Prdm14", "Dazl", "Ddx4")

TSC_all=c("Cdx2","Gata3","Tfap2c","Rhox5","Fgfr2","Gata2","Krt18","Ddah1","Rin3","Plet1","Hand1","Elf5",
          "Tead4","Eomes","Ets2","Yap1","Erf","Ly6a")

ESC_all=c("Klf4","Tbx3","Zfp42","Klf2","Tfcp2l1","Dppa5a","Esrrb","Prdm14","Nanog","Nr5a2",
          "Gbx2","Pbx1","Tcf15","Tcf3","Tgif1")

monocle=AddModuleScore(monocle,features = list(ESC_all),
                       name = "ESC_all")
monocle=AddModuleScore(monocle,features = list(TSC_all),
                       name = "TSC_all")
monocle=AddModuleScore(monocle,features = list(c(PGC_all)),
                       name = "PGC_all")
monocle=AddModuleScore(monocle,features = list(c(esc_pgc_relative)),
                       name = "ESC_and_PGCLC_AUC")



####2.Klf5 active ETP gene####

k5_a=read.table("K5_D2_L_bind_active_anno.bed")
k5_a=k5_a[k5_a$V5=="Bind","V4"]


monocle=AddModuleScore(monocle,features = list(c(k5_a)),
                       name = "K5_ac_2_bind")


####3.Vlnplot in EpiSC D1,2####
monocle_sub=monocle[,monocle$leiden%in%c("Epi-13","Epi-4",
                                         40,24,21,28)]


c(20,25,39)
c(10,8,14,18,2,3,12,31,29,1,33,
17,22,27,26,30,36)

p1=FeaturePlot(monocle,features = c("pseudotime"),order = T)
p2=DimPlot(monocle,reduction = "umap", group.by = "sample", label = TRUE,pt.size = .1)
p1+p2

VlnPlot(monocle, features = c("pseudotime"),group.by = "sample",pt.size=0)
# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)


monocle_raw=monocle
monocle=monocle_sub
monocle_sub=AddModuleScore(monocle_sub,features = list(c(deg_res[deg_res$category=="Cluster 13 up","gene"])),
                       name = "Cluster 13 up")

# 设置参数
data_type <- "Cluster 13 up1"  # 可以改为 "TSC_AUC", "Nanog", "Klf5" 等

# 获取数据 - 区分基因表达和AUC值
if (data_type %in% colnames(monocle@meta.data)) {
  values <- monocle_sub[[data_type]][,1]
  y_label <- paste(data_type, "Score")
} else {
  # 获取基因表达数据
  if (data_type %in% rownames(monocle_sub)) {
    values <- GetAssayData(monocle_sub, slot = "data")[data_type, ]
  } else {
    stop(paste("未找到基因", data_type, "的表达数据"))
  }
  y_label <- paste(data_type, "Expression Level")
}

# 创建绘图数据框
plot_data <- data.frame(
  value = values,
  leiden_original = as.character(monocle_sub$leiden),
  stringsAsFactors = FALSE
)

# 重新标记leiden聚类
plot_data$leiden <- plot_data$leiden_original
plot_data$leiden[plot_data$leiden_original %in% c("20", "25", "39")] <- "TSC"
plot_data$leiden[plot_data$leiden_original %in% c("10", "8", "14", "18", "2", "3", "12", "31", "29", "1", "33", "17", "22", "27", "26", "30", "36")] <- "NR"

# 排除35
plot_data <- plot_data[plot_data$leiden != "35", ]

# 定义分组规则（排除35后的分组）
group_mapping <- list(
  "Group1" = c("Epi-13", "Epi-4"),
  "Group2" = c("40", "24"),
  "Group3" = c("21", "28")
)

# 创建分组信息
plot_data$group <- NA
plot_data$x_position <- NA

# 应用分组规则
for (i in seq_along(group_mapping)) {
  group_name <- names(group_mapping)[i]
  clusters <- group_mapping[[group_name]]
  
  for (cluster in clusters) {
    plot_data$group[plot_data$leiden == cluster] <- group_name
    plot_data$x_position[plot_data$leiden == cluster] <- i
  }
}

# 移除未分组的细胞
plot_data <- plot_data[!is.na(plot_data$group), ]

# 设置因子水平以确保正确的顺序
plot_data$group <- factor(plot_data$group, levels = names(group_mapping))
plot_data$leiden <- factor(plot_data$leiden, levels = unique(unlist(group_mapping)))
plot_data$x_position <- as.factor(plot_data$x_position)

# 创建自定义的x轴标签
x_labels <- c("EpiSCs", "Day1", "Day2")

# 手动定义颜色
manual_colors <- c(
  "Epi-13" = "#A1D99B",
  "Epi-4" = "#FDAE6B", 
  "40" = "#41AB5D",
  "24" = "#F16913",
  "21" = "#00441B",
  "28" = "#A50F15"
)

# 只保留实际存在的leiden颜色
leiden_levels <- levels(plot_data$leiden)
fill_colors <- manual_colors[names(manual_colors) %in% leiden_levels]

# 创建对应的深色边框颜色
border_colors <- sapply(fill_colors, function(color) {
  rgb_val <- col2rgb(color)
  dark_rgb <- rgb_val * 0.6
  rgb(dark_rgb[1], dark_rgb[2], dark_rgb[3], maxColorValue = 255)
})

# 修复：使用基础R方法计算平均值
mean_data <- data.frame()
for (pos in levels(plot_data$x_position)) {
  for (cluster in levels(plot_data$leiden)) {
    subset_data <- plot_data[plot_data$x_position == pos & plot_data$leiden == cluster, ]
    if (nrow(subset_data) > 0) {
      mean_val <- mean(subset_data$value, na.rm = TRUE)
      mean_data <- rbind(mean_data, 
                         data.frame(x_position = pos, 
                                    leiden = cluster, 
                                    mean_value = mean_val))
    }
  }
}

# 确保mean_data的x_position是因子
mean_data$x_position <- factor(mean_data$x_position, levels = levels(plot_data$x_position))

# 修复：手动计算每组内的p值
pvalue_results <- data.frame()
for (pos in levels(plot_data$x_position)) {
  subset_data <- plot_data[plot_data$x_position == pos, ]
  clusters <- unique(subset_data$leiden)
  
  if (length(clusters) == 2) {
    values1 <- subset_data$value[subset_data$leiden == clusters[1]]
    values2 <- subset_data$value[subset_data$leiden == clusters[2]]
    
    if (length(values1) > 1 & length(values2) > 1) {
      p_val <- wilcox.test(values1, values2)$p.value
      p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))
    } else {
      p_val <- NA
      p_label <- ""
    }
  } else {
    p_val <- NA
    p_label <- ""
  }
  
  pvalue_results <- rbind(pvalue_results, 
                          data.frame(x_position = pos, 
                                     p_value = p_val, 
                                     p_label = p_label))
}

# 计算y轴位置
y_max <- max(plot_data$value, na.rm = TRUE)

# 检查数据结构
print("修复后的 mean_data 结构:")
print(str(mean_data))
print("mean_data 内容:")
print(mean_data)

# 绘制小提琴图
p_optimized <- ggplot(plot_data, aes(x = x_position, y = value)) +
  # 小提琴图
  geom_violin(aes(fill = leiden), scale = "width", width = 0.7, 
              alpha = 0.4, position = position_dodge(width = 0.8), size = 0) +
  geom_violin(aes(color = leiden), scale = "width", width = 0.7, 
              alpha = 0, position = position_dodge(width = 0.8), size = 0.8) +
  # 添加平均值点
  geom_point(data = mean_data, 
             aes(x = x_position, y = mean_value, color = leiden),
             shape = 18, size = 3.5, 
             position = position_dodge(width = 0.8)) +
  # 添加leiden cluster名字标签
  geom_text(data = mean_data,
            aes(x = x_position, y = mean_value, label = leiden, color = leiden),
            size = 3.5, fontface = "bold", vjust = -1.2,
            position = position_dodge(width = 0.8), show.legend = FALSE) +
  # 添加p值标注
  geom_segment(data = pvalue_results,
               aes(x = as.numeric(x_position) - 0.4, 
                   xend = as.numeric(x_position) + 0.4,
                   y = y_max * 0.9, yend = y_max * 0.9),
               color = "black", size = 0.8) +
  geom_text(data = pvalue_results,
            aes(x = as.numeric(x_position), y = y_max * 1, label = p_label),
            size = 4, fontface = "bold", color = "red") +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = fill_colors, name = "Leiden Cluster") +
  scale_color_manual(values = border_colors, name = "Leiden Cluster") +
  labs(
    title = paste(data_type, "by Cluster Groups"),
    x = "Cluster Groups",
    y = y_label
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

print(p_optimized)
#6x3.5

# 假设比较cluster A和cluster B的某个基因表达
cluster_A_cells <- WhichCells(monocle, idents = "Epi-13")
cluster_B_cells <- WhichCells(monocle, idents = "Epi-4")

# 获取基因表达值
gene_expression <- monocle$PGC_all1

# Wilcoxon检验（推荐）
p_value <- wilcox.test(
  gene_expression[cluster_A_cells],
  gene_expression[cluster_B_cells]
)$p.value
p_value
# 保存图形
# ggsave(paste0("violin_", data_type, ".png"), p_optimized, 
#        width = 12, height = 7, dpi = 300)

# 显示统计信息
cat("Data summary by group and leiden:\n")
data_summary <- plot_data %>%
  group_by(group, leiden) %>%
  summarise(
    n_cells = n(),
    mean_value = mean(value),
    median_value = median(value),
    sd_value = sd(value),
    .groups = 'drop'
  ) %>%
  arrange(group, leiden)

print(data_summary)

# 显示重新标记的汇总信息
cat("\nReclustering summary:\n")
recluster_summary <- plot_data %>%
  group_by(leiden_original, leiden) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  arrange(leiden, leiden_original)

print(recluster_summary)









###Monocle3####

monocle_meta=monocle@meta.data

cds <- new_cell_data_set(
  expression_data = GetAssayData(monocle, assay = "RNA", slot = "counts"),
  cell_metadata = monocle@meta.data,
  gene_metadata = data.frame(
    gene_short_name = rownames(monocle),
    row.names = rownames(monocle)
  )
)
# 1. 数据预处理
cds <- preprocess_cds(cds, num_dim = 50)

# 2. 降维和可视化（使用UMAP）
cds <- reduce_dimension(cds,
                        umap.min_dist = 1,
                        umap.n_neighbors = 50)
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "sample")
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "leiden")
plot_cells(cds, 
           genes = c("Nanog", "Klf5"),
           label_cell_groups = FALSE,
           show_trajectory_graph = TRUE)

# 3. 细胞聚类
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")


# 4. 学习细胞轨迹
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "sample", 
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE, 
           label_branch_points = FALSE)


# 方法1：通过交互式选择起点
#cds <- order_cells(cds)  # 会弹出交互界面选择根节点

# 方法2：编程指定起点（推荐）
get_earliest_principal_node <- function(cds, cell_type="D1_K5"){
  cell_ids <- which(colData(cds)$sample == cell_type)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}

root_node <- get_earliest_principal_node(cds, "D1_K5")
cds <- order_cells(cds, root_pr_nodes = root_node)
pseudotime_values1=as.data.frame(pseudotime(cds))

root_node <- get_earliest_principal_node(cds, "EpiSCs")
cds <- order_cells(cds, root_pr_nodes = root_node)
pseudotime_values2=as.data.frame(pseudotime(cds))




pseudotime_values2[,1]=pseudotime_values2[,1]*0.65
pv2_max=max(pseudotime_values2[!is.infinite(pseudotime_values2[,1]),])
pseudotime_values1[,1]=pseudotime_values1[,1]+pv2_max



# 假设您的两个data.frame叫做pseudotime_values1和pseudotime_values2

# 1. 确保两个data.frame有相同的行名
common_cells <- intersect(rownames(pseudotime_values1), rownames(pseudotime_values2))
pseudotime_values1 <- pseudotime_values1[common_cells, , drop = FALSE]
pseudotime_values2 <- pseudotime_values2[common_cells, , drop = FALSE]

# 2. 创建一个新的data.frame来存储合并后的结果
merged_pseudotime <- data.frame(
  cell_id = common_cells,
  pseudotime_final = NA,
  stringsAsFactors = FALSE
)

# 3. 定义处理函数
replace_inf_values <- function(df1, df2) {
  # 提取数值向量（假设每列只有一个数值列）
  values1 <- df1[, 1]
  values2 <- df2[, 1]
  
  # 计算正常值的中位数（排除Inf）
  normal_values1 <- values1[!is.infinite(values1)]
  normal_values2 <- values2[!is.infinite(values2)]
  all_normal_values <- c(normal_values1, normal_values2)
  
  if (length(all_normal_values) > 0) {
    median_value <- median(all_normal_values, na.rm = TRUE)
  } else {
    median_value <- 0  # 如果所有值都是Inf，使用0
  }
  
  # 创建结果向量
  result <- rep(NA, nrow(df1))
  
  for (i in 1:nrow(df1)) {
    val1 <- values1[i]
    val2 <- values2[i]
    
    if (!is.infinite(val1)) {
      # 如果df1有正常值，优先使用
      result[i] <- val1
    } else if (!is.infinite(val2)) {
      # 如果df2有正常值，使用df2的值
      result[i] <- val2
    } else {
      # 如果两边都是Inf，使用中位数
      result[i] <- median_value
    }
  }
  
  return(result)
}

# 4. 应用处理函数
final_pseudotime <- replace_inf_values(pseudotime_values1, pseudotime_values2)

# 5. 创建最终的data.frame
result_df <- data.frame(
  cell_id = common_cells,
  pseudotime_final = final_pseudotime,
  stringsAsFactors = FALSE
)

# 设置行名
rownames(result_df) <- common_cells

# 6. 查看结果
head(result_df)
tail(result_df)
summary(result_df)


result_df=readRDS("result_df.rds")

saveRDS(result_df,"result_df.rds")
save_monocle_objects(cds,"/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/")


# 可视化伪时间
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5)





# 7. 统计替换情况
inf_count_df1 <- sum(is.infinite(pseudotime_values1[, 1]))
inf_count_df2 <- sum(is.infinite(pseudotime_values2[, 1]))
replaced_count <- sum(is.infinite(pseudotime_values1[, 1]) & is.infinite(pseudotime_values2[, 1]))

cat("数据统计:\n")
cat("DF1中Inf值的数量:", inf_count_df1, "\n")
cat("DF2中Inf值的数量:", inf_count_df2, "\n")
cat("被中位数替换的细胞数量:", replaced_count, "\n")
cat("最终有效伪时间值的数量:", sum(!is.infinite(result_df$pseudotime_final)), "\n")



# 将伪时间添加回Seurat对象
seurat_obj$pseudotime <- pseudotime_values


#Seurat to h5ad####
library(data.table)
library(Seurat)
library(loomR)
library(SeuratObject)
library(ggplot2)
library(plyr)
library(limma)
library(AUCell)
#library(RcisTarget)
#library(GENIE3)
#library(SCENIC)
library(doRNG)
library(SeuratData)
library(patchwork)
#library(monocle)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)
library(harmony)
library(cowplot)
library(Seurat)
#library(tidyverse)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
#library(scCustomize)
library(harmony)
library(corrplot)
library("SeuratDisk")
setwd("/home/data/t020406/iK5_SGC_PNT_D0-5_iK24_D5_scRNAseq/")
seurat=readRDS("monocle_for_PAGA.rds")
episc=readRDS("EpiSCs_QC.rds")

episc=NormalizeData(episc)
episc=ScaleData(episc)
episc=FindVariableFeatures(episc,nfeatures = 2000)
episc <- RunPCA(episc)

episc <- FindNeighbors(episc, dims = 1:30, reduction = "pca")
episc <- FindClusters(episc,resolution = 0.1)

Klfs=c("Klf2","Klf5")
ESC_Marker=c("Nanog","Tfcp2l1","Zfp42","Klf4","Tbx3")
PGC_Marker=c("Dppa3","Prdm1","Nanos3","Sox15","Id3")
TSC_and_ESC=c("Esrrb","Socs3","Tcf7l1","Zbtb18")
ESC_and_PGC=c("Pou5f1","Dppa5a","Nr5a2","Prdm14")
PGC_and_TSC=c("Tfap2c","Rhox6","Rhox9","Msx2","Msx1")
TSC_Marker=c("Gata3","Cdx2","Krt7","Gata2","Rin3","Plet1","Hand1","Elf5","Bmp7","Krt19","Ascl2")

FeaturePlot(episc,features = c(Klfs,ESC_Marker,PGC_Marker,TSC_and_ESC,
                               ESC_and_PGC,PGC_and_TSC,TSC_Marker),order = T)
FeaturePlot(episc,features = c("Nanog","Bmp7","Krt19","Ascl2"),order = T)
Idents(episc)="seurat_clusters"
episc_deg <- FindAllMarkers(episc, only.pos = TRUE, logfc.threshold = 0.25)


episc_leiden=as.numeric(episc$seurat_clusters)
table(episc_leiden)
episc$leiden=episc_leiden+41
table(episc$leiden)
filter_gene=
  function(deg_res){
    deg_res_filter=deg_res[deg_res$p_val_adj<0.01,]
    deg_res_filter=deg_res_filter[deg_res_filter$avg_log2FC>1,]
    deg_res_filter=deg_res_filter[order(deg_res_filter$avg_log2FC,decreasing = T),]
    return(deg_res_filter)
  }
episc_deg=filter_gene(episc_deg)

write.csv(episc_deg,"episc_deg.csv")

seurat_episc=merge(seurat,episc)

seurat_episc=NormalizeData(seurat_episc)
seurat_episc=ScaleData(seurat_episc)
seurat_episc=FindVariableFeatures(seurat_episc,nfeatures = 2000)
seurat_episc <- RunPCA(seurat_episc)



table(seurat_episc$leiden)

seurat_episc@meta.data=seurat_episc@meta.data[,c(1:10,13,21:24,89:95,117)]
colnames(seurat_episc@meta.data)
table(seurat_episc$leiden)

SaveH5Seurat(seurat_episc, filename = "K5_PAGA_leiden_add_EpiSC.h5Seurat")
Convert("K5_PAGA_leiden_add_EpiSC.h5Seurat", dest = "h5ad")
##Stringency_data_strat_plot####






#Other function####
Plot_Pie=
  function(df,title="Pieplot",color=brewer.pal(2,"Set2"),
           show_legend=F,Only_show_precentage=F){
    df$condition=factor(df$condition,levels = rev(df$condition))
    
    myLabel = as.vector(df$condition)   
    myLabel = paste(myLabel, "(", round(df$num / sum(df$num) * 100, 2), "%)", sep = "")   
    
    if(Only_show_precentage==T){
      myLabel = as.vector(df$condition)   
      myLabel = paste(round(df$num / sum(df$num) * 100, 2), "%", sep = "")   
      
    }
    
    #if(length(df$condition)>2){color=brewer.pal(length(df$condition),"Set2")[c(1,3,2)]}
    
    title=paste(title,"(n=",sum(df$num),")",sep = "")
    
    p=ggplot(df,aes(x="",y=num,fill=condition))+
      geom_bar(stat = "identity",width = 1)+
      coord_polar(theta = "y",start = 0)+
      scale_fill_manual(values = color)+
      theme_bw()+
      labs(x = "", y = "", title = title) + 
      theme(panel.grid=element_blank(),## 去掉白色外框
            axis.text=element_blank(),## 把图旁边的标签去掉
            axis.ticks=element_blank(),## 去掉左上角的坐标刻度线
            panel.border=element_blank()## 去掉最外层的正方形边框
      ) +
    #geom_text(aes(y = num/2 + c(0, cumsum(num)[-length(num)]), x = 1, label = myLabel), size = 4)
    geom_label_repel(aes(y = num/2+ c(0, cumsum(num)[-length(num)]), x = 1,label = myLabel), 
                     size=4, show.legend = F,fill=NA,label.padding = 0.7, nudge_x = 1)
    if(show_legend==F){
      p=p+guides(fill=FALSE)
    }
    return(p)
  }




library(Rcpp)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  IntegerMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
' )

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}



#normalize pseudotime
max(l10[,1])-min(l10[,1])

study_gene=c("Prdm1","Nanog","Klf4","Gata3","Sox15")

make_lineage_long_df=
  function(pancreas_sub,l_df,study_gene,lineage_name){
l_exp=as.data.frame(t(pancreas_sub@assays$RNA$data[study_gene,colnames(pancreas_sub)%in%rownames(l_df)]))
l_exp$cellname=rownames(l_exp)

l_meta=data.frame(pstime=l_df[,1],
              cellname=rownames(l_df),
              order=paste(lineage_name,1:dim(l_df)[1],sep = "_"),
              lineage=lineage_name)
l_meta$pstime=c(l_meta$pstime-min(l_meta$pstime))/c(max(l_meta$pstime)-min(l_meta$pstime))
meta_data=pancreas_sub@meta.data[,c("sample","seurat_clusters")]
meta_data$cellname=rownames(meta_data)
l_fin=plyr::join(l_meta,meta_data)
l_fin=plyr::join(l_fin,l_exp)

return(l_fin)
}




l_df=pancreas_sub@tools$DynamicFeatures_Lineage4$upr_matrix[,1:2]
l4=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l4")

l_df=pancreas_sub@tools$DynamicFeatures_Lineage5$upr_matrix[,1:2]
l5=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l5")

l_df=pancreas_sub@tools$DynamicFeatures_Lineage15$upr_matrix[,1:2]
l15=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l15")

l_df=pancreas_sub@tools$DynamicFeatures_Lineage12$upr_matrix[,1:2]
l12=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l12")

l_df=pancreas_sub@tools$DynamicFeatures_Lineage10$upr_matrix[,1:2]
l10=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l10")

l_df=pancreas_sub@tools$DynamicFeatures_Lineage7$upr_matrix[,1:2]
l7=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l7")

l_df=pancreas_sub@tools$DynamicFeatures_Lineage8$upr_matrix[,1:2]
l8=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l8")

l_df=pancreas_sub@tools$DynamicFeatures_Lineage2$upr_matrix[,1:2]
l2=make_lineage_long_df(pancreas_sub,l_df,study_gene,"l2")




df=rbind(l4,l5,l15,l12,l10,l7,l8,l2)
df[,7:dim(df)[2]]=scale(df[,7:dim(df)[2]])
head(df)

smooth_plot=function(df){}
ggplot(df) +
  #geom_point(aes(x = pstime, y = Prdm1,color=sample),alpha = .5,shape=20) +
  geom_smooth(aes(pstime, Nanog,color=lineage), method="loess",se = T)+
  scale_color_manual(values=c(#palette_list$Paired[1:6],
                              c( "#6699CC","#661100",
                                                          "#CC6677","#332288" ,"#117733" ,
                                                          "#999933","#AA4499","#888888")))+
  theme_bw()+
  coord_cartesian(ylim=c(-2.5,2.5))+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10))+
  ggtitle("Prdm1")
  #scale_y_continuous(limits=c(0,2))+

options(repr.plot.width =8, repr.plot.height =3, repr.plot.res = 100)
p4



pseudotime_heatmap=function(df){}
pseudotime<-cds$Pseudotime
newdata<-data.frame(Pseudotime = seq(min(pseudotime), 
                                     max(pseudotime), length.out = 100))
gene1<-row.names(subset(BEAM_res, qval==0))
genSmoothCurves_mat1<-genSmoothCurves(cds[gene1,],
                                      new_data = newdata,
                                      cores = 1)
data1<-log10(genSmoothCurves_mat1+1)
data1

hmcols=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
bks <- seq(-3.1, 3.1, length.out = length(hmcols))
p<-pheatmap(data1,show_rownames = T, scale = "row", 
            cluster_cols = F, cluster_rows = T,gaps_col=100,
            border_color=NA,#annotation_row = annotation_row ,
            clustering_distance_rows = "correlation",#"correlation，euclidean，maximum，manhattan，canberra，binary，minkowski"
            clustering_method =  'mcquitty', #'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'
            cutree_row=3,
            breaks=bks,
            border=T,
            color = hmcols)





Btable_insert_into_Atable=
  function(A,B,index_A="rownames",index_B="rownames",C="Merge"){
    
    
    A=as.data.frame(A)
    B=as.data.frame(B)
    
    ###先生成A中给B表插入的位置
    newcol_1=data.frame(c(rep(NA,length(A[,1]))))
    if(length(B)!=1){
      for(col in 1:(length(B)-1)){
        newcol=data.frame(c(rep(NA,length(A[,1]))))
        newcol_1=cbind(newcol_1,newcol)
      }
    }
    ###把生成的表合并到A中
    colnames(newcol_1)=colnames(B)
    A_B=cbind(A,newcol_1)
    if(index_A=="rownames"){index_A=rownames(A)}else{index_A=A_B[,index_A]}
    if(index_B=="rownames"){index_B=rownames(B)}else{index_B=B[,index_B]}
    ###开始把B表的数据插入刚才的空白位置中去
    for(row in 1:length(B[,1])){
      if(sum(as.character(index_B)[row]%in%as.character(index_A))==1){
        A_B[which(as.character(index_A)==as.character(index_B)[row]),(length(A)+1):(length(A)+length(B))]=as.character(B[row,])}
    }
    B=A_B[,(length(A[1,])+1):length(A_B)]
    if(C=="Merge"){
      return(A_B)
    }
    if(C=="B"){
      return(B)
    }
    
  }




DynamicHeatmap_modify=
  function (srt, lineages, features = NULL, study_gene_set=NULL,use_fitted = FALSE, 
          border = TRUE, flip = FALSE, min_expcells = 20, r.sq = 0.2, 
          dev.expl = 0.2, padjust = 0.05, num_intersections = NULL, 
          cell_density = 1, cell_bins = 100, order_by = c("peaktime", 
                                                          "valleytime"), slot = "counts", assay = NULL, exp_method = c("zscore", 
                                                                                                                       "raw", "fc", "log2fc", "log1p"), exp_legend_title = NULL, 
          limits = NULL, lib_normalize = identical(slot, "counts"), 
          libsize = NULL, family = NULL, cluster_features_by = NULL, 
          cluster_rows = FALSE, cluster_row_slices = FALSE, cluster_columns = FALSE, 
          cluster_column_slices = FALSE, show_row_names = FALSE, show_column_names = FALSE, 
          row_names_side = ifelse(flip, "left", "right"), column_names_side = ifelse(flip, 
                                                                                     "bottom", "top"), row_names_rot = 0, column_names_rot = 90, 
          row_title = NULL, column_title = NULL, row_title_side = "left", 
          column_title_side = "top", row_title_rot = 0, column_title_rot = ifelse(flip, 
                                                                                  90, 0), feature_split = NULL, feature_split_by = NULL, 
          n_split = NULL, split_order = NULL, split_method = c("mfuzz", 
                                                               "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime"), 
          decreasing = FALSE, fuzzification = NULL, anno_terms = FALSE, 
          anno_keys = FALSE, anno_features = FALSE, terms_width = unit(4, 
                                                                       "in"), terms_fontsize = 8, keys_width = unit(2, "in"), 
          keys_fontsize = c(6, 10), features_width = unit(2, "in"), 
          features_fontsize = c(6, 10), IDtype = "symbol", species = "Homo_sapiens", 
          db_update = FALSE, db_version = "latest", db_combine = FALSE, 
          convert_species = FALSE, Ensembl_version = 103, mirror = NULL, 
          db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, 
          maxGSSize = 500, GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", 
          simplify_method = "Wang", simplify_similarityCutoff = 0.7, 
          pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, 
          show_termid = FALSE, topWord = 20, words_excluded = NULL, 
          nlabel = 20, features_label = NULL, label_size = 10, label_color = "black", 
          pseudotime_label = NULL, pseudotime_label_color = "black", 
          pseudotime_label_linetype = 2, pseudotime_label_linewidth = 3, 
          heatmap_palette = "viridis", heatmap_palcolor = NULL, pseudotime_palette = "cividis", 
          pseudotime_palcolor = NULL, feature_split_palette = "simspec", 
          feature_split_palcolor = NULL, cell_annotation = NULL, cell_annotation_palette = "Paired", 
          cell_annotation_palcolor = NULL, cell_annotation_params = if (flip) list(width = unit(5, 
                                                                                                "mm")) else list(height = unit(5, "mm")), feature_annotation = NULL, 
          feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, 
          feature_annotation_params = if (flip) list(height = unit(5, 
                                                                   "mm")) else list(width = unit(5, "mm")), separate_annotation = NULL, 
          separate_annotation_palette = "Paired", separate_annotation_palcolor = NULL, 
          separate_annotation_params = if (flip) list(width = unit(10, 
                                                                   "mm")) else list(height = unit(10, "mm")), reverse_ht = NULL, 
          use_raster = NULL, raster_device = "png", raster_by_magick = FALSE, 
          height = NULL, width = NULL, units = "inch", seed = 11, 
          ht_params = list()) {
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_R("magick")
  }
  split_method <- match.arg(split_method)
  order_by <- match.arg(order_by)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", 
                      ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(as.character(x = formals()$exp_method), 
                       "(", data_nm, ")")
  }
  else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }
  exp_name <- exp_legend_title %||% exp_name
  assay <- assay %||% DefaultAssay(srt)
  if (any(!lineages %in% colnames(srt@meta.data))) {
    lineages_missing <- lineages[!lineages %in% colnames(srt@meta.data)]
    for (l in lineages_missing) {
      if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
        pseudotime <- srt@tools[[paste0("DynamicFeatures_", 
                                        l)]][["lineages"]]
        srt@meta.data[[l]] <- srt@meta.data[[pseudotime]]
      }
      else {
        stop("lineages: ", l, " is not in the meta data of the Seurat object")
      }
    }
  }
  if (is.null(feature_split_by)) {
    feature_split_by <- lineages
  }
  if (any(!feature_split_by %in% lineages)) {
    stop("'feature_split_by' must be a subset of lineages.")
  }
  if (!split_method %in% c("mfuzz", "kmeans", "kmeans-peaktime", 
                           "hclust", "hclust-peaktime")) {
    stop("'split_method' must be one of 'mfuzz', 'kmeans', 'kmeans-peaktime', 'hclust', 'hclust-peaktime'.")
  }
  if (!is.null(feature_split) && is.null(names(feature_split))) {
    stop("'feature_split' must be named.")
  }
  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (is.numeric(pseudotime_label)) {
    if (length(pseudotime_label_color) == 1) {
      pseudotime_label_color <- rep(pseudotime_label_color, 
                                    length(pseudotime_label))
    }
    if (length(pseudotime_label_linetype) == 1) {
      pseudotime_label_linetype <- rep(pseudotime_label_linetype, 
                                       length(pseudotime_label))
    }
    if (length(pseudotime_label_linewidth) == 1) {
      pseudotime_label_linewidth <- rep(pseudotime_label_linewidth, 
                                        length(pseudotime_label))
    }
    npal <- unique(c(length(pseudotime_label), length(pseudotime_label_color), 
                     length(pseudotime_label_linetype), length(pseudotime_label_linewidth)))
    if (length(npal[npal != 0]) > 1) {
      stop("Parameters for the pseudotime_label must be the same length!")
    }
  }
  if (!is.null(cell_annotation)) {
    if (length(cell_annotation_palette) == 1) {
      cell_annotation_palette <- rep(cell_annotation_palette, 
                                     length(cell_annotation))
    }
    if (length(cell_annotation_palcolor) == 1) {
      cell_annotation_palcolor <- rep(cell_annotation_palcolor, 
                                      length(cell_annotation))
    }
    npal <- unique(c(length(cell_annotation_palette), length(cell_annotation_palcolor), 
                     length(cell_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("cell_annotation_palette and cell_annotation_palcolor must be the same length as cell_annotation")
    }
    if (any(!cell_annotation %in% c(colnames(srt@meta.data), 
                                    rownames(srt@assays[[assay]])))) {
      stop("cell_annotation: ", paste0(cell_annotation[!cell_annotation %in% 
                                                         c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))], 
                                       collapse = ","), " is not in the Seurat object.")
    }
  }
  if (!is.null(feature_annotation)) {
    if (length(feature_annotation_palette) == 1) {
      feature_annotation_palette <- rep(feature_annotation_palette, 
                                        length(feature_annotation))
    }
    if (length(feature_annotation_palcolor) == 1) {
      feature_annotation_palcolor <- rep(feature_annotation_palcolor, 
                                         length(feature_annotation))
    }
    npal <- unique(c(length(feature_annotation_palette), 
                     length(feature_annotation_palcolor), length(feature_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("feature_annotation_palette and feature_annotation_palcolor must be the same length as feature_annotation")
    }
    if (any(!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features))) {
      stop("feature_annotation: ", paste0(feature_annotation[!feature_annotation %in% 
                                                               colnames(srt@assays[[assay]]@meta.features)], 
                                          collapse = ","), " is not in the meta data of the ", 
           assay, " assay in the Seurat object.")
    }
  }
  if (!is.null(separate_annotation)) {
    if (length(separate_annotation_palette) == 1) {
      separate_annotation_palette <- rep(separate_annotation_palette, 
                                         length(separate_annotation))
    }
    if (length(separate_annotation_palcolor) == 1) {
      separate_annotation_palcolor <- rep(separate_annotation_palcolor, 
                                          length(separate_annotation))
    }
    npal <- unique(c(length(separate_annotation_palette), 
                     length(separate_annotation_palcolor), length(separate_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("separate_annotation_palette and separate_annotation_palcolor must be the same length as separate_annotation")
    }
    if (any(!unique(unlist(separate_annotation)) %in% c(colnames(srt@meta.data), 
                                                        rownames(srt@assays[[assay]])))) {
      stop("separate_annotation: ", paste0(unique(unlist(separate_annotation))[!unique(unlist(separate_annotation)) %in% 
                                                                                 c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))], 
                                           collapse = ","), " is not in the Seurat object.")
    }
  }
  if (length(width) == 1) {
    width <- rep(width, length(lineages))
  }
  if (length(height) == 1) {
    height <- rep(height, length(lineages))
  }
  if (length(width) >= 1) {
    names(width) <- lineages
  }
  if (length(height) >= 1) {
    names(height) <- lineages
  }
  column_split <- NULL
  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_row_slices_raw <- cluster_row_slices
    cluster_column_slices_raw <- cluster_column_slices
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
    cluster_row_slices <- cluster_column_slices_raw
    cluster_column_slices <- cluster_row_slices_raw
  }
  cell_union <- unique(colnames(srt@assays[[1]])[apply(srt@meta.data[, 
                                                                     lineages, drop = FALSE], 1, function(x) !all(is.na(x)))])
  Pseudotime_assign <- rowMeans(srt@meta.data[cell_union, 
                                              lineages, drop = FALSE], na.rm = TRUE)
  cell_metadata <- cbind.data.frame(data.frame(row.names = cell_union, 
                                               cells = cell_union), Pseudotime_assign = Pseudotime_assign, 
                                    srt@meta.data[cell_union, lineages, drop = FALSE])
  if (cell_density != 1) {
    bins <- cut(Pseudotime_assign, breaks = seq(min(Pseudotime_assign, 
                                                    na.rm = TRUE), max(Pseudotime_assign, na.rm = TRUE), 
                                                length.out = cell_bins), include.lowest = TRUE)
    ncells <- ceiling(max(table(bins), na.rm = TRUE) * cell_density)
    message("ncell/bin=", ncells, "(", cell_bins, "bins)")
    cell_keep <- unlist(sapply(levels(bins), function(x) {
      cells <- names(Pseudotime_assign)[bins == x]
      out <- sample(cells, size = min(length(cells), ncells))
      return(out)
    }))
    cell_metadata <- cell_metadata[cell_keep, , drop = FALSE]
  }
  cell_order_list <- list()
  for (l in lineages) {
    if (!is.null(reverse_ht) && (l %in% lineages[reverse_ht] || 
                                 l %in% reverse_ht)) {
      cell_metadata_sub <- na.omit(cell_metadata[, l, 
                                                 drop = FALSE])
      cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], 
                                                   decreasing = TRUE), , drop = FALSE]
    }
    else {
      cell_metadata_sub <- na.omit(cell_metadata[, l, 
                                                 drop = FALSE])
      cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], 
                                                   decreasing = FALSE), , drop = FALSE]
    }
    cell_order_list[[l]] <- paste0(rownames(cell_metadata_sub), 
                                   l)
  }
  if (!is.null(cell_annotation)) {
    cell_metadata <- cbind.data.frame(cell_metadata, cbind.data.frame(srt@meta.data[rownames(cell_metadata), 
                                                                                    c(intersect(cell_annotation, colnames(srt@meta.data))), 
                                                                                    drop = FALSE], t(srt@assays[[assay]]@data[intersect(cell_annotation, 
                                                                                                                                        rownames(srt@assays[[assay]])) %||% integer(), rownames(cell_metadata), 
                                                                                                                              drop = FALSE])))
  }
  dynamic <- list()
  if (is.null(features)) {
    for (l in lineages) {
      DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", 
                                           l)]][["DynamicFeatures"]]
      if (is.null(DynamicFeatures)) {
        stop("DynamicFeatures result for ", l, " found in the srt object. Should perform RunDynamicFeatures first!")
      }
      DynamicFeatures <- DynamicFeatures[DynamicFeatures$exp_ncells > 
                                           min_expcells & DynamicFeatures$r.sq > r.sq & 
                                           DynamicFeatures$dev.expl > dev.expl & DynamicFeatures$padjust < 
                                           padjust, , drop = FALSE]
      dynamic[[l]] <- DynamicFeatures
      features <- c(features, DynamicFeatures[["features"]])
    }
    message(length(unique(features)), " features from ", 
            paste0(lineages, collapse = ","), " passed the threshold (exp_ncells>", 
            min_expcells, " & r.sq>", r.sq, " & dev.expl>", 
            dev.expl, " & padjust<", padjust, "): \n", paste0(head(features, 
                                                                   10), collapse = ","), "...")
  }
  else {
    for (l in lineages) {
      DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", 
                                           l)]][["DynamicFeatures"]]
      if (is.null(DynamicFeatures)) {
        srt <- RunDynamicFeatures(srt, lineages = l, 
                                  features = features, assay = assay, slot = slot, 
                                  family = family, libsize = libsize)
        DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", 
                                             l)]][["DynamicFeatures"]]
      }
      if (any(!features %in% rownames(DynamicFeatures))) {
        srt <- RunDynamicFeatures(srt, lineages = l, 
                                  features = features[!features %in% rownames(DynamicFeatures)], 
                                  assay = assay, slot = slot, family = family, 
                                  libsize = libsize)
        DynamicFeatures <- rbind(DynamicFeatures, srt@tools[[paste0("DynamicFeatures_", 
                                                                    l)]][["DynamicFeatures"]])
      }
      DynamicFeatures <- DynamicFeatures[features, , drop = FALSE]
      dynamic[[l]] <- DynamicFeatures
    }
  }
  all_calculated <- Reduce(union, lapply(lineages, function(l) rownames(srt@tools[[paste0("DynamicFeatures_", 
                                                                                              l)]][["DynamicFeatures"]])))
  features_tab <- table(features)
  features <- names(features_tab)[which(features_tab %in% 
                                          (num_intersections %||% seq_along(lineages)))]
  if (!all(features %in% all_calculated)) {
    message("Some features were missing in at least one lineage: \n", 
            paste0(head(features[!features %in% all_calculated], 
                        10), collapse = ","), "...")
    features <- union(features, all_calculated)
  }
  if(!is.null(study_gene_set)){
    features=intersect(features,study_gene_set)
    message("Detect study_gene_set, features now length:",length(features),"(",head(features),"...)")
  }
  gene <- features[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  if (length(gene) == 0 && length(meta) == 0) {
    stop("No dynamic features found in the meta.data or in the assay: ", 
         assay)
  }
  feature_metadata <- data.frame(row.names = features, features = features)
  for (l in lineages) {
    feature_metadata[rownames(dynamic[[l]]), paste0(l, order_by)] <- dynamic[[l]][, 
                                                                                  order_by]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "exp_ncells")] <- dynamic[[l]][, 
                                                                                      "exp_ncells"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "r.sq")] <- dynamic[[l]][, 
                                                                                "r.sq"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "dev.expl")] <- dynamic[[l]][, 
                                                                                    "dev.expl"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "padjust")] <- dynamic[[l]][, 
                                                                                   "padjust"]
  }
  feature_metadata[, order_by] <- apply(feature_metadata[, 
                                                         paste0(lineages, order_by), drop = FALSE], 1, max, na.rm = TRUE)
  feature_metadata <- feature_metadata[order(feature_metadata[, 
                                                              order_by], decreasing = decreasing), , drop = FALSE]
  feature_metadata <- feature_metadata[rownames(feature_metadata) %in% 
                                         features, , drop = FALSE]
  features <- rownames(feature_metadata)
  if (!is.null(feature_annotation)) {
    feature_metadata <- cbind.data.frame(feature_metadata, 
                                         srt@assays[[assay]]@meta.features[rownames(feature_metadata), 
                                                                           feature_annotation, drop = FALSE])
  }
  if (isTRUE(use_fitted)) {
    mat_list <- list()
    for (l in lineages) {
      fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", 
                                         l)]][["fitted_matrix"]][, -1]
      rownames(fitted_matrix) <- paste0(rownames(fitted_matrix), 
                                        l)
      mat_list[[l]] <- t(fitted_matrix[, features])
    }
    mat_raw <- do.call(cbind, mat_list)
  }
  else {
    mat_list <- list()
    Y_libsize <- colSums(slot(srt@assays[[assay]], "counts"))
    for (l in lineages) {
      cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
      mat_tmp <- as_matrix(rbind(slot(srt@assays[[assay]], 
                                      slot)[gene, cells, drop = FALSE], t(srt@meta.data[cells, 
                                                                                        meta, drop = FALSE])))[features, , drop = FALSE]
      if (isTRUE(lib_normalize) && min(mat_tmp, na.rm = TRUE) >= 
          0) {
        if (!is.null(libsize)) {
          libsize_use <- libsize
        }
        else {
          libsize_use <- Y_libsize[colnames(mat_tmp)]
          isfloat <- any(libsize_use%%1 != 0, na.rm = TRUE)
          if (isTRUE(isfloat)) {
            libsize_use <- rep(1, length(libsize_use))
            warning("The values in the 'counts' slot are non-integer. Set the library size to 1.", 
                    immediate. = TRUE)
          }
        }
        mat_tmp[gene, ] <- t(t(mat_tmp[gene, , drop = FALSE])/libsize_use * 
                               median(Y_libsize))
      }
      colnames(mat_tmp) <- paste0(colnames(mat_tmp), l)
      mat_list[[l]] <- mat_tmp
    }
    mat_raw <- do.call(cbind, mat_list)
  }
  mat <- matrix_process(mat_raw, method = exp_method)
  mat[is.infinite(mat)] <- max(abs(mat[!is.infinite(mat)]), 
                               na.rm = TRUE) * ifelse(mat[is.infinite(mat)] > 0, 1, 
                                                      -1)
  mat[is.na(mat)] <- mean(mat, na.rm = TRUE)
  mat_split <- mat[, unlist(cell_order_list[feature_split_by]), 
                   drop = FALSE]
  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", 
                                                      "log2fc")) {
      b <- ceiling(min(abs(quantile(mat, c(0.01, 0.99), 
                                    na.rm = TRUE)), na.rm = TRUE) * 2)/2
      colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette, 
                                                                 palcolor = heatmap_palcolor))
    }
    else {
      b <- quantile(mat, c(0.01, 0.99), na.rm = TRUE)
      colors <- colorRamp2(seq(b[1], b[2], length = 100), 
                           palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    }
  }
  else {
    colors <- colorRamp2(seq(limits[1], limits[2], length = 100), 
                         palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
  }
  lgd <- list()
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, 
                        border = TRUE)
  ha_top_list <- list()
  pseudotime <- na.omit(unlist(cell_metadata[, lineages]))
  pseudotime_col <- colorRamp2(breaks = seq(min(pseudotime, 
                                                na.rm = TRUE), max(pseudotime, na.rm = TRUE), length = 100), 
                               colors = palette_scp(palette = pseudotime_palette, palcolor = pseudotime_palcolor))
  for (l in lineages) {
    ha_top_list[[l]] <- HeatmapAnnotation(Pseudotime = anno_simple(x = cell_metadata[gsub(pattern = l, 
                                                                                          replacement = "", x = cell_order_list[[l]]), l], 
                                                                   col = pseudotime_col, which = ifelse(flip, "row", 
                                                                                                        "column"), na_col = "transparent", border = TRUE), 
                                          which = ifelse(flip, "row", "column"), show_annotation_name = l == 
                                            lineages[1], annotation_name_side = ifelse(flip, 
                                                                                       "top", "left"))
  }
  lgd[["pseudotime"]] <- Legend(title = "Pseudotime", col_fun = pseudotime_col, 
                                border = TRUE)
  if (!is.null(cell_annotation)) {
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_annotation_palette[i]
      palcolor <- cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, cellan]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, 
                                                    FALSE))
        }
        else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", 
                                x = cell_order_list[[l]])
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(x = as.character(cell_anno[lineage_cells]), 
                                           col = palette_scp(cell_anno, palette = palette, 
                                                             palcolor = palcolor), which = ifelse(flip, 
                                                                                                  "row", "column"), na_col = "transparent", 
                                           border = TRUE)
          anno_args <- c(ha_cell, which = ifelse(flip, 
                                                 "row", "column"), show_annotation_name = l == 
                           lineages[1], annotation_name_side = ifelse(flip, 
                                                                      "top", "left"))
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), 
                                                                   names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          }
          else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], 
                                  ha_top)
          }
        }
        lgd[[cellan]] <- Legend(title = cellan, labels = levels(cell_anno), 
                                legend_gp = gpar(fill = palette_scp(cell_anno, 
                                                                    palette = palette, palcolor = palcolor)), 
                                border = TRUE)
      }
      else {
        col_fun <- colorRamp2(breaks = seq(min(cell_anno, 
                                               na.rm = TRUE), max(cell_anno, na.rm = TRUE), 
                                           length = 100), colors = palette_scp(palette = palette, 
                                                                               palcolor = palcolor))
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", 
                                x = cell_order_list[[l]])
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(x = cell_anno[lineage_cells], 
                                           col = col_fun, which = ifelse(flip, "row", 
                                                                         "column"), na_col = "transparent", border = TRUE)
          anno_args <- c(ha_cell, which = ifelse(flip, 
                                                 "row", "column"), show_annotation_name = l == 
                           lineages[1], annotation_name_side = ifelse(flip, 
                                                                      "top", "left"))
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), 
                                                                   names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          }
          else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], 
                                  ha_top)
          }
        }
        lgd[[cellan]] <- Legend(title = cellan, col_fun = col_fun, 
                                border = TRUE)
      }
    }
  }
  if (!is.null(separate_annotation)) {
    subplots_list <- list()
    for (i in seq_along(separate_annotation)) {
      cellan <- separate_annotation[[i]]
      palette <- separate_annotation_palette[i]
      palcolor <- separate_annotation_palcolor[[i]]
      if (length(cellan) == 1 && cellan %in% colnames(srt@meta.data)) {
        cell_anno <- srt@meta.data[[cellan]]
      }
      else {
        cell_anno <- numeric()
      }
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, 
                                                    FALSE))
        }
        else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", 
                                x = cell_order_list[[l]])
          subplots <- CellDensityPlot(srt = srt, cells = lineage_cells, 
                                      group.by = cellan, features = l, decreasing = TRUE, 
                                      x_order = "rank", palette = palette, palcolor = palcolor, 
                                      flip = flip, reverse = l %in% lineages[reverse_ht] || 
                                        l %in% reverse_ht) + theme_void()
          subplots_list[[paste0(cellan, ":", l)]] <- subplots
          graphics <- list()
          nm <- paste0(cellan, ":", l)
          funbody <- paste0("\n            g <- as_grob(subplots_list[['", 
                            nm, "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));\n            g$name <- '", 
                            nm, "';\n            grid.draw(g);\n            grid.rect(gp = gpar(fill = 'transparent', col = 'black'));\n            ")
          funbody <- gsub(pattern = "\n", replacement = "", 
                          x = funbody)
          eval(parse(text = paste("block_graphics <- function(index, levels) {", 
                                  funbody, "}", sep = "")), envir = environment())
          ha_cell <- list()
          ha_cell[[paste0(cellan, "\n(separate)")]] <- anno_block(panel_fun = block_graphics, 
                                                                  which = ifelse(flip, "row", "column"), show_name = l == 
                                                                    lineages[1])
          anno_args <- c(ha_cell, which = ifelse(flip, 
                                                 "row", "column"), show_annotation_name = l == 
                           lineages[1], annotation_name_side = ifelse(flip, 
                                                                      "top", "left"))
          anno_args <- c(anno_args, separate_annotation_params[setdiff(names(separate_annotation_params), 
                                                                       names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          }
          else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], 
                                  ha_top)
          }
        }
        lgd[[paste0("separate:", cellan)]] <- Legend(title = paste0(cellan, 
                                                                    "\n(separate)"), labels = levels(cell_anno), 
                                                     legend_gp = gpar(fill = palette_scp(cell_anno, 
                                                                                         palette = palette, palcolor = palcolor)), 
                                                     border = TRUE)
      }
      else {
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", 
                                x = cell_order_list[[l]])
          subplots <- DynamicPlot(srt = srt, cells = lineage_cells, 
                                  lineages = l, group.by = NULL, features = cellan, 
                                  line_palette = palette, line_palcolor = palcolor, 
                                  add_rug = FALSE, legend.position = "none", 
                                  compare_features = TRUE, x_order = "rank", 
                                  flip = flip, reverse = l %in% lineages[reverse_ht] || 
                                    l %in% reverse_ht) + theme_void()
          subplots_list[[paste0(paste0(cellan, collapse = ","), 
                                ":", l)]] <- subplots
          graphics <- list()
          nm <- paste0(paste0(cellan, collapse = ","), 
                       ":", l)
          funbody <- paste0("\n            g <- as_grob(subplots_list[['", 
                            nm, "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));\n            g$name <- '", 
                            nm, "';\n            grid.draw(g);\n            grid.rect(gp = gpar(fill = 'transparent', col = 'black'));\n            ")
          funbody <- gsub(pattern = "\n", replacement = "", 
                          x = funbody)
          eval(parse(text = paste("block_graphics <- function(index, levels) {", 
                                  funbody, "}", sep = "")), envir = environment())
          ha_cell <- list()
          ha_cell[[paste0(paste0(cellan, collapse = ","), 
                          "\n(separate)")]] <- anno_block(panel_fun = block_graphics, 
                                                          which = ifelse(flip, "row", "column"), show_name = l == 
                                                            lineages[1])
          anno_args <- c(ha_cell, which = ifelse(flip, 
                                                 "row", "column"), show_annotation_name = l == 
                           lineages[1], annotation_name_side = ifelse(flip, 
                                                                      "top", "left"))
          anno_args <- c(anno_args, separate_annotation_params[setdiff(names(separate_annotation_params), 
                                                                       names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          }
          else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], 
                                  ha_top)
          }
        }
        lgd[[paste0("separate:", paste0(cellan, collapse = ","))]] <- Legend(title = "Features\n(separate)", 
                                                                             labels = cellan, legend_gp = gpar(fill = palette_scp(cellan, 
                                                                                                                                  palette = palette, palcolor = palcolor)), 
                                                                             border = TRUE)
      }
    }
  }
  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat) <= n_split)) {
      row_split_raw <- row_split <- feature_split <- NULL
    }
    else {
      if (n_split == 1) {
        row_split_raw <- row_split <- feature_split <- setNames(rep(1, 
                                                                    nrow(mat_split)), rownames(mat_split))
      }
      else {
        if (split_method == "mfuzz") {
          status <- tryCatch(check_R("e1071"), error = identity)
          if (inherits(status, "error")) {
            warning("The e1071 package was not found. Switch split_method to 'kmeans'", 
                    immediate. = TRUE)
            split_method <- "kmeans"
          }
          else {
            mat_split_tmp <- mat_split
            colnames(mat_split_tmp) <- make.unique(colnames(mat_split_tmp))
            mat_split_tmp <- standardise(mat_split_tmp)
            min_fuzzification <- mestimate(mat_split_tmp)
            if (is.null(fuzzification)) {
              fuzzification <- min_fuzzification + 0.1
            }
            else {
              if (fuzzification <= min_fuzzification) {
                warning("fuzzification value is samller than estimated:", 
                        round(min_fuzzification, 2), immediate. = TRUE)
              }
            }
            cl <- e1071::cmeans(mat_split_tmp, centers = n_split, 
                                method = "cmeans", m = fuzzification)
            if (length(cl$cluster) == 0) {
              stop("Clustering with mfuzz failed (fuzzification=", 
                   round(fuzzification, 2), "). Please set a larger fuzzification parameter manually.")
            }
            row_split <- feature_split <- cl$cluster
          }
        }
        if (split_method == "kmeans") {
          km <- kmeans(mat_split, centers = n_split, 
                       iter.max = 10000, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "kmeans-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), 
                                        order_by]
          names(feature_y) <- rownames(mat_split)
          km <- kmeans(feature_y, centers = n_split, 
                       iter.max = 10000, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "hclust") {
          hc <- hclust(as.dist(dist(mat_split)))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
        if (split_method == "hclust-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), 
                                        order_by]
          names(feature_y) <- rownames(mat_split)
          hc <- hclust(as.dist(dist(feature_y)))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      df <- data.frame(row_split = row_split, order_by = feature_metadata[names(row_split), 
                                                                          order_by])
      df_order <- aggregate(df, by = list(row_split), 
                            FUN = mean)
      df_order <- df_order[order(df_order[["order_by"]], 
                                 decreasing = decreasing), , drop = FALSE]
      if (!is.null(split_order)) {
        df_order <- df_order[split_order, , drop = FALSE]
      }
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", 
                                                         i)
        level <- paste0("C", i, "(", sum(row_split == 
                                           raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split_raw <- row_split <- factor(row_split, 
                                           levels = split_levels)
      feature_split <- factor(feature_split, levels = paste0("C", 
                                                             seq_len(nrow(df_order))))
    }
  }
  else {
    row_split_raw <- row_split <- feature_split <- feature_split[row.names(mat)]
  }
  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split
  }
  else {
    feature_metadata[["feature_split"]] <- NA
  }
  ha_left <- NULL
  if (!is.null(row_split)) {
    if (isTRUE(cluster_row_slices)) {
      if (!isTRUE(cluster_rows)) {
        dend <- cluster_within_group(t(mat_split), row_split_raw)
        cluster_rows <- dend
        row_split <- length(unique(row_split_raw))
      }
    }
    funbody <- paste0("\n      grid.rect(gp = gpar(fill = palette_scp(", 
                      paste0("c('", paste0(levels(row_split_raw), collapse = "','"), 
                             "')"), ",palette = '", feature_split_palette, 
                      "',palcolor=c(", paste0("'", paste0(unlist(feature_split_palcolor), 
                                                          collapse = "','"), "'"), "))[nm]))\n    ")
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("panel_fun <- function(index, nm) {", 
                            funbody, "}", sep = "")), envir = environment())
    ha_clusters <- HeatmapAnnotation(features_split = anno_block(align_to = split(seq_along(row_split_raw), 
                                                                                  row_split_raw), panel_fun = getFunction("panel_fun", 
                                                                                                                          where = environment()), width = unit(0.1, "in"), 
                                                                 height = unit(0.1, "in"), show_name = FALSE, which = ifelse(flip, 
                                                                                                                             "column", "row")), which = ifelse(flip, "column", 
                                                                                                                                                               "row"), border = TRUE)
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    }
    else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- Legend(title = "Cluster", labels = intersect(levels(row_split_raw), 
                                                                     row_split_raw), legend_gp = gpar(fill = palette_scp(intersect(levels(row_split_raw), 
                                                                                                                                   row_split_raw), type = "discrete", palette = feature_split_palette, 
                                                                                                                         palcolor = feature_split_palcolor, matched = TRUE)), 
                               border = TRUE)
  }
  if (isTRUE(cluster_rows) && !is.null(cluster_features_by)) {
    mat_cluster <- mat[, unlist(cell_order_list[cluster_features_by]), 
                       drop = FALSE]
    if (is.null(row_split)) {
      dend <- as.dendrogram(hclust(as.dist(dist(mat_cluster))))
      dend_ordered <- reorder(dend, wts = colMeans(mat_cluster), 
                              agglo.FUN = mean)
      cluster_rows <- dend_ordered
    }
    else {
      row_split <- length(unique(row_split_raw))
      dend <- cluster_within_group2(t(mat_cluster), row_split_raw)
      cluster_rows <- dend
    }
  }
  l <- lineages[1]
  ht_args <- list(matrix = mat[, cell_order_list[[l]], drop = FALSE], 
                  col = colors, row_split = row_split, cluster_rows = cluster_rows, 
                  cluster_row_slices = cluster_row_slices, column_split = column_split, 
                  cluster_columns = cluster_columns, cluster_column_slices = cluster_column_slices, 
                  use_raster = TRUE)
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), 
                                          names(ht_args))])
  ht_list <- do.call(Heatmap, args = ht_args)
  features_ordered <- rownames(mat)[unlist(suppressWarnings(row_order(ht_list)))]
  feature_metadata[["index"]] <- setNames(object = seq_along(features_ordered), 
                                          nm = features_ordered)[rownames(feature_metadata)]
  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(features) > nlabel) {
        index_from <- ceiling((length(features_ordered)/nlabel)/2)
        index_to <- length(features_ordered)
        index <- unique(round(seq(from = index_from, 
                                  to = index_to, length.out = nlabel)))
      }
      else {
        index <- seq_along(features_ordered)
      }
    }
    else {
      index <- NULL
    }
  }
  else {
    index <- which(features_ordered %in% features_label)
    drop <- setdiff(features_label, features_ordered)
    if (length(drop) > 0) {
      warning(paste0(paste0(drop, collapse = ","), "was not found in the features"), 
              immediate. = TRUE)
    }
  }
  if (length(index) > 0) {
    ha_mark <- HeatmapAnnotation(gene = anno_mark(at = which(rownames(feature_metadata) %in% 
                                                               features_ordered[index]), labels = feature_metadata[which(rownames(feature_metadata) %in% 
                                                                                                                           features_ordered[index]), "features"], side = ifelse(flip, 
                                                                                                                                                                                "top", "left"), labels_gp = gpar(fontsize = label_size, 
                                                                                                                                                                                                                 col = label_color), link_gp = gpar(fontsize = label_size, 
                                                                                                                                                                                                                                                    col = label_color), which = ifelse(flip, "column", 
                                                                                                                                                                                                                                                                                       "row")), which = ifelse(flip, "column", "row"), 
                                 show_annotation_name = FALSE)
    if (is.null(ha_left)) {
      ha_left <- ha_mark
    }
    else {
      ha_left <- c(ha_mark, ha_left)
    }
  }
  ha_right <- NULL
  if (length(lineages) > 1) {
    ha_list <- list()
    for (l in lineages) {
      ha_list[[l]] <- anno_simple(x = is.na(feature_metadata[, 
                                                             paste0(l, order_by)]) + 0, col = c(`0` = "#181830", 
                                                                                                `1` = "transparent"), width = unit(5, "mm"), 
                                  height = unit(5, "mm"), which = ifelse(flip, 
                                                                         "column", "row"))
    }
    ha_lineage <- do.call("HeatmapAnnotation", args = c(ha_list, 
                                                        which = ifelse(flip, "column", "row"), annotation_name_side = ifelse(flip, 
                                                                                                                             "left", "top"), border = TRUE))
    if (is.null(ha_right)) {
      ha_right <- ha_lineage
    }
    else {
      ha_right <- c(ha_right, ha_lineage)
    }
  }
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featan <- feature_annotation[i]
      palette <- feature_annotation_palette[i]
      palcolor <- feature_annotation_palcolor[[i]]
      featan_values <- feature_metadata[, featan]
      if (!is.numeric(featan_values)) {
        if (is.logical(featan_values)) {
          featan_values <- factor(featan_values, levels = c(TRUE, 
                                                            FALSE))
        }
        else if (!is.factor(featan_values)) {
          featan_values <- factor(featan_values, levels = unique(featan_values))
        }
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(x = as.character(featan_values), 
                                            col = palette_scp(featan_values, palette = palette, 
                                                              palcolor = palcolor), which = ifelse(flip, 
                                                                                                   "column", "row"), na_col = "transparent", 
                                            border = TRUE)
        anno_args <- c(ha_feature, which = ifelse(flip, 
                                                  "column", "row"), show_annotation_name = TRUE, 
                       annotation_name_side = ifelse(flip, "left", 
                                                     "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), 
                                                                    names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        }
        else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(title = featan, labels = levels(featan_values), 
                                legend_gp = gpar(fill = palette_scp(featan_values, 
                                                                    palette = palette, palcolor = palcolor)), 
                                border = TRUE)
      }
      else {
        col_fun <- colorRamp2(breaks = seq(min(featan_values, 
                                               na.rm = TRUE), max(featan_values, na.rm = TRUE), 
                                           length = 100), colors = palette_scp(palette = palette, 
                                                                               palcolor = palcolor))
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(x = featan_values, 
                                            col = col_fun, which = ifelse(flip, "column", 
                                                                          "row"), na_col = "transparent", border = TRUE)
        anno_args <- c(ha_feature, which = ifelse(flip, 
                                                  "column", "row"), show_annotation_name = TRUE, 
                       annotation_name_side = ifelse(flip, "left", 
                                                     "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), 
                                                                    names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        }
        else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(title = featan, col_fun = col_fun, 
                                border = TRUE)
      }
    }
  }
  enrichment <- heatmap_enrichment(geneID = feature_metadata[["features"]], 
                                   geneID_groups = feature_metadata[["feature_split"]], 
                                   feature_split_palette = feature_split_palette, feature_split_palcolor = feature_split_palcolor, 
                                   ha_right = ha_right, flip = flip, anno_terms = anno_terms, 
                                   anno_keys = anno_keys, anno_features = anno_features, 
                                   terms_width = terms_width, terms_fontsize = terms_fontsize, 
                                   keys_width = keys_width, keys_fontsize = keys_fontsize, 
                                   features_width = features_width, features_fontsize = features_fontsize, 
                                   IDtype = IDtype, species = species, db_update = db_update, 
                                   db_version = db_version, db_combine = db_combine, convert_species = convert_species, 
                                   Ensembl_version = Ensembl_version, mirror = mirror, 
                                   db = db, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, 
                                   minGSSize = minGSSize, maxGSSize = maxGSSize, GO_simplify = GO_simplify, 
                                   GO_simplify_cutoff = GO_simplify_cutoff, simplify_method = simplify_method, 
                                   simplify_similarityCutoff = simplify_similarityCutoff, 
                                   pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, 
                                   topTerm = topTerm, show_termid = show_termid, topWord = topWord, 
                                   words_excluded = words_excluded)
  res <- enrichment$res
  ha_right <- enrichment$ha_right
  ht_list <- NULL
  for (l in lineages) {
    if (l == lineages[1]) {
      left_annotation <- ha_left
    }
    else {
      left_annotation <- NULL
    }
    if (l == lineages[length(lineages)]) {
      right_annotation <- ha_right
    }
    else {
      right_annotation <- NULL
    }
    ht_args <- list(name = l, matrix = if (flip) t(mat[, 
                                                       cell_order_list[[l]], drop = FALSE]) else mat[, 
                                                                                                     cell_order_list[[l]], drop = FALSE], col = colors, 
                    row_title = row_title %||% if (flip) l else character(0), 
                    row_title_side = row_title_side, column_title = column_title %||% 
                      if (flip) character(0) else l, column_title_side = if (flip) "top" else column_title_side, 
                    row_title_rot = row_title_rot, column_title_rot = column_title_rot, 
                    row_split = if (flip) column_split else row_split, 
                    column_split = if (flip) row_split else column_split, 
                    cluster_rows = if (flip) cluster_columns else cluster_rows, 
                    cluster_columns = if (flip) cluster_rows else cluster_columns, 
                    cluster_row_slices = if (flip) cluster_column_slices else cluster_row_slices, 
                    cluster_column_slices = if (flip) cluster_row_slices else cluster_column_slices, 
                    show_row_names = show_row_names, show_column_names = show_column_names, 
                    row_names_side = row_names_side, column_names_side = column_names_side, 
                    row_names_rot = row_names_rot, column_names_rot = column_names_rot, 
                    top_annotation = if (flip) left_annotation else ha_top_list[[l]], 
                    left_annotation = if (flip) ha_top_list[[l]] else left_annotation, 
                    bottom_annotation = if (flip) right_annotation else NULL, 
                    right_annotation = if (flip) NULL else right_annotation, 
                    show_heatmap_legend = FALSE, border = border, use_raster = use_raster, 
                    raster_device = raster_device, raster_by_magick = raster_by_magick, 
                    width = if (is.numeric(width[l])) unit(width[l], 
                                                           units = units) else NULL, height = if (is.numeric(height[l])) unit(height[l], 
                                                                                                                              units = units) else NULL)
    if (any(names(ht_params) %in% names(ht_args))) {
      warning("ht_params: ", paste0(intersect(names(ht_params), 
                                              names(ht_args)), collapse = ","), " were duplicated and will not be used.", 
              immediate. = TRUE)
    }
    ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), 
                                            names(ht_args))])
    if (isTRUE(flip)) {
      if (is.null(ht_list)) {
        ht_list <- do.call(Heatmap, args = ht_args)
      }
      else {
        ht_list <- ht_list %v% do.call(Heatmap, args = ht_args)
      }
    }
    else {
      ht_list <- ht_list + do.call(Heatmap, args = ht_args)
    }
  }
  if ((!is.null(row_split) && length(index) > 0) || any(c(anno_terms, 
                                                          anno_keys, anno_features)) || !is.null(width) || !is.null(height)) {
    fix <- TRUE
    if (is.null(width) || is.null(height)) {
      message("\nThe size of the heatmap is fixed because certain elements are not scalable.\nThe width and height of the heatmap are determined by the size of the current viewport.\nIf you want to have more control over the size, you can manually set the parameters 'width' and 'height'.\n")
    }
  }
  else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(width = width, height = height, 
                                   units = units, ha_top_list = ha_top_list, ha_left = ha_left, 
                                   ha_right = ha_right, ht_list = ht_list, legend_list = lgd, 
                                   flip = flip)
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]
  if (isTRUE(fix)) {
    fixsize <- heatmap_fixsize(width = width, width_sum = width_sum, 
                               height = height, height_sum = height_sum, units = units, 
                               ht_list = ht_list, legend_list = lgd)
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]
    gTree <- grid.grabExpr({
      draw(ht_list, annotation_legend_list = lgd)
      for (enrich in db) {
        enrich_anno <- names(ha_right)[grep(paste0("_split_", 
                                                   enrich), names(ha_right))]
        if (length(enrich_anno) > 0) {
          for (enrich_anno_element in enrich_anno) {
            enrich_obj <- strsplit(enrich_anno_element, 
                                   "_split_")[[1]][1]
            decorate_annotation(enrich_anno_element, 
                                slice = 1, {
                                  grid.text(paste0(enrich, " (", enrich_obj, 
                                                   ")"), x = unit(1, "npc"), y = unit(1, 
                                                                                      "npc") + unit(2.5, "mm"), just = c("left", 
                                                                                                                         "bottom"))
                                })
          }
        }
      }
      if (is.numeric(pseudotime_label)) {
        for (n in seq_along(pseudotime_label)) {
          pse <- pseudotime_label[n]
          col <- pseudotime_label_color[n]
          lty <- pseudotime_label_linetype[n]
          lwd <- pseudotime_label_linewidth[n]
          for (l in lineages) {
            for (slice in 1:max(nlevels(row_split), 
                                1)) {
              decorate_heatmap_body(l, {
                pseudotime <- cell_metadata[gsub(pattern = l, 
                                                 replacement = "", x = cell_order_list[[l]]), 
                                            l]
                i <- which.min(abs(pseudotime - pse))
                if (flip) {
                  x <- 1 - (i/length(pseudotime))
                  grid.lines(c(0, 1), c(x, x), gp = gpar(lty = lty, 
                                                         lwd = lwd, col = col))
                }
                else {
                  i <- which.min(abs(pseudotime - pse))
                  x <- i/length(pseudotime)
                  grid.lines(c(x, x), c(0, 1), gp = gpar(lty = lty, 
                                                         lwd = lwd, col = col))
                }
              }, row_slice = ifelse(flip, 1, slice), 
              column_slice = ifelse(flip, slice, 1))
            }
          }
        }
      }
    }, width = ht_width, height = ht_height, wrap = TRUE, 
    wrap.grobs = TRUE)
  }
  else {
    ht_width <- unit(width_sum, units = units)
    ht_height <- unit(height_sum, units = units)
    gTree <- grid.grabExpr({
      draw(ht_list, annotation_legend_list = lgd)
      for (enrich in db) {
        enrich_anno <- names(ha_right)[grep(paste0("_split_", 
                                                   enrich), names(ha_right))]
        if (length(enrich_anno) > 0) {
          for (enrich_anno_element in enrich_anno) {
            enrich_obj <- strsplit(enrich_anno_element, 
                                   "_split_")[[1]][1]
            decorate_annotation(enrich_anno_element, 
                                slice = 1, {
                                  grid.text(paste0(enrich, " (", enrich_obj, 
                                                   ")"), x = unit(1, "npc"), y = unit(1, 
                                                                                      "npc") + unit(2.5, "mm"), just = c("left", 
                                                                                                                         "bottom"))
                                })
          }
        }
      }
      if (is.numeric(pseudotime_label)) {
        for (n in seq_along(pseudotime_label)) {
          pse <- pseudotime_label[n]
          col <- pseudotime_label_color[n]
          lty <- pseudotime_label_linetype[n]
          lwd <- pseudotime_label_linewidth[n]
          for (l in lineages) {
            for (slice in 1:max(nlevels(row_split), 
                                1)) {
              decorate_heatmap_body(l, {
                pseudotime <- cell_metadata[gsub(pattern = l, 
                                                 replacement = "", x = cell_order_list[[l]]), 
                                            l]
                i <- which.min(abs(pseudotime - pse))
                if (flip) {
                  x <- 1 - (i/length(pseudotime))
                  grid.lines(c(0, 1), c(x, x), gp = gpar(lty = lty, 
                                                         lwd = lwd, col = col))
                }
                else {
                  i <- which.min(abs(pseudotime - pse))
                  x <- i/length(pseudotime)
                  grid.lines(c(x, x), c(0, 1), gp = gpar(lty = lty, 
                                                         lwd = lwd, col = col))
                }
              }, row_slice = ifelse(flip, 1, slice), 
              column_slice = ifelse(flip, slice, 1))
            }
          }
        }
      }
    }, width = ht_width, height = ht_height, wrap = TRUE, 
    wrap.grobs = TRUE)
  }
  if (isTRUE(fix)) {
    p <- panel_fix_overall(gTree, width = as.numeric(ht_width), 
                           height = as.numeric(ht_height), units = units)
  }
  else {
    p <- wrap_plots(gTree)
  }
  return(list(plot = p, matrix = mat, cell_order = cell_order_list,
              feature_split = feature_split, cell_metadata = cell_metadata, 
              feature_metadata = feature_metadata, enrichment = res))
}


source("/home/data/t020406/R/SCP-plot.R")




modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin =unit(c(0, 0, 0, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(0, 0, 0, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


Compare_Cluster_GO_results=
  function(cluster_list,
           pvalueCutoff = 0.1,
           qvalueCutoff = 0.1,
           pvalue_max_limit=1,
           ontology = "BP",
           path=getwd(),
           top_go_term=5,
           custom_design_go_term=1,
           custom_add_go_term=NULL,
           input_go_cluster_resluts_pathlist=1,
           return_plot_or_matrix="matrix",
           study_genes="Klf5",
           term_select_by_study_genes=F,
           convert_id_to_ENTREZID=T,
           simplify=T,#
           simplify_ver2=F,
           term_reoder=F
  ){
    
    top_go_term_row=top_go_term
    total_go_result_mat=data.frame(ID=1,
                                   Description="NA",
                                   GeneRatio=1,
                                   BgRatio=1,
                                   pvalue=1,
                                   p.adjust=1,
                                   qvalue=1,
                                   Count=1,
                                   group="NA")
    
    ignore=c()
    for(i in 1:length(cluster_list)){
      if(length(cluster_list[[i]])<5){#排除基因数小于5的类
        cat(names(cluster_list)[i],"gene number is <5! Ignore","\n")}
      ignore=c(ignore,i)
    }
    
    
    
    total_display_top_go_term=c()
    for (ls in 1:length(cluster_list)) {
      top_go_term=top_go_term_row
      name=paste(names(cluster_list)[ls]," ","(",length(cluster_list[[ls]]),")",sep = "")
      #编辑group名字
      
      if(length(input_go_cluster_resluts_pathlist)==1){
        
        if(convert_id_to_ENTREZID==T){
          cat("convert_id = TURE !","\n","Convert SYMBOL to ENTREZID !","\n")
          ENTREZID=bitr(cluster_list[[ls]],fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
          cluster_gene=ENTREZID$ENTREZID
        }else{
          cluster_gene=cluster_list[[ls]]
        }
        
        cat(paste("Doing",name,"GO...","\n"))
        # xx.mus.go=compareCluster(genes[[3]], OrgDb='org.Mm.eg.db', 
        #                    pvalueCutoff = 0.1, 
        #                    qvalueCutoff = 0.1, 
        #                    fun='enrichGO',
        #                    ont = "BP", 
        #                    readable=T)
        xx.mus.go=enrichGO(cluster_gene, OrgDb='org.Mm.eg.db', 
                           pvalueCutoff = pvalueCutoff, 
                           qvalueCutoff = qvalueCutoff, 
                           ont = ontology, readable=T)
        
        pvalueCutoff_change=pvalueCutoff
        qvalueCutoff_change=qvalueCutoff
        while(dim(xx.mus.go)[1]<top_go_term){
          pvalueCutoff_change=pvalueCutoff_change+0.1
          qvalueCutoff_change=qvalueCutoff_change+0.1
          xx.mus.go <- enrichGO(cluster_gene, OrgDb='org.Mm.eg.db', 
                                pvalueCutoff = pvalueCutoff_change,
                                qvalueCutoff = qvalueCutoff_change, 
                                ont = ontology, readable=T)
          cat(paste("pvalueCutoff and qvalueCutoff both adjust to:",pvalueCutoff_change,"\n"))
          if(pvalueCutoff_change==0.9){break}
        }
        
        
        
        #output=gsub(c("\\|"),c("+"),output)
        #output=gsub(c("\\="),c("_"),output)
        #output=gsub(c("\\("),c("_"),output)
        #output=gsub(c("\\)"),c("_"),output)
        
        
        if(simplify==T){
          cat("Simplify = TURE !","\n")
          xx.mus.go.f=clusterProfiler::simplify(xx.mus.go,by = "p.adjust",cutoff=0.7)
          output=paste(path,"/",name,"_",ontology,"_","Simplify_go.csv",sep="")
          write.csv(xx.mus.go.f,file = output)
          print(paste("\n","Export GO_Simplify_table file to :",output))
        }else{
          xx.mus.go.f=xx.mus.go
          output=paste(path,"/",name,"_",ontology,"_","Row_go.csv",sep="")
          write.csv(xx.mus.go.f,file = output)
          print(paste("\n","Export GO_row_table file to :",output))
        }
        
        
        if(simplify_ver2==T){
          cat("Simplify_ver2 = TURE !","\n")
          term_gene_list=list()
          for(n in 1:length(xx.mus.go.f[,1])){
            term_gene_list[[n]]=str_split(xx.mus.go.f[n,"geneID"],pattern = "/")[[1]]
          }
          
          if(length(xx.mus.go.f[,1])<=50){
            max_simplify_term=length(xx.mus.go.f[,1])
          }else{max_simplify_term=50}
          
          delete=c()
          for(n in max_simplify_term:2){
            for(i in 1:(n-1)){
              if(length(intersect(term_gene_list[[n]],term_gene_list[[i]]))>=length(term_gene_list[[n]])*0.8){delete=c(delete,n)}
            }
          }
          delete=rev(delete[!duplicated(delete)])
          cat("Remove GO term:","\n")
          print(head(xx.mus.go.f[delete,1:2]))
          if(length(delete)>6){cat("......",paste(length(delete)-6,"GO terms not show"),"\n")}
          xx.mus.go.f=xx.mus.go.f[-delete,]
          
          output=paste(path,"/",name,"_",ontology,"_","Simplify_ver2_go.csv",sep="")
          write.csv(xx.mus.go.f,file = output)
          print(paste("\n","Export GO_Simplify_ver2_table file to :",output))
          
        }
        
      }else{
        xx.mus.go.f=read.csv(input_go_cluster_resluts_pathlist[ls],row.names = 1)
        xx.mus.go.f$Cluster=NULL}
      
      if(term_select_by_study_genes==T){
        cat("Only display GO term have study_genes!","\n")
        xx.mus.go.f.s=xx.mus.go.f[1,]
        for(n in 1:length(xx.mus.go.f[,1])){
          term_gene=str_split(xx.mus.go.f[n,"geneID"],pattern = "/")[[1]]
          hitgene=term_gene[term_gene%in%study_genes]
          if(length(hitgene)>=1){
            xx.mus.go.f.s=rbind(xx.mus.go.f.s,xx.mus.go.f[n,])
          }
        }
        if(length(xx.mus.go.f.s[,1])<=1){
          paste("There are no match with study_genes!!! use oringial GO results")
        }else{
          xx.mus.go.f.s=xx.mus.go.f.s[-1,]
          xx.mus.go.f=xx.mus.go.f.s}
      }
      
      if(top_go_term>length(xx.mus.go.f[,1])){
        term_len=length(xx.mus.go.f[,1])
        cat(name," GO term is ",term_len," smaller than ",top_go_term," !!","\n")
        top_go_term=length(xx.mus.go.f[,1])
      }
      
      if(length(xx.mus.go.f[,1])>0){
        display_top_go_term=xx.mus.go.f[1:top_go_term,2]
        xx.mus.go.f=xx.mus.go.f[,c(1:7,9)]
        xx.mus.go.f[,9]=rep(name,length(xx.mus.go.f[,1]))
        colnames(xx.mus.go.f)=colnames(total_go_result_mat)
        total_go_result_mat=rbind(total_go_result_mat,xx.mus.go.f)
        total_display_top_go_term=c(total_display_top_go_term,display_top_go_term)
      }
      
      if(length(xx.mus.go.f[,1])==0){cat(name," GO term is ",term_len,"skip !!","\n")}
      
    }
    
    
    total_go_result_mat=total_go_result_mat[-1,]
    
    total_go_plot_mat=total_go_result_mat
    total_generatio=c()
    for (ratio in total_go_plot_mat$GeneRatio) {
      counts=as.numeric(strsplit(ratio,split="/",fixed=T)[[1]][1])
      total_counts=as.numeric(strsplit(ratio,split="/",fixed=T)[[1]][2])
      generatio=round(counts/total_counts,2)
      total_generatio=c(total_generatio,generatio)
    }
    total_go_plot_mat$GeneRatio=total_generatio
    
    
    if( !is.null(custom_add_go_term)){
      total_display_top_go_term=c(total_display_top_go_term,custom_add_go_term)
    }
    
    
    total_display_top_go_mat=subset(total_go_plot_mat,total_go_plot_mat$Description%in%total_display_top_go_term)
    
    order_mat=total_display_top_go_mat
    
    order=order_mat$group
    order=order[!duplicated(order)]
    
    order_mat$group=factor(order_mat$group,levels = order)
    order_mat=order_mat[order(order_mat$pvalue),]
    
    all_term=as.character(order_mat$Description)
    dup_term=order_mat$Description[duplicated(as.character(order_mat$Description))]
    
    term_in_cluster=1:length(all_term)
    
    dup_term_index=which(duplicated(as.character(order_mat$Description)))
    
    need_term=c()
    for(i in dup_term){
      index=which(all_term==i)
      need_term=c(need_term,index[which.min(term_in_cluster[index])])
    }
    
    need_rm_term_index=setdiff(dup_term_index,need_term)
    
    order_mat=order_mat[-dup_term_index,]
    
    order=c()
    for(i in levels(order_mat$group)){
      order=c(order,order_mat[order_mat$group==i,"Description"])
    }
    
    #site order
    # 
    # all_term=as.character(total_display_top_go_mat$Description)
    # dup_term=total_display_top_go_mat$Description[duplicated(as.character(total_display_top_go_mat$Description))]
    # 
    # dup_term_index=which(duplicated(as.character(total_display_top_go_mat$Description)))
    # 
    # group_term_num=as.numeric(table(total_display_top_go_mat$group))
    # 
    # term_in_cluster=c()
    # for(i in group_term_num){
    #   term_in_cluster=c(term_in_cluster,1:i)
    # }
    # 
    # need_term=c()
    # for(i in dup_term){
    #   index=which(all_term==i)
    #   need_term=c(need_term,index[which.min(term_in_cluster[index])])
    #   
    # }
    # need_rm_term_index=setdiff(dup_term_index,need_term)
    # 
    # 
    # term_order=total_display_top_go_mat[-dup_term_index,"Description"]
    # 
    
    # level=rev(total_go_result_mat$Description[!duplicated(as.character(total_go_result_mat$Description))])
    # total_go_result_mat$Description=factor(total_go_result_mat$Description,
    #                                        levels=level)
    
    plot_go_plot_mat=total_display_top_go_mat
    plot_go_plot_mat$Description=factor(plot_go_plot_mat$Description,
                                        levels = rev(order))
    
    
    if(length(custom_design_go_term)>1){
      plot_go_plot_mat=subset(total_go_plot_mat,total_go_plot_mat$Description%in%custom_design_go_term)
      plot_go_plot_mat$Description=factor(plot_go_plot_mat$Description,
                                          levels = c(custom_design_go_term))
    }
    
    plot_go_plot_mat$pvalue=as.numeric(plot_go_plot_mat$pvalue)
    
    ##这里大概根据聚类情况重新排序
    if(term_reoder==T){
      cl_mat=reshape2::dcast(plot_go_plot_mat[,c(2,3,5,6,7,8,9)],Description~group,value.var = "pvalue",fill = 1)
      rownames(cl_mat)=cl_mat$Description
      cl_mat$Description=NULL
      cl_mat=-log10(cl_mat)
      set.seed(2)
      k=kmeans(cl_mat,3)
      k_mat=as.data.frame(k$cluster)
      k_mat$gene=rownames(k_mat)
      k_mat=k_mat[order(k_mat$`k$cluster`),]
      
      term_order=k_mat$gene
      plot_go_plot_mat$Description=factor(plot_go_plot_mat$Description,
                                          levels = c(term_order))
    }
    
    
    
    ##plot function##
    colorbar = guide_colorbar(
      #title = 'Colorbar',
      title.position = 'top',  #left,right,top,bottom
      #title.theme = element_text(size = 15,face = "italic",colour = "red",angle = 45),
      #title.hjust =  .5,
      #title.vjust =  0,
      label = TRUE,
      label.position = 'top',
      #label.theme = element_text(size = 15,face = "italic",colour = "red",angle = 45),
      #label.hjust = .5,
      #label.vjust = .5,
      #barwidth = unit(6,"cm"),  #Default value is legend.key.width or legend.key.siz
      #barheight = unit(1.2,"cm"),  #Default value is legend.key.height or legend.key.size
      #nbin = 20,
      ticks = T,
      raster = F,
      draw.ulim = T,
      draw.llim = T,
      direction = "horizontal",  #"horizontal" or "vertical."
      reverse = T,
      order = 0)
    sizebar = guide_legend(
      #title = 'legend',
      title.position = "left",
      #title.hjust =  .5,
      #title.vjust = .5,
      #title.theme = element_text(size = 15,face = "italic",colour = "red",angle = 45),
      #keywidth = 2.5,
      #keyheight = 2,
      label.position = "right",
      direction = "vertical.",
      label = TRUE,
      #label.hjust = 0.5,
      #label.vjust = 0.5,
      #label.theme = element_text(size = 15,face = "italic",colour = "red",angle = 45),
      #nrow = 2,
      #byrow = TRUE,
      reverse = TRUE,
      order = 0)
    top_go_term=top_go_term_row
    scale.by="radius"
    scale.func <- switch(EXPR = scale.by, size = scale_size, #该函数可以改变圆圈大小，让圆圈之间的区别更大！！！
                         radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    
    order=plot_go_plot_mat$group
    order=order[!duplicated(order)]
    
    plot_go_plot_mat$group=factor(plot_go_plot_mat$group,levels = order)
    plot_go_plot_mat$pvalue=c(-log10(plot_go_plot_mat$pvalue))
    
    pmax=max(plot_go_plot_mat$pvalue)
    if(pmax<12){pmax=12}
    if(pvalue_max_limit!=1){
      pmax=pvalue_max_limit
      plot_go_plot_mat[plot_go_plot_mat$pvalue>pmax,"pvalue"]=pmax
    }
    
    red_ratio=(pmax-2)*5
    grey_ratio=10
    
    generatiomax=max(plot_go_plot_mat$GeneRatio)
    
    p=ggplot(plot_go_plot_mat,aes(x=group,y=Description,fill=pvalue))+
      #geom_point(aes(size=GeneRatio))+
      geom_point(aes_string(size="GeneRatio"),color="black",shape=21)+
      scale_size_continuous(breaks = seq(0,generatiomax,by=round(generatiomax/4,2)))+
      #scale.func(#range = c(1,6),
      #           breaks = seq(0,generatiomax,by=generatiomax/4))+  
      labs(fill="-Log10(Pvalue)",title = paste("GO top",top_go_term,ontology,"terms"),size="GeneRatio")+
      scale_fill_gradientn(colors = rev(c(colorRampPalette(colors =  rev(c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026")))(red_ratio),
                                          colorRampPalette(colors = c("#EEEEEE","#DEDEDE","#CECECE","#BEBEBE"))(grey_ratio))),
                           limits=c(0,pmax))+ #设置颜色梯度，2为3梯度
      theme_bw() +
      theme(#清空一切没必要的线
        legend.background=element_blank(), 
        legend.key=element_blank(),
        #legend.title = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust = 0.5),#标题居中
        panel.grid.major = element_line(color = "grey",size=0.3),
        panel.grid.minor = element_blank()
      ) +
      #scale_y_discrete(labels=function(x)str_wrap(x,width = 65))+#自动换行
      guides(color=colorbar,size=sizebar)
    #scale_x_discrete(labels=function(x)str_wrap(x,width = 3))#遇到空格换行
    
    print(p)
    
    if(return_plot_or_matrix=="matrix"){return(plot_go_plot_mat)}
    if(return_plot_or_matrix=="plot"){return(p)}
    if(return_plot_or_matrix=="total_matrix"){return(total_go_plot_mat)}
    
  }



Study_Gene_Distribution_in_gene_list=
  function(study_genes,gene_list,print=T){
    no_bind_index=length(gene_list)+1
    sub_gene_list=gene_list
    for(i in 1:length(gene_list)){
      part_gene=as.character(unlist(gene_list[c(i)]))
      sub_gene_list[[i]]=study_genes[study_genes%in%part_gene]
    }
    #names(sub_gene_list)=gene_list
    total=length(study_genes)
    outer=setdiff(study_genes,as.character(unlist(sub_gene_list)))
    if(print==T){
      print(paste(length(outer),"/",total,round((length(outer)/total)*100),"%"," Genes not in list"))}
    
    sub_gene_list[[no_bind_index]]=outer
    names(sub_gene_list)[no_bind_index]="No_binding"
    return(sub_gene_list)
  }


Plot_GO_result=
  function(xx.mus.go.f,
           display_top=10,
           name="GO-TERM",
           bar_color="#FFB6C1",
           xlab="Biological.Process",
           study_gene=c("Klf5"),
           sub_study_gene=c(),
           force_show_gene_in_term=F,
           force_show_gene_in_term_number=10,
           geom_text_size=6,
           save=F){
    go_plot=as.data.frame(xx.mus.go.f)
    go_plot=go_plot[c(1:display_top),]
    ##set plot info
    x=go_plot$Description
    z=c(-log10(go_plot$pvalue))
    
    
    display_gene=c()
    for(n in 1:length(go_plot[,1])){
      term_gene=str_split(go_plot[n,"geneID"],pattern="/")[[1]]
      bingo_gene=intersect(term_gene,study_gene)
      
      if(length(bingo_gene)<4 & length(bingo_gene)>=1){
        bingo_gene=c(bingo_gene,intersect(term_gene,sub_study_gene))}
      
      for(i in bingo_gene){
        if(i%in%sub_study_gene & length(bingo_gene)>7){
          bingo_gene=setdiff(bingo_gene,i)
        }
      }
      
      if(length(intersect(term_gene,study_gene))>=7){
        bingo_gene=intersect(term_gene,study_gene)
      }
      if(length(bingo_gene)<1){bingo_gene=" "}
      
      bingo_gene=paste(bingo_gene,collapse = ",")
      display_gene=c(display_gene,bingo_gene)
    }
    
    
    if(force_show_gene_in_term==T){
      num=force_show_gene_in_term_number
      
      for(n in 1:length(display_gene)){
        term_gene=str_split(go_plot[n,"geneID"],pattern="/")[[1]]
        term_num=paste("(n=",length(term_gene),")",sep="")
        bingo_gene=str_split(display_gene[n],pattern=",")[[1]]
        
        if(length(bingo_gene)<=num & length(term_gene)>num){
          without_gene=setdiff(term_gene,bingo_gene)
          bingo_gene=c(bingo_gene,without_gene[1:(num-length(bingo_gene))])
        }
        bingo_gene=setdiff(bingo_gene," ")
        if(length(term_gene)<=num){
          display_gene[n]=paste(term_gene,collapse = ",")
          display_gene[n]=paste(display_gene[n],term_num,collapse = ",")
        }else{
          display_gene[n]=paste(bingo_gene,collapse = ",")
          display_gene[n]=paste(display_gene[n],term_num,collapse = ",")
        }
      }
      geom_text_size=geom_text_size
    }
    
    #brewer.pal(n=2,palette = "Set1")
    p=ggplot(go_plot, aes(x=reorder(x,z),y=z)) +   
      geom_bar(stat="identity",fill=bar_color,width = 0.8) + 
      geom_text(aes(label=display_gene,y=0),color="black",
                position = "identity",
                vjust="center",#??ֱ???????뷽ʽ
                hjust=0,#ˮƽ???????뷽ʽ,0??????,1?Ҷ???
                size=geom_text_size)+
      coord_flip() +
      geom_abline(intercept = c(-log10(0.05)),slope = 0,color="#666666",size=0.8, linetype="dotted")+
      #scale_x_reverse()+
      labs(x=xlab,y="-Log10(P-value)",title=name)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      theme_gray()+
      theme(axis.line = element_line(size=1.5), #??????????ϸ
            panel.background = element_rect(fill=NA,size=rel(20)),
            panel.grid.minor = element_line(colour = NA),
            plot.title = element_text(size=18,face ="bold",colour = "black"),
            axis.text = element_text(size=17,face ="bold",colour = "black"),
            axis.title = element_text(size=18,face ="bold",colour = "black"),
            axis.ticks = element_line(size=1.5)#?Ŀ̶ȴ?ϸ
      )+scale_x_discrete(labels=function(x)str_wrap(x,width = 50))#?Զ?????
    p
    if(save==T){
      ggsave(filename =paste(name,"_go.pdf",sep=""), width=12, height=5,device = "pdf")}
    return(p)
  }


project2MST=
  function (cds, Projection_Method){
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)
  cds <- findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  tip_leaves <- names(which(degree(dp_mst) == 1))
  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  }
  else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(V(dp_mst)[suppressWarnings(.nei(closest_vertex_names[i], 
                                                         mode = "all"))])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], 
                                            neighbor)])
        }
        else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], 
                                              neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, dist(rbind(Z_i, tmp)))
      }
      if (!is(projection, "matrix")) 
        projection <- as.matrix(projection)
      P[, i] <- projection[which(distance == min(distance))[1], 
      ]
    }
  }
  colnames(P) <- colnames(Z)
  dp <- as.matrix(dist(t(P)))
  min_dist = min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0
  cellPairwiseDistances(cds) <- dp
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
  }
orderCells2=
function (cds, root_state = NULL, num_paths = NULL, reverse = NULL) 
{
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 
            0))) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  root_cell <- select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  if (cds@dim_reduce_type == "ICA") {
    if (is.null(num_paths)) {
      num_paths = 1
    }
    adjusted_S <- t(cds@reducedDimS)
    dp <- as.matrix(dist(adjusted_S))
    cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    order_list <- extract_good_branched_ordering(res$subtree, 
                                                 res$root, cellPairwiseDistances(cds), num_paths, 
                                                 FALSE)
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name
    minSpanningTree(cds) <- as.undirected(order_list$cell_ordering_tree)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
    ]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)), 
    ]$cell_state
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                        2)]$name
    minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- as.undirected(order_list$cell_ordering_tree)
  }
  else if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
    ]$pseudo_time
    K_old <- reducedDimK(cds)
    old_dp <- cellPairwiseDistances(cds)
    old_mst <- minSpanningTree(cds)
    old_A <- reducedDimA(cds)
    old_W <- reducedDimW(cds)
    cds <- project2MST(cds, project_point_to_line_segment)
    minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(V(old_mst)$name == root_cell, 
                           arr.ind = T)
    cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == 
                                          root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 
                                1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% 
                                              tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- select_root_cell(cds, root_state, reverse)
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, 
                                                           root_cell)
    pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)), 
    ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      pData(cds)$State <- cc_ordering[closest_vertex[, 
                                                     1], ]$cell_state
    }
    reducedDimK(cds) <- K_old
    cellPairwiseDistances(cds) <- old_dp
    minSpanningTree(cds) <- old_mst
    reducedDimA(cds) <- old_A
    reducedDimW(cds) <- old_W
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                        2)]$name
  }
  else if (cds@dim_reduce_type == "SimplePPT") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
    ]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)), 
    ]$cell_state
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                        2)]$name
  }
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  cds
}