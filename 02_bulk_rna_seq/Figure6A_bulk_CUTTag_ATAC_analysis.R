source(("F:/Backup_20240703/秦组/Data processing/WXY_R_function/WXY_R_function.R"))
Load_library()
setwd("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation")


# '-------------------------------- ####
# CUTTag Part ####
# 1. load peak ####
setwd("E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3")

CUTTag_bed=list.files("F:/Backup_20240703/秦组/Data processing/#Klf5 Figure#/All_CUTTag_bedfile/q0_sorted/",full.names = T)

K2_L=readPeakFile(CUTTag_bed[9])
K4_L=readPeakFile(CUTTag_bed[11])
K5_L=readPeakFile(CUTTag_bed[13])


K2_ESC=readPeakFile(CUTTag_bed[5])
K4_ESC=readPeakFile(CUTTag_bed[6])
K5_ESC=readPeakFile(CUTTag_bed[4])



# 2. K245 overkap ####

Make_3_Sample_Overlaps_Clustered_Bedfile(peak1 = K2_L,
                                         peak2 = K4_L,
                                         peak3 = K5_L,
                                         name1 = "K2_L",name2 = "K4_L",name3 = "K5_L",maxgap = 400
                                         )
Make_3_Sample_Overlaps_Clustered_Bedfile(peak1 = K2_ESC,
                                         peak2 = K4_ESC,
                                         peak3 = K5_ESC,
                                         name1 = "K2_ESC",name2 = "K4_ESC",name3 = "K5_ESC",maxgap = 400
)

Split_Bed_by_Groups("title_K2_L_vs_K4_L_vs_K5_L_400.bed")
Split_Bed_by_Groups("title_K2_ESC_vs_K4_ESC_vs_K5_ESC_400.bed")

# 3. K245 split peak annotation ####

peakfile=list.files("title_K2_L_vs_K4_L_vs_K5_L_400/",full.names = T)

genelist=list()
for(p in 1:length(peakfile)){
  name=gsub("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_","",peakfile[p])
  name=gsub("_.bed","",name)
  peak=readPeakFile(peakfile[p])
  anno=Peak_Anno(peak,up = 20000,down = 20000)
  genelist[[p]]=anno
  names(genelist)[p]=name
}


Study_Gene_Distribution_in_gene_list(studygene,genelist)


K2_L_anno=Peak_Anno(K2_L,up = 20000,down = 20000,return = "promoter_matrix")
K4_L_anno=Peak_Anno(K4_L,up = 20000,down = 20000,return = "promoter_matrix")
K5_L_anno=Peak_Anno(K5_L,up = 20000,down = 20000,return = "promoter_matrix")



olpeak=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400.bed")

olpeak_anno=Peak_Anno(olpeak,up = 20000,down = 20000,return = "promoter_matrix")

TSCr=ETP1[grepl("TSC",ETP1$Cluster),1]
PGCr=ETP1[grepl("PGCLC",ETP1$Cluster),1]
ESCr=ETP1[grepl("ESC",ETP1$Cluster),1]

TSCr_mat=data.frame(TSCr=TSCr,row.names = TSCr)
PGCr_mat=data.frame(PGCr=PGCr,row.names = PGCr)
ESCr_mat=data.frame(ESCr=ESCr,row.names = ESCr)

olpeak_anno_mat=Btable_insert_into_Atable(olpeak_anno,TSCr_mat,index_A =17 )
olpeak_anno_mat=Btable_insert_into_Atable(olpeak_anno_mat,PGCr_mat,index_A =17 )
olpeak_anno_mat=Btable_insert_into_Atable(olpeak_anno_mat,ESCr_mat,index_A =17 )

anno_mat=olpeak_anno_mat[,19:21]
anno_mat[!is.na(anno_mat)]=1
anno_mat[is.na(anno_mat)]=0
anno_mat=apply(anno_mat,2,as.numeric)

anno=data.frame(row.names = rownames(anno_mat),kmean=olpeak_anno_mat$V4)

pheatmap::pheatmap(anno_mat,annotation_row = anno,
                   cluster_row = F,
                   cluster_cols = F,show_rownames = F)

table(olpeak_anno_mat$V4)
sub_olpeak_anno_mat=olpeak_anno_mat[olpeak_anno_mat$V4=="K4_L(n=6937)",]
sub_olpeak_anno_mat$position=paste0(sub_olpeak_anno_mat$seqnames,":",sub_olpeak_anno_mat$start,"-",sub_olpeak_anno_mat$end)
sub_olpeak_anno_mat[sub_olpeak_anno_mat$SYMBOL%in%studygene,]


K245_split_peak_distribustion=function(){}

common=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K2_L+K4_L+K5_L_n_20920_.bed")
k2=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K2_L_n_11657_.bed")
k4=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K4_L_n_6937_.bed")
k5=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K5_L_n_5810_.bed")
k2_k4=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K2_L+K4_L_n_7164_.bed")
k2_k5=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K2_L+K5_L_n_3208_.bed")
k4_k5=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K4_L+K5_L_n_2016_.bed")


filelist=list(#ESC_K5=ESC_K5,
  #TSC_K5=TSC_K5,
  common=common,
  k2=k2,
  k4=k4,
  k5=k5,
  k2_k4=k2_k4,
  k2_k5=k2_k5,
  k4_k5=k4_k5)

row_list=lapply(filelist, annotatePeak,TxDb=txdb,tssRegion=c(-2000, 2000),annoDb="org.Mm.eg.db",verbose=F)
list=row_list

plotAnnoBar(list)


K245_unique_peak_annotation=function(){}


k2_cuttag_u=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K2_L_n_11657_.bed")
k4_cuttag_u=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K4_L_n_6937_.bed")
k5_cuttag_u=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K5_L_n_5810_.bed")


usegene=PGCu
k2_cuttag_u_anno=intersect(Peak_Anno(k2_cuttag_u,up = 20000,down = 20000),usegene)
k4_cuttag_u_anno=intersect(Peak_Anno(k4_cuttag_u,up = 20000,down = 20000),usegene)
k5_cuttag_u_anno=intersect(Peak_Anno(k5_cuttag_u,up = 20000,down = 20000),usegene)

Venner_plot(list(k2_cuttag_u_anno=k2_cuttag_u_anno,
                 k4_cuttag_u_anno=k4_cuttag_u_anno,
                 k5_cuttag_u_anno=k5_cuttag_u_anno))

venn=get.venn.partitions(list(k2_cuttag_u_anno=k2_cuttag_u_anno,
                         k4_cuttag_u_anno=k4_cuttag_u_anno,
                         k5_cuttag_u_anno=k5_cuttag_u_anno))

names(venn$..values..)=venn$..count..
venn$..values..$`61`


# 4. K245 bind ETP gene ####

ETP1=read.csv("F:/Backup_20240703/秦组/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/ETP_gene_WGCNA_ver3_4431_.csv")
ETP1=ETP1[ETP1$Cluster!="Primed",]
ETP1$Cluster=factor(ETP1$Cluster,levels = c("ESC","TSC","PGCLC","ESC&TSC","ESC&PGCLC",
                                            "TSC&PGCLC","ESC&TSC&PGCLC"))
ETP1=ETP1[order(ETP1$Cluster),]

K2_L_anno=na.omit(Peak_Anno(K2_L,up = 20000,down = 20000,return = "promoter_matrix"))
K4_L_anno=na.omit(Peak_Anno(K4_L,up = 20000,down = 20000,return = "promoter_matrix"))
K5_L_anno=na.omit(Peak_Anno(K5_L,up = 20000,down = 20000,return = "promoter_matrix"))



K2_L_mat=data.frame(K2_L=K2_L_anno,row.names = K2_L_anno)
K4_L_mat=data.frame(K4_L=K4_L_anno,row.names = K4_L_anno)
K5_L_mat=data.frame(K5_L=K5_L_anno,row.names = K5_L_anno)


ETP1_mat=Btable_insert_into_Atable(ETP1,K2_L_mat,index_A =1 )
ETP1_mat=Btable_insert_into_Atable(ETP1_mat,K4_L_mat,index_A =1 )
ETP1_mat=Btable_insert_into_Atable(ETP1_mat,K5_L_mat,index_A =1 )

anno_mat=ETP1_mat[,7:9]
anno_mat[!is.na(anno_mat)]=1
anno_mat[is.na(anno_mat)]=0
anno_mat=as.data.frame(apply(anno_mat,2,as.numeric))

anno=na.omit(data.frame(row.names = rownames(anno_mat),kmean=ETP1$Cluster))

plot_mat.class=anno$kmean
gaps.row=c()
gaps=0
for(i in 1:c(length(table(plot_mat.class))-1)){
  gaps=gaps+as.numeric(table(plot_mat.class))[i]
  gaps.row=c(gaps.row,gaps)
}

# 创建二值颜色映射
col_fun <- colorRampPalette(c("grey", "red"))(2)

# 绘制热图
Heatmap(
  anno_mat,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # 行间隔设置
  row_split = rep(1:(length(gaps.row) + 1), 
                  diff(c(0, gaps.row, nrow(anno_mat)))),
  col = col_fun,
  show_row_names = FALSE,
  use_raster = TRUE,
  raster_quality = 5,
  rect_gp = gpar(col = NA), # 无边框
  show_heatmap_legend = FALSE, # 通常注释矩阵不需要图例
)


length(intersect(K4_L_anno,ETP1$X))/dim(ETP1)[1]


VennPlot(list(K2=na.omit(ETP1_mat$K2_L),
                 K4=na.omit(ETP1_mat$K4_L),
                 K5=na.omit(ETP1_mat$K5_L)))
vennpart=get.venn.partitions(list(K2=na.omit(ETP1_mat$K2_L),
                         K4=na.omit(ETP1_mat$K4_L),
                         K5=na.omit(ETP1_mat$K5_L)))
vennpart$..count..
intersect(vennpart$..values..$`1`,studygene)

K245_unique_bind_ETP_gene=function(){}



return_mat=ETP1_mat[1,6:9]
for(i in 1:dim(ETP1_mat)[1]){
  if(sum(is.na(ETP1_mat[i,7:9]))==2 & !is.na(ETP1_mat[i,"K2_L"])){
    return_mat=rbind(return_mat,ETP1_mat[i,c(6:9)])
  }
}
dim(return_mat[-1,])



K245_bind_ETP_gene_peak_distribustion=function(){}


ETP1=read.csv("F:/Backup_20240703/秦组/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/ETP_gene_WGCNA_ver3_4431_.csv")

TSCr=ETP1[grepl("TSC",ETP1$Cluster),1]
PGCr=ETP1[grepl("PGCLC",ETP1$Cluster),1]
ESCr=ETP1[grepl("ESC",ETP1$Cluster),1]

TSCu=ETP1[ETP1$Cluster=="TSC",1]
PGCu=ETP1[ETP1$Cluster=="PGCLC",1]
ESCu=ETP1[ETP1$Cluster=="ESC",1]

use=ETP1[ETP1$Cluster=="ESC&PGCLC",1]

K2_L_anno=Peak_Anno(K2_L,up = 20000,down = 20000,return = "promoter_matrix")
K4_L_anno=Peak_Anno(K4_L,up = 20000,down = 20000,return = "promoter_matrix")
K5_L_anno=Peak_Anno(K5_L,up = 20000,down = 20000,return = "promoter_matrix")



esc_K2_L_anno=toGRanges(K2_L_anno[K2_L_anno$SYMBOL%in%use,c(1:3,23)])
esc_K4_L_anno=toGRanges(K4_L_anno[K4_L_anno$SYMBOL%in%use,c(1:3,23)])
esc_K5_L_anno=toGRanges(K5_L_anno[K5_L_anno$SYMBOL%in%use,c(1:3,23)])


filelist=list(#ESC_K5=ESC_K5,
  #TSC_K5=TSC_K5,
  tsc_K2=esc_K2_L_anno,
  tsc_K4=esc_K4_L_anno,
  tsc_K5=esc_K5_L_anno)

row_list=lapply(filelist, annotatePeak,TxDb=txdb,tssRegion=c(-2000, 2000),annoDb="org.Mm.eg.db",verbose=F)
list=row_list

plotAnnoBar(list)

for(i in 1:length(list)){
  sub_peak=list[[i]]
  sub_annomat=sub_peak@annoStat
  
  annoStat_new=data.frame(Feature=c("Promoter (1-2kb)","Exon","Intron"),
                          Frequency=c(sum(sub_annomat[2,2]),
                                      sum(sub_annomat[grep(pattern = "Exon",sub_annomat[,1],fixed = T),2]),
                                      sum(sub_annomat[grep(pattern = "Intron",sub_annomat[,1],fixed = T),2]))
  )
  annoStat_new=rbind(sub_annomat[grep(pattern = "Promoter (<=1kb)",sub_annomat[,1],fixed = T),],
                     annoStat_new,
                     sub_annomat[grep(pattern = "UTR",sub_annomat[,1],fixed = T),],
                     sub_annomat[grep(pattern = "Downstream (<=300)",sub_annomat[,1],fixed = T),],
                     sub_annomat[grep(pattern = "Distal Intergenic",sub_annomat[,1],fixed = T),]
  )
  annoStat_new$Feature=factor(annoStat_new$Feature,levels = annoStat_new$Feature)
  list[[i]]@annoStat=annoStat_new
}
pdf("Feature_Distribution.pdf",width =6,height = 3)
p=plotAnnoBar(list)
pp=p + guides(fill=guide_legend(reverse = T))+
  scale_fill_manual(values = colorRampPalette(c("#F08080","grey","#8DD3C7","#BEBADA","#2E8B57","#B3DE69","#80B1D3","#00BFFF"))(8))
pp
dev.off()


# 5. K245 bind TSC PGC 20 marker ####
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(pheatmap)

PGC_all=c("Dppa3", "Prdm1", "Tfap2c", "Nanos3","Id3", "Rhox6", "Irf1", "Sox15", "Rhox9", "Zic3", "Kit", 
          "Rbm47", "Phf19", "Msx1",  "Msx2","Dmrt1", "Dnd1", "Prdm14", "Dazl", "Ddx4")

TSC_all=c("Cdx2","Gata3","Tfap2c","Rhox5","Fgfr2","Gata2","Krt7","Krt18","Krt19","Ddah1","Rin3","Plet1","Hand1","Elf5",
          "Tead4","Eomes","Ets2","Yap1","Erf","Ly6a")


ESC_all=c("Klf4","Tbx3","Zfp42","Klf2","Tfcp2l1","Dppa5a","Esrrb","Prdm14","Nanog","Nr5a2",
          "Gbx2","Pbx1","Tcf15","Tcf3","Tgif1")


setwd("E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3")

CUTTag_bed=list.files("F:/Backup_20240703/秦组/Data processing/#Klf5 Figure#/All_CUTTag_bedfile/q0_sorted/",full.names = T)

K2_L=readPeakFile(CUTTag_bed[9])
K4_L=readPeakFile(CUTTag_bed[11])
K5_L=readPeakFile(CUTTag_bed[13])
K2_L_anno=Peak_Anno(K2_L,up = 20000,down = 20000)
K4_L_anno=Peak_Anno(K4_L,up = 20000,down = 20000)
K5_L_anno=Peak_Anno(K5_L,up = 20000,down = 20000)

k2_cuttag_u=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K2_L_n_11657_.bed")
k4_cuttag_u=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K4_L_n_6937_.bed")
k5_cuttag_u=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K5_L_n_5810_.bed")

k2_cuttag_u_anno=Peak_Anno(k2_cuttag_u,up = 20000,down = 20000)
k4_cuttag_u_anno=Peak_Anno(k4_cuttag_u,up = 20000,down = 20000)
k5_cuttag_u_anno=Peak_Anno(k5_cuttag_u,up = 20000,down = 20000)

study_gene=PGC_all

using_mat=data.frame(KLF5=rep(0,length(study_gene)),
                     KLF2=rep(0,length(study_gene)),
                     KLF4=rep(0,length(study_gene)),
                     row.names = study_gene)

using_mat_bind=using_mat
using_mat_bind[rownames(using_mat_bind)%in%K5_L_anno,1]=1
using_mat_bind[rownames(using_mat_bind)%in%K2_L_anno,2]=1
using_mat_bind[rownames(using_mat_bind)%in%K4_L_anno,3]=1

pheatmap(using_mat_bind,cluster_rows = F,cluster_cols = F,border_color = "black",
         cellwidth = 20,cellheight = 20,#gaps_row = c(3,7,13,17),
         color = c("#D4D4D4","#F54040"))


using_mat=data.frame(KLF5=rep(0,length(study_gene)),
                     KLF2=rep(0,length(study_gene)),
                     KLF4=rep(0,length(study_gene)),
                     row.names = study_gene)

using_mat_bind=using_mat
using_mat_bind[rownames(using_mat_bind)%in%k5_cuttag_u_anno,1]=1
using_mat_bind[rownames(using_mat_bind)%in%k2_cuttag_u_anno,2]=1
using_mat_bind[rownames(using_mat_bind)%in%k4_cuttag_u_anno,3]=1

pheatmap(using_mat_bind,cluster_rows = F,cluster_cols = F,border_color = "black",
         cellwidth = 20,cellheight = 20,#gaps_row = c(3,7,13,17),
         color = c("#D4D4D4","#F54040"))




study_gene=PGC_all


tpm_mat_mean=read.csv("F:/Backup_20240703/秦组/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_mat_mean)


using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%study_gene,c(16,12,4,8)]
colnames(using_mat)=c("iEmpty","iKlf5","iKlf2","iKlf4")


column_split_mat=data.frame(sample=colnames(using_mat),
                            split=rep("PNT_D2",4))

tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))
tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,1]+0.1))

Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,remove_na = F,by_step = 2,cut_max_min = "Yes",
                            study_gene = study_gene,bubble_scale_by_new_function = F,limit_exp = c(0,9),
                            limit_fold = c(-2,2),column_split_mat = column_split_mat,
                            title = "K25_PNT_Day25",z_score =F)#10


using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%study_gene,c(#1,116,
  19,20,22,23,106)]
using_mat=data.frame(TSC=apply(using_mat[,3:4],1,mean),
                     PGCLC=using_mat[,5],
                     ESC=apply(using_mat[,1:2],1,mean))


tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))

column_split_mat=data.frame(sample=colnames(using_mat),
                            split=rep("Control",3))

Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,by_step = 0.75,
                            study_gene =  study_gene,limit_exp = c(2,12),
                            limit_fold = c(-1.5,1.5),column_split_mat = column_split_mat,
                            col_vector =c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF","#FFFFFF","#D7B5D8","#DF65B0","#DD1C77","#980043"),
                            title = "Control",z_score =T)#10





# 6. GREAT GO####

study_term=c("vesicle-mediated transport","cellular response to DNA damage stimulus","cytoskeleton organization",
             "chromatin organization","regulation of stem cell population maintenance","meiotic cell cycle process",
             "immune system development","ossification","regulation of cell adhesion","regulation of apoptotic signaling pathway",
             "regulation of cell cycle","pattern specification process","cell fate commitment",
             "embryonic placenta morphogenesis","placenta blood vessel development")


plot_mat=Plot_GREAT_GO_results(path = "E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3/GREAT GO CUTTag/",
                               study_term =study_term,
                               pmax = 30,cluster_rows = F,
                               data_type ="greatExportAll",top_go_term = 0,
                               #col_reoder = c(3,1,2),
                               return_plot_or_matrix = "total_matrix")
plot_mat=plot_mat[plot_mat$Desc%in%study_term,]

brewer.pal(5,"Blues")
display_mat=plot_mat
display_mat[sapply(display_mat,function(x){x<2})]<-0
display_mat[sapply(display_mat,function(x){x>=2})]<-"*"
display_mat[sapply(display_mat,function(x){x==0})]<-""

bk=seq(0,30,0.1)
pheatmap::pheatmap(plot_mat,
                   cluster_row = F,
                   cluster_cols = F,
                   border=T,
                   breaks = bk,
                   display_numbers =as.matrix(display_mat),
                   border_color = "black",
                   show_rownames = T,
                   show_colnames = T,
                   cellwidth = 10,
                   cellheight = 10,
                   fontsize = 10,
                   #gaps_row = c(6,11),
                   legend_labels = "-log10 P-value",
                   color = c(colorRampPalette(colors = c("#FFFFFF"))((length(bk))*0.05),
                             colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.95))
)





# 7. Motif dicovery####
setwd("E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3/motifdiscover_ver3/")
show_motif=c("KLF5(Zf)","KLF6(Zf)","KLF3(Zf)","KLF1(Zf)","Sox2(HMG)","TEAD4(TEA)",
             "OCT4-SOX2-TCF-NANOG(POU,Homeobox,HMG)","Oct4(POU,Homeobox)","Esrrb(NR)","LEF1(HMG)","Nr5a2(NR)",
             "Jun-AP1(bZIP)","Fos(bZIP)","Eomes(T-box)","ELF5(ETS)","AP-2gamma(AP2)","Cdx2(Homeobox)"
             #"AP-1(bZIP)","JunB(bZIP)","Atf3(bZIP)", 
)
total_mat=Compare_Homer_known_motif_results(select_sites = c(".bed_400"),
                                            select_condition = "",max_value = 200,
                                            #return="total_matrix",
                                            rename_vector = c("KLF2","KLF2&KLF4",
                                                              "KLF2&KLF4&KLF5","KLF2&KLF5",
                                                              "KLF4","KLF4&KLF5","KLF5"),
                                            show_top_sig._motif = 20)


plot_mat=total_mat[,c(1:3,5,4,6,7)]

plot_mat$class=as.data.frame(strsplit2(rownames(plot_mat),split = "\\("))[,2]
plot_mat[plot_mat$class=="POU,Homeobox,HMG)","class"]="POU,Homeobox)"
plot_mat[plot_mat$class=="TEA)","class"]="TEAD)"
plot_mat[plot_mat$class=="bHLH),near_PU.1","class"]="bHLH)"

#?Ӹ?TF????
plot_mat$max=apply(plot_mat[,1:c(dim(plot_mat)[2]-1)],1,max)
delete_class=c("CCAAT)","Forkhead,bHLH)","NR)","?)","CP2)","CTF)","Homeobox)")#class member>1
plot_mat_new=t(data.frame(rep(1,dim(plot_mat)[2])))
colnames(plot_mat_new)=colnames(plot_mat)

for(c in levels(factor(plot_mat$class))){
  if(length(plot_mat[plot_mat$class==c,1])>10){
    sub_plot_mat=plot_mat[plot_mat$class==c,]
    sub_plot_mat=sub_plot_mat[order(sub_plot_mat$max,decreasing = T),]
    sub_plot_mat=sub_plot_mat[1:10,]
    plot_mat_new=rbind(plot_mat_new,sub_plot_mat)
  }else{
    sub_plot_mat=plot_mat[plot_mat$class==c,]
    sub_plot_mat=sub_plot_mat[order(sub_plot_mat$max,decreasing = T),]
    plot_mat_new=rbind(plot_mat_new,sub_plot_mat)}
}
plot_mat_new=plot_mat_new[-1,]
plot_mat_new=plot_mat_new[!plot_mat_new$class%in%delete_class,]


plot_mat=plot_mat_new
plot_mat$class=factor(plot_mat$class,levels = c("bHLH)","Zf)",
                                                "bZIP)","E2F)",
                                                "POU,Homeobox)","Forkhead)","IRF)"))
c("Homeobox)","bHLH)","HMG)","AP2)","HTH)",
  "T-box)","Forkhead)","POU,Homeobox)",
  "TEAD)","bZIP)","Zf)")

plot_mat=plot_mat[order(plot_mat$class),]
plot_mat.class=plot_mat$class

anno=data.frame(Motif_family=plot_mat$class,row.names =rownames(plot_mat) )

plot_mat$class=NULL
plot_mat$max=NULL


gaps.row=c()
gaps=0
for(i in 1:c(length(table(plot_mat.class))-1)){
  gaps=gaps+as.numeric(table(plot_mat.class))[i]
  gaps.row=c(gaps.row,gaps)
}
####

###??????
#plot_mat=plot_mat[,c(1,4,7,2,5,8,3,6,9)]
#colnames(plot_mat)
####

rownames(plot_mat)=as.data.frame(strsplit2(str_to_upper(rownames(plot_mat)),split = "\\("))[,1]
plot_mat=plot_mat[c(1:9,11,9,12:17),]

write.csv(plot_mat,"Heatmap_plot_mat2.csv")
#plot_mat=read.csv("Heatmap_plot_mat2.csv",header = T,row.names = 1)



display_mat=plot_mat
display_mat[sapply(display_mat,function(x){x<2})]<-0
display_mat[sapply(display_mat,function(x){x>=2})]<-"*"
display_mat[sapply(display_mat,function(x){x==0})]<-""

pheatmap::pheatmap(plot_mat,
                   cluster_row =T,
                   cluster_cols = F,
                   border=T,
                   border_color = "black",
                   show_rownames = T,
                   show_colnames = T,
                   clustering_method = "ward.D2",
                   cellwidth =10,
                   cellheight = 10,
                   fontsize = 10,
                   cutree_rows = 5,
                   #gaps_col = seq(2,2,by=2),
                   #gaps_row = gaps.row,
                   legend_labels = "-log10 P-value",
                   display_numbers =as.matrix(display_mat),
                   #main = title,
                   #filename = pdf,
                   color=c(colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))(200)))

# 8. overlap Klf5 unique Motif ####

ap2=readPeakFile("motifdiscover_ver3/AP2.sites.mm10.bed")
tead4=readPeakFile("motifdiscover_ver3/Tead4.sites.mm10.bed")

cuttagall=readPeakFile("E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3/title_K2_L_vs_K4_L_vs_K5_L_400.bed")
atacall=toGRanges(rpkm_mat[rownames(rpkm_mat)%in%rownames(cluster_results),2:4])

Make_3_Sample_Overlaps_Clustered_Bedfile(name1 = "AP2",
                                         name2 = "CUTTagall",
                                         name3 = "ATACall",
                                         peak1 = ap2,
                                         peak2 = cuttagall,
                                         peak3 = atacall,maxgap = 400)

Make_3_Sample_Overlaps_Clustered_Bedfile(name1 = "Tead4",
                                         name2 = "CUTTagall",
                                         name3 = "ATACall",
                                         peak1 = tead4,
                                         peak2 = cuttagall,
                                         peak3 = atacall,maxgap = 400)


ap2=readPeakFile("motifdiscover_ver3/title_AP2_vs_CUTTagall_vs_ATACall_400.bed")
tead4=readPeakFile("motifdiscover_ver3/title_Tead4_vs_CUTTagall_vs_ATACall_400.bed")

ap2_ues=ap2[!ap2$V4%in%c("AP2(n=544844)","CUTTagall|ATACall(n=5174)",
                         "ATACall(n=14312)",
                         "CUTTagall(n=17580)"),]
tead4_ues=tead4[!tead4$V4%in%c("Tead4(n=940506)",
                         "ATACall(n=13041)",
                         "CUTTagall(n=22866)"),]

Output_bedfile(as.data.frame(ap2_ues),name = "AP2_ol_CUTTAG_ATAC")
Output_bedfile(as.data.frame(tead4_ues),name = "Tead4_ol_CUTTAG_ATAC")






# '-------------------------------- ####
# ATAC Part ####
# 1. Load atac rpkm matrix####

rpkm_mat=read.csv("F:/Backup_20240703/秦组/Data processing/20221201_K245_bind_opening_expression/rpkm_mat.csv")
rpkm_mat=rpkm_mat[,-1]
rpkm_mat_rpkm=rpkm_mat[,c(6,8,10,12,14)]

rpkm_mat_rpkm_cv=rpkm_mat_rpkm[apply(rpkm_mat_rpkm,1,function(x) (sd(x)/(mean(x)))^2>0.2),]
dim(rpkm_mat_rpkm_cv)
head(rpkm_mat_rpkm_cv)

# 2. find unique open region####
library(Mfuzz)

rpkm_mat_rpkm_cv <- rpkm_mat_rpkm_cv[apply(rpkm_mat_rpkm_cv, 1, sd) > 0, ]
rpkm_mat_rpkm_cv=na.omit(rpkm_mat_rpkm_cv)

sample_info <- data.frame(
  row.names = colnames(rpkm_mat_rpkm_cv),
  sample =  colnames(rpkm_mat_rpkm_cv)# 根据你的实验设计修改
)

# 创建特征信息（featureData）
# 这是一个数据框，行数等于基因数，行名必须与表达矩阵的行名匹配。
# 通常可以先创建一个空的或包含基因ID的信息
feature_info <- data.frame(
  row.names = rownames(rpkm_mat_rpkm_cv),
  GeneID = rownames(rpkm_mat_rpkm_cv)
)

# 创建 ExpressionSet 对象
eset <- ExpressionSet(assayData = as.matrix(rpkm_mat_rpkm_cv),
                      phenoData = new("AnnotatedDataFrame", data = sample_info),
                      featureData = new("AnnotatedDataFrame", data = feature_info))

# 示例：过滤掉在所有样本中标准差为0的基因（即没有变化的基因）


# 方法1：使用 scale 函数手动标准化
exprs_matrix <- exprs(eset)
exprs_standardized <- t(scale(t(exprs_matrix))) # 对行进行标准化

# 创建新的 ExpressionSet
eset_s <- eset
exprs(eset_s) <- exprs_standardized


# # 设置一个可能的聚类数量范围，例如从 2 到 20
# c_range <- 2:20
# 
# # 计算不同 c 值下的最小距离
# # 此步骤可能较慢，取决于数据量和c_range的大小
# m <- mestimate(eset_s) # 先估算模糊化参数 m
# Dmin_vals <- c()
# for (i in c_range) {
#   cl <- mfuzz(eset_s, c = i, m = m)
#   Dmin_vals <- c(Dmin_vals, cl$Dmin)
# }

# 绘制 “肘部法则” 图
#plot(c_range, Dmin_vals, type = "b", pch = 19, xlab = "Number of clusters (c)", ylab = "Dmin")
# 选择图中“拐点”（elbow point）对应的 c 值

# 估算模糊加权指数 m (通常使用默认估算值即可，结果对m不敏感)
m <- mestimate(eset_s)
print(paste("The estimated m is:", m))

# 设置随机种子以保证结果可重现
set.seed(123)

# 执行 Mfuzz 聚类
cl <- mfuzz(eset_s, c = 8, m = m) # 将 c 替换为你选择的值


# 获取标准化后的表达矩阵
expression_data <- rpkm_mat_rpkm_cv

# 获取聚类结果
cluster_results <- data.frame(
  GeneID = rownames(expression_data),
  Cluster = cl$cluster,
  Max_Membership = apply(cl$membership, 1, max)
)

# 按聚类号和成员值排序
cluster_results <- cluster_results %>%
  arrange(Cluster, desc(Max_Membership))

# 重新排序表达矩阵以匹配聚类顺序
expression_ordered <- expression_data[cluster_results$GeneID, ]

# 创建行注释数据框
row_annotation <- data.frame(
  Cluster = factor(paste0("Cluster_", cluster_results$Cluster)),
  Membership = cluster_results$Max_Membership
)
rownames(row_annotation) <- cluster_results$GeneID

# 为不同聚类设置颜色
cluster_colors <- brewer.pal(n = length(unique(cl$cluster)), name = "Set3")
names(cluster_colors) <- unique(paste0("Cluster_", unique(cl$cluster)))

# 成员值颜色梯度
mem_colors <- colorRampPalette(c("white", "darkred"))(100)

annotation_colors <- list(
  Cluster = cluster_colors,
  Membership = mem_colors
)


#d=Cluster_mat_by_Kmeans(rpkm_mat_rpkm_cv,kmeans = 8,set.seed = 1)#auto 8
#anno=data.frame(row.names = rownames(d),kmean=d$kmean)

#d$kmean=factor(d$kmean,levels = c(3:1))
#d=d[order(d$kmean),]

expression_ordered$cluster=strsplit2(row_annotation$Cluster,"_")[,2]
plotmat=expression_ordered[expression_ordered$cluster%in%c(4,2,7,1),]
plotmat$cluster=factor(plotmat$cluster,levels=c(4,2,7,1))
plotmat=plotmat[order(plotmat$cluster),]

row_annotation=data.frame(cluster=plotmat$cluster,row.names = rownames(plotmat))

bk=seq(2,6,by=0.1)
pheatmap::pheatmap(#d[,1:c(dim(d)[2]-1)],
                   plotmat[,1:5]*2,
                   cluster_cols = F,cluster_rows = F,border=F,
                   breaks = bk,show_rownames = F,
                   annotation_row = row_annotation,
                   cellwidth = 30,fontsize = 20,raster =T,
                   #display_numbers = round(as.matrix(using_tpm),digits = 1),
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0),
                              colorRampPalette(colors = viridis(8, option = "B"))((length(bk))*1)))


peaks=rpkm_mat[rownames(rpkm_mat)%in%rownames(d),c(2:4)]
peaks$index=rownames(peaks)
peaks=toGRanges(peaks)

peaks_anno=Peak_Anno(peaks,up = 20000,down = 20000)

intersect(studygene,peaks_anno)
ha_row <- NULL
if (!is.null(row_annotation)) {
  ha_row <- rowAnnotation(
    df = row_annotation,
    show_annotation_name = TRUE,
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 14)
    )
  )
}

col_fun <- colorRamp2(
  breaks = bk,space = "RGB",
  colors = c(colorRampPalette(colors = viridis::viridis(8, option = "B"))(length(bk) * 1)
  )
)
 Heatmap(
  plotmat[,1:5] * 2,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = NA),
  col = col_fun,
  show_row_names = FALSE,
  
  # 精确控制每个单元格的尺寸
  #width = ncol(plotmat[, 1:5]) * unit(30, "pt"),
  #height = nrow(plotmat) * unit(30, "pt"),
  left_annotation = ha_row,
  use_raster = TRUE,
  raster_quality = 6
)
 


d[rownames(d)%in%pluripotency,]


epi_upeak=rpkm_mat[rownames(rpkm_mat)%in%rownames(cluster_results[cluster_results$Cluster==3,]),2:4]
ctrl_upeak=rpkm_mat[rownames(rpkm_mat)%in%rownames(cluster_results[cluster_results$Cluster==4,]),2:4]
k2_upeak=rpkm_mat[rownames(rpkm_mat)%in%rownames(cluster_results[cluster_results$Cluster==2,]),2:4]
k4_upeak=rpkm_mat[rownames(rpkm_mat)%in%rownames(cluster_results[cluster_results$Cluster==7,]),2:4]
k5_upeak=rpkm_mat[rownames(rpkm_mat)%in%rownames(cluster_results[cluster_results$Cluster==1,]),2:4]

Output_bedfile(epi_upeak,name = "epi_upeak_ver2")
Output_bedfile(ctrl_upeak,name = "ctrl_upeak_ver2")
Output_bedfile(k2_upeak,name = "k2_upeak_ver2")
Output_bedfile(k4_upeak,name = "k4_upeak_ver2")
Output_bedfile(k5_upeak,name = "k5_upeak_ver2")

k5_upeak=readPeakFile("k5_upeak_ver2.bed")

peaks_anno=Peak_Anno(k5_upeak,up = 20000,down = 20000,return = "promoter_matrix")
peaks_anno[peaks_anno$SYMBOL%in%studygene,]


rpkm_mat$position=paste0(rpkm_mat$Chr,":",rpkm_mat$Start,"-",rpkm_mat$End)
peaks_anno$position=paste0(peaks_anno$seqnames,":",peaks_anno$start,"-",peaks_anno$end)
peaks_anno[peaks_anno$SYMBOL%in%TSCu,"SYMBOL"]

rpkm_mat[rpkm_mat$position=="chr1:13132095-13132466",]

# 3. GREAT GO####

study_term=c("regulation of developmental growth","neuron differentiation","heart development",
             "DNA hypermethylation","histone modification","covalent chromatin modification","cellular response to leukemia inhibitory factor","regulation of stem cell population maintenance",
             "regulation of cell adhesion","animal organ morphogenesis","epithelial cell differentiation",
             "positive regulation of apoptotic process","skin development",
             "pattern specification process","blood vessel morphogenesis",
             "placenta development","cell fate commitment","meiotic cell cycle process")



plot_mat=Plot_GREAT_GO_results(path = "E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3/GREAT GO/",
                               study_term =study_term,pmax = 25,cluster_rows = F,
                               data_type ="greatExportAll",top_go_term = 0,
                               #col_reoder = c(3,1,2),
                               return_plot_or_matrix = "heatmap_matrix")
display_mat=plot_mat
display_mat[sapply(display_mat,function(x){x<2})]<-0
display_mat[sapply(display_mat,function(x){x>=2})]<-"*"
display_mat[sapply(display_mat,function(x){x==0})]<-""

bk=seq(0,25,0.1)
pheatmap::pheatmap(plot_mat,
                   cluster_row = F,
                   cluster_cols = F,
                   border=T,
                   breaks = bk,
                   display_numbers =as.matrix(display_mat),
                   border_color = "black",
                   show_rownames = T,
                   show_colnames = T,
                   cellwidth = 10,
                   cellheight = 10,
                   fontsize = 10,
                   #gaps_row = c(6,11),
                   legend_labels = "-log10 P-value",
                   color = c(colorRampPalette(colors = c("#FFFFFF"))((length(bk))*0.),
                             colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*1))
)

# 4. Motif dicovery####
setwd("E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3/motifdiscover_ver2/")
show_motif=c("KLF5(Zf)","KLF6(Zf)","KLF3(Zf)","KLF1(Zf)","Sox2(HMG)","TEAD4(TEA)",
             "OCT4-SOX2-TCF-NANOG(POU,Homeobox,HMG)","Oct4(POU,Homeobox)","Esrrb(NR)","LEF1(HMG)","Nr5a2(NR)",
             "Jun-AP1(bZIP)","Fos(bZIP)","Eomes(T-box)","ELF5(ETS)","AP-2gamma(AP2)","Cdx2(Homeobox)"
             #"AP-1(bZIP)","JunB(bZIP)","Atf3(bZIP)", 
)
total_mat=Compare_Homer_known_motif_results(select_sites = c("200"),
                                            select_condition = "",max_value = 200,
                                            #return="total_matrix",
                                            #rename_vector = c("ESC_unique","Common","TSC_unique"),
                                            show_top_sig._motif = 6)


plot_mat=total_mat

plot_mat$class=as.data.frame(strsplit2(rownames(plot_mat),split = "\\("))[,2]
plot_mat[plot_mat$class=="POU,Homeobox,HMG)","class"]="POU,Homeobox)"
plot_mat[plot_mat$class=="TEA)","class"]="TEAD)"
plot_mat[plot_mat$class=="bHLH),near_PU.1","class"]="bHLH)"

#?Ӹ?TF????
plot_mat$max=apply(plot_mat[,1:c(dim(plot_mat)[2]-1)],1,max)
delete_class=c("CCAAT)","Forkhead,bHLH)","NR)","?)","CP2)","CTF)","Homeobox)")#class member>1
plot_mat_new=t(data.frame(rep(1,dim(plot_mat)[2])))
colnames(plot_mat_new)=colnames(plot_mat)

for(c in levels(factor(plot_mat$class))){
  if(length(plot_mat[plot_mat$class==c,1])>10){
    sub_plot_mat=plot_mat[plot_mat$class==c,]
    sub_plot_mat=sub_plot_mat[order(sub_plot_mat$max,decreasing = T),]
    sub_plot_mat=sub_plot_mat[1:10,]
    plot_mat_new=rbind(plot_mat_new,sub_plot_mat)
  }else{
    sub_plot_mat=plot_mat[plot_mat$class==c,]
    sub_plot_mat=sub_plot_mat[order(sub_plot_mat$max,decreasing = T),]
    plot_mat_new=rbind(plot_mat_new,sub_plot_mat)}
}
plot_mat_new=plot_mat_new[-1,]
plot_mat_new=plot_mat_new[!plot_mat_new$class%in%delete_class,]


plot_mat=plot_mat_new
plot_mat$class=factor(plot_mat$class,levels = c("bHLH)","Zf)",
                                                "bZIP)","E2F)",
                                                "POU,Homeobox)","Forkhead)","IRF)"))
c("Homeobox)","bHLH)","HMG)","AP2)","HTH)",
  "T-box)","Forkhead)","POU,Homeobox)",
  "TEAD)","bZIP)","Zf)")

plot_mat=plot_mat[order(plot_mat$class),]
plot_mat.class=plot_mat$class

anno=data.frame(Motif_family=plot_mat$class,row.names =rownames(plot_mat) )

plot_mat$class=NULL
plot_mat$max=NULL


gaps.row=c()
gaps=0
for(i in 1:c(length(table(plot_mat.class))-1)){
  gaps=gaps+as.numeric(table(plot_mat.class))[i]
  gaps.row=c(gaps.row,gaps)
}
####

###??????
#plot_mat=plot_mat[,c(1,4,7,2,5,8,3,6,9)]
#colnames(plot_mat)
####

rownames(plot_mat)=as.data.frame(strsplit2(str_to_upper(rownames(plot_mat)),split = "\\("))[,1]
plot_mat=plot_mat[c(1:9,11,9,12:17),]

write.csv(plot_mat,"Heatmap_plot_mat2.csv")
#plot_mat=read.csv("Heatmap_plot_mat2.csv",header = T,row.names = 1)


rownames(plot_mat)=as.data.frame(strsplit2(rownames(plot_mat),split = "\\("))[,1]
display_mat=plot_mat
display_mat[sapply(display_mat,function(x){x<25})]<-0
display_mat[sapply(display_mat,function(x){x>=25})]<-"*"
display_mat[sapply(display_mat,function(x){x==0})]<-""

pheatmap::pheatmap(plot_mat,
                     cluster_row =T,
                     cluster_cols = F,
                     border=T,
                     border_color = "black",
                     show_rownames = T,
                     show_colnames = T,
                     clustering_method = "ward.D2",
                     cellwidth =10,
                     cellheight = 10,
                     fontsize = 10,
                     cutree_rows = 5,
                     #gaps_col = seq(2,2,by=2),
                     #gaps_row = gaps.row,
                     legend_labels = "-log10 P-value",
                     display_numbers =as.matrix(display_mat),
                     #main = title,
                     #filename = pdf,
                     color=c(colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))(200)))


# 5. Overlap with CUTTag ####

ctrl_upeak=readPeakFile("ctrl_upeak_ver2.bed")
k2_upeak=readPeakFile("k2_upeak_ver2.bed")
k4_upeak=readPeakFile("k4_upeak_ver2.bed")
k5_upeak=readPeakFile("k5_upeak_ver2.bed")

cuttag_k2_upeak=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K2_L_n_11657_.bed")
cuttag_k4_upeak=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K4_L_n_6937_.bed")
cuttag_k5_upeak=readPeakFile("title_K2_L_vs_K4_L_vs_K5_L_400/title_K2_L_vs_K4_L_vs_K5_L_400_K5_L_n_5810_.bed")

Make_5_Sample_Overlaps_Clustered_Bedfile(name1 = "cuttag_k2_upeak",
                                         name2 = "ctrl_upeak",
                                         name3 = "k2_upeak",
                                         name4 = "k4_upeak",
                                         name5 = "k5_upeak",
                                         peak1 = cuttag_k2_upeak,
                                         peak2 = ctrl_upeak,
                                         peak3 = k2_upeak,
                                         peak4 = k4_upeak,
                                         peak5 = k5_upeak,maxgap = 400)

Make_5_Sample_Overlaps_Clustered_Bedfile(name1 = "cuttag_k4_upeak",
                                         name2 = "ctrl_upeak",
                                         name3 = "k2_upeak",
                                         name4 = "k4_upeak",
                                         name5 = "k5_upeak",
                                         peak1 = cuttag_k4_upeak,
                                         peak2 = ctrl_upeak,
                                         peak3 = k2_upeak,
                                         peak4 = k4_upeak,
                                         peak5 = k5_upeak,maxgap = 400)

Make_5_Sample_Overlaps_Clustered_Bedfile(name1 = "cuttag_k5_upeak",
                                         name2 = "ctrl_upeak",
                                         name3 = "k2_upeak",
                                         name4 = "k4_upeak",
                                         name5 = "k5_upeak",
                                         peak1 = cuttag_k5_upeak,
                                         peak2 = ctrl_upeak,
                                         peak3 = k2_upeak,
                                         peak4 = k4_upeak,
                                         peak5 = k5_upeak,maxgap = 400)


# '-------------------------------- ####
# Part1 ESC TSC CUTTag&K245DEG ####
# 1.peak anno ####
K2_L=readPeakFile(CUTTag_bed[9])
K4_L=readPeakFile(CUTTag_bed[11])
K5_L=readPeakFile(CUTTag_bed[13])

K2_ESC=readPeakFile(CUTTag_bed[5])
K4_ESC=readPeakFile(CUTTag_bed[6])
K5_ESC=readPeakFile(CUTTag_bed[4])

K2_L_anno=Peak_Anno(K2_L,up = 20000,down = 20000,return = "promoter_matrix")
K4_L_anno=Peak_Anno(K4_L,up = 20000,down = 20000,return = "promoter_matrix")
K5_L_anno=Peak_Anno(K5_L,up = 20000,down = 20000,return = "promoter_matrix")

# 2.DEG list ####
deg_path="D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis"

# 2.1 load TPM matrix ####
tpm_total_mat=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new.csv",header = T,row.names = 1)
head(tpm_total_mat)
colnames(tpm_total_mat)

tpm_mat_mean=Combine_Repeat(tpm_total_mat)
as.data.frame(colnames(tpm_mat_mean))

early_mat=read.csv("D:/????/Data processing/#mES_to_2c_data#/Eatly-Embryo-GSE45719_scRNAseq/Eatly-Embryo-GSE45719_scRNAseq_tpm_anno.csv",header = T,row.names = 1)
early_mat=early_mat[,c(1:256)]
#boxplot(log2(early_mat+0.1))
colnames(early_mat)
head(early_mat)
early_mat_mean=Combine_Repeat(early_mat)
head(early_mat_mean)

tpm_mat_mean=cbind(tpm_mat_mean,early_mat_mean)


k2_d2_deg=read.table(paste(deg_path,"Day2_L_iEmpty _vs_ Day2_L_iKlf2 _2/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",sep = "/"))
k4_d2_deg=read.table(paste(deg_path,"Day2_L_iEmpty _vs_ Day2_L_iKlf4 _2/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",sep = "/"))
k5_d2_deg=read.table(paste(deg_path,"Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",sep = "/"))
k2_d5_deg=read.table(paste(deg_path,"Day5_L_iEmpty _vs_ Day5_L_iKlf2 _2/Day5_L_iEmpty vs Day5_L_iKlf2_degene_list.txt",sep = "/"))
k4_d5_deg=read.table(paste(deg_path,"Day5_L_iEmpty _vs_ Day5_L_iKlf4 _2/Day5_L_iEmpty vs Day5_L_iKlf4_degene_list.txt",sep = "/"))
k5_d5_deg=read.table(paste(deg_path,"Day5_L_iEmpty _vs_ Day5_L_iKlf5 _2/Day5_L_iEmpty vs Day5_L_iKlf5_degene_list.txt",sep = "/"))

c(5,9,6,10)c(5,13,6,14)c(5,17,6,18)#D25
k2_exp_gene=rownames(tpm_mean_total_mat[apply(tpm_mean_total_mat[,c(5,9)], 1, function(x) max(x)>2),]) #D2
k4_exp_gene=rownames(tpm_mean_total_mat[apply(tpm_mean_total_mat[,c(5,13)], 1, function(x) max(x)>2),])#D2
k5_exp_gene=rownames(tpm_mean_total_mat[apply(tpm_mean_total_mat[,c(5,17)], 1, function(x) max(x)>2),])#D2

k2_d2_up=intersect(na.omit(k2_d2_deg$V2),k2_exp_gene)
k2_d2_dn=intersect(na.omit(k2_d2_deg$V1),k2_exp_gene)
k4_d2_up=intersect(na.omit(k4_d2_deg$V2),k4_exp_gene)
k4_d2_dn=intersect(na.omit(k4_d2_deg$V1),k4_exp_gene)
k5_d2_up=intersect(na.omit(k5_d2_deg$V2),k5_exp_gene)
k5_d2_dn=intersect(na.omit(k5_d2_deg$V1),k5_exp_gene)

k2_d5_up=intersect(na.omit(k2_d5_deg$V2),k2_exp_gene)
k2_d5_dn=intersect(na.omit(k2_d5_deg$V1),k2_exp_gene)
k4_d5_up=intersect(na.omit(k4_d5_deg$V2),k4_exp_gene)
k4_d5_dn=intersect(na.omit(k4_d5_deg$V1),k4_exp_gene)
k5_d5_up=intersect(na.omit(k5_d5_deg$V2),k5_exp_gene)
k5_d5_dn=intersect(na.omit(k5_d5_deg$V1),k5_exp_gene)

# 3.K245 DEG overlap ####

list=list(k2_d2_up=k2_d2_up,
          k4_d2_up=k4_d2_up,
          k5_d2_up=k5_d2_up)
Venner_plot(list)

k245_up=c(k2_d2_up,k4_d2_up,k5_d2_up)
k245_up=k245_up[!duplicated(k245_up)]
length(k245_up)

colnames(tpm_mat_mean)
c(1,5:6,9:10,13:14,17:18,19,21,28,29)#D25
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k245_up,c(1,5:6,9:10,13:14,17:18,19,21,28,29,75,76)]
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k245_up,c(1,5,9,13,17,19,21,28,29,75,76)]
dim(using_mat)
colnames(using_mat)

d=Cluster_mat_by_Kmeans(using_mat,use_col = c(2,4,6,8),set.seed = 3,kmeans =6)
anno=data.frame(kmean=d$kmean,row.names = rownames(d))
d$kmean=factor(d$kmean,levels = c(3,1,5,6,2,4))
d=d[order(d$kmean),]
pheatmap::pheatmap(d[,1:c(dim(d)[2]-1)],cluster_rows = F,cluster_cols = F,
                   scale = "row",show_rownames = F,annotation_row = anno)

Plot_clustered_heatmap_R4.2.2(d[,1:c(dim(d)[2]-1)],
                              colum_split_fun = c(1,2,2,2,2,3,3,4,4,5,5),
                              kmeans = 7,
                              name="K245_D2_up_overlap_TPM2",
                              clustering_method="complete",
                              clustered_and_order =d$kmean,
                              scale = "row",
                              set_seed = 169,
                              study_genes=c(tsc,"Rhox6","Rhox9","Dppa3","Prdm1"),
                              dont_clustered = T,
                              show_study_genes = T,
                              show_row_names = F,
                              breaks=c(-2,2),#Default
                              height = 5,#Default
                              width = 4,#Default
                              anno_font_size=8,#use_raster = T,raster_quality = 25
)

d=read.csv("K245_D2_up_overlap_TPM2_2641_.csv",header = T,row.names = 1)
d$gene=rownames(d)
colnames(d)
cluster_gene_list=Convert_Long_dataframe_to_list(d,group_index = 12,variable_index = 13)

# 3.1 GO KEGG K245 DEG overlap ####
go_bp_mat=Compare_Cluster_GO_results(cluster_gene_list,top_go_term = 5,return_plot_or_matrix = "total_matrix")
write.csv(go_bp_mat,"K245_D2_up_overlap_TPM2_2641_GO_BP.csv")
go_bp_mat=read.csv("K245_D2_up_overlap_TPM2_2641_GO_BP.csv",header = T,row.names = 1)
term_selected=c( "response to leukemia inhibitory factor","DNA methylation or demethylation","alpha-amino acid metabolic process","meiotic cell cycle","stem cell population maintenance",
                 "skin development","epidermis development","muscle system process","leukocyte migration","cell-substrate adhesion",
                 "placenta development","pattern specification process","cell fate commitment","epithelial cell development","hormone metabolic process",
                 "cellular carbohydrate metabolic process","autophagy","positive regulation of tumor necrosis factor production","regulation of neurogenesis","positive regulation of apoptotic signaling pathway",
                 "regulation of DNA-binding transcription factor activity","stem cell division","BMP signaling pathway","cell fate specification","regulation of heterochromatin formation",
                 "cell junction assembly","regulation of actin filament-based process","monosaccharide metabolic process","cell-substrate adhesion","epithelial cell development")
Plot_Compare_Cluster_results(go_bp_mat,selected_term = term_selected)


kegg_mat=Compare_Cluster_KEGG_results(cluster_gene_list,universe_gene =rownames(tpm_mean_total_mat) ,top_kegg_term = 5,qvalueCutoff = 0.5)
term_selected=c("Biosynthesis of cofactors","Thiamine metabolism","Glycine, serine and threonine metabolism",
                "MAPK signaling pathway","Calcium signaling pathway","Toxoplasmosis","Osteoclast differentiation","Ras signaling pathway",
                "Glioma","Mineral absorption","Signaling pathways regulating pluripotency of stem cells","Hippo signaling pathway","Wnt signaling pathway",
                "Metabolism of xenobiotics by cytochrome P450",
                "Cortisol synthesis and secretion","Chemical carcinogenesis","TGF-beta signaling pathway","Pantothenate and CoA biosynthesis","GnRH secretion",
                "Cell adhesion molecules","Fluid shear stress and atherosclerosis","Amino sugar and nucleotide sugar metabolism","Sphingolipid metabolism","Galactose metabolism")
Plot_Compare_Cluster_results(kegg_mat,results_type = "KEGG",selected_term = term_selected)


# 4.binding gene overlap DEG ####

k2_list=list(K2_L_anno=K2_L_anno$SYMBOL,
             k2_d25_up=c(k2_d2_up))
k2_venn=get.venn.partitions(k2_list)
k2_bind_up=k2_venn$..values..$`1`

k4_list=list(K4_L_anno=K4_L_anno$SYMBOL,
             k4_d25_up=c(k4_d2_up))
k4_venn=get.venn.partitions(k4_list)
k4_bind_up=k4_venn$..values..$`1`

k5_list=list(K5_L_anno=K5_L_anno$SYMBOL,
             k5_d25_up=c(k5_d2_up))
k5_venn=get.venn.partitions(k5_list)
k5_bind_up=k5_venn$..values..$`1`


k245_bind_up_list=list(k2_bind_up=k2_bind_up,
                       k4_bind_up=k4_bind_up,
                       k5_bind_up=k5_bind_up)
Venner_plot(k245_bind_up_list)

k245_bind_deg=c(k2_bind_up,k4_bind_up,k5_bind_up)
k245_bind_deg=k245_bind_deg[!duplicated(k245_bind_deg)]
length(k245_bind_deg)

colnames(tpm_mat_mean)
c(1,5:6,9:10,13:14,17:18,19,21,28,29)#D25
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k245_bind_deg,c(1,5,9,13,17,19,21,28,29,75,76)]
dim(using_mat)
colnames(using_mat)

d=Cluster_mat_by_Kmeans(using_mat,use_col = c(2:5),set.seed = 3,kmeans = 6)
anno=data.frame(kmean=d$kmean,row.names = rownames(d))
d$kmean=factor(d$kmean,levels = c(2,3,5,4,6,1))
d=d[order(d$kmean),]
pheatmap::pheatmap(d[,1:c(dim(d)[2]-1)],cluster_rows = F,cluster_cols = F,
                   scale = "row",show_rownames = F,annotation_row = anno)

Plot_clustered_heatmap_R4.2.2(d[,1:c(dim(d)[2]-1)],
                              colum_split_fun = c(1,2,2,2,2,3,3,4,4,5,5),
                              kmeans = 7,
                              name="K245_bind_D2_up_overlap_TPM2",
                              clustering_method="complete",
                              clustered_and_order =d$kmean,
                              scale = "row",
                              set_seed = 169,
                              study_genes=c(tsc,"Rhox6","Rhox9","Dppa3","Prdm1"),
                              dont_clustered = T,
                              show_study_genes = T,
                              show_row_names = F,
                              breaks=c(-2,2),#Default
                              height = 5,#Default
                              width = 4,#Default
                              anno_font_size=8,#use_raster = T,raster_quality = 25
)

d=read.csv("K245_bind_D2_up_overlap_TPM2_2401_.csv",header = T,row.names = 1)
d$gene=rownames(d)
cluster_gene_list=Convert_Long_dataframe_to_list(d,group_index =12,variable_index = 13)

# 4.1 GO KEGG K245 DEG overlap ####
go_bp_mat=Compare_Cluster_GO_results(cluster_gene_list,top_go_term = 5,return_plot_or_matrix = "total_matrix")
go_bp_mat=read.csv("K245_D2_up_overlap_TPM2_2641_GO_BP.csv",header = T,row.names = 1)
term_selected=c( "response to leukemia inhibitory factor","DNA methylation or demethylation","alpha-amino acid metabolic process","meiotic cell cycle","stem cell population maintenance",
                 "skin development","epidermis development","muscle system process","leukocyte migration","cell-substrate adhesion",
                 "placenta development","pattern specification process","cell fate commitment","epithelial cell development","hormone metabolic process",
                 "positive regulation of tumor necrosis factor production","autophagy","cellular carbohydrate metabolic process","small GTPase mediated signal transduction","regulation of neurogenesis",
                 "stem cell division","BMP signaling pathway","cell fate specification","regulation of DNA-binding transcription factor activity","regulation of oxidative phosphorylation",
                 "cell junction assembly","regulation of actin filament-based process","monosaccharide metabolic process","cell-substrate adhesion","epithelial cell development")
Plot_Compare_Cluster_results(go_bp_mat,selected_term =term_selected )

kegg_mat=Compare_Cluster_KEGG_results(cluster_gene_list,universe_gene =rownames(tpm_mean_total_mat) ,top_kegg_term = 5,qvalueCutoff = 0.5)
term_selected=c("Biosynthesis of cofactors","Thiamine metabolism","Glycine, serine and threonine metabolism","ABC transporters","Cysteine and methionine metabolism",
                "MAPK signaling pathway","Calcium signaling pathway","Ras signaling pathway","Osteoclast differentiation","Parathyroid hormone synthesis, secretion and action",
                "Mineral absorption","Signaling pathways regulating pluripotency of stem cells","Glioma ","Vitamin digestion and absorption","Hepatocellular carcinoma",
                "Metabolism of xenobiotics by cytochrome P450","Fatty acid degradation",
                "Chemical carcinogenesis","Cortisol synthesis and secretion","TGF-beta signaling pathway","Pantothenate and CoA biosynthesis","GnRH secretion",
                "Cell adhesion molecules","Fluid shear stress and atherosclerosis","Amino sugar and nucleotide sugar metabolism","Leukocyte transendothelial migration","Hippo signaling pathway")
Plot_Compare_Cluster_results(kegg_mat,results_type = "KEGG",selected_term = term_selected)


# 5.Plot Marker gene ####
c(1,5:6,9:10,13:14,17:18,19,21,28,29)#D25
c(5,3,9,7,13,11,17,15)#D2 L CL
c("PNT_D2_iKlfs")#D2
c("iEmpty","iEmpty","iK2","iK2","iK4","iK4","iK5","iK5")#D25
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%c(pluripotency,PGC,TSC),c(16:17,4:5,8:9,12:13)]

tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,2]+0.1))
tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))

column_split_mat=data.frame(sample=colnames(using_mat),
                            split=c("iEmpty","iEmpty","iK2","iK2","iK4","iK4","iK5","iK5"))

fin1=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,study_gene = c(pluripotency),limit_fold = c(-3,3),column_split_mat = column_split_mat,
                            title = "Pluripotency_Marker",z_score =T,return = "fin_mat")#10
#fin2=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,study_gene = c(TSC),limit_fold = c(-3,3),column_split_mat = column_split_mat,
#                                  title = "TSC_stemness_Marker",z_score =T,return = "fin_mat")#10
fin3=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,study_gene = c(PGC),limit_fold = c(-3,3),column_split_mat = column_split_mat,
                                 title = "PGC_Marker",z_score =T,return = "fin_mat")#10

using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%c(pluripotency,PGC,TSC),c(1,108,19,20,103:104,109,105:107)]
tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,2]+0.1))
tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))

column_split_mat=data.frame(sample=colnames(using_mat),
                            split=c("EpiSCs","EpiSCs","ESC","ESC","TSC","TSC","PGCLC","PGCLC"))
column_split_mat=data.frame(sample=colnames(using_mat),
                            split=c("EpiSCs","EpiSCs","ESC","ESC","PGCLC","PGCLC","BV_SC","PGC","PGC","PGC"))

fin4=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,study_gene = c(pluripotency),limit_fold = c(-3,3),column_split_mat = column_split_mat,
                                 title = "Pluripotency_Marker",z_score =T,return = "fin_mat")#10
#fin5=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,study_gene = c(TSC),limit_fold = c(-3,3),column_split_mat = column_split_mat,
#                                  title = "TSC_stemness_Marker",z_score =T,return = "fin_mat")#10
fin6=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,study_gene = c(PGC),limit_fold = c(-3,3),column_split_mat = column_split_mat,
                                 title = "PGC_Marker",z_score =T,return = "fin_mat")#10

# fin7=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,study_gene = c(K245),limit_fold = c(-3,3),column_split_mat = column_split_mat,
#                                  title = "Klf2,4,5",z_score =T)#10

fin_mat=rbind(fin1,fin3,fin4,fin6)
fin_mat$Group=factor(fin_mat$Group,levels =c("Pluripotency_Marker","PGC_Marker"))

P1=do_plot(fin_mat,title = "K245_PNT_Day2_Day5",breaks = c(-3,3))
P1=P1+facet_grid(vars(Group),vars(column_split),scales = "free",space = "free")+
  theme(strip.text.y = element_text(size=11),
        strip.background = element_rect(fill = "#EEEEEE"))
gt1=ggplotGrob(P1)
panel_x <- unique(gt1$layout$l[grepl("panel", gt1$layout$name)])[1]
# Set that width based on the number of x-axis variables, plus 0.2 because
# of the expand arguments in the scales
#gt1$widths[panel_x] <- unit(nlevels(droplevels(fin_mat$Sample)) + 0.2, "null")
# Respect needs to be true to have 'null' units match in x- and y-direction
gt1$respect <- TRUE
grid.newpage()
grid.draw(gt1)


# '-------------------------------- ####
# Part2 2C in_vivo CUTTag&K245DEG ####
# Boxplot ####
zga_gene=as.data.frame(read_excel("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation_TEs_and_ZGA_CREs/CUTTag_ZGA_CREs/early_embryos_multiomics_defined_ZGA_gene.xls"))
two_cell_specific=as.data.frame(read_excel("D:/????/Data processing/#mES_to_2c_data#/2C-Specific Genes macfarlan 2012 .xlsx"))
head(zga_gene)
head(two_cell_specific)

k245_up=c(k2_d2_up,k4_d2_up,k5_d2_up)
k245_up=k245_up[!duplicated(k245_up)]
length(k245_up)

k245_up_ol_2c=intersect(k245_up,two_cell_specific[,1])
k245_up_ol_zga=intersect(k245_up,zga_gene$symbol)

colnames(tpm_mat_mean)
c(1,5:6,9:10,13:14,17:18,19,21,28,29)#D25
c(1,5:6,9:10,13:14,17:18,19,21,28,29,75,76,
  42,67,68,92,69,93,43,70,44,71,96,45,72,97,46,
  47,73,48,
  79,77,78)#all
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k245_up,c(1,5,9,13,17)]#D25]
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k245_up,c(1,5,9,13,17,19,21,28,29,75,76)]
dim(using_mat)
colnames(using_mat)

d=Cluster_mat_by_Kmeans(using_mat,use_col = c(2:5),set.seed = 3,kmeans =6)
anno=data.frame(kmean=d$kmean,row.names = rownames(d))
d$kmean=factor(d$kmean,levels = c(3,1,5,6,2,4))
d=d[order(d$kmean),]
pheatmap::pheatmap(d[,1:c(dim(d)[2]-1)],cluster_rows = F,cluster_cols = F,
                   scale = "row",show_rownames = F,annotation_row = anno)

d1=d[d$kmean%in%c(1),]
d1$kmean=NULL
d1=Cluster_mat_by_Kmeans(d1,use_col = c(10:36),kmeans = 8)
anno=data.frame(kmean=d1$kmean,row.names = rownames(d1))
#d$kmean=factor(d$kmean,levels = c(2,3,5,4,6,1))
#d=d[order(d$kmean),]
pheatmap::pheatmap(d1[,1:c(dim(d1)[2]-1)],cluster_rows = F,cluster_cols = F,
                   scale = "row",show_rownames = F,annotation_row = anno)

rownames(d1[d1$kmean%in%c(4,5),])
intersect(rownames(d1),zga_gene$symbol)
intersect(rownames(d1[d1$kmean%in%c(2,4),]),zga_gene$symbol)
intersect(rownames(d1),two_cell_specific[,1])
intersect(rownames(d1[d1$kmean%in%c(4,7),]),two_cell_specific[,1])


using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%rownames(d[d$kmean%in%c(3),]),c(92:101,67:73,42:47)]
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k2_d2_up,c(92:101,67:73,42:47)]
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%rownames(d[d$kmean%in%c(1),]),c(92:101,67:73,42:47)]
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k4_d2_up,c(92:101,67:73,42:47)]
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%rownames(d[d$kmean%in%c(5),]),c(92:101,67:73,42:47)]
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%k5_d2_up,c(92:101,67:73,42:47)]
boxplot(log2(using_mat+0.1))

using_mat1=tpm_mat_mean[rownames(tpm_mat_mean)%in%rownames(d[d$kmean%in%c(2),]),c(92:101,67:73,42:47)]
using_mat2=tpm_mat_mean[rownames(tpm_mat_mean)%in%k2_bind_up,c(92:101,67:73,42:47)]
using_mat3=tpm_mat_mean[rownames(tpm_mat_mean)%in%rownames(d[d$kmean%in%c(3),]),c(92:101,67:73,42:47)]
using_mat4=tpm_mat_mean[rownames(tpm_mat_mean)%in%k4_bind_up,c(92:101,67:73,42:47)]
using_mat5=tpm_mat_mean[rownames(tpm_mat_mean)%in%rownames(d[d$kmean%in%c(5),]),c(92:101,67:73,42:47)]
using_mat6=tpm_mat_mean[rownames(tpm_mat_mean)%in%k5_bind_up,c(92:101,67:73,42:47)]
boxplot(log2(using_mat+0.1))



col_fun=c(colorRampPalette(brewer.pal(9,"Purples"))(10),
brewer.pal(7,"YlOrBr"),
brewer.pal(9,"Greens"))
column_split_mat=data.frame(sample=colnames(using_mat),
                               split=c(rep("in_vivo_data_SC",10),
                                 rep("in_vivo_data_XW",7),
                                       rep("in_vivo_data_GSR",6)))
p1=Plot_Box_plot(using_mat1,return = "melt_data",ylim = c(-3,8),col_fun = col_fun,title = "K2_unique",column_split_mat = column_split_mat,show_mean_median_line = F,show_median_line = T)
p2=Plot_Box_plot(using_mat2,return = "melt_data",ylim = c(-3,8),col_fun = col_fun,title = "K2_D2_UP",column_split_mat = column_split_mat,show_mean_median_line = F,show_median_line = T)
p3=Plot_Box_plot(using_mat3,return = "melt_data",ylim = c(-3,8),col_fun = col_fun,title = "K4_unique",column_split_mat = column_split_mat,show_mean_median_line = F,show_median_line = T)
p4=Plot_Box_plot(using_mat4,return = "melt_data",ylim = c(-3,8),col_fun = col_fun,title = "K4_D2_UP",column_split_mat = column_split_mat,show_mean_median_line = F,show_median_line = T)
p5=Plot_Box_plot(using_mat5,return = "melt_data",ylim = c(-3,8),col_fun = col_fun,title = "K5_unique",column_split_mat = column_split_mat,show_mean_median_line = F,show_median_line = T)
p6=Plot_Box_plot(using_mat6,return = "melt_data",ylim = c(-3,8),col_fun = col_fun,title = "K5_D2_UP",column_split_mat = column_split_mat,show_mean_median_line = F,show_median_line = T)
p1$group="K2_unique(n=302)"
p2$group="K2_D2_UP(n=1073)"
p3$group="K4_unique(n=1131)"
p4$group="K4_D2_UP(n=1850)"
p5$group="K5_unique(n=232)"
p6$group="K5_D2_UP(n=775)"
fin_p=rbind(p1,p2,p3,p4,p5,p6)
fin_p=rbind(p1,p3,p5)
title="K245_bind_D2_up_overlap_TPM2_2401_in_vivo_Boxplot_ver2"
ggplot(fin_p,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(width=0.7,position = position_dodge(0.9),alpha=0.6,
               fatten=2.5,#??λ?ߴ?ϸ
               outlier.size = 0.1,outlier.color = "gray",outlier.alpha = 0.5)+
  #geom_jitter(position = position_jitter(width = 0.2))+
  coord_cartesian(ylim=c(-3,5))+
  #coord_fixed()+
  theme_classic()+
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5),
    #panel.border = element_blank(),
    #axis.line.x.top = element_line(),
    #axis.line.y.right = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background =element_blank())+
  guides(fill=FALSE)+#????????legend
  ggtitle(title) +
  xlab(label = NULL)+
  ylab(label = "Log2(TPM+0.1)")+
  scale_fill_manual(values =col_fun)+
  geom_line(aes(x=variable,y=V4,group=1),
            size=0.7,color="red",alpha=1,
            linetype="longdash")+
  facet_grid(cols = vars(column_split),rows = vars(group),scales = "free_x",space = "free")





d=read.csv("K245_D2_up_overlap_TPM2_2641_.csv",header = T,row.names = 1)
d$gene=rownames(d)
tc_mat=d
tc_mat[tc_mat$Cluster=="C_1","Cluster"]="K2_unique"
tc_mat[tc_mat$Cluster=="C_2","Cluster"]="K4_unique"
tc_mat[tc_mat$Cluster=="C_3","Cluster"]="K5_unique"


# '-------------------------------- ####
# Study gene ####
pluripotency=c("Pou5f1","Nanog","Dppa5a","Zfp42","Nr5a2","Tfcp2l1","Tbx3",
               "Fgf5","Otx2")
PGC=c("Rhox6","Rhox9","Tfap2c","Kit","Dppa3","Prdm1","Blimp1","Nanos3","Prdm14","Ddx4")
TSC=c("Krt7","Gata2","Plet1","Hand1","Elf5","Gata3","Cdx2","Rin3")
K245=c("Klf2","Klf4","Klf5")

studygene=c(pluripotency,PGC,TSC)

ETP1=read.csv("F:/Backup_20240703/秦组/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/ETP_gene_WGCNA_ver3_4431_.csv")

TSCu=ETP1[ETP1$Cluster=="TSC",1]
PGCu=ETP1[ETP1$Cluster=="PGCLC",1]
ESCu=ETP1[ETP1$Cluster=="ESC",1]

TSCr=ETP1[grepl("TSC",ETP1$Cluster),1]
PGCr=ETP1[grepl("PGCLC",ETP1$Cluster),1]
ESCr=ETP1[grepl("ESC",ETP1$Cluster),1]





# '-------------------------------- ####
# PCA ####
tpm_total=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_total)

tpm_total=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new.csv",row.names = 1)
colnames(tpm_total)


select_sample=c(1,2,40,41,54,55,48,49)

tpm=tpm_total[,select_sample]
colnames(tpm)
using_tpm=tpm[apply(tpm, 1, function(x) max(x)>2),]


topVarGenes = head(order(rowVars(as.matrix(using_tpm)), decreasing = TRUE),2000)
using_tpm_cv=using_tpm[apply(using_tpm,1,function(x) (sd(x)/(mean(x)))^2>0.5),]
using_tpm_cv=using_tpm[topVarGenes,]

using_tpm_cv=na.omit(using_tpm_cv)
#using_tpm_cv=using_tpm_cv[apply(using_tpm_cv, 1, function(x) !sd(x)==0),]
dim(using_tpm_cv)
#head(using_tpm_cv)
#using_tpm_cv=expr_batch_mean_mat[rownames(expr_batch_mean_mat)%in%geneset,select_sample]


# gene.pca <- prcomp(expr_batch_mean_PCA_mat, scale. = T)
# pca_result=data.frame(condition=colnames(expr_batch_mean_PCA_mat),
#                                  PC1=gene.pca$rotation[,1],
#                                  PC2=gene.pca$rotation[,2],
#                                  PC3=gene.pca$rotation[,3])

gene.pca <- FactoMineR::PCA(t(using_tpm_cv), ncp = 10, scale.unit = T, graph = FALSE)
pca_result=data.frame(PC1_percentVar=round(gene.pca$eig[1,2],1),
                      PC2_percentVar=round(gene.pca$eig[2,2],1),
                      PC3_percentVar=round(gene.pca$eig[3,2],1),
                      #condition=strsplit2(colnames(using_tpm_cv),split = "_rep")[,1],
                      condition=colnames(using_tpm_cv),
                      PC1=gene.pca$ind$coord[,1],
                      PC2=gene.pca$ind$coord[,2],
                      PC3=gene.pca$ind$coord[,3])
pca_result$condition=factor(pca_result$condition,levels = pca_result$condition)
fig=plotly::plot_ly(pca_result,x=pca_result$PC1,y=pca_result$PC2,z=pca_result$PC3,color = pca_result$condition)
fig




set_df=data.frame(
  col_fun=c("grey","grey","#33A02C","#33A02C",
            "#1F78B4","#1F78B4","#FF7F00", "#FF7F00"
  ),
  shape_fun=c(
              rep(21,8)
              #21,21,
              #11,21,
              #21,21,4
  ),
  condition=pca_result$condition)


pca_result=readRDS("E:/#Data processing#/8CLC_DDi/e4CL_zygote-like/pca_result.rds")

set_df=data.frame(
  col_fun=rep("black",21),
  shape_fun=c(
    rep(21,21)
    #21,21,
    #11,21,
    #21,21,4
  ),
  condition=pca_result$condition)

ggplot(pca_result, aes(PC1, PC2,fill=condition,shape=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",pca_result$PC1_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",pca_result$PC2_percentVar[1],"% variance")) + 
  coord_fixed()+
  #scale_x_continuous(limits = c(min,max))+ #?޸?????
  #scale_y_continuous(limits = c(min,max))+
  #labs(title = name)+
  scale_fill_manual(values =set_df$col_fun )+
  scale_shape_manual(values = set_df$shape_fun)+
  theme_bw() +
  theme(#????һ??û??Ҫ????
    legend.background=element_blank(), 
    legend.key=element_blank(),
    #legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),#????????
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+geom_text_repel(size=2.5, 
                    nudge_y = 0,
                    aes(label=pca_result$condition,x=pca_result$PC1,y= pca_result$PC2),
                    colour= "#000000" ,
                    #arrow = arrow(length=unit(0.005, "npc")),
                    force = 1,
                    max.overlaps = 100,
                    force_pull=0,
                    box.padding = 0.2, #?ֵ????ľ???
                    point.padding = 0.33, #?ֵ????ľ??룬????Χ?Ŀհ׿???
                    min.segment.length = 0)



# '-------------------------------- ####
# K245_D25_activate have bind ####
tpm_total=read.csv("F:/Backup_20240703/秦组/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_total)

d2_K5_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d2_K2_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _2/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",header = T)
d2_K4_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _2/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",header = T)


d2_K5_l_up_2=intersect(d2_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K2_l_up_2=intersect(d2_K2_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K4_l_up_2=intersect(d2_K4_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))

d2_K5_l_dn_2=intersect(d2_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K2_l_dn_2=intersect(d2_K2_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K4_l_dn_2=intersect(d2_K4_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))





d5_K5_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf5 _2/Day5_L_iEmpty vs Day5_L_iKlf5_degene_list.txt",header = T)
d5_K2_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf2 _2/Day5_L_iEmpty vs Day5_L_iKlf2_degene_list.txt",header = T)
d5_K4_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf4 _2/Day5_L_iEmpty vs Day5_L_iKlf4_degene_list.txt",header = T)

d5_K5_l_up_2=intersect(d5_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))
d5_K5_l_dn_2=intersect(d5_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))
d5_K2_l_up_2=intersect(d5_K2_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))
d5_K2_l_dn_2=intersect(d5_K2_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))
d5_K4_l_up_2=intersect(d5_K4_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))
d5_K4_l_dn_2=intersect(d5_K4_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))


#####have bind####
CUTTag_bed=list.files("D:/????/Data processing/#Klf5 Figure#/All_CUTTag_bedfile/q0_sorted/",full.names = T)

K2_L=readPeakFile(CUTTag_bed[9])
K4_L=readPeakFile(CUTTag_bed[11])
K5_L=readPeakFile(CUTTag_bed[13])

K2_L_anno=Peak_Anno(K2_L,up = 20000,down =20000)
K4_L_anno=Peak_Anno(K4_L,up = 20000,down =20000)
K5_L_anno=Peak_Anno(K5_L,up = 20000,down =20000)

K5_L_anno


Make_Gene_bed=
  function(deg_list,anno_peak,name="gene_list"){

homer_mm10_TSS_2k_bed=read.table("D:/????/Data processing/Data_base/homer_mm10_TSS_2k.bed")

K5_active_mat=bitr(deg_list,fromType = "SYMBOL",toType = "REFSEQ",OrgDb = org.Mm.eg.db)
K5_active_matNM=K5_active_mat[grepl("NM",K5_active_mat$REFSEQ),]
K5_active_matNR=K5_active_mat[grepl("NR",K5_active_mat$REFSEQ),]
K5_active_mat=rbind(K5_active_matNM,K5_active_matNR)
K5_active_bed=Btable_insert_into_Atable(K5_active_mat,homer_mm10_TSS_2k_bed,2,5)
K5_active_bed=na.omit(K5_active_bed)
K5_active_bed=K5_active_bed[!duplicated(K5_active_bed$SYMBOL),]

K5_L_anno_bind_bed=K5_active_bed[K5_active_bed$SYMBOL%in%intersect(anno_peak,K5_active_bed$SYMBOL),]
K5_L_anno_bind_bed$deepTools_group="Bind"

K5_L_anno_nobind_bed=K5_active_bed[K5_active_bed$SYMBOL%in%setdiff(K5_active_bed$SYMBOL,anno_peak),]
K5_L_anno_nobind_bed$deepTools_group="No_Bind"

K5_L_anno_bed_fin=rbind(K5_L_anno_bind_bed,K5_L_anno_nobind_bed)
K5_L_anno_bed_fin=K5_L_anno_bed_fin[,c(3:5,1,8)]
colnames(K5_L_anno_bed_fin)=c("#Chr","Start","End","Genes","deepTools_group")

out=setdiff(deg_list,K5_L_anno_bed_fin$Genes)
out_mat=data.frame(a="chr",b=1,end=1,c=out,d="No_Bind")
colnames(out_mat)=c("#Chr","Start","End","Genes","deepTools_group")

K5_L_anno_bed_fin=rbind(K5_L_anno_bed_fin,out_mat)
dim(K5_L_anno_bed_fin)
head(K5_L_anno_bed_fin)

Output_bedfile(K5_L_anno_bed_fin,name=name,colname = T)

return(K5_L_anno_bed_fin)
}

K5_L_anno_bed_fin=Make_Gene_bed(d2_K5_l_up_2,K5_L_anno,"K5_L_anno_bed_fin")
K2_L_anno_bed_fin=Make_Gene_bed(d2_K2_l_up_2,K2_L_anno,"K2_L_anno_bed_fin")
K4_L_anno_bed_fin=Make_Gene_bed(d2_K4_l_up_2,K4_L_anno,"K4_L_anno_bed_fin")

table(K2_L_anno_bed_fin$deepTools_group)
table(K4_L_anno_bed_fin$deepTools_group)



#####Mfuzz Clusters up#####
library(Mfuzz)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)
library("stringr")
library(ggplot2)

setwd("E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3/")


tpm_total=read.csv("F:/Backup_20240703/秦组/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_total)

d2_K5_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d2_K2_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _2/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",header = T)
d2_K4_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _2/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",header = T)


d2_K5_l_up_2=intersect(d2_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K2_l_up_2=intersect(d2_K2_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K4_l_up_2=intersect(d2_K4_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))

d2_K5_l_dn_2=intersect(d2_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K2_l_dn_2=intersect(d2_K2_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K4_l_dn_2=intersect(d2_K4_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))

using_mat=tpm_total[rownames(tpm_total)%in%c(d2_K5_l_up_2,
                                             d2_K2_l_up_2,
                                             d2_K4_l_up_2),c(16,4,8,12)]



using_mat=read.csv("E:/Klf5_paper_submission/NC_revision/Analysis part/R1_comment3/K245_D2_up_overlap_TPM2_2641_.csv",row.names = 1)
using_mat=using_mat[using_mat$Cluster%in%c("C_2","C_3"),]

rpkm_mat_rpkm_cv=using_mat[,1:4]


rpkm_mat_rpkm_cv <- rpkm_mat_rpkm_cv[apply(rpkm_mat_rpkm_cv, 1, sd) > 0, ]
rpkm_mat_rpkm_cv=na.omit(rpkm_mat_rpkm_cv)

sample_info <- data.frame(
  row.names = colnames(rpkm_mat_rpkm_cv),
  sample =  colnames(rpkm_mat_rpkm_cv)# 根据你的实验设计修改
)

# 创建特征信息（featureData）
# 这是一个数据框，行数等于基因数，行名必须与表达矩阵的行名匹配。
# 通常可以先创建一个空的或包含基因ID的信息
feature_info <- data.frame(
  row.names = rownames(rpkm_mat_rpkm_cv),
  GeneID = rownames(rpkm_mat_rpkm_cv)
)

# 创建 ExpressionSet 对象
eset <- ExpressionSet(assayData = as.matrix(rpkm_mat_rpkm_cv),
                      phenoData = new("AnnotatedDataFrame", data = sample_info),
                      featureData = new("AnnotatedDataFrame", data = feature_info))

# 示例：过滤掉在所有样本中标准差为0的基因（即没有变化的基因）


# 方法1：使用 scale 函数手动标准化
exprs_matrix <- exprs(eset)
exprs_standardized <- t(scale(t(exprs_matrix))) # 对行进行标准化

# 创建新的 ExpressionSet
eset_s <- eset
exprs(eset_s) <- exprs_standardized


# # 设置一个可能的聚类数量范围，例如从 2 到 20
# c_range <- 2:20
# 
# # 计算不同 c 值下的最小距离
# # 此步骤可能较慢，取决于数据量和c_range的大小
# m <- mestimate(eset_s) # 先估算模糊化参数 m
# Dmin_vals <- c()
# for (i in c_range) {
#   cl <- mfuzz(eset_s, c = i, m = m)
#   Dmin_vals <- c(Dmin_vals, cl$Dmin)
# }

# 绘制 “肘部法则” 图
#plot(c_range, Dmin_vals, type = "b", pch = 19, xlab = "Number of clusters (c)", ylab = "Dmin")
# 选择图中“拐点”（elbow point）对应的 c 值

# 估算模糊加权指数 m (通常使用默认估算值即可，结果对m不敏感)
m <- mestimate(eset_s)
print(paste("The estimated m is:", m))

# 设置随机种子以保证结果可重现
set.seed(123)

# 执行 Mfuzz 聚类
cl <- mfuzz(eset_s, c = 30, m = m) # 将 c 替换为你选择的值 #8


# 获取标准化后的表达矩阵
expression_data <- rpkm_mat_rpkm_cv

# 获取聚类结果
cluster_results <- data.frame(
  GeneID = rownames(expression_data),
  Cluster = cl$cluster,
  Max_Membership = apply(cl$membership, 1, max)
)

# 按聚类号和成员值排序
cluster_results <- cluster_results %>%
  arrange(Cluster, desc(Max_Membership))

# 重新排序表达矩阵以匹配聚类顺序
expression_ordered <- expression_data[cluster_results$GeneID, ]

# 创建行注释数据框
row_annotation <- data.frame(
  Cluster = factor(paste0("Cluster_", cluster_results$Cluster)),
  Membership = cluster_results$Max_Membership
)
rownames(row_annotation) <- cluster_results$GeneID

# 为不同聚类设置颜色
cluster_colors <- brewer.pal(n = length(unique(cl$cluster)), name = "Set3")
names(cluster_colors) <- unique(paste0("Cluster_", unique(cl$cluster)))

# 成员值颜色梯度
mem_colors <- colorRampPalette(c("white", "darkred"))(100)

annotation_colors <- list(
  Cluster = cluster_colors,
  Membership = mem_colors
)


#d=Cluster_mat_by_Kmeans(rpkm_mat_rpkm_cv,kmeans = 8,set.seed = 1)#auto 8
#anno=data.frame(row.names = rownames(d),kmean=d$kmean)

#d$kmean=factor(d$kmean,levels = c(3:1))
#d=d[order(d$kmean),]

expression_ordered$cluster=strsplit2(row_annotation$Cluster,"_")[,2]
plotmat=expression_ordered
plotmat$cluster=factor(plotmat$cluster,levels=c(6,8,4,1,5,2,7,3))
plotmat$cluster=factor(plotmat$cluster,levels=c(5,7,13,21,
                                                4,15,17,18,25,
                                                28,11,
                                                3,9,23,
                                                14,22,26,
                                                1,2,8,10,12,19,20,24,29,27,
                                                6,16,30))

plotmat=plotmat[order(plotmat$cluster),]


plotmat$cluster[plotmat$cluster %in% c(5,7,13,21)] <- "C_1"
plotmat$cluster[plotmat$cluster %in% c( 4,15,17,18,25)] <- "C_2"
plotmat$cluster[plotmat$cluster %in% c(28,11)] <- "C_3"
plotmat$cluster[plotmat$cluster %in% c(3,9,23)] <- "C_4"
plotmat$cluster[plotmat$cluster %in% c(14,22,26)] <- "C_5"
plotmat$cluster[plotmat$cluster %in% c(1,2,8,10,12,19,20,24,29,27)] <- "C_6"
plotmat$cluster[plotmat$cluster %in% c(6,16,30)] <- "C_7"
plotmat$cluster=factor(plotmat$cluster,levels=c("C_1","C_2","C_3","C_4",
                                                "C_5","C_6","C_7"))
plotmat=plotmat[order(plotmat$cluster),]

row_annotation=data.frame(cluster=plotmat$cluster,row.names = rownames(plotmat))


plotmat=order_fun(plotmat,1:4)


bk=seq(-1.5,1.5,by=0.1)
pheatmap::pheatmap(#d[,1:c(dim(d)[2]-1)],
  plotmat[,1:4],scale = "row",
  cluster_cols = F,cluster_rows = F,border=F,
  breaks = bk,show_rownames = F,
  annotation_row = row_annotation,
  cellwidth = 30,fontsize = 20,raster =T,
  #display_numbers = round(as.matrix(using_tpm),digits = 1),
  color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
             colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))



Plot_clustered_heatmap_R4.2.2(plotmat[,1:4],
                              colum_split_fun = c(1,1,1,1),
                              kmeans = 7,
                              name="K245_D2_up_overlap_TPM2",
                              clustering_method="complete",
                              clustered_and_order =plotmat$cluster,
                              scale = "row",
                              set_seed = 169,
                              #study_genes=c(tsc,"Rhox6","Rhox9","Dppa3","Prdm1"),
                              dont_clustered = T,
                              show_study_genes = F,
                              show_row_names = F,
                              breaks=c(-1.5,1.5),#Default
                              height = 5,#Default
                              width = 2.2,#Default
                              anno_font_size=8,use_raster = T,raster_quality = 4
)


plotmat=read.csv("K245_D2_up_overlap_TPM2_2641_.csv",row.names = 1)
deglist=split(rownames(plotmat), plotmat$Cluster)

#GO_function(rownames(plotmat[plotmat$cluster==6,]),simplify_ver2 = T,name = "K5_u")

Study_Gene_Distribution_in_gene_list(c(study_gene,PGC_all,TSC_all,"Eomes","Notch1","Map2",
                                       "Wasf2","Anxa2","Arpc1b","Krt7","Gata2","Bmp4","Stat3",
                                       "Rin3","Tead4","Esrrb","Id1","Id2","Dppa2",
                                       "Rec8","Stra8","Hormad2",
                                       "Krt14","Krt18","Sppr2a2","Itga6","Sox9","Fgf8"),deglist)



term=rev(c("placenta development","cell fate commitment","somatic stem cell population maintenance",
           "positive regulation of meiotic nuclear division",
           "organic anion transport","epithelial cell development","monosaccharide metabolic process",
           "stem cell division","endodermal cell differentiation","gastrulation","regulation of DNA-binding transcription factor activity",
           "BMP signaling pathway",
  "cellular response to leukemia inhibitory factor","response to retinoic acid","lipid catabolic process",
  "meiotic cell cycle",
           "regulation of neurogenesis","regulation of cellular component size","developmental maturation",
  "skin development","muscle system process","actin filament organization",
  "mesenchyme development","cell junction assembly","Wnt signaling pathway"
))#K2,5 up d5 1236 597 689  

term=rev(c("hexose metabolic process","organophosphate catabolic process","regulation of autophagy",
           "developmental growth involved in morphogenesis","regulation of neurogenesis","calcium ion transport",
           "cell fate specification","Tgf-β receptor superfamily signaling pathway","stem cell division",
           "cell junction assembly","epidermis development","skin development",
           "cellular response to leukemia inhibitory factor","cellular response to retinoic acid","meiotic cell cycle",
           "skin development","epidermal cell differentiation","cell-substrate adhesion",
           "placenta development","response to BMP","positive regulation of meiotic nuclear division"))
term=term[!duplicated(term)]

Compare_Cluster_GO_results(deglist,simplify_ver2 = T,
                           custom_design_go_term = term,
                           pvalue_max_limit = 10,
                           input_go_cluster_resluts_pathlist = list.files(path = "Klfs_up_GO_ver2/",pattern = "BP_Simplify_go.csv",full.names = T),
                           top_go_term = 5)
#6,


#####Mfuzz Clusters dn#####
library(Mfuzz)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(limma)
library(ComplexHeatmap)
library(circlize)

tpm_total=read.csv("F:/Backup_20240703/秦组/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_total)

d2_K5_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d2_K2_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _2/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",header = T)
d2_K4_l_2=read.table("F:/Backup_20240703/秦组/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _2/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",header = T)


d2_K5_l_up_2=intersect(d2_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K2_l_up_2=intersect(d2_K2_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K4_l_up_2=intersect(d2_K4_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))

d2_K5_l_dn_2=intersect(d2_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K2_l_dn_2=intersect(d2_K2_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K4_l_dn_2=intersect(d2_K4_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))

using_mat=tpm_total[rownames(tpm_total)%in%c(d2_K5_l_dn_2,
                                             d2_K2_l_dn_2,
                                             d2_K4_l_dn_2),c(16,4,8,12)]
dim(using_mat)

rpkm_mat_rpkm_cv=using_mat


rpkm_mat_rpkm_cv <- rpkm_mat_rpkm_cv[apply(rpkm_mat_rpkm_cv, 1, sd) > 0, ]
rpkm_mat_rpkm_cv=na.omit(rpkm_mat_rpkm_cv)

sample_info <- data.frame(
  row.names = colnames(rpkm_mat_rpkm_cv),
  sample =  colnames(rpkm_mat_rpkm_cv)# 根据你的实验设计修改
)

# 创建特征信息（featureData）
# 这是一个数据框，行数等于基因数，行名必须与表达矩阵的行名匹配。
# 通常可以先创建一个空的或包含基因ID的信息
feature_info <- data.frame(
  row.names = rownames(rpkm_mat_rpkm_cv),
  GeneID = rownames(rpkm_mat_rpkm_cv)
)

# 创建 ExpressionSet 对象
eset <- ExpressionSet(assayData = as.matrix(rpkm_mat_rpkm_cv),
                      phenoData = new("AnnotatedDataFrame", data = sample_info),
                      featureData = new("AnnotatedDataFrame", data = feature_info))

# 示例：过滤掉在所有样本中标准差为0的基因（即没有变化的基因）


# 方法1：使用 scale 函数手动标准化
exprs_matrix <- exprs(eset)
exprs_standardized <- t(scale(t(exprs_matrix))) # 对行进行标准化

# 创建新的 ExpressionSet
eset_s <- eset
exprs(eset_s) <- exprs_standardized


# # 设置一个可能的聚类数量范围，例如从 2 到 20
# c_range <- 2:20
# 
# # 计算不同 c 值下的最小距离
# # 此步骤可能较慢，取决于数据量和c_range的大小
# m <- mestimate(eset_s) # 先估算模糊化参数 m
# Dmin_vals <- c()
# for (i in c_range) {
#   cl <- mfuzz(eset_s, c = i, m = m)
#   Dmin_vals <- c(Dmin_vals, cl$Dmin)
# }

# 绘制 “肘部法则” 图
#plot(c_range, Dmin_vals, type = "b", pch = 19, xlab = "Number of clusters (c)", ylab = "Dmin")
# 选择图中“拐点”（elbow point）对应的 c 值

# 估算模糊加权指数 m (通常使用默认估算值即可，结果对m不敏感)
m <- mestimate(eset_s)
print(paste("The estimated m is:", m))

# 设置随机种子以保证结果可重现
set.seed(123)

# 执行 Mfuzz 聚类
cl <- mfuzz(eset_s, c = 30, m = m) # 将 c 替换为你选择的值


# 获取标准化后的表达矩阵
expression_data <- rpkm_mat_rpkm_cv

# 获取聚类结果
cluster_results <- data.frame(
  GeneID = rownames(expression_data),
  Cluster = cl$cluster,
  Max_Membership = apply(cl$membership, 1, max)
)

# 按聚类号和成员值排序
cluster_results <- cluster_results %>%
  arrange(Cluster, desc(Max_Membership))

# 重新排序表达矩阵以匹配聚类顺序
expression_ordered <- expression_data[cluster_results$GeneID, ]

# 创建行注释数据框
row_annotation <- data.frame(
  Cluster = factor(paste0("Cluster_", cluster_results$Cluster)),
  Membership = cluster_results$Max_Membership
)
rownames(row_annotation) <- cluster_results$GeneID

# 为不同聚类设置颜色
cluster_colors <- brewer.pal(n = length(unique(cl$cluster)), name = "Set3")
names(cluster_colors) <- unique(paste0("Cluster_", unique(cl$cluster)))

# 成员值颜色梯度
mem_colors <- colorRampPalette(c("white", "darkred"))(100)

annotation_colors <- list(
  Cluster = cluster_colors,
  Membership = mem_colors
)


#d=Cluster_mat_by_Kmeans(rpkm_mat_rpkm_cv,kmeans = 8,set.seed = 1)#auto 8
#anno=data.frame(row.names = rownames(d),kmean=d$kmean)

#d$kmean=factor(d$kmean,levels = c(3:1))
#d=d[order(d$kmean),]

expression_ordered$cluster=strsplit2(row_annotation$Cluster,"_")[,2]
plotmat=expression_ordered
plotmat$cluster=factor(plotmat$cluster,levels=c(3,5,1,2,4,6))
plotmat=plotmat[order(plotmat$cluster),]
plotmat$cluster=factor(plotmat$cluster,levels=c(1,3,8,11,16,20,29,
                                                5,10,12,15,24,27,
                                                6,7,17,18,22,23,25,28,
                                                26,
                                                2,13,14,
                                                9,19,
                                                4))

plotmat=plotmat[order(plotmat$cluster),]


plotmat$cluster[plotmat$cluster %in% c(1,8,11,16,20,29)] <- "C_1"
plotmat$cluster[plotmat$cluster %in% c(5,10,12,15,24,27,3)] <- "C_2"
plotmat$cluster[plotmat$cluster %in% c(6,7,17,18,22,23,25,28,30)] <- "C_3"
plotmat$cluster[plotmat$cluster %in% c(26)] <- "C_4"
plotmat$cluster[plotmat$cluster %in% c(2,13,14)] <- "C_5"
plotmat$cluster[plotmat$cluster %in% c(9,19)] <- "C_6"
plotmat$cluster[plotmat$cluster %in% c(4,21)] <- "C_7"
plotmat$cluster=factor(plotmat$cluster,levels=c("C_1","C_2","C_3","C_4",
                                                "C_5","C_6","C_7"))
plotmat=plotmat[order(plotmat$cluster),]

row_annotation=data.frame(cluster=plotmat$cluster,row.names = rownames(plotmat))


plotmat=order_fun(plotmat,1:4)


bk=seq(-1.5,1.5,by=0.1)
pheatmap::pheatmap(#d[,1:c(dim(d)[2]-1)],
  plotmat[,1:4],scale = "row",
  cluster_cols = F,cluster_rows = F,border=F,
  breaks = bk,show_rownames = F,
  annotation_row = row_annotation,
  cellwidth = 30,fontsize = 20,raster =T,
  #display_numbers = round(as.matrix(using_tpm),digits = 1),
  color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
             colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))


Plot_clustered_heatmap_R4.2.2(plotmat[,1:4],
                              colum_split_fun = c(1,1,1,1),
                              kmeans = 7,
                              name="K245_D2_up_overlap_TPM2",
                              clustering_method="complete",
                              clustered_and_order =plotmat$cluster,
                              scale = "row",
                              set_seed = 169,
                              #study_genes=c(tsc,"Rhox6","Rhox9","Dppa3","Prdm1"),
                              dont_clustered = T,
                              show_study_genes = F,
                              show_row_names = F,
                              breaks=c(-1.5,1.5),#Default
                              height = 5,#Default
                              width = 2.2,#Default
                              anno_font_size=8,use_raster = T,raster_quality = 4
)



plotmat=read.csv("K245_D2_up_overlap_TPM2_1382_.csv",row.names = 1)
deglist=split(rownames(plotmat), plotmat$Cluster)

#GO_function(rownames(plotmat[plotmat$cluster==6,]),simplify_ver2 = T,name = "K5_u")

Study_Gene_Distribution_in_gene_list(c(wnt,"Dab1","Hes3","Hoxa1","Rspo1","Frzb","Sfrp1Btg2","Sox4","Flt1",
                                       "Cxcl12","Col7a1","Tgfb2Kcnh2","Scn3b","Cdh2",
                                       "Irx5","Twist1","Ets1Itga5","Cald1","Tmsb4x","Fzd3",
                                       "Dusp10","Bcar1","Ksr1Neurod2","Zic1","Rora","Id1","Id2","Id3","Id4","Hoxb1","Gata4",
                                       "Fzd7","Sfrp2","Mdk","Igfbp4","Gdf11","Tgfbr3"),deglist)


term=rev(c("regulation of epithelial cell proliferation","embryonic skeletal system morphogenesis","ossification",
           "regulation of canonical Wnt signaling pathway","SMAD protein signal transduction",
           "positive regulation of nervous system development","appendage development","limb development","cardiac chamber development",
           "regulation of cell-substrate adhesion","chemotaxis","negative regulation of phosphorus metabolic process",
           "regulation of protein catabolic ","regulation of ERK1 and ERK2 cascade",
           "regulation of membrane potential","cell junction assembly","heart morphogenesis",
           "mesenchyme development","endothelium development",
           
           "regulation of neurogenesis","regulation of vasculature development","regulation of wound healing",
           "skin development",
           
           "forebrain development","pattern specification process","gland development",
           'embryonic organ morphogenesis',"Wnt signaling pathway"))#K2,5 up d5 1236 597 689  

term=rev(c("regulation of synapse structure or activity","pattern specification process","gland development",
           "pattern specification process","ossification","cell-cell signaling by wnt",
          "regulation of cell-substrate adhesion","hormone secretion","regulation of neurogenesis",
           "microtubule-based movement","neuron migration","forebrain development",
           "muscle organ development","positive regulation of epithelial cell migration","regulation of Wnt signaling pathway",
           "regulation of epithelial cell proliferation","endothelial cell proliferation","kidney morphogenesis",
           "regulation of neurogenesis","cellular response to vascular endothelial growth factor stimulus","monocyte chemotaxis"))


Compare_Cluster_GO_results(deglist,simplify_ver2 = T,
                           custom_design_go_term = term,
                           pvalue_max_limit = 12,
                           input_go_cluster_resluts_pathlist = list.files(path = "Klfs_down_GO_ver2/",pattern = "",full.names = T),
                           top_go_term = 5)


 #####Heatmap####
tpm_total_rep=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new.csv",row.names = 1)
colnames(tpm_total_rep)


K5_L_anno_bed_fin=read.table('K5_L_anno_fin.bed',header = F)
K2_L_anno_bed_fin=read.table('K2_L_anno_fin.bed',header = F)
K4_L_anno_bed_fin=read.table('K4_L_anno_fin.bed',header = F)

k5_active_tpm_mat=Btable_insert_into_Atable(K5_L_anno_bed_fin,tpm_total_rep[,c(34,35,26,27)],4)
k5_active_tpm_mat=as.data.frame(apply(k5_active_tpm_mat[,6:9], 2, as.numeric))
rownames(k5_active_tpm_mat)=K5_L_anno_bed_fin$Genes
head(k5_active_tpm_mat)

k2_active_tpm_mat=Btable_insert_into_Atable(K2_L_anno_bed_fin,tpm_total_rep[,c(34,35,10,11)],4)
k2_active_tpm_mat=as.data.frame(apply(k2_active_tpm_mat[,6:9], 2, as.numeric))
rownames(k2_active_tpm_mat)=K2_L_anno_bed_fin$Genes
head(k2_active_tpm_mat)

k4_active_tpm_mat=Btable_insert_into_Atable(K4_L_anno_bed_fin,tpm_total_rep[,c(34,35,18,19)],4)
k4_active_tpm_mat=as.data.frame(apply(k4_active_tpm_mat[,6:9], 2, as.numeric))
rownames(k4_active_tpm_mat)=K4_L_anno_bed_fin$Genes
head(k4_active_tpm_mat)


bk=seq(-1.5,1.5,by=0.01)
p1=pheatmap::pheatmap(k4_active_tpm_mat,show_rownames = F,cellwidth = 25,cellheight = 0.2,main = paste0("iKlf4_D2_Active_Genes(n=",dim(K4_L_anno_bed_fin)[1],")"),
                      cluster_cols = F,cluster_rows = F,breaks = bk,scale = "row",gaps_row = table(K4_L_anno_bed_fin$V5)[1],
                      color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                                 colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))



#####GO####

deglist=list(d5_K5_l_up_2=d5_K5_l_up_2,
             d5_K2_l_up_2=d5_K2_l_up_2)
venn=get.venn.partitions(deglist)
deglist=venn$..values..
names(deglist)=c("Overlap","iKlf2_D5_unique","iKlf5_D5_unique")

deglist=list(d5_K5_l_dn_2=d5_K5_l_dn_2,
             d5_K2_l_dn_2=d5_K2_l_dn_2)
venn=get.venn.partitions(deglist)
deglist=venn$..values..
names(deglist)=c("Overlap","iKlf2_D5_unique","iKlf5_D5_unique")


deglist=list(d2_K5_l_up_2=intersect(d2_K5_l_up_2,K5_L_anno),
             d5_K5_l_up_2=d5_K5_l_up_2)
venn=get.venn.partitions(deglist)
deglist=venn$..values..
names(deglist)=c("Overlap","iKlf5_D5_unique","iKlf5_D2_unique")




deglist=list(d2_K2_l_dn_2=d2_K2_l_dn_2,
             d2_K5_l_dn_2=d2_K5_l_dn_2,
             d2_rfp_up=d2_rfp_v_epi_l_up_2)

deglist=list(d5_K2_l_dn_2=d5_K2_l_dn_2,
             d5_K5_l_dn_2=d5_K5_l_dn_2)

deglist=list(d5_K2_l_up_2=d5_K2_l_up_2,
             d5_K5_l_up_2=d5_K5_l_up_2)


deglist=list(d2_rfp_up=d2_rfp_v_epi_l_up_2,
             d5_rfp_up=d5_rfp_v_epi_l_up_2)

deglist=list(d2_rfp_dn=d2_rfp_v_epi_l_dn_2,
             d5_rfp_dn=d5_rfp_v_epi_l_dn_2)

deglist=list(d2_rfp_dn=d2_rfp_v_epi_l_dn_2,
             d5_rfp_dn=d5_rfp_v_epi_l_dn_2)

deglist=list(iKlf5_D2=K5_L_anno_bed_fin[K5_L_anno_bed_fin$V5=="Bind",4],
             iKlf2_D2=K2_L_anno_bed_fin[K2_L_anno_bed_fin$V5=="Bind",4])

deglist=list(d2_K5_l_dn_2=intersect(d2_K5_l_dn_2,K5_L_anno),
             d2_K2_l_dn_2=intersect(d2_K2_l_dn_2,K2_L_anno),
             d2_rfp_up=d2_rfp_v_epi_l_up_2)

deglist=list(d2_K5_l_dn_2=intersect(d2_K5_l_dn_2,K5_L_anno),
             d2_K2_l_dn_2=intersect(d2_K2_l_dn_2,K2_L_anno))

deglist=list(d2_K5_l_up_2=intersect(d2_K5_l_up_2,K5_L_anno),
             d2_K2_l_up_2=intersect(d2_K2_l_up_2,K2_L_anno))


venn=get.venn.partitions(deglist)
deglist=venn$..values..
names(deglist)=venn$..count..
Venner_plot(deglist)

c=GO_function(d2_K5_l_up_2,simplify = F)
Plot_GO_result(c,display_top = 10)

c556=deglist$`556`

deglist=list(overlap=deglist$`395`,
  K5_unique=deglist$`170`,
  K2_unqiue=deglist$`520`,
  none_bind=c556
  )

deglist=list(D2_UP=1,
             D5_UP=1,
             D2_DN=1,
             D5_DN=1)


Compare_Cluster_GO_results(deglist,top_go_term = 10,simplify_ver2 = T)


term=rev(c("placenta development","epidermis development","pattern specification process",
           "mesenchyme development","response to BMP","muscle organ development",
           "regulation of hemopoiesis","regulation of DNA-binding transcription factor activity",
           "regulation of stem cell proliferation","positive regulation of meiotic nuclear division",
           "axonogenesis","cell junction assembly",
           "cellular response to leukemia inhibitory factor",
           "DNA methylation or demethylation",
           "regulation of neurogenesis","stem cell differentiation"
))

term=rev(c("urogenital system development","ossification",
           "muscle cell differentiation","Wnt signaling pathway",
           "axonogenesis","forebrain development","pattern specification process",'sensory organ morphogenesis',
           "mesenchyme development"
))#RFP up D2,5 #32,39

term=rev(c("urogenital system development","ossification",
           "muscle cell differentiation","Wnt signaling pathway",
           "regulation of ERK1 and ERK2 cascade","response to BMP",
           "axonogenesis","forebrain development","pattern specification process",'sensory organ morphogenesis',
           "mesenchyme development",
           "hexose metabolic process",
           "generation of precursor metabolites and energy",
           "response to fibroblast growth factor",
           "cellular amino acid metabolic process",
           "stem cell differentiation","mitotic cell cycle phase transition",
           "stem cell population maintenance"
))#RFP dn D2,5

term=rev(c("axonogenesis","cell-substrate adhesion","urogenital system development","ossification",
           "pattern specification process","muscle cell differentiation",
           "mesenchymal cell differentiation","Wnt signaling pathway","regulation of epithelial cell proliferation",
           "regulation of ERK1 and ERK2 cascade","response to BMP",
           "hexose metabolic process","generation of precursor metabolites and energy",
           "nucleoside diphosphate phosphorylation",
           "stem cell differentiation","small molecule catabolic process",
           "cellular response to leukemia inhibitory factor",
           "response to fibroblast growth factor"
))#RFP up dn D2


term=rev(c("axonogenesis","ossification","cell-substrate adhesion",
           "urogenital system development","pattern specification process","muscle cell differentiation",
           "heart morphogenesis","epithelial tube morphogenesis","Wnt signaling pathway",
           "muscle tissue development","regulation of ERK1 and ERK2 cascade",
           "regulation of cellular response to growth factor stimulus","gland development",
           "regulation of endothelial cell apoptotic process"
))#K2,5 dn d2

term=rev(c(
  "renal system development","ossification",
  "pattern specification process",
  "sensory organ morphogenesis","cell-cell signaling by wnt",
  "axonogenesis","synapse organization","cell-substrate adhesion","cell junction assembly",
           "regulation of neurogenesis","gland development",
           "muscle tissue development"
))#K2,5 dn d2 520 395  170 

term=rev(c("pattern specification process","glial cell differentiation",
           "myeloid cell differentiation",
           "axonogenesis","forebrain development",
           "cell-cell adhesion via plasma-membrane adhesion molecules",
           "regulation of developmental growth",
           "urogenital system development",
           "outflow tract morphogenesis","neural tube development",
           "cell junction assembly","skeletal system morphogenesis"
))#K2,5 dn d5 1181 539  159 

term=rev(c("axonogenesis","cell junction assembly","carboxylic acid biosynthetic process",
           "placenta development","epidermis development",
           "skin development","BMP signaling pathway","stem cell differentiation",
           "cellular response to leukemia inhibitory factor",
           "embryonic placenta development","cell fate commitment","muscle organ development",
           "mesenchyme development"
))#K2,5 up d2 637 436 339  

term=rev(c("cellular response to leukemia inhibitory factor","germ cell development",
           "response to oxidative stress","negative regulation of cell cycle","histone modification",
           "epidermis development","stem cell differentiation",
           "morphogenesis of a branching structure","regulation of apoptotic signaling pathway",
           "placenta development","cell fate commitment","gland development","response to BMP","urogenital system development"
           ))#K2,5 up d5 1236 597 689  


Compare_Cluster_GO_results(deglist,
                           top_go_term = 5,#pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = list.files(pattern = "BP_Simplify_go.csv",full.names = T)[c(35,42,34,41)],
                           custom_design_go_term = term)

Compare_Cluster_GO_results(deglist,
                           top_go_term = 5,#pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = list.files(pattern = "BP_Simplify_ver2_go.csv",full.names = T)[c(43,35,32)],
                           custom_design_go_term = term)

xx=read.csv("K2,5_D2_up_go/iKlf5_D2 (775)_BP_Simplify_go.csv")
Plot_GO_result(xx[c(1,3,6,9,11,18,23,45,73,256),],display_top = 10)

xx.mus.go.f=xx[c(1,3,6,9,11,18,23,45,73,256),]
display_top=10
name="GO-TERM"
bar_color="#F27C4B"
xlab="Biological.Process"
study_gene=c("Klf5")
sub_study_gene=c()
force_show_gene_in_term=F
force_show_gene_in_term_number=10
geom_text_size=6
save=F


d2_K5_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _1.5/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d5_K5_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf5 _1.5/Day5_L_iEmpty vs Day5_L_iKlf5_degene_list.txt",header = T)

d2_K5_l_up_1.5=intersect(d2_K5_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K5_l_dn_1.5=intersect(d2_K5_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d5_K5_l_up_1.5=intersect(d5_K5_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))
d5_K5_l_dn_1.5=intersect(d5_K5_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))

d2_K2_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _1.5/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",header = T)
d5_K2_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf2 _1.5/Day5_L_iEmpty vs Day5_L_iKlf2_degene_list.txt",header = T)

d2_K2_l_up_1.5=intersect(d2_K2_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K2_l_dn_1.5=intersect(d2_K2_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d5_K2_l_up_1.5=intersect(d5_K2_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))
d5_K2_l_dn_1.5=intersect(d5_K2_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))

d2_K4_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _1.5/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",header = T)
d5_K4_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf4 _1.5/Day5_L_iEmpty vs Day5_L_iKlf4_degene_list.txt",header = T)

d2_K4_l_up_1.5=intersect(d2_K4_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))
d2_K4_l_dn_1.5=intersect(d2_K4_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))
d5_K4_l_up_1.5=intersect(d5_K4_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))
d5_K4_l_dn_1.5=intersect(d5_K4_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))



k2_d5_vs_d2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/K2_L_D2_vs_K2_L_D5_1.5/K2_L_D2_vs_K2_L_D5_1.5_degene_list.txt",header = T)
k4_d5_vs_d2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/K4_L_D2_vs_K4_L_D5_1.5/K4_L_D2_vs_K4_L_D5_1.5_degene_list.txt",header = T)
k5_d5_vs_d2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/K5_L_D2_vs_K5_L_D5_1.5/K5_L_D2_vs_K5_L_D5_1.5_degene_list.txt",header = T)

k2_d5_vs_d2_up=intersect(k2_d5_vs_d2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(4,5)],1,function(x) max(x)>2),]))
k4_d5_vs_d2_up=intersect(k4_d5_vs_d2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(8,9)],1,function(x) max(x)>2),]))
k5_d5_vs_d2_up=intersect(k5_d5_vs_d2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(12,13)],1,function(x) max(x)>2),]))


k2_d5_vs_d2_up_clean=intersect(d5_K2_l_up_1.5,k2_d5_vs_d2_up)
k4_d5_vs_d2_up_clean=intersect(d5_K4_l_up_1.5,k4_d5_vs_d2_up)
k5_d5_vs_d2_up_clean=intersect(d5_K5_l_up_1.5,k5_d5_vs_d2_up)


using_mat=tpm_total[rownames(tpm_total)%in%k4_d5_vs_d2_up_clean,c(1,19,20,22,23,106)]
colnames(using_mat)=c("EpiSC","ESC-1","ESC-2","TSC-1","TSC-2","PGCLC_D6")


Plot_Box_Violin_plot(using_mat,show_median_absolute_line = T,ylab = "Log2(TPM+0.1)",violin_width = 1,
                     title = "K4",ylim = c(-2,8),
                     # show_p_value = T,pvalue_label_y_position = 5,
                     # compare_list = list(c("PGC_E14_mPGCLC_d6","ESC_DMSO_2iL"))
                     )



d=Cluster_mat_by_Kmeans(using_mat,kmeans = 4)

anno=data.frame(kmean=d$kmean,row.names = rownames(d))

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(d[,1:c(dim(d)[2]-1)],cluster_rows = F,cluster_cols = F,
                   scale = "row",show_rownames = F,annotation_row = anno,
                   breaks = bk,main="d2",
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))


# '-------------------------------- ####
# D2_K5_vs_RFP_DEG have bind ####
setwd("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation")

d2_K5_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _1.5/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d2_K5_l_up_1.5=intersect(d2_K5_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K5_l_dn_1.5=intersect(d2_K5_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))

d2_K5_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d2_K5_l_up_2=intersect(d2_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K5_l_dn_2=intersect(d2_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))


CUTTag_bed=list.files("D:/????/Data processing/#Klf5 Figure#/All_CUTTag_bedfile/q0_sorted/",full.names = T)

K5_L=readPeakFile(CUTTag_bed[13])
K5_L_anno=Peak_Anno(K5_L,up = 20000,down =20000)

K5_L_anno

####1.K5 active have bind####
#####Version1####
K5_L_anno
d2_K5_l_up_2
length(d2_K5_l_up_2)
length(intersect(K5_L_anno,d2_K5_l_up_2))

k5_bind_up=data.frame(condition=factor(c("Bind","No_Bind"),levels = c("No_Bind","Bind")),
              num=c(1252,146))

Plot_Pie(k5_bind_up,"Klf5_Active_and_Binding_Genes",show_legend = F,Only_show_precentage = F,
         color=c("#33AEEE","#F08080"))

#####Version2####
homer_mm10_TSS_2k_bed=read.table("D:/????/Data processing/Data_base/homer_mm10_TSS_2k.bed")


K5_active_mat=bitr(d2_K5_l_dn_2,fromType = "SYMBOL",toType = "REFSEQ",OrgDb = org.Mm.eg.db)
K5_active_matNM=K5_active_mat[grepl("NM",K5_active_mat$REFSEQ),]
K5_active_matNR=K5_active_mat[grepl("NR",K5_active_mat$REFSEQ),]

K5_active_mat=rbind(K5_active_matNM,K5_active_matNR)

K5_active_bed=Btable_insert_into_Atable(K5_active_mat,homer_mm10_TSS_2k_bed,2,5)
K5_active_bed=na.omit(K5_active_bed)
K5_active_bed=K5_active_bed[!duplicated(K5_active_bed$SYMBOL),]


K5_L_anno_bind_bed=K5_active_bed[K5_active_bed$SYMBOL%in%intersect(K5_L_anno,K5_active_bed$SYMBOL),]
K5_L_anno_bind_bed$deepTools_group="Bind"

K5_L_anno_nobind_bed=K5_active_bed[K5_active_bed$SYMBOL%in%setdiff(K5_active_bed$SYMBOL,K5_L_anno),]
K5_L_anno_nobind_bed$deepTools_group="No_Bind"

K5_L_anno_bed_fin=rbind(K5_L_anno_bind_bed,K5_L_anno_nobind_bed)
K5_L_anno_bed_fin=K5_L_anno_bed_fin[,c(3:5,1,8)]
colnames(K5_L_anno_bed_fin)=c("#Chr","Start","End","Genes","deepTools_group")
dim(K5_L_anno_bed_fin)
head(K5_L_anno_bed_fin)

Output_bedfile(K5_L_anno_bed_fin,name="K5_L_bind_active_anno.bed",colname = T)


Output_bedfile(K5_active_bed[,c(3:5,1)],name="K5_active.bed",colname = F)


K5_active=readPeakFile("K5_active.bed")

CUTTag_bed=list.files("D:/????/Data processing/#Klf5 Figure#/All_CUTTag_bedfile/q0_sorted/",full.names = T)

K5_L=readPeakFile(CUTTag_bed[13])

ol=findOverlapsOfPeaks(K5_active,K5_L,maxgap = 200)

K5_active_have_bind=ol$overlappingPeaks$`K5_active///K5_L`
K5_active_have_bind=K5_active_have_bind[!duplicated(K5_active_have_bind$V4),]

uniquePeaks=rownames(as.data.frame(ol$uniquePeaks))
K5_active_no_bind=K5_active[as.numeric(strsplit2(uniquePeaks[grepl("K5_active__",uniquePeaks)],split = "\\__")[,2]),]


K5_active_have_bind=as.data.frame(K5_active_have_bind[,c(2:4,7)])
K5_active_have_bind$deepTools_group="Bind"

K5_active_no_bind=as.data.frame(K5_active_no_bind)
K5_active_no_bind=K5_active_no_bind[,c(1:3,6)]
K5_active_no_bind$deepTools_group="No_Bind"

K5_L_anno_bed_fin2=rbind(K5_active_have_bind,K5_active_no_bind)
colnames(K5_L_anno_bed_fin2)=c("#Chr","Start","End","Genes","deepTools_group")
dim(K5_L_anno_bed_fin2)
head(K5_L_anno_bed_fin2)

Output_bedfile(K5_L_anno_bed_fin2,name="K5_active2.bed",colname = T)


#####Heatmap####
tpm_total_rep=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new.csv",row.names = 1)
colnames(tpm_total_rep)
# d2_K5_l_2_mat=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_lfcmatrix.txt",header = T)
# 
# fold_mat=Btable_insert_into_Atable(K5_L_anno_bed_fin2,d2_K5_l_1.5_mat,4)
# fold_mat=fold_mat[,c(4,7)]
# fold_mat$log2FoldChange=as.numeric(fold_mat$log2FoldChange)*c(-1)
# rownames(fold_mat)=fold_mat$Genes
# 
# anno=Btable_insert_into_Atable(fold_mat,log2(tpm_total[,c(16,12)]))
# anno=data.frame(`Log2(TPM)`=as.numeric(anno$PNT_K5_L_D2),row.names = anno$Genes)
# fold_mat=data.frame(`Log2(FoldChange)`=as.numeric(fold_mat$log2FoldChange),row.names = fold_mat$Genes)

K5_L_anno_bed_fin=read.table('K5_D5_L_bind_active_anno.bed',header = F)

k5_active_tpm_mat=Btable_insert_into_Atable(K5_L_anno_bed_fin,tpm_total_rep[,c(34,35,26,27)],4)#D2
k5_active_tpm_mat=Btable_insert_into_Atable(K5_L_anno_bed_fin,tpm_total_rep[,c(36,37,28,29)],4)#D5

k5_active_tpm_mat=as.data.frame(apply(k5_active_tpm_mat[,6:9], 2, as.numeric))
rownames(k5_active_tpm_mat)=K5_L_anno_bed_fin$Genes
head(k5_active_tpm_mat)

k5_active_tpm_mat=rbind(tpm_total_rep[rownames(tpm_total_rep)%in%intersect(d2_K5_l_dn_2,K5_L_anno),c(34,35,26,27)],
                        tpm_total_rep[rownames(tpm_total_rep)%in%setdiff(d2_K5_l_dn_2,K5_L_anno),c(34,35,26,27)])


bk=seq(-1.5,1.5,by=0.01)
p1=pheatmap::pheatmap(k5_active_tpm_mat,show_rownames = F,cellwidth = 25,cellheight = 0.2,main = paste0("iKlf5_D5_Active_Genes(n=",dim(K5_L_anno_bed_fin)[1],")"),
                   cluster_cols = F,cluster_rows = F,breaks = bk,scale = "row",gaps_row = length(intersect(d2_K5_l_dn_2,K5_L_anno)),
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))
bk=seq(0,10,by=1)
p2=pheatmap::pheatmap(anno,gaps_row = dim(K5_active_have_bind)[1],cellwidth = 25,cellheight = 0.2,
                   cluster_cols = F,cluster_rows = F,breaks = bk,show_rownames = F,
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0),
                              colorRampPalette(colors = brewer.pal(n=5,name = "Greens"))((length(bk))*1)))

write.csv(cbind(k5_active_tpm_mat,anno),"iK5_D2_active_gene_heatmap_1389.csv")

#####GO####


k5_bind_ac=read.table("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/K5_L_bind_active_anno.bed")
k5_bind_ac=k5_bind_ac[k5_bind_ac$V5=="Bind","V4"]

d25_bind_up=intersect(d5_K5_l_up_1.5,k5_bind_ac)

d2_K5_l_up_1.5
d25_up=intersect(d2_K5_l_up_1.5,d5_K5_l_up_1.5)

xx=GO_function(d25_up,simplify_ver2 = T,name = "d25_ac")

xx=GO_function(k5_bind_ac,simplify_ver2 = T,name = "d2_bind_ac")

xx=GO_function(d25_bind_up,simplify_ver2 = T,name = "d25_bind_ac")

d2_K5_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d5_K5_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf5 _2/Day5_L_iEmpty vs Day5_L_iKlf5_degene_list.txt",header = T)

d2_K5_l_up_2=intersect(d2_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d5_K5_l_up_2=intersect(d5_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))

ETP2=read.csv("ETP_gene_WGCNA_ver4_4240_.csv")

ETP2list=Convert_Long_dataframe_to_list(ETP2,group_index = 6,variable_index = 1)

Study_Gene_Distribution_in_gene_list(d2_K5_l_up_2,gene_list =ETP2list)


d2_rpf_v_k5_l_up_2=intersect(d2_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d5_rpf_v_k5_l_up_2=intersect(d5_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))


d2_rfp_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Epi2.1_FA_X_vs_PNT_RFP_L_D2_2/Epi2.1_FA_X_vs_PNT_RFP_L_D2_2_degene_list.txt",header = T)
d5_rfp_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Epi2.1_FA_X_vs_PNT_RFP_L_D5_2/Epi2.1_FA_X_vs_PNT_RFP_L_D5_2_degene_list.txt",header = T)

d2_rfp_v_epi_l_up_2=intersect(d2_rfp_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,1)],1,function(x) max(x)>2),]))
d5_rfp_v_epi_l_up_2=intersect(d5_rfp_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,1)],1,function(x) max(x)>2),]))
d2_rfp_v_epi_l_dn_2=intersect(d2_rfp_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,1)],1,function(x) max(x)>2),]))
d5_rfp_v_epi_l_dn_2=intersect(d5_rfp_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,1)],1,function(x) max(x)>2),]))



#####RFP vs Epi D2,5 bubble plot####
primed=c("Pou5f1","Nanog","Nodal","Fgf5","Otx2","Zic2","Fgf8")
using_mat=tpm_total[rownames(tpm_total)%in%primed,c(1,16,17)]
tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,1]+0.1))
Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_fold,#by_step = 0.75,
                            study_gene = primed,bubble_scale_by_new_function = F,limit_exp = c(0,9),
                            limit_fold = c(-2,2),#column_split_mat = column_split_mat,
                            title = "RFP_PNT_Day25",z_score =T)


rfp_v_epi_d25_deglist=list(d2_rfp_up=d2_rfp_v_epi_l_up_2,
                           d5_rfp_up=d5_rfp_v_epi_l_up_2)
rfp_v_epi_d25_deglist=list(d2_rfp_dn=d2_rfp_v_epi_l_dn_2,
                           d5_rfp_dn=d5_rfp_v_epi_l_dn_2)

Compare_Cluster_GO_results(rfp_v_epi_d25_deglist,top_go_term = 10,simplify_ver2 = T)

deglist=list(d2_K5_l_up_2=intersect(d2_K5_l_up_2,K5_L_anno),
             d5_K5_l_up_2=intersect(d5_K5_l_up_2,K5_L_anno))
VennPlot(deglist)


Venner_plot(list(d2_K5_l_up_2=d2_K5_l_up_2,
                 d5_K5_l_up_2=d5_K5_l_up_2,
                 d2_rfp_v_epi_l_up_2=d2_rfp_v_epi_l_up_2,
                 d5_rfp_v_epi_l_up_2=d5_rfp_v_epi_l_up_2))
            
Venner_plot(list(d5_rpf_v_k5_l_up_2=d5_rpf_v_k5_l_up_2,
                 d5_rfp_v_epi_l_up_2=d5_rfp_v_epi_l_up_2,
                 d5_rpf_v_k5_l_dn_2=d5_K5_l_up_2))

path="D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/K5_D2,5_up_go/"
Compare_Cluster_GO_results(deglist,top_go_term = 10,simplify_ver2 = T)



list=list(D2_K5_LIF_up=d2_K5_l_up_2,
     D5_K5_LIF_up=d5_K5_l_up_2)
degene_table=as.data.frame(sapply(list, "[", i = 1:max(sapply(list, length))))
write.csv(degene_table,"K5_D2,5_up_gene.csv")

term=rev(c("placenta development","epithelial cell development","mesenchyme development",
           "pattern specification process","response to BMP","muscle organ development",
  "regulation of fat cell differentiation","regulation of stem cell proliferation",
  "somatic stem cell population maintenance",
  "positive regulation of meiotic nuclear division",
  "gland development","cell-substrate adhesion","stem cell differentiation","neuron fate commitment",
  "regulation of neurogenesis"
  ))

term=rev(c("placenta development","epithelial cell development","pattern specification process",
           "mesenchyme development","response to BMP","muscle organ development","regulation of hemopoiesis",
           "regulation of stem cell proliferation","response to leukemia inhibitory factor",
           "positive regulation of meiotic nuclear division",
           "gland development","neuron fate commitment","urogenital system development","endothelium development"
           #"axonogenesis","urogenital system development","regulation of vasculature development",
           #"regulation of neurogenesis","Wnt signaling pathway",
           #"forebrain development","regulation of developmental growth","cell-cell signaling by wnt"
))


Compare_Cluster_GO_results(list(d2_K5_l_up_2=intersect(d2_K5_l_up_2,K5_L_anno),
                                d5_K5_l_up_2=intersect(d5_K5_l_up_2,K5_L_anno)),
                                #d2_rfp_v_epi_l_up_2=d2_rfp_v_epi_l_up_2,d5_rfp_v_epi_l_up_2=d5_rfp_v_epi_l_up_2),
                                top_go_term = 10,#pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = list.files(path = "K5_D2,5_up_go",
                                                                          pattern = "_BP_Simplify_go.csv",full.names = T)[c(1,4)],
                           custom_design_go_term = term)



kegg=KEGG_function(tpm_total_rep,k5_bind_ac)
Plot_KEGG_result(kegg)



#PAGA DEG
setwd("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/PAGA_DEG_GO_heatmap/")

deg_res=readRDS("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/PAGA_DEG_GO_heatmap/PAGA_big_group_DEG_resluts_minpct0.05.rds")
deg_res_filter=deg_res[deg_res$p_val_adj<0.01,]
deg_res_filter=deg_res_filter[deg_res_filter$avg_log2FC>1,]
table(deg_res_filter$cluster)

deg_list=Convert_Long_dataframe_to_list(deg_res_filter,group_index = 6,variable_index = 7)

term=rev(c("cellular response to leukemia inhibitory factor",
           "germ cell development",
           "DNA methylation or demethylation",
           "stem cell population maintenance",
           "histone modification",
           "regulation of neuron differentiation",
           "muscle cell differentiation",
           "response to BMP",
           "regulation of epithelial cell differentiation",
           "endothelial cell differentiation",
           "epidermis development",
           "placenta development",
           "regulation of vasculature development",
           "cell-cell junction organization",
           "blood vessel diameter maintenance"))

Compare_Cluster_GO_results(deg_list[c(3,1,2)],custom_design_go_term = term,pvalue_max_limit = 30,
                           input_go_cluster_resluts_pathlist = list.files(pattern = "BP_Simplify_ver2_go.csv")[c(2,1,3)],
                           simplify_ver2 = T,top_go_term = 5)



####2.ESC_TSC_PGC_gene K5 active bind####
esc=c("Pou5f1","Zfp42","Nanog","Dppa5a","Tfcp2l1","Nr5a2","Tbx3","Nanog","Klf2","Klf4")
tsc=c("Cdx2","Tfap2c","Elf5","Gata3","Hand1","Eomes","Tead4","Fgfr2","Krt7",
      "Ascl2","Plet1","Rin3","Tinagl1")
pgc=c("Prdm14","Prdm1","Dppa3","Nanos3","Kit","Tfap2c","Ddx4","Rhox6","Rhox9",
      "Dazl","Msx1","Prdm1","Tfap2c","Msx2","Dlx5","Id3","Sox15")

condition=c("iKLF5_LIF_D2","iKLF5_LIF_D2","iEmpty_LIF_D2","iEmpty_LIF_D2")
Plot_comparation_of_two_sample_transcriptomic_profiles(total_counts_mat[rownames(total_counts_mat)%in%NEpi_up,c(26,27,34,35)],
                                                       condition_1 ="iKLF5_LIF_D2",condition_2 ="iEmpty_LIF_D2",condition_list =condition,hide_nosig_gene = F,
                                                       study_genes = esc,DESeq2_normalized = T,title = "ESC Genes(n=1111)")
Plot_comparation_of_two_sample_transcriptomic_profiles(total_counts_mat[rownames(total_counts_mat)%in%TEpi_up,c(26,27,34,35)],
                                                       condition_1 ="iKLF5_LIF_D2",condition_2 ="iEmpty_LIF_D2",condition_list =condition,hide_nosig_gene = F,
                                                       study_genes = tsc,DESeq2_normalized = T,title = "TSC Genes(n=979)")
Plot_comparation_of_two_sample_transcriptomic_profiles(total_counts_mat[rownames(total_counts_mat)%in%PEpi_up,c(26,27,34,35)],
                                                       condition_1 ="iKLF5_LIF_D2",condition_2 ="iEmpty_LIF_D2",condition_list =condition,hide_nosig_gene = F,
                                                       study_genes = tsc,DESeq2_normalized = T,title = "PGC Genes(n=1607)")



####3.GO####

xx=GO_function(intersect(d2_K5_l_up_1.5,K5_L_anno),simplify_ver2 = T)
xx=GO_function(K5_active_have_bind$V4,simplify_ver2 = T)
xxx=as.data.frame(xx)
term_select=c("placenta development","epithelial cell development","pattern specification process",
              "response to BMP","maintenance of cell number","trophectodermal cell differentiation",
              "cellular response to leukemia inhibitory factor",
              "regulation of DNA-binding transcription factor activity","germ cell migration",
              "axonogenesis","cell-substrate adhesion","muscle cell development",
              "heart morphogenesis","regulation of endothelial cell apoptotic process",
              "response to transforming growth factor beta","Wnt signaling pathway",
              "fat cell differentiation","mesenchymal cell differentiation")
xxx=xxx[c(1,5,10,28,119,151,431),]
xxxx=xxx[xxx$Description%in%term_select,]
Plot_GO_result(xxxx,display_top = 6,bar_color = "#F08080",force_show_gene_in_term = T,force_show_gene_in_term_number = 6)

Compare_Cluster_GO_results(list(UP=intersect(d2_K5_l_up_2,K5_L_anno),
                                DN=intersect(d2_K5_l_dn_2,K5_L_anno)))

Compare_Cluster_GO_results(input_go_cluster_resluts_pathlist = list.files(path = "K5_D2,5_up_go",pattern = "y_go.csv",full.names = T)[c(9,8)],
                           list(UP=intersect(d2_K5_l_up_2,K5_L_anno),
                                DN=intersect(d2_K5_l_dn_2,K5_L_anno)),pvalue_max_limit = 25,
                            term_select_by_study_genes = F,custom_add_go_term  = rev(term_select),
                           top_go_term = 0,simplify_ver2 = T)




####4.K5 activate ETP gene fiter by tf in vivo regulon####

ETP=read.csv("ETP_gene_WGCNA_ver3_4431_.csv")
ETP=ETP[ETP$Cluster!="Primed",]
#table(ETP$Cluster)
k5_ac_bind=read.table("K5_L_bind_active_anno.bed")
#head(k5_ac_bind)
k5_expr_mat=read.csv("iK5_active_gene_heatmap_1389_have_k5.csv")
#head(k5_expr_mat)
k5_ac_ETP=ETP[ETP$X%in%k5_ac_bind[k5_ac_bind$V5=="Bind",4],]
k5_ac_ETP=Btable_insert_into_Atable(k5_ac_ETP,k5_expr_mat,index_B = 1,index_A = 1)
dim(k5_ac_ETP)
#head(k5_ac_ETP)
tf=read.table("D:/????/Data processing/Data_base/mm_mgi_tfs.txt")#form SCENIC


k5_ac_ETP_tf=k5_ac_ETP[k5_ac_ETP$X%in%tf$V1,]
#write.csv(k5_ac_ETP_tf,"iK5_Day2_activate_TFs.csv")


filter_regulon_mat=read.csv("ICM_TE2_sorted_regulon.csv")
filter_regulon_mat2=read.csv("PGC_sorted_regulon.csv")
invivo_regulon=c(filter_regulon_mat$E3.5_ICM[1:30],
                 filter_regulon_mat$E3.5_TE[1:30],
                 filter_regulon_mat2$E9.5_PGC[1:30])

k5_ac_ETP_tf_invivo=k5_ac_ETP_tf[k5_ac_ETP_tf$X%in%invivo_regulon,]
#write.csv(k5_ac_ETP_tf_invivo,"iK5_Day2_activate_TFs_invivo.csv")


# '-------------------------------- ####
# ESC TSC PGCLC vs EpiSC ####


tpm_mat_mean=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_mat_mean)
using_mat=tpm_mat_mean[,c(1,116,
  19,20,22,23)]

esc_vs_tsc_up=read.table("Naive_vs_TSC_2/Naive_vs_TSC_2_degene_list.txt",header = T)
esc_vs_tsc_up=na.omit(esc_vs_tsc_up$Up_regulate)
esc_vs_tsc_up=intersect(esc_vs_tsc_up,rownames(using_mat[apply(using_mat[,c(3,4)],1,function(x) mean(x) >10),]))

pid_vs_tsc_up=read.table("TSC_vs_Primed_2/TSC_vs_Primed_2_degene_list.txt",header = T)
pid_vs_tsc_up=na.omit(pid_vs_tsc_up$Down_regulate)
pid_vs_tsc_up=intersect(pid_vs_tsc_up,rownames(using_mat[apply(using_mat[,c(1:2)],1,function(x) mean(x) >10),]))

intersect(esc_vs_tsc_up,pid_vs_tsc_up)











tpm_total=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_total)


use_counts=tpm_total[,c(1,46,53,54,56,57,113:115)]

total_counts_mat=read.csv("D:/????/Data processing/20230820_Klf5_figure_all_RNAseq_count_new.csv",row.names = 1,header = T)
colnames(total_counts_mat)

use_counts=total_counts_mat[,c(40:43,81:83,
                               48:49,54:55,
                               90:92,104:105)]
colnames(use_counts)
head(use_counts)
condition=c(rep("Naive",7),rep("TSC",4),rep("PGCLC",5))

k=levels(factor(condition))##??????Ҫ????ȡ??��?е?Ԫ??????
com=t(combn(k,2))
for(n in 1:length(com[,1])){
  c1=com[n,1]
  c2=com[n,2]
Auto_DEseq2(use_counts,DESeq2_normalization = T,
            condition = condition,
            celltype = condition,
            condition1 = c1,
            condition2 = c2,
            path = paste(c1,"vs",c2,"p.adj0.1",sep = "_"),p.adj = 0.1,
            lfc = log2(4))
}


NEpi=read.table("Naive_vs_Primed_2/Naive_vs_Primed_2_degene_list.txt",header = T)
TEpi=read.table("TSC_vs_Primed_2/TSC_vs_Primed_2_degene_list.txt",header = T)
PEpi=read.table("PGCLC_vs_Primed_2/PGCLC_vs_Primed_2_degene_list.txt",header = T)

NEpi_up=intersect(NEpi$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(19,20,1)],1,function(x) max(x)>10),]))
NEpi_dn=intersect(NEpi$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(19,20,1)],1,function(x) max(x)>10),]))

TEpi_up=intersect(TEpi$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(22,23,1)],1,function(x) max(x)>10),]))
TEpi_dn=intersect(TEpi$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(22,23,1)],1,function(x) max(x)>10),]))

PEpi_up=intersect(PEpi$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(106,117,1)],1,function(x) max(x)>10),]))
PEpi_dn=intersect(PEpi$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(106,117,1)],1,function(x) max(x)>10),]))

Epiu=intersect(NEpi_dn,intersect(TEpi_dn,PEpi_dn))
ETP=intersect(NEpi_up,intersect(TEpi_up,PEpi_up))





NP=read.table("Naive_vs_PGCLC_2/Naive_vs_PGCLC_2_degene_list.txt",header = T)
NT=read.table("Naive_vs_TSC_2/Naive_vs_TSC_2_degene_list.txt",header = T)
PT=read.table("PGCLC_vs_TSC_2/PGCLC_vs_TSC_2_degene_list.txt",header = T)

NP_up=intersect(NP$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(19,20,103,106,117)],1,function(x) max(x)>10),]))
NP_dn=intersect(NP$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(19,20,103,106,117)],1,function(x) max(x)>10),]))
NP_no=intersect(NP$No_sig,rownames(tpm_total[apply(tpm_total[,c(19,20,103,106,117)],1,function(x) max(x)>10),]))

NT_up=intersect(NT$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(19,20,103,22,23)],1,function(x) max(x)>10),]))
NT_dn=intersect(NT$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(19,20,103,22,23)],1,function(x) max(x)>10),]))
NT_no=intersect(NT$No_sig,rownames(tpm_total[apply(tpm_total[,c(19,20,103,22,23)],1,function(x) max(x)>10),]))


PT_up=intersect(PT$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(106,117,22,23)],1,function(x) max(x)>10),]))
PT_dn=intersect(PT$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(106,117,22,23)],1,function(x) max(x)>10),]))
PT_no=intersect(PT$No_sig,rownames(tpm_total[apply(tpm_total[,c(106,117,22,23)],1,function(x) max(x)>10),]))


c(d2_K2_l_up_1.5,d2_K4_l_up_1.5,d2_K5_l_up_1.5)
using_mat=tpm_total[,c(16,4,8,12,19,20,103,22,23,106,117)]
dim(using_mat)
colnames(using_mat)

plot_mat_sub=plot_mat[plot_mat$group=="TP",]
d=Cluster_mat_by_Kmeans(plot_mat_sub[rownames(plot_mat_sub)%in%c(d2_K5_l_up_1.5,d2_K2_l_up_1.5),1:10],
                        kmeans = 3,use_col = c(1:3))


anno=data.frame(kmean=d$kmean,row.names = rownames(d))

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(d[,1:c(dim(d)[2]-1)],cluster_rows = F,cluster_cols = F,
                   scale = "row", 
                   breaks = bk,
                   show_rownames = F,annotation_row = anno,
                   main="d2",
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))
rownames(d[d$kmean==3,])


d$gene=rownames(d)
test=Convert_Long_dataframe_to_list(d,4,5)

Tu=intersect(PT_dn,NT_dn)
Nu=intersect(NP_up,NT_up)
Pu=intersect(PT_up,NP_dn)
TP=intersect(intersect(NT_dn,NP_dn),PT_no)
NP=intersect(intersect(NT_up,PT_up),NP_no)
NT=intersect(intersect(NP_up,PT_dn),NT_no)


plotlist=list(Tu=Tu,Nu=Nu,Pu=Pu,NT=NT,NP=NP,TP=TP)
length(unlist(plotlist))

Study_Gene_Distribution_in_gene_list(study_gene,plotlist)


#write.csv(plot_mat,"ETP_genelist.csv")
plotlist_mat=read.csv("ETP_genelist.csv",row.names = 1)
plotlist_mat$gene=rownames(plotlist_mat)
plotlist=Convert_Long_dataframe_to_list(plotlist_mat,group_index = 11,variable_index = 12)
plotlist_mat$gene=NULL

plot_mat=plotlist_mat
plot_mat$group=factor(plot_mat$group,levels = c("Pu","Eu","Tu","ET","EP","PT"))
plot_mat=plot_mat[order(plot_mat$group),]
anno=data.frame(group=plot_mat$group,row.names = rownames(plot_mat))
bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(plot_mat[,c(4:5,9:10,7:8)],cluster_rows = F,cluster_cols = F,
                   scale = "row", 
                   breaks = bk,
                   show_rownames = F,annotation_row = anno,
                   main="d2",
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))
Plot_clustered_heatmap_R4.2.2(plot_mat[,c(4:5,9:10,7:8)],
                              colum_split_fun = c(1,1,2,2,3,3),
                              kmeans = 7,
                              name="ETP_genelist",
                              clustering_method="complete",
                              clustered_and_order =plot_mat$group,rename_cluster = F,
                              scale = "row",
                              set_seed = 169,
                              study_genes=study_gene,
                              dont_clustered = T,
                              show_study_genes = T,
                              show_row_names = F,
                              breaks=c(-2,2),#Default
                              height = 8,#Default
                              width = 4,#Default
                              anno_font_size=8,#use_raster = T,raster_quality = 25
)




using_mat_scale=t(scale(t(using_mat)))
using_mat_scale=using_mat
using_mat_scale=data.frame(ESC=as.numeric(apply(using_mat_scale[,1:3], 1, mean)),
                           TSC=as.numeric(apply(using_mat_scale[,4:5], 1, mean)),
                           PGC=as.numeric(apply(using_mat_scale[,6:7], 1, mean)),
                           row.names = rownames(using_mat))
using_mat_scale=using_mat_scale[1:200,]
p=ggtern(data = using_mat_scale,aes(x = ESC,y = TSC,z = PGC))+geom
  stat_density_tern(geom="polygon",aes(fill=..level..,alpha=..level..),bdl.val=NA)



gene.pca <- FactoMineR::PCA(t(using_mat_scale), ncp = 10, scale.unit = T, graph = FALSE)
dim(gene.pca$var$contrib)




#overlap with K2,4,5 d2,5 up gene ####
d2_K5_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _1.5/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d5_K5_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf5 _1.5/Day5_L_iEmpty vs Day5_L_iKlf5_degene_list.txt",header = T)

d2_K5_l_up_1.5=intersect(d2_K5_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K5_l_dn_1.5=intersect(d2_K5_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d5_K5_l_up_1.5=intersect(d5_K5_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))
d5_K5_l_dn_1.5=intersect(d5_K5_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))

d2_K2_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _1.5/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",header = T)
d5_K2_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf2 _1.5/Day5_L_iEmpty vs Day5_L_iKlf2_degene_list.txt",header = T)

d2_K2_l_up_1.5=intersect(d2_K2_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K2_l_dn_1.5=intersect(d2_K2_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d5_K2_l_up_1.5=intersect(d5_K2_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))
d5_K2_l_dn_1.5=intersect(d5_K2_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))

d2_K4_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _1.5/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",header = T)
d5_K4_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf4 _1.5/Day5_L_iEmpty vs Day5_L_iKlf4_degene_list.txt",header = T)

d2_K4_l_up_1.5=intersect(d2_K4_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))
d2_K4_l_dn_1.5=intersect(d2_K4_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))
d5_K4_l_up_1.5=intersect(d5_K4_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))
d5_K4_l_dn_1.5=intersect(d5_K4_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))



d2_K5_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)
d5_K5_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf5 _2/Day5_L_iEmpty vs Day5_L_iKlf5_degene_list.txt",header = T)

d2_K5_l_up_2=intersect(d2_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K5_l_dn_2=intersect(d2_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d5_K5_l_up_2=intersect(d5_K5_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))
d5_K5_l_dn_2=intersect(d5_K5_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,13)],1,function(x) max(x)>2),]))

d2_K2_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _2/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",header = T)
d5_K2_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf2 _2/Day5_L_iEmpty vs Day5_L_iKlf2_degene_list.txt",header = T)

d2_K2_l_up_2=intersect(d2_K2_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K2_l_dn_2=intersect(d2_K2_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d5_K2_l_up_2=intersect(d5_K2_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))
d5_K2_l_dn_2=intersect(d5_K2_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,5)],1,function(x) max(x)>2),]))

d2_K4_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _2/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",header = T)
d5_K4_l_2=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day5_L_iEmpty _vs_ Day5_L_iKlf4 _2/Day5_L_iEmpty vs Day5_L_iKlf4_degene_list.txt",header = T)

d2_K4_l_up_2=intersect(d2_K4_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))
d2_K4_l_dn_2=intersect(d2_K4_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))
d5_K4_l_up_2=intersect(d5_K4_l_2$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))
d5_K4_l_dn_2=intersect(d5_K4_l_2$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(17,9)],1,function(x) max(x)>2),]))


####1.radarchart####

using_mat=tpm_total[,c(16,4,8,12,19,20,103,22,23,106,117)]
for(g in 1:length(plotlist)){
  using_mat[rownames(using_mat)%in%plotlist[[g]],"group"]=names(plotlist)[g]
}
plot_mat=na.omit(using_mat)


library(fmsb)
df=plot_mat[rownames(plot_mat)%in%d2_K5_l_up_1.5,]
#df=plot_mat
df=merge_mat
df$group=factor(df$group,levels = c("Eu","EP","Pu","PT","Tu","ET"))
group=c()
K5=c()
for(g in levels(factor(df$group))){
  #K5=c(K5,median(log2(df[df$group%in%g,"PNT_K5_L_D2"]+0.1)))
  K5=c(K5,-mean(df[df$group%in%g,4]))
  group=c(group,g)
}

df2=plot_mat[rownames(plot_mat)%in%d2_K2_l_up_1.5,]
#df2=plot_mat
df2=merge_mat
df2$group=factor(df2$group,levels = c("Eu","EP","Pu","PT","Tu","ET"))
K2=c()
for(g in levels(factor(df$group))){
  #K2=c(K2,median(log2(df[df$group%in%g,"PNT_K2_L_D2"]+0.1)))
  K2=c(K2,-mean(df2[df2$group%in%g,10]))
}

#df3=plot_mat
df3=merge_mat
df3=plot_mat[rownames(plot_mat)%in%d2_K4_l_up_1.5,]
df3$group=factor(df3$group,levels = c("Eu","EP","Pu","PT","Tu","ET"))
K4=c()
for(g in levels(factor(df3$group))){
  #K4=c(K4,median(log2(df3[df3$group%in%g,"PNT_K4_L_D2"]+0.1)))
  K4=c(K4,-mean(df3[df3$group%in%g,16]))
}

df0=plot_mat
df0$group=factor(df0$group,levels = c("Eu","EP","Pu","PT","Tu","ET"))
RFP=c()
for(g in levels(factor(df0$group))){
  RFP=c(RFP,median(log2(df0[df0$group%in%g,"PNT_RFP_L_D2"]+0.1)))
}

#Precentage
p_df=data.frame(t(data.frame(max=rep(0.35,6),min=rep(0.05,6),
                             K2=as.numeric(table(df2$group)/table(plot_mat$group)),
                             #K4=as.numeric(table(df3$group)/table(plot_mat$group)),
                             K5=as.numeric(table(df$group)/table(plot_mat$group))
                             )))
p_df=p_df[,c(3,1,5,4,6,2)]*100
#Number
p_df=data.frame(t(data.frame(max=rep(212,6),min=rep(0,6),
                             K2=as.numeric(table(df2$group)),
                             #K4=as.numeric(table(df3$group)),
                             K5=as.numeric(table(df$group)))))
p_df=p_df[,c(3,1,5,4,6,2)]
#mean&median
p_df=data.frame(t(data.frame(max=rep(5,6),min=rep(1,6),
                             RFP=RFP,
                             K2=K2,
                             #K4=K4,
                             K5=K5)))

#LFC
p_df=data.frame(t(data.frame(max=rep(0.1,6),min=rep(-0.1,6),
                             K2=K2,
                             K4=K4,
                             K5=K5)))




colnames(p_df)=group
color1=c("black","#00AFBB","#E7B800","red")
color2=c("black","#00AFBB","red")
color3=c("#00AFBB","#E7B800","red")
color4=c("#00AFBB","red")

color=color4
radarchart(p_df,
           axistype = 1,
           # Customize the polygon
           pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
           # Customize the grid
           cglcol = "grey", cglty = 1, cglwd = 0.8,
           # Customize the axis
           axislabcol = "black", 
           caxislabels = c("7%", "14%", "21%", "28%","35%"),
           #caxislabels = c("0%", "7%", "14%", "21%","28%"),
           #caxislabels = c(0,70, 140, 210,280),
           #caxislabels = c(0,53, 106, 159,212),
           #caxislabels = c(1,2,3,4,5),
           # Variable labels
           vlcex = 0.7, vlabels = colnames(p_df))
legend(
  x = "bottom", legend = rownames(p_df[-c(1,2),]), horiz = T,
  bty = "n", pch = 20 , col =color,
  text.col = "black", cex = 1, pt.cex = 1.5
)



xx=as.data.frame(GO_function(d2_K5_l_up_1.5))


####2.scatterplot####

d2_K5_l_2_mat=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _2/Day2_L_iEmpty vs Day2_L_iKlf5_degene_lfcmatrix.txt",header = T)
d2_K2_l_2_mat=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _2/Day2_L_iEmpty vs Day2_L_iKlf2_degene_lfcmatrix.txt",header = T)
d2_K4_l_2_mat=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _2/Day2_L_iEmpty vs Day2_L_iKlf4_degene_lfcmatrix.txt",header = T)

d2_K5_l_2_mat=d2_K5_l_2_mat[rownames(d2_K5_l_2_mat)%in%rownames(plotlist_mat),]
d2_K2_l_2_mat=d2_K2_l_2_mat[rownames(d2_K2_l_2_mat)%in%rownames(plotlist_mat),]
d2_K4_l_2_mat=d2_K4_l_2_mat[rownames(d2_K4_l_2_mat)%in%rownames(plotlist_mat),]

merge_mat=plotlist_mat[,c(1,11)]
merge_mat=Btable_insert_into_Atable(merge_mat,d2_K5_l_2_mat)
merge_mat=Btable_insert_into_Atable(merge_mat,d2_K2_l_2_mat)
merge_mat=Btable_insert_into_Atable(merge_mat,d2_K4_l_2_mat)
merge_mat=as.data.frame(merge_mat)
merge_mat[is.na(merge_mat[,4]),4]=0
merge_mat[is.na(merge_mat[,10]),10]=0
merge_mat[is.na(merge_mat[,8]),8]=1
merge_mat[is.na(merge_mat[,14]),14]=1
merge_mat[is.na(merge_mat[,16]),16]=0
merge_mat[is.na(merge_mat[,20]),20]=1
merge_mat[,c(4,10,8,14,16,20)]=apply(merge_mat[,c(4,10,8,14,16,20)], 1, as.numeric)
tail(merge_mat)
boxplot(merge_mat[merge_mat$group=="PT",c(4,10,16)])

plot_mat=data.frame(gene=merge_mat$X,
                    K5_active_lfc=-as.numeric(merge_mat[,4]),
                    K2_active_lfc=-as.numeric(merge_mat[,10]),
                    padj=apply(merge_mat[,c(8,14)], 1, min),
                    top_gene=-as.numeric(merge_mat[,4])-c(-as.numeric(merge_mat[,10])),
                    group=merge_mat[,2])
for(g in 1:length(plotlist)){
  plot_mat[plot_mat$gene%in%plotlist[[g]],"group"]=names(plotlist)[g]
}
plot_mat=plot_mat[order(plot_mat$top_gene),]
plot_mat$padj_sig="No.Sig"
plot_mat[plot_mat$padj<0.1,"padj_sig"]="P.adj<0.1"

dim(plot_mat)
top_gene=c(plot_mat[1:10,1],plot_mat[3308:3318,1])


scaterrplot=
  function(plot_mat,study_gene,gene_type,color,title){
    
    plot_mat[plot_mat[,2]>12,2]=12
    plot_mat[plot_mat[,2]<c(-4),2]=c(-4)
    plot_mat[plot_mat[,3]>12,3]=12
    plot_mat[plot_mat[,2]<c(-4),3]=c(-4)
    
plot_mat_s=plot_mat[plot_mat$group%in%c(gene_type),]

marker_genes_mat=plot_mat_s[plot_mat_s$gene%in%c(study_gene),]
#marker_genes_mat=marker_genes_mat[abs(marker_genes_mat$K5_active_lfc-marker_genes_mat$K2_active_lfc)>1,]

plot_mat_sp=plot_mat_s[plot_mat_s$padj_sig=="P.adj<0.1",]
k2up=plot_mat_sp[plot_mat_sp$top_gene<c(-1),]
k2up=k2up[k2up$K2_active_lfc>1,"gene"]
k5up=plot_mat_sp[plot_mat_sp$top_gene>1,]
k5up=k5up[k5up$K5_active_lfc>1,"gene"]


p=ggplot(plot_mat_s, aes(x=K5_active_lfc,y=K2_active_lfc,color=padj_sig)) +
  geom_point(alpha=0.5,shape=16,size=2.1) +
  #show_col(brewer.pal(n = 6,name = "Set2"))
  annotate("segment",x=12,y=11,xend=1,yend = 0, linetype="dotted",color="#666666",size=0.8)+
  annotate("segment",x=11,y=12,xend=0,yend = 1, linetype="dotted",color="#666666",size=0.8)+
  annotate("segment",x=1,y=0,xend=1,yend = -4, linetype="dotted",color="#666666",size=0.8)+
  annotate("segment",x=0,y=1,xend=-4,yend = 1, linetype="dotted",color="#666666",size=0.8)+
  geom_point(size = 1.5, shape = 1, data = marker_genes_mat,color="black",stroke=1.7) +#marker?????Ӹ???Ȧ
  geom_point(size = 1.5, data=marker_genes_mat,color="#DC143C" )+
  theme_bw() +
  theme(#????һ??û??Ҫ????
    legend.background=element_blank(), 
    legend.key=element_blank(),
    #legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),#????????
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(limits = c(-4,12),breaks = seq(-4,12,by=4))+ #?޸?????
  scale_y_continuous(limits = c(-4,12),breaks = seq(-4,12,by=4))+
  coord_fixed()+#??ͼƬʼ?ձ???ԭ?г?????????????
  labs(x="Day2_iKlf5_vs_iEmpty_LFC",y="Day2_iKlf2_vs_iEmpty_LFC",
       title=paste(title,"(n=",dim(plot_mat_s)[1],")",sep = ""))+ #?޸???????????
  #guides(color=guide_legend(title = "Groups"))
annotate("text", label =paste("P.adj<0.1 n=",length(k5up),sep=""), x =9.2, y = c(-4), size = 4, colour = "black")+
annotate("text", label =paste("P.adj<0.1 n=",length(k2up),sep=""), x =c(-1), y = 12, size = 4, colour = "black")+
scale_color_manual(values = c("P.adj<0.1"=color,
                                "No.Sig"="#CCCCCC"))



#??ͼ?б??ǳ?marker gene
pp=p+ggrepel::geom_text_repel(data=marker_genes_mat,aes(label = gene,x = K5_active_lfc, y = K2_active_lfc),
                              fontface="bold",#?????Ӵ? 
                              size = 4,
                              box.padding = 0.3, #?ֵ????ľ???
                              point.padding = 0, #?ֵ????ľ??룬????Χ?Ŀհ׿???
                              min.segment.length = 0, #???߶ο???ʡ??
                              segment.color = "black", #segment.colour = NA, ????ʾ?߶?
                              colour = "black",inherit.aes = F,show.legend=F)
return(pp)
  }


scaterrplot(plot_mat,study_gene,"Nu",brewer.pal(n = 6,name = "Set2")[1],"ESC_Specific_Genes")
scaterrplot(plot_mat,study_gene,"Pu",brewer.pal(n = 6,name = "Set2")[3],"PGC_Specific_Genes")
scaterrplot(plot_mat,study_gene,"Tu",brewer.pal(n = 6,name = "Set2")[2],"TSC_Specific_Genes")
scaterrplot(plot_mat,study_gene,"NP",brewer.pal(n = 6,name = "Set2")[5],"ESC&PGC_Genes")
scaterrplot(plot_mat,study_gene,"NT",brewer.pal(n = 6,name = "Set2")[6],"ESC&TSC_Genes")
scaterrplot(plot_mat,study_gene,"TP",brewer.pal(n = 6,name = "Set2")[4],"PGC&TSC_Genes")




####3.Bubbleplot####
wnt=c("Axin2","Wnt3a","Wnt5b","Wnt6","Wnt10a","Fzd9","Fzd10","Tcf7",
      "Hoxc3","Hoxb4","Hoxa5")
x2c <- c(
  "Abcb5", "Ddr2", "Usp17lc", "Ddit41", "Chrm3", "Zfp560", "Usp17le",
  "Zfp809", "Pik2", "Snai1", "Gpx3", "Usp17la", "Cdkn1a", "Hspa1a",
  "Btg2", "Mdm2", "Emp3", "Mafk", "Sp110", "Rab43", "Ctsb", "Igfbp7",
  "Usp17ld", "Dcc", "Btg1", "Olfr376", "Zfp352", "Pnp", "Gm40910",
  "Pid1", "Lbp", "Bmpr2","Cpd", "Hmgn3", "Pdrg1", "Egfr", "Duoxa2", "Tmem132c", "Rfpl4b",
  "Obox3", "Trp53inp2", "Gm11487", "Zscan4e", "Tcstv3", "Zscan4a",
  "Zscan4f", "Zscan4c", "Zscan4b", "Zscan4d", "Mmp19", "Prnp",
  "4933440M02Rik", "4931428L18Rik", "Tgfbr3", "Tcstv1", "Zfp760",
  "Nrp1", "Fgf7", "Ptger4", "Cdkn2c"
)


PGC_all=c("Tfap2c", "Id3", "Rhox6", "Irf1", "Sox15", "Rhox9", "Dppa3", "Zic3", "Kit", 
          "Prdm1", "Rbm47", "Phf19", "Msx1", "Nanos3", "Msx2", "Prdm14", "Dazl", "Ddx4",
          "Dmrt1", "Dnd1")

TSC_all=c("Gata3","Gata2","Cdx2","Krt7","Krt18","Krt19","Ddah1","Rin3","Plet1","Hand1","Elf5",
          "Tead4","Eomes","Ets2","Fgfr2","Yap1","Tfap2c","Erf","Ly6a","Rhox5")


PGC_all=c("Prdm14","Prdm1","Dppa3","Nanos3","Kit","Tfap2c","Ddx4","Rhox6","Rhox9",
          "Sik1","Phf19","Wnt7b",
          "Irf1","Rbm47","Dazl","Ki67",
          "Msx1","Prdm1","Tfap2c","Msx2","Zic3","Dlx5",
          "Id3","Sox15")

Klfs=c("Klf2","Klf5")
ESC_Marker=c("Nanog","Tfcp2l1","Zfp42","Klf4","Tbx3")
PGC_Marker=c("Dppa3","Prdm1","Nanos3","Sox15","Id3")
TSC_and_ESC=c("Esrrb","Socs3","Tcf7l1","Zbtb18")
ESC_and_PGC=c("Pou5f1","Dppa5a","Nr5a2","Prdm14")
PGC_and_TSC=c("Tfap2c","Rhox6","Rhox9","Msx2")
TSC_Marker=c("Gata3","Cdx2","Krt7","Ddah1","Rin3","Plet1","Hand1","Elf5")

ETS=rev(c("Esrrb","Dppa4","Sox15"))
TSC_and_ESC=rev(c("Sox2","Ly6a","Tfeb"))
ESC_and_PGC=rev(c("Pou5f1","Nanog","Dppa5a","Nr5a2"))
PGC_and_TSC=rev(c("Tfap2c","Rhox6","Rhox9","Msx2"))
ESC_Marker=rev(c("Tfcp2l1","Tbx3","Zfp42"))
TSC_Marker=rev(c("Gata3","Cdx2","Fgfr2","Krt19","Hand1","Elf5"))
PGC_Marker=rev(c("Prdm1","Dppa3","Sox15","Dmrt1","Nanos3","Dazl"))
primde=rev(c("Fgf5","Otx2","Dnmt3a","Dnmt3b","Lefty1","T","Lin28a"))


en=c("Gdf3","Nckap1","Nodal","Ext1","Arc","Hhex","Gsc","Hnf1b","Sox17","Gata4",
     "Gata6","Foxa2","Cfc1","Sox7","Gata5",
     "Claudin6","Dab2","Foxa2","Fxyd3","Tm4sf20","Ttf1","Amn","Cxcr4","Emp2","Sdc4" #K4,K5,KD NC????
)
me=c("Srf","Bmpr2","Eomes","T","Msx2","Wnt3","Cebpa","Bmp6","Cdx4","Mixl1","Runx1","Tal1",
     "Flk1","Myt6","Myod1","Spef2","Opn4","Tbx20","Vegfr1","Bmp2","Hand1" #K4,K5,KD NC????
)
ex=c("Tcf15","Pax6","Dsp","Nes","Otx2","Gja1","Krt18","Tgm3","Col5a1","Sox18",
     "Runx3","Fgf5","Gfap",
     "Msi2","Cdh2","Neurod1"#K4,K5,KD NC????
)

diff=c("Col1a1","Col4a5","Col5a1","Tubb3","Nrcam","Thy1","Cdh2",
       "Hoxa1","Hoxa2","Hoxd3","Junb","Runx1","Tal1","Cxcr4","Hnf1b",
       "Pax6","Sox1","Neurod1","Neurod4"
       )

study_gene=c("Klf5","Prdm1","Dppa3","Nanos3","Fgf5","Otx2")
study_gene=x2c

study_list=list(ESC=ESC_Marker,
                TSC=TSC_Marker,
                PGC=PGC_Marker,
                `E&T`=TSC_and_ESC,
                `E&P`=ESC_and_PGC,
                `T&P`=PGC_and_TSC,
                `E&T&P`=ETS)

study_list=list(ESC=ESC_Marker,
                `E&P`=ESC_and_PGC,
                PGC=PGC_Marker,
                `T&P`=PGC_and_TSC,
                TSC=TSC_Marker)

study_list=study_gene


study_gene=as.character(unlist(study_list))

tpm_mat_mean=read.csv("20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_mat_mean)

c("iEmpty","iEmpty","iK2","iK2","iK4","iK4","iK5","iK5")#D25
c("iEmpty","iK2","iK4","iK5")#D2
c(16,4,8,12)#D2
c(16:17,4:5,8:9,12:13)#D25
c(14,2,6,10)#D2
c(14:15,2:3,6:7,10:11)#D25
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%study_gene,c(16,4,8,12)]
colnames(using_mat)=c("iEmpty","iK2","iK4","iK5")
colnames(using_mat)=c("iEmpty_LIF_D2","iEmpty_LIF_D5",#"iEmpty_CL_D2","iEmpty_CL_D5",
                      "iKlf2_LIF_D2","iKlf2_LIF_D5",#"iKlf2_CL_D2","iKlf2_CL_D5",
                      #"iKlf4_LIF_D2","iKlf4_LIF_D5",#"iKlf4_CL_D2","iKlf4_CL_D5",
                      "iKlf5_LIF_D2","iKlf5_LIF_D5")#"iKlf5_CL_D2","iKlf5_CL_D5")
colnames(using_mat)=c("iEmpty_LIF_D2","iEmpty_CL_D2","iKlf2_LIF_D2","iKlf2_CL_D2",
                      "iKlf4_LIF_D2","iKlf4_CL_D2","1GGGGgGGGiKlf5_LIF_D2","1GGGGgGGGiKlf5_CL_D2")
colnames(using_mat)=c("iEmpty_CL_D2","iEmpty_CL_D5","iKlf2_CL_D2","iKlf2_CL_D5",
                      "iKlf4_CL_D2","iKlf4_CL_D5","iKlf5_CL_D2","1GGGGgGGGiKlf5_CL_D5")

# using_mat[rownames(using_mat)=="Klf2",3:4]=0
# using_mat[rownames(using_mat)=="Klf4",5:6]=0
# using_mat[rownames(using_mat)=="Klf5",5:6]=0

tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,1]+0.1))
tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))

column_split_mat=data.frame(sample=colnames(using_mat),
                            split=rep("PNT_D2",dim(using_mat)[2]))
column_split_mat=data.frame(sample=colnames(using_mat),
                            split=c(rep("iEmpty",2),rep("iKlf2",2),
                                    rep("iKlf4",2),
                                    rep("iKlf5",2))
)

# tpm_mean_scale[rownames(tpm_mean_scale)=="Klf2",3:4]=NA
# tpm_mean_scale[rownames(tpm_mean_scale)=="Klf4",5:6]=NA
# tpm_mean_scale[rownames(tpm_mean_scale)=="Klf5",5:6]=NA


finmat=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_fold,remove_na = T,by_step = 2,cut_max_min = "Yes",
                               study_gene = study_gene,bubble_scale_by_new_function = F,limit_exp = c(0,9),
                               limit_fold = c(-2,2),column_split_mat = column_split_mat,
                               title = "K25_PNT_Day25",z_score =F,return = )#10
#3.6x5.5 ver3
#3.6x6.5 ver4


using_mat_bind=using_mat
using_mat_bind[,]=0
using_mat_bind[rownames(using_mat_bind)%in%K5_L_anno,1]=1
using_mat_bind[rownames(using_mat_bind)%in%K2_L_anno,2]=1
using_mat_bind$gene=rownames(using_mat_bind)
using_mat_bind$gene=factor(rownames(using_mat_bind),levels = study_gene)
using_mat_bind=using_mat_bind[order(using_mat_bind$gene),]
using_mat_bind$gene=NULL

pheatmap(using_mat_bind,cluster_rows = F,cluster_cols = F,border_color = "black",
         cellwidth = 20,cellheight = 20,#gaps_row = c(3,7,13,17),
         color = c("#D4D4D4","#F54040"))




tpm_mat_mean=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_mat_mean)

klfs_family=c()
for(i in 1:17){
  klf=paste("Klf",i,sep = "")
  klfs_family=c(klfs_family,klf)}

using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%klfs_family,c(19,25,1,67)]
colnames(using_mat)=c("ESC_1","ESC_2","EpiSC_1","EpiSC_2")


# using_mat[rownames(using_mat)=="Klf2",3:4]=0
# using_mat[rownames(using_mat)=="Klf4",5:6]=0
# using_mat[rownames(using_mat)=="Klf5",5:6]=0

tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,1]+0.1))
tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))

column_split_mat=data.frame(sample=colnames(using_mat),
                            split=c(rep("ESC",2),rep("Epi",2)))

# tpm_mean_scale[rownames(tpm_mean_scale)=="Klf2",3:4]=NA
# tpm_mean_scale[rownames(tpm_mean_scale)=="Klf4",5:6]=NA
# tpm_mean_scale[rownames(tpm_mean_scale)=="Klf5",5:6]=NA

Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,remove_na = F,by_step = 1,#cut_max_min = "Yes",
                            study_gene = klfs_family,bubble_scale_by_new_function =T,limit_exp = c(1,9),
                            limit_fold = c(-2,2),column_split_mat = column_split_mat,
                            title = "ESC_Epi",z_score =T)#10




using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%study_gene,c(#1,116,
                                                              19,20,22,23,106)]
using_mat=data.frame(TSC=apply(using_mat[,3:4],1,mean),
                     PGCLC=using_mat[,5],
                     ESC=apply(using_mat[,1:2],1,mean))
colnames(using_mat)=c(#"EpiSC_1","EpiSC_2",
                      "ESC_1","ESC_2","TSC_1","TSC_2",
                      "PGC_diff_PGCLC_d6"#,"BV+SC+"
                      )
tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,2]+0.1))
tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))

column_split_mat=data.frame(sample=colnames(using_mat),
                            split=rep("Control",3))
column_split_mat=data.frame(sample=colnames(using_mat),
                            split=c(#"EpiSC","EpiSC",
                                    "ESC","ESC","TSC","TSC",
                                    rep("PGCLC",2)))

Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,by_step = 0.75,
                               study_gene =  study_list,limit_exp = c(2,12),
                               limit_fold = c(-1.5,1.5),column_split_mat = column_split_mat,
                               col_vector =c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF","#FFFFFF","#D7B5D8","#DF65B0","#DD1C77","#980043"),
                               title = "Control",z_score =T)#10
#3.6x5.5 ver3
#3.6x6.2 ver4




####4.DEG vs Epi Motif and ESC TSC ATAC motif####

esc_atac1=readPeakFile("GSE110950_ESC_ATAC_merge_rmdupchrm_q0_sorted_peaks.narrowPeak")
tsc_atac1=readPeakFile("GSE110950_TSC_ATAC_merge_rmdupchrm_q0_sorted_peaks.narrowPeak")
esc_atac2=readPeakFile("GSE175632_ESC_ATAC_merge_rmdupchrm_q0_sorted_peaks.narrowPeak")
tsc_atac2=readPeakFile("GSE175632_TSC_ATAC_merge_rmdupchrm_q0_sorted_peaks.narrowPeak")


Make_2_Sample_Overlaps_Clustered_Bedfile(peak1 = esc_atac1,peak2 = esc_atac2,name1 = "esc_atac1",name2 = "esc_atac2")
Make_2_Sample_Overlaps_Clustered_Bedfile(peak1 = tsc_atac1,peak2 = tsc_atac2,name1 = "tsc_atac1",name2 = "tsc_atac2")

Split_Bed_by_Groups(list.files(pattern = "title"),filetype = ".bed")



setwd("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/ESC_TSC_ATAC_RNAseq_Motif")

show_motif=c("KLF5(Zf)","Sp5(Zf)",
             "Tcfcp2l1(CP2)","TFE3(bHLH)","Usf2(bHLH)",
             "NRF1(NRF)","Klf9(Zf)",
             "Fos(bZIP)","Jun-AP1(bZIP)","TEAD1(TEAD)","TEAD3(TEA)",
             "Maz(Zf)","Eomes(T-box)","Sox2(HMG)","Zfp281(Zf)","ZNF467(Zf)")

total_mat=Compare_Homer_known_motif_results(select_sites = c(".txt"),
                                            select_condition = "1000_100",
                                            study_motif = show_motif,
                                            #return="total_matrix",
                                            #rename_vector = c("ESC_unique","Common","TSC_unique"),
                                            show_top_sig._motif = 0)

plot_mat=total_mat[,c(2,4,6,1,3,5)]
colnames(plot_mat)=c("Day5","Day10")

plot_mat$order=factor(rownames(plot_mat),levels = show_motif)
plot_mat=plot_mat[order(plot_mat$order),]
plot_mat$order=NULL

write.csv(plot_mat,"Heatmap_plot_mat2.csv")
#plot_mat=read.csv("Heatmap_plot_mat2.csv",header = T,row.names = 1)


display_mat=plot_mat
display_mat[sapply(display_mat,function(x){x<=2})]<-0
display_mat[sapply(display_mat,function(x){x>2})]<-"*"
display_mat[sapply(display_mat,function(x){x==0})]<-""

colnames(plot_mat)=c("ESC_vs_Epi_up","PGCLC_vs_Epi_up","TSC_vs_Epi_up",
                     "ESC_vs_Epi_dn","PGCLC_vs_Epi_dn","TSC_vs_Epi_dn")

bk=seq(0,6,by=0.1)
p=pheatmap::pheatmap(plot_mat,
                     cluster_row =T,
                     cluster_cols = F,
                     border=T,
                     #breaks=bk,
                     border_color = "black",
                     show_rownames = T,
                     show_colnames = T,
                     clustering_method = "ward.D2",
                     cellwidth =15,
                     cellheight = 15,
                     fontsize = 15,
                     cutree_rows = 5,
                     gaps_col = 3,
                     #gaps_row = c(1,5),
                     legend_labels = "-log10 P-value",
                     display_numbers =as.matrix(display_mat),
                     main = "Motif Discovery",
                     #filename = pdf,
                     color=c(colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))(length(bk))))





####5.Lineage Heatmap in bulk RNAseq####

formative=c("Utf1","Pou3f1","Fgf5","Sox3","Lef1","Dnmt3a","Dnmt3b","Dnmt3l",
            "Tead2","Otx2","Sox12","Pim2","Sall2","Hes6","Zic2","Zic3","Zic5","Lin28a",
            "Lin28b","Zyx","Arid2","Foxp1","Cebpb","Mapk12","Fzd2","Rspo1","Smo")
study_gene=read.csv("D:/????/Data processing/#Klf5 Figure#/Figure.6 iK5 SGC PNT D0-5 iK2,4 D5 scRNAseq/marker_anno.csv")


colnames(tpm_mat_mean)
dynamic=read.csv("D:/????/Data processing/#Klf5 Figure#/Figure.6 iK5 SGC PNT D0-5 iK2,4 D5 scRNAseq/Lineage12_ht_feature.csv")
head(dynamic$X)
study_gene=study_gene$X
study_gene=c(formative,"Prdm1","Dppa3","Nanos3")
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%study_gene,c(103:106)]
using_mat$gene=factor(rownames(using_mat),levels = study_gene)
using_mat=using_mat[order(using_mat$gene),]
using_mat$gene=NULL

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(using_mat,cluster_rows = F,cluster_cols = F,
                   scale = "row",show_rownames = T,#annotation_row = anno,
                   breaks = bk,main="d2",gaps_col = 4,
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))

#overlap with K2,4,5 d2,5 up gene ####
d2_K5_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf5 _1.5/Day2_L_iEmpty vs Day2_L_iKlf5_degene_list.txt",header = T)

d2_K5_l_up_1.5=intersect(d2_K5_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))
d2_K5_l_dn_1.5=intersect(d2_K5_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,12)],1,function(x) max(x)>2),]))

d2_K2_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf2 _1.5/Day2_L_iEmpty vs Day2_L_iKlf2_degene_list.txt",header = T)

d2_K2_l_up_1.5=intersect(d2_K2_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))
d2_K2_l_dn_1.5=intersect(d2_K2_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,4)],1,function(x) max(x)>2),]))

d2_K4_l_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/Day2_L_iEmpty _vs_ Day2_L_iKlf4 _1.5/Day2_L_iEmpty vs Day2_L_iKlf4_degene_list.txt",header = T)

d2_K4_l_up_1.5=intersect(d2_K4_l_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))
d2_K4_l_dn_1.5=intersect(d2_K4_l_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(16,8)],1,function(x) max(x)>2),]))


d25_K2_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/K2_L_D2_vs_K2_L_D5_1.5/K2_L_D2_vs_K2_L_D5_1.5_degene_list.txt",header = T)
d25_K4_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/K4_L_D2_vs_K4_L_D5_1.5/K4_L_D2_vs_K4_L_D5_1.5_degene_list.txt",header = T)
d25_K5_1.5=read.table("D:/????/Data processing/20220205 R,K245 LIF CHIR-L Day2,5 RNA-seq/DEG_analysis/K5_L_D2_vs_K5_L_D5_1.5/K5_L_D2_vs_K5_L_D5_1.5_degene_list.txt",header = T)


d25_K2_5up_1.5=intersect(d25_K2_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(4,5)],1,function(x) max(x)>2),]))
d25_K4_5up_1.5=intersect(d25_K4_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(8,9)],1,function(x) max(x)>2),]))
d25_K5_5up_1.5=intersect(d25_K5_1.5$Down_regulate,rownames(tpm_total[apply(tpm_total[,c(12,13)],1,function(x) max(x)>2),]))

d25_K2_5dn_1.5=intersect(d25_K2_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(4,5)],1,function(x) max(x)>2),]))
d25_K4_5dn_1.5=intersect(d25_K4_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(8,9)],1,function(x) max(x)>2),]))
d25_K5_5dn_1.5=intersect(d25_K5_1.5$Up_regulate,rownames(tpm_total[apply(tpm_total[,c(12,13)],1,function(x) max(x)>2),]))



d25_K2_25up_1.5=setdiff(d2_K2_l_up_1.5,d25_K2_5dn_1.5)
d25_K4_25up_1.5=setdiff(d2_K4_l_up_1.5,d25_K4_5dn_1.5)
d25_K5_25up_1.5=setdiff(d2_K5_l_up_1.5,d25_K5_5dn_1.5)


Venner_plot(list = list(d25_K2_25up_1.5=d25_K2_25up_1.5,
                        d25_K4_25up_1.5=d25_K4_25up_1.5,
                        d25_K5_25up_1.5=d25_K5_25up_1.5))
venn=get.venn.partitions(list(d25_K2_25up_1.5=d25_K2_25up_1.5,
                         d25_K4_25up_1.5=d25_K4_25up_1.5,
                         d25_K5_25up_1.5=d25_K5_25up_1.5))
names(venn$..values..)=venn$..count..

Study_Gene_Distribution_in_gene_list(study_gene,venn$..values..)

Compare_Cluster_GO_results(venn$..values..,simplify_ver2 = T)



using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%c(d25_K2_25up_1.5,d25_K4_25up_1.5,d25_K5_25up_1.5),c(4,5,8,9,12,13)]

d=Cluster_mat_by_Kmeans(using_mat[,c(2,4,6)],kmeans = 7)


anno=data.frame(kmean=d$kmean,row.names = rownames(d))

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(d[,1:c(dim(d)[2]-1)],cluster_rows = F,cluster_cols = F,
                   scale = "row", 
                   breaks = bk,
                   show_rownames = F,annotation_row = anno,
                   main="d2",
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))
rownames(d[d$kmean==3,])

# '-------------------------------- ####
# WGCNA and cytoscape GRNCOP2 hub genes ####

library(WGCNA)


tpm_total_raw=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new.csv",row.names = 1)
colnames(tpm_total_raw)

use_tpm=tpm_total_raw[,c(1:2,122:123,
                         40:43,#50:53,206:208,
                         46:49,215:217#,246:247
                         )]
use_tpm=use_tpm[apply(use_tpm,1,function(x) max(x)>10),]
colnames(use_tpm)
# batch=c(1,1,2,2,2,3,3,4,4,5,5,3,3,6,6,6,7,7)
# batch=c(1,1,2,2,2,rep(1,13))
# #rep(8,7)
# #mod=model.matrix(~tissue)
# data.frame(batch=batch,m=colnames(use_tpm))
# use_tpm_rm_batch=ComBat_seq(counts = as.matrix(use_tpm),batch =batch)
# use_tpm=use_tpm_rm_batch

sampleTree = hclust(dist(use_tpm), method = "average") 
par(cex = 0.6) 
par(mar = c(0,6,6,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1.5, cex.main = 2) 
abline(h = 10000, col = "red")##???и?
cl=as.data.frame(cutree(sampleTree,h = 25000))##???и?
table(cl)#???ݷ????????жϲü??߶?
cl$gene=rownames(cl)
need_gene=cl[cl[,1]==1,2]
length(need_gene)
setdiff(study_gene,need_gene)

#remove outline gene
use_tpm=use_tpm[rownames(use_tpm)%in%need_gene,]

use_tpm_log2=log2(use_tpm+0.1)

condition=data.frame(row.names = colnames(use_tpm_log2),
                     condition=c(rep("EpiSC",4),rep("ESC",4),
                                 rep("TSC",4),rep("PGC",5)))
condition=data.frame(row.names = colnames(use_tpm_log2),
                     condition=c(rep("EpiSC",4),rep("ESC",4),
                                 rep("TSC",4),rep("PGC",3)))

grouptype <- data.frame(condition=sort(unique(condition$condition)),
                        groupNo=1:length(unique(condition$condition)))

condition$conditionNo = "NA"
for(i in 1:nrow(grouptype)){
  condition[which(condition$condition == grouptype$condition[i]),'conditionNo'] <- grouptype$groupNo[i]}
condition


use_tpm_log2_t <- as.data.frame(t(use_tpm_log2))

gsg <- goodSamplesGenes(use_tpm_log2_t,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(use_tpm_log2_t,verbose = 3)
gsg$allOK


#?????ʿ?
  sampleTree <- hclust(dist(use_tpm_log2_t), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## ??????????״?????ͣ????????Ӷ?Ӧ??ɫ???鿴?Ƿ?????????
  sample_colors <- numbers2colors(as.numeric(factor(condition$condition)), 
                                  colors = rainbow(length(table(condition$condition))), 
                                  signed = FALSE)
  ## ??????Ʒ??ϵͳ??????????Ӧ??״
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )

  
  group_list <- condition$condition
  dat.pca <- PCA(use_tpm_log2_t, graph = F) 
  pca <- fviz_pca_ind(dat.pca,
                      title = "Principal Component Analysis",
                      legend.title = "Groups",
                      geom.ind = c("point","text"), #"point","text"
                      pointsize = 2,
                      labelsize = 4,
                      repel = TRUE, #??ǩ???ص?
                      col.ind = group_list, # ??????ɫ
                      axes.linetype=NA,  # remove axeslines
                      mean.point=F#ȥ?????????ĵ?
  ) +
    theme(legend.position = "none")+  # "none" REMOVE legend
    coord_fixed(ratio = 1) #?????????ݺ???
  pca
  
  datExpr <- use_tpm_log2_t
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  R.sq_cutoff = 0.9  #????R^2 cut-off??Ĭ??Ϊ0.85

  # Call the network topology analysis function
  #????power????ѡ????Χ
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "signed",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ?? -1
  #pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: Ѱ?ҹյ㣬ȷ??????powerȡֵ
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  #dev.off()

  
  sft$powerEstimate  #?鿴?��Ƶ?????power
  power = sft$powerEstimate
  #power=9
  # ????????????powerС??15??????????powerС??30?ڣ?û??һ??powerֵʹ
  # ?ޱ???????ͼ?׽ṹR^2?ﵽ0.8??ƽ??��?Ӷ???100???£???????????
  # ??????Ʒ????????Ʒ????̫????????????????ЧӦ????Ʒ?????Ի?ʵ????????
  # ????Ӱ??̫???????ɡ?????ͨ????????Ʒ?????鿴??????Ϣ???????쳣??Ʒ??
  # ??????ȷʵ?????????????????仯?????ģ?Ҳ????ʹ???????ľ???powerֵ??
  if(is.na(power)){
    # ?ٷ??Ƽ? "signed" ?? "signed hybrid"
    # Ϊ??ԭ?ĵ?һ?£???δ?޸?
    type = "signed"
    nSamples=nrow(datExpr)
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))      
                   )
    )
  }
  
  power=18
  
  #һ??????????Ȩ?????????磬ʶ??????ģ??
    net <- blockwiseModules(
      datExpr,
      power = power,
      maxBlockSize = ncol(datExpr),
      corType = "pearson", #Ĭ??Ϊ"pearson","bicor"?????ܿ?????Ⱥ????Ӱ??
      networkType = "signed",
      TOMType = "signed", 
      minModuleSize = 200,    ##Խ??ģ??Խ??
      mergeCutHeight = 0, ##Խ??ģ??Խ??
      numericLabels = TRUE, 
      saveTOMs = F,
      verbose = 3
    )
    
    table(net$colors) 
    # power: ??һ????????????ֵ
    # maxBlockSize:???????ܴ?????????ģ???Ļ?????��(Ĭ??5000),16G?ڴ????Դ???2??????
    # ??????Դ???????????????÷???һ??block???档
    # corType???????????Եķ???????ѡpearson(Ĭ??)??bicor?????߸??ܿ?????Ⱥ????Ӱ?졣
    # networkType:?????ڽӾ???ʱ???Ƿ??????????????ԣ?Ĭ??Ϊ"unsigned",??ѡ"signed", "signed hybrid"
    # TOMType??????TOM????ʱ???Ƿ??????????????ԣ?Ĭ??Ϊ"signed",??ѡ"unsigned"?????Ǹ???????ת?????ڽӾ???(Ȩ??)?ķǸ??ԣ???????Ϊ????ѡ??"signed"Ҳû??̫???????塣
    # numericLabels: ???????ֶ???????ɫ??Ϊģ???????֣???????????ת??Ϊ??ɫ
    # saveTOMs?????ķ?ʱ???ļ??㣬?ɴ洢??��??????ʹ?ã?
    # mergeCutHeight: ?ϲ?ģ??????ֵ??Խ??ģ??Խ??,һ??Ϊ0.25
    # minModuleSize: ÿ??ģ???????ٷŶ??ٸ????????趨Խ??ģ??Խ??
    # ????????????ģ???л?????Ŀ?Ķ??٣????????У??��α???Ϊ `1-????ģ????`??
    # **0 (grey)**??ʾ**δ**?????κ?ģ???Ļ?????
# 
#     #?ֲ??????????繹????һ?㲻??
# 
#       ## ??????Ȩ????????????Ϊ��????
#       ## 1. ?????ڽ?ֵ??Ҳ?Ǿ???��???????ڲ?ͬ??Ʒ?б???��?ı???????ϵ??(pearson correlation rho)??
#       ## 2. ????topology overlap similarity (TOM)?? ??TOM??ʾ��?????????????ṹ?ϵ??????ԣ???��?????????????????Ƶ??ڽ?????????��?????????????????໥???á?
#       
#       ###(1)???繹?? Co-expression similarity and adjacency 
#       adjacency = adjacency(datExpr, power = power) 
#       
#       ###(2) ?ڽ??????????˾?????ת????Turn adjacency into topological overlap
#       TOM = TOMsimilarity(adjacency)
#       dissTOM = 1-TOM
#       
#       ###(3) ???????˾??? Clustering using TOM
#       # Call the hierarchical clustering function
#       geneTree = hclust(as.dist(dissTOM), method = "average");
#       # Plot the resulting clustering tree (dendrogram)
#       ## ????ʱ????geneTree??һ?????? net$dendrograms[[1]] ???????ƣ????ǻ???Ҫ???н?һ??????
#       #plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#       #     labels = FALSE, hang = 0.04)
#       
#       ###(4) ??????֧?????? dynamicTreeCut 
#       ################# set the minimum module size ##############################
#       minModuleSize = 200
#       ####
#       # Module identification using dynamic tree cut:
#       dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
#                                   deepSplit = 2, pamRespectsDendro = FALSE,
#                                   minClusterSize = minModuleSize)
#       table(dynamicMods)
#       # Convert numeric lables into colors
#       dynamicColors = labels2colors(dynamicMods)
#       table(dynamicColors)
#       # Plot the dendrogram and colors underneath
#       plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#                           dendroLabels = FALSE, hang = 0.03,
#                           addGuide = TRUE, guideHang = 0.05,
#                           main = "Gene dendrogram and module colors")
#       
#       ###(5) ????????????ģ?????ں? Merging of modules whose expression profiles are very similar
#       # Calculate eigengenes
#       MEList = moduleEigengenes(datExpr, colors = dynamicColors)
#       MEs = MEList$eigengenes
#       # Calculate dissimilarity of module eigengenes
#       MEDiss = 1-cor(MEs)
#       # Cluster module eigengenes
#       METree = hclust(as.dist(MEDiss), method = "average")
#       #һ??ѡ?? height cut Ϊ0.25,??Ӧ????75%?????ԣ??????ں?
#       ###################### set  Merging height cut  ################################
#       MEDissThres = 0.0
#       ####
#       # Plot the result
#       plot(METree, main = "Clustering of module eigengenes",
#            xlab = "", sub = "")
#       # Plot the cut line into the dendrogram
#       abline(h=0.5, col = "red")
#       # Call an automatic merging function
#       merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
#       # The merged module colors
#       mergedColors = merge$colors
#       # Eigengenes of the new merged modules:
#       mergedMEs = merge$newMEs
#       # ͳ??mergedmodule
#       
#       table(mergedColors)
#       
#       ### (6) plot the gene dendrogram 
#       plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#                           c("Dynamic Tree Cut", "Merged dynamic"),
#                           dendroLabels = FALSE, hang = 0.03,
#                           addGuide = TRUE, guideHang = 0.05)
# 
#       
#       ### ????????
#       # Rename to moduleColors
#       moduleColors = mergedColors
#       # Construct numerical labels corresponding to the colors
#       colorOrder = c("grey", standardColors(50))
#       moduleLabels = match(moduleColors, colorOrder)-1
#       MEs = mergedMEs
#       # Save module colors and labels for use in subsequent parts

#     
# #WGCNA???ӻ???TOMplot  Eigengene-adjacency-heatmap 
#   TOM=TOMsimilarityFromExpr(datExpr,power=power)
#   dissTOM=1-TOM
#   ## draw all genes 
# 
#     # geneTree = net$dendrograms[[1]]
#     # plotTOM = dissTOM^7
#     # diag(plotTOM)=NA
#     # TOMplot(plotTOM,geneTree,moduleColors,
#     #         col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
#     #         main="Network heapmap plot")
# 
#   ### draw selected genes to save time...just for test...
# 
#     # nSelect =0.1*nGenes
#     # set.seed(123)
#     # select=sample(nGenes,size = nSelect)
#     # selectTOM = dissTOM[select,select]
#     # selectTree = hclust(as.dist(selectTOM),method = "average")
#     # selectColors = moduleColors[select]
#     # plotDiss=selectTOM^7
#     # diag(plotDiss)=NA
#     # TOMplot(plotDiss,selectTree,selectColors,
#     #         col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
#     #         main="Network heapmap plot of selected gene")

    
    
    ## ģ?????ӻ????㼶??????չʾ????ģ??
    # Convert labels to colors for plotting
    moduleColors <- labels2colors(net$colors)
    table(moduleColors)
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    condition$condition <- as.factor(condition$condition)
    design <- model.matrix(~0+condition$condition)
    colnames(design) <- levels(condition$condition) #get the group
    MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
    MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
    moduleTraitCor <- cor(MEs,design,use = "p")
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
    textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                         signif(moduleTraitPvalue,1),")")
    dim(textMatrix) <- dim(moduleTraitCor)
    
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(design),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = F,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = F,
                   cex.text = 0.5,
                   zlim = c(-1,1), 
                   main = "Module-trait relationships")
    
### ģ????????չʾ Eigengene-adjacency-heatmap
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # ?????ӱ???????

    ## ��??????״
    #MET = orderMEs(cbind(MEs,condition$conditionNo))
    ## ??��??????״???轫?Ƿ????????????ͽ???0,1??ֵ?????Ѵ???design??
    design
    choose_celltype = as.data.frame(design[,1:4])
    names(choose_celltype) = colnames(design)
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, choose_celltype))

  plotEigengeneNetworks(MET, setLabels="",#plotHeatmaps = F,
                        marDendro = c(0,4,1,4),  # ???ף?????????
                        marHeatmap = c(5,5,1,2), # ???ף?????????
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
 
  

    check=as.data.frame(t(datExpr))
    check$moduleColors=moduleColors
    check$gene=rownames(check)
    checklist=Convert_Long_dataframe_to_list(check,group_index = dim(check)[2]-1,variable_index = dim(check)[2])
    
    table(check$moduleColors)
    Study_Gene_Distribution_in_gene_list(study_genes,checklist)
    
    module =c("magenta")
    module =c("lightcyan","grey60","tan")
    module =c("turquoise","lightgreen")
    module =c("red","blue","green")
    module=c("magenta","lightcyan","grey60","tan","turquoise",
                      "lightgreen","red","blue","green")
    module=c("lightgreen","magenta","lightyellow","darkred","midnightblue","purple","salmon","tan","lightyellow",
                         "blue","black","darkgreen","darkturquoise","green","royalblue")
    module=c("lightcyan","yellow","red","pink","purple","turquoise",
                        "green","cyan","greenyellow","lightyellow",
                        "royalblue","salmon")
    module2=c("lightcyan","green","greenyellow","purple","pink","cyan",
             "red","lightyellow","royalblue","salmon","yellow")
    module3=c("brown","darkred","lightyellow","pink","royalblue","yellow","lightgreen",
              "salmon","purple","blue")
      ### ????Ȥģ?黭??ͼ 

      dat=as.data.frame(t(datExpr[,moduleColors%in%module3]))
      dat=as.data.frame(t(datExpr[,moduleColors%in%"blue"]))
      dat=as.data.frame(t(datExpr[,colnames(datExpr)%in%c("Msx1")]))
      library(pheatmap)
      n=t(scale(t(dat))) #?Ի?????scale????ת?ñ???????Ϊ??Ϊ????????Ϊ??????ʽ
      # n[n>2]=2 
      # n[n< -2]= -2
      # n[1:4,1:4]
      
      group_list=condition$condition
      ac=data.frame(g=group_list)
      rownames(ac)=colnames(n)
      
      p=pheatmap::pheatmap(n[,],
                         fontsize = 8,
                         show_colnames =T,
                         show_rownames = F,
                         cluster_cols = F,
                         annotation_col =ac,
                         width = 8,
                         height = 6,
                         angle_col=45)
      setdiff(study_genes,rownames(n))

      
      # ????׼???????????????м???
      # ????KMEֵ
      datKME=signedKME(datExpr, MEs, outputColumnName="")
      # ??ȡ????Ȥ??module????????lightyellowΪ??
      datKME_select=datKME[,colnames(datKME)%in%module3]
      FilterGenes =rownames(datKME_select[apply(datKME_select,1, function(x) abs(max(x))> 0.7),])
      length(FilterGenes)
      # ???ؽ???
      table(FilterGenes)
      
      FilterGenes=as.data.frame(FilterGenes)
      FilterGenes$sum=apply(FilterGenes$kME_MM.lightyellow,1,sum)
      FilterGenes=FilterGenes[FilterGenes$sum>0,]
      FilterGenes=FilterGenes[FilterGenes$kME_MM.lightyellow==T,]
      # ????41??gene??????????
      hubgene <- rownames(FilterGenes)
      hubgene_lightyellow <- dimnames(data.frame(datExpr))[[2]][FilterGenes]
      setdiff(study_gene,hubgene)
      
      #select hubgene
      select_module_gene=rownames(as.data.frame(t(datExpr[,moduleColors%in%module3])))
      datKME=signedKME(datExpr, MEs, outputColumnName="")
      datKME_select=datKME[rownames(datKME)%in%select_module_gene,colnames(datKME)%in%module3]
      # ??ȡ????Ȥ??module????????lightyellowΪ??
      hubgene=c()
      for(i in 1:dim(datKME_select)[1]){
        if(max(abs(datKME_select[i,]))>0.65){
          hubgene=c(hubgene,rownames(datKME_select[i,]))
        }
      }
      length(hubgene)
      setdiff(study_genes,hubgene)

      
      turquoise_dat=as.data.frame(t(datExpr[,moduleColors%in%"turquoise"]))
      ETS=read.csv("ETP_gene_WGCNA_ver2_5254_.csv",row.names = 1)
      ETS=ETS[!rownames(ETS)%in%rownames(turquoise_dat),]
      
colnames(dat)

plot_mat=data.frame(EpiSCs=apply(dat[,c(1:2)], 1, mean),
                    ESC=apply(dat[,c(3:6)], 1, mean),
                    TSC=apply(dat[,7:10],1,mean),
                    PGCLC=apply(dat[,c(11:15)],1,mean))
plot_mat=data.frame(EpiSCs=apply(dat[,c(1:4)], 1, mean),
                    ESC=apply(dat[,c(5:8)], 1, mean),
                    TSC=apply(dat[,9:12],1,mean),
                    PGCLC=apply(dat[,c(13:15)],1,mean))

#plot_mat=plot_mat[rownames(plot_mat)%in%c(rownames(plotlist_mat),Epiu,ETP),]
plot_mat=plot_mat[rownames(plot_mat)%in%c(hubgene),]

d=Cluster_mat_by_Kmeans(plot_mat,scale_data = T,kmeans = 35,set.seed =1)

anno=data.frame(kmean=d$kmean,row.names = rownames(d))

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(d[,1:c(dim(d)[2]-1)],cluster_rows =F,cluster_cols = F,
                   scale = "row", 
                   breaks = bk,
                   show_rownames = F,annotation_row = anno,
                   main="d2",
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))
#rownames(d[d$kmean==4,])
dd=d
dd$gene=rownames(d)
checklist=Convert_Long_dataframe_to_list(dd,group_index = dim(dd)[2]-1,variable_index = dim(dd)[2])
table(dd$kmean)
Study_Gene_Distribution_in_gene_list(study_gene,checklist)

dat[rownames(dat)=="Nanos3",]



epi=c(3,4,8,13,32,34)
e=c(14,27,28)
t=c(2,11,12,16,19,21,22,23,25,26,35)
p=c(10,17,18,20,33,31)
ep=c(9,15,24)
et=c(1,5,6)
tp=c(7)
ept=c(29,30)

epi=d[d$kmean%in%epi,]
epi$kmean="Primed"
e=d[d$kmean%in%e,]
e$kmean="ESC"
t=d[d$kmean%in%t,]
t$kmean="TSC"
p=d[d$kmean%in%p,]
p$kmean="PGCLC"
ep=d[d$kmean%in%ep,]
ep$kmean="ESC&PGCLC"
et=d[d$kmean%in%et,]
et$kmean="ESC&TSC"
tp=d[d$kmean%in%tp,]
tp$kmean="TSC&PGCLC"
ept=d[d$kmean%in%ept,]
ept$kmean="ESC&TSC&PGCLC"

order_fun=function(plot_mat,coln,z_score=F){
  
  cluster=c()
  for(i in 1:length(levels(factor(plot_mat[,dim(plot_mat)[2]])))){
    sub_cl=levels(factor(plot_mat[,dim(plot_mat)[2]]))[i]
    cluster=c(cluster,rep(sub_cl,table(plot_mat[,dim(plot_mat)[2]])[i]))
  }
  
  row_plot_mat=plot_mat
  cl=plot_mat[,dim(plot_mat)[2]]
  cl=cl[!duplicated(cl)]
  
  plot_mat_new=plot_mat[1,]
  rownames(plot_mat_new)=1
  if(z_score==T){
    plot_mat=as.data.frame(t(scale(t(plot_mat[,1:c(dim(plot_mat)[2]-1)]))))
    plot_mat[,colnames(row_plot_mat)[dim(row_plot_mat)[2]]]=row_plot_mat[,dim(row_plot_mat)[2]]
  }
  for(c in cl){
    if(length(coln)==1){
      sub_plot_mat=as.data.frame(subset(plot_mat,plot_mat[,dim(plot_mat)[2]]%in%c))
      plot_mat_new=rbind(plot_mat_new,sub_plot_mat[order(sub_plot_mat[,coln],decreasing=T),])
    }else{
      sub_plot_mat=as.data.frame(subset(plot_mat,plot_mat[,dim(plot_mat)[2]]%in%c))
      sub_plot_mat$tmp_mean=apply(sub_plot_mat[,coln], 1, mean)
      sub_plot_mat=sub_plot_mat[order(sub_plot_mat$tmp_mean,decreasing=T),]
      plot_mat_new=rbind(plot_mat_new,sub_plot_mat[,1:dim(plot_mat)[2]])
    }
  }
  plot_mat_new=plot_mat_new[-1,]
  
  gene=rownames(plot_mat_new)
  plot_mat_new=Btable_insert_into_Atable(plot_mat_new[,1:2],row_plot_mat)
  plot_mat_new[,1:c(dim(plot_mat_new)[2]-1)]=as.data.frame(apply(plot_mat_new[,1:c(dim(plot_mat_new)[2]-1)], 2, as.numeric))
  rownames(plot_mat_new)=gene
  
  plot_mat_new=plot_mat_new[,3:c(dim(plot_mat_new)[2])]
  #plot_mat_new[,dim(plot_mat_new)[2]]=cluster
  return(plot_mat_new)
}

merge_mat=ETS
merge_mat=rbind(e,t,p,et,ep,tp,ept,epi)
#merge_mat=merge_mat[rownames(merge_mat)%in%c(rownames(plotlist_mat),Epiu,ETP),]
#merge_mat=merge_mat[rownames(merge_mat)%in%hubgene,]
merge_mat=order_fun(merge_mat,coln = c(1:4),z_score = T)
merge_mat$kmean=factor(merge_mat[,dim(merge_mat)[2]],levels=c("ESC","TSC","PGCLC","ESC&TSC","ESC&PGCLC","TSC&PGCLC","ESC&TSC&PGCLC",
                                                "Primed"))

rownames(merge_mat[merge_mat$kmean=="ESC&TSC",])

setdiff(study_genes,rownames(merge_mat))
Plot_clustered_heatmap_R4.2.2(merge_mat[,1:4],
                              colum_split_fun = c(1,1,1,1),
                              kmeans = 7,
                              name="ETP_gene_WGCNA_ver4",
                              clustering_method="complete",
                              clustered_and_order =merge_mat[,dim(merge_mat)[2]],
                              rename_cluster = F,
                              scale = "row",
                              set_seed = 169,
                              study_genes=c(study_gene,"Klf5","Klf2","Klf4",
                                            "Fgf5","Sox2","Otx2","Dnmt3a"),
                              dont_clustered = T,
                              show_study_genes = T,
                              show_row_names = F,
                              breaks=c(-2,2),#Default
                              height =7,#Default
                              width = 3.5,#Default
                              anno_font_size=8#use_raster = T,raster_quality = 25
)


anno=data.frame(kmean=merge_mat$kmean,row.names = rownames(merge_mat))

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(merge_mat[,1:c(dim(merge_mat)[2]-1)],cluster_rows =F,cluster_cols = F,
                   scale = "row",
                   breaks = bk,
                   show_rownames = F,annotation_row = anno,
                   main="d2",
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))



ETP1=read.csv("ETP_gene_WGCNA_ver3_4431_.csv")
ETP2=read.csv("ETP_gene_WGCNA_ver4_4240_.csv")

pgclc_relative1=ETP1[grepl("PGCLC",ETP1$Cluster),]
pgclc_relative2=ETP2[grepl("PGCLC",ETP2$Cluster),]


setdiff(pgclc_relative1$X,pgclc_relative2$X)

Venner_plot(list(a=pgclc_relative1$X,
            b=pgclc_relative2$X))
venn=get.venn.partitions(list(a=pgclc_relative1$X,
                         b=pgclc_relative2$X))
venn$..values..$`3`
venn$..values..$`2`




### ??ȡ????Ȥģ????????
datExpr_sub=datExpr[c(1:2,5:8,18:21,22:26),]
gene <- colnames(datExpr) 

celltype="ETP_2"
select_mat=merge_mat[grepl(celltype,merge_mat$kmean),]
select_mat=select_mat[!select_mat$kmean%in%"ESC&TSC&PGCLC",]
inModule <- gene%in%intersect(modgene,ETP)
inModule <- gene%in%rownames(merge_mat[merge_mat$kmean=="TSC",])
modgene <- gene[inModule]


gene <- colnames(datExpr) 
inModule <- moduleColors=="lightyellow"
modgene <- gene[inModule]
inModule <- gene%in%intersect(modgene,ETP)
modgene <- gene[inModule]
### ģ????Ӧ?Ļ?????ϵ????
TOM <- TOMsimilarityFromExpr(datExpr,
                             power=18,
                             corType = "pearson",
                             networkType = "signed",
                             TOMType = "signed", 
                             verbose = 3)
modTOM <- TOM[inModule,inModule]
dimnames(modTOM) <- list(modgene,modgene)

### ɸѡ��?Ӷ???????top100????
nTop = 100
IMConn = softConnectivity(datExpr[, modgene]) #????��?Ӷ?
top = (rank(-IMConn) <= nTop) #ѡȡ��?Ӷ???????top100
filter_modTOM <- modTOM[top, top]


# for cytoscape
cyt <- exportNetworkToCytoscape(filter_modTOM,
                                edgeFile = paste("step8_CytoscapeInput-edges-", paste(celltype, collapse="-"), ".txt", sep=""),
                                nodeFile = paste("step8_CytoscapeInput-nodes-", paste(celltype, collapse="-"), ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.15,  #weightedȨ??ɸѡ??ֵ???ɵ???
                                nodeNames = modgene[top])


edge=as.data.frame(fread("step8_CytoscapeInput-edges-ESC.txt"))
edge=edge[edge$fromNode%in%c("Klf5","Pou5f1","Nanog","Tbx3","Tfcp2l1"),1:3]
edge$target_color="#FFD700"
edge$target_size=60
#edge$node_color=GO_resluts
edge$node_size=40
write.table(edge,"edges-ESC.txt",sep = "\t",row.names = F,quote = F)

edge=as.data.frame(fread("step8_CytoscapeInput-edges-TSC.txt"))
edge=edge[edge$fromNode%in%c("Klf5","Cdx2","Elf5","Eomes","Gata3"),1:3]
edge$target_color="#FFD700"
edge$target_size=60
#edge$node_color=GO_resluts
edge$node_size=40
write.table(edge,"edges-TSC.txt",sep = "\t",row.names = F,quote = F)

edge=as.data.frame(fread("step8_CytoscapeInput-edges-PGCLC.txt"))
edge=edge[edge$fromNode%in%c("Klf5","Prdm14","Prdm1","Dppa3","Nanog"),1:3]
edge$target_color="#FFD700"
edge$target_size=60
#edge$node_color=GO_resluts
edge$node_size=40
write.table(edge,"edges-PGCLC.txt",sep = "\t",row.names = F,quote = F)


# '-------------------------------- ####
# K5_L binding distance to ETP gene ####

ETP1=read.csv("ETP_gene_WGCNA_ver3_4431_.csv")

CUTTag_bed=list.files("D:/????/Data processing/#Klf5 Figure#/All_CUTTag_bedfile/q0_sorted/",full.names = T)
K5_L=readPeakFile(CUTTag_bed[13])
K5_L_anno=Peak_Anno(K5_L,up = 20000,down = 20000,return = "total_matrix")


mm10_gene_region=read.table("D:/????/Data processing/Data_base/mm10_gene_region.bed")
all_gene=bitr(mm10_gene_region$V4,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Mm.eg.db")
all_gene=all_gene[!duplicated(all_gene$SYMBOL),]


K5_L_anno[K5_L_anno$SYMBOL%in%ETP1[grep("ESC",ETP1$Cluster),1],"ESC"]=1
K5_L_anno[K5_L_anno$SYMBOL%in%ETP1[grep("TSC",ETP1$Cluster),1],"TSC"]=1
K5_L_anno[K5_L_anno$SYMBOL%in%ETP1[grep("PGCLC",ETP1$Cluster),1],"PGCLC"]=1
K5_L_anno[is.na(K5_L_anno)]=0
K5_L_anno[apply(K5_L_anno[,25:27],1,function(x) sum(x)==0),"NO_ETP"]=1
K5_L_anno[is.na(K5_L_anno)]=0
tail(K5_L_anno)


# K5_L iEmpty_L ATAC and ETP gene ####

setwd("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/Epi_vs_iEmpty_iK5_ATAC_MAnorm")

iem_upbed=readPeakFile("Epi_FA_X_rep1_vs_iEmpty_L_rep1_M_below_-0.5849625_biased_peaks.bed")
ik5_upbed=readPeakFile("Epi_FA_X_rep1_vs_iKlf5_L_rep1_M_below_-0.5849625_biased_peaks.bed")

Make_2_Sample_Overlaps_Clustered_Bedfile(name1 = "iEmpty_o",
                                         name2 = "iK5_o",
                                         peak1 = iem_upbed,
                                         peak2 = ik5_upbed)


iK5_gain=readPeakFile("title_iEmpty_o_vs_iK5_o_100__iEmpty_o_vs_iK5_o_iK5_o_n_3389_.bed")
iK5_loss=readPeakFile("title_iEmpty_o_vs_iK5_o_100__iEmpty_o_vs_iK5_o_iEmpty_o_n_6961_.bed")


Make_2_Sample_Overlaps_Clustered_Bedfile(name1 = "iK5_gain",
                                         name2 = "K5_L",
                                         peak1 = iK5_gain,
                                         peak2 = K5_L)

Make_2_Sample_Overlaps_Clustered_Bedfile(name1 = "iK5_loss",
                                         name2 = "K5_L",
                                         peak1 = iK5_loss,
                                         peak2 = K5_L)


iK5_gain_anno=Peak_Anno(iK5_gain,up = 20000,down = 20000,return = "promoter_matrix")
iK5_loss_anno=Peak_Anno(iK5_loss,up = 20000,down = 20000,return = "promoter_matrix")


dim(iK5_gain_anno[iK5_gain_anno$SYMBOL%in%intersect(iK5_gain_anno$SYMBOL,d2_K5_l_dn_2),])
dim(iK5_gain_anno[iK5_gain_anno$SYMBOL%in%intersect(iK5_gain_anno$SYMBOL,d2_K5_l_up_2),])

dim(iK5_loss_anno[iK5_loss_anno$SYMBOL%in%intersect(iK5_loss_anno$SYMBOL,d2_K5_l_dn_2),])
dim(iK5_loss_anno[iK5_loss_anno$SYMBOL%in%intersect(iK5_loss_anno$SYMBOL,d2_K5_l_up_2),])




GREAT_GO=function(){}

study_term=c("keratinocyte differentiation","placenta development",
             "meiotic DNA recombinase assembly",
             "pattern specification process","central nervous system development",
             "epithelial cell differentiation","regulation of cell cycle",
             "stem cell population maintenance",
  "positive regulation of cell differentiation",
  "blood vessel development","skeletal system development","muscle structure development")
study_term=c("DNA replication-dependent nucleosome assembly","DNA repair","regulation of translation","cell cycle process","mitochondrial transmembrane transport","heterochromatin assembly",
             "maintenance of cell number","regulation of stem cell population maintenance","histone H3-K27 methylation","cellular response to leukemia inhibitory factor","regulation of cellular localization","fatty acid biosynthetic process",
             "cytoskeleton organization","embryonic morphogenesis","regulation of cell adhesion","embryonic organ development","Wnt signaling pathway","embryonic placenta morphogenesis")

plot_mat=Plot_GREAT_GO_results(path = "D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/Epi_vs_iEmpty_iK5_ATAC_MAnorm/GREAT_GO/",
                               study_term =study_term,data_type ="greatExportAll",top_go_term = 0,col_reoder = c(3,1,2),return_plot_or_matrix = "heatmap_matrix")

brewer.pal(5,"Blues")

bk=seq(0,25,0.1)
pheatmap::pheatmap(plot_mat,
                   cluster_row = F,
                   cluster_cols = F,
                   border=T,
                   breaks = bk,
                   border_color = "black",
                   show_rownames = T,
                   show_colnames = T,
                   cellwidth = 20,
                   cellheight = 20,
                   fontsize = 20,
                   #gaps_row = c(6,11),
                   legend_labels = "-log10 P-value",
                   color = c(colorRampPalette(colors = c("#FFFFFF"))((length(bk))*0.05),
                             colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.95))
)


# '-------------------------------- ####
# pySCENIC downstream ####
library(SCENIC)
library(AUCell)

setwd("D:/????/Data processing/#Klf5 Figure#/Figure.6.A K245_PNT_LIF_Comparation/ICM_TE_PGC_pySCENIC")

cellInfo=read.csv("GSE84892_cellInfo.csv",row.names = 1)

cellInfo=read.table("GSE45719_rawcounts_metadata.txt",row.names = 1)
                  
##1.AUCscore in different celltype####
AUC=as.data.frame(fread("Morula_ICM_TE_AUC.txt"))
rownames(AUC)=AUC$Cell
AUC$Cell=NULL
AUC=AUC[rownames(AUC)%in%rownames(cellInfo),]

motif=as.data.frame(fread("ITP_motifs.csv"))


##2.regulon in different celltype####
bin_mtx=read.csv("Zy_Morula_bin_mtx.csv",row.names = 1)
bin_mtx=bin_mtx[rownames(bin_mtx)%in%rownames(cellInfo),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Celltype),
                                     function(cells) rowMeans(plot_mat[cells,]))



plot_mat=t(bin_mtx)
colnames(plot_mat)=strsplit2(colnames(plot_mat),split = "\\.")[,1]
plot_mat_sub=as.data.frame(t(plot_mat[#cellInfo$Celltype%in%c("Morula")
                    ,colnames(plot_mat)%in%study_gene
                    ]))
plot_mat_sub=t(plot_mat[cellInfo$Celltype%in%c("E6.5_E7.5_SPGCs"),
  #colnames(plot_mat)%in%c("Prdm1","Klf5","Nanog","Klf4","Pou5f1","Otx2","Cdx2")
])
#plot_mat_sub=t(plot_mat)

anno=data.frame(celltype=cellInfo[,c(1)],row.names = rownames(cellInfo))

plot_mat_sub=plot_mat_sub[apply(plot_mat_sub, 1, function(x) sd(x)!=0),]
plot_mat_sub=plot_mat_sub[apply(plot_mat_sub, 1, function(x) mean(x)>0.05),]
plot_mat_sub$mean=apply(plot_mat_sub,1,mean)
plot_mat_sub=plot_mat_sub[order(plot_mat_sub$mean,decreasing = T),]

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(plot_mat_sub,#scale = "column",
                  show_colnames = F,cluster_cols = T,
                  #breaks = bk,
                  #color = c("white","black"),
                  #color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"))((length(bk))*0.5),
                  #           colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)),
                  annotation_col =anno)
mean(plot_mat_sub[rownames(plot_mat_sub)%in%c("Klf5"),])


a=plot_mat_sub[rownames(plot_mat_sub)%in%c("Klf5","Cdx2"),]
b=a[,a[2,]>0]

table(cellInfo[rownames(cellInfo)%in%colnames(b),"Celltype"])

adjacencies=as.data.frame(fread("ITP_adjacencies.txt"))
k5_adjacencies=adjacencies[adjacencies$TF=="Id1",]
k5_adjacencies[k5_adjacencies$target%in%study_gene,]



aucfile=list.files(pattern = "AUC.txt")[1:6]
names=strsplit2(aucfile,split = ".txt")[,1]
regulon_all=c()

for(n in 1:length(aucfile)){
  auc=aucfile[n]
  AUC=read.table(auc,row.names = 1,header = T,sep = "\t")
  regulon=strsplit2(colnames(AUC),split = "\\.")[,1]
  regulon_all=c(regulon_all,regulon)
}
regulon_freq=as.data.frame(table(regulon_all))
regulon_freq_filter=regulon_freq[regulon_freq$Freq>1,]
setdiff(study_gene,regulon_freq_filter$regulon_all)

newauc=as.data.frame(t(matrix(data=0,nrow=1,ncol=dim(regulon_freq_filter)[1])))
rownames(newauc)=regulon_freq_filter$regulon_all

for(n in 1:c(length(aucfile))){
  auc=aucfile[n]
  AUC=t(read.table(auc,row.names = 1,header = T,sep = "\t"))
  rownames(AUC)=strsplit2(rownames(AUC),split = "\\.")[,1]
  newauc=Btable_insert_into_Atable(newauc,AUC)
}
newauc[is.na(newauc)]=0
newauc$V1=NULL
newauc=as.data.frame(t(newauc))
rowna=rownames(newauc)
newauc=apply(newauc, 2, as.numeric)
rownames(newauc)=rowna

dim(newauc)
write.csv(newauc,"newauc.csv")


head(cellInfo)
cellInfo=Btable_insert_into_Atable(newauc[,1:2],cellInfo,index_B = 6)
cellInfo=cellInfo[,3:7]
write.csv(cellInfo,"newauc_cellInfo.csv")



##3.check regulon expr in different celltype####
cellInfo=read.csv("GSE148032_cellInfo.csv",row.names = 1)
regulon_expr=read.csv("PGC_regulonexpr.csv",row.names = 1)
regulon_AUC=read.csv("PGC_regulonAUC.csv",row.names = 1)

cellInfo=read.csv("GSE123046_cellInfo.csv",row.names = 1,check.names = T)
regulon_expr=read.csv("ICM_TE2_regulonexpr.csv",row.names = 1)
regulon_AUC=read.csv("ICM_TE2_regulonAUC.csv",row.names = 1)

cellInfo=as.data.frame(read.csv("GSE136714_cellInfo.csv",header = T))
rownames(cellInfo)=cellInfo$cellid
regulon_expr=read.csv("Zy_32c_regulonexpr.csv",row.names = 1)
regulon_AUC=read.csv("Zy_32c_regulonAUC.csv",row.names = 1)


dim(cellInfo)
dim(regulon_expr)
dim(regulon_AUC)

regulon_expr[regulon_expr<0.1]=0
regulon_expr[regulon_expr>0.1]=1


regulon_checkexpr_mat=as.data.frame(matrix(data=0,nrow = dim(regulon_expr)[2],ncol = length(levels(factor(cellInfo$Celltype)))))
colnames(regulon_checkexpr_mat)=levels(factor(cellInfo$Celltype))
rownames(regulon_checkexpr_mat)=colnames(regulon_expr)

for(i in 1:length(levels(factor(cellInfo$Celltype)))){
  celltype=levels(factor(cellInfo$Celltype))[i]
  cell=rownames(cellInfo[cellInfo$Celltype==celltype,])
  
  for(n in 1:dim(regulon_expr)[2]){

  regulon=colnames(regulon_expr)[n]
  sub_regulon_expr=regulon_expr[rownames(regulon_expr)%in%cell,]
  expr=mean(sub_regulon_expr[regulon_AUC[rownames(regulon_AUC)%in%cell,regulon]==1,regulon])
  if(is.na(expr)){expr=0}
  regulon_checkexpr_mat[n,i]=expr

  }
  print(paste("Finish...",celltype))
}


pheatmap::pheatmap(regulon_checkexpr_mat[rownames(regulon_checkexpr_mat)%in%study_gene,],
                   display_numbers = T,cluster_cols = F)

#write.csv(regulon_checkexpr_mat,"ICM_TE2_regulon_checkexpr_mat.csv")
#write.csv(regulon_checkexpr_mat,"PGC_regulon_checkexpr_mat.csv")

regulon_checkexpr_mat=read.csv("ICM_TE2_regulon_checkexpr_mat.csv",row.names = 1)

##4.defined top regulon in different celltype####
bin_mtx=t(read.csv("PGC_bin_mtx.csv",row.names = 1))
bin_mtx=t(read.csv("ICM_TE2_bin_mtx.csv",row.names = 1))
bin_mtx=t(read.csv("Zy_32c_bin_mtx.csv",row.names = 1))



regulonActivity <- sapply(split(rownames(cellInfo), cellInfo$Celltype),
                                     function(cells) rowMeans(bin_mtx[,cells]))
rownames(regulonActivity)=strsplit2(rownames(regulonActivity),split = "\\.")[,1]
regulonActivity=as.data.frame(regulonActivity)

test2=c("Msx2","Prdm1","Klf2","Isl2","Pole3","Sox15","Etv4",
        "Esrrb","Sap30","Jdp2","Cebpg","Pou3f2","Tfap2c","Tfe3",    
        "Zfp64","Hdac1","Snai3","Sox2","Stat6","Gcm1","Klf5")  
pheatmap::pheatmap(regulonActivity[rownames(regulonActivity)%in%study_gene,],cluster_cols = F,display_numbers = T)


filter_regulon_by_checkexpr=
function(regulonActivity,regulon_checkexpr_mat,col,top=dim(regulonActivity)[1]){
regulonActivity=regulonActivity[order(regulonActivity[,col],decreasing = T),]
sub_regulonActivity=regulonActivity[regulonActivity[,col]>0,]
regulon=intersect(rownames(sub_regulonActivity),
                  rownames(regulon_checkexpr_mat[regulon_checkexpr_mat[,col]>0.6,]))[1:top]
print(cat("\n",colnames(regulonActivity)[col],"\n",
  intersect(rownames(sub_regulonActivity),rownames(regulon_checkexpr_mat[regulon_checkexpr_mat[,col]<0.6,]))))
return(regulon)
}
#PGC>0.8
#ICM_TE>0.6
#Zy-Morula>0.8

filter_regulon_list=list(E3.5_ICM=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,1),
                      E3.5_TE=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,2),
                      E4.5_EPI=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,3),
                      E4.5_TE=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,4),
                      E5.5_EPI=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,5),
                      E5.5_ExE=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,6),
                      E6.5_EPI=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,7),
                      E6.5_ExE=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,8))
filter_regulon=unique(as.character(na.omit(unlist(filter_regulon_list))))

filter_regulon_list2=list(#E6.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,6),
                         E7.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,7),
                         #E8.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,8),
                         E9.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,9),
                         #E10.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,1),
                         E11.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,2),
                         #E12.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,3),
                         E13.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,4)
                         #E14.5_PGC=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,5),
                         )
filter_regulon2=unique(as.character(na.omit(unlist(filter_regulon_list2))))

filter_regulon_list3=list(
  Zy=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,7),
  Early_2cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,4),
  Mid_2cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,6),
  Late_2cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,5),
  a4cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,2),
  a8cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,3),
  Morula=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,1)
)
filter_regulon3=unique(as.character(na.omit(unlist(filter_regulon_list3))))

filter_regulon_list3=list(
  Zygote=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,7),
  a2cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,2),
  a4cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,4),
  Late4cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,6),
  a8cell=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,5),
  Morula=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,1),
  blast=filter_regulon_by_checkexpr(regulonActivity,regulon_checkexpr_mat,3)
)
filter_regulon3=unique(as.character(na.omit(unlist(filter_regulon_list3))))


##5.regulon activity ranking plot####

#filter_regulon_mat=as.data.frame(sapply(filter_regulon_list, "[", i = 1:max(sapply(filter_regulon_list, length))))
#write.csv(filter_regulon_mat,"ICM_TE2_sorted_regulon.csv")

#filter_regulon_mat2=as.data.frame(sapply(filter_regulon_list2, "[", i = 1:max(sapply(filter_regulon_list2, length))))
#write.csv(filter_regulon_mat2,"PGC_sorted_regulon.csv")

filtered_regulonActivity=regulonActivity[rownames(regulonActivity)%in%filter_regulon,]


write.csv(filtered_regulonActivity,"ICM_TE_filtered_regulonActivity.csv")
filtered_regulonActivity=read.csv("ICM_TE_filtered_regulonActivity.csv",row.names = 1)

Plot_Ranking=
  function(filtered_regulonActivity,regulonActivity,col,show_top=15,
           regulon_pool=NULL,return="plot",coord_fixed=100){
    filtered_regulonActivity$Regulon=rownames(filtered_regulonActivity)
    if(is.null(regulon_pool)){
    sub_filtered_regulonActivity=filtered_regulonActivity
    title=paste(colnames(sub_filtered_regulonActivity)[col],"(n=",length(na.omit( filtered_regulonActivity$Regulon)),")")
    }else{
      sub_filtered_regulonActivity=filtered_regulonActivity[filtered_regulonActivity$Regulon%in%regulon_pool,]
      title=paste(colnames(sub_filtered_regulonActivity)[col],"(n=",length(na.omit(regulon_pool)),")")
    }
    
    
    plot_mat=sub_filtered_regulonActivity[,c(col,dim(sub_filtered_regulonActivity)[2])]
    plot_mat=plot_mat[order(plot_mat[,1],decreasing = T),]
    plot_mat$Regulon=factor(plot_mat$Regulon,levels = plot_mat$Regulon)
    
    top_regulon=rownames(plot_mat[1:show_top,])
    regulon_highac=rownames(regulonActivity[regulonActivity[,col]>0.9,])
             
    colnames(plot_mat)=c("celltype","Regulon")
    top_regulon_mat=plot_mat[rownames(plot_mat)%in%c(top_regulon,regulon_highac),]
    
    
    p=ggplot(plot_mat,aes(x=Regulon,y=celltype))+
      geom_point(color='#666666')+
      ylim(0,1)+
      geom_point(size = 2.2, shape=16,data=top_regulon_mat,color="#DC143C")+
      geom_abline(intercept = 0,slope = 0,color="#444444",size=0.5, linetype="dotted")+
      geom_abline(intercept = 0.9,slope = 0,color="#DC143C",size=0.5, linetype="dotted")+
      coord_fixed(coord_fixed)+
      theme_bw()+
      theme(#axis.text.x = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background =element_blank())+
      ggtitle(title) +
      ylab(label = "Activity")
    
    pp=p+ggrepel::geom_text_repel(data=top_regulon_mat,
                                aes(label =Regulon ,x=Regulon,y=celltype),
                                fontface="italic",#?????Ӵ? 
                                size = 3.5,nudge_y = c(-0.05),force = T,max.overlaps = 100,
                                #direction = "x",
                                box.padding = 0.4, #?ֵ????ľ???
                                point.padding = 0, #?ֵ????ľ??룬????Χ?Ŀհ׿???
                                min.segment.length = 0, #???߶ο???ʡ??
                                segment.color = "black", #segment.colour = NA, ????ʾ?߶?
                                colour = "black",inherit.aes = F,show.legend=F)
    if(return=="top_gene"){
      return(top_regulon_mat$Regulon)
    }
    return(pp)
  }

p1=Plot_Ranking(filtered_regulonActivity,regulonActivity,1,regulon_pool = filter_regulon_list$E3.5_ICM)
p2=Plot_Ranking(filtered_regulonActivity,regulonActivity,3,regulon_pool = filter_regulon_list$E4.5_EPI)
p3=Plot_Ranking(filtered_regulonActivity,regulonActivity,5,regulon_pool = filter_regulon_list$E5.5_EPI)
p4=Plot_Ranking(filtered_regulonActivity,regulonActivity,7,regulon_pool = filter_regulon_list$E6.5_EPI)
p5=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list$E3.5_TE,col = 2,coord_fixed=90)
p6=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list$E4.5_TE,col = 4)
p7=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list$E5.5_ExE,col = 6)
p8=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list$E6.5_ExE,col = 8)

p9=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list2$E7.5_PGC,col = 7)
p10=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list2$E9.5_PGC,col = 9)
p11=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list2$E11.5_PGC,col = 2)
p12=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list2$E13.5_PGC,col = 4)

p13=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Zy,col = 7)
p14=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Early_2cell,col = 4)
p15=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Mid_2cell,col = 6)
p16=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Late_2cell,col = 5)
p17=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$a4cell,col = 2)
p18=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$a8cell,col = 3)
p19=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Morula,col = 1)

p13=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Zygote,col = 7)
p14=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$a2cell,col = 2)
p15=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$a4cell,col = 4)
p16=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Late4cell,col = 6)
p17=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$a8cell,col = 5)
p18=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$Morula,col = 1)
p19=Plot_Ranking(filtered_regulonActivity,regulonActivity,regulon_pool = filter_regulon_list3$blast,col = 3)

#save pdf 3.5x3.5

##5.regulon expr bubble plot####
library(data.table)
library(Seurat)
library(loomR)
library(SeuratObject)
library(ggplot2)
library(plyr)
library(limma)

seurat=readRDS("accu_E3.5_4.5_5.5_6.5_ICM_TE.Rds")
class(seurat)
seurat_ICM=seurat[,seurat@meta.data$Celltype%in%c("E3.5_ICM","E4.5_EPI","E5.5_EPI","E6.5_EPI")]
seurat_TE=seurat[,seurat@meta.data$Celltype%in%c("E3.5_TE","E4.5_TE","E5.5_ExE","E6.5_ExE")]

seurat=readRDS("PGC_E6.5-14.5.Rds")
# seurat_PGC=NormalizeData(seurat_PGC)
# all.genes <- rownames(seurat_PGC)
# seurat_PGC <- ScaleData(seurat_PGC, features = all.genes)
# saveRDS(seurat_PGC,"PGC_E6.5-14.5.Rds")
seurat_PGC=seurat[,seurat@meta.data$Celltype%in%c("E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC")]


color = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#FFFFFF"   ))((20)),
          colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026" ))((80)))

ICM_gene=rev(c("Klf5","Klf4","Nanog","Klf2","Pou5f1"))
TE_gene=rev(c("Klf5","Cdx2","Gata3","Erf","Tfap2a"))
PGC_gene=rev(c("Klf5","Prdm1","Nanog","Tfap2c","Klf2"))

gene=TE_gene
p <- DotPlot(seurat_TE, features=gene,
             scale.min = 0,
             scale.max = 100,
             scale = F,
             group.by='Celltype',
             cols =  c("RdYlBu"))
p+ggplot2::coord_flip()
p_data=p$data
p_data$features.plot=factor(gene,levels = gene)
#p_data$id=factor(p_data$id,levels =c("E7.5_PGC","E9.5_PGC","E11.5_PGC","E13.5_PGC"))

scale.by="radius"
scale.func <- switch(EXPR = scale.by, size = scale_size, #?ú??????Ըı?ԲȦ??С????ԲȦ֮???????????󣡣???
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
#limit_max=2
limit_max=2.5

p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
#?ú?????????ȫ????scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  scale.func(range = c(0,15))+#????С
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(0,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#???ӱ߿?
  #facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size = 20),
        axis.text.y = element_text(face=c("italic"),size = 20),
        strip.text.x = element_text(size=8))
#4.2x6.1


##6.regulon network####

filtered_regulonActivity1=read.csv("ICM_TE_filtered_regulonActivity.csv",row.names = 1)
filtered_regulonActivity2=read.csv("PGC_filtered_regulonActivity.csv",row.names = 1)

adjacencies=as.data.frame(fread("ICM_TE2_adjacencies.txt"))

adjacencies2=as.data.frame(fread("PGC_adjacencies.txt"))

adjacencies3=as.data.frame(fread("Zy_32c_adjacencies.txt"))

filter_regulon_list1=read.csv("ICM_TE2_sorted_regulon.csv",row.names = 1)
filter_regulon_list2=read.csv("PGC_sorted_regulon.csv",row.names = 1)

E3.5_ICM_top=filter_regulon_list1$E3.5_ICM[1:15]
E3.5_TE_top=filter_regulon_list1$E3.5_TE[1:15]
E9.5_PGC_top=filter_regulon_list2$E9.5_PGC[1:15]


ICM_gene=rev(c("Klf5","Klf4","Nanog","Klf2","Pou5f1"))
TE_gene=rev(c("Klf5","Cdx2","Gata3","Erf","Tfap2a"))
PGC_gene=rev(c("Klf5","Prdm1","Nanog","Tfap2c","Klf2"))

make_net=
function(adjacencies,top_regulon,marker_regulon,name){
net=adjacencies[adjacencies$TF%in%top_regulon,]
net=net[net$target%in%top_regulon,]
net$highlight="normal"
net[net$TF%in%marker_regulon,"highlight"]="highlight"
write.csv(x = net,file = name,row.names = F,quote = F)

}




make_net2=
  function(adjacencies,top_regulon,marker_regulon,
           filtered_regulonActivity,col,
           center_regulon="Klf5",name){
    net=adjacencies[adjacencies$TF%in%top_regulon,]
    net=net[net$target%in%top_regulon,]
    net$highlight="normal"
    net[net$TF%in%marker_regulon,"highlight"]="highlight"
    
    net[net$TF%in%center_regulon,"highlight"]="center"
    
    colname=colnames(net)
    
    filtered_regulonActivity$regulon=rownames(filtered_regulonActivity)
    
    net=Bmat_join_Amat(net,filtered_regulonActivity[,c(col,dim(filtered_regulonActivity)[2])],1,2)
    colnames(net)=c(colname,"activity")
    
    write.csv(x = net,file = name,row.names = F,quote = F)
    
  }

make_net2(adjacencies,E3.5_ICM_top,ICM_gene,filtered_regulonActivity1,1,name="E3.5_ICM_network2.csv")
make_net2(adjacencies,E3.5_TE_top,TE_gene,filtered_regulonActivity1,2,name="E3.5_TE_network2.csv")
make_net2(adjacencies2,E9.5_PGC_top,PGC_gene,filtered_regulonActivity2,2,name="E9.5_PGC_network2.csv")



top_regulon_mat=as.data.frame(matrix(data=0,nrow = length(unique(unlist(top_regulon_list))),ncol = length(top_regulon_list)))
colnames(top_regulon_mat)=names(top_regulon_list)
rownames(top_regulon_mat)=unique(unlist(top_regulon_list))

for(i in 1:length(top_regulon_list)){
  names(top_regulon_list)[i]
  for(gene in unique(unlist(top_regulon_list))){
    top_regulon_mat[rownames(top_regulon_mat)==gene,i]=sum(top_regulon_list[[i]]%in%gene)
  }
}
head(top_regulon_mat)


top_regulon_mat$group=apply(top_regulon_mat,1,function(x) paste(x,collapse = ","))
top_regulon_mat$sum=rowSums(top_regulon_mat[,1:6])


test=as.data.frame(table(top_regulon_mat$group))
top_regulon_mat[top_regulon_mat$group%in%test[test$Freq==1,1],]
top_regulon_mat[top_regulon_mat$sum==1,]


top_regulonActivity=Btable_insert_into_Atable(top_regulon_mat[,1:2],regulonActivity)
top_regulonActivity=top_regulonActivity[,3:dim(top_regulonActivity)[2]]
top_regulonActivity=apply(top_regulonActivity,2,as.numeric)

for(i in 1:dim(top_regulon_mat)[1]){
  top_regulon_mat[i,"activity_mean"]=mean(top_regulonActivity[i,which(top_regulon_mat[i,1:6]==1)])
}

top_regulon_mat_0.2=top_regulon_mat[top_regulon_mat$activity_mean>0.2,]

top_adjacencies=adjacencies[adjacencies$TF%in%rownames(top_regulon_mat_0.2),]
top_adjacencies=top_adjacencies[top_adjacencies$target%in%rownames(top_regulon_mat_0.2),]
top_adjacencies=top_adjacencies[top_adjacencies$importance>1,]

make_network=
  function(top_regulon_mat_0.2){
    
  }




make_omic=
  function(top_regulon_mat_0.2){
    omic_mat=data.frame(regulon="a",celltype="a")
    
    for(i in 1:dim(top_regulon_mat_0.2)[1]){
      regulon=rownames(top_regulon_mat_0.2)[i]
      celltype=rep("blank",6)
      celltype[which(top_regulon_mat_0.2[i,1:6]==1)]=colnames(top_regulon_mat_0.2)[which(top_regulon_mat_0.2[i,1:6]==1)]
      sub_omic_mat=data.frame(regulon=rep(regulon,6),celltype=celltype)
      omic_mat=rbind(omic_mat,sub_omic_mat)
    }
    omic_mat=omic_mat[-1,]
    return(omic_mat)
  }

omic_mat=make_omic(top_regulon_mat_0.2)





# '-------------------------------- ####
# Mouse_Organogenesis_E9.5_GSE87038 ####

library(Seurat)
library(loomR)
library(SeuratObject)
library(ggplot2)
library(plyr)
library(data.table)


setwd("D:/????/Data processing/#TSC_Epi_and_ES_gene_and_epigenetic#/All_Data/Mouse_Organogenesis_E9.5_GSE87038")


data=read.csv("GSE87038_Mouse_Organogenesis_expression_matrix_log2_based_TPM_divided_by_10.csv",header = T,row.names=1)
head(data[,1:10])

meta=as.data.frame(strsplit2(colnames(data),"_e"))
colnames(meta)=c("Celltype","cell")
head(meta)


seurat=CreateSeuratObject(data)
seurat@meta.data=meta


p <- DotPlot(seurat, features=c("Klf5","Cdx2"),
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
scale.func <- switch(EXPR = scale.by, size = scale_size, #?ú??????Ըı?ԲȦ??С????ԲȦ֮???????????󣡣???
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
limit_max=2
limit_max=2.5

name="Figure.1.A_Stringency_Single_cell_data_new_ver1"
pdf(paste(name,".pdf",sep=""),width=5.2,height =5)

p_data[p_data$avg.exp.scaled>=limit_max,"avg.exp.scaled"]=as.numeric(limit_max)
#?ú?????????ȫ????scale_size_...!!!!!
ggplot(p_data,aes(x=id,y=features.plot,fill=avg.exp.scaled))+
  geom_point(aes_string(size="pct.exp"),color="black",shape=21)+
  scale.func(range = c(0,6.5))+#????С
  guides(size = guide_legend(title = "Percent Expressed"))+
  scale_fill_gradientn(colors = color,limits=c(0,limit_max),
                       name = "Average Expression")+
  ylab("")+xlab("")+theme_bw()+
  # scale_y_continuous(breaks = 1:length(levels(p_data$features.plot)),
  #                    labels = levels(p_data$features.plot),
  #                    sec.axis = dup_axis())+#???ӱ߿?
  facet_grid(~label,scales = "free_x",space = "free")+theme_classic()+
  theme(panel.spacing=unit(2, "mm"),
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(face=c("italic")),
        strip.text.x = element_text(size=6.3))
dev.off()
dev.new()
























signal=read.csv("D:/????/Data processing/#Klf5 Figure#/Figure.6 iK5 SGC PNT D0-5 iK2,4 D5 scRNAseq/D5_K245_tSNE_cl/mouse.signaling.cgz.csv")

tpm_mat_mean=read.csv("D:/????/Data processing/20230104_Klf5_figure_all_RNAseq_tpm_new_mean_mat.csv",row.names = 1)
colnames(tpm_mat_mean)

c(16,4,8,12)#D2
c(16:17,4:5,8:9,12:13)#D25
c(14,2,6,10)#D2
c(14:15,2:3,6:7,10:11)#D25
using_mat=tpm_mat_mean[rownames(tpm_mat_mean)%in%signal$gene,c(16:17,12:13,14:15,10,11)]
using_mat$gene=rownames(using_mat)
using_mat$gene=factor(using_mat$gene,levels=signal$gene)
using_mat=using_mat[order(using_mat$gene),]
using_mat=Btable_insert_into_Atable(using_mat,signal,index_B = 1)

anno=data.frame(signal=using_mat$signaling,row.names = using_mat$gene)

bk=seq(-2,2,by=0.1)
pheatmap::pheatmap(using_mat[,1:8],cluster_rows = F,cluster_cols = F,display_numbers = round(using_mat[,1:8]),
                   scale = "row",show_rownames = T,annotation_row = anno,
                   breaks = bk,main="d2_5",gaps_col = c(4),
                   color  = c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#FFFFFF"  ))((length(bk))*0.5),
                              colorRampPalette(colors = c("#FFFFFF","#FDAE61","#F46D43","#D73027","#A50026"))((length(bk))*0.5)))



signal=signal[signal$gene%in% rownames(using_mat),]

signal_list=Convert_Long_dataframe_to_list(signal,2,1)

using_mat=using_mat[,1:4]
tpm_mean_fold=log2((using_mat+0.1)/(using_mat[,1]+0.1))
tpm_mean_scale=as.data.frame(t(scale(t(using_mat))))

tpm_mean_scale[tpm_mean_scale>2]=2
tpm_mean_scale[tpm_mean_scale<c(-2)]=-2

tpm_mean_fold[tpm_mean_fold>3]=3
tpm_mean_fold[tpm_mean_fold<c(-3)]=-3

p1=Plot_Flodchange_Dot_Heatmap(using_mat,tpm_mean_scale,remove_na = F,
                               study_gene =signal_list,
                               bubble_scale_by_new_function = F,
                               limit_exp = c(0,9),
                               limit_fold = c(-2,2),#column_split_mat = column_split_mat,
                               title = "K25_PNT_Day25",z_score =F)#10

