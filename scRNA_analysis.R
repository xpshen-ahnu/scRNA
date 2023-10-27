library(Seurat)
samples=list.files("data/")
samples
dir <- file.path('./data',samples)
names(dir) <- samples

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i], gene.column = 2)
  
    
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples[i],
                                       min.cells=3, min.features = 200)
  
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-") 
  }
  
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^Rp[sl]")
  }
 
}


#  -------------------------------------------------------------
# dir.create("Integrate")
setwd("F:/scRNA_analysis/Integrate/")

names(scRNAlist) <- samples

system.time(saveRDS(scRNAlist, file = "scRNAlist0.rds"))

scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA

table(scRNA$orig.ident)

save(scRNA,file = 'scRNA_orig.Rdata')
load('scRNA_orig.Rdata')


# ----------------------------------------------------------------------


feats <- c("nFeature_RNA", "nCount_RNA")
library(patchwork)

p_filt_b_1=VlnPlot(scRNA, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()

feats <- c("percent.mt","percent.rb")

p_filt_b_2=VlnPlot(scRNA, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()
p_filt_b_1 | p_filt_b_2
ggplot2::ggsave('quality_control_before.png')


#  ------------------------------------------------------------


retained_c_umi_low <- scRNA$nFeature_RNA > 200
retained_c_umi_high <- scRNA$nFeature_RNA < 8000

retained_c_mito <- scRNA$percent.mt < 5
retained_c_ribo <- scRNA$percent.rb > 10
table(retained_c_mito & retained_c_ribo & retained_c_umi_low & retained_c_umi_high)

table(scRNA$orig.ident[retained_c_mito & retained_c_ribo & retained_c_umi_low & retained_c_umi_high])


scRNA_filt = scRNA[, retained_c_mito & retained_c_ribo & retained_c_umi_low & retained_c_umi_high]
dim(scRNA_filt)

table(scRNA_filt$orig.ident)


feats <- c("nFeature_RNA", "nCount_RNA")
p_filt_a_1=VlnPlot(scRNA_filt, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()
feats <- c("percent.mt","percent.rb")
p_filt_a_2=VlnPlot(scRNA_filt, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()
p_filt_a_1 | p_filt_a_2
ggplot2::ggsave('quality_control_after.png')
saveRDS(scRNA_filt, file = "sce_filt.rds")

scRNA = readRDS("sce_filt.rds")

# --------------------------------------------


sce.all.filt <- scRNA
library(DoubletFinder)
library(tidyverse)
table(sce.all.filt@meta.data$orig.ident)
sce.all.list <- SplitObject(sce.all.filt, split.by = "orig.ident")
sce.all.list
dataL <- sce.all.list
for (idx in 1:length(dataL)){
  print(paste0('第',idx,'样品'))
  ## preprocessing
  data <- NormalizeData(dataL[[idx]])
  data <- ScaleData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, verbose = FALSE)
  data <- RunPCA(data, npcs = 30, verbose = FALSE)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
  data <- FindClusters(data, resolution = 0.5)
  sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- data@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)          
  nExp_poi <- round(0.04*nrow(data@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  ## save results
  dataL[[idx]]$doubFind_res = data@meta.data %>% select(contains('DF.classifications'))
  dataL[[idx]]$doubFind_score = data@meta.data %>% select(contains('pANN'))
}
head(dataL[[1]]@meta.data)
table(dataL[[1]]@meta.data$doubFind_res)  

#dir.create('F:/scRNA_analysis/2-int-doublet_removal')
setwd('F:/scRNA_analysis/2-int-doublet_removal')

dataL1 <- dataL
for (i in 1:length(dataL1)) {
  dataL1[[i]] <- dataL1[[i]][,dataL1[[i]]@meta.data$doubFind_res=='Singlet']
}
sce.all.list <- dataL1
save(sce.all.list,file = "F:/tmp/sce.RData")

# -------------------------------------------------------------------------


library(harmony)
load('sce.RData')
library(Seurat)
library(tidyverse)
sce <- merge(sce.all.list[[1]],sce.all.list[[2]])
sce@meta.data
scRNA <- sce
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", plot_convergence = TRUE) 
harmony_embeddings <- Embeddings(scRNA,"harmony")
harmony_embeddings[1:4,1:4]

scRNA <- scRNA %>% RunUMAP(reduction = "harmony",dims = 1:40) %>%
  FindNeighbors(reduction = "harmony",dims = 1:40)
scRNA <- FindClusters(scRNA, resolution = seq(from=0,by=.2,length=10))

library(clustree)
p_cutree = clustree(scRNA@meta.data, prefix = "RNA_snn_res.");p_cutree 


save(scRNA, file = "scRNA20230629.RData")


# -------------------------------------------------------------------------

load('scRNA20230629.RData')
library(Seurat)
library(tidyverse)

# -----------------------------------------------------------

retained_c_umi_low <- scRNA$nFeature_RNA > 1200

scRNA <- scRNA[,retained_c_umi_low]
p_dim = DimPlot(scRNA, group.by = "orig.ident");p_dim
p_dim1 = DimPlot(scRNA);p_dim1
scRNA <- scRNA %>% RunUMAP(reduction = "harmony",dims = 1:40) %>%
  FindNeighbors(reduction = "harmony",dims = 1:40)
scRNA <- FindClusters(scRNA, resolution = seq(from=0,by=.2,length=10))



Idents(scRNA) = scRNA@meta.data$RNA_snn_res.0.6
table(scRNA@active.ident)

table(scRNA@meta.data$orig.ident, scRNA@active.ident)



# ------------------------------------------------------
# dir.create('./cluster_marker')

setwd('./cluster_marker')
for (i in 0:8) {
  marker <- FindMarkers(scRNA, ident.1 = i)
  marker <- marker %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))
  write.csv(marker, file = paste0("cluster_",i,".csv"))
}
DimPlot(scRNA, label = TRUE, label.size = 6)
ggsave("UMAP_cluster.png")
ggsave("UMAP_cluster.pdf")

DimPlot(scRNA, label = TRUE, label.size = 6, split.by = "orig.ident")
ggsave("UMAP_cluster_split.png", width = 8)
ggsave("UMAP_cluster_split.pdf", width = 8)


save(scRNA,file = "scRNA20230703.RData")


# -----------------------------------------------------------------
scRNA$group=ifelse(scRNA@active.ident==6,0,
                   ifelse(scRNA@active.ident==0,1,
                          ifelse(scRNA@active.ident==1,2,
                                 ifelse(scRNA@active.ident==2,3,
                                        ifelse(scRNA@active.ident==3,4,
                                               ifelse(scRNA@active.ident==4,5,
                                                      ifelse(scRNA@active.ident==5,7,
                                                             ifelse(scRNA@active.ident==7,8,6))))))))
table(scRNA@active.ident,scRNA@meta.data$group)

DotPlot(scRNA, features = c("Pou5f1","Nanog","Sox2","Zfp42", "T","Mixl1","Mesp1","Eomes","Gsc",
                            "Btg2","Crip2","Cd24a","Nkx2-9","Snai3"), group.by = "group")
genes <- list(pluripotent=c("Pou5f1","Nanog","Sox2","Zfp42"),
              Mesendoderm=c("T","Mixl1","Mesp1","Eomes","Gsc"),
              ecotderm=c("Btg2","Crip2","Cd24a","Nkx2-9","Snai3"))

DotPlot(scRNA,features = genes,group.by = "group")
Idents(scRNA) = scRNA@meta.data$group
ggsave("marker_dotplot.png",width = 11)
ggsave("marker_dotplot.pdf",width = 11)


# annotation --------------------------------------------------------------

scRNA$celltype <- ifelse(scRNA@meta.data$group==0,"Undiff. ES",
                         ifelse(scRNA@meta.data$group==7,"Mesendoderm",
                                ifelse(scRNA@meta.data$group==8,"Ecotderm","Epiblast")))
table(scRNA$group,scRNA$celltype)

Idents(scRNA) = scRNA@meta.data$celltype
DimPlot(scRNA, label = TRUE, label.size = 5,pt.size = 1.5)+NoLegend()
ggsave("UMAP_anno.png")
ggsave("UMAP_anno.pdf")
genes <- list(pluripotent=c("Pou5f1","Nanog","Sox2","Zfp42"),
              Mesendoderm=c("T","Mixl1","Mesp1","Eomes","Gsc"),
              Ecotderm=c("Btg2","Crip2","Cd24a","Nkx2-9","Snai3"))

DotPlot(scRNA,features = genes,group.by = "celltype")
ggsave("marker_dotplot_anno.png",width = 11)
ggsave("marker_dotplot_anno.pdf",width = 11)


# ---------------------------------------------------

table(scRNA@meta.data$orig.ident,scRNA$celltype)

table(scRNA@meta.data$orig.ident)

a1 <- as.data.frame(table(scRNA@meta.data$orig.ident,scRNA$celltype))


a1$ratio <- ifelse(a1$Var1=="KO",a1$Freq/7207, a1$Freq/8841)
a1 <- a1 %>% group_by(Var1) %>% arrange(desc(Var2)) %>% mutate(lab.pos=cumsum(ratio)-0.5*ratio)
a1$ratio1 <- a1$ratio*100
a1$ratio1 <- round(a1$ratio1,digits = 2)
a1$ratio1 <- paste0(a1$ratio1,"%")
a1$Var1 <- factor(a1$Var1,levels = c("WT","KO"))
mycols <- c("#0073C2FF","#EFC000FF","#868686FF","#CD534CFF")
ggplot(a1, aes(x="",y= ratio, fill=Var2))+
  geom_bar(width = 1, stat = "identity",color="white")+
  coord_polar("y",start = 0)+
  geom_text(aes(y=lab.pos,label=ratio1),nudge_x = 0.6, size=7)+scale_fill_manual(values = mycols)+
  facet_grid(~Var1)+theme_void()+theme(strip.text.x  = element_text(size=20),
                                       legend.title = element_blank(),
                                       legend.text = element_text(size=15))
ggsave("ratio_subgroup.png",width = 12)
ggsave("ratio_subgroup.pdf",width = 12)


# --------------------------------------------------------

matrix <- as.data.frame(scRNA[["RNA"]]@counts)["T",]
mat.t <- t(matrix)
phe <- as.data.frame(scRNA@meta.data)
identical(rownames(mat.t),rownames(phe)) 
phe1 <- cbind(phe,mat.t)
fivenum(phe1$T)

phe1$T.pos <- ifelse(phe1$T>0,"T+","T-")
phe1$orig.ident <- factor(phe1$orig.ident, levels = c("WT","KO"))
a4 <- phe1 %>% group_by(orig.ident) %>% summarise(T_pos=mean(T.pos=="T+"),
                                                  T_neg=mean(T.pos=="T-"))
library(reshape2)
a5 <- melt(a4)
colnames(a5)[1:2] <- c("Group","T.pos")
a5$T.pos <- ifelse(a5$T.pos=="T_pos","T+","T-")
a5$label <- paste0(round(a5$value*100,2),"%")
library(ggpubr)
library(ggthemes)
ggbarplot(a5,
          x="Group",
          y="value",
          fill="T.pos",
          palette =  c("#0073C2FF","#EFC000FF"),
          ylab = "",xlab = "",
          label = a5$label,lab.pos = "in",lab.col = "white",lab.size = 5)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_base()+
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 15))
ggsave("T_pos_ratio.png")
ggsave("T_pos_ratio.pdf")


#  ------------------------------------------------------

epi <- subset(scRNA, subset = celltype=="Epiblast")
epi@meta.data
table(epi$orig.ident)

marker_epi <- FindMarkers(epi,ident.1 = "KO", group.by = "orig.ident")
write.csv(marker_epi,file = "epi_DEG_KO_vs_WT.csv")



# -----------------------------------------

epi.dat <- as.data.frame(epi[["RNA"]]@counts)
epi.phe <- as.data.frame(epi@meta.data)
epi.phe1 <- epi.phe %>% arrange(desc(orig.ident))
epi.dat1 <- epi.dat[,rownames(epi.phe1)]



library(edgeR)
epi.dat2 <- cpm(epi.dat1)
epi.merge <- cbind(epi.phe1,t(epi.dat2))


# ----------------------------------------------------------


epi$orig.ident <- factor(epi$orig.ident,levels = c("WT","KO"))

VlnPlot(epi,c("Lef1","Axin2","Wnt3","Wnt8a","Fzd10","Sp5","Tcf7","Ctnnb1","Tcf7l1"), group.by = "orig.ident",pt.size = 0,cols = c("limegreen", "navy"))


singlecell_gene_test <- function(SerautObj, 
                                 genes.use, 
                                 group.by=NULL, 
                                 assay = "RNA", 
                                 comp = NULL, 
                                 alpha_start = .05, 
                                 Bonferroni = T,
                                 only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
  else{
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname]
      group1_exp <- group1_exp[which(group1_exp>0)] 
      
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      group2_exp <- group2_exp[which(group2_exp>0)] 
      
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
    
  }
  
  if (Bonferroni == T){
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
    
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
    
  }
  
  else{
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
  }
  
  return(dfOUT)
}



A <- singlecell_gene_test(epi, 
                          genes.use = c("Wnt3","Wnt8a","Fzd10","Ctnnb1","Tcf7","Lef1","Sp5","Axin2","Tcf7l1"),
                          group.by = 'orig.ident', 
                          comp = c("WT", "KO"))

A1 <- singlecell_gene_test(epi,
                           genes.use = c("Wnt3","Wnt8a","Fzd10","Ctnnb1","Tcf7","Lef1","Sp5","Axin2","Tcf7l1"),
                           group.by = 'orig.ident', 
                           comp = c("WT", "KO"),
                           only_postive = T)

anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(epi, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "orig.ident",
                         features = c("Wnt3","Wnt8a","Fzd10","Ctnnb1","Tcf7","Lef1","Sp5","Axin2","Tcf7l1"), 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("WT","KO"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)
ggsave("wnt_can_vln.png")
ggsave("wnt_can_vln.pdf")



# --------------------------------------------------------------------

library(tidyverse)
library(data.table)
all_df <- cbind(Name = rownames(epi.dat1), DESCRIPTION = "na", epi.dat1[,1:ncol(epi.dat1)])
a <- c("#1.2",rep("",ncol(all_df)-1))
b <- c(nrow(all_df),ncol(all_df)-2,rep("",ncol(all_df)-2))
tmp <- rbind(a,b, colnames(all_df),all_df)

group <- group_list
group <- paste(group, collapse = " ")
group <- c(paste(c(length(group_list), 2, 1), collapse = " "), "# WT KO", group)

write.table(file = "mat.gct", tmp, sep = "\t", col.names = F, row.names = F, quote = F)
write.table(file = "group.cls", group, col.names = F, row.names = F, quote = F)



library(data.table)
library(ggthemes)
gsea1 <- fread('F:/tmp1/epi/gsea/my_analysis.Gsea.1688695903191_development/gsea_report_for_KO_1688695903191.tsv')
gsea1 <- gsea1 %>% filter(`NOM p-val`<0.05) %>% arrange(desc(NES))
gsea2 <- gsea1[c(grep('DEVELOP',gsea1$NAME),3,9),]
gsea3 <- gsea2[-c(4:6),]
gsea3 <- gsea3 %>% arrange(desc(NES)) %>% mutate(`-log(pvale)`=-log(`NOM p-val`), sig=ifelse(gsea3$`NOM p-val`<0.05, "Sigificant","Not Significant"))
gsea3$NAME <- str_remove(gsea3$NAME,"GOBP_") %>% str_replace_all("_"," ")
gsea3$NAME <- factor(gsea3$NAME,levels = rev(c("CELL FATE COMMITMENT","GASTRULATION","MESENDODERM DEVELOPMENT","MESODERM DEVELOPMENT","ENDODERM DEVELOPMENT","ECTODERM DEVELOPMENT")))
ggbarplot(gsea3,x="NAME",y="NES", fill = "sig",
          palette = c("grey", "#FC4E07"),
          #sorting = "descending",
          #label = round(gsea3$NES,1),
          dot.size = 8,
          font.label = list(color = "white", size = 11, 
                            vjust = 0.5),               # Adjust label parameters
          ggtheme = theme_pubr(),  
          add = "segments",
          xlab = "Terms",
          rotate = TRUE)+theme(legend.title = element_text(size = 15),
                               legend.text = element_text(size = 12))+labs(fill="")
ggsave("gsea.png")
ggsave("gsea.pdf")


gsea4 <- gsea1[c(grep('SIGN',gsea1$NAME),grep('ERK',gsea1$NAME)),]
gsea5 <- gsea4[!grep('REGU',gsea4$NAME),]

gsea6 <- gsea4[1:10,]
gsea7 <- fread('./epi/gsea/bmp.Gsea.1688956524349/gsea_report_for_KO_1688956524349.tsv')
tmp8 <- rbind(gsea4,gsea7)
tmp8$pvalue <- ifelse(tmp8$`NOM p-val`==0,-log10(0.0001),-log10(tmp8$`NOM p-val`))
tmp8$term <- tmp8$NAME%>% str_remove('GOBP_') %>% str_replace_all('_',' ')
tmp8 <- tmp8[-c(1,3,8,11,15,19,21,23),]
tmp8$group <- ifelse(grepl("NON CANONICAL WNT",tmp8$term),"NON CANONICAL WNT",
                     ifelse(grepl("NODAL",tmp8$term),"NODAL",
                            ifelse(grepl("NOTCH",tmp8$term),"NOTCH",
                                   ifelse(grepl("ERK1",tmp8$term),"FGF",
                                          ifelse(grepl("BMP",tmp8$term),"BMP" ,"CANONICAL WNT")))))
tmp8$group <- factor(tmp8$group, levels = c("CANONICAL WNT","NON CANONICAL WNT","BMP","FGF","NODAL","NOTCH"))
tmp8$label <- ifelse(tmp8$group=="CANONICAL WNT",tmp8$term,"")
ggplot(tmp8, aes(x=pvalue,y=NES))+geom_point(aes(fill=pvalue,size=NES,color=group), shape=21,stroke=2)+
  geom_vline(xintercept = -log10(0.001),lty=2, size=1.5, color='darkgrey')+
  geom_hline(yintercept = 1.5,lty=2, size=1.5, color='darkgrey')+
  ggrepel::geom_label_repel(aes(x=pvalue,y=NES, label=label),segment.color="grey50",size=4,nudge_y = 0.05)+
  scale_fill_gradient2(low = "lightgreen", high = "pink", 
                       limits = c(0, max(tmp8$pvalue)),
                       midpoint = mean(tmp8$pvalue))+
  scale_size(range = c(2,8))+labs(color="Signaling pathways", fill="-log10(P-value)", x="-log10(P-value)")+
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("gsea_pathway.png",width = 9, height = 7)
ggsave("gsea_pathway.pdf",width = 9, height = 7)


epi$orig.ident <- factor(epi$orig.ident,levels = c("WT","KO"))

VlnPlot(epi,c("Pgk1","Ldha","Aldoa"), group.by = "orig.ident",pt.size = 0,cols = c("limegreen", "navy"))


A <- singlecell_gene_test(epi, 
                          genes.use = c("Pgk1","Ldha","Aldoa"),
                          group.by = 'orig.ident', 
                          comp = c("WT", "KO"))

A1 <- singlecell_gene_test(epi,
                           genes.use = c("Pgk1","Ldha","Aldoa"),
                           group.by = 'orig.ident', 
                           comp = c("WT", "KO"),
                           only_postive = T)

anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(epi, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "orig.ident",
                         features = c("Pgk1","Ldha","Aldoa"), 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("WT","KO"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)
patchwork::wrap_plots(plots_violins)
ggsave("../hif1a_vln.png",width = 9, height = 3.5)
ggsave("../hif1a_vln.pdf",width = 9, height = 3.5)


save(scRNA,file = "scRNA20230710.RData")


#dir.create('F:/tmp1/mes')
setwd('F:/tmp1/mes')
mes <- subset(scRNA, subset = celltype=="Mesendoderm")
mes@meta.data

marker_mes <- FindMarkers(mes,ident.1 = "KO", group.by = "orig.ident")
write.csv(marker_mes,file = "mes_DEG_KO_vs_WT.csv")
mes_up <- marker_mes %>% filter(avg_log2FC>0,p_val_adj<0.05)
mes_down <- marker_mes %>% filter(avg_log2FC<0,p_val_adj<0.05)

library(clusterProfiler)
library(org.Mm.eg.db) 
gene_u <- bitr(rownames(mes_up),fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Mm.eg.db)

head(gene_u)
nrow(gene_u)

gene_d <- bitr(rownames(mes_down),fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Mm.eg.db)

head(gene_d)
nrow(gene_d)

ego <- enrichGO(gene= gene_u$ENTREZID,
                keyType = "ENTREZID",OrgDb = org.Mm.eg.db,
                ont= "ALL",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,
                readable      = TRUE)

write.csv(ego,file="GO_enrichment.csv")
save(ego,file = "ego_up.Rdata")

ego <- enrichGO(gene= gene_d$ENTREZID,
                keyType = "ENTREZID",OrgDb = org.Mm.eg.db,
                ont= "ALL",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,
                readable      = TRUE)

write.csv(ego,file="GO_down.csv")
save(ego,file = "ego_down.Rdata")



#  ------------------------------------------

epi.dat <- as.data.frame(epi[["RNA"]]@counts)
epi.phe <- as.data.frame(epi@meta.data)
epi.phe1 <- epi.phe %>% arrange(desc(orig.ident))
epi.dat1 <- epi.dat[,rownames(epi.phe1)]
identical(colnames(epi.dat1),rownames(epi.phe1))

table(epi.phe1$orig.ident) #总共有WT 7779个细胞，KO 5190个细胞


library(edgeR)
epi.dat2 <- cpm(epi.dat1)
epi.merge <- cbind(epi.phe1,t(epi.dat2))


library(ggpubr)
library(ggthemes)
library(ggsci)


data <- read.csv("GO_down1.csv")
if(F){
  library(tidyverse)
  data_BP <- data%>% filter(ONTOLOGY=="BP")%>% arrange(desc(Count))
  data_CC <- data%>% filter(ONTOLOGY=="CC")%>% arrange(desc(Count))
  data_MF <- data%>% filter(ONTOLOGY=="MF")%>% arrange(desc(Count))
  data_BP <- data_BP[1:10,]
  data_CC <- data_CC[1:10,]
  data_MF <- data_MF[1:10,]
  data <- na.omit(rbind(data_BP,data_CC, data_MF))
}





a <- read.csv("GO_enrichment1.csv")
data$Count <- -data$Count
data <- rbind(a,data)
data <- na.omit(data)
data$group <- ifelse(data$Count>0,"Upregulated","Downregulated")
data$P_log = -log10(data$p.adjust)

data$Description<- sapply(as.character(data$Description),function(string) {ifelse (nchar(string)>40, paste(substr(string,1,40),"...",sep=""),string)})

p1 <- ggbarplot(data,x="Description", y="Count", fill = "P_log", rotate=T,
                sort.val = "asc", 
                xlab = "Terms")+
  scale_y_continuous(expand = c(0,0),limits = c(-13,21))+
  scale_fill_continuous(low="blue", high="red", space='rgb')+labs(fill="-log10(P-value)")+
  geom_hline(yintercept = 0,size=1)+
  theme_pubr()+
  theme(axis.text = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size=20),
        axis.title = element_text(size=20));p1
ggsave("GO_MERGE_Bar.pdf",p1, width =12,height = 10, units = "in")


# -------------------------------------------------------------------------

load('./scRNA20230710.RData')
library(monocle3)
library(Seurat)
expression <- GetAssayData(scRNA, assay = 'RNA', slot = 'counts')

cell_metadata <- scRNA@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression))
rownames(gene_annotation) <- rownames(expression)
cds <- new_cell_data_set(expression, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
head(colData(cds))
fData(cds)
rownames(fData(cds))[1:10]
head(counts(cds))


# Retrieve clustering information from Surat object -----------------------


recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- scRNA@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- scRNA@reductions$umap@cell.embeddings
library(ggplot2)
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj


# Learn Trajectory --------------------------------------------------------

cds <- learn_graph(cds, use_partition=F, close_loop = T, learn_graph_control = list(minimal_branch_len=5))
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = T, label_leaves = T,graph_label_size = 3,
           group_label_size = 5)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = T, label_leaves = T,graph_label_size = 2,
           group_label_size = 5)
ggsave('./traj_UMAP_cluster.png')
ggsave('./traj_UMAP_cluster.pdf')
plot_cells(cds, color_cells_by = "pseudotime")
ggsave('./traj_pseu.png')
ggsave('./traj_pseu.pdf')


# --------------------------------------------------------------

save(cds, file = "./monocle3_cds.RData")


# ------------------------------------------------------------
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
library(ggpubr)
ggplot(data.pseudo, aes(monocle3_pseudotime, celltype, fill = orig.ident)) + geom_boxplot()
data.pseudo$celltype <- factor(data.pseudo$celltype, levels = c("Mesendoderm","Ecotderm","Epiblast","Undiff. ES"))
p1 <- ggboxplot(data.pseudo,
                x="celltype",
                y="monocle3_pseudotime",
                fill = "orig.ident",
                palette = "npg")+stat_compare_means(aes(group=celltype),
                                                    label = "p.signif")+
  labs(x="",y="Pseudotime",fill="")+
  coord_flip()
ggsave('./boxplot_pseu_celltype.png')
ggsave('./boxplot_pseu_celltype.pdf')


# ------------------------------------------------------

scRNA$pseudotime <- pseudotime(cds)
FeaturePlot(scRNA, "pseudotime")
#scRNA$orig.ident <- factor(scRNA$orig.ident, levels = c("Mesendoderm","Ecotderm","Epiblast","Undiff. ES"))
p <- RidgePlot(scRNA, "pseudotime" )
p
tmp1 <- p$data
tmp1$ident <- factor(tmp1$ident, levels = c("Mesendoderm","Ecotderm","Epiblast","Undiff. ES"))
library(ggridges)
ggplot(tmp1, aes_string("pseudotime","ident"))+
  geom_density_ridges_gradient(aes(fill=ident))+theme_ridges()
ggsave("celltype_ridgeplot_all.png")
ggsave("celltype_ridgeplot_all.pdf")
identical(rownames(tmp1),rownames(scRNA@meta.data))
#TRUE
tmp1$group <- scRNA@meta.data$orig.ident

p2 <- ggplot(tmp1, aes_string("pseudotime","ident"))+
  geom_density_ridges_gradient(aes(fill=group))+theme_ridges();p2
ggsave("celltype_ridgeplot.png")
ggsave("celltype_ridgeplot.pdf")

p3 <- p1+theme_void()
p4 <- p2+NoLegend()
library(patchwork)
p4+p3+plot_layout(widths = c(2,1))+scale_fill_npg()
ggsave("celltype_ridgeplot_merge.png",width = 9)
ggsave("celltype_ridgeplot_merge.pdf",width = 9)



# -------------------------------------------------------------
pluripotent=c("Pou5f1","Nanog","Sox2","Zfp42")
Mesendoderm=c("T","Mixl1","Mesp1","Eomes","Gsc")
Ecotderm=c("Btg2","Crip2","Cd24a","Nkx2-9","Snai3")



marker_pseudo <- function(gene.select){
  my_genes <- row.names(subset(fData(cds), gene_short_name %in% gene.select)) 
  cds_subset <- cds[my_genes,]
  plot_genes_in_pseudotime(cds_subset, color_cells_by = "celltype" )
  ggsave(paste0(substitute(gene.select),"_marker_pseduo.png"))
  ggsave(paste0(substitute(gene.select),"_marker_pseduo.pdf"))
}


marker_pseudo(pluripotent)
marker_pseudo(Mesendoderm)
marker_pseudo(Ecotderm)

