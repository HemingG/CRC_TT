library(Seurat)
library(data.table)

A1_004<-Load10X_Spatial('/data/20230713-V42L22-004-A1-1492164-1-outs/')
A1_310<-Load10X_Spatial('/data/20230713-V42L22-310-A1-1476353-2-outs/')
D1_004<-Load10X_Spatial('/data/20230713-V42L22-004-D1-1521037-5-outs/')
D1_310<-Load10X_Spatial('/data/20230713-V42L22-310-D1-1524789-2-outs/')



plot1 <- VlnPlot(A1_004, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(A1_004, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##################空转标准分析
##################空转标准分析

setwd('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-9-27/01.spatial_plot')

A1_004@project.name<-'A1_004'

spatial_fun<-function(A1){
  plot1 <- VlnPlot(A1_004, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(A1_004, features = "nCount_Spatial") + theme(legend.position = "right")
  p3<-wrap_plots(plot1, plot2)

  dir.create(A1@project.name)
  setwd(A1@project.name)
  ggsave(p3,filename = paste0("QC.pdf"),width = 15,height = 8)
  A1 <- SCTransform(A1, assay = "Spatial", verbose = FALSE)

  A1 <- RunPCA(A1, assay = "SCT", verbose = FALSE)

  A1 <- FindNeighbors(A1, reduction = "pca", dims = 1:15)
  A1 <- FindClusters(A1, verbose = FALSE) # resolution = 0.2-1.2

  A1 <- RunUMAP(A1, reduction = "pca", dims = 1:15)
  setwd('../')
  A1

}


plot1 <- VlnPlot(D1_310, features = "nCount_Spatial", pt.size = 0.1,group.by = 'orig.ident') + NoLegend()
plot2 <- SpatialFeaturePlot(D1_310, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


ElbowPlot(A1_004_1)
A1_310@project.name<-"A1_310"
A1_004_1<-spatial_fun(A1_004)
A1_310_1<-spatial_fun(A1_310)
ElbowPlot(A1_310_1)
D1_004@project.name<-'D1_004'
D1_004_1<-spatial_fun(D1_004)
ElbowPlot(D1_004_1)

D1_310@project.name<-'D1_310'
D1_310_1<-spatial_fun(D1_310)
ElbowPlot(D1_310_1)

plot1 <- VlnPlot(D1_310_1, features = "nCount_Spatial", pt.size = 0.1,group.by = 'orig.ident') + NoLegend()
plot2 <- SpatialFeaturePlot(D1_310_1, features = "nCount_Spatial") + theme(legend.position = "right")
p3<-wrap_plots(plot1, plot2)


#######################################################region mapping
A1_310_coor<-fread('/share/home/guests01/bhzhang/project/peiqian/cloupe/1476353.csv',header = T)
A1_004_coor<-fread('/share/home/guests01/bhzhang/project/peiqian/cloupe/1492164.csv',header = T)
D1_004_coor<-fread('/share/home/guests01/bhzhang/project/peiqian/cloupe/1521037.csv',header = T)
D1_310_coor<-fread('/share/home/guests01/bhzhang/project/peiqian/cloupe/1524789.csv',header = T)


A1_004_1$type<-mapvalues(colnames(A1_004_1),A1_004_coor$Barcode,A1_004_coor$`1492164`)
A1_004_1$type[grepl("-1",A1_004_1$type)]="other"
A1_310_1$type<-mapvalues(colnames(A1_310_1),A1_310_coor$Barcode,A1_310_coor$`1476353`)
A1_310_1$type[grepl("-1",A1_310_1$type)]="other"
D1_004_1$type<-mapvalues(colnames(D1_004_1),D1_004_coor$Barcode,D1_004_coor$`1521037`)
D1_004_1$type[grepl("-1",D1_004_1$type)]="other"
D1_310_1$type<-mapvalues(colnames(D1_310_1),D1_310_coor$Barcode,D1_310_coor$`1524789`)
D1_310_1$type[grepl("-1",D1_310_1$type)]="other"


#####################################################diff

A1_004_1_c_vs_n<-subset(FindMarkers(A1_004_1,group.by = 'type',ident.1 = 'Carcinothrombosis','Tumor_Noncarcinothrombosis'),p_val<0.05)
fwrite(A1_004_1_c_vs_n,file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2023-9-27/01.spatial_plot/A1_004/Carcinothrombosis_vs_Tumor_Noncarcinothrombosis.xls',sep = '\t',col.names = T,row.names = T)
A1_310_1_c_vs_n<-subset(FindMarkers(A1_310_1,group.by = 'type',ident.1 = 'Carcinothrombosis','Tumor_Noncarcinothrombosis'),p_val<0.05)
fwrite(A1_310_1_c_vs_n,file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2023-9-27/01.spatial_plot/A1_310/Carcinothrombosis_vs_Tumor_Noncarcinothrombosis.xls',sep = '\t',col.names = T,row.names = T)
D1_004_1_c_vs_n<-subset(FindMarkers(D1_004_1,group.by = 'type',ident.1 = 'Carcinothrombosis','Tumor_Noncarcinothrombosis'),p_val<0.05)
fwrite(D1_004_1_c_vs_n,file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2023-9-27/01.spatial_plot/D1_004/Carcinothrombosis_vs_Tumor_Noncarcinothrombosis.xls',sep = '\t',col.names = T,row.names = T)
D1_310_1_c_vs_n<-subset(FindMarkers(D1_310_1,group.by = 'type',ident.1 = 'Carcinothrombosis','Tumor_Noncarcinothrombosis'),p_val<0.05)
fwrite(D1_310_1_c_vs_n,file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2023-9-27/01.spatial_plot/D1_310/Carcinothrombosis_vs_Tumor_Noncarcinothrombosis.xls',sep = '\t',col.names = T,row.names = T)


##########################################################差异富集分析
diff_list=list(A1_004=A1_004_1_c_vs_n,A1_310=A1_310_1_c_vs_n,D1_004=D1_004_1_c_vs_n,
               D1_310=D1_310_1_c_vs_n
               )

dim(D1_004_1_c_vs_n)

setwd('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-9-27/01.spatial_plot')
for(nm in names(diff_list)){
  print(nm)
  setwd(nm)
  try({
    print(head(rownames(diff_list[[nm]])))

  ego <- enrichGO(gene = rownames(diff_list[[nm]]),
                  keyType = "SYMBOL",
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 1,
                  readable = F) ###改成F的keyType才能使用symbol
  head(ego@result)
  write.table(ego, file=paste0("GO.xls"), sep="\t", quote=F, row.names=F)
  p1<-barplot(ego,showCategory = 10,color='p.adjust',split="ONTOLOGY") +facet_grid(ONTOLOGY~.,scales = 'free')
  ggsave(p1,filename = paste0('GO.pdf'),width = 20,height = 20)
  ggsave(p1,filename = paste0('GO.png'),width = 20,height = 20)
  })
  try({
    # kegg

    common_gene_id<-bitr(rownames(diff_list[[nm]]), fromType = c('SYMBOL'), toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')


    enrich_kegg <- enrichKEGG(gene = common_gene_id$ENTREZID,organism = "hsa", pvalueCutoff =0.05,use_internal_data = T)

    write.table(enrich_kegg, file=paste0("kegg.txt"), sep="\t", quote=F, row.names=F)
    p1<-barplot(enrich_kegg ,color = "p.adjust", showCategory = 20, size = NULL, split = NULL,font.size = 12, title = "")
    ggsave(p1,filename = paste0('KEGG.pdf'),width = 20,height = 13)
    ggsave(p1,filename = paste0('KEGG.png'),width = 20,height = 13)
  })

  setwd('../')
}

####################################################TF

###################################################提取 矩阵 出来进行整合
A1_004_1_seurat <-CreateSeuratObject(A1_004_1@assays$Spatial@counts,project = 'A1-004')
A1_310_1_seurat <-CreateSeuratObject(A1_310_1@assays$Spatial@counts,project = 'A1-310')
D1_004_1_seurat <-CreateSeuratObject(D1_004_1@assays$Spatial@counts,project = 'D1-004')
D1_310_1_seurat <-CreateSeuratObject(D1_310_1@assays$Spatial@counts,project = 'D1-310')

###################################################合并
merge1<-merge(x=A1_004_1_seurat,y=c(A1_310_1_seurat,D1_004_1_seurat,D1_310_1_seurat),add.cell.ids=c('A1-004',"A1-310","D1-004","D1-310"))
# saveRDS(merge1,file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2023-10-11/data/merge1.rds')
###################################################
###################################################

# 先整合
inte_data<-readRDS('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-10-11/01.marker/result/backup/g_combined.rds')


setwd('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-10-11/02.transcript_factor')
# fwrite(data.frame(inte_data@assays$SCT@data,check.names = F),file = 'mat.xls',sep = ',',col.names = T,row.names = T)
# fwrite(data.frame(row.names=rownames(inte_data@meta.data),cluster=inte_data$type,check.names = F),file = 'meta.xls',sep = ',',col.names = T,row.names = T)


########################################################cellchat

# rm(data.input)
data.input = GetAssayData(A1_004_1, slot = "data", assay = "SCT") # normalized data matrix
meta = data.frame(labels = A1_004_1$type, row.names = names(Idents(A1_004_1))) # manually create a dataframe consisting of the cell labels
spatial.locs = GetTissueCoordinates(A1_004_1, scale = NULL, cols = c("imagerow", "imagecol"))
scale_factor <-jsonlite::fromJSON('/spatial/scalefactors_json.json')
View(scale_factor)
scale_factor = list(spot.diameter = 65, spot = scale_factor$spot_diameter_fullres, # these two information are required
                     fiducial = scale_factor$fiducial_diameter_fullres, hires = scale_factor$tissue_hires_scalef, lowres = scale_factor$tissue_lowres_scalef # these three information are not required
)


cellchat_A1_004 <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs,scale.factors = scale_factor)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_A1_004@DB <- CellChatDB.use
cellchat <- subsetData(cellchat_A1_004)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- c("CXCL")
# VISFATIN
pathways.show <- c("VISFATIN")
# CSF
pathways.show <- c("CSF")
cellchat@netP$pathways
grep("CXCL",cellchat@netP$pathways,value = T)
# Circle plot
par(mfrow=c(1,1))
p1<-netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
p1

#############################
# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

############################
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

############################

par(mfrow=c(1,1))
p1<-netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 3, vertex.label.cex = 3.5)
print(p1)


######################################################merge cellchat




setwd('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-10-16/01.cellchat/')
A1_004_1<-readRDS('../data/A1_004_1_cellchat.rds')
A1_310_1<-readRDS('../data/A1_310_1_cellchat.rds')
D1_004_1<-readRDS('../data/D1_004_1_cellchat.rds')
D1_310_1<-readRDS('../data/D1_310_1_cellchat.rds')
A1_004_1 <- updateCellChat(A1_004_1)
A1_310_1 <- updateCellChat(A1_310_1)
D1_004_1 <- updateCellChat(D1_004_1)
D1_310_1 <- updateCellChat(D1_310_1)

object.list <- list(A1_004_1 = A1_004_1, A1_310_1 = A1_310_1,D1_004_1=D1_004_1,D1_310_1=D1_310_1)

cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = T)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

levels(cellchat@idents$A1_004_1)

?netAnalysis_signalingChanges_scatter

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Carcinothrombosis",comparison = c(3, 4))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tumor_Noncarcinothrombosis",comparison = c(3, 4))
patchwork::wrap_plots(plots = list(gg1,gg2))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional",comparison=c(1,2))
cellchat <- netEmbedding(cellchat, type = "functional",comparison = c(1,2),umap.method=c('uwot'))
cellchat <- netClustering(cellchat, type = "functional",comparison = c(1,2))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,comparison = c(1,2))



rankSimilarity(cellchat, type = "functional",comparison1 = c(1,2))
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional",comparison=c(3,4))
cellchat <- netEmbedding(cellchat, type = "functional",comparison = c(3,4),umap.method=c('uwot'))
cellchat <- netClustering(cellchat, type = "functional",comparison = c(3,4))

rankSimilarity(cellchat, type = "functional",comparison1 = c(3,4))


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(3, 4))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(3, 4))
gg1 + gg2


netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:3),  comparison = c(1, 2,3,4), angle.x = 45)



gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:3),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in A1_004", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:3),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in A1_004", angle.x = 45, remove.isolate = T)
gg1 + gg2


gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:3),  comparison = c(3, 4), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:3),  comparison = c(3, 4), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
gg1 + gg2

##################################################
##################################################TF  rss

# 重新绘制rss 热图

matrice<-t(data.frame(fread('/02.transcript_factor/outdir/PySCENIC/00.data/spatial_total_rss.csv'),row.names = 1,check.names = F))

temp_matrix=matrice[,order(colnames(matrice))]
pheatmap::pheatmap(matrice, border = NA,
              angle_col = 45, cluster_cols = F,
              fontsize_row = 3)
# ?plotRSS
rssPlot <-plotRSS(matrice, zThreshold = 0.02)

auc1="/02.transcript_factor/outdir/PySCENIC/00.data/spatial_auc_mtx.csv"
meta="/02.transcript_factor/meta.xls"
auc1 <- t(data.frame(fread(auc1,header = T),check.names = F,row.names = 1))
anno <- data.frame(fread(meta,header = T),check.names = F,row.names = 1)

anno <- anno[order(anno$cluster),,drop = FALSE]
auc1 <- auc1[,order(match(colnames(auc1),rownames(anno)))]
pheatmap::pheatmap(auc1, show_colnames = F,
                        annotation_col = anno,
                        treeheight_row = 20, treeheight_col = 20,
                        cluster_cols = F,show_rownames = T,
                        color = colorRampPalette(c("blue","white","red"))(100),
                        breaks = seq(-3, 3, length.out = 100),
                        fontsize_row = 5,
                        border_color = NA,
                        cutree_rows = 4)


clusters_means_auc <- sapply(split(colnames(auc1),anno$cluster),function(cells){
  rowMeans(auc[,cells])
})
clusters_means_auc <- scale(clusters_means_auc, center = T, scale = T)
p1<-pheatmap::pheatmap(clusters_means_auc, fontsize_row = 5,
                       cutree_rows = 4,
                       color=colorRampPalette(c("blue","white","red"))(100),
                       breaks=seq(-3, 3, length.out = 100),
                       treeheight_row = 15,
                       treeheight_col = 15,
                       border_color = NA,
                       angle_col = 45)


######################################################################### deconvolution singlecell preprocess

anno<-fread("/data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
mat<-fread("/data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",data.table = F)

rownames(mat)<-mat$Index
mat$Index<-NULL

GSE132465<-CreateSeuratObject(mat,project = 'GSE132465')
GSE132465[["percent.mt"]] <- PercentageFeatureSet(GSE132465, pattern = "^MT-")

GSE132465$sample<-"sample"

VlnPlot(GSE132465,group.by = 'sample', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

GSE132465$anno <- mapvalues(colnames(GSE132465),anno$Index,anno$Cell_type)

# integrated
GSE132465_1<-harmony_func(GSE132465)

Idents(GSE132465_1)<-GSE132465_1$anno
levels(GSE132465_1)<-names(sort(table(Idents(GSE132465_1)),decreasing = T))


###############################################################################  spatial deconvolution


plot1 <- VlnPlot(A1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(A1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


# TT_1  D1_310_1524789  D1_310
# TT_2  A1_004_1492164  A1_004
# TT_3  A1_310_1476353  A1_310
# TT_4  D1_004_1521037  D1_004
# NTT_1     GSM7089855  CR_22_15980_FP_R_JSV_1
# NTT_2     GSM7089857  CR_22_07526_TS_R_SSV_1
# NTT_3     GSM7089856  CR_22_14422_FP_R_JSV_1
# NTT_4     GSM7089858  CR_22_15979_FP_R_JSV_1
# GSM7089857	P3	CR_22_07526_TS_R_SSV_1
# GSM7089856	P2	CR_22_14422_FP_R_JSV_1
# GSM7089858	P4	CR_22_15979_FP_R_JSV_1
# GSM7089855	P1	CR_22_15980_FP_R_JSV_1
A1_004<-Load10X_Spatial('/data/20230713-V42L22-004-A1-1492164-1-outs/')
A1_310<-Load10X_Spatial('/data/20230713-V42L22-310-A1-1476353-2-outs/')
D1_004<-Load10X_Spatial('/data/20230713-V42L22-004-D1-1521037-5-outs/')
D1_310<-Load10X_Spatial('/data/20230713-V42L22-310-D1-1524789-2-outs/')

CR_22_07526_TS_R_SSV_1<-Load10X_Spatial("/data/ncbi_data/CR_22_07526_TS_R_SSV_1/outs/")
CR_22_14422_FP_R_JSV_1<-Load10X_Spatial("/data/ncbi_data/CR_22_14422_FP_R_JSV_1/outs/")
CR_22_15979_FP_R_JSV_1<-Load10X_Spatial("/data/ncbi_data/CR_22_15979_FP_R_JSV_1/outs/")
CR_22_15980_FP_R_JSV_1<-Load10X_Spatial("/data/ncbi_data/CR_22_15980_FP_R_JSV_1/outs/")
SpatialFeaturePlot(CR_22_07526_TS_R_SSV_1, features = "nCount_Spatial")
SpatialFeaturePlot(CR_22_14422_FP_R_JSV_1, features = "nCount_Spatial")
SpatialFeaturePlot(CR_22_15979_FP_R_JSV_1, features = "nCount_Spatial")


merge1<-merge(x=D1_310,y=c(A1_004,A1_310,D1_004,CR_22_15980_FP_R_JSV_1,CR_22_07526_TS_R_SSV_1,CR_22_14422_FP_R_JSV_1,CR_22_15979_FP_R_JSV_1),
               add.cell.ids=c("TT-1","TT-2","TT-3","TT-4","NTT-1","NTT-2","NTT-3","NTT-4"))





######################################  integrated by shell R script in background
######################################  read the integrated rds into rstudio
integrated_spatial<-readRDS('/01.spatial_integrate/result/backup/g_combined.rds')


# 输出QC结果
library(cowplot)
?SpatialFeaturePlot
integrated_spatial$orig.ident<-sapply(strsplit(rownames(integrated_spatial@meta.data),"_"),'[[',1)
integrated_spatial$orig.ident <- gsub("-",".",integrated_spatial$orig.ident)

integrated_spatial$orig.ident <- factor(integrated_spatial$orig.ident,levels = c("TT.1","TT.2","TT.3","TT.4","NTT.1","NTT.2","NTT.3","NTT.4"))

SpatialFeaturePlot(integrated_spatial,features = 'nFeature_Spatial',ncol = 4)


options(scipen = 999)
count =1
# QC
for(samp in unique(integrated_spatial$orig.ident)){
  tmp<-subset(integrated_spatial,orig.ident==samp)
  tmp@images <- tmp@images[count]
  p1 <- SpatialFeaturePlot(tmp,features = 'nFeature_Spatial')
  p2 <- VlnPlot(tmp,features = 'nFeature_Spatial',group.by = 'orig.ident',pt.size = 0) +ggtitle("")

  p1_1 <- SpatialFeaturePlot(tmp,features = 'nCount_Spatial')
  p2_1 <- VlnPlot(tmp,features = 'nCount_Spatial',group.by = 'orig.ident',pt.size = 0) +ggtitle("")

  p1_2 <- SpatialFeaturePlot(tmp,features = 'Mito.percent')
  p2_2 <- VlnPlot(tmp,features = 'Mito.percent',group.by = 'orig.ident',pt.size = 0) +ggtitle("")

  p3<-plot_grid(p2, p1, labels = c('A', 'B'), ncol = 2, rel_widths = c(1, 2))
  p3_1<-plot_grid(p2_1, p1_1, labels = c('A', 'B'), ncol = 2, rel_widths = c(1, 2))
  p3_2<-plot_grid(p2_2, p1_2, labels = c('A', 'B'), ncol = 2, rel_widths = c(1, 2))
  p4<-p3_1|p3|p3_2
  ggsave(p4,filename = paste0(samp,"_QC.pdf"),width = 28,height = 10)
  ggsave(p4,filename = paste0(samp,"_QC.png"),width = 28,height = 10)
  count = count +1
}
gc()



####################################### read the region from 10x cloupe

class(integrated_spatial$orig.ident)
TT_3<-fread("/data/coordinate/cloupe-1476353.csv",data.table=F)
TT_2<-fread("/data/coordinate/cloupe-1492164.csv",data.table=F)
TT_4<-fread("/data/coordinate/cloupe-1521037.csv",data.table=F)
TT_1<-fread("/data/coordinate/cloupe-1524789.csv",data.table=F)
NTT_2<-fread("/data/coordinate/CR_22_07526_TS_R_SSV_1.csv",data.table=F)
NTT_3<-fread("/data/coordinate/CR_22_14422_FP_R_JSV_1.csv",data.table=F)
NTT_4<-fread("/data/coordinate/CR_22_15979_FP_R_JSV_1.csv",data.table=F)
NTT_1<-fread("/data/coordinate/CR_22_15980_FP_R_JSV_1.csv",data.table=F)

########################################################
########################################################

TT_3$Barcode<-paste0("TT-3_",TT_3$Barcode)
colnames(TT_3)<-c('Barcode',"area")
table(TT_3$area)

TT_2$Barcode<-paste0("TT-2_",TT_2$Barcode)
colnames(TT_2)<-c('Barcode',"area")
table(TT_2$area)

TT_1$Barcode<-paste0("TT-1_",TT_1$Barcode)
colnames(TT_1)<-c('Barcode',"area")
table(TT_1$area)

TT_4$Barcode<-paste0("TT-4_",TT_4$Barcode)
colnames(TT_4)<-c('Barcode',"area")
table(TT_4$area)
####################

NTT_2
NTT_2$Barcode<-paste0("NTT-2_",NTT_2$Barcode)
colnames(NTT_2)<-c('Barcode',"area")
table(NTT_2$area)


NTT_1$Barcode<-paste0("NTT-1_",NTT_1$Barcode)
colnames(NTT_1)<-c('Barcode',"area")
table(NTT_1$area)


NTT_3$Barcode<-paste0("NTT-3_",NTT_3$Barcode)
colnames(NTT_3)<-c('Barcode',"area")
table(NTT_3$area)

NTT_4$Barcode<-paste0("NTT-4_",NTT_4$Barcode)
colnames(NTT_4)<-c('Barcode',"area")

all_df<-do.call(rbind,list(TT_1,TT_2,TT_3,TT_4,NTT_1,NTT_2,NTT_3,NTT_4))



integrated_spatial$area <- mapvalues(rownames(integrated_spatial@meta.data),all_df$Barcode,all_df$area)





##################################################################  inferCNV  cancer recognize
# running in the background  in shell R script
setwd('/data/02.spatial_CNV')
df<-data.frame(integrated_spatial@assays$Spatial@counts,check.names = F)


anno<-data.frame(barcode=colnames(integrated_spatial),anno=integrated_spatial$area,check.rows = F)
fwrite(anno,file = 'anno.tsv',sep = '\t',col.names = F,row.names = F)

setwd("/data/02.spatial_CNV/no_other")
# 没有other的
remove_other<-subset(integrated_spatial,idents = 'other',invert = T)
df<-data.frame(remove_other@assays$Spatial@counts,check.names = F)
fwrite(df,file="count.tsv",sep="\t",col.names = T,row.names = T)

anno<-data.frame(barcode=colnames(remove_other),anno=remove_other$area,check.rows = F)
fwrite(anno,file = 'anno.tsv',sep = '\t',col.names = F,row.names = F)


################################################################################ deconvolution single cell preprocess
singlecell_anno<-readRDS('/data/singlecell_anno.rds')
# Macrophage, CD4 T cells, CD8 T cells, B cells, Endothelial cells, Fibroblast cells, SMC, Mast cells, cDC


DimPlot(adjust_annotation_singlecell,label = T,group.by = 'subtype1')

meta_data<-fread('/data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt',data.table = F)

#############
setwd('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-12-4/02.spatial_CNV/single_cell_cnv')

df<-data.frame(Epithelial_cells@assays$RNA@counts,check.names = F)
fwrite(df,file="count.tsv",sep="\t",col.names = T,row.names = T)

anno<-data.frame(barcode=colnames(Epithelial_cells),anno= paste0(Epithelial_cells$sample,"_",Epithelial_cells$subtype1),check.rows = F)
View(anno)
table(anno$anno)

anno$anno<-gsub(".*Normal_Epithelial_cells","Normal_Epithelial_cells",anno$anno)

fwrite(anno,file = 'anno.tsv',sep = '\t',col.names = F,row.names = F)
gc()


#####################################################
##################################################### QC

DimPlot(integrated_spatial,reduction = 'tsne')

options(scipen = 999)
count =1
for(samp in unique(integrated_spatial$orig.ident)){
  tmp<-subset(integrated_spatial,orig.ident==samp)
  tmp@images <- tmp@images[count]
  p1 <- SpatialFeaturePlot(tmp,features = 'nFeature_Spatial')
  p2 <- VlnPlot(tmp,features = 'nFeature_Spatial',group.by = 'orig.ident',pt.size = 0) +ggtitle("")

  p1_1 <- SpatialFeaturePlot(tmp,features = 'nCount_Spatial')
  p2_1 <- VlnPlot(tmp,features = 'nCount_Spatial',group.by = 'orig.ident',pt.size = 0) +ggtitle("")

  p1_2 <- SpatialFeaturePlot(tmp,features = 'Mito.percent')
  p2_2 <- VlnPlot(tmp,features = 'Mito.percent',group.by = 'orig.ident',pt.size = 0) +ggtitle("")

  p3<-plot_grid(p2, p1, labels = c('A', 'B'), ncol = 2, rel_widths = c(1, 2))
  p3_1<-plot_grid(p2_1, p1_1, labels = c('A', 'B'), ncol = 2, rel_widths = c(1, 2))
  p3_2<-plot_grid(p2_2, p1_2, labels = c('A', 'B'), ncol = 2, rel_widths = c(1, 2))
  p4<-p3_1|p3|p3_2
  ggsave(p4,filename = paste0(samp,"_QC.pdf"),width = 28,height = 10)
  ggsave(p4,filename = paste0(samp,"_QC.png"),width = 28,height = 10)
  count = count +1
}


count =1
for(samp in unique(integrated_spatial$orig.ident)){
  tmp<-subset(integrated_spatial,orig.ident==samp)
  tmp@images <- tmp@images[count]
  p1 <- SpatialFeaturePlot(tmp,features = 'nFeature_Spatial')
  p2 <- VlnPlot(tmp,features = 'nFeature_Spatial',group.by = 'orig.ident',pt.size = 0) +ggtitle("")+NoLegend()

  p1_1 <- SpatialFeaturePlot(tmp,features = 'nCount_Spatial')
  p2_1 <- VlnPlot(tmp,features = 'nCount_Spatial',group.by = 'orig.ident',pt.size = 0) +ggtitle("")+NoLegend()

  p1_2 <- SpatialFeaturePlot(tmp,features = 'Mito.percent')
  p2_2 <- VlnPlot(tmp,features = 'Mito.percent',group.by = 'orig.ident',pt.size = 0) +ggtitle("")+NoLegend()

  ggsave(p1,filename = paste0(samp,"_spatial_Feature_QC.pdf"),width = 10,height = 10)
  ggsave(p1,filename = paste0(samp,"_spatial_Feature_QC.png"),width = 10,height = 10)
  ggsave(p2,filename = paste0(samp,"_Feature_QC.pdf"),width = 10,height = 10)
  ggsave(p2,filename = paste0(samp,"_Feature_QC.png"),width = 10,height = 10)

  ggsave(p1_1,filename = paste0(samp,"_spatial_nCount_QC.pdf"),width = 10,height = 10)
  ggsave(p1_1,filename = paste0(samp,"_spatial_nCount_QC.png"),width = 10,height = 10)
  ggsave(p2_1,filename = paste0(samp,"_nCount_QC.pdf"),width = 10,height = 10)
  ggsave(p2_1,filename = paste0(samp,"_nCount_QC.png"),width = 10,height = 10)

  ggsave(p1_2,filename = paste0(samp,"_spatial_Mito.percent_QC.pdf"),width = 10,height = 10)
  ggsave(p1_2,filename = paste0(samp,"_spatial_Mito.percent_QC.png"),width = 10,height = 10)
  ggsave(p2_2,filename = paste0(samp,"_Mito.percent_QC.pdf"),width = 10,height = 10)
  ggsave(p2_2,filename = paste0(samp,"_Mito.percent_QC.png"),width = 10,height = 10)
  count = count +1
}


########################################################### infercnv result process read from the shell background 

single_cnv<-readRDS('/01.RunCNV/run.final.infercnv_obj')


tree1 <- read.tree(file = "/01.RunCNV/infercnv.preliminary.observations_dendrogram.txt")

plot(tree1, edge.color="blue", edge.width=2, font=2)

plot(tree1, edge.color="blue")
# show.tip.label = FALSE
plot(tree1, edge.color="blue",show.tip.label = FALSE)
hc <- as.hclust(tree1)
k <- 6  
clusters <- cutree(hc, k)
cluster_labels <- paste0("Cluster", 1:k)

dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = k, groupLabels = cluster_labels)
plot(dend)

############################################################################ cluster4 was cancer
df <- data.frame(clusters)
non_maligant_cluster <- subset(df,clusters=='4')

Epithelial_cells$maligan_type <-mapvalues(rownames(Epithelial_cells@meta.data),rownames(df),df$clusters)

DimPlot(Epithelial_cells,group.by = "maligan_type")
table(Epithelial_cells$maligan_type)
Epithelial_cells$maligan_type <- gsub(".*-N_.*","",Epithelial_cells$maligan_type)
Epithelial_cells$maligan_type[Epithelial_cells$maligan_type==""]<-'Normal_epithelial'
Epithelial_cells$maligan_type[Epithelial_cells$maligan_type=="Non_malignant"]<-'Normal_epithelial'
Epithelial_cells$maligan_type[Epithelial_cells$maligan_type=="4"]<-'Non_malignant_epithelial'
Epithelial_cells$maligan_type[Epithelial_cells$maligan_type%in%c("1","2","3","5","6")]<-'Malignant_epithelial'

###########################################################单细胞处理

maligant_Epithelial_cells<-subset(Epithelial_cells,maligan_type=="Normal_epithelial",invert = T)
maligant_Epithelial_cells$anno1<-maligant_Epithelial_cells$maligan_type
table(maligant_Epithelial_cells$anno1)
maligant_Epithelial_cells$anno1<-mapvalues(rownames(maligant_Epithelial_cells@meta.data),rownames(Epithelial_cells@meta.data),paste0("Malignant_epithelial",Epithelial_cells$maligan_type))
maligant_Epithelial_cells$anno1[maligant_Epithelial_cells$anno1=="Malignant_epithelial4"]="Non_malignant_epithelial"


# 读入 空转infercnv  获取infercnv树的cluster
library(ape)

spatial_cnv<-readRDS('/01.RunCNV/run.final.infercnv_obj')
tree2<-read.tree(file = "/01.RunCNV/infercnv.observations_dendrogram.txt")
plot(tree2, edge.color="blue",show.tip.label = FALSE)
hc <- as.hclust(tree2)
k <- 6  
clusters <- cutree(hc, k)
cluster_labels <- paste0("Cluster", 1:k)
class(clusters)
print(clusters)
View(data.frame(clusters))
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = k, groupLabels = cluster_labels)
plot(dend)

df1<-data.frame(clusters)
integrated_spatial_carcinoma_regions$maligant_type<-mapvalues(rownames(integrated_spatial_carcinoma_regions@meta.data),rownames(df1),df1$clusters)
integrated_spatial_carcinoma_regions@meta.data$maligant_type[grepl("_",integrated_spatial_carcinoma_regions@meta.data$maligant_type)]="other"
integrated_spatial_carcinoma_regions@meta.data$maligant_type[integrated_spatial_carcinoma_regions@meta.data$maligant_type=='4']="Malignant_epithelial"
integrated_spatial_carcinoma_regions@meta.data$maligant_type[integrated_spatial_carcinoma_regions@meta.data$maligant_type=='Malignant_epithelial']="Non_malignant_epithelial"
integrated_spatial_carcinoma_regions@meta.data$maligant_type[integrated_spatial_carcinoma_regions@meta.data$maligant_type=='1']="Malignant_epithelial1"
integrated_spatial_carcinoma_regions@meta.data$maligant_type[integrated_spatial_carcinoma_regions@meta.data$maligant_type=='2']="Malignant_epithelial2"
integrated_spatial_carcinoma_regions@meta.data$maligant_type[integrated_spatial_carcinoma_regions@meta.data$maligant_type=='3']="Malignant_epithelial3"
integrated_spatial_carcinoma_regions@meta.data$maligant_type[integrated_spatial_carcinoma_regions@meta.data$maligant_type=='5']="Malignant_epithelial5"
integrated_spatial_carcinoma_regions@meta.data$maligant_type[integrated_spatial_carcinoma_regions@meta.data$maligant_type=='6']="Malignant_epithelial6"


integrated_spatial_carcinoma_regions<-readRDS(file = '/data/integrated_spatial_carcinoma_regions_add_maligant_type.rds')
integrated_spatial<-readRDS(file = '/data/integrated_spatial_add_area.rds')

integrated_spatial
integrated_spatial<-ScaleData(integrated_spatial,features = rownames(integrated_spatial))
sig_marker<-fread("/01.spatial_integrate/result/sig_marker.xls",data.table =F)
top10 <- sig_marker %>% dplyr::filter(!grepl("^RPL|^RPS",gene,perl = T),!grepl("^MT-|^mt-|^Mt-",gene,perl = T))%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
DefaultAssay(integrated_spatial)<-'SCT'
dim(integrated_spatial)[2]
DoHeatmap( subset(integrated_spatial, downsample = 500), features = top10$gene, size = 3)

#######################################################################
NTT.1<-readRDS('/03.spatial_deconvolution/deconvolution_NTT/NTT.1_df.rds')
NTT.2<-readRDS('/03.spatial_deconvolution/deconvolution_NTT/NTT.2_df.rds')
NTT.3<-readRDS('/03.spatial_deconvolution/deconvolution_NTT/NTT.3_df.rds')
NTT.4<-readRDS('/03.spatial_deconvolution/deconvolution_NTT/NTT.4_df.rds')


NTT<-do.call(rbind,list(NTT.1,NTT.2,NTT.3,NTT.4))

SpatialDimPlot(integrated_spatial_carcinoma_regions,ncol = 4)

integrated_spatial_carcinoma_regions$maligant_type1<-mapvalues(rownames(integrated_spatial_carcinoma_regions@meta.data),NTT$Cell,NTT$CT)




adjust_annotation_singlecell<-readRDS('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-12-4/data/adjust_annotation_singlecell.rds')
DimPlot(adjust_annotation_singlecell,label = T)
table(adjust_annotation_singlecell$orig.ident)

View(adjust_annotation_singlecell@meta.data)
unique(adjust_annotation_singlecell$orig.ident)
str(adjust_annotation_singlecell@meta.data)

adjust_annotation_singlecell$sample_type<-sapply(strsplit(as.character(adjust_annotation_singlecell$orig.ident),'-'),'[[',2)

View(adjust_annotation_singlecell@meta.data)
levels(adjust_annotation_singlecell)
all_Epithelial_cells<-WhichCells(adjust_annotation_singlecell,idents="Epithelial_cells")
table(adjust_annotation_singlecell$sample_type)
n_barcode=rownames(adjust_annotation_singlecell@meta.data)[adjust_annotation_singlecell$sample_type=='N']
n_epi_barcode=intersect(n_barcode,all_Epithelial_cells)
length(n_epi_barcode)


########################################################## deconvolution

adjust_annotation_singlecell_not_normal_epi<-subset(adjust_annotation_singlecell,cells=n_epi_barcode,invert =T)

TT.1<-readRDS('/03.spatial_deconvolution/CARD/TT.1_CARD_visium.rds')
TT.2<-readRDS('/03.spatial_deconvolution/CARD/TT.2_CARD_visium.rds')
TT.3<-readRDS('/03.spatial_deconvolution/CARD/TT.3_CARD_visium.rds')
TT.4<-readRDS('/03.spatial_deconvolution/CARD/TT.4_CARD_visium.rds')

# TT.1_proportion<-TT.1@Proportion_CARD
View(TT.1_proportion)
class(TT.1_proportion)
TT.1_proportion<-TT.1@Proportion_CARD%>%data.frame()
TT.2_proportion<-TT.2@Proportion_CARD%>%data.frame()
TT.3_proportion<-TT.3@Proportion_CARD%>%data.frame()
TT.4_proportion<-TT.4@Proportion_CARD%>%data.frame()

####################################

TT.1_proportion$dominant_celltype <- apply(TT.1_proportion, 1, function(x) names(TT.1_proportion)[which.max(x)])
TT.2_proportion$dominant_celltype <- apply(TT.2_proportion, 1, function(x) names(TT.2_proportion)[which.max(x)])
TT.3_proportion$dominant_celltype <- apply(TT.3_proportion, 1, function(x) names(TT.3_proportion)[which.max(x)])
TT.4_proportion$dominant_celltype <- apply(TT.4_proportion, 1, function(x) names(TT.4_proportion)[which.max(x)])

##################################NTT
library(plyr)
library(dplyr)

NTT.1<-readRDS('/03.spatial_deconvolution/CARD_NTT/NTT.1_CARD_visium.rds')
NTT.2<-readRDS('/03.spatial_deconvolution/CARD_NTT/NTT.2_CARD_visium.rds')
NTT.3<-readRDS('/03.spatial_deconvolution/CARD_NTT/NTT.3_CARD_visium.rds')
NTT.4<-readRDS('/03.spatial_deconvolution/CARD_NTT/NTT.4_CARD_visium.rds')

NTT.1_proportion<-NTT.1@Proportion_CARD%>%data.frame()
NTT.2_proportion<-NTT.2@Proportion_CARD%>%data.frame()
NTT.3_proportion<-NTT.3@Proportion_CARD%>%data.frame()
NTT.4_proportion<-NTT.4@Proportion_CARD%>%data.frame()



NTT.1_proportion$dominant_celltype <- apply(NTT.1_proportion, 1, function(x) names(NTT.1_proportion)[which.max(x)])
NTT.2_proportion$dominant_celltype <- apply(NTT.2_proportion, 1, function(x) names(NTT.2_proportion)[which.max(x)])
NTT.3_proportion$dominant_celltype <- apply(NTT.3_proportion, 1, function(x) names(NTT.3_proportion)[which.max(x)])
NTT.4_proportion$dominant_celltype <- apply(NTT.4_proportion, 1, function(x) names(NTT.4_proportion)[which.max(x)])


total_proportion<-do.call(rbind,list(TT.1_proportion,TT.2_proportion,TT.3_proportion,TT.4_proportion,NTT.1_proportion,NTT.2_proportion,NTT.3_proportion,NTT.4_proportion))


carcinoma_regions_integrated1$cell_type <-mapvalues(colnames(carcinoma_regions_integrated1),rownames(total_proportion),total_proportion$dominant_celltype)


Epithelial_cells1$maligant_type<- ifelse(Epithelial_cells1$Epithelial_cells>0.6,'Malignant',"Non_malignant")

carcinoma_regions_integrated1$decon_type<-carcinoma_regions_integrated1$cell_type
Epithelial_cells1$maligant_type

carcinoma_regions_integrated1$decon_type[colnames(carcinoma_regions_integrated1)%in%rownames(Epithelial_cells1)[Epithelial_cells1$maligant_type=="Malignant"]]="Malignant"
table(carcinoma_regions_integrated1$decon_type)
carcinoma_regions_integrated1$decon_type=ifelse(carcinoma_regions_integrated1$decon_type=="Malignant","Malignant","Non_malignant")
DimPlot(carcinoma_regions_integrated1,group.by = 'decon_type',cols = c('red','grey'))

#######################################################

Malignant_epithelial_barcode<-rownames(carcinoma_regions_integrated1@meta.data)[carcinoma_regions_integrated1$maligant_type=="Malignant_epithelial"]
Malignant_epithelial_barcode1<-rownames(carcinoma_regions_integrated1@meta.data)[carcinoma_regions_integrated1$decon_type=="Malignant"]
maligant_barcode<-unique(c(Malignant_epithelial_barcode1,Malignant_epithelial_barcode))
carcinoma_regions_integrated1$combine_cnv_decon<-ifelse(rownames(carcinoma_regions_integrated1@meta.data)%in%maligant_barcode,"Malignant_epithelial","Non_malignant_epithelial")



DimPlot(carcinoma_regions_integrated1,group.by ="combine_cnv_decon",cols = c('red','grey') ) +ggtitle("")
Idents(carcinoma_regions_integrated1)<-carcinoma_regions_integrated1$seurat_clusters
DimPlot(carcinoma_regions_integrated1,cols = color_panel[1:5])

DimPlot(carcinoma_regions_integrated1,cols = color_panel[1:5],split.by = 'orig.ident',ncol=4)
colnames(carcinoma_regions_integrated1@meta.data)
DimPlot(carcinoma_regions_integrated1,group.by = 'maligant_type',cols = c('red','grey'))



DimPlot(integrated_spatial_carcinoma_regions)
View(integrated_spatial_carcinoma_regions@meta.data)
dim(integrated_spatial_carcinoma_regions)
# dim(integrated_spatial)
View(integrated_spatial@meta.data)
# saveRDS(total_proportion,file = "~/personality/2023-8-8_jingdi_wuhanxiangya/2023-12-11/data/total_proportion.rds",compress = F)
# saveRDS(integrated_spatial_carcinoma_regions,file = "/data/integrated_spatial_carcinoma_regions_2023-12-12.rds",compress = F)
# saveRDS(integrated_spatial_carcinoma_regions,file = "/data/integrated_spatial_carcinoma_regions_2023-12-12.rds",compress = F)
integrated_spatial_carcinoma_regions<-readRDS("/data/integrated_spatial_carcinoma_regions_2023-12-12.rds")

# saveRDS(integrated_spatial_carcinoma_regions,file = "/data/integrated_spatial_carcinoma_regions_2023-12-12_test.rds",compress = T)

integrated_spatial_carcinoma_regions<-readRDS("/data/integrated_spatial_carcinoma_regions_2023-12-12.rds")

carcinoma_regions_integrated1<-readRDS("/data/carcinoma_regions_integrated1.rds")

SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4)

############################################################

# 读入 所有反卷积的比例的数据框图
total_proportion<-readRDS("/data/total_proportion.rds")

Epithelial_cells1<-subset(total_proportion,dominant_celltype=="Epithelial_cells")


Epithelial_cells1$maligant_type<- ifelse(Epithelial_cells1$Epithelial_cells>0.5,'Malignant',"Non_malignant")
DimPlot(carcinoma_regions_integrated1,group.by = 'maligant_type',cols = c('red','grey'))

carcinoma_regions_integrated1$decon_type[colnames(carcinoma_regions_integrated1)%in%rownames(Epithelial_cells1)[Epithelial_cells1$maligant_type=="Malignant"]]="Malignant"
carcinoma_regions_integrated1$decon_type[colnames(carcinoma_regions_integrated1)%in%rownames(Epithelial_cells1)[Epithelial_cells1$maligant_type=="Non_malignant"]]="Non_malignant"
DimPlot(carcinoma_regions_integrated1,group.by = 'decon_type',cols = c('red','grey'))


##################################################
##################################################
# total_proportion
carcinoma_regions_integrated1$decon_type[colnames(carcinoma_regions_integrated1)%in%rownames(Epithelial_cells1)]="Malignant"

DimPlot(carcinoma_regions_integrated1,group.by = 'decon_type',cols = c('red','grey'))
DimPlot(carcinoma_regions_integrated1,group.by = 'maligant_type',cols = c('red','grey'))



DimPlot(carcinoma_regions_integrated1)
SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4)

############################################################
############################################################

Malignant_epithelial_barcode_cnv<-rownames(carcinoma_regions_integrated1@meta.data)[carcinoma_regions_integrated1$maligant_type_cnv=="Malignant_epithelial"]
Malignant_epithelial_barcode_decon<-rownames(carcinoma_regions_integrated1@meta.data)[carcinoma_regions_integrated1$decon_type=="Malignant"]
combined_maligant_barcode<-unique(c(Malignant_epithelial_barcode_decon,Malignant_epithelial_barcode_cnv))
carcinoma_regions_integrated1$combine_cnv_decon<-ifelse(rownames(carcinoma_regions_integrated1@meta.data)%in%combined_maligant_barcode,"Malignant cell","Nonmalignant cell")
DimPlot(carcinoma_regions_integrated1,group.by ="combine_cnv_decon",cols = c('red','grey') ) +ggtitle("")



integrated_spatial$B_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$B_cells))
integrated_spatial$CD4._T_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$CD4._T_cells))
integrated_spatial$CD8._T_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$CD8._T_cells))
integrated_spatial$cDC<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$cDC))
integrated_spatial$Endothelial_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$Endothelial_cells))
integrated_spatial$Nonmalignant_epithelial_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$Nonmalignant_epithelial_cells))
integrated_spatial$Fibroblast_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$Fibroblast_cells))
integrated_spatial$Macrophage<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$Macrophage))
integrated_spatial$Mast_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$Mast_cells))
integrated_spatial$SMC<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$SMC))
integrated_spatial$Malignant_epithelial_cells<-mapvalues(colnames(integrated_spatial),total_proportion1_total1$barcode,as.numeric(total_proportion1_total1$Malignant_epithelial_cells))




integrated_spatial$B_cells<-as.numeric(integrated_spatial$B_cells)
integrated_spatial$CD4._T_cells<-as.numeric(integrated_spatial$CD4._T_cells)
integrated_spatial$CD8._T_cells<-as.numeric(integrated_spatial$CD8._T_cells)
integrated_spatial$cDC<-as.numeric(integrated_spatial$cDC)
integrated_spatial$Endothelial_cells<-as.numeric(integrated_spatial$Endothelial_cells)
integrated_spatial$Nonmalignant_epithelial_cells<-as.numeric(integrated_spatial$Nonmalignant_epithelial_cells)
integrated_spatial$Fibroblast_cells<-as.numeric(integrated_spatial$Fibroblast_cells)
integrated_spatial$Macrophage<-as.numeric(integrated_spatial$Macrophage)
integrated_spatial$Mast_cells<-as.numeric(integrated_spatial$Mast_cells)
integrated_spatial$SMC<-as.numeric(integrated_spatial$SMC)
integrated_spatial$Malignant_epithelial_cells<-as.numeric(integrated_spatial$Malignant_epithelial_cells)




#########################################################################  malignant color
setwd('/01.combined_decon_cnv_umap/each_sample')

count=1
class(tmp@images[[1]])
DefaultAssay(carcinoma_regions_integrated1)<-"Spatial"
cols = c("#DE2D26","grey50")
names(cols)<-levels(carcinoma_regions_integrated1)
for(sample in unique(carcinoma_regions_integrated1$orig.ident)){
  tmp<-subset(carcinoma_regions_integrated1,orig.ident==sample)
  tmp1<-subset(carcinoma_regions_integrated1,orig.ident==sample)

  Idents(tmp)<-tmp$combine_cnv_decon
  names(cols)<-levels(tmp)

  tmp@images<-tmp1@images[sprintf("slice1_%s",sample)]
  p1<-SpatialDimPlot(tmp,group.by = 'combine_cnv_decon',cols = cols)
  ggsave(p1,filename = paste0(sample,"_spatial_malignant.pdf"),width = 15,height = 8)


}
Idents(carcinoma_regions_integrated1)<-carcinoma_regions_integrated1$combine_cnv_decon
cols = c( "grey50","#DE2D26")
names(cols)<-levels(carcinoma_regions_integrated1)
p1<-SpatialDimPlot(carcinoma_regions_integrated1,group.by = 'combine_cnv_decon',cols = cols,ncol = 4)
ggsave(p1,filename = paste0("total_spatial_malignant.pdf"),width = 15,height = 8)




#############################cluster heatmap

sig_marker<-fread('/01.carcinoma_integrate/result/sig_marker.xls',data.table = F)

top10<-sig_marker%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)

Idents(carcinoma_regions_integrated1)<-carcinoma_regions_integrated1$seurat_clusters
DefaultAssay(carcinoma_regions_integrated1)<-'SCT'
DoHeatmap(subset(carcinoma_regions_integrated1,downsample=200),features = top10$gene,size = 1.5)+theme(
  axis.text.y = element_text(size = 4))



#########################################
#########################################

colnames(carcinoma_regions_integrated1@meta.data)
cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,decon_celltype)%>%group_by(decon_celltype)%>%table()%>%prop.table(margin = 1)%>%data.frame()
cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(decon_celltype,orig.ident)%>%table()%>%prop.table(margin = 1)%>%data.frame()

View(cell_type_bar)

levels(carcinoma_regions_integrated1)
cell_type_bar$orig.ident<-factor(cell_type_bar$orig.ident,levels=c(paste0("TT.",c(1:4)),paste0("NTT.",c(1:4))))

cell_type_bar$decon_celltype<-factor(cell_type_bar$decon_celltype,levels=levels(carcinoma_regions_integrated1))
colnames(cell_type_bar)[1]<-"celltype"
colnames(cell_type_bar)
color_select
unique(cell_type_bar$orig.ident)
names(color_select)<-unique(cell_type_bar$orig.ident)
ggplot(data = cell_type_bar,aes(x=celltype,y=Freq,fill=orig.ident))+
  geom_bar(position = 'stack',stat="identity")+
  scale_fill_manual(values = color_select)+
  xlab(label = "samples")+ylab(label = 'proportion(%)')+theme(axis.text.x = element_text(angle = 45))+
  theme_classic()




# p1<-SpatialFeaturePlot(integrated_spatial,features = cls,ncol = 4) & scale_fill_gradientn(colours = cols)
cols=c("#43387F","#288D88","#8AD24C","#FBE625")


cell_type <-fread("./cluster.txt",data.table = F,header = F)
SpatialFeaturePlot(carcinoma_regions_integrated1,features = "B_cells",ncol = 4,crop = F) & scale_fill_gradientn(colours = cols)

setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2023-12-22/02.spatial_celltype_distribution")
for(cls in cell_type$V1){
  p1<-SpatialFeaturePlot(carcinoma_regions_integrated1,features = cls,ncol = 4,crop = F) & scale_fill_gradientn(colours = cols)
  ggsave(p1,filename = paste0(cls,"_featureplot.pdf"),width = 15,height = 8)
}


###############################################################barplot

names(color_select)<-c("Malignant_epithelial_cells","Nonmalignant_epithelial_cells","CD4._T_cells","Fibroblast_cells","SMC",
                       "Endothelial_cells","B_cells","Mast_cells","Macrophage","CD8._T_cells")

colnames(carcinoma_regions_integrated1@meta.data)

colnames(carcinoma_regions_integrated1@meta.data)
cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,dominant_celltype_guozhiyan)%>%table()%>%prop.table(margin = 1)%>%data.frame()

View(cell_type_bar)

levels(carcinoma_regions_integrated1)
cell_type_bar$orig.ident<-factor(cell_type_bar$orig.ident,levels=c(paste0("TT.",c(1:4)),paste0("NTT.",c(1:4))))

cell_type_bar$dominant_celltype_guozhiyan<-factor(cell_type_bar$dominant_celltype_guozhiyan,levels=levels(carcinoma_regions_integrated1))
colnames(cell_type_bar)[2]<-"celltype"
colnames(cell_type_bar)
ggplot(data = cell_type_bar,aes(x=orig.ident,y=Freq,fill=celltype))+
  geom_bar(position = 'stack',stat="identity")+
  scale_fill_manual(values = color_select)+
  xlab(label = "samples")+ylab(label = 'proportion(%)')+theme(axis.text.x = element_text(angle = 45))+
  theme_classic()







cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(dominant_celltype_guozhiyan,orig.ident)%>%table()%>%prop.table(margin = 1)%>%data.frame()



levels(carcinoma_regions_integrated1)
cell_type_bar$orig.ident<-factor(cell_type_bar$orig.ident,levels=c(paste0("TT.",c(1:4)),paste0("NTT.",c(1:4))))

cell_type_bar$dominant_celltype_guozhiyan<-factor(cell_type_bar$dominant_celltype_guozhiyan,levels=levels(carcinoma_regions_integrated1))
colnames(cell_type_bar)[1]<-"celltype"
colnames(cell_type_bar)
color_select
unique(cell_type_bar$orig.ident)
names(color_select)<-unique(cell_type_bar$orig.ident)
ggplot(data = cell_type_bar,aes(x=celltype,y=Freq,fill=orig.ident))+
  geom_bar(position = 'stack',stat="identity")+
  scale_fill_manual(values = color_select)+
  xlab(label = "cluster")+ylab(label = 'proportion(%)')+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))

######################################################################
######################################################################

#################################################################2023-12-22


levels(carcinoma_regions_integrated1)
carcinoma_regions_integrated1$group<-gsub("\\..*","",carcinoma_regions_integrated1$orig.ident)

celltype=c("CD4._T_cells","CD8._T_cells","Endothelial_cells","Macrophage")
for(cls in celltype){
  diff<-FindMarkers(carcinoma_regions_integrated1,group.by = "group",ident.1 = "TT",ident.2 = "NTT",subset.ident = cls)
  fwrite(diff,file = paste0(cls,"-TT_vs_NTT.xls"),sep = '\t',col.names = T,row.names = T)
}



##################################################################2023-12-25



carcinoma_regions_integrated1<-readRDS('~/personality/2023-8-8_jingdi_wuhanxiangya/2023-12-25/data/carcinoma_regions_integrated1.rds')

View(carcinoma_regions_integrated1@meta.data)
colnames(carcinoma_regions_integrated1@meta.data)
DimPlot(carcinoma_regions_integrated1,group.by = 'combine_cnv_decon')

malignant_infercnv_barcode<-rownames(subset(carcinoma_regions_integrated1@meta.data,combine_cnv_decon=="Malignant cell"))
colnames(cna_predict)

cna_predict_barcode<-subset(cna_predict,copykat.pred=="aneuploid")$cell.names

table(malignant_infercnv_barcode%in%cna_predict_barcode)
DimPlot(carcinoma_regions_integrated1,cells.highlight = cna_predict_barcode,sizes.highlight = 0.1)

carcinoma_regions_integrated1$cna<-mapvalues(colnames(carcinoma_regions_integrated1),cna_predict$cell.names,cna_predict$copykat.pred)


######################################### cna process


carcinoma_regions_integrated1$cna[grepl("-1",carcinoma_regions_integrated1$cna)]="other"

carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,cna)%>%group_by(orig.ident)%>%table()

carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,combine_cnv_decon)%>%group_by(orig.ident)%>%table()


cnv_rds<-readRDS('01.RunCNV/run.final.infercnv_obj')


colmean_normal_cnv<-colMeans(cnas_normal)

################################################################ 处理tumor cna
cnas_tumors <- cnv_rds@expr.data[,cnv_rds@observation_grouped_cell_indices$`carcinoma regions`] %>%data.frame(check.rows = F,check.names = F)


cnas_tumors$chr<-cnv_rds@gene_order$chr
cnas_tumors$chr<-NULL
# 加到第一列的方法
cnas_tumors<-cbind(chr=cnv_rds@gene_order$chr,cnas_tumors)



chrom_average<-cnas_tumors%>%dplyr::group_by(chr)%>%summarise(across(everything(),mean,na.rm=TRUE))


chrom_average_df<-t(chrom_average)%>%data.frame(check.names = F,check.rows = F)


chrom_average_df<-setNames(chrom_average_df[-1,],chrom_average_df[1,])

chrom_average_df1<-chrom_average_df%>%mutate_all(~as.numeric(as.character(.)))

calc_frac<-function(x){
  print(x)
  if(mean(x)>1){
    frac = (x-1)/(max(x)-1)
  }else{
    frac = (1-x)/(1-min(x))
  }
  frac
}


chrom_average_df1_stat<-apply(chrom_average_df1,1,calc_frac)

cell_frac = apply(chrom_average_df1_stat,2,max)


######################################## normal cna
cnas_normal <- cnv_rds@expr.data[,cnv_rds@reference_grouped_cell_indices$`noncarcinoma regions`]  %>%data.frame(check.rows = F,check.names = F)
identical(rownames(cnv_rds@gene_order),rownames(cnas_normal)) # TRUE

# 加到第一列的方法
cnas_normal<-cbind(chr=cnv_rds@gene_order$chr,cnas_normal)

chrom_average_normal<- cnas_normal%>%dplyr::group_by(chr)%>%summarise(across(everything(),mean,na.rm=T))


chrom_average_normal<-t(chrom_average_normal)%>%data.frame(check.rows = F,check.names = F)

chrom_average_normal<-setNames(chrom_average_normal[-1,],chrom_average_normal[1,])


# chrom_average_normal<-data.frame(lapply(chrom_average_normal, function(x){as.numeric(x)})) # 这种方法，行名不见了
chrom_average_normal1 <- chrom_average_normal%>%mutate_all(~as.numeric(as.character(.)))

chrom_average_normal_spotstat<-apply(chrom_average_normal1,2,median)

######################################################judgment normal or non-normal


# chrom_average_tumor_spotstat<-chrom_average_df1 - chrom_average_normal_spotstat
result <- sweep(chrom_average_df1, 2, chrom_average_normal_spotstat, FUN = "-")

calc_frac<-function(x){
  max(abs(x))
}


table(apply(result,1,calc_frac) >0.12)
total_barcode<-rownames(result)

sum(maligant_score$chr1)/dim(maligant_score)[1]
sum(result$chr1)/dim(result)[1]
sum(result$chr2)/dim(result)[1]
sum(result$chr3)/dim(result)[1]
sum(result$chr4)/dim(result)[1]
sum(result$chr5)/dim(result)[1]
sum(result$chr6)/dim(result)[1]
sum(result$chr7)/dim(result)[1]
sum(result$chr8)/dim(result)[1]
sum(result$chr9)/dim(result)[1]
sum(result$chr10)/dim(result)[1]
sum(result$chr11)/dim(result)[1]
sum(result$chr12)/dim(result)[1]
sum(result$chr13)/dim(result)[1]
sum(result$chr14)/dim(result)[1]
sum(result$chr15)/dim(result)[1]
sum(result$chr16)/dim(result)[1]
sum(result$chr17)/dim(result)[1]
sum(result$chr18)/dim(result)[1]
sum(result$chr19)/dim(result)[1]
sum(result$chr20)/dim(result)[1]
sum(result$chr21)/dim(result)[1]
sum(result$chr22)/dim(result)[1]



# 0.1
result <- sweep(chrom_average_df1, 2, chrom_average_normal_spotstat, FUN = "-")
View(result)
calc_frac<-function(x){
  max(abs(x))
}
total_barcode<-rownames(result)
malignant_barcode<-total_barcode[apply(result,1,calc_frac) >=0.1]
nonmalignant_barcode<-total_barcode[apply(result,1,calc_frac) <0.1]
# malignant_barcode
sample_stat<-sapply(strsplit(malignant_barcode,"_"),'[[',1)
sample_normal_stat<-sapply(strsplit(nonmalignant_barcode,"_"),'[[',1)



#################################################################cnas union with deconvolution

carcinoma_regions_integrated1$dominant_celltype_combined<-carcinoma_regions_integrated1$dominant_celltype_guozhiyan
carcinoma_regions_integrated1$dominant_celltype_combined <- ifelse(rownames(carcinoma_regions_integrated1@meta.data) %in% names(malignant_barcode_mapping),
                                                                   malignant_barcode_mapping[rownames(carcinoma_regions_integrated1@meta.data)],
                                                                   carcinoma_regions_integrated1$dominant_celltype_guozhiyan)



Idents(carcinoma_regions_integrated1)<-carcinoma_regions_integrated1$dominant_celltype_combined

color_select=c("red",color_panel[1:9])
names(color_select)<-levels(carcinoma_regions_integrated1)
names(color_select)<-c("Malignant_epithelial_cells","Nonmalignant_epithelial_cells","CD4._T_cells","Fibroblast_cells","SMC",
                       "Endothelial_cells","B_cells","Mast_cells","Macrophage","CD8._T_cells")
SpatialDimPlot(carcinoma_regions_integrated1,ncol = 2,crop = F,cols = color_select)


####################################################################
carcinoma_regions_integrated1$dominant_celltype_combined_0.1<-carcinoma_regions_integrated1$dominant_celltype_combined



df<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,dominant_celltype_combined_0.1)%>%group_by(orig.ident)%>%table()%>%data.frame()%>%tidyr::spread(key=dominant_celltype_combined_0.1,value=Freq)



############################################ 0.1

colnames(carcinoma_regions_integrated1@meta.data)
carcinoma_regions_integrated1$malignant_type<-carcinoma_regions_integrated1$dominant_celltype_combined_0.1
celltype<-unique(carcinoma_regions_integrated1$malignant_type)
celltype1<-setdiff(celltype,c("Malignant_epithelial_cells"))

carcinoma_regions_integrated1$malignant_type[carcinoma_regions_integrated1$malignant_type%in%celltype1]="Nonmalignant cells"
carcinoma_regions_integrated1$malignant_type[carcinoma_regions_integrated1$malignant_type=="Malignant_epithelial_cells"]="Malignant cells"

DimPlot(carcinoma_regions_integrated1,group.by = 'malignant_type',cols = c("red","grey"))+ggtitle("")

colnames(carcinoma_regions_integrated1@meta.data)
colnames(carcinoma_regions_integrated1@meta.data)[30]<-"decon_cnas_combined_malignant_type"
DimPlot(carcinoma_regions_integrated1,group.by = 'decon_cnas_combined_malignant_type',cols = c("red","grey"))+ggtitle("")
#######################################
#######################################

colnames(carcinoma_regions_integrated1@meta.data)
DimPlot(carcinoma_regions_integrated1,group.by = 'decon_celltype')+ggtitle("")
#######################################
#######################################
# cnas only
malignant_barcode<-total_barcode[apply(result,1,calc_frac) >=0.1]
nonmalignant_barcode<-total_barcode[apply(result,1,calc_frac) <0.1]

length(malignant_barcode)
length(nonmalignant_barcode)
cell_type_mapping<-data.frame(barcode=c(malignant_barcode,nonmalignant_barcode),type=c(rep("Malignant_epithelial_cells",length(malignant_barcode)),rep("Nonmalignant_epithelial_cells",length(nonmalignant_barcode))))


carcinoma_regions_integrated1$only_cnas<-carcinoma_regions_integrated1$malignant_type

carcinoma_regions_integrated1$only_cnas<-mapvalues(colnames(carcinoma_regions_integrated1),cell_type_mapping$barcode,cell_type_mapping$type)


DimPlot(carcinoma_regions_integrated1,group.by = 'only_cnas',cols = c("red","grey"))+ggtitle("")


######################################################## 2023-12-28
colnames(carcinoma_regions_integrated1@meta.data)
table(carcinoma_regions_integrated1$dominant_celltype_guozhiyan)
carcinoma_regions_integrated1$only_decon<-carcinoma_regions_integrated1$dominant_celltype_guozhiyan

yy1<-unique(carcinoma_regions_integrated1$only_decon)
yy2<-setdiff(yy1,c("Malignant_epithelial_cells"))

carcinoma_regions_integrated1$only_decon

carcinoma_regions_integrated1$only_decon[carcinoma_regions_integrated1$only_decon%in%yy2]="Nonmalignant_epithelial_cells"
carcinoma_regions_integrated1$only_decon[carcinoma_regions_integrated1$only_decon=="Malignant_epithelial_cells"]="Malignant_epithelial_cells"
table(carcinoma_regions_integrated1$only_decon)
DimPlot(carcinoma_regions_integrated1,group.by = 'only_decon',cols = c("red","grey"))+ggtitle("")



##############################################################
# 修改图例

table(carcinoma_regions_integrated1$decon_cnas_combined_malignant_type)
carcinoma_regions_integrated1$decon_cnas_combined_malignant_type[carcinoma_regions_integrated1$decon_cnas_combined_malignant_type=="Malignant cells"]="Malignant_epithelial_cells"
carcinoma_regions_integrated1$decon_cnas_combined_malignant_type[carcinoma_regions_integrated1$decon_cnas_combined_malignant_type=="Nonmalignant cells"]="Nonmalignant_epithelial_cells"
table(carcinoma_regions_integrated1$only_cnas)
DimPlot(carcinoma_regions_integrated1,group.by = 'decon_cnas_combined_malignant_type',cols = c("red","grey"))+ggtitle("")


###################celltype color selected
color_panel <-c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
                '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
                '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
                '#CCC9E6','#3A6963','#68A180','#476D87','#9FA3A8','#968175')
# 用这种颜色
color_select=c("red",color_panel[1:3],color_panel[25],color_panel[22],color_panel[7],color_panel[17],color_panel[19:20])

# color_select=c("red",color_panel[1:3],color_panel[28],color_panel[22],color_panel[16:17],color_panel[19:20])


names(color_select)<-levels(carcinoma_regions_integrated1)
names(color_select)<-c("Malignant_epithelial_cells","Nonmalignant_epithelial_cells","CD4._T_cells","Fibroblast_cells","SMC",
                       "Endothelial_cells","B_cells","Mast_cells","Macrophage","CD8._T_cells")

DimPlot(carcinoma_regions_integrated1,cols = color_select,pt.size = 1.5)

DimPlot(carcinoma_regions_integrated1,cols = color_select)

SpatialDimPlot(carcinoma_regions_integrated1,ncol = 2,crop = F,cols = color_select)


########################################
######################################## 做组装图统计


# barplot
names(color_select)<-c("Malignant_epithelial_cells","Nonmalignant_epithelial_cells","CD4._T_cells","Fibroblast_cells","SMC",
                       "Endothelial_cells","B_cells","Mast_cells","Macrophage","CD8._T_cells")

colnames(carcinoma_regions_integrated1@meta.data)
# 柱状图
colnames(carcinoma_regions_integrated1@meta.data)
table(carcinoma_regions_integrated1$decon_cnas_combined_malignant_type)
table(carcinoma_regions_integrated1$dominant_celltype_combined_0.1)
table(carcinoma_regions_integrated1$group)
cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(group,dominant_celltype_combined_0.1)%>%table()%>%prop.table(margin = 1)%>%data.frame()



levels(carcinoma_regions_integrated1)
# cell_type_bar$orig.ident<-factor(cell_type_bar$orig.ident,levels=c(paste0("TT.",c(1:4)),paste0("NTT.",c(1:4))))
cell_type_bar$group<-factor(cell_type_bar$group,levels=c("TT","NTT"))

cell_type_bar$dominant_celltype_combined_0.1<-factor(cell_type_bar$dominant_celltype_combined_0.1,levels=levels(carcinoma_regions_integrated1))
colnames(cell_type_bar)[2]<-"celltype"
colnames(cell_type_bar)

ggplot(data = cell_type_bar,aes(x=group,y=Freq,fill=celltype))+
  geom_bar(position = 'stack',stat="identity")+
  scale_fill_manual(values = color_select)+
  xlab(label = "samples")+ylab(label = 'proportion(%)')+theme(axis.text.x = element_text(angle = 45))+
  theme_classic()






colnames(carcinoma_regions_integrated1@meta.data)
cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(dominant_celltype_combined_0.1,group)%>%table()%>%prop.table(margin = 1)%>%data.frame()




cell_type_bar$group<-factor(cell_type_bar$group,levels=c("TT","NTT"))

cell_type_bar$dominant_celltype_combined_0.1<-factor(cell_type_bar$dominant_celltype_combined_0.1,levels=levels(carcinoma_regions_integrated1))
# cell_type_bar$celltype<-factor(cell_type_bar$celltype,levels=levels(carcinoma_regions_integrated1))
# cell_type_bar$dominant_celltype_combined_0.1<-factor(cell_type_bar$dominant_celltype_combined_0.1,levels=levels(carcinoma_regions_integrated1))
colnames(cell_type_bar)[1]<-"celltype"
# colnames(cell_type_bar)[2]<-"Samples"
colnames(cell_type_bar)[2]<-"group"
colnames(cell_type_bar)
View(cell_type_bar)
color_select
unique(cell_type_bar$orig.ident)
names(color_select)<-c("#C1E6F3","#23452F")
color_select<-c("#58A4C3","#AB3282")
color_select<-c("#AB3282","#58A4C3")
ggplot(data = cell_type_bar,aes(x=celltype,y=Freq,fill=group))+
  geom_bar(position = 'stack',stat="identity")+
  scale_fill_manual(values = color_select)+
  xlab(label = "cluster")+ylab(label = 'proportion(%)')+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))



###########################################################
###########################################################

DimPlot(carcinoma_regions_integrated1)
unique(Idents(carcinoma_regions_integrated1))
Nonmalignant_epithelial_cells<-subset(carcinoma_regions_integrated1,idents = "Malignant_epithelial_cells",invert =T)

SpatialDimPlot(Nonmalignant_epithelial_cells,crop = F,ncol = 4)

DimPlot(carcinoma_regions_integrated1)
table(Idents(carcinoma_regions_integrated1))


total_proportion<-total_proportion%>%mutate(orig.ident=sapply(strsplit(rownames(total_proportion),"_"),'[[',1))

df1<-total_proportion%>%dplyr::mutate(barcode=rownames(total_proportion))%>%dplyr::select(dominant_celltype,orig.ident)

colnames(carcinoma_regions_integrated1@meta.data)
table(carcinoma_regions_integrated1$dominant_celltype_combined_0.1)
df2<-carcinoma_regions_integrated1@meta.data%>%
  dplyr::filter(dominant_celltype_combined_0.1=="Malignant_epithelial_cells")%>%dplyr::rename(dominant_celltype=dominant_celltype_combined_0.1)%>%
  dplyr::select(dominant_celltype,orig.ident)


df<-rbind(df1,df2)


df_res<-df%>%group_by(orig.ident)%>%table()%>%data.frame()%>%tidyr::spread(key=dominant_celltype,value=Freq)



########################################################2024-1-3

setwd("2024-1-2/03.celltype_feature")


cell_type<-colnames(carcinoma_regions_integrated1@meta.data)[16:25]


cols=c("#43387F","#288D88","#8AD24C","#FBE625")
View(carcinoma_regions_integrated1@meta.data)

setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-2/03.celltype_feature")
unique(Idents(carcinoma_regions_integrated1))
for(cls in cell_type){
  p1<-SpatialFeaturePlot(carcinoma_regions_integrated1,features = cls,ncol = 4,crop = F) & scale_fill_gradientn(colours = cols)
  ggsave(p1,filename = paste0(cls,"_featureplot.pdf"),width = 15,height = 8)
}

cls="B_cells"

SpatialFeaturePlot(carcinoma_regions_integrated1,features = cls,ncol = 4,crop = F) & scale_fill_gradientn(colours = cols)

SpatialFeaturePlot(carcinoma_regions_integrated1,features = cls,ncol = 4,crop = F) & scale_fill_gradientn(colours = cols)
SpatialFeaturePlot(carcinoma_regions_integrated1,features = cls,ncol = 4,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols)

tt.3<-subset(carcinoma_regions_integrated1,orig.ident=='TT.3')
sample="TT.3"
tt.3@images<-tt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F,max.cutoff = 0.3) & scale_fill_gradientn(colours = cols)

SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.3))

sample="NTT.3"
ntt.3<-subset(carcinoma_regions_integrated1,orig.ident=='NTT.3')
ntt.3@images<-ntt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols)
max(ntt.3$B_cells)
# na.value = "grey50"
SpatialFeaturePlot(ntt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4))
SpatialFeaturePlot(ntt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "grey50")


#######################################################


str(integrated_spatial@meta.data)
cols=c("#43387F","#288D88","#8AD24C","#FBE625")
SpatialFeaturePlot(integrated_spatial,features = "B_cells",crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,limits=c(0.0,0.4),na.value = "lightgrey")

setwd("/2024-1-2/04.celltype_feature_adjust")
View(integrated_spatial@meta.data)

SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,na.value = "lightgrey",oob=scales::squish)

for(cls in cell_type){
  # p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,limits=c(0.0,1.0),na.value = "lightgrey",oob=scales::squish)
  p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,na.value = "lightgrey",oob=scales::squish)

  ggsave(p1,filename = paste0(cls,"_featureplot.pdf"),width = 18,height = 8)
}

############################ plot
setwd("/2024-1-2/05.celltype_feature_adjust_nogrey")
cls="B_cells"

tt.3<-subset(integrated_spatial,orig.ident=='TT.3')
sample="TT.3"
tt.3@images<-tt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.35),na.value = "lightgrey",,oob=scales::squish)


ntt.3<-subset(integrated_spatial,orig.ident=='NTT.3')
sample="NTT.3"
ntt.3@images<-ntt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.5),na.value = "lightgrey")


cls="CD4._T_cells"
tt.2<-subset(integrated_spatial,orig.ident=='TT.2')
sample="TT.2"
tt.2@images<-tt.2@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.3),na.value = "lightgrey")


cls="CD8._T_cells"
tt.3<-subset(integrated_spatial,orig.ident=='TT.3')
sample="TT.3"
tt.3@images<-tt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.25),na.value = "lightgrey")



cls="CD8._T_cells"
ntt.2<-subset(integrated_spatial,orig.ident=='NTT.2')
sample="NTT.2"
ntt.2@images<-ntt.2@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey")



cls="CD8._T_cells"
ntt.3<-subset(integrated_spatial,orig.ident=='NTT.3')
sample="NTT.3"
ntt.3@images<-ntt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey")



cls="CD8._T_cells"
ntt.4<-subset(integrated_spatial,orig.ident=='NTT.4')
sample="NTT.4"
ntt.4@images<-ntt.4@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.4,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey")


cls="Macrophage"
tt.3<-subset(integrated_spatial,orig.ident=='TT.3')
sample="TT.3"
tt.3@images<-tt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.25),na.value = "lightgrey")




cls="Malignant_epithelial_cells"
tt.2<-subset(integrated_spatial,orig.ident=='TT.2')
sample="TT.2"
tt.2@images<-tt.2@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey",oob=scales::squish)
# SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.35),na.value = "lightgrey")
####################


cls="Malignant_epithelial_cells"
ntt.1<-subset(integrated_spatial,orig.ident=='NTT.1')
sample="NTT.1"
ntt.1@images<-ntt.1@images[sprintf("slice1_%s",sample)]
# SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey",oob=scales::squish)
SpatialFeaturePlot(ntt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey",oob=scales::squish)
####################



table(integrated_spatial$CD4._T_cells)

test11<-subset(integrated_spatial,orig.ident=='NTT.1')
View(test11@meta.data)
sample="NTT.1"
test11@images<-test11@images[sprintf("slice1_%s",sample)]
table(test11$CD4._T_cells=="NA")
table(is.na(test11$CD4._T_cells))

max(test11$CD4._T_cells[!is.na(test11$CD4._T_cells)])
min(test11$CD4._T_cells[!is.na(test11$CD4._T_cells)])


#

SpatialFeaturePlot(test11,features = "CD4._T_cells",crop = F,alpha = 1) & scale_fill_gradientn(colours = cols,limits=c(0,1),na.value = "lightgrey")

SpatialFeaturePlot(test11,features = "CD4._T_cells",crop = F,alpha = 1) & scale_fill_gradientn(colours = cols,na.value = "lightgrey")




setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-2/05.celltype_feature_adjust_nogrey")

for(cls in cell_type){
  p1<-SpatialFeaturePlot(carcinoma_regions_integrated1,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols)
  ggsave(p1,filename = paste0(cls,"_featureplot.pdf"),width = 15,height = 8)
}




#############################################################################
#############################################################################
#############################################################################nogrey

setwd("/05.celltype_feature_adjust_nogrey")
cls="B_cells"

tt.3<-subset(carcinoma_regions_integrated1,orig.ident=='TT.3')
sample="TT.3"
tt.3@images<-tt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.35),oob=scales::squish)


ntt.3<-subset(carcinoma_regions_integrated1,orig.ident=='NTT.3')
sample="NTT.3"
ntt.3@images<-ntt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.5),na.value = "lightgrey")


cls="CD4._T_cells"
tt.2<-subset(carcinoma_regions_integrated1,orig.ident=='TT.2')
sample="TT.2"
tt.2@images<-tt.2@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.3),na.value = "lightgrey")


cls="CD8._T_cells"
tt.3<-subset(carcinoma_regions_integrated1,orig.ident=='TT.3')
sample="TT.3"
tt.3@images<-tt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.25),na.value = "lightgrey")



cls="CD8._T_cells"
ntt.2<-subset(carcinoma_regions_integrated1,orig.ident=='NTT.2')
sample="NTT.2"
ntt.2@images<-ntt.2@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey")



cls="CD8._T_cells"
ntt.3<-subset(carcinoma_regions_integrated1,orig.ident=='NTT.3')
sample="NTT.3"
ntt.3@images<-ntt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey")



cls="CD8._T_cells"
ntt.4<-subset(carcinoma_regions_integrated1,orig.ident=='NTT.4')
sample="NTT.4"
ntt.4@images<-ntt.4@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.4,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey")


# Macrophage


cls="Macrophage"
tt.3<-subset(carcinoma_regions_integrated1,orig.ident=='TT.3')
sample="TT.3"
tt.3@images<-tt.3@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.3,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.25),na.value = "lightgrey")


# Malignant_epithelial_cells


cls="Malignant_epithelial_cells"
tt.2<-subset(carcinoma_regions_integrated1,orig.ident=='TT.2')
sample="TT.2"
tt.2@images<-tt.2@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),oob=scales::squish)
# SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.35),na.value = "lightgrey")
####################


cls="Malignant_epithelial_cells"
ntt.1<-subset(carcinoma_regions_integrated1,orig.ident=='NTT.1')
sample="NTT.1"
ntt.1@images<-ntt.1@images[sprintf("slice1_%s",sample)]
# SpatialFeaturePlot(tt.2,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),na.value = "lightgrey",oob=scales::squish)
SpatialFeaturePlot(ntt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.4),oob=scales::squish)







SpatialFeaturePlot(test11,features = "CD4._T_cells",crop = F,alpha = 1) & scale_fill_gradientn(colours = cols,limits=c(0,1),na.value = "lightgrey")

SpatialFeaturePlot(test11,features = "CD4._T_cells",crop = F,alpha = 1) & scale_fill_gradientn(colours = cols,na.value = "lightgrey")



##########################################################2024-1-5
# 空间分布图 明天还是需要帮我修改2张图。CD4._T_cells：NTT1 梯度到1.0的图
# 。CD8._T_cells：TT1 试一下梯度到0.3和0.2的图。
setwd('/2024-1-5/01.adjust_celltype_feature')

# no_grey
cls="CD4._T_cells"
ntt.1<-subset(carcinoma_regions_integrated1,orig.ident=='NTT.1')
sample="NTT.1"
ntt.1@images<-ntt.1@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,1.0))

# grey
cls="CD4._T_cells"
ntt.1<-subset(integrated_spatial,orig.ident=='NTT.1')
sample="NTT.1"
ntt.1@images<-ntt.1@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(ntt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,1.0),na.value = "lightgrey")




####################################
####################################

cls="CD8._T_cells"
tt.1<-subset(carcinoma_regions_integrated1,orig.ident=='TT.1')
sample="TT.1"
tt.1@images<-tt.1@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.3),oob=scales::squish)

cls="CD8._T_cells"
tt.1<-subset(carcinoma_regions_integrated1,orig.ident=='TT.1')
sample="TT.1"
tt.1@images<-tt.1@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.2),oob=scales::squish)


##################################
# grey
cls="CD8._T_cells"
tt.1<-subset(integrated_spatial,orig.ident=='TT.1')
sample="TT.1"
tt.1@images<-tt.1@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.3),oob=scales::squish,na.value = "lightgrey")

cls="CD8._T_cells"
tt.1<-subset(integrated_spatial,orig.ident=='TT.1')
sample="TT.1"
tt.1@images<-tt.1@images[sprintf("slice1_%s",sample)]
SpatialFeaturePlot(tt.1,features = cls,crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colours = cols,limits=c(0,0.2),oob=scales::squish,na.value = "lightgrey")


##################################################
##################################################
# 1.7
# 1.8
SpatialDimPlot(integrated_spatial,ncol = 4)


View(integrated_spatial@meta.data)
View(carcinoma_regions_integrated1@meta.data)
carcinoma_regions_integrated1$decon_cnas_combined_malignant_type
table(carcinoma_regions_integrated1$decon_cnas_combined_malignant_type)


integrated_spatial$decon_cnas_combined_malignant_type<-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$decon_cnas_combined_malignant_type)

# celltype
integrated_spatial$decon_cnas_combined_malignant_type[grepl("-1",integrated_spatial$decon_cnas_combined_malignant_type)]="Nonmalignant cells"
integrated_spatial$decon_cnas_combined_malignant_type[grepl("Nonmalignant_epithelial_cells",integrated_spatial$decon_cnas_combined_malignant_type)]="Nonmalignant cells"
integrated_spatial$decon_cnas_combined_malignant_type[grepl("Malignant_epithelial_cells",integrated_spatial$decon_cnas_combined_malignant_type)]="Malignant cells"

table(integrated_spatial$decon_cnas_combined_malignant_type)

###################################


setwd("/1.7_spatial_spot_malignant_nonmalignant")
color_select=c("red","lightgrey")

names(color_select)<-c("Malignant cells","Nonmalignant cells")
Idents(integrated_spatial)<-integrated_spatial$decon_cnas_combined_malignant_type
levels(integrated_spatial)<-c("Malignant cells","Nonmalignant cells")
SpatialDimPlot(integrated_spatial,cols = color_select,ncol = 4)

####################################

SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4)

####################################百分比柱状图
# carcinoma_regions_integrated1$decon_cnas_combined_malignant_type
cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,decon_cnas_combined_malignant_type)%>%table()%>%prop.table(margin = 1)%>%data.frame()


cell_type_bar$decon_cnas_combined_malignant_type <- as.character(cell_type_bar$decon_cnas_combined_malignant_type)

cell_type_bar$orig.ident<-factor(cell_type_bar$orig.ident,levels=c(paste0("TT.",c(1:4)),paste0("NTT.",c(1:4))))

cell_type_bar$decon_cnas_combined_malignant_type[cell_type_bar$decon_cnas_combined_malignant_type=="Malignant_epithelial_cells"]="Malignant cells"
cell_type_bar$decon_cnas_combined_malignant_type[cell_type_bar$decon_cnas_combined_malignant_type=="Nonmalignant_epithelial_cells"]="Nonmalignant cells"
# cell_type_bar$group<-factor(cell_type_bar$group,levels=c("TT","NTT"))

cell_type_bar$decon_cnas_combined_malignant_type<-factor(cell_type_bar$decon_cnas_combined_malignant_type,levels=c("Nonmalignant cells","Malignant cells"))


colnames(cell_type_bar)[1]<-"samples"
# colnames(cell_type_bar)[2]<-"Samples"
colnames(cell_type_bar)[2]<-"Malignant status"
cell_type_bar$`Malignant status`<-factor(cell_type_bar$`Malignant status`,levels=c("Nonmalignant cells","Malignant cells"))


unique(cell_type_bar$orig.ident)

color_select<-c("red","lightgrey")
ggplot(data = cell_type_bar,aes(x=samples,y=Freq,fill=`Malignant status`))+
  geom_bar(position = 'stack',stat="identity")+


  xlab(label = "cluster")+ylab(label = 'proportion(%)')+

  theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_discrete(type = color_select,limits = c("Malignant cells","Nonmalignant cells"))




Idents(carcinoma_regions_integrated1) <-carcinoma_regions_integrated1$dominant_celltype_combined_0.1
levels(carcinoma_regions_integrated1)<-names(sort(table(Idents(carcinoma_regions_integrated1)),decreasing = T))

color_panel <-c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
                '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
                '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
                '#CCC9E6','#3A6963','#68A180','#476D87','#9FA3A8','#968175')


# 用这种颜色
color_select=c("red",color_panel[1:3],color_panel[25],color_panel[22],color_panel[7],color_panel[17],color_panel[19:20])



names(color_select)<-levels(carcinoma_regions_integrated1)
names(color_select)<-c("Malignant_epithelial_cells","Nonmalignant_epithelial_cells","CD4._T_cells","Fibroblast_cells","SMC",
                       "Endothelial_cells","B_cells","Mast_cells","Macrophage","CD8._T_cells")
DimPlot(carcinoma_regions_integrated1,cols = color_select)


sig_marker <- subset(FindAllMarkers(carcinoma_regions_integrated1,only.pos = T),p_val_adj <0.05)

top10<-sig_marker%>%group_by(cluster)%>%dplyr::filter(avg_log2FC>0.8)%>%slice_head(n=10)

dim(carcinoma_regions_integrated1@assays$SCT@scale.data)

DoHeatmap(subset(carcinoma_regions_integrated1,downsample=1000),features = top10$gene)

carcinoma_regions_integrated1<-SCTransform(carcinoma_regions_integrated1,assay = "Spatial")
# PTGS2
FeaturePlot(carcinoma_regions_integrated1,features = 'PTGS2')
VlnPlot(carcinoma_regions_integrated1,split.by = 'orig.ident',features = 'PTGS2',pt.size = 0)



#################################################2024-1-8
#################################################2024-1-8
################################################# marker基因分析
View(sig_marker)
color_heatmap=c("#2E7DCC","white","#CD3435")


DoHeatmap(subset(carcinoma_regions_integrated1,downsample=500),features = top10$gene,group.colors=color_select,raster=F) +
  scale_fill_gradientn(colors=color_heatmap)+theme(axis.text.y = element_text(size = 3))




fwrite(sig_marker,file = "01.heatmap_plot/sig_marker.xls",sep="\t",col.names = T,row.names = F)


############################################################ enrichment
test1<-read.table('/02.diff/B_cells-TT_vs_NTT_p0.05.xls',sep = '\t',header = T)



hs_msigdbr <- msigdbr(species="Homo sapiens")


GO_db <- msigdbr(species="Homo sapiens",category="H") %>%
  dplyr::select(gs_name, entrez_gene, gene_symbol)

DEgenelist<-read.table('02.diff/Endothelial_cells-TT_vs_NTT_p0.05.xls',sep = '\t',header = T)

res <- enricher(DEgenelist[,1],TERM2GENE=GO_db[,c(1,3)],pvalueCutoff =0.05)
View(res)
res@result
barplot(res ,color = "p.adjust",showCategory = 10,size = NULL,split = NULL,font.size = 12,title = "")

############################################################
############################################################
#
mycol8 <-c('#FF0000','#FF8C00','#FFD700','#9ACD32','#008000', '#20B2AA','#6674A5','#0000FF','#8000FF','#800080')

# 麻烦帮我先把senic的代码跑上 然后帮我把这个堆积柱状图的颜色改一下。谢谢


names(color_select)<-c("Malignant_epithelial_cells","Nonmalignant_epithelial_cells","CD4._T_cells","Fibroblast_cells","SMC",
                       "Endothelial_cells","B_cells","Mast_cells","Macrophage","CD8._T_cells")



cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,dominant_celltype_combined_0.1)%>%table()%>%prop.table(margin = 1)%>%data.frame()




cell_type_bar$orig.ident<-factor(cell_type_bar$orig.ident,levels=c(paste0("TT.",c(1:4)),paste0("NTT.",c(1:4))))

cell_type_bar$dominant_celltype_combined_0.1<-factor(cell_type_bar$dominant_celltype_combined_0.1,levels=levels(carcinoma_regions_integrated1))

colnames(cell_type_bar)[2]<-"celltype"

color_select<-mycol8
names(color_select)<-levels(carcinoma_regions_integrated1)
ggplot(data = cell_type_bar,aes(x=orig.ident,y=Freq,fill=celltype))+
  geom_bar(position = 'stack',stat="identity")+
  scale_fill_manual(values = color_select)+
  xlab(label = "samples")+ylab(label = 'proportion(%)')+theme(axis.text.x = element_text(angle = 45))+
  theme_classic()

###################################
names(color_select)<-levels(carcinoma_regions_integrated1)

DimPlot(carcinoma_regions_integrated1,cols = color_select)
################################### 2024-1-11


SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4,crop = F,pt.size.factor = 1)
SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4,crop = F,pt.size.factor = 1.6)

setwd("/01.spatial_point_size_adjust")

for(i in c(1.0,1.1,1.2,1.3,1.4,1.5,1.6)){
  p1<-SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4,crop = F,pt.size.factor = i,cols = color_select)
  ggsave(p1,filename = paste0("point_size_",i,"_spatial.pdf"),width = 20,height = 15)
}

for(i in c(1.0,1.1,1.2,1.3,1.4,1.5)){
  p1<-DimPlot(carcinoma_regions_integrated1,pt.size = i, cols = color_select)
  ggsave(p1,filename = paste0("point_size_",i,"_umap.pdf"),width = 20,height = 15)
}



Malignant_epithelial_cells<-subset(carcinoma_regions_integrated1,idents=c("Malignant_epithelial_cells"))





average_data<-AverageExpression(Malignant_epithelial_cells_cms,assays = 'SCT',slot = 'scale.data',group.by="orig.ident")


df<-average_data$SCT%>%data.frame()
View(df)
colnames(df)
df1<-df[,c("TT.1","TT.2","TT.3","TT.4","NTT.1","NTT.2","NTT.3","NTT.4")]
cor_res<-cor(df1,method = c("pearson"))


annotation_col<-data.frame(row.names = c("TT.1","TT.2","TT.3","TT.4","NTT.1","NTT.2","NTT.3","NTT.4"),type=c(rep("TT",4),rep("NTT",4)))
annotation_row<-data.frame(row.names = c("TT.1","TT.2","TT.3","TT.4","NTT.1","NTT.2","NTT.3","NTT.4"),type=c(rep("TT",4),rep("NTT",4)))

ann_colors = list(
  type = c(TT="#E1472B",NTT="#73C8DD")

)



ann_colors = list(
  type = c(TT="#AB3282",NTT="#58A4C3")

)

color<-colorRampPalette(rev(c('#FF0000', '#FFD700', '#FFFFFF', '#7B68EE', '#4169E1')))(500)
color <- colorRampPalette(rev(c('#B22222', '#FFD700', '#FFFFFF', '#7B68EE', '#4169E1')))(500)
color <- colorRampPalette(rev(c('#FF6347', '#FFD700', '#FFFFFF', '#7B68EE', '#4169E1')))(500)
color <- colorRampPalette(rev(c('#FF6347', '#FFD700', '#FFFFFF', '#7B68EE', '#6a7fdb')))(500)

color <- colorRampPalette(rev(c('#FF0000', '#FFFF00', '#FFD700', '#FFFFFF', '#7B68EE', '#4169E1')))(500)
# 引入柔和的过渡色，例如浅黄色 (#FFFFCC) 和浅蓝色 (#CCCCFF)
color <- colorRampPalette(rev(c('#FF0000', '#FFFF00', '#FFD700', '#FFFFCC', '#FFFFFF', '#CCCCFF', '#7B68EE', '#4169E1')))(500)
# 定义一个更平滑的颜色渐变
color <- colorRampPalette(rev(c(
  '#FF0000',   # 红色
  '#FF4500',   # 橙红色
  '#FFA500',   # 橙色
  '#FFD700',   # 金色
  '#FFFF00',   # 黄色
  '#FFFFE0',   # 浅黄色
  '#F0FFFF',   # 亚兹特兰的蓝色
  '#B0E0E6',   # 粉蓝色
  '#87CEFA',   # 亮天蓝色
  '#1E90FF',   # 道奇蓝
  '#4169E1'    # 皇家蓝
)))(500)


color <- colorRampPalette(rev(c(
  '#FF0000',   # 红色
  '#FFFF00',   # 黄色
  '#FFD700',   # 金色
  '#FFFFCC',   # 浅黄色
  '#FFFFFF',   # 白色
  '#CCCCFF',   # 淡蓝色
  '#7B68EE',   # 亮紫色
  '#27408B'    # 海军蓝，替代原来的皇家蓝，使其更暗
)))(500)

#############################aaaa
#############################aaaa
color <- colorRampPalette(rev(c(
  '#FF0000',   # 红色
  '#FF4500',   # 橙红色
  '#FFA500',   # 橙色
  '#FFD700',   # 金色
  '#FFFF00',   # 黄色
  '#FFFFE0',   # 浅黄色
  '#F0FFFF',   # 亚兹特兰的蓝色
  '#B0E0E6',   # 粉蓝色
  '#87CEFA',   # 亮天蓝色
  '#1E90FF',   # 道奇蓝
  '#27408B'    # 海军蓝
)))(500)

################################

color <- colorRampPalette(rev(c(
  '#FF0000',   # 红色
  '#FF4500',   # 橙红色
  '#FFA500',   # 橙色
  '#FFD700',   # 金色
  '#FFFF00',   # 黄色
  '#FFFFE0',   # 浅黄色
  '#F0FFFF',   # 亚兹特兰的蓝色
  '#B0E0E6',   # 粉蓝色
  '#87CEFA',   # 亮天蓝色
  '#87CEEB',   # 道奇蓝
  '#89CFF0'    # 婴儿蓝
)))(500)
#
color <- colorRampPalette(rev(c(
  '#FF0000',   # 红色
  '#FF4500',   # 橙红色
  '#FFA500',   # 橙色
  '#FFD700',   # 金色
  '#FFFF00',   # 黄色
  '#FFFFE0',   # 浅黄色
  '#F0FFFF',   # 亚兹特兰的蓝色
  '#B0E0E6',   # 粉蓝色
  '#87CEFA',   # 亮天蓝色
  '#87CEEB',   # 道奇蓝
  '#89CFF0'    # 婴儿蓝
)))(500)
#
color <- colorRampPalette(rev(c(
  '#FF0000',   # 红色
  '#FF4500',   # 橙红色
  '#FFA500',   # 橙色
  '#FFD700',   # 金色
  '#FFFF00',   # 黄色
  '#FFFFE0',   # 浅黄色
  '#F0FFFF',   # 亚兹特兰的蓝色
  '#B0E0E6',   # 粉蓝色
  '#87CEFA',   # 亮天蓝色
  '#89CFF0',   # 道奇蓝
  '#1E90FF'    # 婴儿蓝
)))(500)

pheatmap::pheatmap(cor_res,color = color,border_color = NA,
                   cluster_rows = F,cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   gaps_row = 4,gaps_col = 4,
                   annotation_colors=ann_colors
)





average_data<-AverageExpression(Malignant_epithelial_cells,assays = 'SCT',slot = 'scale.data',group.by="orig.ident")

average_data

View(average_data$SCT
df<-average_data$SCT%>%data.frame()
View(df)
colnames(df)
df1<-df[,c("TT.1","TT.2","TT.3","TT.4","NTT.1","NTT.2","NTT.3","NTT.4")]
cor_res<-cor(df1,method = c("pearson"))
# color hight #A11138   #A31A3B
# color low #6651B4
# color middle white
View(cor_res)
color<-colorRampPalette(c("#6651B4","#B0ECAE","white","#FEF8C6","#FC914B","#A31A3B"))(500)
color<-colorRampPalette(c("#6651B4","white","#A31A3B"))(500)
annotation_col<-data.frame(row.names = c("TT.1","TT.2","TT.3","TT.4","NTT.1","NTT.2","NTT.3","NTT.4"),type=c(rep("TT",4),rep("NTT",4)))
annotation_row<-data.frame(row.names = c("TT.1","TT.2","TT.3","TT.4","NTT.1","NTT.2","NTT.3","NTT.4"),type=c(rep("TT",4),rep("NTT",4)))

ann_colors = list(
  type = c(TT="#E1472B",NTT="#73C8DD")

)



ann_colors = list(
  type = c(TT="#AB3282",NTT="#58A4C3")

)
color<-colorRampPalette(c("#6651B4","#B0ECAE","#FEF8C6","#FC914B","#A31A3B"))(500)
#
color<-colorRampPalette(c("#6651B4","#B0E0E0","#FEF8C6","#FC8B3B","#A31A3B"))(500)
pheatmap::pheatmap(cor_res,color = color,border_color = NA,
                   cluster_rows = F,cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   gaps_row = 4,gaps_col = 4,
                   annotation_colors=ann_colors
                   )



################################################


################################################
################################################
# 默认的颜色，空间分布图，默认的点大小


p1<-SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4,crop = F,cols = color_select)



#################################################

carcinoma_regions_integrated1$orig.ident <- factor(carcinoma_regions_integrated1$orig.ident,levels=unique(carcinoma_regions_integrated1$orig.ident))
DimPlot(carcinoma_regions_integrated1,group.by = "orig.ident") +ggtitle("")




cell_type_bar<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(dominant_celltype_combined_0.1,orig.ident)%>%table()%>%prop.table(margin = 1)%>%data.frame()



cell_type_bar$orig.ident<-factor(cell_type_bar$orig.ident,levels=c(paste0("TT.",c(1:4)),paste0("NTT.",c(1:4))))

cell_type_bar$dominant_celltype_combined_0.1<-factor(cell_type_bar$dominant_celltype_combined_0.1,levels=levels(carcinoma_regions_integrated1))
colnames(cell_type_bar)[1]<-"celltype"
colnames(cell_type_bar)[2]<-"samples"
colnames(cell_type_bar)

unique(cell_type_bar$orig.ident)
names(cols_sample)<-unique(cell_type_bar$orig.ident)


ggplot(data = cell_type_bar,aes(x=celltype,y=Freq,fill=samples))+
  geom_bar(position = 'stack',stat="identity")+
  scale_fill_manual(values = cols_sample)+
  xlab(label = "cluster")+ylab(label = 'proportion(%)')+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))


##################################################### plot


SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4,crop = F,cols = color_select)

integrated_spatial$dominant_celltype_combined_0.1 <-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$dominant_celltype_combined_0.1)

integrated_spatial$dominant_celltype_combined_0.1[grepl("-1",integrated_spatial$dominant_celltype_combined_0.1)]=NA
integrated_spatial$dominant_celltype_combined_0.1[is.na(integrated_spatial$dominant_celltype_combined_0.1)]="Normal_region"

Idents(integrated_spatial)<-integrated_spatial$dominant_celltype_combined_0.1

levels(integrated_spatial)<-c(levels(carcinoma_regions_integrated1),"Normal_region")

color_select1<-c(color_select,"lightgrey")
names(color_select1)[11]<-"Normal_region"

integrated_spatial$cancer_region_and_noncancer_region<-integrated_spatial$dominant_celltype_combined_0.1


integrated_spatial$cancer_region_and_noncancer_region[integrated_spatial$cancer_region_and_noncancer_region!="Normal_region"]="Carcinoma regions"

integrated_spatial$cancer_region_and_noncancer_region[integrated_spatial$cancer_region_and_noncancer_region=="Normal_region"]="Other"


Idents(integrated_spatial)<-integrated_spatial$cancer_region_and_noncancer_region

levels(integrated_spatial)<-c("Carcinoma regions","Other")
color_select2<-c("#833c66","lightgrey")
color_select2<-c("#b32e5b","lightgrey")
names(color_select2)<-levels(integrated_spatial)
SpatialDimPlot(integrated_spatial,ncol = 4,cols = color_select2)



#####################################################
##################################################### 2024-1-15 tf plot






tf_binary<-fread("/result/spatial_auc_mtx_binary.csv")%>%data.frame(row.names = 1,check.names = F)

cols=c("white","black")
pheatmap::pheatmap(tf_binary,color = cols)



pheatmap::pheatmap(auc, fontsize_row = 3,
                   annotation_col = anno,
                   scale = 'none',
                   color=c('white','black'),
                   treeheight_row=20, treeheight_col=20,
                   border_color=NA,cluster_cols = T,
                   show_rownames = T,show_colnames = F,
                   cutree_rows = 4)



tf_auc<-fread("02.pyscenic_analysis/result/spatial_auc_mtx.csv")%>%data.frame(row.names = 1,check.names = F)

malignant_tf<-subset(carcinoma_regions_integrated1,cells=rownames(tf_auc))
SpatialDimPlot(malignant_tf,crop = F,ncol = 4)
########################################################
########################################################
tf_auc_t<-t(tf_auc)

malignant_tf[["TF"]] <- CreateAssayObject(counts = tf_auc_t)

DefaultAssay(malignant_tf)<-'TF'
SpatialFeaturePlot(malignant_tf,features = 'AHR (78g)',crop = F,ncol = 4)
#######################################
#######################################

diff_tt_vs_ntt<-subset(FindMarkers(malignant_tf,group.by = 'group',ident.1 = "TT",ident.2 = 'NTT',logfc.threshold = 0.2),p_val_adj <0.05)

fwrite(diff_tt_vs_ntt,file = "diff_tt_vs_ntt.xls",sep = '\t',col.names = T,row.names = T)


#####################################################2024-1-15
# 热图更改颜色
# 后面更改了细胞类型的颜色 这个热图的颜色忘记更改了


DoHeatmap(subset(carcinoma_regions_integrated1,downsample=500),features = top10$gene,group.colors=color_select,raster=F) +
  scale_fill_gradientn(colors=color_heatmap)+theme(axis.text.y = element_text(size = 3))


# Malignant_epithelial_cells
Malignant_epithelial_cells <- c("CEACAM6", "FABP1", "MUC12")
# B_cells
B_cells <- c("MZB1", "IGHM", "CD79A")
# CD4_T_cells
CD4_T_cells <- c("EIF3K", "HMGB1", "H2AFV")
# CD8_T_cells
CD8_T_cells <- c("CCL5", "CCR7", "GZMA")
# Endothelial_cells
Endothelial_cells <- c("ACKR1", "APLNR", "PECAM1")
# Fibroblast_cells
Fibroblast_cells <- c("EFEMP1", "COL1A1", "DCN")
# Macrophage
Macrophage <- c("CYTH4", "ITGAM", "GBP5")
# Nonmalignant_epithelial_cells
Nonmalignant_epithelial_cells <- c("PLAC8", "PIGR", "CTSE")
# SMC
SMC <- c("DES", "ACTG2", "CNN1")


top10_select=c(Malignant_epithelial_cells,Nonmalignant_epithelial_cells,Fibroblast_cells,
               CD4_T_cells,SMC,Endothelial_cells,B_cells,Macrophage,CD8_T_cells)





colors1=c("#567DBB","white","#DC4A44")
levels(dot_carcinoma_regions_integrated1)<-rev(levels(dot_carcinoma_regions_integrated1))
DotPlot(dot_carcinoma_regions_integrated1, features = unique(top10_select),
        cols = colors1)+scale_colour_gradientn(colours=colors1)+
  theme(axis.text.x = element_text(hjust=1,vjust=1,size = 8))+
  RotatedAxis()



############################################################测试kegg为什么没有富集 2024-1-9
diff_table="02.diff/Nonmalignant_epithelial_cells-TT_vs_NTT_p0.05.xls"
diff_table <- fread(diff_table,data.table=F)


common_gene_id <-bitr(diff_table$V1,fromType = c('SYMBOL'),toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
enrich_kegg <-enrichKEGG(gene = common_gene_id$ENTREZID,organism = "hsa",pvalueCutoff = 0.05,qvalueCutoff = 1,use_internal_data = T)


############################################################ cms  


Malignant_epithelial_cells<-subset(carcinoma_regions_integrated1,idents=c("Malignant_epithelial_cells"))
Malignant_epithelial_cells@commands$FindIntegrationAnchors
DefaultAssay(Malignant_epithelial_cells)<-'integrated'
Malignant_epithelial_cells<-RunUMAP(Malignant_epithelial_cells,dims = 1:30)
Malignant_epithelial_cells<-FindNeighbors(Malignant_epithelial_cells,dims = 1:30)
Malignant_epithelial_cells<-FindClusters(Malignant_epithelial_cells,resolution = 0.2)
DefaultAssay(Malignant_epithelial_cells)<-'integrated'
DefaultAssay(Malignant_epithelial_cells) <- 'SCT'
Malignant_epithelial_cells <- PrepSCTFindMarkers(object = Malignant_epithelial_cells)
sig <- subset(FindAllMarkers(Malignant_epithelial_cells, only.pos = T,assay="SCT"), p_val_adj < 0.05)

cms<-fread('/2024-1-17/01.cms_marker_gene/cms.txt',header = T,data.table = F) %>%data.frame()
cms_list<-as.list(cms)

cms_list<-lapply(cms_list,function(x)x[x!=""])

###############################
Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,cms_list)
View(Malignant_epithelial_cells@meta.data)
colnames(Malignant_epithelial_cells@meta.data)
VlnPlot(Malignant_epithelial_cells,features = 'Cluster1',pt.size = 0)
VlnPlot(Malignant_epithelial_cells,features = 'Cluster2',pt.size = 0)
VlnPlot(Malignant_epithelial_cells,features = 'Cluster3',pt.size = 0)

VlnPlot(Malignant_epithelial_cells,features = 'KRAS',pt.size = 0)
VlnPlot(Malignant_epithelial_cells,features = 'EPCAM',pt.size = 0)
# TGFB1
VlnPlot(Malignant_epithelial_cells,features = 'TGFB1',pt.size = 0)
VlnPlot(Malignant_epithelial_cells,features = 'IL2',pt.size = 0)
VlnPlot(Malignant_epithelial_cells,features = 'PRF1',pt.size = 0) #PRF1
# IFNG
VlnPlot(Malignant_epithelial_cells,features = 'IFNG',pt.size = 0)

########################################


library(xlsx)
sheet_name<-getSheets(loadWorkbook('./geneset_cms.xlsx'))
class(sheet_name)
sheet_name_test<-names(sheet_name)
setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-17/01.cms_marker_gene")


cms_genesets <- lapply(sheet_name_test, function(sheet) {
  print(sheet)
  read.xlsx('./geneset_cms.xlsx', sheetName = sheet,header = F)
})


cms_genesets[[16]]<-NULL





#################################################################
#################################################################
sheet_name_test<-sheet_name_test[-16]
names(cms_genesets)<-sheet_name_test

##################################################################打分
Malignant_epithelial_cells$Cluster1<-NULL
Malignant_epithelial_cells$Cluster2<-NULL
Malignant_epithelial_cells$Cluster3<-NULL
# ?AddModuleScore
View(Malignant_epithelial_cells@meta.data)
View(cms_genesets$STROMAL_INFILTRATION)
cms_genesets1<-lapply(cms_genesets,function(x){as.character(x$X1)})
names(cms_genesets1)<-sheet_name_test
# cms_genesets1

# Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,cms_list,name = sheet_name_test)
Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,cms_genesets1,name = sheet_name_test)

# names(cms_genesets1)
View(Malignant_epithelial_cells@meta.data)
colnames(Malignant_epithelial_cells@meta.data)

p1<-VlnPlot(Malignant_epithelial_cells,features = 'STROMAL_INFILTRATION1',pt.size = 0) # 3群最高  CMS4
p2<-VlnPlot(Malignant_epithelial_cells,features = 'EMT_ACTIVATION_SIGNATURE2',pt.size = 0)  # 3群最高 #
p3<-VlnPlot(Malignant_epithelial_cells,features = 'COMPLEMENT_ACTIVATION3',pt.size = 0) #  3群第二高，1 群最高

(p1|p2)/p3



p1<-VlnPlot(Malignant_epithelial_cells,features = 'CRYPT_TOP_SIGNATURE4',pt.size = 0)  # 0群最高   # CMS3
p2<-VlnPlot(Malignant_epithelial_cells,features = 'SUGAR_AMINO_ACID_NUCLEOTIDE_MET6',pt.size = 0)  # 0群高
VlnPlot(Malignant_epithelial_cells,features = 'TH17_ACTIVATION5',pt.size = 0)  # 0群高  3群高
p3<-VlnPlot(Malignant_epithelial_cells,features = 'FRUCTOSE_MANNOSE_METABOLISM7',pt.size = 0)  # 0群高
p4<-VlnPlot(Malignant_epithelial_cells,features = 'STARCH_SUCROSE_METABOLISM8',pt.size = 0)  # 0群高

p5<-VlnPlot(Malignant_epithelial_cells,features = 'GALACTOSE_METABOLISM9',pt.size = 0)  # 0群高
VlnPlot(Malignant_epithelial_cells,features = 'GOBP_GLUTAMINE_METABOLIC_PROCES10',pt.size = 0)  # 0群高 (没什么区别)
VlnPlot(Malignant_epithelial_cells,features = 'GLUTATHIONE_METABOLISM11',pt.size = 0)  # 0群高 3 群高
p6<-VlnPlot(Malignant_epithelial_cells,features = 'NITROGEN_METABOLISM12',pt.size = 0)  # 0群高

p7<-VlnPlot(Malignant_epithelial_cells,features = 'TYROSINE_METABOLISM13',pt.size = 0)  # 0群高
p8<-VlnPlot(Malignant_epithelial_cells,features = 'FATTY_ACID_METABOLISM14',pt.size = 0)  # 0群高
VlnPlot(Malignant_epithelial_cells,features = 'LINOLEIC_ACID_METABOLISM15',pt.size = 0)  # 0群高  高的不明显

(p1|p2|p3|p4)/(p5|p6|p7|p8)


###############################################
###############################################

cms_genesets2<-fread("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-17/01.cms_marker_gene/geneset_cms1.txt",data.table = F)

cms_genesets2_list<-as.list(cms_genesets2)
cms_genesets2_list<-lapply(cms_genesets2_list,function(x){x[x!=""]})
colnames(Malignant_epithelial_cells@meta.data)
Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,cms_genesets2_list,name = names(cms_genesets2_list))

p1<-VlnPlot(Malignant_epithelial_cells,features = 'NFILTRATION1',pt.size = 0) # 没有最高的  1群cms1
p2<-VlnPlot(Malignant_epithelial_cells,features = 'TFH_INFILTRATION2',pt.size = 0)  # 没有最高的
p3<-VlnPlot(Malignant_epithelial_cells,features = 'NKC_INFILTRATION3',pt.size = 0) #  3群最高
p4<-VlnPlot(Malignant_epithelial_cells,features = 'JAK_STAT_ACTIVATION4',pt.size = 0) #  1和 3群最高
p5<-VlnPlot(Malignant_epithelial_cells,features = 'CYTOTOXIC_T_CELLS_INFILTRATION6',pt.size = 0) #  1群高
VlnPlot(Malignant_epithelial_cells,features = 'activation_of_immune_response10',pt.size = 0) #
p1
(p1|p2)/(p3|p4|p5)


p1<-VlnPlot(Malignant_epithelial_cells,features = 'WNT_ACTIVATION5',pt.size = 0) #  0和2 比较高 2群 cms2
p2<-VlnPlot(Malignant_epithelial_cells,features = 'WNT_ACTIVATION5',pt.size = 0) #  0和2 比较高
p3<-VlnPlot(Malignant_epithelial_cells,features = 'TRANSLATION_ACTIVATION_REACTOME8',pt.size = 0) #  0和2 比较高
p4<-VlnPlot(Malignant_epithelial_cells,features = 'SRC_ACTIVATION_BIOCARTA7',pt.size = 0) #  0和2 比较高
VlnPlot(Malignant_epithelial_cells,features = 'MYC_ACTIVATION9',pt.size = 0) #  0和2 比较高
# MYC_ACTIVATION9
(p1|p2)/(p3|p4)
################################################################################
VlnPlot(Malignant_epithelial_cells,features = 'APC',pt.size = 0)
VlnPlot(Malignant_epithelial_cells,features = 'KRAS',pt.size = 0)

################################################################################
################################################################################ 2024-1-18
################################################################################ 2024-1-18

View(Malignant_epithelial_cells@meta.data)
colnames(Malignant_epithelial_cells@meta.data)

colnames(Malignant_epithelial_cells@meta.data)
df_score<-Malignant_epithelial_cells@meta.data[35:58]
df_score<-Malignant_epithelial_cells@meta.data[35:60]
df_score$cluster<-Malignant_epithelial_cells$seurat_clusters
class(df_score$cluster)

df_score<-df_score%>%dplyr::arrange(cluster)

View(df_score)
anno_row<-data.frame(row.names=rownames(df_score), type=df_score$cluster)
df_score$cluster<-NULL
pheatmap::pheatmap(df_score,scale = 'row',cluster_rows = F,show_rownames = F,annotation_row = anno_row,angle_col = 45)

########################################################
########################################################

########################################################
########################################################
DefaultAssay(Malignant_epithelial_cells)<-'integrated'

score_name<-colnames(Malignant_epithelial_cells@meta.data[35:60])
setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-17/02.4_cluster")
setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-17/jiaofu/05.score_boxplot")
# CMS4
View(Malignant_epithelial_cells@meta.data)
?VlnPlot
colnames(Malignant_epithelial_cells@meta.data)
VlnPlot(Malignant_epithelial_cells,features = "activation_of_immune_response10",pt.size = 0)+geom_boxplot()
dir.create('CMS4')
for(i in score_name[1:3]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)+geom_boxplot()
  ggsave(p1,filename = paste0("CMS4/",i,"_vlnplot.png"),width = 12,height = 6)
}
# CMS3
dir.create('CMS3')
for(i in score_name[4:15]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)
  ggsave(p1,filename = paste0("CMS3/",i,"_vlnplot.png"),width = 12,height = 6)
}

# CMS1
dir.create('CMS1')
for(i in score_name[c(17,18,19,20,22,26)]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)
  ggsave(p1,filename = paste0("CMS1/",i,"_vlnplot.png"),width = 12,height = 6)
}

# CMS2
dir.create('CMS2')
for(i in score_name[c(21,23,24,25)]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)
  ggsave(p1,filename = paste0("CMS2/",i,"_vlnplot.png"),width = 12,height = 6)
}


##########################################################################

Malignant_epithelial_cells<-FindClusters(Malignant_epithelial_cells,resolution = 0.4)
DimPlot(Malignant_epithelial_cells)
score_name<-colnames(Malignant_epithelial_cells@meta.data[35:60])
setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-17/03.7_cluster")

# CMS4
dir.create('CMS4')
for(i in score_name[1:3]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)
  ggsave(p1,filename = paste0("CMS4/",i,"_vlnplot.png"),width = 12,height = 6)
}
# CMS3
dir.create('CMS3')
for(i in score_name[4:15]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)
  ggsave(p1,filename = paste0("CMS3/",i,"_vlnplot.png"),width = 12,height = 6)
}

# CMS1
dir.create('CMS1')
for(i in score_name[c(17,18,19,20,22,26)]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)
  ggsave(p1,filename = paste0("CMS1/",i,"_vlnplot.png"),width = 12,height = 6)
}

# CMS2
dir.create('CMS2')
for(i in score_name[c(21,23,24,25)]){
  p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)
  ggsave(p1,filename = paste0("CMS2/",i,"_vlnplot.png"),width = 12,height = 6)
}

DimPlot(Malignant_epithelial_cells)
Malignant_epithelial_cells<-FindClusters(Malignant_epithelial_cells,resolution = 0.2)
DefaultAssay(Malignant_epithelial_cells)<-'SCT'
######################################################
###################################################### 2024-1-18
###################################################### 2024-1-18

p1<-FeaturePlot(obj,features = "MS4A2")

##########################################################
##########################################################


GO_db <- msigdbr(species="Homo sapiens",category="H") %>%
  dplyr::select(gs_name, entrez_gene, gene_symbol)


########################################
# cms_genesets1
# cms_genesets2_list


combined_cms_genesets_list <-append(cms_genesets2_list,cms_genesets1)


combined_cms_genesets_list$TFH_INFILTRATION<-NULL

result_df <- do.call(rbind, Map(function(term, genes) data.frame(term = rep(term, length(genes)), gene_name = genes),
                                names(combined_cms_genesets_list),
                                combined_cms_genesets_list))

#################################################################
#################################################################
combined_cms_genesets_list$TH1_INFILTRATION
combined_cms_genesets_list$TFH_INFILTRATION
combined_cms_genesets_list$NKC_INFILTRATION
combined_cms_genesets_list$JAK_STAT_ACTIVATION
colnames(result_df)<-c("gs_name","gene_symbol")



rownames(result_df)<- c(1:dim(result_df)[1])


##############################################################
##############################################################
##############################################################

setwd("/05.score_boxplot")

dir.create('CMS4')
for(i in score_name[1:3]){
  # p1<-VlnPlot(Malignant_epithelial_cells,features = i,pt.size = 0)+geom_boxplot()
  data_to_plot <-data.frame(cluster=as.character(Malignant_epithelial_cells@meta.data[,'seurat_clusters']),score=Malignant_epithelial_cells@meta.data[,i])
  p1<-ggplot(data_to_plot, aes(x=cluster,y = score,fill=cluster)) +
    geom_boxplot()

  ggsave(p1,filename = paste0("CMS4/",i,"_vlnplot.png"),width = 12,height = 6)
}
# CMS3
dir.create('CMS3')
for(i in score_name[4:15]){
  data_to_plot <-data.frame(cluster=as.character(Malignant_epithelial_cells@meta.data[,'seurat_clusters']),score=Malignant_epithelial_cells@meta.data[,i])
  p1<-ggplot(data_to_plot, aes(x=cluster,y = score,fill=cluster)) +
    geom_boxplot()
  ggsave(p1,filename = paste0("CMS3/",i,"_vlnplot.png"),width = 12,height = 6)
}

# CMS1
dir.create('CMS1')
for(i in score_name[c(17,18,19,20,22,26)]){
  data_to_plot <-data.frame(cluster=as.character(Malignant_epithelial_cells@meta.data[,'seurat_clusters']),score=Malignant_epithelial_cells@meta.data[,i])
  p1<-ggplot(data_to_plot, aes(x=cluster,y = score,fill=cluster)) +
    geom_boxplot()
  ggsave(p1,filename = paste0("CMS1/",i,"_vlnplot.png"),width = 12,height = 6)
}

# CMS2
dir.create('CMS2')
for(i in score_name[c(21,23,24,25)]){
  data_to_plot <-data.frame(cluster=as.character(Malignant_epithelial_cells@meta.data[,'seurat_clusters']),score=Malignant_epithelial_cells@meta.data[,i])
  p1<-ggplot(data_to_plot, aes(x=cluster,y = score,fill=cluster)) +
    geom_boxplot()
  ggsave(p1,filename = paste0("CMS2/",i,"_vlnplot.png"),width = 12,height = 6)
}

#############################################################
#############################################################
#############################################################
df_stat<-Malignant_epithelial_cells@meta.data%>%dplyr::select(seurat_clusters,group)%>%table()%>%prop.table(margin = 2)%>%data.frame()



ggplot(df_stat,aes(x=group,y=Freq,fill=seurat_clusters))+geom_bar(stat = 'identity')


df_stat<-Malignant_epithelial_cells@meta.data%>%dplyr::select(seurat_clusters,orig.ident)%>%table()%>%prop.table(margin = 2)%>%data.frame()


ggplot(df_stat,aes(x=orig.ident,y=Freq,fill=seurat_clusters))+geom_bar(stat = 'identity')


###########################################2024-1-25
# carcinoma_regions_integrated1
SpatialDimPlot(carcinoma_regions_integrated1,ncol = 4)

sig_marker_spatial<-subset(FindAllMarkers(carcinoma_regions_integrated1,only.pos = T),p_val_adj <0.05)


###########################################

sig_marker_spatial_subset<-subset(sig_marker_spatial,cluster%in%c("Malignant_epithelial_cells","Macrophage"))

sig_marker_spatial_subset_top20<-sig_marker_spatial_subset%>%group_by(cluster)%>%top_n(n=20,wt=avg_log2FC)


sig_marker_spatial_subset_top20_genes<-unique(sig_marker_spatial_subset_top20$gene)

sig_marker_spatial_subset_top20_genes_list<-list(sig_marker_spatial_subset_top20_genes=sig_marker_spatial_subset_top20_genes)

carcinoma_regions_integrated1<-AddModuleScore(carcinoma_regions_integrated1,features = sig_marker_spatial_subset_top20_genes_list,name = "sig_marker_spatial_subset_top20_genes_list")

SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'sig_marker_spatial_subset_top20_genes_list1',crop = F,ncol = 4)




integrated_spatial<-readRDS('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-11/data/integrated_spatial.rds')

integrated_spatial$sig_marker_spatial_subset_top20_genes_list <- mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$sig_marker_spatial_subset_top20_genes_list1)


View(integrated_spatial@meta.data)
integrated_spatial$sig_marker_spatial_subset_top20_genes_list<-as.numeric(integrated_spatial$sig_marker_spatial_subset_top20_genes_list)

theme_set(theme_minimal(base_size = 14, base_family = "Helvetica", na.rm = TRUE))

p1<-SpatialFeaturePlot(integrated_spatial,features = 'sig_marker_spatial_subset_top20_genes_list',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")


library(scales)
library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
FeaturePalettes <- list(
  'Spatial' = SpatialColors(n = 100),
  'Seurat' = c('lightgrey', 'blue')
)


################################################# 2024-1-29
# 

geneset<-read.xlsx('/share/home/hhhong/personality/2023-8-8_jingdi_wuhanxiangya/2024-1-29/data/gene_list1.xlsx',sheetIndex = 1,header = F)

geneset$X2<-toupper(geneset$X2)

####################################
geneset1<-geneset%>%dplyr::select(c(1,2,5:515))
rownames(geneset1)<-geneset1$X2
duplicated(geneset1$X2)
geneset1$X1<-paste0(geneset1$X1,"_",geneset1$X2)
rownames(geneset1)<-geneset1$X1
geneset1$X2<-NULL
geneset2<-t(geneset1)
geneset2<-geneset2[-1,]

geneset2<-geneset2%>%data.frame(check.names = F,check.rows = F)
geneset2_list=as.list(geneset2)
geneset2_list1<-lapply(geneset2_list,function(x){x[!is.na(x)]})
SpatialDimPlot(Malignant_epithelial_cells,crop = F,ncol=4)

Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,features = geneset2_list1,name = "geneset")

SpatialFeaturePlot(Malignant_epithelial_cells,features = 'sig_marker_spatial_subset_top20_genes_list1',crop = F,ncol = 4)


colnames(Malignant_epithelial_cells@meta.data)[62:97]=names(geneset2_list1)
Malignant_epithelial_cells$seurat_clusters
Malignant_epithelial_cells$cms_type<-factor(Malignant_epithelial_cells$cms_type,levels=c("CMS1","CMS2","CMS3","CMS4"))
DimPlot(Malignant_epithelial_cells,group.by = 'cms_type',label = T)+ggtitle("")
mapping_df <- c("0","1","2","3")
maping_df_type<-c("CMS2","CMS1","CMS3","CMS4")
Malignant_epithelial_cells$cms_type<-mapvalues(Malignant_epithelial_cells$seurat_clusters,mapping_df,maping_df_type)
score_df<-Malignant_epithelial_cells@meta.data%>%dplyr::select(c(62:97))%>%data.frame(check.names = F,check.rows = F)


# DimPlot(Malignant_epithelial_cells,group.by = "cms_type",label = T)


score_df_stat<-Malignant_epithelial_cells@meta.data%>%dplyr::select(c(62:98))

score_df_stat1<-score_df_stat%>%dplyr::group_by(cms_type)%>%summarise(across(everything(),mean,na.rm=TRUE))%>%data.frame(check.rows = F,check.names = F)
rownames(score_df_stat1)<-score_df_stat1$cms_type
score_df_stat1$cms_type<-NULL
score_df_stat1<-score_df_stat1%>%t()%>%data.frame(check.rows = F,check.names = F)
score_df_stat1<-score_df_stat1[,c("CMS1","CMS2","CMS3","CMS4")]
rownames(score_df_stat1)[10]<-"CMS2_RESPIRATORY_ELECTRON_TRANSPORT-ATP_SYNTHESIS"
pheatmap::pheatmap(score_df_stat1,cluster_rows = F,cluster_cols = F,scale = 'row',border_color = NA)

#####################################################
#####################################################
cms2<-c("EPCAM","CDH1","KRT19","GGH")
cms1<-c("IGKC","IGHG3","IGHA1","IGHG1","IGHM","IGLC1","IFI6")
cms3<-c("ATP5IF1","ATP2C2","ATPSCKMT")
cms4<-c("COL1A1","LUM","SPARC","COL1A2","MMP2")
Idents(Malignant_epithelial_cells)<-Malignant_epithelial_cells$cms_type

gene_aver<-AverageExpression(Malignant_epithelial_cells)

gene_aver1<-gene_aver$SCT
View(gene_aver1)

gene_aver2<-gene_aver1[c(cms1,cms2,cms3,cms4),]



pheatmap::pheatmap(gene_aver2,cluster_rows = F,cluster_cols = F,scale = 'row',border_color = NA)


color_cms<-c("#FFD700","#53A85F","#57C3F3","#AB3282")
levels(Malignant_epithelial_cells)
DimPlot(Malignant_epithelial_cells,cols = color_cms,label = T) +
  guides(colour=guide_legend(title="CMS cluster",override.aes=list(size=3)))


sig_marker<-subset(FindAllMarkers(Malignant_epithelial_cells,only.pos = T),p_val_adj <0.05)
View(sig_marker)
top10<-sig_marker%>%dplyr::group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
#################################################################

Malignant_epithelial_cells<-SCTransform(Malignant_epithelial_cells,assay = "Spatial",variable.features.n=8000)
DoHeatmap(subset(Malignant_epithelial_cells,downsample=500),group.colors = color_cms,features = top10$gene)


################################################################绘制空间图片

SpatialDimPlot(integrated_spatial,ncol = 4)

integrated_spatial$cms_type<-mapvalues(colnames(integrated_spatial),colnames(Malignant_epithelial_cells),as.character(Malignant_epithelial_cells$cms_type))


integrated_spatial$cms_type[grepl("-1",integrated_spatial$cms_type)]=NA


Idents(integrated_spatial)<-integrated_spatial$cms_type
color_cms1<-c("#FFD700","#53A85F","#57C3F3","#AB3282","lightgrey")
color_cms1<-c("#FFD700","#53A85F","#57C3F3","#AB3282")

Idents(integrated_spatial)<-integrated_spatial$cms_type
levels(integrated_spatial)<-c("CMS1","CMS2","CMS3","CMS4","NA")
SpatialDimPlot(integrated_spatial,ncol = 4,crop = F,cols = color_cms1) & scale_fill_manual(values = color_cms1,breaks = c("CMS1","CMS2","CMS3","CMS4"),na.value = 'lightgrey')
########################################


statistic1<-Malignant_epithelial_cells@meta.data%>%dplyr::select(orig.ident,cms_type)%>%group_by(orig.ident)%>%
  table()%>%prop.table(margin = 1)%>%data.frame(check.names = F,check.rows = F)
rm(statistic1)
View(statistic1)
##########################
sum(statistic1$Freq[statistic1$orig.ident=='TT.2'])
colnames(statistic1)<-c("Samples","CMS_type","Proportion")
color_cms1<-c("#FFD700","#53A85F","#57C3F3","#AB3282")
ggplot(data = statistic1,aes(x=Samples,y=Proportion,fill=CMS_type))+
  geom_bar(stat = 'identity',position = 'stack') +
  scale_fill_manual(values = color_cms1) +
  theme_classic() +guides(fill=guide_legend(title = "CMS type"))+
  theme(legend.position = 'bottom') +
  xlab(label = "")

########################## enrichment


geneset<-read.xlsx('/data/gene_list_update.xlsx',sheetIndex = 1,header = F)

geneset$X2<-toupper(geneset$X2)

####################################
geneset1<-geneset%>%dplyr::select(c(1,2,3,5:515))
dim(geneset)


geneset1$X2[duplicated(geneset1$X2)]

geneset1$X1<-paste0(geneset1$X3,"_",geneset1$X2)


rownames(geneset1)<-geneset1$X1
geneset1$X1<-NULL
geneset1$X2<-NULL
geneset1$X3<-NULL
geneset2<-t(geneset1)

class(geneset2)
geneset2<-geneset2%>%data.frame(check.names = F,check.rows = F)
geneset2_list=as.list(geneset2)
geneset2_list1<-lapply(geneset2_list,function(x){x[!is.na(x)]})


Malignant_epithelial_cells@meta.data<-Malignant_epithelial_cells@meta.data[,-c(35:(ncol(Malignant_epithelial_cells@meta.data)-1))]
SpatialDimPlot(Malignant_epithelial_cells,crop = F,ncol=4)

Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,features = geneset2_list1,name = "geneset")



mapping_df <- c("0","1","2","3")
maping_df_type<-c("CMS2","CMS1","CMS3","CMS4")
Malignant_epithelial_cells$cms_type<-mapvalues(Malignant_epithelial_cells$seurat_clusters,mapping_df,maping_df_type)


colnames(Malignant_epithelial_cells@meta.data)[36:66]=names(geneset2_list1)
Idents(Malignant_epithelial_cells)<-Malignant_epithelial_cells$cms_type
Malignant_epithelial_cells$cms_type<-factor(Malignant_epithelial_cells$cms_type,levels=c("CMS1","CMS2","CMS3","CMS4"))
DimPlot(Malignant_epithelial_cells,group.by = 'cms_type',label = T)+ggtitle("")
mapping_df <- c("0","1","2","3")
maping_df_type<-c("CMS2","CMS1","CMS3","CMS4")
Malignant_epithelial_cells$cms_type<-mapvalues(Malignant_epithelial_cells$seurat_clusters,mapping_df,maping_df_type)
score_df<-Malignant_epithelial_cells@meta.data%>%dplyr::select(c(35:66))%>%data.frame(check.names = F,check.rows = F)




score_df_stat1<-score_df%>%dplyr::group_by(cms_type)%>%summarise(across(everything(),mean,na.rm=TRUE))%>%data.frame(check.rows = F,check.names = F)
rownames(score_df_stat1)<-score_df_stat1$cms_type
score_df_stat1$cms_type<-NULL
score_df_stat1<-score_df_stat1%>%t()%>%data.frame(check.rows = F,check.names = F)
score_df_stat1<-score_df_stat1[,c("CMS1","CMS2","CMS3","CMS4")]

rownames(score_df_stat1)[19]<-"REACTOME_RESPIRATORY_ELECTRON_TRANSPORT-ATP_SYNTHESIS"
#
score_df_stat1))
pheatmap::pheatmap(score_df_stat1,cluster_rows = F,cluster_cols = F,scale = 'row',border_color = NA)


########################################################2024-1-31
library("Nebulosa")
genes <- c("EPCAM", "CDH1", "IGHG3", "IGHM", "ATP5IF1", "ATPSCKMT", "COL1A1", "LUM","KRT19","GGH")

DimPlot(Malignant_epithelial_cells,label = T)
# p4 <- plot_density(Malignant_epithelial_cells, genes, joint = TRUE)
# p4 + plot_layout(ncol = 4)

options(scipen = -1)
options(digits = 5)
p4 <- plot_density(Malignant_epithelial_cells, genes, joint = F,size = 2,reduction = 'umap')
p4 + plot_layout(ncol = 4)
Idents(Malignant_epithelial_cells)<-Malignant_epithelial_cells$cms_type
DimPlot(Malignant_epithelial_cells,pt.size = 2,cols = color_cms1)

p4 <- plot_density(Malignant_epithelial_cells, genes, joint = F,size = 0.2,reduction = 'umap')
p4 + plot_layout(ncol = 4)
########################################堆积柱状体


statistic1<-Malignant_epithelial_cells@meta.data%>%dplyr::select(group,cms_type)%>%group_by(group)%>%
  table()%>%prop.table(margin = 1)%>%data.frame(check.names = F,check.rows = F)
rm(statistic1)
View(statistic1)
##########################
sum(statistic1$Freq[statistic1$orig.ident=='TT.2'])
statistic1$Group<-factor(statistic1$Group,levels = c('TT','NTT'))
colnames(statistic1)<-c("Group","CMS_type","Proportion")
color_cms1<-c("#FFD700","#53A85F","#57C3F3","#AB3282")
ggplot(data = statistic1,aes(x=Group,y=Proportion,fill=CMS_type))+
  geom_bar(stat = 'identity',position = 'stack') +
  scale_fill_manual(values = color_cms1) +
  theme_classic() +guides(fill=guide_legend(title = "CMS type"))+
  theme(legend.position = 'bottom') +
  xlab(label = "")



########################## spatial expression
setwd('/2024-1-31/03.spatial_expr')
for(i in genes){
  p1<-SpatialFeaturePlot(Malignant_epithelial_cells,crop = F,features = i,ncol = 4)
ggsave(p1,filename = paste0(i,"_spatial_expr.pdf"),width = 30,height = 25)
}



######################################### correlation

levels(Malignant_epithelial_cells$cms_type)
levels(Malignant_epithelial_cells$orig.ident)
Malignant_epithelial_cells$correlation<-paste0(Malignant_epithelial_cells$cms_type,"_",Malignant_epithelial_cells$orig.ident)
Idents(Malignant_epithelial_cells)<-Malignant_epithelial_cells$correlation
levels(Malignant_epithelial_cells)<-c("CMS1_TT.1", "CMS2_TT.1", "CMS3_TT.1", "CMS4_TT.1",
                                      "CMS1_TT.2", "CMS2_TT.2", "CMS3_TT.2", "CMS4_TT.2",
                                      "CMS1_TT.3", "CMS2_TT.3", "CMS3_TT.3", "CMS4_TT.3",
                                      "CMS1_TT.4", "CMS2_TT.4", "CMS3_TT.4", "CMS4_TT.4",
                                      "CMS1_NTT.1", "CMS2_NTT.1", "CMS3_NTT.1", "CMS4_NTT.1",
                                      "CMS1_NTT.2", "CMS2_NTT.2", "CMS3_NTT.2", "CMS4_NTT.2",
                                      "CMS1_NTT.3", "CMS2_NTT.3", "CMS3_NTT.3", "CMS4_NTT.3",
                                      "CMS1_NTT.4", "CMS2_NTT.4", "CMS3_NTT.4", "CMS4_NTT.4")


# cms_type
DimPlot(integrated_spatial,group.by = 'cms_type')
Malignant_epithelial_cells_cms$cms_type<-mapvalues(colnames(Malignant_epithelial_cells_cms),rownames(integrated_spatial@meta.data),integrated_spatial$cms_type)
DimPlot(Malignant_epithelial_cells_cms,group.by = 'cms_type')

Idents(Malignant_epithelial_cells_cms)<-Malignant_epithelial_cells_cms$orig.ident

###########################开始分析
corr_expr<-AverageExpression(Malignant_epithelial_cells_cms)
corr_expr_sct<-corr_expr$SCT
df<-corr_expr_sct%>%data.frame()


cor_res<-cor(df,method = c("pearson"))
#####################################



cor_res<-cor(df,method = c("pearson"))

# 第一种图
color<-colorRampPalette(c("#6651B4","#B0ECAE","white","#FEF8C6","#FC914B","#A31A3B"))(500)
color<-colorRampPalette(c("#6651B4","white","#A31A3B"))(500)
annotation_col<-data.frame(row.names = col_names,cms_type=gsub("_.*","",col_names),type=group)
annotation_row<-data.frame(row.names = col_names,cms_type=gsub("_.*","",col_names),type=group)

ann_colors = list(
  type = c(TT="#E1472B",NTT="#73C8DD")

)



ann_colors = list(
  type = c(TT="#AB3282",NTT="#58A4C3")

)
color<-colorRampPalette(c("#6651B4","#B0ECAE","#FEF8C6","#FC914B","#A31A3B"))(500)
#
color<-colorRampPalette(c("#6651B4","#B0E0E0","#FEF8C6","#FC8B3B","#A31A3B"))(500)
pheatmap::pheatmap(cor_res,color = color,border_color = NA,
                   cluster_rows = F,cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   gaps_row = 4,gaps_col = 4,
                   annotation_colors=ann_colors
)










# second plot

ann_colors = list(
  type = c(TT="#E1472B",NTT="#73C8DD")

)
View(cor_res)
order_col_name <-c("CMS1_TT.1", "CMS1_TT.2", "CMS1_TT.3", "CMS1_TT.4",
                   "CMS2_TT.1", "CMS2_TT.2", "CMS2_TT.3", "CMS2_TT.4",
                   "CMS3_TT.1", "CMS3_TT.2", "CMS3_TT.3", "CMS3_TT.4",
                   "CMS4_TT.1", "CMS4_TT.2", "CMS4_TT.3", "CMS4_TT.4",
                   "CMS1_NTT.1", "CMS1_NTT.2", "CMS1_NTT.3", "CMS1_NTT.4",
                   "CMS2_NTT.1", "CMS2_NTT.2", "CMS2_NTT.3", "CMS2_NTT.4",
                   "CMS3_NTT.1", "CMS3_NTT.2", "CMS3_NTT.3", "CMS3_NTT.4",
                   "CMS4_NTT.1", "CMS4_NTT.2", "CMS4_NTT.3", "CMS4_NTT.4")

cor_res1<-cor_res[,order_col_name]
cor_res1<-cor_res1[order_col_name,]
group<-gsub("CMS[1-4]+_","",order_col_name)
group<-gsub("\\..*","",group)
color<-colorRampPalette(c("#6651B4","#B0ECAE","white","#FEF8C6","#FC914B","#A31A3B"))(500)
color<-colorRampPalette(c("#6651B4","white","#A31A3B"))(500)
annotation_col<-data.frame(row.names = order_col_name,cms_type=gsub("_.*","",order_col_name),group=group)
annotation_row<-data.frame(row.names = order_col_name,cms_type=gsub("_.*","",order_col_name),group=group)



ann_colors = list(
  group = c(TT="#AB3282",NTT="#58A4C3"),
  cms_type=c(CMS1=color_cms1[1],CMS2=color_cms1[2],CMS3=color_cms1[3],CMS4=color_cms1[4])

)
color<-colorRampPalette(c("#6651B4","#B0ECAE","#FEF8C6","#FC914B","#A31A3B"))(500)
#
color<-colorRampPalette(c("#6651B4","#B0E0E0","#FEF8C6","#FC8B3B","#A31A3B"))(500)
pheatmap::pheatmap(cor_res1,color = color,border_color = NA,
                   cluster_rows = F,cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   gaps_row = 16,gaps_col = 16,
                   annotation_colors=ann_colors
)

# third plot



ann_colors = list(
  type = c(TT="#E1472B",NTT="#73C8DD")

)
View(cor_res)
order_col_name <-c("CMS1_TT.1", "CMS1_TT.2", "CMS1_TT.3", "CMS1_TT.4",
                   "CMS1_NTT.1", "CMS1_NTT.2", "CMS1_NTT.3", "CMS1_NTT.4",
                   "CMS2_TT.1", "CMS2_TT.2", "CMS2_TT.3", "CMS2_TT.4",
                   "CMS2_NTT.1", "CMS2_NTT.2", "CMS2_NTT.3", "CMS2_NTT.4",
                   "CMS3_TT.1", "CMS3_TT.2", "CMS3_TT.3", "CMS3_TT.4",
                   "CMS3_NTT.1", "CMS3_NTT.2", "CMS3_NTT.3", "CMS3_NTT.4",
                   "CMS4_TT.1", "CMS4_TT.2", "CMS4_TT.3", "CMS4_TT.4",
                   "CMS4_NTT.1", "CMS4_NTT.2", "CMS4_NTT.3", "CMS4_NTT.4")


cor_res1<-cor_res[,order_col_name]
cor_res1<-cor_res1[order_col_name,]
group<-gsub("CMS[1-4]+_","",order_col_name)
group<-gsub("\\..*","",group)
color<-colorRampPalette(c("#6651B4","#B0ECAE","white","#FEF8C6","#FC914B","#A31A3B"))(500)
color<-colorRampPalette(c("#6651B4","white","#A31A3B"))(500)
annotation_col<-data.frame(row.names = order_col_name,cms_type=gsub("_.*","",order_col_name),group=group)
annotation_row<-data.frame(row.names = order_col_name,cms_type=gsub("_.*","",order_col_name),group=group)



ann_colors = list(
  type = c(TT="#AB3282",NTT="#58A4C3")

)
color<-colorRampPalette(c("#6651B4","#B0ECAE","#FEF8C6","#FC914B","#A31A3B"))(500)
#
color<-colorRampPalette(c("#6651B4","#B0E0E0","#FEF8C6","#FC8B3B","#A31A3B"))(500)
pheatmap::pheatmap(cor_res1,color = color,border_color = NA,
                   cluster_rows = F,cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   gaps_row = 16,gaps_col = 16,
                   annotation_colors=ann_colors
)



##################################################### spatial plot



genes_df<-Malignant_epithelial_cells@assays$SCT[genes,]
genes_df1<-genes_df%>%t()%>%data.frame(check.names = F,check.rows = F)
genes_df1$barcode <- rownames(genes_df1)
colnames(genes_df1)[ncol(genes_df1)]<-"barcode"


colnames(genes_df1)<-paste0(colnames(genes_df1),"_expr")

integrated_spatial@meta.data<-integrated_spatial@meta.data[,-c(23:33)]
integrated_spatial$barcode <- colnames(integrated_spatial)
integrated_spatial@meta.data<-left_join(integrated_spatial@meta.data,genes_df1,by=c("barcode"))

rownames(integrated_spatial@meta.data)<-integrated_spatial$barcode




setwd('/2024-1-31/03.spatial_expr')
genes_expr<-paste0(genes,"_expr")
for(i in genes_expr){
  p1<-SpatialFeaturePlot(integrated_spatial,crop = F,features = i,ncol = 4)  & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
  ggsave(p1,filename = paste0(i,"_spatial_expr_adjust.pdf"),width = 30,height = 25)
}



###########################################geneset score 2024-2-1
carcinoma_regions_integrated1<-readRDS('/data/carcinoma_regions_integrated1.rds')

sig_marker_spatial<-subset(FindAllMarkers(carcinoma_regions_integrated1,only.pos = T),p_val_adj <0.05)


sig_marker_spatial_subset<-subset(sig_marker_spatial,cluster%in%c("Endothelial_cells","Macrophage"))

sig_marker_spatial_subset_top20<-sig_marker_spatial_subset%>%group_by(cluster)%>%top_n(n=20,wt=avg_log2FC)


sig_marker_spatial_subset_top20_genes<-unique(sig_marker_spatial_subset_top20$gene)

sig_marker_spatial_subset_top20_genes_list<-list(sig_marker_spatial_subset_top20_genes=sig_marker_spatial_subset_top20_genes)

carcinoma_regions_integrated1<-AddModuleScore(carcinoma_regions_integrated1,features = sig_marker_spatial_subset_top20_genes_list,name = "sig_marker_spatial_subset_top20_genes_list_endo_cell_macrophage")

SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'sig_marker_spatial_subset_top20_genes_list1',crop = F,ncol = 4)



integrated_spatial$sig_marker_spatial_subset_top20_genes_list_endo_cell_macrophage <- mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$sig_marker_spatial_subset_top20_genes_list_endo_cell_macrophage1)



integrated_spatial$sig_marker_spatial_subset_top20_genes_list_endo_cell_macrophage<-as.numeric(integrated_spatial$sig_marker_spatial_subset_top20_genes_list_endo_cell_macrophage)


colnames(integrated_spatial@meta.data)[34]="Macrophage_and_Endothelial_cells_top20_score"
p1<-SpatialFeaturePlot(integrated_spatial,features = 'Macrophage_and_Endothelial_cells_top20_score',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")


#################################
################################# diff heatmap
# macrophage
gene_set1<-c("COMP/RAP2B/ACTG1/CSRP1/DGKH/SERPINE1/PPIA/SERPING1/SYK/LMAN1/GNAQ/JAK2/F13A1/PRDX2/THBS1/FAP/FCER1G/DGKE/AP3B1/SERPINA1/JMJD1C/CD9/PIK3CA/LCK/PLCG2/LNPK/DGKZ/STXBP1/FERMT3/PTPN6/VKORC1/PLEK")
gene_set1<-c("JUNB/BTG1/CEBPB/DUSP1/CCND1/BHLHE40/TRIP10/RELB/CCL5/BTG2/SAT1/KLF2/CD44/PNRC1/SERPINE1/TNFAIP3/CFLAR/SPSB1/ZFP36/MSC/MAP3K8/LITAF/MCL1/RHOB/CD83/LDLR/VEGFA/FOSL2/PER1/GPR183/TANK/IL1B/CXCL1/NINJ1/ETS2/HES1/MARCKS/EGR1/REL/NFIL3/IER5/B4GALT5/ACKR3/PLPP3/ABCA1/FOS/KLF6/EFNA1/TNIP2/JAG1/IER2/PLEK")
gene_set1<-c("IGHG3，TMSB10，IGKC，IGHG1，CEACAM6，IGLC1，MGP，GGH，JCHAIN，SPP1")
# endo
gene_set1 <-c("IGHG3，IGKC，TMSB10，IGHG1，IGHA1，CEACAM6，MGP，PIM2，CFL1，HIST1H1C，CXCL8")
gene_set1<-c("PARK7/COX5A/COX5B/NDUFA13/UQCRQ/NDUFS5/AKR7A3/COX6A1/NDUFB8/ETFRF1/NDUFA4/GPD2/AKR1B1/COX7C/CYCS/COX6B1/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX11/COX4I1/PLEC")
gene_set1<-c("PARK7/COX5A/COX5B/UQCRQ/NDUFS5/COX6A1/NDUFB8/ETFRF1/NDUFA4/GPD2/COX7C/CYCS/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX4I1/PLEC")
gene_set1<-c("PARK7/COX5A/COX5B/UQCRQ/NDUFS5/COX6A1/NDUFB8/NDUFA4/COX7C/CYCS/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX4I1")
gene_set1<-strsplit(gene_set1,"/")[[1]]
gene_set1<-strsplit(gene_set1,"，")[[1]]

####################################################



Idents(Malignant_epithelial_cells)<-Malignant_epithelial_cells$orig.ident
expr1<-AverageExpression(Malignant_epithelial_cells)
expr <-expr1$SCT
expr<-expr[gene_set1,]
# pheatmap::pheatmap(expr,scale = 'column')
color1<-colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)
pheatmap::pheatmap(expr,scale = 'row',
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-2,2,length.out=101))



########################################################2024-2-13
Malignant_epithelial_cells<-readRDS('/data/Malignant_epithelial_cells.rds')

########################################################
cancer_marker<-fread("/data/cancer_marker.xls",header = F)

colnames(cancer_marker)<-"CRC_Cancer_marker"
cancer_marker<-as.list(cancer_marker)

##############scoring
Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,features = cancer_marker,name = "CRC_Cancer_marker")


##########
cancer_marker_score_df <- Malignant_epithelial_cells@meta.data%>%
  dplyr::select(orig.ident,group,CRC_Cancer_marker1)%>%group_by(orig.ident)%>%
  dplyr::mutate(mean_var=mean(CRC_Cancer_marker1,na.rm = TRUE),sd_var=sd(CRC_Cancer_marker1,na.rm = TRUE))%>%
  ungroup()


cancer_marker_score_df$group<-factor(cancer_marker_score_df$group,levels=c("TT",'NTT'))
colnames(cancer_marker_score_df)[1]<-'samples'
colnames(cancer_marker_score_df)[4]<-'CRC_signature_score'
# scatter
ggplot(data = cancer_marker_score_df,aes(x=group,y=CRC_signature_score,color=samples))+
  geom_point(stat="identity",position = position_dodge(width = 0.75))+
  geom_errorbar(aes(ymin=CRC_signature_score -sd_var,ymax=CRC_signature_score+sd_var),width=0.2,position = position_dodge(0.75)) +
  theme_classic2()+ylab(label = "CRC signature score")

################
################

spatial_list_object<-readRDS('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-2-13/data/spatial_list_object.rds')
integrated_spatial<-spatial_list_object$integrated_spatial


###########################################
integrated_spatial$CRC_score<-mapvalues(colnames(integrated_spatial),colnames(Malignant_epithelial_cells),Malignant_epithelial_cells$CRC_Cancer_marker1)


View(integrated_spatial@meta.data)
integrated_spatial$CRC_score <-as.numeric(integrated_spatial$CRC_score)

SpatialFeaturePlot(integrated_spatial,features = 'CRC_score',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")

View(integrated_spatial@meta.data)


################################
# 巨噬细胞中下列基因的热图展示
genelist = c("FSTL3、CD274、IL1B、IL10RA、SOCS3、STAT3、CEBPB、BHLHE40、A2M、C1QA、LGMN、C3、THBS1、TIMP1")
genelist <-strsplit(genelist,"、")[[1]]
# 内皮细胞
genelist1 = c("CXCL8、CXCL1、CXCL10、CXCL14、CXCL6、IL11、IL4I1、ICAM1、VCAM1")
genelist1 <-strsplit(genelist1,"、")[[1]]

carcinoma_regions_integrated1<-spatial_list_object$carcinoma_regions_integrated1

carcinoma_regions_integrated1_Endothelial_cells<-subset(carcinoma_regions_integrated1,idents = "Endothelial_cells")
carcinoma_regions_integrated1_Macrophage<-subset(carcinoma_regions_integrated1,idents = "Macrophage")
Idents(carcinoma_regions_integrated1_Macrophage)<-carcinoma_regions_integrated1_Macrophage$orig.ident

Idents(carcinoma_regions_integrated1_Endothelial_cells)<-carcinoma_regions_integrated1_Endothelial_cells$orig.ident
expr1<-AverageExpression(carcinoma_regions_integrated1_Macrophage)
expr1<-AverageExpression(carcinoma_regions_integrated1_Endothelial_cells)
expr <-expr1$SCT
expr<-expr[genelist,]
expr<-expr[genelist1,]
# pheatmap::pheatmap(expr,scale = 'column')
color1<-colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)
pheatmap::pheatmap(expr,scale = 'row',
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-2,2,length.out=101))



########################################2024-2-20
# 这个结果图里面TT和NTT之间有无差异？ 上方加个显著性标记
library(ggsignif)

p1<-ggplot(data = cancer_marker_score_df,aes(x=group,y=CRC_signature_score,color=samples))+
  geom_point(stat="identity",position = position_dodge(width = 0.75))+
  geom_errorbar(aes(ymin=CRC_signature_score -sd_var,ymax=CRC_signature_score+sd_var),width=0.2,position = position_dodge(0.75)) +
  theme_classic2()+ylab(label = "CRC signature score")+
  # geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = T)
geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)
p1

p_value=2.2e-16
p1 + annotate("text", x = 1.8, y = 0.25,
                  label = paste("p =", formatC(p_value, format = "e", digits = 2)),
                  color = "red", size = 4, vjust = -2)



tt_score<-subset(cancer_marker_score_df,group=='TT')
ntt_score<-subset(cancer_marker_score_df,group=='NTT')
tt_vs_ntt<-t.test(tt_score$CRC_signature_score,ntt_score$CRC_signature_score)



############################################


SpatialFeaturePlot(integrated_spatial,features = 'CRC_score',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey",limits=c(NA,0.1))



DefaultAssay(integrated_spatial)<-'SCT'
SpatialFeaturePlot(integrated_spatial,features = 'PROM1',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")

####################################################
####################################################
integrated_spatial@assays$SCT@scale.data[1:10,1:10]
integrated_spatial@assays$SCT@data[1:10,1:10]

integrated_spatial$PROM1 <-integrated_spatial@assays$SCT@scale.data['PROM1',]
integrated_spatial$CD44 <-integrated_spatial@assays$SCT@scale.data['CD44',]
integrated_spatial$LGR5 <-integrated_spatial@assays$SCT@scale.data['LGR5',]
integrated_spatial$ALDH1A1 <-integrated_spatial@assays$SCT@scale.data['ALDH1A1',]
integrated_spatial$EPCAM <-integrated_spatial@assays$SCT@scale.data['EPCAM',]

View(integrated_spatial@meta.data)
integrated_spatial$PROM1[which(is.na(integrated_spatial$CRC_score))]=NA
integrated_spatial$CD44[which(is.na(integrated_spatial$CRC_score))]=NA
integrated_spatial$LGR5[which(is.na(integrated_spatial$CRC_score))]=NA
integrated_spatial$ALDH1A1[which(is.na(integrated_spatial$CRC_score))]=NA
integrated_spatial$EPCAM[which(is.na(integrated_spatial$CRC_score))]=NA

for(i in c('PROM1','CD44','LGR5','ALDH1A1','EPCAM')){
 p1<- SpatialFeaturePlot(integrated_spatial,features = i,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
  ggsave(p1,filename = paste0(i,"_spatial.pdf"),width = 30,height = 35)

}

############################EMT GENESETS：

HALLMARK_EPITHELIAL_MESENCHYMAL<-clusterProfiler::read.gmt('/data/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gmt')

HALLMARK_EPITHELIAL_MESENCHYMAL1<-lapply(unique(HALLMARK_EPITHELIAL_MESENCHYMAL$term),function(x){HALLMARK_EPITHELIAL_MESENCHYMAL$gene[HALLMARK_EPITHELIAL_MESENCHYMAL$term==x]})


names(HALLMARK_EPITHELIAL_MESENCHYMAL1)<-unique(HALLMARK_EPITHELIAL_MESENCHYMAL$term)

Malignant_epithelial_cells<-AddModuleScore(Malignant_epithelial_cells,HALLMARK_EPITHELIAL_MESENCHYMAL1,name = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
#####################################################


integrated_spatial$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION<-mapvalues(colnames(integrated_spatial),colnames(Malignant_epithelial_cells),Malignant_epithelial_cells$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1)

integrated_spatial$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <-as.numeric(integrated_spatial$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)


SpatialFeaturePlot(integrated_spatial,features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")




integrated_spatial$CDH1 <-integrated_spatial@assays$SCT@scale.data['CDH1',]
integrated_spatial$CDH2 <-integrated_spatial@assays$SCT@scale.data['CDH2',]
integrated_spatial$VIM <-integrated_spatial@assays$SCT@scale.data['VIM',]
integrated_spatial$SNAI1 <-integrated_spatial@assays$SCT@scale.data['SNAI1',]
integrated_spatial$TWIST1 <-integrated_spatial@assays$SCT@scale.data['TWIST1',]
integrated_spatial$ZEB1 <-integrated_spatial@assays$SCT@scale.data['ZEB1',]
integrated_spatial$ZEB2 <-integrated_spatial@assays$SCT@scale.data['ZEB2',]



integrated_spatial$CDH1[which(is.na(integrated_spatial$CRC_score))]=NA
integrated_spatial$VIM[which(is.na(integrated_spatial$CRC_score))]=NA
integrated_spatial$ZEB1[which(is.na(integrated_spatial$CRC_score))]=NA
integrated_spatial$ZEB2[which(is.na(integrated_spatial$CRC_score))]=NA

setwd("/2024-2-20/04.EMT")

for(i in c('CDH1','VIM','ZEB1','ZEB2')){
  p1<- SpatialFeaturePlot(integrated_spatial,features = i,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
  ggsave(p1,filename = paste0(i,"_spatial.pdf"),width = 30,height = 35)

}





cancer_marker_score_df1 <- Malignant_epithelial_cells@meta.data%>%
  dplyr::select(orig.ident,group,HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1)%>%group_by(orig.ident)%>%
  dplyr::mutate(mean_var=mean(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1,na.rm = TRUE),sd_var=sd(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1,na.rm = TRUE))%>%
  ungroup()

View(cancer_marker_score_df1)
View(cancer_marker_score_df)
str(cancer_marker_score_df1)
cancer_marker_score_df1$group<-factor(cancer_marker_score_df1$group,levels=c("TT",'NTT'))
colnames(cancer_marker_score_df1)[1]<-'samples'
colnames(cancer_marker_score_df1)[4]<-'EMT_score'
# 开始绘制 点图
ggplot(data = cancer_marker_score_df1,aes(x=group,y=EMT_score,color=samples))+
  geom_point(stat="identity",position = position_dodge(width = 0.75))+
  geom_errorbar(aes(ymin=EMT_score -sd_var,ymax=EMT_score+sd_var),width=0.2,position = position_dodge(0.75)) +
  theme_classic2()+ylab(label = "EMT signature score") +
  geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)


tt_score<-subset(cancer_marker_score_df1,group=='TT')
ntt_score<-subset(cancer_marker_score_df1,group=='NTT')
tt_vs_ntt<-t.test(tt_score$EMT_score,ntt_score$EMT_score)



#################################差异基因热图
# macrophage
gene_set1<-c("COMP/RAP2B/ACTG1/CSRP1/DGKH/SERPINE1/PPIA/SERPING1/SYK/LMAN1/GNAQ/JAK2/F13A1/PRDX2/THBS1/FAP/FCER1G/DGKE/AP3B1/SERPINA1/JMJD1C/CD9/PIK3CA/LCK/PLCG2/LNPK/DGKZ/STXBP1/FERMT3/PTPN6/VKORC1/PLEK")
gene_set1<-c("JUNB/BTG1/CEBPB/DUSP1/CCND1/BHLHE40/TRIP10/RELB/CCL5/BTG2/SAT1/KLF2/CD44/PNRC1/SERPINE1/TNFAIP3/CFLAR/SPSB1/ZFP36/MSC/MAP3K8/LITAF/MCL1/RHOB/CD83/LDLR/VEGFA/FOSL2/PER1/GPR183/TANK/IL1B/CXCL1/NINJ1/ETS2/HES1/MARCKS/EGR1/REL/NFIL3/IER5/B4GALT5/ACKR3/PLPP3/ABCA1/FOS/KLF6/EFNA1/TNIP2/JAG1/IER2/PLEK")
gene_set1<-c("IGHG3，TMSB10，IGKC，IGHG1，CEACAM6，IGLC1，MGP，GGH，JCHAIN，SPP1")
# endo
gene_set1 <-c("IGHG3，IGKC，TMSB10，IGHG1，IGHA1，CEACAM6，MGP，PIM2，CFL1，HIST1H1C，CXCL8")
gene_set1<-c("PARK7/COX5A/COX5B/NDUFA13/UQCRQ/NDUFS5/AKR7A3/COX6A1/NDUFB8/ETFRF1/NDUFA4/GPD2/AKR1B1/COX7C/CYCS/COX6B1/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX11/COX4I1/PLEC")
gene_set1<-c("PARK7/COX5A/COX5B/UQCRQ/NDUFS5/COX6A1/NDUFB8/ETFRF1/NDUFA4/GPD2/COX7C/CYCS/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX4I1/PLEC")
gene_set1<-c("PARK7/COX5A/COX5B/UQCRQ/NDUFS5/COX6A1/NDUFB8/NDUFA4/COX7C/CYCS/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX4I1")
gene_set1<-strsplit(gene_set1,"/")[[1]]
gene_set1<-strsplit(gene_set1,"，")[[1]]

####################################################


carcinoma_regions_integrated1_Macrophage<-subset(carcinoma_regions_integrated1,idents=c("Macrophage"))
carcinoma_regions_integrated1_Endothelial_cells<-subset(carcinoma_regions_integrated1,idents=c("Endothelial_cells"))
# 巨噬细胞
Idents(carcinoma_regions_integrated1_Macrophage)<-carcinoma_regions_integrated1_Macrophage$orig.ident
expr1<-AverageExpression(carcinoma_regions_integrated1_Macrophage)
# 内皮细胞
Idents(carcinoma_regions_integrated1_Endothelial_cells)<-carcinoma_regions_integrated1_Endothelial_cells$orig.ident
expr1<-AverageExpression(carcinoma_regions_integrated1_Endothelial_cells)

expr <-expr1$SCT
expr<-expr[gene_set1,]
# pheatmap::pheatmap(expr,scale = 'column')
color1<-colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)
pheatmap::pheatmap(expr,scale = 'row',
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-2,2,length.out=101))



################################ 热图 2024-3-4，重新调整热图为 TT和NTT
# macrophage
gene_set1<-c("COMP/RAP2B/ACTG1/CSRP1/DGKH/SERPINE1/PPIA/SERPING1/SYK/LMAN1/GNAQ/JAK2/F13A1/PRDX2/THBS1/FAP/FCER1G/DGKE/AP3B1/SERPINA1/JMJD1C/CD9/PIK3CA/LCK/PLCG2/LNPK/DGKZ/STXBP1/FERMT3/PTPN6/VKORC1/PLEK")
gene_set1<-c("JUNB/BTG1/CEBPB/DUSP1/CCND1/BHLHE40/TRIP10/RELB/CCL5/BTG2/SAT1/KLF2/CD44/PNRC1/SERPINE1/TNFAIP3/CFLAR/SPSB1/ZFP36/MSC/MAP3K8/LITAF/MCL1/RHOB/CD83/LDLR/VEGFA/FOSL2/PER1/GPR183/TANK/IL1B/CXCL1/NINJ1/ETS2/HES1/MARCKS/EGR1/REL/NFIL3/IER5/B4GALT5/ACKR3/PLPP3/ABCA1/FOS/KLF6/EFNA1/TNIP2/JAG1/IER2/PLEK")
gene_set1<-c("IGHG3，TMSB10，IGKC，IGHG1，CEACAM6，IGLC1，MGP，GGH，JCHAIN，SPP1")
# endo
gene_set1 <-c("IGHG3，IGKC，TMSB10，IGHG1，IGHA1，CEACAM6，MGP，PIM2，CFL1，HIST1H1C，CXCL8") # 前10个基因
gene_set1<-c("PARK7/COX5A/COX5B/NDUFA13/UQCRQ/NDUFS5/AKR7A3/COX6A1/NDUFB8/ETFRF1/NDUFA4/GPD2/AKR1B1/COX7C/CYCS/COX6B1/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX11/COX4I1/PLEC") # 电子呼吸链传递异常基因
gene_set1<-c("PARK7/COX5A/COX5B/UQCRQ/NDUFS5/COX6A1/NDUFB8/ETFRF1/NDUFA4/GPD2/COX7C/CYCS/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX4I1/PLEC")  #   电子传递链基因
gene_set1<-c("PARK7/COX5A/COX5B/UQCRQ/NDUFS5/COX6A1/NDUFB8/NDUFA4/COX7C/CYCS/UQCR10/NDUFB5/UQCRH/UQCR11/NDUFS6/COX4I1") # 消除氧自由基的相关蛋白基因
gene_set1<-strsplit(gene_set1,"/")[[1]]
gene_set1<-strsplit(gene_set1,"，")[[1]]

####################################################


carcinoma_regions_integrated1_Macrophage<-subset(carcinoma_regions_integrated1,idents=c("Macrophage"))
carcinoma_regions_integrated1_Endothelial_cells<-subset(carcinoma_regions_integrated1,idents=c("Endothelial_cells"))
# 巨噬细胞
# Idents(carcinoma_regions_integrated1_Macrophage)<-carcinoma_regions_integrated1_Macrophage$orig.ident
Idents(carcinoma_regions_integrated1_Macrophage)<-carcinoma_regions_integrated1_Macrophage$group
expr1<-AverageExpression(carcinoma_regions_integrated1_Macrophage)
# 内皮细胞
# Idents(carcinoma_regions_integrated1_Endothelial_cells)<-carcinoma_regions_integrated1_Endothelial_cells$orig.ident
Idents(carcinoma_regions_integrated1_Endothelial_cells)<-carcinoma_regions_integrated1_Endothelial_cells$group
expr1<-AverageExpression(carcinoma_regions_integrated1_Endothelial_cells)

expr <-expr1$SCT
expr<-expr[gene_set1,]
# pheatmap::pheatmap(expr,scale = 'column')
color1<-colorRampPalette(rev(brewer.pal(n = 7, name =
                                         "RdYlBu")))(100)

expr_scale<-t(scale(t(expr),center = T,scale = F))


pheatmap::pheatmap(expr_scale,scale = 'none',
                   cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out=101))




pheatmap::pheatmap(expr,scale = 'row',
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-6,0,length.out=101))

pheatmap::pheatmap(expr,scale = 'none',
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-6,6,length.out=101))





# 提取TT
tt_carcinoma_regions_integrated1 <- subset(carcinoma_regions_integrated1,group=='TT')
# 提取 NTT
ntt_carcinoma_regions_integrated1 <- subset(carcinoma_regions_integrated1,group=='NTT')


tt_carcinoma_regions_integrated1_aver<-AverageExpression(tt_carcinoma_regions_integrated1)
ntt_carcinoma_regions_integrated1_aver<-AverageExpression(ntt_carcinoma_regions_integrated1)




library(psych)
library(pheatmap)
library(corrplot)
library(export)
#
# st_p5 <- tt_carcinoma_regions_integrated1_aver$SCT
# st_p5 <- data.frame(ntt_carcinoma_regions_integrated1_aver$SCT,check.names = F)
st_p5 <- data.frame(tt_carcinoma_regions_integrated1_aver$integrated,check.names = F)
# # order_cluster<-colnames(aa$corr)
order_cluster<-c("CD4._T_cells","Mast_cells","Macrophage","CD8._T_cells","B_cells",
                 "Nonmalignant_epithelial_cells","Malignant_epithelial_cells","Fibroblast_cells",
                 "SMC","Endothelial_cells")
st_p5<-st_p5[,order_cluster]


st_p5 <- data.frame(ntt_carcinoma_regions_integrated1_aver$integrated,check.names = F)
# order_cluster<-colnames(aa$corr)
st_p5<-st_p5[,order_cluster]


# tt <- corr.test(t(st_p5),method="pearson",adjust="holm")
tt <- corr.test(st_p5,method="pearson",adjust="holm")

# tt$r[is.na(tt$r)] <- 0
# tt$p[is.na(tt$p)] <- 1
# TT和NTT的这两张图更改下颜色-1：#00008B，-0.5：#7B68EE，0：#FFFFFF，0.5：#FFD700，1：#FF0000
# color_palette<-c("#00008B","#7B68EE","#FFFFFF","#FFD700","#FF0000")
color_palette<-c("#1E90FF","#7B68EE","#FFFFFF","#FFD700","#FF0000")

# 显示相关性数字
corrplot(as.matrix(tt$r),col = colorRampPalette(color_palette)(100),
         method = "color",tl.col="black",
         tl.cex = 0.7,cl.pos = "r",cl.lim=c(-1,1),
         cl.ratio = 0.1,cl.cex=0.8,cl.length=5,p.mat = as.matrix(tt$p),
         sig.level = 1,outline="white",insig = "blank",type="full",addCoef.col = "black")


# 不显示相关性数字
corrplot(as.matrix(tt$r),col = colorRampPalette(color_palette)(100),
         method = "color",tl.col="black",
         tl.cex = 0.7,cl.pos = "r",cl.lim=c(-1,1),
         cl.ratio = 0.1,cl.cex=0.8,cl.length=5,p.mat = as.matrix(tt$p),
         sig.level = 0.05,outline="white",insig = "blank",type="full")





##################################热图调整绘制 2024-3-5

# 4.2）热图展示巨噬细胞中以下基因在TT和NTT表达差异（分成两个热图分别展示）
# PDL1/CD274，HAVCR2，VEGFA，IL1B，IL10RA，SOCS3，STAT3，CEBPB，CD163
#
# SERPINE1，SERPINE2，C1qa，LGMN，C3，THBS1，TIMP1，CSRP1，SPARC，CTSB，A2M

genes <- c("CD274", "HAVCR2", "VEGFA", "IL1B", "IL10RA", "SOCS3", "STAT3", "CEBPB", "CD163",
           "SERPINE1", "SERPINE2", "C1QA", "LGMN", "C3", "THBS1", "TIMP1", "CSRP1", "SPARC", "CTSB", "A2M")



carcinoma_regions_integrated1_Macrophage<-subset(carcinoma_regions_integrated1,idents=c("Macrophage"))
carcinoma_regions_integrated1_Endothelial_cells<-subset(carcinoma_regions_integrated1,idents=c("Endothelial_cells"))
# 巨噬细胞
# Idents(carcinoma_regions_integrated1_Macrophage)<-carcinoma_regions_integrated1_Macrophage$orig.ident
Idents(carcinoma_regions_integrated1_Macrophage)<-carcinoma_regions_integrated1_Macrophage$group
expr1<-AverageExpression(carcinoma_regions_integrated1_Macrophage)
# 内皮细胞
# Idents(carcinoma_regions_integrated1_Endothelial_cells)<-carcinoma_regions_integrated1_Endothelial_cells$orig.ident
Idents(carcinoma_regions_integrated1_Endothelial_cells)<-carcinoma_regions_integrated1_Endothelial_cells$group
expr1<-AverageExpression(carcinoma_regions_integrated1_Endothelial_cells)

expr <-expr1$SCT
expr<-expr[genes,]
# pheatmap::pheatmap(expr,scale = 'column')
color1<-colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)



expr_scale<-t(scale(t(expr),center = T,scale = F))


pheatmap::pheatmap(expr_scale,scale = 'none',
                   cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out=101))




pheatmap::pheatmap(expr,scale = 'row',
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-6,0,length.out=101))

pheatmap::pheatmap(expr,scale = 'none',
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-6,6,length.out=101))







genes <- c("HLA-E", "HLA-DMA", "HLA-DRA", "HLA-DPB1", "HLA-F", "HLA-DPA1", "HLA-A", "IGHG2", "HMGB1", "IL6ST")


expr <-expr1$SCT
expr<-expr[genes,]

color1<-colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)


expr_scale<-t(scale(t(expr),center = T,scale = F))


pheatmap::pheatmap(expr_scale,scale = 'none',
                   cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out=101))





genes <- c("CXCL8", "CXCL1", "CXCL10", "CXCL14", "CXCL6", "CXCL3", "IL11", "IL4I1", "ICAM1", "VCAM1")



expr <-expr1$SCT
expr<-expr[genes,]

color1<-colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)


expr_scale<-t(scale(t(expr),center = T,scale = F))

pheatmap::pheatmap(expr_scale,scale = 'none',
                   cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out=101))


#######################################################

cell_type<-levels(carcinoma_regions_integrated1)

df1<-integrated_spatial@meta.data[,cell_type]

SpatialFeaturePlot(integrated_spatial,features = 'B_cells',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey",oob=scales::squish)
setwd('/share/home/hhhong/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-5/02.heatmap_adjust/4.5_plot')
for(cls in cell_type){
  # p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,limits=c(0.0,1.0),na.value = "lightgrey",oob=scales::squish)
  p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey",oob=scales::squish)

  ggsave(p1,filename = paste0(cls,"_featureplot.pdf"),width = 18,height = 18)
}



####################################################### 2024-3-6




# 调整空间相关性热图

######################################################2024-3-6

aa<-fread('/01.correlation_adjust/FigS4_ST_P5.transfer.xls',header = T,data.table = F)


DimPlot(carcinoma_regions_integrated1)

tt_sct_df <- tt_carcinoma_regions_integrated1@meta.data%>%dplyr::select(c(16:25))  %>%data.frame(check.names = F,check.rows = F)






# 提取TT
tt_carcinoma_regions_integrated1 <- subset(carcinoma_regions_integrated1,group=='TT')
tt_sct_df <- tt_carcinoma_regions_integrated1@meta.data%>%dplyr::select(c(16:25))  %>%data.frame(check.names = F,check.rows = F)

# 提取 NTT
ntt_carcinoma_regions_integrated1 <- subset(carcinoma_regions_integrated1,group=='NTT')
ntt_sct_df <- ntt_carcinoma_regions_integrated1@meta.data%>%dplyr::select(c(16:25))  %>%data.frame(check.names = F,check.rows = F)


##########################################总的相关性热图，分TT和NTT
##########################################总的相关性热图，分TT和NTT 2024-5-8 重新修改图的颜色

tt_ntt_sct_df <- ntt_carcinoma_regions_integrated1@meta.data%>%dplyr::select(c(16:25))  %>%data.frame(check.names = F,check.rows = F)

library(psych)
library(pheatmap)
library(corrplot)
library(export)
#

# st_p5 <- data.frame(tt_sct_df,check.names = F)
st_p5 <- data.frame(tt_ntt_sct_df,check.names = F)


order_cluster1<-c("Nonmalignant_epithelial_cells","CD4._T_cells","Fibroblast_cells",
                  "SMC","Endothelial_cells","B_cells","Mast_cells","Macrophage","CD8._T_cells",
                  "Malignant_epithelial_cells")
st_p5<-st_p5[,order_cluster1]

tt <- corr.test(st_p5,method="pearson",adjust="holm")
color_palette<-c("#1E90FF","#7B68EE","#FFFFFF","#FFD700","#FF0000")

dev.off()

corrplot(as.matrix(tt$r),col = colorRampPalette(color_palette)(100),
         method = "color",tl.col="black",
         tl.cex = 0.7,cl.pos = "r",cl.lim=c(-1,1),
         mar = c(0, 0, 0, 0),
         cl.ratio = 0.1,cl.cex=0.8,cl.length=5,p.mat = as.matrix(tt$p),
         sig.level = 1,outline="white",insig = "blank",type="full",addCoef.col = "black")


# 
corrplot(as.matrix(tt$r),col = colorRampPalette(color_palette)(100),
         method = "color",tl.col="black",
         tl.cex = 0.7,cl.pos = "r",cl.lim=c(-1,1),
         cl.ratio = 0.1,cl.cex=0.8,cl.length=5,p.mat = as.matrix(tt$p),
         sig.level = 0.05,outline="white",insig = "blank",type="full")



aa<-corrplot(as.matrix(tt$r), order = "hclust",col = colorRampPalette(c("darkblue","white","darkred"))(100),
             method = "color",tl.col="black",
             tl.cex = 0.7,cl.pos = "r",cl.lim=c(-1,1),
             cl.ratio = 0.1,cl.cex=0.8,cl.length=5,p.mat = as.matrix(tt$p),
             sig.level = 0.05,outline="white",insig = "blank",type="full")


order_cluster<-colnames(aa$corr)



min(as.matrix(tt$r))
max(as.matrix(tt$r))
corrplot(as.matrix(tt$r), order = "hclust",col = colorRampPalette(c("darkblue","darkred"))(100),
         method = "color",tl.col="black",
         tl.cex = 0.7,cl.pos = "r",cl.lim=c(-1,1),
         cl.ratio = 0.1,cl.cex=0.8,cl.length=5,p.mat = as.matrix(tt$p),
         sig.level = 0.05,outline="white",insig = "blank",type="full")



#################################################################### 绘制特定分群的相关性



compare_cluster<-list(a=c("Fibroblast_cells","Malignant_epithelial_cells"),
                      b=c("CD4._T_cells","Malignant_epithelial_cells"),
                      d=c("CD4._T_cells","Nonmalignant_epithelial_cells"),
                      e=c("CD4._T_cells","Fibroblast_cells"))

for(i in name(compare_cluster)){
  dat<-tt_sct_df[,compare_cluster[[i]]]
}
aa<-tt_ntt_sct_df[,c("Fibroblast_cells","Malignant_epithelial_cells")]


aa_cor<-cor(t(aa),method = 'spearman')


aa<-tt_ntt_sct_df[,c("Fibroblast_cells","Malignant_epithelial_cells")]

aa1<-aa[apply(aa,1,sum)>0,]
aa2<-subset(aa1,!(Fibroblast_cells==0 | Malignant_epithelial_cells ==0 ))



g <- ggplot(aa2, aes(Malignant_epithelial_cells,Fibroblast_cells))


df.lm<-lm(`Malignant_epithelial_cells`~
            `Fibroblast_cells`,
          data=aa2)
p_value <- summary(df.lm)$coefficients["Fibroblast_cells", "Pr(>|t|)"]
# 提取R-squared
r_squared <- summary(df.lm)$r.squared
p_value <- summary(df.lm)$coefficients
# 获取R值
R_value <- sqrt(summary(df.lm)$r.squared)

# 可以
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

summary_output <- capture.output(summary(df.lm))

p_value_line <- grep("p-value:", summary_output, value = TRUE)
p_value <- sub(".*p-value: (.*)$", "\\1", p_value_line)
matches <-stringr::str_match(p_value, "< (-?\\d+\\.?\\d*)e(-?\\d+)")

# Scatterplot
g + geom_point(colour='red') +
  geom_smooth(method="lm", se=F,color='black') +
  labs(subtitle="correlation",
       y="Malignant_epithelial_cells",
       x="Fibroblast_cells",
       title="",
       caption="")+
  theme_bw()+
  annotate(geom = "text",

           x=median(aa2$Malignant_epithelial_cells), y=max(aa2$Fibroblast_cells),
           label = paste0("R=",round(r_squared,2),",","P ",matches[1,1]))





############################################################# 将单细胞注释结果整合到seurat object 2024-3-8



crc_dat_singlecell <- readRDS('/data/adjust_annotation_singlecell.rds')
colnames(crc_dat_singlecell@meta.data)
DimPlot(crc_dat_singlecell,label = T,group.by = 'subtype',repel = T)
DimPlot(crc_dat_singlecell,label = T,group.by = 'subtype',repel = T)
Idents(crc_dat_singlecell)<-crc_dat_singlecell$subtype
DimPlot(crc_dat_singlecell,cells.highlight = WhichCells(crc_dat_singlecell,idents = "Regulatory T cells"))

FeaturePlot(crc_dat_singlecell,features = 'FOXP3',order = T)
FeaturePlot(crc_dat_singlecell,features = 'CTLA4',order = T)



crc_dat_singlecell$subtype2 <- crc_dat_singlecell$subtype

#############################################
DimPlot(crc_dat_singlecell,label = T,group.by = 'subtype2',repel = F)
DimPlot(crc_dat_singlecell,label = T,group.by = 'seurat_clusters',repel = F)
# 巨噬细胞
crc_dat_singlecell$subtype2[crc_dat_singlecell$subtype2=="SPP1+"]="M2_macrophage"
#
CAF <- c('ACTA2', 'FAP', 'PDGFRA', 'PDGFRB', 'PDPN', 'S100A4', 'TNC', 'VIM', 'COL1A1','COL1A3','COL1A2','DCN')
FeaturePlot(crc_dat_singlecell,features = CAF,order = T)
DimPlot(crc_dat_singlecell,label = T,group.by = 'subtype1',repel = F)
unique(crc_dat_singlecell$subtype1)
crc_dat_singlecell$subtype2[crc_dat_singlecell$subtype1=="Fibroblast_cells"]="CAF"


T_Ex <- c('PDCD1', 'CTLA4', 'LAG3', 'HAVCR2','CD38')

FeaturePlot(crc_dat_singlecell,features = T_Ex,order = T)

crc_dat_singlecell$subtype2[crc_dat_singlecell$subtype1=="Fibroblast_cells"]="CAF"
Idents(crc_dat_singlecell)<-crc_dat_singlecell$subtype2
crc_dat_singlecell$subtype2[rownames(crc_dat_singlecell@meta.data)%in%WhichCells(crc_dat_singlecell,expression = PDCD1 >0 & Idents(crc_dat_singlecell)=="CD8+ T cells")]="Exhausted_CD8+T_cells"



# STAT4,IFNG,IL12,IL2,RB2,TBX21   Th1
# GATA3,STAT6,IL4,IL10,TGFB1 Th2
th1<-strsplit("STAT4,IFNG,IL12A,IL2,RBL2,TBX21",",")[[1]]
th2<-strsplit("GATA3,STAT6,IL4,IL10,TGFB1",",")[[1]]


FeaturePlot(crc_dat_singlecell,features = th1,order = T)

FeaturePlot(crc_dat_singlecell,features = th2,order = T)
FeaturePlot(crc_dat_singlecell,features = "CD3E",order = T)


Th1 <- c('CCR1','CCR5','CD3','CD4','CXCR3','TBX21','GATA1','GATA4','MS4AB','CCL5','IFNG')
FeaturePlot(crc_dat_singlecell,features = Th1,raster = F)
crc_dat_singlecell$subtype2[crc_dat_singlecell$seurat_clusters=="11"]="Th1"


Th2 <- c('IL1RL1',"GATA3", 'IL13', 'IL5', 'IL4','IL17RB', 'LTB4R1','CCR8')
FeaturePlot(crc_dat_singlecell,features = Th2,raster = F)

M_MDSCs <- c('CXCR1','WFDC17','CD14','CD274')
FeaturePlot(crc_dat_singlecell,features = M_MDSCs,raster = F)

PMN_MDSCs <- c('ITGAM','CXCL1','CXCL2','CXCR1','CXCR2','LILRA3','TREM1','OLR1','CD84','ARG1','ARG2')
FeaturePlot(crc_dat_singlecell,features = PMN_MDSCs,raster = F,order = T)

FeaturePlot(crc_dat_singlecell,features = "CD14",raster = F,order = T)
Idents(crc_dat_singlecell)<-crc_dat_singlecell$seurat_clusters

WhichCells(crc_dat_singlecell,expression = ITGAM >0 & CD14 == 0 & Idents(crc_dat_singlecell)=="2")
PMN_MDSCs=WhichCells(crc_dat_singlecell,expression = ITGAM >0 & CD14 == 0 & Idents(crc_dat_singlecell)=="2")
crc_dat_singlecell$subtype2[rownames(crc_dat_singlecell@meta.data)%in%PMN_MDSCs]="PMN-MDSCs"

table(Idents(crc_dat_singlecell)=="2")
# FUT4 CD15
M_MDSCs=WhichCells(crc_dat_singlecell,expression = ITGAM >0 & CD14 > 0 & CD33 >0&  Idents(crc_dat_singlecell)=="2")
FeaturePlot(crc_dat_singlecell,features = "CD33",raster = F)
crc_dat_singlecell$subtype2[rownames(crc_dat_singlecell@meta.data)%in%M_MDSCs]="M-MDSCs"
DimPlot(crc_dat_singlecell,label = T,group.by = 'subtype2',repel = F)


tDC <- c('IDO1', 'CCR7', 'LGALS3', 'LGALS9', 'NECTIN2','CCL19','CCL17','BIRC3','CXCL13','CRIP1')
FeaturePlot(crc_dat_singlecell,features = tDC,raster = F)

Idents(crc_dat_singlecell)<-crc_dat_singlecell$subtype2
levels(crc_dat_singlecell)
TDC=WhichCells(crc_dat_singlecell,expression = CCR7 >0 & IDO1 > 0 & CRIP1 >0&  Idents(crc_dat_singlecell)=="cDC")
crc_dat_singlecell$subtype2[rownames(crc_dat_singlecell@meta.data)%in%TDC]="tDCs"

DimPlot(crc_dat_singlecell,label = T,group.by = 'subtype2')
crc_dat_singlecell$subtype2[crc_dat_singlecell$subtype1=="SMC"]="SMC"
crc_dat_singlecell$subtype2[crc_dat_singlecell$subtype1=="Epithelial_cells"]="Epithelial_cells"
crc_dat_singlecell$subtype2[crc_dat_singlecell$subtype1=="B_cells"]="B_cells"
crc_dat_singlecell$subtype2[crc_dat_singlecell$subtype1=="Endothelial_cells"]="Endothelial_cells"
crc_dat_singlecell1<-subset(crc_dat_singlecell,idents=c("CD19+CD20+ B","Regulatory T cells","CMS1","CMS2","CMS3","Unknown","T helper 17 cells"),invert=T)

DimPlot(crc_dat_singlecell1,label = T,group.by = 'subtype2')
DimPlot(crc_dat_singlecell1,label = T,group.by = 'subtype1')

crc_dat_singlecell1$subtype2[crc_dat_singlecell1$subtype2=="Pro-inflammatory"]="M1_macrophage"


CAF <- c('ACTA2', 'FAP', 'PDGFRA', 'PDGFRB', 'PDPN', 'S100A4', 'TNC', 'VIM', 'COL1A1','COL1A3','COL1A2','DCN')
FeaturePlot(crc_dat_singlecell1,features = CAF,order = T)


crc_dat_singlecell1_subset<-subset(crc_dat_singlecell1,idents=c("Th1","Exhausted_CD8+T_cells","M-MDSCs","M2_macrophage","tDCs","cDC"))



############################################################2024-3-10

integrated_spatial_carcinoma_regions<-readRDS('/data/integrated_spatial_carcinoma_regions.rds')



for(i in 46:62){
  nm = colnames(integrated_spatial@meta.data)[i]

  p1<-SpatialFeaturePlot(integrated_spatial,features = nm,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")

  ggsave(p1,filename = paste0(nm,"_featureplot.pdf"),width = 18,height = 18)
}


SpatialFeaturePlot(integrated_spatial,features = "Epithelial_cells",crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
tt3<-subset(integrated_spatial,orig.ident=='TT.3')

SpatialFeaturePlot(tt3,features = "Epithelial_cells",crop = F,alpha = 1,interactive=F,ncol = 4)


VlnPlot(tt3,features = 'Epithelial_cells',pt.size = 0,group.by = 'orig.ident')


tt3_2<-subset(tt3,cells=rownames(tt3@meta.data)[!is.na(tt3@meta.data$Epithelial_cells)])


tt3_1<-all_portion$TT.3%>%data.frame()

boxplot(tt3_1$Epithelial_cells)

###############################################################################2024-3-11

cols=c("#43387F","#288D88","#8AD24C","#FBE625")
setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-11/02.color_adjust")


cell_type = colnames(integrated_spatial@meta.data)[46:62]
for(cls in cell_type){
  # p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,limits=c(0.0,1.0),na.value = "lightgrey",oob=scales::squish)
  p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,na.value = "lightgrey",oob=scales::squish)

  ggsave(p1,filename = paste0(cls,"_featureplot.pdf"),width = 18,height = 18)
}

############细胞数统计表


NTT4_proportions$dominant_celltype <- apply(NTT4_proportions, 1, function(x) names(NTT4_proportions)[which.max(x)])


df<-integrated_spatial@meta.data[46:62]

dominant_celltype<-apply(df,1,function(x)names(df)[which.max(x)])
dominant_celltype <- as.character(dominant_celltype)


df$dominant_celltype<-dominant_celltype


df$dominant_celltype[df$dominant_celltype=="character(0)"]=NA
# df$dominant_celltype[df$dominant_celltype=="NA"]<-NA
df1<-na_if(df, 'NA')

#############################################
integrated_spatial$dominant_celltype <- df$dominant_celltype


df2<-df1%>%dplyr::select(dominant_celltype,sample)%>%dplyr::group_by(sample)%>%table()
df3<-df2%>%data.frame()%>%tidyr::spread(key=dominant_celltype,value=Freq)

fwrite(df3,file = "~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-11/03.cell_stat/stat.xls",sep = '\t',col.names = T)


# 4.4）热图展示基因（TT和NTT整个病理癌区域所有细胞类型都放进去比较）
geneset1<-"ISG15，LGALS9，FGL2，CXCL10，CCL2，CCL12，ARG1，CD160，CD244，CXCL12，SERPINB9，CD274，CTLA4，TGFB1，TGFB2，IDO1，LAG3，HAVCR2， VSIR，PDCD1，CD276，VTCN1"
geneset1<-"HLA-A，HLA-B，HLA-C，HLA-DPB1，HLA-DQB1，HLA-DRA，HLA-DMA，HLA-DOA，HLA-E，HLA-G，CD1A，CD1B，CD1C，CD1D，CD1E"
geneset1<-"CTLA4，PDCD1，CD274，PDCD1LG2，LAG3，TIGIT，HAVCR2，IDO1，IL4I1，ENTPD1，NT5E，KLRC1，SIRPA，ICOS，TLR3，TLR4，TLR7，TLR8，TLR9，VSIR，SIGLEC15，PVRIG，KIR2DL1，KIR3DL1，KIR2DL3，SIGLEC7，SIGLEC9，HLA-G，CD47，STAB1，CD24，SIGLEC10"
# 细胞因子基因
geneset1<-"TGFB1，TGFB2，TGFB3，IL10，IL6，EBI3，IL12A，SOCS1，SOCS2，VEGFA，VEGFB，VEGFC，VEGFD，IL4，IL13，IL2，IFNA1，IFNA2，IFNG，IL12A，IL12B，TNF，IL7，IL15，IL21，IL18"
# 趋化因子基因
geneset1<-"CXCL8，CXCL12，CCL2，CCL5，CCL20，CXCL9，CXCL10，CXCL1，CCL18，CX3CL1，CCL17，CCL21，CCL22，CXCL11"
geneset<-strsplit(geneset1,"，")[[1]]




Idents(carcinoma_regions_integrated1)<-carcinoma_regions_integrated1$orig.ident
expr<-AverageExpression(carcinoma_regions_integrated1)
geneset[!geneset%in% rownames(expr$SCT)]

expr1<-expr$SCT[geneset[geneset%in% rownames(expr$SCT)],]

pheatmap::pheatmap(expr1,scale = 'row',
                   cluster_rows = F,cluster_cols = F,breaks = seq(-2,2,length.out=101))


##########################################2024-3-11
setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-11/05.TT_vs_NTT")
diff_tt_vs_ntt<-subset(FindMarkers(carcinoma_regions_integrated1,group.by = 'group',ident.1 = 'TT',ident.2 = 'NTT'),p_val<0.05)

diff_tt_vs_ntt$genes<-rownames(diff_tt_vs_ntt)

diff_tt_vs_ntt$regulated<-ifelse(diff_tt_vs_ntt$avg_log2FC>0,"up-regulated","down-regulated")

setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-11/05.TT_vs_NTT")
fwrite(diff_tt_vs_ntt,file = "TT_vs_NTT_p0.05.xls",sep="\t",col.names = T)


###################################################################2024-3-12
###################################################################2024-3-12
# cDC：将全部8张图 图例上限调整为0.015
# Exhausted_CD8+T_cells：将全部8张图 图例上限调整为0.15
# M1_macrophage：将全部8张图 图例上限调整为0.01

# M-MDSCs：TT1图例上限调整为0.15，NTT3图例上限调整为0.02
# PMN-MDSCs：TT1图例上限调整为0.05，NTT1图例上限调整为0.03


# tDCs:将全部8张图 图例上限调整为0.02出一张图，0.01出一张图，0.008出一张图

#############################################################################
#############################################################################
#############################################################################
#############################################################################

cols=c("#43387F","#288D88","#8AD24C","#FBE625")
setwd("/2024-3-11/02.color_adjust")


cell_type = colnames(integrated_spatial@meta.data)[46:62]
for(cls in cell_type){
  # p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,limits=c(0.0,1.0),na.value = "lightgrey",oob=scales::squish)
  p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,na.value = "lightgrey",oob=scales::squish)

  ggsave(p1,filename = paste0(cls,"_featureplot.pdf"),width = 18,height = 18)
}

setwd("/2024-3-12/01.legend_adjust")
cls="cDC"
cls="Exhausted_CD8+T_cells"
cls="M1_macrophage"
cls="tDCs"



p1<-SpatialFeaturePlot(integrated_spatial,features = cls,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = cols,limits=c(0,0.008),na.value = "lightgrey",oob=scales::squish)

ggsave(p1,filename = paste0(cls,"_0.008_featureplot.pdf"),width = 18,height = 18)




##########################################################################TIDE


carcinoma_regions_integrated1<-SCTransform(carcinoma_regions_integrated1,assay = "Spatial"
            , method = "poisson", verbose = FALSE,variable.features.n=18643)


Idents(carcinoma_regions_integrated1)<-carcinoma_regions_integrated1$orig.ident
expr1<-AverageExpression(carcinoma_regions_integrated1,slot = 'scale.data')

expr1<-expr1$SCT


fwrite(data.frame(carcinoma_regions_integrated1@assays$SCT@scale.data,check.names=F,check.rows = F),file = "/02.tide_analysis/total_splot_expr_scale.xls",sep = '\t',col.names = T,row.names = T)



################################################ adjust legend
setwd('2024-3-13/4.1_legend_adjust')


# 热图展示基因（下面都调整了每个类别里面的基因，删除了效果不佳的基因）
# 抗原呈递分子


# 4.4）热图展示基因（TT和NTT整个病理癌区域所有细胞类型都放进去比较）
# 抗原呈递分子
geneset1<-"HLA-A，HLA-C，HLA-DPB1，HLA-DQB1，HLA-DRA，HLA-DMA，HLA-DOA，HLA-E，HLA-G，CD1A，CD1B，CD1C，CD1D，CD1E"
# 免疫检查点基因
geneset1<-"CTLA4，PDCD1，CD274，PDCD1LG2，LAG3，TIGIT，HAVCR2，IDO1，IL4I1，ENTPD1，NT5E，KLRC1，SIRPA，ICOS，TLR3，TLR4，TLR7，TLR8，TLR9，VSIR，SIGLEC15，SIGLEC7，SIGLEC9，CD47，STAB1，CD24，SIGLEC10"
# 免疫抑制基因
geneset1<-"ISG15，LGALS9，FGL2，ARG1，CD160，CD244，SERPINB9，VTCN1"
# 趋化因子
geneset1<-"CXCL8，CXCL12，CCL2，CCL5，CCL20，CXCL9，CXCL10，CXCL1，CCL17，CCL22"
# 细胞因子基因
geneset1<-"TGFB2，TGFB3，IL10，IL6，EBI3，IL12A，SOCS1，SOCS2，VEGFA，VEGFB，VEGFC，VEGFD，IL4，IL13"
geneset<-strsplit(geneset1,"，")[[1]]
geneset_temp<-strsplit(geneset1,"，")[[1]]

###############整合的
total_geneset1<-"CTLA4，PDCD1，CD274，PDCD1LG2，LAG3，TIGIT，HAVCR2，IDO1，IL4I1，ENTPD1，NT5E，KLRC1，SIRPA，ICOS，TLR3，TLR4，TLR7，TLR8，TLR9，VSIR，SIGLEC15，SIGLEC7，SIGLEC9，CD47，STAB1，CD24，SIGLEC10，ISG15，LGALS9，FGL2，ARG1，CD160，CD244，SERPINB9，VTCN1，TGFB2，TGFB3，IL10，IL6，EBI3，IL12A，SOCS1，SOCS2，VEGFA，VEGFB，VEGFC，VEGFD，IL4，IL13，CXCL8，CXCL12，CCL2，CCL5，CCL20，CXCL9，CXCL10，CXCL1，CCL17，CCL22，HLA-A，HLA-C，HLA-DPB1，HLA-DQB1，HLA-DRA，HLA-DMA，HLA-DOA，HLA-E，HLA-G，CD1A，CD1B，CD1C，CD1D，CD1E"
geneset<-strsplit(total_geneset1,"，")[[1]]

expr<-AverageExpression(carcinoma_regions_integrated1)

expr1<-expr$SCT[geneset[geneset%in% rownames(expr$SCT)],]
row_anno<-data.frame(genes=geneset,type=c(rep("Immune Checkpoints",27),
                                rep("Immunosuppressive Genes",8),
                                rep("Cytokines",14),rep("Chemokines",10),
                                rep("Antigen Presentation Molecules",14)),
                     colors=c(rep("#C77CFF",27),
                              rep("#00BFC4",8),
                              rep("#9ACD32",14),rep("#FFD92F",10),
                              rep("#b15928",14))
                     ,row.names = 1)


row_anno<-data.frame(genes=geneset,type=c(rep("Immune Checkpoints",27),
                                          rep("Immunosuppressive Genes",8),
                                          rep("Cytokines",14),rep("Chemokines",10),
                                          rep("Antigen Presentation Molecules",14))

                     ,row.names = 1)

ann_colors = list(
  type = c(`Immune Checkpoints`="#C77CFF",`Immunosuppressive Genes`="#00BFC4",
           `Cytokines`="#9ACD32",`Chemokines`="#FFD92F",`Antigen Presentation Molecules`="#b15928")

)

pheatmap::pheatmap(expr1,scale = 'row',annotation_row = row_anno,
                   annotation_colors = ann_colors,
                   cluster_rows = F,cluster_cols = F,
                   breaks = seq(-2,2,length.out=101),
                   border_color = NA,
                   gaps_row = cumsum(c(27,8,14,10,14)))



###################################################################2024-3-13
"抗原呈递分子"<-"HLA-A，HLA-C，HLA-DPB1，HLA-G，CD1A，CD1C，CD1D"
ff=list("抗原呈递分子"=strsplit("HLA-A，HLA-C，HLA-DPB1，HLA-G，CD1A，CD1C，CD1D","，")[[1]],
        "免疫检查点基因"=strsplit("CTLA4，PDCD1，CD274，PDCD1LG2","，")[[1]],
        "免疫抑制基因"=strsplit("LGALS9，FGL2，ARG1","，")[[1]],
        "趋化因子"=strsplit("CXCL12，CCL2，CCL5，CXCL9","，")[[1]],
        "细胞因子基因"=strsplit("TGFB2，TGFB3，IL10，VEGFA，VEGFB","，")[[1]])
ff
#####################################################
#####################################################绘制基因
#####################################################绘制基因

View(integrated_spatial@meta.data)
setwd("~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-13/4.3_gene_spatial")
for(nm in names(ff)){
  dir.create(nm)
  setwd(nm)
  for(j in ff[[nm]]){
    tmp_expr=as.numeric(integrated_spatial@assays$SCT@data[j,])
    tmp_expr[is.na(integrated_spatial$Th1)]<-NA
    new_feature<-paste0("gene_",j)

    integrated_spatial@meta.data[[new_feature]]<-tmp_expr

    p1<-SpatialFeaturePlot(integrated_spatial,features = new_feature,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colours = FeaturePalettes[['Spatial']],na.value = "lightgrey",oob=scales::squish)

    ggsave(p1,filename = paste0(j,"_featureplot.pdf"),width = 18,height = 18)
    integrated_spatial@meta.data[[new_feature]]<-NULL
  }
  setwd("../")

}

View(integrated_spatial@meta.data)
integrated_spatial@assays$SCT@scale.data<-matrix()
gc()

table(is.na(integrated_spatial$Th1))
integrated_spatial<-SCTransform(integrated_spatial,assay = "Spatial",method="poisson",variable.features.n = 18645)
dim(integrated_spatial@assays$SCT@data)
dim(integrated_spatial@assays$SCT@scale.data)



# tide
tide_splot<-fread("/share/home/hhhong/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-12/02.tide_analysis/total_predicted_result.xls",data.table = F)

tide_splot$sample<-sapply(strsplit(tide_splot$V1,"_"),'[[',1)
tide_splot$sample<-gsub("-","\\.",tide_splot$sample)
tide_splot$group<-gsub("\\..*","",tide_splot$sample)
tide_splot<-as.data.frame(tide_splot)


tide_splot$group<-factor(tide_splot$group,levels=c('TT','NTT'))
tide_splot$sample<-factor(tide_splot$sample,levels=c('TT.1','TT.2','TT.3','TT.4','NTT.1',"NTT.2","NTT.3","NTT.4"))
cancer_marker_score_df1_bak<-cancer_marker_score_df1



ggplot(data = tide_splot,aes(x=group,y=TIDE,fill=sample))+geom_boxplot(position = position_dodge(0.8))+
  theme_classic2()+geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)
# 显示点

ggplot(data = tide_splot,aes(x=group,y=TIDE,fill=sample))+geom_boxplot(position = position_dodge(0.8))+
  # geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),alpha=0.5)+
  geom_point(position = position_jitterdodge(jitter.width = 0.7,dodge.width = 0.8))+
  theme_classic2()+geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)


ggplot(data = tide_splot,aes(x=group,y=Dysfunction,fill=sample))+geom_boxplot(position = position_dodge(0.8))+
  theme_classic2()+geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)
# 点图
ggplot(data = tide_splot,aes(x=group,y=Dysfunction,fill=sample))+geom_boxplot(position = position_dodge(0.8))+
  geom_point(position = position_jitterdodge(jitter.width = 0.7,dodge.width = 0.8))+
  theme_classic2()+geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)


ggplot(data = tide_splot,aes(x=group,y=Exclusion,fill=sample))+geom_boxplot(position = position_dodge(0.8))+
  theme_classic2()+geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)

# 点图
ggplot(data = tide_splot,aes(x=group,y=Exclusion,fill=sample))+geom_boxplot(position = position_dodge(0.8))+
  geom_point(position = position_jitterdodge(jitter.width = 0.7,dodge.width = 0.8))+
  theme_classic2()+geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)



group1<-subset(tide_splot,group=="TT")$TIDE
group2<-subset(tide_splot,group=="NTT")$TIDE





#############################################################
library(Seurat)
library(Signac)
library(CellChat)
rm(carcinoma_regions_integrated1_Macrophage)

# 5.1）配受体对相互作用表（4个癌栓样本合并，4个非癌栓样本合并后展示细胞通讯，不需要展示单个样本的细胞通讯）

carcinoma_list=list(
  Malignant_epithelial_cells=Malignant_epithelial_cells,
  integrated_spatial=integrated_spatial,
  integrated_spatial_carcinoma_regions=integrated_spatial_carcinoma_regions,
  carcinoma_regions_integrated1=carcinoma_regions_integrated1

)


total_rds<-readRDS(file = "~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-13/data/total_rds_list.rds")
integrated_spatial<-total_rds$integrated_spatial
carcinoma_regions_integrated1<-total_rds$carcinoma_regions_integrated1


###################################################################################2024-3-13
TT<-subset(carcinoma_regions_integrated1,orig.ident%in%paste0("TT.",seq(1,4)))
NTT<-subset(carcinoma_regions_integrated1,orig.ident%in%paste0("NTT.",seq(1,4)))

carcinoma_regions_integrated1@images
data.input_list = list()
coor_list = list()
meta_list = list()
for(nm in unique(carcinoma_regions_integrated1$orig.ident)){
  tmp<-subset(carcinoma_regions_integrated1,orig.ident==nm)
  Idents(tmp)<-tmp$dominant_celltype_combined_0.1
  levels(tmp)<-names(sort(table(Idents(tmp)),decreasing = T))
  slice_name = paste0("slice1_",nm)
  tmp_image<-tmp@images[[slice_name]]
  tmp@images<-list()
  tmp@images[[slice_name]]<-tmp_image
  meta = data.frame(labels = Idents(tmp), samples = nm)
  meta_list[[nm]]<-meta
  data.input1 = Seurat::GetAssayData(tmp, slot = "data", assay = "SCT")
  data.input_list[[nm]]=data.input1
  spatial.locs1 = Seurat::GetTissueCoordinates(tmp, scale = NULL, cols = c("imagerow", "imagecol"))
  coor_list[[nm]]=spatial.locs1
}


data.input1 = Seurat::GetAssayData(TT, slot = "data", assay = "SCT")

data.input2 = Seurat::GetAssayData(NTT, slot = "data", assay = "SCT")

genes.common <- intersect(rownames(data.input1), rownames(data.input2))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ])


Idents(TT)<-TT$dominant_celltype_combined_0.1
levels(TT)<-names(sort(table(Idents(TT)),decreasing = T))

Idents(NTT)<-NTT$dominant_celltype_combined_0.1
levels(NTT)<-names(sort(table(Idents(NTT)),decreasing = T))



meta1 = data.frame(labels = Idents(TT), samples = "A1") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(NTT), samples = "A2")
meta1$samples<-TT$orig.ident
meta2$samples<-NTT$orig.ident

meta <- rbind(meta1, meta2)

rownames(meta) <- colnames(data.input)

meta$labels <- factor(meta$labels, levels = levels(TT))
meta$samples <- factor(meta$samples, levels = c(paste0("TT.",seq(1,4)),
                                                paste0("NTT.",seq(1,4))))



spatial.locs1 = Seurat::GetTissueCoordinates(TT, scale = NULL, cols = c("imagerow", "imagecol"))

spatial.locs2 = Seurat::GetTissueCoordinates(NTT, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs <- rbind(spatial.locs1, spatial.locs2)


rownames(spatial.locs) <- colnames(data.input)



total_coor=do.call(rbind,coor_list)
total_data_input <- do.call(cbind,data.input_list)

total_meta = do.call(rbind,meta_list)

total_meta$group<-gsub("\\..*","",total_meta$samples)
rownames(total_coor)<-colnames(total_data_input)



total_meta$group<-factor(total_meta$group,levels=c("TT","NTT"))
total_meta$samples<-factor(total_meta$samples,levels=c(paste0("TT.",seq(1,4)),
                                                       paste0("NTT.",seq(1,4))))


##################################################读入scalefactor.json


NTT_2 = jsonlite::fromJSON(txt = file.path("/data/ncbi_data/CR_22_07526_TS_R_SSV_1/outs/spatial", 'scalefactors_json.json'))
NTT_4 = jsonlite::fromJSON(txt = file.path("/data/ncbi_data/CR_22_15979_FP_R_JSV_1/outs/spatial", 'scalefactors_json.json'))
NTT_3 = jsonlite::fromJSON(txt = file.path("/data/ncbi_data/CR_22_14422_FP_R_JSV_1/outs/spatial", 'scalefactors_json.json'))
NTT_1 = jsonlite::fromJSON(txt = file.path("/data/ncbi_data/CR_22_15980_FP_R_JSV_1/outs/spatial", 'scalefactors_json.json'))

TT2 = jsonlite::fromJSON(txt = file.path("/data/our_data", '004-A1_scalefactors_json.json'))
TT3 = jsonlite::fromJSON(txt = file.path("/data/our_data", '310-A1_scalefactors_json.json'))
TT4 = jsonlite::fromJSON(txt = file.path("/data/our_data", '004-D1_scalefactors_json.json'))
TT1 = jsonlite::fromJSON(txt = file.path("/data/our_data", '310-D1_scalefactors_json.json'))

scale_factor_list = list(TT1=TT1,TT2=TT2,TT3=TT3,TT4=TT4,
                         NTT1=NTT_1,NTT2=NTT_2,NTT3=NTT_3,NTT4=NTT_4)

spot.size = 65
spatial.factors1_list =list()
for(nm in names(scale_factor_list)){
  tmp = scale_factor_list[[nm]]
  conversion.factor1 = spot.size/tmp$spot_diameter_fullres
  spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)
  spatial.factors1_list[[nm]]=spatial.factors1
}

# the theoretical spot size (um) in 10X Visium


spatial.factors1_df<-do.call(rbind,spatial.factors1_list)



scalefactors2 = jsonlite::fromJSON(txt = file.path("/Users/suoqinjin/Library/CloudStorage/OneDrive-Personal/works/CellChat/tutorial/spatial_imaging_data-intestinalA2", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)

spatial.factors <- rbind(spatial.factors1, spatial.factors2)
rownames(spatial.factors) <- c("A1", "A2")

rownames(spatial.factors1_df)<-c(paste0("TT.",seq(1,4)),paste0("NTT.",seq(1,4)))

###################################################################
###################################################################
###################################################################
###################################################################

rownames(total_meta)<-colnames(total_data_input)
cellchat <- createCellChat(object = total_data_input, meta = total_meta, group.by = "labels",
                           datatype = "spatial", coordinates = total_coor, spatial.factors = spatial.factors1_df)


##################################################################


CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use



cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 422

execution.time = Sys.time() - ptm


cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)



groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

pathways.show <- c("EGF")
par(mfrow=c(2,1))


netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = c("TT.2"),
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1,
                    alpha.image = 0.2, vertex.label.cex = 0)


################################################ 类似单细胞分析


setwd('/2023-10-16/01.cellchat/')
A1_004_1<-readRDS('../data/A1_004_1_cellchat.rds')
A1_310_1<-readRDS('../data/A1_310_1_cellchat.rds')
D1_004_1<-readRDS('../data/D1_004_1_cellchat.rds')
D1_310_1<-readRDS('../data/D1_310_1_cellchat.rds')
A1_004_1@meta$labels

cellchat_saptial<-function(obj){

  data.input = GetAssayData(A1_004_1, slot = "data", assay = "SCT") # normalized data matrix
  meta = data.frame(labels = A1_004_1$type, row.names = names(Idents(A1_004_1))) # manually create a dataframe consisting of the cell labels
  spatial.locs = GetTissueCoordinates(A1_004_1, scale = NULL, cols = c("imagerow", "imagecol"))
  scale_factor <-jsonlite::fromJSON('/data/20230713-V42L22-004-A1-1492164-1-outs/spatial/scalefactors_json.json')
  View(scale_factor)
  scale_factor = list(spot.diameter = 65, spot = scale_factor$spot_diameter_fullres, # these two information are required
                      fiducial = scale_factor$fiducial_diameter_fullres, hires = scale_factor$tissue_hires_scalef, lowres = scale_factor$tissue_lowres_scalef # these three information are not required
  )


  cellchat_A1_004 <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                                    datatype = "spatial", coordinates = spatial.locs,scale.factors = scale_factor)


  
  cellchat_A1_004
  CellChatDB <- CellChatDB.human

 
  CellChatDB.use <- CellChatDB
  cellchat_A1_004@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat_A1_004)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)

  cellchat
}





A1_004_1 <- updateCellChat(A1_004_1)
A1_310_1 <- updateCellChat(A1_310_1)
D1_004_1 <- updateCellChat(D1_004_1)
D1_310_1 <- updateCellChat(D1_310_1)


object.list <- list(A1_004_1 = A1_004_1, A1_310_1 = A1_310_1,D1_004_1=D1_004_1,D1_310_1=D1_310_1)

cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = T)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")



gg1 <- netVisual_heatmap(cellchat)

gg2 <- netVisual_heatmap(cellchat, measure = "weight")

gg1 + gg2



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


cellchat_saptial1<-function(obj){

  data.input = GetAssayData(obj, slot = "data", assay = "SCT") # normalized data matrix

  labels = Idents(obj)
  meta <- data.frame(labels = labels, row.names = names(labels))
  print(meta)
  print("**********************************")
  print(str(meta))
  meta$labels<-as.character(meta$labels)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")



  CellChatDB <- CellChatDB.human

  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(cellchat,type = "triMean")
 
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")



  cellchat
}
edit(identifyOverExpressedGenes)

dim(tt_group@assays$SCT@data)
dim(ntt_group@assays$SCT@data)
tt_group<-subset(carcinoma_regions_integrated1,group=="TT")
ntt_group<-subset(carcinoma_regions_integrated1,group=="NTT")
Idents(tt_group)<-tt_group$dominant_celltype_combined_0.1
Idents(ntt_group)<-ntt_group$dominant_celltype_combined_0.1
table(ntt_group$orig.ident)
levels(tt_group)<-names(sort(table(Idents(tt_group)),decreasing = T))
levels(ntt_group)<-names(sort(table(Idents(ntt_group)),decreasing = T))
Idents(tt_group)
tt_cellchat<-cellchat_saptial1(tt_group)
ntt_cellchat<-cellchat_saptial1(ntt_group)

object.list <- list(TT = tt_cellchat, NTT = ntt_cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



######################################################## 开始分析

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2



par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


cellchat@var.features$TT$features
tt_sig_lr<-cellchat@LR$TT$LRsig
ntt_sig_lr<-cellchat@LR$NTT$LRsig


# 绘制5.2图
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method='uwot')
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
rankSimilarity(cellchat, type = "functional")
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2



library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


###################################################################
###################################################################

i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- object.list[[i]]@netP$pathways
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


i = 2
# combining all the identified signaling pathways from different datasets
pathway.union <- object.list[[i]]@netP$pathways
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# 这里返回  互作的source 和target表

p1<-netVisual_bubble(cellchat, sources.use = c(7:8), targets.use = c(7:8),  comparison = c(1, 2), angle.x = 45,return.data = T)
# ?netVisual_bubble

View(data.frame(p1$communication,check.names = F,check.rows = F))
library(data.table)
fwrite(data.frame(p1$communication,check.names = F,check.rows = F),file="~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-13/5_cellchat/5.1/communication_interaction.xls",sep = '\t',col.names = T,row.names = F)


gg1 <- compareInteractions(cellchat, show.legend = F, group=c(1:10),group.facet=carcinoma_regions_integrated1$dominant_celltype_combined_0.1)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

tt_count<-cellchat@net$TT$count
ntt_count<-cellchat@net$NTT$count





total_count<-tt_count + ntt_count

total_count<-data.frame(total_count,check.names = F,check.rows = F)
total_count$celltype<-rownames(total_count)
###################################################
###################################################

aa<-total_count%>%data.frame(check.names = F,check.rows = F)%>%tidyr::gather(-c("celltype"))


library(tidyverse)


long_df <- total_count %>%
  pivot_longer(
    cols = -celltype, #
    names_to = "interaction", # 
    values_to = "value" # 
  )



long_df <- long_df %>%
  group_by(celltype) %>%
  mutate(mean_value = mean(value)) %>%
  ungroup() %>%
  arrange(desc(mean_value), desc(value))


ggplot(long_df,aes(x=reorder(celltype,-as.numeric(mean_value)),y= reorder(value,-as.numeric(value)),fill=interaction))+
  geom_bar(stat = 'identity',position = position_dodge(0.75))

color_select=c("red",color_panel[1:3],color_panel[25],color_panel[22],color_panel[7],color_panel[17],color_panel[19:20])
ggplot(long_df,aes(x=reorder(celltype,-as.numeric(mean_value)),y= value,fill=interaction))+
  geom_bar(stat = 'identity',position = position_dodge(0.75))+
  scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                             Macrophage = color_select[2], Fibroblast_cells = color_select[3] ,
                             CD4._T_cells = color_select[4], Endothelial_cells = color_select[5],
                             B_cells=color_select[6],Nonmalignant_epithelial_cells=color_select[7],
                             SMC=color_select[8],CD8._T_cells=color_select[9],Mast_cells=color_select[10]))+
  theme(axis.text.x =  element_text(angle = 45))+
  theme_classic()+xlab(label = "celltype")+

  ylab(label = "interaction number")





tt_weight$celltype<-NULL
tt_weight11<-tt_weight[upper.tri(tt_weight)]


tt_weight<-cellchat@net$TT$weight
ntt_weight<-cellchat@net$NTT$weight
total_weight<-tt_weight+ntt_weight
total_weight<-data.frame(total_weight,check.names = F,check.rows = F)
total_weight$celltype<-rownames(total_weight)


long_df1 <- total_weight %>%
  pivot_longer(
    cols = -celltype, # 
    names_to = "interaction", # 
    values_to = "value" # 
  )


ggplot(long_df1,aes(x=reorder(celltype,-as.numeric(value)),y= value,fill=interaction))+
  geom_bar(stat = 'identity',position = position_dodge(0.75))+
  scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                             Macrophage = color_select[2], Fibroblast_cells = color_select[3] ,
                             CD4._T_cells = color_select[4], Endothelial_cells = color_select[5],
                             B_cells=color_select[6],Nonmalignant_epithelial_cells=color_select[7],
                             SMC=color_select[8],CD8._T_cells=color_select[9],Mast_cells=color_select[10]))+
  theme(axis.text.x =  element_text(angle = 45))+
  theme_classic()+xlab(label = "celltype")+

  ylab(label = "interaction weight")



##################################################################
##################################################################



tt_weight<-tt_weight%>%data.frame(check.names = F,check.rows = F)%>%mutate(celltype=row.names(.))
ntt_weight<-ntt_weight%>%data.frame(check.names = F,check.rows = F)%>%mutate(celltype=row.names(.))

tt_weight_long_df <- tt_weight%>%pivot_longer(
  cols = -celltype,
  names_to = "interaction",
  values_to = "value"

)

ntt_weight_long_df <- ntt_weight%>%pivot_longer(
  cols = -celltype,
  names_to = "interaction",
  values_to = "value"

)
# TT weight
color_select<-c("red","#F3B1A0","#58A4C3","#6778AE","#BD956A","#8C549C","#53A85F","#AB3282","#F1BB72","#57C3F3")
setwd("5_cellchat/5.5_adjust/TT")
for(i in unique(tt_weight_long_df$celltype)){

  tmp<-subset(tt_weight_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  tmp$interaction<-factor(tmp$interaction,levels=tmp$interaction)
  p1<-ggplot(tmp,aes(x=celltype,y= value,fill=interaction))+
    geom_bar(stat = 'identity',position = position_dodge(0.95))+

    scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                               Macrophage = color_select[2],
                               Fibroblast_cells = color_select[3] ,
                               CD4._T_cells = color_select[4],
                               Endothelial_cells = color_select[5],
                               B_cells=color_select[6],
                               Nonmalignant_epithelial_cells=color_select[7],
                               SMC=color_select[8],
                               CD8._T_cells=color_select[9],
                               Mast_cells=color_select[10]))+

    theme_classic()+xlab(label = "celltype")+
    ylab(label = "interaction weight")+
    theme(axis.text.x =  element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    geom_text(aes(label=paste(celltype, "-", interaction),y=0),
              position=position_dodge(0.95), angle=90,hjust=0)

  ggsave(p1,filename = paste0(i,"_weight.pdf"),width = 12,height = 6)

}

# NTT weight
color_select<-c("red","#F3B1A0","#58A4C3","#6778AE","#BD956A","#8C549C","#53A85F","#AB3282","#F1BB72","#57C3F3")
for(i in unique(ntt_weight_long_df$celltype)){

  tmp<-subset(ntt_weight_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  tmp$interaction<-factor(tmp$interaction,levels=tmp$interaction)
  p1<-ggplot(tmp,aes(x=celltype,y= value,fill=interaction))+
    geom_bar(stat = 'identity',position = position_dodge(0.95))+

    scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                               Macrophage = color_select[2],
                               Fibroblast_cells = color_select[3] ,
                               CD4._T_cells = color_select[4],
                               Endothelial_cells = color_select[5],
                               B_cells=color_select[6],
                               Nonmalignant_epithelial_cells=color_select[7],
                               SMC=color_select[8],
                               CD8._T_cells=color_select[9],
                               Mast_cells=color_select[10]))+

    theme_classic()+xlab(label = "celltype")+
    ylab(label = "interaction weight")+
    theme(axis.text.x =  element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    geom_text(aes(label=paste(celltype, "-", interaction),y=0),
              position=position_dodge(0.95), angle=90,hjust=0)

  ggsave(p1,filename = paste0(i,"_weight.pdf"),width = 12,height = 6)

}



############################################################################
tt_count<-tt_count%>%data.frame(check.names = F,check.rows = F)%>%mutate(celltype=row.names(.))
ntt_count<-ntt_count%>%data.frame(check.names = F,check.rows = F)%>%mutate(celltype=row.names(.))

tt_count_long_df <- tt_count%>%pivot_longer(
  cols = -celltype,
  names_to = "interaction",
  values_to = "value"

)

ntt_count_long_df <- ntt_count%>%pivot_longer(
  cols = -celltype,
  names_to = "interaction",
  values_to = "value"

)

# TT number

############################################################################
color_select<-c("red","#F3B1A0","#58A4C3","#6778AE","#BD956A","#8C549C","#53A85F","#AB3282","#F1BB72","#57C3F3")
setwd("5_cellchat/5.5_adjust/TT")
for(i in unique(tt_count_long_df$celltype)){

  tmp<-subset(tt_count_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  tmp$interaction<-factor(tmp$interaction,levels=tmp$interaction)
  p1<-ggplot(tmp,aes(x=celltype,y= value,fill=interaction))+
    geom_bar(stat = 'identity',position = position_dodge(0.95))+

    scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                               Macrophage = color_select[2],
                               Fibroblast_cells = color_select[3] ,
                               CD4._T_cells = color_select[4],
                               Endothelial_cells = color_select[5],
                               B_cells=color_select[6],
                               Nonmalignant_epithelial_cells=color_select[7],
                               SMC=color_select[8],
                               CD8._T_cells=color_select[9],
                               Mast_cells=color_select[10]))+

    theme_classic()+xlab(label = "celltype")+
    ylab(label = "interaction number")+
    theme(axis.text.x =  element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    geom_text(aes(label=paste(celltype, "-", interaction),y=0),
              position=position_dodge(0.95), angle=90,hjust=0)

  ggsave(p1,filename = paste0(i,"_number.pdf"),width = 12,height = 6)

}




##################################
setwd("/5_cellchat/5.5_adjust/NTT")
for(i in unique(ntt_count_long_df$celltype)){

  tmp<-subset(ntt_count_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  tmp$interaction<-factor(tmp$interaction,levels=tmp$interaction)
  p1<-ggplot(tmp,aes(x=celltype,y= value,fill=interaction))+
    geom_bar(stat = 'identity',position = position_dodge(0.95))+

    scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                               Macrophage = color_select[2],
                               Fibroblast_cells = color_select[3] ,
                               CD4._T_cells = color_select[4],
                               Endothelial_cells = color_select[5],
                               B_cells=color_select[6],
                               Nonmalignant_epithelial_cells=color_select[7],
                               SMC=color_select[8],
                               CD8._T_cells=color_select[9],
                               Mast_cells=color_select[10]))+

    theme_classic()+xlab(label = "celltype")+
    ylab(label = "interaction number")+
    theme(axis.text.x =  element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    geom_text(aes(label=paste(celltype, "-", interaction),y=0),
              position=position_dodge(0.95), angle=90,hjust=0)

  ggsave(p1,filename = paste0(i,"_number.pdf"),width = 12,height = 6)

}




############################################################# 绘制到一张图上
# ntt_weight

tt_weight_long_df_order<-data.frame()
for(i in unique(tt_weight_long_df$celltype)){

  tmp<-subset(tt_weight_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  tt_weight_long_df_order<-rbind(tt_weight_long_df_order,tmp)
}

ntt_weight_long_df_order<-data.frame()
for(i in unique(ntt_weight_long_df$celltype)){

  tmp<-subset(ntt_weight_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  ntt_weight_long_df_order<-rbind(ntt_weight_long_df_order,tmp)
}



tt_count_long_df_order<-data.frame()
for(i in unique(tt_count_long_df$celltype)){

  tmp<-subset(tt_count_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  tt_count_long_df_order<-rbind(tt_count_long_df_order,tmp)
}



ntt_count_long_df_order<-data.frame()
for(i in unique(ntt_count_long_df$celltype)){

  tmp<-subset(ntt_count_long_df,celltype==i)
  tmp<-tmp[order(tmp$value,decreasing = T),]
  ntt_count_long_df_order<-rbind(ntt_count_long_df_order,tmp)
}




# 绘制到一张图上


p1<-ggplot(ntt_weight_long_df_order,aes(x=reorder(celltype,-as.numeric(value)),y= value,fill=interaction))+
  geom_bar(stat = 'identity',position = position_dodge(0.95))+

  scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                             Macrophage = color_select[2],
                             Fibroblast_cells = color_select[3] ,
                             CD4._T_cells = color_select[4],
                             Endothelial_cells = color_select[5],
                             B_cells=color_select[6],
                             Nonmalignant_epithelial_cells=color_select[7],
                             SMC=color_select[8],
                             CD8._T_cells=color_select[9],
                             Mast_cells=color_select[10]))+

  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction weight")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)

##############################################################################################
##############################################################################################



p1<-ggplot(ntt_count_long_df_order,aes(x=reorder(celltype,-as.numeric(value)),y= value,fill=interaction))+
  geom_bar(stat = 'identity',position = position_dodge(0.95))+

  scale_fill_manual(values=c(Malignant_epithelial_cells = "red",
                             Macrophage = color_select[2],
                             Fibroblast_cells = color_select[3] ,
                             CD4._T_cells = color_select[4],
                             Endothelial_cells = color_select[5],
                             B_cells=color_select[6],
                             Nonmalignant_epithelial_cells=color_select[7],
                             SMC=color_select[8],
                             CD8._T_cells=color_select[9],
                             Mast_cells=color_select[10]))+

  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction number")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)












tt_group<-subset(carcinoma_regions_integrated1,group=="TT")
ntt_group<-subset(carcinoma_regions_integrated1,group=="NTT")
Idents(tt_group)<-tt_group$group
Idents(ntt_group)<-ntt_group$group

levels(tt_group)<-names(sort(table(Idents(tt_group)),decreasing = T))
levels(ntt_group)<-names(sort(table(Idents(ntt_group)),decreasing = T))

tt_cellchat<-cellchat_saptial(tt_group)
ntt_cellchat<-cellchat_saptial(ntt_group)


i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- object.list[[i]]@netP$pathways
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

?netAnalysis_signalingRole_heatmap




##############################################################################cellchat 2024-3-16
setwd('2024-3-16/5.5_adjust/TT')


tt_count$celltype<-NULL
ntt_count$celltype<-NULL

tt_count_t<-t(tt_count)
tt_count_process<-tt_count_t+tt_count
tt_count_process[upper.tri(tt_count_process)] <- 0



tt_count_process<-as.matrix(tt_count_process)
diag_value<-diag(as.matrix(tt_count_process))

diag(tt_count_process) <- diag_value / 2
diag(tt_count_process)


tt_count_process1<-data.frame(tt_count_process,check.rows = F,check.names = F)
tt_count_process1$celltype<-rownames(tt_count_process1)


########################################
tt_count_process1_long <- tt_count_process1 %>%
  pivot_longer(
    cols = -celltype, # 
    names_to = "interaction", # 
    values_to = "value" # 
  )

tt_count_process1_long1<-subset(tt_count_process1_long,value!=0)


tt_count_process1_long1$label<-paste0(tt_count_process1_long1$celltype,"-",tt_count_process1_long1$interaction)

ggplot(tt_count_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+
  scale_fill_viridis_d()+
  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction number")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)


# 第二种配色
ggplot(tt_count_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+

  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction number")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)

fwrite(data.frame(label=tt_count_process1_long1$label,value=tt_count_process1_long1$value),file="/TT_count.xls",sep = '\t')


##########################ntt
##########################ntt
##########################ntt
ntt_count$celltype<-NULL


ntt_count_t<-t(ntt_count)
ntt_count_process<-ntt_count_t+ntt_count


ntt_count_process[upper.tri(ntt_count_process)] <- 0



ntt_count_process<-as.matrix(ntt_count_process)
diag_value<-diag(as.matrix(ntt_count_process))

diag(ntt_count_process) <- diag_value / 2



ntt_count_process1<-data.frame(ntt_count_process,check.rows = F,check.names = F)
ntt_count_process1$celltype<-rownames(ntt_count_process1)


ntt_count_process1_long <- ntt_count_process1 %>%
  pivot_longer(
    cols = -celltype, # 
    names_to = "interaction", # 
    values_to = "value" # 
  )

ntt_count_process1_long1<-subset(ntt_count_process1_long,value!=0)



ntt_count_process1_long1$label<-paste0(ntt_count_process1_long1$celltype,"-",ntt_count_process1_long1$interaction)


# 第一种配色
ggplot(ntt_count_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+
  scale_fill_viridis_d()+
  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction number")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)


# 第二种配色
ggplot(ntt_count_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+

  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction number")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)



fwrite(data.frame(label=ntt_count_process1_long1$label,value=ntt_count_process1_long1$value),file="~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-18/01.stat/NTT_count.xls",sep = '\t')


##############################################################
############################################################## weight

tt_weight$celltype<-NULL


tt_weight_t<-t(tt_weight)
tt_weight_process<-tt_weight_t+tt_weight

tt_weight_process[upper.tri(tt_weight_process)] <- 0


tt_weight_process<-as.matrix(tt_weight_process)
diag_value<-diag(as.matrix(tt_weight_process))

diag(tt_weight_process) <- diag_value / 2

tt_weight_process1<-data.frame(tt_weight_process,check.rows = F,check.names = F)
tt_weight_process1$celltype<-rownames(tt_weight_process1)

########################################
########################################
########################################
########################################
tt_weight_process1_long <- tt_weight_process1 %>%
  pivot_longer(
    cols = -celltype, # 
    names_to = "interaction", # 
    values_to = "value" #
  )

tt_weight_process1_long1<-subset(tt_weight_process1_long,value!=0)


tt_weight_process1_long1$label<-paste0(tt_weight_process1_long1$celltype,"-",tt_weight_process1_long1$interaction)
# 第一种配色
ggplot(tt_weight_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+
  scale_fill_viridis_d()+
  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction strength")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)


# 第二种配色
ggplot(tt_weight_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+

  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction strength")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)

fwrite(data.frame(label=tt_weight_process1_long1$label,value=tt_weight_process1_long1$value),file="~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-18/01.stat/TT_weight.xls",sep = '\t')


#####################################################################ntt weight
#####################################################################ntt weight


ntt_weight$celltype<-NULL




ntt_weight_t<-t(ntt_weight)
ntt_weight_process<-ntt_weight_t+ntt_weight

ntt_weight_process[upper.tri(ntt_weight_process)] <- 0

ntt_weight_process<-as.matrix(ntt_weight_process)
diag_value<-diag(as.matrix(ntt_weight_process))

diag(ntt_weight_process) <- diag_value / 2



ntt_weight_process1<-data.frame(ntt_weight_process,check.rows = F,check.names = F)
ntt_weight_process1$celltype<-rownames(ntt_weight_process1)

########################################
########################################
########################################
########################################


ntt_weight_process1_long <- ntt_weight_process1 %>%
  pivot_longer(
    cols = -celltype, # 
    names_to = "interaction", # 
    values_to = "value" # 
  )

ntt_weight_process1_long1<-subset(ntt_weight_process1_long,value!=0)


ntt_weight_process1_long1$label<-paste0(ntt_weight_process1_long1$celltype,"-",ntt_weight_process1_long1$interaction)

# 第一种配色
ggplot(ntt_weight_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+
  scale_fill_viridis_d()+
  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction strength")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)


fwrite(data.frame(label=ntt_weight_process1_long1$label,value=ntt_weight_process1_long1$value),file="~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-18/01.stat/NTT_weight.xls",sep = '\t')

# 第二种配色
ggplot(ntt_weight_process1_long1,aes(x=reorder(label,-as.numeric(value)),y= value,fill=label))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+

  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction strength")+
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0)













###################################################################2024-3-19
##################################################### 统计互作对（每一个cluster的数量和强度）

################################################ barplot stat
procee_weight2df<-function(weight_df){
  ntt_weight<-weight_df
  ntt_weight_t<-t(ntt_weight)
  ntt_weight_process<-ntt_weight_t+ntt_weight

  ntt_weight_process[upper.tri(ntt_weight_process)] <- 0
  ntt_weight_process<-as.matrix(ntt_weight_process)
  diag_value<-diag(as.matrix(ntt_weight_process))

  diag(ntt_weight_process) <- diag_value / 2

  ntt_weight_process1<-data.frame(ntt_weight_process,check.rows = F,check.names = F)
  ntt_weight_process1$celltype<-rownames(ntt_weight_process1)

#####################################################################

  ntt_weight_process1_long <- ntt_weight_process1 %>%
    pivot_longer(
      cols = -celltype, # 
      names_to = "interaction", # 
      values_to = "value" # 
    )

  ntt_weight_process1_long1<-subset(ntt_weight_process1_long,value!=0)


  ntt_weight_process1_long1$label<-paste0(ntt_weight_process1_long1$celltype,"-",ntt_weight_process1_long1$interaction)
  ntt_weight_process1_long1
}


procee_count2df<-function(count_df){
  ntt_weight<-count_df
  ntt_weight_t<-t(ntt_weight)
  ntt_weight_process<-ntt_weight_t+ntt_weight

  ntt_weight_process[upper.tri(ntt_weight_process)] <- 0
  ntt_weight_process<-as.matrix(ntt_weight_process)
  diag_value<-diag(as.matrix(ntt_weight_process))

  diag(ntt_weight_process) <- diag_value / 2

  ntt_weight_process1<-data.frame(ntt_weight_process,check.rows = F,check.names = F)
  ntt_weight_process1$celltype<-rownames(ntt_weight_process1)

  #####################################################################

  ntt_weight_process1_long <- ntt_weight_process1 %>%
    pivot_longer(
      cols = -celltype, # 保留celltype列为标识符，转换其他所有列
      names_to = "interaction", # 新列，用于保存原来的列名（细胞类型）
      values_to = "value" # 新列，用于保存原来列中的数值
    )

  ntt_weight_process1_long1<-subset(ntt_weight_process1_long,value!=0)


  ntt_weight_process1_long1$label<-paste0(ntt_weight_process1_long1$celltype,"-",ntt_weight_process1_long1$interaction)
  ntt_weight_process1_long1
}



tt_weight<-cellchat@net$TT$weight
ntt_weight<-cellchat@net$NTT$weight
tt_weight_process1_long1<-procee_weight2df(tt_weight)
ntt_weight_process1_long1<-procee_weight2df(ntt_weight)



###################################weight
combined_weight_long_order$source<-factor(combined_weight_long_order$source,levels=c('TT','NTT'))

combined_weight_long_order<-process_after_stat(tt_weight_process1_long1,ntt_weight_process1_long1)

# combined_weight_long_order_total<-combined_weight_long_order
fwrite(data.frame(value=combined_weight_long_order_total$value,label=combined_weight_long_order_total$label,source=combined_weight_long_order_total$source),file = "combined_weight.xls",sep = '\t')


# 第一种配色
p<-ggplot(combined_weight_long_order_total,aes(x=label,y= value,fill=source))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+
  scale_fill_manual(values = c(TT="#AB3282",NTT="#58A4C3"))+
  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction weight")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        )+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0,size=2.5)

# 加边框
p + theme(plot.margin=margin(t=1, r=5, b=1, l=4, unit="cm"))



##########################################count
cellchat<-readRDS('/data/total_db_cellchat.rds')


tt_count<-cellchat@net$TT$count
ntt_count<-cellchat@net$NTT$count
tt_count_process1_long1<-procee_count2df(tt_count)
ntt_count_process1_long1<-procee_count2df(ntt_count)

View(tt_count_process1_long1)


process_after_stat<-function(tt_count_process1_long1,ntt_count_process1_long1){
  tt_weight_process1_long1<-tt_count_process1_long1
  ntt_weight_process1_long1<- ntt_count_process1_long1
  tt_weight_process1_long1$label1<-paste0(tt_weight_process1_long1$label,"@TT")
  ntt_weight_process1_long1$label1<-paste0(ntt_weight_process1_long1$label,"@NTT")

  tt_weight_process1_long1_order<-tt_weight_process1_long1[order(tt_weight_process1_long1$value,decreasing = T),]

  tt_weight_process1_long1_order$source<-'df1'
  ntt_weight_process1_long1$source<-'df2'

  # 正好反过来
  ntt_weight_process1_long1_order<-ntt_weight_process1_long1[match(tt_weight_process1_long1_order$label,ntt_weight_process1_long1$label),]

  combined_weight_long_order<-data.frame(check.names = F,check.rows = F)



  for(i in 1:nrow(tt_weight_process1_long1_order)){
    combined_weight_long_order<-rbind(combined_weight_long_order,rbind(tt_weight_process1_long1_order[i,],ntt_weight_process1_long1_order[i,]))
  }



  combined_weight_long_order$source<-sapply(strsplit(combined_weight_long_order$label1,"@"),'[[',2)
  combined_weight_long_order$label<-factor(combined_weight_long_order$label,levels=unique(combined_weight_long_order$label))
  ###################################weight
  combined_weight_long_order$source<-factor(combined_weight_long_order$source,levels=c('TT','NTT'))
  combined_weight_long_order
}



combined_count_long_order<-process_after_stat(tt_count_process1_long1,ntt_count_process1_long1)



setwd('/2024-3-21/5.5/total_db_stat')
fwrite(data.frame(value=combined_count_long_order$value,label=combined_count_long_order$label,source=combined_count_long_order$source),file = "combined_count.xls",sep = '\t')


p<-ggplot(combined_count_long_order,aes(x=label,y= value,fill=source))+
  geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+
  scale_fill_manual(values = c(TT="#AB3282",NTT="#58A4C3"))+
  theme_classic()+xlab(label = "celltype")+
  ylab(label = "interaction count")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
  )+
  geom_text(aes(label=paste(celltype, "-", interaction),y=0),
            position=position_dodge(0.95), angle=90,hjust=0,size=2.5)



# 加边框
p + theme(plot.margin=margin(t=1, r=5, b=1, l=4, unit="cm"))


#####################################################
##################################################### 做Secreted Signaling
cellchat_saptial2<-function(obj,db="Secreted Signaling"){
  # rm(data.input)
  data.input = GetAssayData(obj, slot = "data", assay = "SCT") # normalized data matrix

  labels = Idents(obj)
  meta <- data.frame(labels = labels, row.names = names(labels))
  print(meta)
  print("**********************************")
  print(str(meta))
  meta$labels<-as.character(meta$labels)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")



  CellChatDB <- CellChatDB.human

  CellChatDB.use <- subsetDB(CellChatDB, search = db)
  # CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(cellchat,type = "triMean")
  # cellchat <- computeCommunProb(cellchat,type = "truncatedMean", trim = 0.1,
  #                               distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
  #                               contact.dependent = TRUE, contact.range = 100)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")



  cellchat
}

#################################################### 做统计
plot_interaction_bar<-function(data_df,title="count"){
  p<-ggplot(data_df,aes(x=label,y= value,fill=source))+
    geom_bar(stat = 'identity',position = position_dodge(0.95),width = 0.8)+
    scale_fill_manual(values = c(TT="#AB3282",NTT="#58A4C3"))+
    theme_classic()+xlab(label = "celltype")+
    ylab(label = paste0("interaction ",title))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
    )+
    geom_text(aes(label=paste(celltype, "-", interaction),y=0),
              position=position_dodge(0.95), angle=90,hjust=0,size=2.5)



  # 加边框
 p<- p + theme(plot.margin=margin(t=1, r=5, b=1, l=4, unit="cm"))
  p
}



setwd('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-21/5.5/Secreted_Signaling')
tt_count<-cellchat_Secreted_Signaling@net$TT$count
ntt_count<-cellchat_Secreted_Signaling@net$NTT$count
tt_count_process1_long1<-procee_count2df(tt_count)
ntt_count_process1_long1<-procee_count2df(ntt_count)


combined_count_long_order_Secreted_Signaling<-process_after_stat(tt_count_process1_long1,ntt_count_process1_long1)




set_total_order_into_split_db_count<-function(count_df){
  combined_count_long_order_Secreted_Signaling <-count_df
  combined_count_long_order_total1<-combined_count_long_order_total[combined_count_long_order_total$label%in%combined_count_long_order_Secreted_Signaling$label,]


  combined_count_long_order_Secreted_Signaling<-combined_count_long_order_Secreted_Signaling%>%data.frame(check.rows = F,check.names = F)
  combined_count_long_order_Secreted_Signaling1<-combined_count_long_order_Secreted_Signaling[match(combined_count_long_order_total1$label1,combined_count_long_order_Secreted_Signaling$label1),]
  combined_count_long_order_Secreted_Signaling1$label<-as.character(combined_count_long_order_Secreted_Signaling1$label)
  combined_count_long_order_Secreted_Signaling1$label<-factor(combined_count_long_order_Secreted_Signaling1$label,levels=unique(as.character(combined_count_long_order_Secreted_Signaling1$label)))
  combined_count_long_order_Secreted_Signaling1
}


combined_count_long_order_Secreted_Signaling1<-set_total_order_into_split_db_count(combined_count_long_order_Secreted_Signaling)

p1<-plot_interaction_bar(data_df = combined_count_long_order_Secreted_Signaling1,title = 'count')


fwrite(data.frame(value=combined_count_long_order$value,label=combined_count_long_order$label,source=combined_count_long_order$source),file = "combined_count.xls",sep = '\t')

###########################weight

tt_weight<-cellchat_Secreted_Signaling@net$TT$weight
ntt_weight<-cellchat_Secreted_Signaling@net$NTT$weight
tt_weight_process1_long1<-procee_count2df(tt_weight)
ntt_weight_process1_long1<-procee_count2df(ntt_weight)


set_total_order_into_split_db_weight<-function(weight_df){
  combined_weight_long_order_Secreted_Signaling <-weight_df
  combined_weight_long_order_total1<-combined_weight_long_order_total[combined_weight_long_order_total$label%in%combined_weight_long_order_Secreted_Signaling$label,]


  combined_weight_long_order_Secreted_Signaling<-combined_weight_long_order_Secreted_Signaling%>%data.frame(check.rows = F,check.names = F)
  combined_weight_long_order_Secreted_Signaling1<-combined_weight_long_order_Secreted_Signaling[match(combined_weight_long_order_total1$label1,combined_weight_long_order_Secreted_Signaling$label1),]
  combined_weight_long_order_Secreted_Signaling1$label<-as.character(combined_weight_long_order_Secreted_Signaling1$label)
  combined_weight_long_order_Secreted_Signaling1$label<-factor(combined_weight_long_order_Secreted_Signaling1$label,levels=unique(as.character(combined_weight_long_order_Secreted_Signaling1$label)))
  combined_weight_long_order_Secreted_Signaling1
}



combined_weight_long_order_Secreted_Signaling<-process_after_stat(tt_weight_process1_long1,ntt_weight_process1_long1)

combined_weight_long_order_Secreted_Signaling1<-set_total_order_into_split_db_weight(combined_weight_long_order_Secreted_Signaling)


p1<-plot_interaction_bar(data_df = combined_weight_long_order_Secreted_Signaling1,title = 'weight')

p1
fwrite(data.frame(value=combined_weight_long_order$value,label=combined_weight_long_order$label,source=combined_weight_long_order$source),file = "combined_weight.xls",sep = '\t')


#####################################################
##################################################### 做ECM-Receptor Interaction
showDatabaseCategory(CellChatDB.human)
tt_cellchat<-cellchat_saptial2(tt_group,db = "ECM-Receptor")
ntt_cellchat<-cellchat_saptial2(ntt_group,db="ECM-Receptor")

object.list <- list(TT = tt_cellchat, NTT = ntt_cellchat)
cellchat_ECM_Receptor <- mergeCellChat(object.list, add.names = names(object.list))

setwd("/2024-3-21/5.5/ECM_Receptor")
####################################################

tt_count<-cellchat_ECM_Receptor@net$TT$count
ntt_count<-cellchat_ECM_Receptor@net$NTT$count
tt_count_process1_long1<-procee_count2df(tt_count)
ntt_count_process1_long1<-procee_count2df(ntt_count)



combined_count_long_order_ECM_Receptor<-process_after_stat(tt_count_process1_long1,ntt_count_process1_long1)

combined_count_long_order_ECM_Receptor1<-set_total_order_into_split_db_count(combined_count_long_order_ECM_Receptor)


p1<-plot_interaction_bar(data_df = combined_count_long_order_ECM_Receptor1,title = 'count')

p1

fwrite(data.frame(value=combined_count_long_order$value,label=combined_count_long_order$label,source=combined_count_long_order$source),file = "combined_count.xls",sep = '\t')
###################################################weight
tt_weight<-cellchat_ECM_Receptor@net$TT$weight
ntt_weight<-cellchat_ECM_Receptor@net$NTT$weight
tt_weight_process1_long1<-procee_count2df(tt_weight)
ntt_weight_process1_long1<-procee_count2df(ntt_weight)


combined_weight_long_order_ECM_Receptor<-process_after_stat(tt_weight_process1_long1,ntt_weight_process1_long1)

combined_weight_long_order_ECM_Receptor1<-set_total_order_into_split_db_weight(combined_weight_long_order_ECM_Receptor)



p1<-plot_interaction_bar(data_df = combined_weight_long_order_ECM_Receptor1,title = 'weight')

p1



fwrite(data.frame(value=combined_weight_long_order$value,label=combined_weight_long_order$label,source=combined_weight_long_order$source),file = "combined_weight.xls",sep = '\t')




####################################################
#################################################### 统计
tt_count<-cellchat_Cell_Cell_Contact@net$TT$count
ntt_count<-cellchat_Cell_Cell_Contact@net$NTT$count
tt_count_process1_long1<-procee_count2df(tt_count)
ntt_count_process1_long1<-procee_count2df(ntt_count)



combined_count_long_order_Cell_Cell_Contact<-process_after_stat(tt_count_process1_long1,ntt_count_process1_long1)


combined_count_long_order_Cell_Cell_Contact1<-set_total_order_into_split_db_count(combined_count_long_order_Cell_Cell_Contact)

p1<-plot_interaction_bar(data_df = combined_count_long_order_Cell_Cell_Contact1,title = 'count')

p1

fwrite(data.frame(value=combined_count_long_order$value,label=combined_count_long_order$label,source=combined_count_long_order$source),file = "combined_count.xls",sep = '\t')
###################################################weight
tt_weight<-cellchat_Cell_Cell_Contact@net$TT$weight
ntt_weight<-cellchat_Cell_Cell_Contact@net$NTT$weight
tt_weight_process1_long1<-procee_count2df(tt_weight)
ntt_weight_process1_long1<-procee_count2df(ntt_weight)


combined_weight_long_order_Cell_Cell_Contact<-process_after_stat(tt_weight_process1_long1,ntt_weight_process1_long1)

combined_weight_long_order_Cell_Cell_Contact1<-set_total_order_into_split_db_weight(combined_weight_long_order_Cell_Cell_Contact)

p1<-plot_interaction_bar(data_df = combined_weight_long_order_Cell_Cell_Contact1,title = 'weight')

p1



fwrite(data.frame(value=combined_weight_long_order$value,label=combined_weight_long_order$label,source=combined_weight_long_order$source),file = "combined_weight.xls",sep = '\t')





########################################将所有的顺序按照total cellchat数据库的顺序进行排序



netVisual_chord_gene(cellchat, sources.use = c(3), targets.use = 2, legend.pos.x = 15)


par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(3), targets.use = c(2),  title.name = paste0("signaling - ", names(object.list)[i]), legend.pos.x = 10)
}



par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(8), targets.use = c(2),  title.name = paste0("signaling - ", names(object.list)[i]), legend.pos.x = 10)
}


par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(5), targets.use = c(4),  title.name = paste0("signaling - ", names(object.list)[i]), legend.pos.x = 10)
}


par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(6), targets.use = c(4),  title.name = paste0("signaling - ", names(object.list)[i]), legend.pos.x = 10)
}


c(c(3,8,5,6),c(2,2,4,4))

gg2 <- netVisual_bubble(cellchat, sources.use = c(3,8,5,6), targets.use = c(2,2,4,4), comparison = c(1,2),  angle.x = 90, remove.isolate = T)

gg2

gg2 <- netVisual_bubble(cellchat, sources.use =c(2,2,4,4), targets.use = c(3,8,5,6), comparison = c(1,2),  angle.x = 90, remove.isolate = T)

gg2



gg2 <- netVisual_bubble(cellchat, sources.use =c(3,2,8,2,5,4,6,4), targets.use = c(2,3,2,8,4,5,4,6), comparison = c(1,2),  angle.x = 90, remove.isolate = T)

gg2




gg2 <- netVisual_bubble(cellchat, sources.use =c(c(3,8,5,6),c(2,2,4,4)), targets.use = c(c(3,8,5,6),c(2,2,4,4)), comparison = c(1,2),  angle.x = 90, remove.isolate = T)

gg2
###############################气泡图
# 以Malignant_epithelial_cells作为作为信号发送者，到所有细胞类型的配受体对（TT和NTT分别展示）
# 以Malignant_epithelial_cells作为作为信号接收者，到所有细胞类型的配受体对（TT和NTT分别展示）
# 上面共4张图

levels(cellchat@idents$TT)
gg2 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:10), comparison = c(1,2),  angle.x = 90, remove.isolate = T)

View(gg2$data)


fwrite(data.frame(gg2$data),file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-27/5.3/aligant_epithelial_为发送者.xls',sep = '\t',col.names = T)

levels(cellchat@idents$TT)
gg3 <- netVisual_bubble(cellchat, sources.use = c(1:10), targets.use = c(7), comparison = c(1,2),  angle.x = 90, remove.isolate = T)

View(gg3$data)
fwrite(data.frame(gg3$data),file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-27/5.3/maligant_epithelial_为接收者.xls',sep = '\t',col.names = T)


ggsave(gg3,filename = 'tt2.pdf',width = 25,height = 45)


# Macrophage
# CD8._T_cells
# CD4._T_cells
# Endothelial_cells

ff1<-list(a=c(6,1:10),b=c(3,1:10),d=c(2,1:10),f=c(4,1:10))
cluster_name = levels(cellchat@idents$TT)
for(nm in names(ff1)){
  source_num = ff1[[nm]][1]
  target_num =ff1[[nm]][2:10]
  cls_name=cluster_name[source_num]
  gg3 <- netVisual_bubble(cellchat, sources.use = source_num, targets.use = target_num, comparison = c(1,2),  angle.x = 90, remove.isolate = T)
  ggsave(gg3,filename = paste0('/2024-3-27/5.3/',cls_name, '_为发送者.pdf'),limitsize = F,width = 20,height = 80)

  fwrite(data.frame(gg3$data),file = paste0('/2024-3-27/5.3/',cls_name, '_为发送者.xls'),sep = '\t',col.names = T)
  #####################################################

  gg4 <- netVisual_bubble(cellchat, sources.use = target_num, targets.use = source_num, comparison = c(1,2),  angle.x = 90, remove.isolate = T)
  ggsave(gg4,filename = paste0('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-27/5.3/',cls_name, '_为接收者.pdf'),limitsize = F,width = 20,height = 80)
  fwrite(data.frame(gg4$data),file = paste0('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-27/5.3/',cls_name, '_为接收者.xls'),sep = '\t',col.names = T)

}

class(gg3)

###########################################################2024-3-28
###########################################################2024-3-28



gg2 <- netVisual_bubble(cellchat, sources.use = c(2,3), targets.use = c(3,2), comparison = c(1,2), angle.x = 45, remove.isolate = T,thresh = 0.01)


fwrite(data.frame(gg2$data),file = paste0('/5.1/CD8_CD4_total.xls'),sep = '\t',col.names = T)

# split

gg2 <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(2), comparison = c(1,2),  angle.x = 45, remove.isolate = T,thresh = 0.01)

View(gg2$data)
fwrite(data.frame(gg2$data),file = paste0('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-28/5.1/split1_CD8_CD4.xls'),sep = '\t',col.names = T)


gg2 <- netVisual_bubble(cellchat, sources.use = c(2), targets.use = c(3), comparison = c(1,2),  angle.x = 45, remove.isolate = T,thresh = 0.01)
fwrite(data.frame(gg2$data),file = paste0('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-28/5.1/split1_CD4_CD8.xls'),sep = '\t',col.names = T)


kk=NULL
my_netVisual_bubble<-function (object, sources.use = NULL, targets.use = NULL, signaling = NULL,
          pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE,
          sort.by.source.priority = TRUE, color.heatmap = c("Spectral",
                                                            "viridis"), n.colors = 10, direction = -1, thresh = 0.05,
          comparison = NULL, group = NULL, remove.isolate = FALSE,
          max.dataset = NULL, min.dataset = NULL, min.quantile = 0,
          max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE,
          color.text = NULL, dot.size.min = NULL, dot.size.max = NULL,
          title.name = NULL, font.size = 10, font.size.title = 10,
          show.legend = TRUE, grid.on = TRUE, color.grid = "grey90",
          angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE)
{
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  }
  else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      pairLR.use$pathway_name <- as.character(pairLR.use$pathway_name)
    }
    else if ("interaction_name" %in% colnames(pairLR.use)) {
      pairLR.use$interaction_name <- as.character(pairLR.use$interaction_name)
    }
  }
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target,
                                  sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)),
                           targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate),
                                             ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target,
                              " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob <
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T) *
                             1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1),
                                               position)]
    }
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in%
                                                                  unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in%
                                                                  unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))),
                         levels(df.net$target), sep = " -> ")
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2),
    ])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2,
                                        levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target,
                                   levels = cells.order)
    df <- df.net
  }
  else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net",
                                      sources.use = sources.use, targets.use = targets.use,
                                      signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target,
                                    sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)),
                             targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate),
                                               ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target,
                                " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in%
                                                                    unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in%
                                                                    unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))),
                           levels(df.net$target), sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]],
                            ")")
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <=
                      0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      }
      else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names),
                                       ncol = 5))
        colnames(df.net) <- c("interaction_name_2", "source.target",
                              "prob", "pval", "prob.original")
        df.net$source.target <- group.names0
      }
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target,
                                     " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob <
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T) *
                             1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1),
                                               position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2,
                                    levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i],
                                             " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
  max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        dataset.na <- c(df.i.j$dataset[is.na(values)],
                        setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (all(!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset]))) {
            df.i.j$prob <- NA
          }
          else if (all((idx.max != idx.min) & !is.null(min.dataset))) {
            if (all(!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset]))) {
              df.i.j$prob <- NA
            }
            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in%
                                                    dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
  }
  # print(paste0("-->",df))

  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  if (!is.null(pairLR.use)) {
    interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name,
    ]$interaction_name_2, unique(df$interaction_name_2))
    df$interaction_name_2 <- factor(df$interaction_name_2,
                                    levels = interaction_name_2.order)
  }
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),
                                                                    unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use,
                                                        df$target))
      df <- with(df, df[order(target, source), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use,
                                                        df$source))
      df <- with(df, df[order(source, target), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use,
                                                        df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use,
                                                          df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target), ])
      }
      else {
        df <- with(df, df[order(target, source), ])
      }
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }


  df<-df%>%dplyr::filter(source !=target)
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2,
                      color = prob, size = pval)) + geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x,
                                     vjust = vjust.x), axis.title.x = element_blank(),
          axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
  values <- c(1, 2, 3)
  names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
  if (is.null(dot.size.max)) {
    dot.size.max = max(df$pval)
  }
  if (is.null(dot.size.min)) {
    dot.size.min = min(df$pval)
  }
  g <- g + scale_radius(range = c(dot.size.min, dot.size.max),
                        breaks = sort(unique(df$pval)), labels = names(values)[values %in%
                                                                                 sort(unique(df$pval))], name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
                                    na.value = "white", limits = c(quantile(df$prob,
                                                                            0, na.rm = T), quantile(df$prob, 1, na.rm = T)),
                                    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob,
                                                                                         1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                    title = "Commun. Prob."))
  }
  else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
                                    na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5,
                                                                                         title = "Commun. Prob."))
  }
  g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) -
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) -
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5 + length(dataset.name[comparison]),
                       length(group.names0) * length(dataset.name[comparison]),
                       by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype = "dashed",
                          color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      }
      else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order,
                                               "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order,
                                             2, stringr::str_length(dataset.name.order) -
                                               1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  }
  else {
    return(g)
  }
}


gg2 <- my_netVisual_bubble(cellchat, sources.use = c(2,3), targets.use = c(3,2), comparison = c(1,2), angle.x = 45, remove.isolate = T,thresh = 0.01,return.data = T)


fwrite(data.frame(gg2$communication),file = paste0('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-3-28/5.1/CD8_CD4_total.xls'),sep = '\t',col.names = T)




############################################### 2024-3-29
setwd('/2024-3-28/5.2')

tt_ntt_decon_celltype=carcinoma_regions_integrated1@meta.data%>%dplyr::select(17:25)%>%data.frame(check.rows = F,check.names = F)


options(digits = 7)
#
options(scipen = 0)
coor_for_decon<-function(tt_ntt_decon_celltype,samp){
 
  tmp<-tt_ntt_decon_celltype[,c("CD4._T_cells","CD8._T_cells")]

  tmp1<-tmp[apply(tmp,1,sum)>0,]
  tmp2<-subset(tmp1,!(CD4._T_cells==0 | CD8._T_cells ==0 ))

  g <- ggplot(tmp2, aes(`CD4._T_cells`,`CD8._T_cells`))


  df.lm<-lm(`CD4._T_cells`~
              `CD8._T_cells`,
            data=tmp2)
  rr=glance(df.lm)

  pvalue=rr$p.value
  pvalue_formatted <- format(pvalue, scientific = TRUE, digits = 2)
 
  rsquared=rr$r.squared


  p1<-g + geom_point(colour='red') +
    geom_smooth(method="lm", se=F,color='black') +
    labs(subtitle=paste0("correlation ",samp),
         y="CD4._T_cells",
         x="CD8._T_cells",
         title="",
         caption="")+
    theme_bw()+
    annotate(geom = "text",
             x=mean(tmp2$CD4._T_cells)+0.05, y=max(tmp2$CD8._T_cells),
             label = paste0("R=",round(rsquared,2),",","P < ",pvalue))

  p1

}



coor_for_decon4<-function(tt_ntt_decon_celltype,samp){
  
  tmp<-tt_ntt_decon_celltype[,c("CD4._T_cells","CD8._T_cells")]

  tmp1<-tmp[apply(tmp,1,sum)>0,]
  tmp2<-subset(tmp1,!(CD4._T_cells==0 | CD8._T_cells ==0 ))

  g <- ggplot(tmp2, aes(`CD4._T_cells`,`CD8._T_cells`))


  df.lm<-lm(`CD4._T_cells`~
              `CD8._T_cells`,
            data=tmp2)
  rr=glance(df.lm)
  pvalue=rr$p.value
  rsquared=rr$r.squared



  p1<-g + geom_point(colour='red') +
    geom_smooth(method="lm", se=F,color='black') +
    labs(subtitle=paste0("correlation ",samp),
         y="CD4._T_cells",
         x="CD8._T_cells",
         title="",
         caption="")+
    theme_bw()+
    annotate(geom = "text",
             x=mean(tmp2$CD4._T_cells)+0.1, y=max(tmp2$CD8._T_cells),
             label = paste0("R=",round(rsquared,2),",","P < ",pvalue))

  p1

}


summary(df.lm)$`p-value`
library(broom) # glance
rr1<-glance(df.lm) # 获取真是的p值




sample_plot=list()
tt_ntt_decon_celltype$sample<-sapply(strsplit(rownames(tt_ntt_decon_celltype),"_"), '[[',1)
for(samp in unique(tt_ntt_decon_celltype$sample)){
  temp_df<-subset(tt_ntt_decon_celltype,sample==samp)
  p1<-coor_for_decon(temp_df,samp)
  sample_plot[[samp]]=p1
}

sample_plot$`TT-1`
sample_plot$`TT-2`
sample_plot$`TT-3`
sample_plot$`TT-4`
sample_plot$`NTT-1`
sample_plot$`NTT-2`
sample_plot$`NTT-3`
sample_plot$`NTT-4`


for(samp in unique(tt_ntt_decon_celltype$sample)){
  temp_df<-subset(tt_ntt_decon_celltype,sample==samp)
  p1<-coor_for_decon4(temp_df,samp)
  sample_plot[[samp]]=p1
}


######################################################

# 构建CD8._T_cells和CD4._T_cells共spots基因集得分图（每个样本出图）

Idents(carcinoma_regions_integrated1)<-carcinoma_regions_integrated1$dominant_celltype_combined_0.1

sig_marker_spatial<-subset(FindAllMarkers(carcinoma_regions_integrated1,only.pos = T),p_val_adj <0.05)


sig_marker_spatial_subset<-subset(sig_marker_spatial,cluster%in%c("CD4._T_cells","CD8._T_cells"))

sig_marker_spatial_subset_top20<-sig_marker_spatial_subset%>%group_by(cluster)%>%top_n(n=20,wt=avg_log2FC)
sig_marker_spatial_subset_top10<-sig_marker_spatial_subset%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)


sig_marker_spatial_subset_top20_genes<-unique(sig_marker_spatial_subset_top20$gene)
sig_marker_spatial_subset_top10_genes<-unique(sig_marker_spatial_subset_top10$gene)

sig_marker_spatial_subset_top20_genes_list<-list(sig_marker_spatial_subset_top20_genes=sig_marker_spatial_subset_top20_genes)
sig_marker_spatial_subset_top10_genes_list<-list(sig_marker_spatial_subset_top10_genes=sig_marker_spatial_subset_top10_genes)
sig_marker_spatial_subset_cd4_cd8_genes_list<-list(cd4_cd8_geneset_score=c('CD4',"IL7R",'CD8A','CD8B'))

carcinoma_regions_integrated1<-AddModuleScore(carcinoma_regions_integrated1,features = sig_marker_spatial_subset_top10_genes_list,name = "sig_marker_spatial_subset_top10_genes_list_cd4_cd8")
carcinoma_regions_integrated1<-AddModuleScore(carcinoma_regions_integrated1,features = sig_marker_spatial_subset_cd4_cd8_genes_list,name = "sig_marker_spatial_subset_cd4_cd8_genes_list")


SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'sig_marker_spatial_subset_top10_genes_list_cd4_cd81',crop = F,ncol = 4)
SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'sig_marker_spatial_subset_cd4_cd8_genes_list1',crop = F,ncol = 4)


integrated_spatial$sig_marker_spatial_subset_top20_genes_list_cd4_cd8 <- mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$sig_marker_spatial_subset_top20_genes_list_cd4_cd81)


integrated_spatial$sig_marker_spatial_subset_top20_genes_list_cd4_cd8<-as.numeric(integrated_spatial$sig_marker_spatial_subset_top20_genes_list_cd4_cd8)


integrated_spatial$CD4_decon_score_check<-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$CD4._T_cells)
names(integrated_spatial@meta.data)[ncol(integrated_spatial@meta.data)]="CD4_decon_score_check"
integrated_spatial$CD8_decon_score_check<-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$CD8._T_cells)


integrated_spatial$CD4_decon_score_check<-as.numeric(integrated_spatial$CD4_decon_score_check)
integrated_spatial$CD8_decon_score_check<-as.numeric(integrated_spatial$CD8_decon_score_check)


integrated_spatial$CD4_CD8_score_decon<-integrated_spatial$CD4_decon_score_check + integrated_spatial$CD8_decon_score_check
integrated_spatial$CD4_CD8_score_decon_mean<-(integrated_spatial$CD4_decon_score_check + integrated_spatial$CD8_decon_score_check)/2


# theme_set(theme_minimal(base_size = 14, base_family = "Helvetica", na.rm = TRUE))


colnames(integrated_spatial@meta.data)[64]="CD8_and_CD4_top20_score"
p1<-SpatialFeaturePlot(integrated_spatial,features = 'CD8_and_CD4_top20_score',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
p1

tt1<-data.frame(cd4=carcinoma_regions_integrated1$CD4._T_cells,cd8=carcinoma_regions_integrated1$CD8._T_cells)


cor1<-cor(t(tt1))


pheatmap(cor1,cluster_rows = F,cluster_cols = F)
cor2<-cor(tt1)

tt <- corr.test(t(tt1),method="pearson",adjust="holm")

color_palette<-c("#00008B","#7B68EE","#FFFFFF","#FFD700","#FF0000")

# 显示相关性数字
corrplot(as.matrix(tt$r),col = colorRampPalette(color_palette)(100),
         method = "color",tl.col="black",
         tl.cex = 0.7,cl.pos = "r",cl.lim=c(-1,1),
         cl.ratio = 0.1,cl.cex=0.8,cl.length=5,p.mat = as.matrix(tt$p),
         sig.level = 1,outline="white",insig = "blank",type="full",addCoef.col = "black")


SpatialFeaturePlot(integrated_spatial,features = 'CD4',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
SpatialFeaturePlot(integrated_spatial,features = 'CD8A',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
SpatialFeaturePlot(integrated_spatial,features = 'CD8B',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")


# 另外一种方式，即CD4和CD8的decon分数加和
p1<-SpatialFeaturePlot(integrated_spatial,features = 'CD4_CD8_score_decon',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
p1
# 均值
# 另外一种方式，即CD4和CD8的decon分数加和
p1<-SpatialFeaturePlot(integrated_spatial,features = 'CD4_CD8_score_decon_mean',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
p1
#################################################################
#
# 5.4）对每个spots计算缺氧基因集得分
# 缺氧基因集：


geneset1<-list(HALLMARK_HYPOXIA=strsplit("ACKR3	ADM	ADORA2B	AK4	AKAP12	ALDOA	ALDOB	ALDOC	AMPD3	ANGPTL4	ANKZF1	ANXA2	ATF3	ATP7A	B3GALT6	B4GALNT2	BCAN	BCL2	BGN	BHLHE40	BNIP3L	BRS3	BTG1	CA12	CASP6	CAV1	CAVIN1	CAVIN3	CCN1	CCN2	CCN5	CCNG2	CDKN1A	CDKN1B	CDKN1C	CHST2	CHST3	CITED2	COL5A1	CP	CSRP2	CXCR4	DCN	DDIT3	DDIT4	DPYSL4	DTNA	DUSP1	EDN2	EFNA1	EFNA3	EGFR	ENO1	ENO2	ENO3	ERO1A	ERRFI1	ETS1	EXT1	F3	FAM162A	FBP1	FOS	FOSL2	FOXO3	GAA	GALK1	GAPDH	GAPDHS	GBE1	GCK	GCNT2	GLRX	GPC1	GPC3	GPC4	GPI	GRHPR	GYS1	HAS1	HDLBP	HEXA	HK1	HK2	HMOX1	HOXB9	HS3ST1	HSPA5	IDS	IER3	IGFBP1	IGFBP3	IL6	ILVBL	INHA	IRS2	ISG20	JMJD6	JUN	KDELR3	KDM3A	KIF5A	KLF6	KLF7	KLHL24	LALBA	LARGE1	LDHA	LDHC	LOX	LXN	MAFF	MAP3K1	MIF	MT1E	MT2A	MXI1	MYH9	NAGK	NCAN	NDRG1	NDST1	NDST2	NEDD4L	NFIL3	NOCT	NR3C1	P4HA1	P4HA2	PAM	PCK1	PDGFB	PDK1	PDK3	PFKFB3	PFKL	PFKP	PGAM2	PGF	PGK1	PGM1	PGM2	PHKG1	PIM1	PKLR	PKP1	PLAC8	PLAUR	PLIN2	PNRC1	PPARGC1A	PPFIA4	PPP1R15A	PPP1R3C	PRDX5	PRKCA	PYGM	RBPJ	RORA	RRAGD	S100A4	SAP30	SCARB1	SDC2	SDC3	SDC4	SELENBP1	SERPINE1	SIAH2	SLC25A1	SLC2A1	SLC2A3	SLC2A5	SLC37A4	SLC6A6	SRPX	STBD1	STC1	STC2	SULT2B1	TES	TGFB3	TGFBI	TGM2	TIPARP	TKTL1	TMEM45A	TNFAIP3	TPBG	TPD52	TPI1	TPST2	UGP2	VEGFA	VHL	VLDLR	WSB1	XPNPEP1	ZFP36	ZNF292","\t")[[1]])




carcinoma_regions_integrated1<-AddModuleScore(carcinoma_regions_integrated1,features = geneset1,name = "HALLMARK_HYPOXIA")


###############################################################
###############################################################



df<-carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,HALLMARK_HYPOXIA1,group)%>%data.frame(check.rows = F,check.names = F)

library(ggsignif)
df$group<-factor(df$group,levels=c('TT','NTT'))
ggplot(data = df,aes(x=group,y=HALLMARK_HYPOXIA1,fill=orig.ident))+
  geom_boxplot(position = position_dodge2(width = 0.8))+ylab(label = 'HALLMARK HYPOXIA Score')+
  guides(fill=guide_legend(title="samples"))+theme_classic()+theme(axis.title.x = element_blank())+
  geom_signif(comparisons = list(c('TT','NTT')),map_signif_level = F)



########################################################33
########################################################33
########################################################33
fwrite(df,file = "/HALLMARK_HYPOXIA_Score.xls",sep = '\t',col.names = T,row.names = F)



d1=data.frame()
tt=list(d=data.frame(a=c(1,2),b=c(3,4)),f=data.frame(a=c('a1','a2'),b=c('b1','b2')))
for(nm in names(tt)){
  d1 = rbind(d1,tt[[nm]])
}
d1


######################################################
######################################################
######################################################
###################################################### 2024-4-3
###################################################### 2024-4-3

setwd('/2024-4-3/5.1_adjust')


# carcinoma_regions_integrated1
# 方案 1
carcinoma_regions_integrated1$CD4_rank<-ifelse(carcinoma_regions_integrated1$CD4._T_cells<=0.02,0,
                                               ifelse(carcinoma_regions_integrated1$CD4._T_cells>0.02 & carcinoma_regions_integrated1$CD4._T_cells<=0.05,1,
                                                      ifelse(carcinoma_regions_integrated1$CD4._T_cells>0.05 & carcinoma_regions_integrated1$CD4._T_cells<=0.1,2,
                                                             ifelse(carcinoma_regions_integrated1$CD4._T_cells>0.1 & carcinoma_regions_integrated1$CD4._T_cells<=0.15,3,4))))




carcinoma_regions_integrated1$CD8_rank<-ifelse(carcinoma_regions_integrated1$CD8._T_cells<=0.02,0,
                                               ifelse(carcinoma_regions_integrated1$CD8._T_cells>0.02 & carcinoma_regions_integrated1$CD8._T_cells<=0.05,1,
                                                      ifelse(carcinoma_regions_integrated1$CD8._T_cells>0.05 & carcinoma_regions_integrated1$CD8._T_cells<=0.1,2,
                                                             ifelse(carcinoma_regions_integrated1$CD8._T_cells>0.1 & carcinoma_regions_integrated1$CD8._T_cells<=0.15,3,4))))





carcinoma_regions_integrated1$CD4_CD8_rank<-carcinoma_regions_integrated1$CD4_rank * carcinoma_regions_integrated1$CD8_rank


SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'CD4_CD8_rank',crop = F,ncol = 4)


# 方案 2
carcinoma_regions_integrated1$CD4_rank<-ifelse(carcinoma_regions_integrated1$CD4._T_cells<=0.05,0,
                                               ifelse(carcinoma_regions_integrated1$CD4._T_cells>0.05 & carcinoma_regions_integrated1$CD4._T_cells<=0.1,1,
                                                      ifelse(carcinoma_regions_integrated1$CD4._T_cells>0.1 & carcinoma_regions_integrated1$CD4._T_cells<=0.15,2,3
                                                             )))


carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,CD4_rank)%>%group_by(CD4_rank)%>%table()
#
carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,CD8_rank)%>%group_by(CD8_rank)%>%table()
#
carcinoma_regions_integrated1@meta.data%>%dplyr::select(orig.ident,CD4_CD8_rank)%>%group_by(CD4_CD8_rank)%>%table()


carcinoma_regions_integrated1$CD8_rank<-ifelse(carcinoma_regions_integrated1$CD8._T_cells<=0.05,0,
                                               ifelse(carcinoma_regions_integrated1$CD8._T_cells>0.05 & carcinoma_regions_integrated1$CD8._T_cells<=0.1,1,
                                                      ifelse(carcinoma_regions_integrated1$CD8._T_cells>0.1 & carcinoma_regions_integrated1$CD8._T_cells<=0.15,2,3
                                                            )))





carcinoma_regions_integrated1$CD4_CD8_rank<-carcinoma_regions_integrated1$CD4_rank * carcinoma_regions_integrated1$CD8_rank


SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'CD4_CD8_rank',crop = F,ncol = 4)



setwd('/2024-4-3/5.3')
malignant_tt_vs_ntt<-fread('../data/Malignant_epithelial_cells-TT_vs_NTT_pval0.05.xls',header = T)


malignant_tt_vs_ntt<-malignant_tt_vs_ntt[order(malignant_tt_vs_ntt$avg_log2FC,decreasing = T),]


combined_polar<-rbind(top10_df,tail10_df)

combined_polar$group<-ifelse(combined_polar$avg_log2FC>0,'up','down')
combined_polar$group<-factor(combined_polar$group,levels=c('up','down'))
combined_polar$V1<-factor(combined_polar$V1,levels=unique(combined_polar$V1))
ggplot(data = combined_polar,aes(x=V1,y=avg_log2FC,fill=group))+geom_bar(stat='identity')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),axis.title.x=element_blank()) +
  ylab(label = "average log2 Fold-Change") +
  scale_fill_manual(values = c(up="#AB3282",down="#58A4C3"))+
  guides(fill=guide_legend(title="Regulated"))

####################################TF


tf_df<-fread('/01.TT_vs_NTT/diff_tt_vs_ntt.xls',header = T)



tf_df<-tf_df[order(tf_df$avg_log2FC,decreasing = T),]
tf_df_top10 <- head(tf_df,10)
tf_df_tail10 <- tail(tf_df,10)

combined_tf_df<-rbind(tf_df_top10,tf_df_tail10)
combined_tf_df$V1 <-factor(combined_tf_df$V1,levels=unique(combined_tf_df$V1))
combined_tf_df$group<-ifelse(combined_tf_df$avg_log2FC>0,'up','down')

combined_tf_df$group<-factor(combined_tf_df$group,levels=c('up','down'))

ggplot(data = combined_tf_df,aes(x=V1,y=avg_log2FC,fill=group))+geom_bar(stat='identity')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),axis.title.x=element_blank()) +
  ylab(label = "average log2 Fold-Change") +
  scale_fill_manual(values = c(up="#AB3282",down="#58A4C3"))+
  guides(fill=guide_legend(title="Regulated"))



##################################################Treg基因染色
##################################################Treg基因染色
##################################################Treg基因染色

# IL2RA
SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'IL2RA',crop = F,ncol = 4)
# # TGFB1
# SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'TGFB1',crop = F,ncol = 4)

SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'ICOS',crop = F,ncol = 4)
SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'CCR4',crop = F,ncol = 4)
SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'IL7R',crop = F,ncol = 4)
SpatialFeaturePlot(carcinoma_regions_integrated1,features = 'FOXP3',crop = F,ncol = 4)

#####################################im()

integrated_spatial$IL2RA_expr <-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1@assays$SCT@data["IL2RA",])
integrated_spatial$ICOS_expr <-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1@assays$SCT@data["ICOS",])
integrated_spatial$CCR4_expr <-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1@assays$SCT@data["CCR4",])
integrated_spatial$IL7R_expr <-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1@assays$SCT@data["IL7R",])
integrated_spatial$FOXP3_expr <-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1@assays$SCT@data["FOXP3",])
View(integrated_spatial@meta.data)
integrated_spatial$IL2RA_expr<-as.numeric(integrated_spatial$IL2RA_expr)
integrated_spatial$ICOS_expr<-as.numeric(integrated_spatial$ICOS_expr)
integrated_spatial$CCR4_expr<-as.numeric(integrated_spatial$CCR4_expr)
integrated_spatial$IL7R_expr<-as.numeric(integrated_spatial$IL7R_expr)
integrated_spatial$FOXP3_expr<-as.numeric(integrated_spatial$FOXP3_expr)


setwd('/2024-4-3/5.6')
for(i in c("IL2RA_expr","ICOS_expr","CCR4_expr","IL7R_expr","FOXP3_expr")){
  p1<-SpatialFeaturePlot(integrated_spatial,features = i,crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")
  ggsave(filename = paste0(i,"_spatial.pdf"),width = 45,height = 45)
}


#############################################################
#############################################################
#############################################################

setwd('/2024-4-7/5.6')
treg_genelist<-c("CCR4","FOXP3","ICOS","IL2RA","IL7R")

treg_genelist=list(treg=treg_genelist)

carcinoma_regions_integrated1<-AddModuleScore(carcinoma_regions_integrated1,treg_genelist,name="treg_score")


########################################################
cols=c("#43387F","#288D88","#8AD24C","#FBE625")
SpatialFeaturePlot(carcinoma_regions_integrated1,features = "treg_score1",crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = cols,na.value = "lightgrey")

integrated_spatial$treg_score1<-mapvalues(colnames(integrated_spatial),colnames(carcinoma_regions_integrated1),carcinoma_regions_integrated1$treg_score1)
integrated_spatial$treg_score1 <- as.numeric(integrated_spatial$treg_score1)
SpatialFeaturePlot(integrated_spatial,features = "treg_score1",crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = cols,na.value = "lightgrey")


########################################################
########################################################
########################################################

DimPlot(integrated_spatial,label=T,label.size = 10)+theme(legend.text = element_text(size=20))

icms_geneset<-fread('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-4-7/5.5/icms_geneset.xls',header = T)


icms_geneset<-as.list(icms_geneset)
###################################

icms_geneset1<-lapply(icms_geneset,function(x){x[x!=""]})


###################################

Malignant_epithelial_cells_cms<-readRDS(file = "/data/Malignant_epithelial_cells.rds")

icms_geneset2<-lapply(icms_geneset1,function(x){x[1:50]})

Malignant_epithelial_cells_cms<-AddModuleScore(Malignant_epithelial_cells_cms,features = icms_geneset2,nbin = 50)


Malignant_epithelial_cells_cms$iCMS2_Up<-NULL
Malignant_epithelial_cells_cms$iCMS2_Down<-NULL
Malignant_epithelial_cells_cms$iCMS3_Up<-NULL
Malignant_epithelial_cells_cms$iCMS3_Down<-NULL

names(Malignant_epithelial_cells_cms@meta.data)[names(Malignant_epithelial_cells_cms@meta.data)=="Cluster1"]<-"iCMS2_Up"
names(Malignant_epithelial_cells_cms@meta.data)[names(Malignant_epithelial_cells_cms@meta.data)=="Cluster2"]<-"iCMS2_Down"
names(Malignant_epithelial_cells_cms@meta.data)[names(Malignant_epithelial_cells_cms@meta.data)=="Cluster3"]<-"iCMS3_Up"
names(Malignant_epithelial_cells_cms@meta.data)[names(Malignant_epithelial_cells_cms@meta.data)=="Cluster4"]<-"iCMS3_Down"

FeaturePlot(Malignant_epithelial_cells_cms,features = "Cluster1",label = T,order = T)
FeaturePlot(Malignant_epithelial_cells_cms,features = "Cluster2",label = T,order = T)
FeaturePlot(Malignant_epithelial_cells_cms,features = "Cluster3",label = T,order = T)
FeaturePlot(Malignant_epithelial_cells_cms,features = "Cluster4",label = T,order = T)

FeaturePlot(Malignant_epithelial_cells_cms,features = "iCMS2_Up",label = T,order = F,label.size = 10)
FeaturePlot(Malignant_epithelial_cells_cms,features = "iCMS2_Down",label = T,order = T,label.size = 10)
FeaturePlot(Malignant_epithelial_cells_cms,features = "iCMS3_Up",label = T,order = T,label.size = 10)
FeaturePlot(Malignant_epithelial_cells_cms,features = "iCMS3_Down",label = T,order = T,label.size = 10)




DefaultAssay(Malignant_epithelial_cells_cms)<-'integrated'

Malignant_epithelial_cells_cms<-FindClusters(Malignant_epithelial_cells_cms,resolution = 0.4)
colnames(Malignant_epithelial_cells_cms@meta.data)
DimPlot(Malignant_epithelial_cells_cms,label = T)
DimPlot(Malignant_epithelial_cells_cms,label = T,group.by = 'integrated_snn_res.0.2')+ggtitle("")
DefaultAssay(Malignant_epithelial_cells_cms)<-'SCT'
sig_marker<-subset(FindAllMarkers(Malignant_epithelial_cells_cms,only.pos = T),p_val_adj <0.05)



fwrite(sig_marker,file = '/5.5/sig_marker.xls',sep = '\t',col.names = T,row.names = F)


###########################################################
###########################################################
###########################################################

wb<-xlsx::loadWorkbook('/share/home/hhhong/personality/2023-8-8_jingdi_wuhanxiangya/2024-4-7/5.5/41588_2022_1100_MOESM3_ESMprocess.xlsx')

sheet_name<-xlsx::getSheets(wb)


geneset1=list()
for(i in names(sheet_name)){
  geneset<-xlsx::read.xlsx('/41588_2022_1100_MOESM3_ESMprocess.xlsx',sheetName = i)
  geneset1[[i]]=geneset
}

View(geneset1$`Supplementary Table 2`)


geneset2<-lapply(geneset1,function(x){
  reshape2::melt(x,id.vars = NULL, variable.name = "term", value.name = "gene")
})


geneset2<-lapply(geneset2,function(x){
  x[!is.na(x$gene),]
})



sig_marker<-subset(FindAllMarkers(Malignant_epithelial_cells_cms),p_val <0.05)
sig_marker1<-subset(FindAllMarkers(Malignant_epithelial_cells_cms,only.pos = T),p_val <0.05)

temp <- subset(sig_marker,cluster==0)
temp <- subset(sig_marker1,cluster==0)
geneList=temp$avg_log2FC
names(geneList)=temp$gene
geneList=sort(geneList,decreasing = T)
geneset_temp<-geneset2$`Supplementary Table 2`

gsea.hsa = GSEA(geneList = geneList, TERM2GENE=geneset_temp,nPerm = 100000, pvalueCutoff = 1)

fwrite(sig_marker,file = '~/personality/2023-8-8_jingdi_wuhanxiangya/2024-4-7/5.5/diff_marker.xls',sep = '\t',col.names = T,row.names = F)



#####################################################3
#####################################################3
#####################################################3
# 新分的这个cluster0,1,6划分为iCMS2, cluster2,3,4,5划分为iCMS3

DimPlot(Malignant_epithelial_cells_cms,label = T)

cluster1<-c(0,1,6,2,3,4,5)
type1<-c("iCMS2","iCMS2","iCMS2","iCMS3","iCMS3","iCMS3","iCMS3")

Malignant_epithelial_cells_cms$icms_type<-mapvalues(Malignant_epithelial_cells_cms$seurat_clusters,cluster1,type1)

DimPlot(Malignant_epithelial_cells_cms,group.by = 'icms_type',label = T)

Idents(Malignant_epithelial_cells_cms)<-Malignant_epithelial_cells_cms$icms_type

icms2_vs_icms3<-subset(FindMarkers(Malignant_epithelial_cells_cms,ident.1 = "iCMS2",ident.2 = 'iCMS3'),p_val<0.05)
setwd('/2024-4-11/5.5')
fwrite(icms2_vs_icms3,file = 'icms2_vs_icms3.xls',sep = '\t',col.names = T,row.names = T)
########################################################

########################################################2024-4-12


# 5.1）


setwd('/2024-4-12/5.1')

DimPlot(Malignant_epithelial_cells_cms,cols = c("#3C5BA8","#B797C6"),label = T)

DimPlot(Malignant_epithelial_cells_cms,pt.size = 2,cols = c("#3C5BA8","#B797C6"),label = T)


# 5.2）
# i2和i3在8个样本中的空间分布
# 颜色选择：i2:#3C5BA8   i3:#B797C6
setwd('/2024-4-12/5.2')


DimPlot(integrated_spatial,label=T,label.size = 10)+theme(legend.text = element_text(size=20))
integrated_spatial$iCMS_type<-mapvalues(colnames(integrated_spatial),colnames(Malignant_epithelial_cells_cms),as.character(Malignant_epithelial_cells_cms$icms_type))


integrated_spatial$iCMS_type[grepl("-1",integrated_spatial$iCMS_type)]=NA

#####################################################
#####################################################

Idents(integrated_spatial)<-integrated_spatial$iCMS_type
color_icms1<-c("#3C5BA8","#B797C6")
names(color_icms1)<-levels(integrated_spatial)
SpatialDimPlot(integrated_spatial,ncol = 4,crop = F,cols = color_icms1) & scale_fill_manual(values = color_icms1,breaks = c("iCMS2","iCMS3"),na.value = 'lightgrey')


setwd('/2024-4-12/5.3')

# 5.3）
# 柱状图i2和i3在各样本中的占比
# 颜色选择：i2:#3C5BA8   i3:#B797C6


stat1<-Malignant_epithelial_cells_cms@meta.data%>%dplyr::select(group,icms_type)%>%
  table()%>%prop.table(margin = 1)%>%data.frame(check.rows = F,check.names = F)

stat1$group<-factor(stat1$group,levels = c('TT','NTT'))

stat1$Freq<-round(stat1$Freq*100,2)

ggplot(data = stat1,aes(x=group,y=Freq,fill=icms_type))+geom_bar(stat="identity",position = 'stack')+
  geom_text(aes(label=Freq), color="white", size=3.5,position=position_stack(0.5))+
  scale_fill_manual(values = c(iCMS2='#3C5BA8',iCMS3='#B797C6'))+
  ylab(label = "Percentage(%)")+
  guides(fill=guide_legend(title="iCMS type")) +
  theme_classic()+theme(axis.title.x = element_blank())


stat2<-Malignant_epithelial_cells_cms@meta.data%>%dplyr::select(orig.ident,icms_type)%>%
  table()%>%prop.table(margin = 1)%>%data.frame(check.rows = F,check.names = F)

stat2$orig.ident<-factor(stat2$orig.ident,levels = c(paste0('TT.',1:4),paste0('NTT.',1:4)))
stat2$Freq<-round(stat2$Freq*100,2)


ggplot(data = stat2,aes(x=orig.ident,y=Freq,fill=icms_type))+geom_bar(stat="identity",position = 'stack')+
  geom_text(aes(label=Freq), color="white", size=3.5,position=position_stack(0.5))+
  scale_fill_manual(values = c(iCMS2='#3C5BA8',iCMS3='#B797C6'))+
  ylab(label = "Percentage(%)")+
  guides(fill=guide_legend(title="iCMS type")) +
  theme_classic()+theme(axis.title.x = element_blank())

#####################################################
#####################################################
#####################################################
#####################################################

# 5.4）
# i2和i3的差异基因和差异通路的做热图
# 展示的差异基因：（一张图）
# 从上到下按下面的顺序展示
# BLCAP, CCL20, CYP2S1, DYNLRB1, EDN1, TSPAN6, COL16A1, DPYSL2, IFI16, IFI6, IGHG4, TYMP

Malignant_epithelial_cells_cms
geneset2<-strsplit('BLCAP, CCL20, CYP2S1, DYNLRB1, EDN1, TSPAN6, COL16A1, DPYSL2, IFI16, IFI6, IGHG4, TYMP',", " )[[1]]
expr<-AverageExpression(Malignant_epithelial_cells_cms)
expr<-expr$SCT[geneset2,]
expr1<-t(scale(t(expr),center = F,scale = T))
pheatmap::pheatmap(expr,scale = 'column',cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(expr,scale = 'row',cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(expr1,scale = 'none',cluster_rows = F,cluster_cols = F)



geneset2<-fread('~/personality/2023-8-8_jingdi_wuhanxiangya/2024-4-12/5.4/pathwaty1.xls',header = T)
geneset2<-geneset2%>%data.frame(check.names = F,check.rows = F)


geneset3<-as.list(geneset2)


geneset3<-lapply(geneset3,function(x){x[x!=""]})


Malignant_epithelial_cells_cms<-AddModuleScore(Malignant_epithelial_cells_cms,features = geneset3)


score_df<-Malignant_epithelial_cells_cms@meta.data%>%dplyr::select(c(69:90))%>%data.frame(check.names = F,check.rows = F)

score_df$icms_type<-Malignant_epithelial_cells_cms$icms_type



# 打分通路
score_df_stat1<-score_df%>%dplyr::group_by(icms_type)%>%summarise(across(everything(),mean,na.rm=TRUE))%>%data.frame(check.rows = F,check.names = F)

rownames(score_df_stat1)<-score_df_stat1$icms_type
score_df_stat1$icms_type<-NULL
colnames(score_df_stat1)<-names(geneset3)

score_df_stat1<-score_df_stat1%>%t()%>%data.frame(check.rows = F,check.names = F)


# rownames(score_df_stat1)<-gsub("^BP_","GO_BP_",rownames(score_df_stat1))
scale_data<-t(scale(t(score_df_stat1),center = F,scale=T))
pheatmap::pheatmap(scale_data,cluster_rows = F,cluster_cols = F)



setwd('/2024-4-12/5.5')

# 5.5）
# 不同亚型标志基因在降维图中表达
# BLCAP, CCL20, CYP2S1, DYNLRB1, EDN1, TSPAN6, COL16A1, DPYSL2, IFI16, IFI6, IGHG4, TYMP


geneset1<-strsplit("BLCAP, CCL20, CYP2S1, DYNLRB1, EDN1, TSPAN6, COL16A1, DPYSL2, IFI16, IFI6, IGHG4, TYMP",", ")[[1]]

p4 <- plot_density(Malignant_epithelial_cells_cms, reduction ='umap',geneset1, joint = TRUE,combine = T)
p4
# p4 + plot_layout(ncol = 4)

setwd('/2024-4-12/5.5')
for(i in geneset1){
  p1 <- plot_density(Malignant_epithelial_cells_cms, size = 2,reduction ='umap',i, joint = F)
  ggsave(p1,filename = paste0(i,"_density.pdf"),width = 12,height = 6)
}

options(scipen = -1)
options(digits = 5)
p4 <- plot_density(Malignant_epithelial_cells, genes, joint = F,size = 2,reduction = 'umap')
p4 + plot_layout(ncol = 4)
Idents(Malignant_epithelial_cells)<-Malignant_epithelial_cells$cms_type
DimPlot(Malignant_epithelial_cells,pt.size = 2,cols = color_cms1)

p4 <- plot_density(Malignant_epithelial_cells, genes, joint = F,size = 0.2,reduction = 'umap')
p4 + plot_layout(ncol = 4)

#############################################################################
#############################################################################
#############################################################################



#############################################################################
#############################################################################



# 免疫检查点基因PDCD1空间表达


SpatialFeaturePlot(integrated_spatial,features = 'PDCD1',crop = F,alpha = 1,interactive=F,ncol = 4) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],na.value = "lightgrey")

colnames(integrated_spatial@meta.data)
colnames(carcinoma_regions_integrated1@meta.data)
table(is.na(integrated_spatial$IL7R_expr))

integrated_spatial$PDCD1_expr <-integrated_spatial@assays$SCT@data['PDCD1',]
integrated_spatial$PDCD1_expr[is.na(integrated_spatial$IL7R_expr)]=NA
# CD274
integrated_spatial$CD274_expr <-integrated_spatial@assays$SCT@data['CD274',]
integrated_spatial$CD274_expr[is.na(integrated_spatial$IL7R_expr)]=NA

# ARG1
integrated_spatial$ARG1_expr <-integrated_spatial@assays$SCT@data['ARG1',]
integrated_spatial$ARG1_expr[is.na(integrated_spatial$IL7R_expr)]=NA
?subset
library(SeuratObject)
Idents(integrated_spatial)<-integrated_spatial$orig.ident
ntt.1<-subset(integrated_spatial,cells=rownames(integrated_spatial@meta.data)[integrated_spatial$orig.ident=='NTT.1'])
ntt.1<-subset(integrated_spatial,ident='NTT.1')
integrated_spatial$orig.ident<-as.character(integrated_spatial$orig.ident)
ntt.1<-subset(integrated_spatial,orig.ident=='NTT.1')






SpatialDimPlot(ntt.1)
temp_img<-ntt.1@images
aa<-temp_img$slice1_NTT.1


ntt.1@images$slice1_TT.1<-NULL
ntt.1@images$slice1_TT.2<-NULL
ntt.1@images$slice1_TT.3<-NULL
ntt.1@images$slice1_TT.4<-NULL
ntt.1@images$slice1_NTT.2<-NULL
ntt.1@images$slice1_NTT.3<-NULL
ntt.1@images$slice1_NTT.4<-NULL

SpatialFeaturePlot(ntt.1,features = 'PDCD1_expr',crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],limits=c(0,1.6),na.value = "lightgrey",oob=scales::squish)
SpatialFeaturePlot(ntt.1,features = 'CD274_expr',crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],limits=c(0,1.6),na.value = "lightgrey",oob=scales::squish)
SpatialFeaturePlot(ntt.1,features = 'ARG1_expr',crop = F,alpha = 1,interactive=F) & scale_fill_gradientn(colors = FeaturePalettes[['Spatial']],limits=c(0,1.0),na.value = "lightgrey",oob=scales::squish)





