
#---------------------------------- parameters
outputDir = paste0('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/DCOA/')
dir.create(outputDir)

inputDir = '/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis/'

pruning_method2 = TRUE
pruning_method1 = FALSE

#---------------------------------- Read Count
library(data.table)
library(dplyr)

setwd(inputDir)

# profiles = c('GSE184950', 'GSE243639', 'Inhouse')
profiles = c('GSE243639')

for(p in 1:length(profiles))
{
  temp = fread(paste0(profiles[p], '/', profiles[p], '_Neurons_Subset_MetaData.txt'), stringsAsFactors = F, header = T)
  temp = as.data.frame(temp)
  temp$DataSet = profiles[p]
  if(p == 1){MetaData = temp}else{MetaData = rbind(MetaData, temp)}
}

table(MetaData$yan)

p=2

for(p in 1:length(profiles))
{
  temp = fread(paste0(profiles[p], '/', profiles[p], '_Neurons_Subset_GeneCount.txt'), stringsAsFactors = F, header = T)
  temp = as.data.frame(temp)
  if(p == 1){GeneCount = temp}else{GeneCount = merge(GeneCount, temp, by = 'GeneSymbol')}
}

row.names(GeneCount) = GeneCount$GeneSymbol
GeneCount = GeneCount[,-1]

GeneCount = GeneCount[apply(GeneCount, 1, function(x) !all(x==0)),]
GeneCount = GeneCount[apply(GeneCount, 2, function(x) !all(x==0)),]

colnames(GeneCount) = gsub('\\.', '-', colnames(GeneCount))

colnames(GeneCount)[grep('ACGTAACCAGAGTGAC', colnames(GeneCount))]
MetaData$Barcode[grep('ACGTAACCAGAGTGAC', MetaData$Barcode)]

identical(MetaData$Barcode, colnames(GeneCount))

#---------------------------------- Seurat analysis
library(Seurat)
target_Data_Seurat <- CreateSeuratObject(counts = GeneCount, project = "TS_SC_Integration",  meta.data = MetaData, assay = "TS_SC")

#---------------------------------- pairwise analysis (select donors which present V and R both)
Selected_Donors = intersect(unique(MetaData$Donor[which(MetaData$yan == "Resistant" )]), unique(MetaData$Donor[which(MetaData$yan == "Vulnerable" )]))
MetaData_Selected = MetaData[which(MetaData$Donor %in% Selected_Donors),]

target_Data_Seurat <- subset(target_Data_Seurat, cells = MetaData_Selected$Barcode)
target_Data_Seurat

#---------------------------------- visualization as tSNE
target_Data_Seurat <- NormalizeData(target_Data_Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
target_Data_Seurat <- FindVariableFeatures(target_Data_Seurat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(target_Data_Seurat)
target_Data_Seurat <- ScaleData(target_Data_Seurat, features = all.genes)
target_Data_Seurat <- RunPCA(target_Data_Seurat, features = VariableFeatures(object = target_Data_Seurat))

target_Data_Seurat <- FindNeighbors(target_Data_Seurat, dims = 1:10)
target_Data_Seurat <- FindClusters(target_Data_Seurat, resolution = 0.02)  # I changed the resolution to get 6 clusters
table(target_Data_Seurat$seurat_clusters)

target_Data_Seurat <- RunUMAP(target_Data_Seurat, dims = 1:10)
#target_Data_Seurat <- RunTSNE(target_Data_Seurat, dims = 1:10, check_duplicates = FALSE)

#---- renames clusters
# Save old identity classes (the cluster labels) for reference.
# target_Data_Seurat[["old.ident"]] <- Idents(object = target_Data_Seurat)

identical(names(Idents(object = target_Data_Seurat)), target_Data_Seurat@meta.data$Barcode)

new.cluster.ids = target_Data_Seurat@meta.data$yan
names(new.cluster.ids) = names(Idents(object = target_Data_Seurat))
Idents(object = target_Data_Seurat) = factor(new.cluster.ids, levels = c("Resistant", "Vulnerable"))
 
DimPlot(target_Data_Seurat, reduction = "umap", label = TRUE, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))

#---- clean UMAP plot (remove intermediate cells)


#---- method-1
if(pruning_method1)
{
  predicted_clusters = target_Data_Seurat@reductions$umap@cell.embeddings
  predicted_clusters = as.data.frame(predicted_clusters)
  
  library(flexmix)
  library(mvtnorm)
  mixed_model = flexmix(cbind(as.numeric(predicted_clusters[,1]),as.numeric(predicted_clusters[,2]))~1,
                        k=2, data=predicted_clusters,
                        model = FLXMCmvnorm(diag = T), #diag = T to get fix results FLXMCmvpois(), FLXMCmvnorm
                        control = list(tolerance = 1e-2, iter.max = 1000))
  
  #--------------------------- assignments comparison
  predicted_clusters$predicted_clusters = clusters(mixed_model)
  real_clusters = target_Data_Seurat@meta.data
  
  table(real_clusters$yan)
  table( clusters(mixed_model) )
  
  identical(row.names(predicted_clusters), row.names(real_clusters))
  
  merged = merge(predicted_clusters, real_clusters[,c('yan', 'Donor')], by = 'row.names')
  
  freq_c1 = data.frame(table(merged[which(merged$predicted_clusters == 1),"yan"]), stringsAsFactors = F)
  freq_c1 = freq_c1[order(freq_c1$Freq, decreasing = T),]
  freq_c1
  
  merged[which(merged$predicted_clusters == 1), "predicted_clusters"] = as.character(freq_c1$Var1[1])
  
  freq_c2 = data.frame(table(merged[which(merged$predicted_clusters == 2), "yan"]), stringsAsFactors = F)
  freq_c2 = freq_c2[order(freq_c2$Freq, decreasing = T),]
  freq_c2
  
  merged[which(merged$predicted_clusters == 2), "predicted_clusters"] = as.character(freq_c2$Var1[1])
  
  table(merged$predicted_clusters == merged$yan)
  
  merged_subset = merged[which(merged$predicted_clusters == merged$yan),]
  table(merged_subset$yan)
  
  # if I have two distinct clusters, otherwise I keep all the cells as they are
  if(length(unique(merged_subset$yan)) > 1)
  {
  
    #--------------------------- keep cells which generate distinct clusters
    target_Data_Seurat = subset(target_Data_Seurat, cells = merged_subset$Row.names)
    DimPlot(target_Data_Seurat, group.by = 'yan', reduction = "umap", label = TRUE, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))
  
    #---------------------------------- pairwise analysis (select donors which present H and L both)
    Donors = data.frame(table(gsub('H|L','',unique(target_Data_Seurat@meta.data$Donor_subtype))))
    Donors_dup = Donors[which(Donors$Freq > 1),]
  
    MetaData_Selected = MetaData[which(MetaData$Donor %in% as.character(Donors_dup$Var1)),]
  
    target_Data_Seurat <- subset(target_Data_Seurat, cells = MetaData_Selected$Barcode)
    target_Data_Seurat
  
    table(target_Data_Seurat@meta.data$yan)
  }
}
# 
# #---- visualization
# # Set color palette
# library(viridis)
# pal <- viridis(n = 10, option = "D")
# 
# palette.colors(palette = "R4")
# 
# setwd(outputDir)
# # label off
# pdf(paste0(outputDir, 'UMAP_', C, '_', S, '_SC.pdf'), height = 8, width = 10, useDingbats = F)
# print(DimPlot(target_Data_Seurat, group.by = 'yan', reduction = "umap", label = F, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c("#F5C710", "#9E9E9E")))
# dev.off() 
# 
# library(scCustomize)
# pdf(paste0(outputDir, 'UMAP_', C, '_', S, '_SC_scCustom.pdf'), height = 8, width = 10, useDingbats = F)
# print(DimPlot_scCustom(seurat_object = target_Data_Seurat, figure_plot = TRUE, repel = T, pt.size = 2, label.size = 8, label.box = T, colors_use = c("#B4DE2CFF", "#FDE725FF")))
# dev.off() 

#---- method-2
if(pruning_method2)
{
  umap = target_Data_Seurat@reductions$umap@cell.embeddings
  umap = as.data.frame(umap)
  
  umap_V = (umap[which(row.names(umap) %in% names(new.cluster.ids[which(new.cluster.ids == "Vulnerable" )]) ),])
  umap_V_subset = umap_V[which(umap_V$umap_1 > 0),]
  
  umap_R = (umap[which(row.names(umap) %in% names(new.cluster.ids[which(new.cluster.ids == "Resistant" )]) ),])
  umap_R_subset = umap_R[which(umap_R$umap_1 < 0),]
  
  Selected_Cells = c(row.names(umap_R_subset), row.names(umap_V_subset))
  
  target_Data_Seurat = subset(target_Data_Seurat, cells = Selected_Cells)
  
  DimPlot(target_Data_Seurat, reduction = "umap", label = TRUE, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))
  table(target_Data_Seurat@meta.data$yan)

  length(unique(target_Data_Seurat@meta.data$Donor_subtype))
}

#------
# library(scCustomize)
# pdf(paste0(outputDir, 'UMAP_VR_Neurons_SC_nonPairWise_scCustom.pdf'), height = 5, width = 6, useDingbats = F)
# DimPlot_scCustom(seurat_object = target_Data_Seurat, figure_plot = TRUE, repel = T, pt.size = 2, label.size = 8, label.box = T, colors_use = c("#B4DE2CFF", "#FDE725FF"))
# dev.off()

#---------------------------------- pairwise analysis (select donors which present V and R both)

Donors = data.frame(table(gsub('R|V','',unique(target_Data_Seurat@meta.data$Donor_subtype))))
Donors_dup = Donors[which(Donors$Freq > 1),]

MetaData_Selected = MetaData[which(MetaData$Donor %in% as.character(Donors_dup$Var1)),]

target_Data_Seurat <- subset(target_Data_Seurat, cells = MetaData_Selected$Barcode)
target_Data_Seurat

table(target_Data_Seurat@meta.data$yan)

#---- visualization
# setwd(outputDir)
# pdf(paste0(outputDir, 'UMAP_VR_Neurons_SC.pdf'), height = 8, width = 10, useDingbats = F)
# DimPlot(target_Data_Seurat, reduction = "umap", label = TRUE, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))
# dev.off() 
# 
# # Set color palette
# library(viridis)
# pal <- viridis(n = 10, option = "D")
# 
# library(scCustomize)
# pdf(paste0(outputDir, 'UMAP_VR_Neurons_SC_scCustom.pdf'), height = 8, width = 10, useDingbats = F)
# DimPlot_scCustom(seurat_object = target_Data_Seurat, figure_plot = TRUE, repel = T, pt.size = 2, label.size = 8, label.box = T, colors_use = c("#B4DE2CFF", "#FDE725FF"))
# dev.off() 

dim(target_Data_Seurat@meta.data)
dim(GeneCount)

GeneCount_subset = GeneCount[,which(colnames(GeneCount) %in% target_Data_Seurat@meta.data$Barcode)]
dim(GeneCount_subset)

#---------------

setwd(outputDir)

#---------------
library(data.table)
library(SAVER)
library(Seurat)
#========================================== transform assigned cells per selected pairs of individuals

#========= subset gene expression profile
GeneCount_subset = GeneCount_subset[apply(GeneCount_subset, 1, function(x) !all(x==0)),] # I did this in the QC step
dim(GeneCount_subset)

#========= run saver for assigned cells per selected pairs of individuals
set.seed(123)# if you do not use seed, you get different values per run
estimate_S1S2 <- saver(GeneCount_subset, ncores = 6, estimates.only = TRUE)
estimate_S1S2 = data.frame(GeneID=row.names(estimate_S1S2), estimate_S1S2)
estimate_S1S2[1:5,1:5]
#---------------

dir.create(paste0(outputDir, '/SAVER/'))
setwd(paste0(outputDir, '/SAVER/'))

fwrite(estimate_S1S2, 'SAVER.txt', row.names = F, quote = F, sep = '\t')
fwrite(target_Data_Seurat@meta.data, 'Metadata_SAVER.txt', row.names = F, quote = F, sep = '\t')

