
#---------------------------------- parameters
outputDir = paste0('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/')
inputDir = '/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis/'
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
length(unique(MetaData$Donor))

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
target_Data_Seurat <- RunTSNE(target_Data_Seurat, dims = 1:10, check_duplicates = FALSE)

#---- renames clusters
# Save old identity classes (the cluster labels) for reference.
# target_Data_Seurat[["old.ident"]] <- Idents(object = target_Data_Seurat)

identical(names(Idents(object = target_Data_Seurat)), target_Data_Seurat@meta.data$Barcode)

new.cluster.ids = target_Data_Seurat@meta.data$yan
names(new.cluster.ids) = names(Idents(object = target_Data_Seurat))
Idents(object = target_Data_Seurat) = factor(new.cluster.ids, levels = c("Resistant", "Vulnerable"))
 
DimPlot(target_Data_Seurat, reduction = "umap", label = F, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))

#-------------------------- clean UMAP plot (remove intermediate cells)

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

identical(row.names(predicted_clusters), row.names(real_clusters))

if(predicted_clusters$predicted_clusters[1] == 2)
{
  predicted_clusters$predicted_clusters[which(predicted_clusters$predicted_clusters == 1)]=0
  predicted_clusters$predicted_clusters[which(predicted_clusters$predicted_clusters == 2)]=1
  predicted_clusters$predicted_clusters[which(predicted_clusters$predicted_clusters == 0)]=2
}

#---
merged = merge(predicted_clusters, real_clusters[,c('yan', 'Donor')], by = 'row.names')

freq_predicted = data.frame(table(merged$predicted_clusters), stringsAsFactors = F)
freq_predicted = freq_predicted[order(freq_predicted$Freq, decreasing = T),]
freq_predicted

freq_real = data.frame(table(merged$yan), stringsAsFactors = F)
freq_real = freq_real[order(freq_real$Freq, decreasing = T),]
freq_real

merged$predicted_clusters_unlabelled = merged$predicted_clusters

merged$predicted_clusters[which(merged$predicted_clusters == as.integer(freq_predicted[1,1]))] = as.character(freq_real[1,1])
merged$predicted_clusters[which(merged$predicted_clusters == as.integer(freq_predicted[2,1]))] = as.character(freq_real[2,1])
table(merged$predicted_clusters)
#---

table = table(merged$predicted_clusters, merged$yan)
# TPrate = (table[1,1] + table[2,2])/length(clusters(mixed_model))
TPrate = (table[1,1])/length(which(clusters(mixed_model) == as.integer(freq_predicted[1,1])))
print(TPrate)

# freq_c1 = data.frame(table(merged[which(merged$predicted_clusters == 1),"yan"]), stringsAsFactors = F)
# freq_c1 = freq_c1[order(freq_c1$Freq, decreasing = T),]
# freq_c1
# 
# freq_c2 = data.frame(table(merged[which(merged$predicted_clusters == 2), "yan"]), stringsAsFactors = F)
# freq_c2 = freq_c2[order(freq_c2$Freq, decreasing = T),]
# freq_c2
# 
# merged[which(merged$predicted_clusters == 1), "predicted_clusters"] = as.character(freq_c1$Var1[1])
# merged[which(merged$predicted_clusters == 2), "predicted_clusters"] = as.character(freq_c2$Var1[1])
# 
# table(merged$predicted_clusters == merged$yan)
# 
merged_subset = merged[which(merged$predicted_clusters == merged$yan),]
table(merged_subset$yan)

# if I have two distinct clusters, otherwise I keep all the cells as they are
if(length(unique(merged_subset$yan)) > 1)
{
  
  #--------------------------- keep cells which generate distinct clusters
  target_Data_Seurat = subset(target_Data_Seurat, cells = merged_subset$Row.names) # if I proceed with predicted cluster I do get 25 DE genes. The actual cell give me around 253.
  #target_Data_Seurat@meta.data$yan = merged$predicted_clusters
  
  DimPlot(target_Data_Seurat, group.by = 'yan', reduction = "umap", label = F, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))
  
  #---------------------------------- pairwise analysis (select donors which present H and L both)
  # Donors = data.frame(table(gsub('H|L','',unique(target_Data_Seurat@meta.data$Donor_subtype))))
  # Donors_dup = Donors[which(Donors$Freq > 1),]
  # 
  # MetaData_Selected = MetaData[which(MetaData$Donor %in% as.character(Donors_dup$Var1)),]
  # 
  # target_Data_Seurat <- subset(target_Data_Seurat, cells = MetaData_Selected$Barcode)
  # target_Data_Seurat
  
  table(target_Data_Seurat@meta.data$yan)
}

#---- clean UMAP plot (remove intermediate cells)
# umap = target_Data_Seurat@reductions$umap@cell.embeddings
# umap = as.data.frame(umap)
# 
# umap_V = (umap[which(row.names(umap) %in% names(new.cluster.ids[which(new.cluster.ids == "Vulnerable" )]) ),])
# umap_V_subset = umap_V[which(umap_V$umap_1 > 0),]
# 
# umap_R = (umap[which(row.names(umap) %in% names(new.cluster.ids[which(new.cluster.ids == "Resistant" )]) ),])
# umap_R_subset = umap_R[which(umap_R$umap_1 < 0),]
# 
# Selected_Cells = c(row.names(umap_R_subset), row.names(umap_V_subset))
# 
# target_Data_Seurat = subset(target_Data_Seurat, cells = Selected_Cells)
# 
# DimPlot(target_Data_Seurat, reduction = "umap", label = TRUE, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))
# table(target_Data_Seurat@meta.data$yan)
# 
# length(unique(target_Data_Seurat@meta.data$Donor_subtype))

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
setwd(outputDir)
pdf(paste0(outputDir, 'UMAP_VR_Neurons_SC.pdf'), height = 8, width = 10, useDingbats = F)
DimPlot(target_Data_Seurat, reduction = "umap", label = F, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))
dev.off() 

# Set color palette
library(viridis)
pal <- viridis(n = 10, option = "D")

library(scCustomize)
pdf(paste0(outputDir, 'UMAP_VR_Neurons_SC_scCustom.pdf'), height = 8, width = 10, useDingbats = F)
DimPlot_scCustom(seurat_object = target_Data_Seurat, figure_plot = TRUE, repel = T, pt.size = 2, label.size = 8, label.box = T, colors_use = c("#B4DE2CFF", "#FDE725FF"))
dev.off() 

#---------------------------------- Perform aggregation of counts and metadata by subject and cell type.
library(scater)
target_Data_sce <- as.SingleCellExperiment(target_Data_Seurat)

table(target_Data_sce$yan)
# Resistant Vulnerable 
# 852         82 

# included 14 donors with PD.
dim(table(target_Data_sce$Donor_subtype))

library(aggregateBioVar)
aggregate_counts = aggregateBioVar(scExp = target_Data_sce, cellVar = "yan", subjectVar = "Donor_subtype")

State = as.character(colData(aggregate_counts$AllCells)$Donor_subtype)
State[grep('V', State)] = 'Vulnerable'
State[grep('R', State)] = 'Resistant'
colData(aggregate_counts$AllCells)$State = State

#---------------------------------- DESeq2Aggregate----------------------------------------------------------
library(DESeq2)

coldata <- DataFrame(State=relevel(factor(colData(aggregate_counts$AllCells)$State), ref = "Resistant"), subject = colData(aggregate_counts$AllCells)$Donor )
#View(data.frame(State=relevel(factor(colData(aggregate_counts$AllCells)$State), ref = "Resistant"), subject = colData(aggregate_counts$AllCells)$Donor ))

subj_dds_dataset <-
  DESeqDataSetFromMatrix(
countData = assay(aggregate_counts$AllCells, "counts"),
colData = coldata,
design = ~ subject + State) 

#~~~~~~~~~
# ~ subject + State # pairwise analysis

# NOTE: if you use "test="LRT", reduced=~1" in DESeq, you get more results if you use "subject + State" in DESeqDataSetFromMatrix
#~~~~~~~~~

subj_dds <- DESeq(subj_dds_dataset)

#~ optemize DESeq
# subj_dds <- DESeq(subj_dds_dataset, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, fitType='local') # , fitType='local'
# res_LRT <- results(subj_dds)
# rld_LRT<- rlogTransformation(subj_dds)
# pathsLRT<-assay(rld_LRT)
# df_pathLRT <- cbind(rownames(res_LRT), data.frame(res_LRT, row.names=NULL))
# topTableLRT <- as.data.frame(df_pathLRT)
# sigGeneListLRT <- subset(topTableLRT,  padj<=0.05)[,1]
# topMatrixLRT <- pathsLRT[which(rownames(pathsLRT) %in% sigGeneListLRT),]
# Use test="LRT" for significance testing when working with single-cell data, over the Wald test.
# ref: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html and https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
# Since our ‘full’ model only has one factor (sampletype), the ‘reduced’ model is just the intercept.
# The LRT is comparing the full model to the reduced model to identify significant genes. The p-values are determined solely by the difference in deviance between the ‘full’ and ‘reduced’ model formula (not log2 fold changes). Essentially the LRT test is testing whether the term(s) removed in the ‘reduced’ model explains a significant amount of variation in the data?
# local - use the locfit package to fit a local regression of log dispersions over log base mean (normal scale means and dispersions are input and output fordispersionFunction). The points are weighted by normalized mean count in the local regression.
# https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/estimateDispersions
#~

subj_dds_results <-
  results( subj_dds, pAdjustMethod = "fdr", contrast = c("State", "Vulnerable", "Resistant") ) # e.g., State PD [Vulnerable] vs Control
# inserting LFC thresholding into the specification of the null hypothesis: lfcThreshold=1
# ref: https://support.bioconductor.org/p/9142755/
subj_dds_transf <- as.data.frame(subj_dds_results) %>%
  bind_cols(feature = rownames(subj_dds_results)) %>%
  mutate(log_padj = -log(subj_dds_results$padj, base = 10))

table(subj_dds_transf$padj < 0.05)
table(subj_dds_transf$log2FoldChange[which(subj_dds_transf$padj < 0.05)] > 0) # up and down regulation

subj_dds_transf[which(subj_dds_transf$feature == 'TUBB'),]

#--------- effect size (Log fold change (LFC) estimates)
# resultsNames(subj_dds) # lists the coefficients
# res <- results(subj_dds, name="State_Vulnerable_vs_Resistant")
# # or to shrink log fold changes association with condition:
# res <- lfcShrink(subj_dds, coef="State_Vulnerable_vs_Resistant", type="apeglm")
# View(as.data.frame(res))
# 
# x <- model.matrix(design(subj_dds), colData(subj_dds))
# param <- dispersions(subj_dds)
# mle <- log(2) * cbind(subj_dds_results$log2FoldChange, subj_dds_results$lfcSE)
# offset <- matrix(log(sizeFactors(subj_dds)),
#                  ncol=ncol(subj_dds),
#                  nrow=nrow(subj_dds),byrow=TRUE)
# 
# res.shr <- lfcShrink(subj_dds, coef=5, type="normal")

# 'COBL' and 'TUBB'

# Count_aggregated = as.data.frame(aggregate_counts$CellType@assays@data$counts)
# dim(Count_aggregated)
# Count_aggregated = data.frame(GeneSymbol = row.names(Count_aggregated), Count_aggregated)

#- check direction of expression
counts_aggregate_Vulnerable = as.data.frame(aggregate_counts$Vulnerable@assays@data$counts)
counts_aggregate_Resistant = as.data.frame(aggregate_counts$Resistant@assays@data$counts)

counts_aggregate_Vulnerable$sum = rowSums(counts_aggregate_Vulnerable)
counts_aggregate_Resistant$sum = rowSums(counts_aggregate_Resistant)

#-- check if direction of expressin in DEGs profile is accurate
counts_aggregate_Vulnerable[which(row.names(counts_aggregate_Vulnerable) == 'TUBB'),]
counts_aggregate_Resistant[which(row.names(counts_aggregate_Resistant) == 'TUBB'),]

dim(counts_aggregate_Vulnerable)
dim(counts_aggregate_Resistant)

subsetting_up = counts_aggregate_Vulnerable
subsetting_down = counts_aggregate_Vulnerable

identical(gsub('V', '', colnames(counts_aggregate_Vulnerable)), gsub('R', '', colnames(counts_aggregate_Resistant)) )

index = 1
for(index in 1:dim(counts_aggregate_Vulnerable)[2])
{
  subsetting_up[,index] = counts_aggregate_Vulnerable[,index] >= counts_aggregate_Resistant[,index]
  subsetting_down[,index] = counts_aggregate_Vulnerable[,index] <= counts_aggregate_Resistant[,index]
}

subsetting_up = subsetting_up[apply(subsetting_up, 1, function(x) all(x==TRUE)),]
subsetting_down = subsetting_down[apply(subsetting_down, 1, function(x) all(x==TRUE)),]

DEGs_up = subj_dds_transf[which(subj_dds_transf$log2FoldChange > 0 & subj_dds_transf$padj < 0.05),]
DEGs_down = subj_dds_transf[which(subj_dds_transf$log2FoldChange < 0 & subj_dds_transf$padj < 0.05),]

dim(DEGs_up)
dim(DEGs_down)

length(which(rownames(DEGs_up) %in% row.names(subsetting_up)))
length(which(rownames(DEGs_down) %in% row.names(subsetting_down)))

#-------------------------------- visualization [heatmap]
library(hgnc)
library(data.table)
library(ComplexHeatmap)

# Transformation
dds_t <- estimateSizeFactors(subj_dds_dataset)
baseMean <- rowMeans(counts(dds_t, normalized=TRUE))
sum(baseMean > 1)
idx <- sample(which(baseMean > 5), 100)
dds_t.sub <- dds_t[idx, ]
dim(dds_t.sub)
dds_t.sub <- estimateDispersions(dds_t.sub)
dispersionFunction(dds_t) <- dispersionFunction(dds_t.sub)
vsd <- varianceStabilizingTransformation(dds_t, blind=FALSE)

ddsTF <- getVarianceStabilizedData(dds_t)
class(ddsTF)
ddsTF <- as.data.frame(ddsTF)

#------ Count
counts_aggregate_upregulated = as.data.frame(aggregate_counts$AllCells@assays@data$counts)
counts_aggregate_downregulated = as.data.frame(aggregate_counts$AllCells@assays@data$counts)

subj_dds_transf = subj_dds_transf[order(subj_dds_transf$padj, decreasing = F),]

counts_aggregate_upregulated = counts_aggregate_upregulated[which(row.names(counts_aggregate_upregulated) %in% subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05 & subj_dds_transf$log2FoldChange > 0)]),]
counts_aggregate_downregulated = counts_aggregate_downregulated[which(row.names(counts_aggregate_downregulated) %in% subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05 & subj_dds_transf$log2FoldChange <= 0)]),]

#~~~~~~~~~ replace raw reads with dipersion
counts_aggregate_upregulated = ddsTF[which(row.names(ddsTF) %in% row.names(counts_aggregate_upregulated)),which(colnames(ddsTF) %in% colnames(counts_aggregate_upregulated))]
counts_aggregate_downregulated = ddsTF[which(row.names(ddsTF) %in% row.names(counts_aggregate_downregulated)),which(colnames(ddsTF) %in% colnames(counts_aggregate_downregulated))]
#~~~~~~~~~

counts_aggregate_upregulated = counts_aggregate_upregulated[match(subj_dds_transf$feature[which(subj_dds_transf$feature %in% row.names(counts_aggregate_upregulated))], row.names(counts_aggregate_upregulated)),]
identical(row.names(counts_aggregate_upregulated), subj_dds_transf$feature[which(subj_dds_transf$feature %in% row.names(counts_aggregate_upregulated))])

counts_aggregate_downregulated = counts_aggregate_downregulated[match(subj_dds_transf$feature[which(subj_dds_transf$feature %in% row.names(counts_aggregate_downregulated))], row.names(counts_aggregate_downregulated)),]
identical(row.names(counts_aggregate_downregulated), subj_dds_transf$feature[which(subj_dds_transf$feature %in% row.names(counts_aggregate_downregulated))])

dim(counts_aggregate_upregulated)
dim(counts_aggregate_downregulated)

#counts_aggregate_downregulated = counts_aggregate_downregulated[which(row.names(counts_aggregate_downregulated) %in% rownames(DEGs_down) ),]

#------ Count subset
CountMatrix = rbindlist(list(counts_aggregate_upregulated[1:10,], counts_aggregate_downregulated[1:10,]))
CountMatrix = as.data.frame(CountMatrix)
rownames(CountMatrix) = c(row.names(counts_aggregate_upregulated)[1:10], row.names(counts_aggregate_downregulated)[1:10])

dim(CountMatrix)
colnames(CountMatrix) 

#------ heatmap
State = as.character(colnames(CountMatrix))
State[grep('V', State)] = 'Vulnerable Neurons'
State[grep('R', State)] = 'Resistant Neurons'

annotation_col = data.frame(
  State = factor(State, levels = c('Vulnerable Neurons', 'Resistant Neurons'))
)
rownames(annotation_col) = as.character(colnames(CountMatrix))
str(annotation_col)

annotation_row = data.frame(
  GeneClass = factor(rep(c('Up-regulated', 'Down-regulated'), c(10, 10)))
)

rownames(annotation_row) = c(row.names(counts_aggregate_upregulated)[1:10], row.names(counts_aggregate_downregulated)[1:10])
identical(row.names(CountMatrix), row.names(annotation_row))

row.names(subsetting_down)

ann_colors = list(
  State = c(`Vulnerable Neurons` = "#1B9E77", `Resistant Neurons` = "#D95F02"),
  GeneClass = c(`Up-regulated` = "#7570B3", `Down-regulated` = "#E7298A")
)

#-- revise colnames
colnames(CountMatrix)[grep('V', colnames(CountMatrix))] = rep(paste0('Donor-', 1:length(grep('V', colnames(CountMatrix))), '-V' ), each = 1)
colnames(CountMatrix)[grep('R', colnames(CountMatrix))] = rep(paste0('Donor-', 1:length(grep('R', colnames(CountMatrix))), '-R' ), each = 1)
#--

p <- pheatmap(CountMatrix, annotation_col = annotation_col, annotation_row = annotation_row, 
              annotation_colors = ann_colors, 
              row_split = annotation_row$GeneClass,
              column_split = annotation_col$State, name = "Expression", scale='row',
               fontsize_number=8, display_numbers = round(as.matrix(CountMatrix), digits = 1), cluster_cols=FALSE  )
p

#~~~~~~~~~
# display_numbers = as.matrix(CountMatrix),
#~~~~~~~~~
# ref: https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/

dir.create('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/', recursive = T)
pdf(paste0('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/Neurons_pheatmap.pdf'), width = 12, height = 8,useDingbats = FALSE)
print(p)
dev.off()

#-------------------------------- annotation
library(dplyr)
library(annotables)
#devtools::install_github("stephenturner/annotables")
grch38_tx2gene
grch38
# convert to the tibble
grch38 = tbl_df(grch38)

Annotated_Genes = grch38 %>% 
  dplyr::filter(symbol %in% as.character(subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05)])) %>% 
  dplyr::select(symbol, description)

subj_dds_transf = subj_dds_transf[order(subj_dds_transf$pvalue, decreasing = F),]
subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05)]

Annotated_Genes = Annotated_Genes[match(subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05)], Annotated_Genes$symbol),]
head(Annotated_Genes)
colnames(Annotated_Genes)[1] = 'feature'

Annotated_Genes = left_join(subj_dds_transf, Annotated_Genes, by = 'feature' )
fwrite(Annotated_Genes, '/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/Annotated_DEGenes.txt', quote = F, sep = '\t', row.names = F)

#--
CountMatrix = rbindlist(list(counts_aggregate_upregulated[1:10,], counts_aggregate_downregulated[1:10,]))
CountMatrix = as.data.frame(CountMatrix)
rownames(CountMatrix) = c(row.names(counts_aggregate_upregulated)[1:10], row.names(counts_aggregate_downregulated)[1:10])

colnames(CountMatrix)[grep('V', colnames(CountMatrix))] = rep(paste0('Donor-', 1:length(grep('V', colnames(CountMatrix))), '-V' ), each = 1)
colnames(CountMatrix)[grep('R', colnames(CountMatrix))] = rep(paste0('Donor-', 1:length(grep('R', colnames(CountMatrix))), '-R' ), each = 1)

Annotated_Genes_sub = Annotated_Genes[which(Annotated_Genes$feature %in% row.names(CountMatrix)),]
dim(Annotated_Genes_sub)
dim(CountMatrix)

CountMatrix_sub = CountMatrix[match(Annotated_Genes_sub$feature, row.names(CountMatrix)),]
identical(row.names(CountMatrix_sub), Annotated_Genes_sub$feature)

row.names(CountMatrix) = paste0(row.names(CountMatrix), ' (', Annotated_Genes_sub$description, ')')

p <- pheatmap(CountMatrix, annotation_col = annotation_col, annotation_row = annotation_row, 
              annotation_colors = ann_colors, 
              row_split = annotation_row$GeneClass,
              column_split = annotation_col$State, name = "Expression", scale='row',
              fontsize_number=8, display_numbers = round(as.matrix(CountMatrix), digits = 1), cluster_cols=FALSE, legend = FALSE, annotation_legend = T  )
p

pdf(paste0('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/Neurons_pheatmap_with_description.pdf'), width = 12, height = 8,useDingbats = FALSE)
print(p)
dev.off()

#-- summary of MetaData
length(unique(MetaData$Donor))
length(unique(MetaData$Donor[which(MetaData$yan == "Resistant" )]))
length(unique(MetaData$Donor[which(MetaData$yan == "Vulnerable" )]))

Selected_Donors = intersect(unique(MetaData$Donor[which(MetaData$yan == "Resistant" )]), unique(MetaData$Donor[which(MetaData$yan == "Vulnerable" )]))
MetaData_Selected = MetaData[which(MetaData$Donor %in% Selected_Donors),]

length(unique(MetaData_Selected$Donor))
length(unique(MetaData_Selected$Donor[which(MetaData_Selected$yan == "Resistant" )]))
length(unique(MetaData_Selected$Donor[which(MetaData_Selected$yan == "Vulnerable" )]))
table(MetaData_Selected$yan)
dim(MetaData_Selected)
table(duplicated(MetaData_Selected$Barcode))
table(MetaData_Selected$Donor)

library(dplyr)
R = MetaData %>%
  group_by(Donor, DataSet, yan) %>%
  summarize(count=n())


#--------- comparison of DEG with GeoMx_DEA
GeoMx_DEA = fread("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/GeoMx_DEA.txt")
GeoMx_DEA = as.data.frame(GeoMx_DEA)
dim(GeoMx_DEA)

length(intersect(subj_dds_transf$feature, GeoMx_DEA$TargetName))
length(intersect(subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05)], GeoMx_DEA$TargetName))

#--------- Comapre with GWAS
# ref https://pubmed.ncbi.nlm.nih.gov/36777638/
subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05 & subj_dds_transf$feature %in% c('SYT1', 'NEFM', 'NEFL', 'SNAP25', 'GAP43', 'GRIA1'))]

GWAS_hits = fread('/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis/GWAS_hits.txt', stringsAsFactors = F, header = F)
subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05 & subj_dds_transf$feature %in% GWAS_hits$V1)]

#-------------------------------- pathway Enrichment

#--- gene enrichment analysis
library('org.Hs.eg.db')
BD = 'c5.all.v2023.1.Hs.entrez.gmt'

library(qusage)
library(clusterProfiler)
gmtfile <- paste0("/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Scripts_USED_visium/msigdb_v2023_Hs/", BD)
c5 <- read.gmt(gmtfile)

# use mapIds method to obtain Entrez IDs
input = subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05)]
de = as.character(mapIds(org.Hs.eg.db, input, 'ENTREZID', 'SYMBOL'))
table(c5$term)
c5$term = sub("^[^_]*_", "", c5$term)

egmt <- enricher(de[!is.na(de)], TERM2GENE=c5)
 
library(cowplot)
library(ggplot2)
p2 <- dotplot(egmt, showCategory=5) + ggtitle("")
p2

pdf(paste0('PathwayEnrichment.pdf') ,width = 10, height = 10, useDingbats = FALSE)
print(plot_grid(p2, ncol=1))
dev.off()

#--- gene enrichment analysis
library("DOSE")
edo <- enrichDGN(de[!is.na(de)])
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

dev.off()
library(enrichplot)
barplot(edo, showCategory=10)
dotplot(edo, showCategory=15) + ggtitle("dotplot for ORA")
cnetplot(edox, circular = TRUE, colorEdge = TRUE)

pdf(paste0('DiseaseEnrichment.pdf') ,width = 8, height = 9, useDingbats = FALSE)
print(dotplot(edo, showCategory=15) + ggtitle("dotplot for ORA"))
dev.off()


#---- pathview - visualization 

library('org.Hs.eg.db')
keytypes(org.Hs.eg.db)

input = subj_dds_transf[which(subj_dds_transf$padj < 0.05),]
input = input[-which(is.na(as.character(mapIds(org.Hs.eg.db, input$feature, 'ENTREZID', 'SYMBOL')))),]
row.names(input) = as.character(mapIds(org.Hs.eg.db, input$feature, 'ENTREZID', 'SYMBOL'))

getwd()
library("pathview")
###################################################
### code chunk number 15: kegg.native
###################################################
# pv.out <- pathview(gene.data = names(input), pathway.id = '05012',
#    species = "hsa", kegg.native = T,
#    same.layer = F, high=list(gene="gold"), out.suffix = 'VR')
# head(pv.out$plot.data.gene)
#  
# intersect(as.character(input), pv.out$plot.data.gene$labels)
# 
# pv.out <- pathview(gene.data = names(input), pathway.id = '05014',
#                    species = "hsa", kegg.native = T,
#                    same.layer = F, high=list(gene="gold"), out.suffix = 'VR')

input_used = input$log2FoldChange
names(input_used) = row.names(input)

pv.out <- pathview(gene.data = input_used, pathway.id = '05012',
   species = "hsa", kegg.native = T,
   same.layer = F, high=list(gene="gold"), out.suffix = 'VR')

list_pathways = c('05010', '05012', '05014', '05016', '05017', '05020', '05022')

for(p in 1:length(list_pathways))
{
  pv.out <- pathview(gene.data = input_used, pathway.id = list_pathways[p],
                     species = "hsa", kegg.native = T,
                     same.layer = F, high=list(gene="gold"), out.suffix = 'VR')
}

# 05010 N
# Alzheimer disease
# 05012 N
# Parkinson disease
# 05014 N
# Amyotrophic lateral sclerosis
# 05016 N
# Huntington disease
# 05017 N
# Spinocerebellar ataxia
# 05020 N
# Prion disease
# 05022 N
# Pathways of neurodegeneration - multiple diseases








