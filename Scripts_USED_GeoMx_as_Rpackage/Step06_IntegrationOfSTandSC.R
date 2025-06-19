####################################################################################################################
################# Prepare input files for integrating single-cell and spatial transcriptomics ######################
####################################################################################################################

#'@import data.table
#'@import dplyr
#'@import Seurat
#'@import flexmix
#'@import mvtnorm
#'@import viridis
#'@import scCustomize
#'@import scater
#'@import aggregateBioVar
#'@import DESeq2
#'@import SingleCellExperiment
#'@import scmap
#'@export
#'@name IntegSTandSC
#'@title Integrating single-cell and spatial transcriptomics
#'@description This function integrats single-cell and spatial transcriptomics
#'@author {Isar Nassiri}
#'@param InputDir
#'A folder including expression matrix of ST and SC.
#'@param OutputDir
#'You can file the output files in this folder.
#'@return You can find the matrix of gene expression in OutputDir folder.
#'@examples
#'library(data.table)
#'library(dplyr)
#'library(Seurat)
#'library(flexmix)
#'library(mvtnorm)
#'library(scater)  
#'library(DESeq2)
#'library(hgnc)
#'library(ComplexHeatmap)
#'library(annotables)
#'library(scCustomize)
#'InputDir = "/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/DataAnalysis_Output/"
#'OutputDir = "/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/DataAnalysis_Output/"
#'IntegSTandSC(PI, InputDir, OutputDir)
#'@export

#----------------------------------
# PI='GSE184950'
# InputDir = "/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/DataAnalysis_Output/"
# OutputDir = "/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/DataAnalysis_Output/"
# IntegSTandSC(PI, InputDir, OutputDir)
#----------------------------------

IntegSTandSC <- NULL
IntegSTandSC <- function(PI, InputDir, OutputDir)
{
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(flexmix)
  library(mvtnorm)
  library(scater)  
  library(DESeq2)
  library(hgnc)
  library(ComplexHeatmap)
  library(annotables)
  library(scCustomize)
  
  MetaData = fread(paste0(InputDir, '/', PI, '_Neurons_Subset_MetaData.txt'), stringsAsFactors = F, header = T)
  MetaData = as.data.frame(MetaData)
  MetaData$DataSet = PI

  GeneCount = fread(paste0(InputDir, '/', PI, '_Neurons_Subset_GeneCount.txt'), stringsAsFactors = F, header = T)
  GeneCount = as.data.frame(GeneCount)

  row.names(GeneCount) = GeneCount$GeneSymbol
  GeneCount = GeneCount[,-1]
  
  GeneCount = GeneCount[apply(GeneCount, 1, function(x) !all(x==0)),]
  GeneCount = GeneCount[apply(GeneCount, 2, function(x) !all(x==0)),]
  
  colnames(GeneCount) = gsub('\\.', '-', colnames(GeneCount))
  
  #---------------------------------- Seurat analysis
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
  identical(names(Idents(object = target_Data_Seurat)), target_Data_Seurat@meta.data$Barcode)
  
  new.cluster.ids = target_Data_Seurat@meta.data$yan
  names(new.cluster.ids) = names(Idents(object = target_Data_Seurat))
  Idents(object = target_Data_Seurat) = factor(new.cluster.ids, levels = c("Resistant", "Vulnerable"))
  
  #-------------------------- clean UMAP plot (remove intermediate cells)
  predicted_clusters = target_Data_Seurat@reductions$umap@cell.embeddings
  predicted_clusters = as.data.frame(predicted_clusters)
  
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
  TPrate = (table[1,1])/length(which(clusters(mixed_model) == as.integer(freq_predicted[1,1])))
  print(TPrate)
  
  merged_subset = merged[which(merged$predicted_clusters == merged$yan),]
  table(merged_subset$yan)
  
  # if I have two distinct clusters, otherwise I keep all the cells as they are
  if(length(unique(merged_subset$yan)) > 1)
  {
    
    #--------------------------- keep cells which generate distinct clusters
    target_Data_Seurat = subset(target_Data_Seurat, cells = merged_subset$Row.names) # if I proceed with predicted cluster I do get 25 DE genes. The actual cell give me around 253.

    DimPlot(target_Data_Seurat, group.by = 'yan', reduction = "umap", label = F, repel = T, pt.size=1.3, label.size = 8, label.box = T, cols = c('#66CCCC', '#6699CC'))

    table(target_Data_Seurat@meta.data$yan)
  }
  
  #---------------------------------- pairwise analysis (select donors which present V and R both)
  Donors = data.frame(table(gsub('R|V','',unique(target_Data_Seurat@meta.data$Donor_subtype))))
  Donors_dup = Donors[which(Donors$Freq > 1),]
  
  MetaData_Selected = MetaData[which(MetaData$Donor %in% as.character(Donors_dup$Var1)),]
  
  target_Data_Seurat <- subset(target_Data_Seurat, cells = MetaData_Selected$Barcode)
  target_Data_Seurat
  
  table(target_Data_Seurat@meta.data$yan)
  
  #---------------------------------- Perform aggregation of counts and metadata by subject and cell type.
  
  target_Data_sce <- as.SingleCellExperiment(target_Data_Seurat)
  
  table(target_Data_sce$yan)
  dim(table(target_Data_sce$Donor_subtype))
  
  library(aggregateBioVar)
  aggregate_counts = aggregateBioVar(scExp = target_Data_sce, cellVar = "yan", subjectVar = "Donor_subtype")
  
  State = as.character(colData(aggregate_counts$AllCells)$Donor_subtype)
  State[grep('V', State)] = 'Vulnerable'
  State[grep('R', State)] = 'Resistant'
  colData(aggregate_counts$AllCells)$State = State
  
  #---------------------------------- DESeq2Aggregate----------------------------------------------------------
  coldata <- DataFrame(State=relevel(factor(colData(aggregate_counts$AllCells)$State), ref = "Resistant"), subject = colData(aggregate_counts$AllCells)$Donor )
  
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
  
  if(length(grep('TRUE', names(table(subj_dds_transf$padj < 0.05)) )) > 0)
  {
  
  #-------------------------------- annotation
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
  
  setwd(OutputDir)
  fwrite(Annotated_Genes, paste0(PI, '_Annotated_DEGenes.txt'), quote = F, sep = '\t', row.names = F)
  }

  cat(
    paste0(
      "\033[0;",
      47,
      "m",
      "You can find the results in: ",
      "\033[0m",
      "\n",
      OutputDir
    )
  )
} 
