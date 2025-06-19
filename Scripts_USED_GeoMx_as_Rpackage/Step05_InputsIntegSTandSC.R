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
#'@name InputsIntegSTandSC
#'@title Prepare input files for integrating single-cell and spatial transcriptomics
#'@description This function prepares input files for integrating single-cell and spatial transcriptomics
#'@author {Isar Nassiri}
#'@param InputDir
#'A folder including expression profiles of SC.
#'@param OutputDir
#'You can file the output files in this folder.
#'@param PI
#'A name for the dataset
#'@param AnnotationDir
#'A folder including MetaData_SC.txt and assignmentTable.txt files
#'@param Region
#'Region of interest in the MetaData_SC.txt file.
#'@param Diagnosis
#'State of interest in the MetaData_SC.txt file.
#'@return You can find the matrix of gene expression in OutputDir folder.
#'@examples
#'library(data.table)
#'library(dplyr)
#'library(Seurat)
#'library(flexmix)
#'library(mvtnorm)
#'library(viridis)
#'library(scCustomize)
#'library(scater)
#'library(aggregateBioVar)
#'library(DESeq2)
#'library(SingleCellExperiment)
#'library(scmap)
#'PI='GSE184950'
#'InputDir = paste0('/Users/isarnassiri/Documents/Parkinson/Data/',PI,'_RAW/')
#'AnnotationDir='/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/annotation/'
#'OutputDir = "/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/DataAnalysis_Output/"
#'Region='SN'
#'Diagnosis='PD'
#'InputsIntegSTandSC(PI, InputDir, OutputDir, AnnotationDir, Region, Diagnosis)
#'@export

#----------------------------------
# PI='GSE184950'
# InputDir = paste0('/Users/isarnassiri/Documents/Parkinson/Data/',PI,'_RAW/')
# AnnotationDir='/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/annotation/'
# OutputDir = "/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/DataAnalysis_Output/"
# Region='SN'
# Diagnosis='PD'
# InputsIntegSTandSC(PI, InputDir, OutputDir, AnnotationDir, Region, Diagnosis)
#----------------------------------

InputsIntegSTandSC <- NULL
InputsIntegSTandSC <- function(PI, InputDir, OutputDir, AnnotationDir, Region, Diagnosis)
{

  library(data.table)
  library(dplyr)
  library(Seurat)
  library(flexmix)
  library(mvtnorm)
  library(viridis)
  library(scCustomize)
  library(scater)
  library(aggregateBioVar)
  library(DESeq2)
  library(SingleCellExperiment)
  library(scmap)
 
  #------------------------- MetaData
  MetaData = fread(paste0(AnnotationDir, '/MetaData_SC.txt'), stringsAsFactors = F, header = T)
  MetaData = as.data.frame(MetaData)
  MetaData = MetaData[which(MetaData$Region == Region & MetaData$Diagnosis %in% Diagnosis),]
  
  #------------------------- CellType Annotations
  CellTypeEnrichment = fread(paste0(AnnotationDir, 'assignmentTable.txt'), stringsAsFactors = F, header = T)
  CellTypeEnrichment = as.data.frame(CellTypeEnrichment)
  
  CellTypeEnrichment_subset = CellTypeEnrichment[which(CellTypeEnrichment$Sample_ID %in% MetaData$Sample_ID),]
  
  #------------------------- read the scRNAseq profile
  setwd(InputDir)
  
  ListOfFolders_All = list.dirs(full.names = F, recursive = F)
  ListOfFolders = ListOfFolders_All[which(ListOfFolders_All %in% MetaData$Sample_ID)]
  
  i=n=p=1
  for(i in 1:length(ListOfFolders))
  {
    
    Diagnosis = MetaData$Diagnosis[which(MetaData$Sample_ID == ListOfFolders[i])]
    
    if(Diagnosis == "PD")
    {
      eislet <- Read10X(paste(ListOfFolders[i],'/filtered_feature_bc_matrix', sep=''))
      eislet = as.data.frame(eislet)
      
      counts.fname = as.matrix(eislet)
      rm(eislet)
      
      print(dim(counts.fname))
      counts.fname = counts.fname[, which(colnames(counts.fname) %in% gsub('\\-.*','-1', CellTypeEnrichment$cell))]
      print(dim(counts.fname))
      
      if(p!=1){
        counts.fname = counts.fname[,-which(colnames(counts.fname) %in% colnames(TP_profile))]
      }
      
      Barcode_Donor = data.frame(Barcode = colnames(counts.fname), Donor = rep(ListOfFolders[i], length(colnames(counts.fname))) )
      
      if(p==1){TP_profile = counts.fname; p=p+1; Barcode_Donor_All = Barcode_Donor}else{TP_profile = merge(TP_profile, counts.fname, by = "row.names"); row.names(TP_profile) = TP_profile$Row.names; TP_profile = TP_profile[,-which(colnames(TP_profile)[1] == "Row.names")]; Barcode_Donor_All = rbind(Barcode_Donor_All, Barcode_Donor)}
      
      rm(counts.fname)
    }
    
    if(Diagnosis == "normal")
    {
      eislet <- Read10X(paste(ListOfFolders[i],'/filtered_feature_bc_matrix', sep=''))
      eislet = as.data.frame(eislet)
      
      counts.fname = as.matrix(eislet)
      rm(eislet)
      
      print(dim(counts.fname))
      counts.fname = counts.fname[, which(colnames(counts.fname) %in% gsub('\\-.*','-1',CellTypeEnrichment$cell))]
      print(dim(counts.fname))
      
      if(n!=1){
        counts.fname = counts.fname[,-which(colnames(counts.fname) %in% colnames(TP_profile_normal))]
      }
      
      if(n==1){TP_profile_normal = counts.fname; n=n+1;}else{TP_profile_normal = merge(TP_profile_normal, counts.fname, by = "row.names"); row.names(TP_profile_normal) = TP_profile_normal$Row.names; TP_profile_normal = TP_profile_normal[,-which(colnames(TP_profile_normal)[1] == "Row.names")]}
      rm(counts.fname)
    }
  }
  
  if(Diagnosis == "PD")
  {
    TP_profile = TP_profile[apply(TP_profile, 1, function(x) !all(x==0)),]
    TP_profile = TP_profile[apply(TP_profile, 2, function(x) !all(x==0)),]
    dim(TP_profile)
  }
  
  if(Diagnosis == "normal")
  {
    TP_profile_normal = TP_profile_normal[apply(TP_profile_normal, 1, function(x) !all(x==0)),]
    TP_profile_normal = TP_profile_normal[apply(TP_profile_normal, 2, function(x) !all(x==0)),]
    dim(TP_profile_normal)
  }
  
  #--- save PD ---
  CellTypeEnrichment_subset = CellTypeEnrichment_subset[which(CellTypeEnrichment_subset$geneSet == 'Neurons'),]
  
  CellTypeEnrichment_subset_PD = CellTypeEnrichment_subset[which(CellTypeEnrichment_subset$cell_orig %in% colnames(TP_profile)),]
  CellTypeEnrichment_subset_PD = CellTypeEnrichment_subset_PD[!duplicated(CellTypeEnrichment_subset$cell_orig),]
  CellTypeEnrichment_subset_PD = CellTypeEnrichment_subset_PD[!is.na(CellTypeEnrichment_subset_PD$cell_orig),]
  
  dim(TP_profile)
  if(length(setdiff(CellTypeEnrichment_subset$cell_orig, colnames(TP_profile) )) > 0)
  {
    TP_profile = TP_profile[,-which(colnames(TP_profile) %in% setdiff(colnames(TP_profile), CellTypeEnrichment_subset_PD$cell_orig) )]
  }
  
  CellTypeEnrichment_subset_PD = CellTypeEnrichment_subset_PD[match(colnames(TP_profile), CellTypeEnrichment_subset_PD$cell_orig),]
  identical(gsub('\\-.*','-1',CellTypeEnrichment_subset_PD$cell), colnames(TP_profile) )
  dim(CellTypeEnrichment_subset_PD)
  dim(TP_profile)
  
  setwd(OutputDir)
  fwrite(cbind(GENES = row.names(TP_profile), TP_profile), paste0(PI, "_PD.txt"), quote = F, row.names = F, sep = '\t')
  fwrite(cbind(SpotID = gsub('\\-.*','-1',CellTypeEnrichment_subset_PD$cell), CellType = CellTypeEnrichment_subset_PD$geneSet), paste0(PI, '_CellTypeLabel_PD.txt'), quote = F, row.names = F, sep = '\t')
  
  #------------------------- GeoMx -------------------------
  MetaData = fread(paste0(OutputDir, 'MetaData_ST.txt'), header = T, stringsAsFactors = F)
  MetaData = as.data.frame(MetaData)
  
  GeoMx_TP = fread(paste0(OutputDir, 'GeoMx_readCount.txt'), header = T, stringsAsFactors = F)
  GeoMx_TP = as.data.frame(GeoMx_TP)
  
  row.names(GeoMx_TP) = make.names(GeoMx_TP$GeneSymbol,unique=T) 
  GeoMx_TP = GeoMx_TP[,-1]
  
  colnames(GeoMx_TP) = gsub('\\.', '-', gsub('.dcc', '', colnames(GeoMx_TP)))
  dim(GeoMx_TP)
  
  GeoMx_TP = GeoMx_TP[apply(GeoMx_TP, 1, function(x) !all(x==0)),]
  GeoMx_TP = GeoMx_TP[apply(GeoMx_TP, 2, function(x) !all(x==0)),]
  dim(GeoMx_TP)
  
  #---- resistant
  GeoMx_TP_resistant = GeoMx_TP[,which(colnames(GeoMx_TP) %in% MetaData$SampleID[which(MetaData$segment == "LB509-TH+")])]
  colnames(GeoMx_TP_resistant) = paste0('Resistant_', gsub('.*\\-', '', colnames(GeoMx_TP_resistant)))
  dim(GeoMx_TP_resistant)
  
  #---- vulnerable
  GeoMx_TP_vulnerable = GeoMx_TP[,which(colnames(GeoMx_TP) %in% MetaData$SampleID[which(MetaData$segment == 'LB509+TH+')])]
  colnames(GeoMx_TP_vulnerable) = paste0('Vulnerable_', gsub('.*\\-', '', colnames(GeoMx_TP_vulnerable)))
  dim(GeoMx_TP_vulnerable)
  
  identical(row.names(GeoMx_TP_vulnerable), row.names(GeoMx_TP_resistant))
  
  #---------------- scmap [projecting cells from a scRNA-seq experiment on to the cell-types or cells identified in a different experiment.]
  prediction_all = data.frame(row.names = c(colnames(GeoMx_TP_vulnerable), colnames(GeoMx_TP_resistant)), cell_type1=c(rep('Vulnerable', length(colnames(GeoMx_TP_vulnerable))), rep('Resistant', length(colnames(GeoMx_TP_resistant)))))
  
  sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(cbind(GeoMx_TP_vulnerable, GeoMx_TP_resistant))), colData = prediction_all)
  sce <- sce[!duplicated(rownames(sce)), ]
  counts(sce) = normcounts(sce)
  sce
  
  # use gene names as feature symbols
  rowData(sce)$feature_symbol <- rownames(sce)
  
  library(standR)
  sce_tmm <- geomxNorm(sce, method = "TMM")
  sce <- selectFeatures(sce_tmm, suppress_plot = FALSE)
  
  # index of a reference dataset 
  sce <- indexCluster(sce)  #writes the scmap_cluster_index item of the metadata slot of the reference dataset.
  
  #5.2Projection
  prediction_annotation = data.frame(row.names = colnames(TP_profile), cell_type1= rep('Neurons', dim(TP_profile)[2]) ) 
  sce_prediction = SingleCellExperiment(assays = list(normcounts = as.matrix(TP_profile)), colData = prediction_annotation)
  rowData(sce_prediction)$feature_symbol <- rownames(sce_prediction)
  logcounts(sce_prediction) <- normcounts(sce_prediction) + 1 # by using raw reads I get more RV neurons
  
  scmapCluster_results <- scmapCluster(
    projection = sce_prediction,  # 'SingleCellExperiment' object to project
    index_list = list(
      yan = metadata(sce)$scmap_cluster_index
    )
  )
  
  sce_prediction@colData$cell_type1 = scmapCluster_results$scmap_cluster_labs
  sce_prediction_assigned <- sce_prediction[, sce_prediction$cell_type1 %in% c("Resistant", "Vulnerable")]
  
  ## ----  Perform aggregation of counts and metadata by subject and cell type.
  
  Barcode_Donor_All_subset = Barcode_Donor_All[match(row.names(colData(sce_prediction_assigned)), Barcode_Donor_All$Barcode),]
  dim(Barcode_Donor_All_subset)
  identical(Barcode_Donor_All_subset$Barcode, row.names(colData(sce_prediction_assigned)))
  
  colData(sce_prediction_assigned)
  colData(sce_prediction_assigned)$subjectVar = Barcode_Donor_All_subset$Donor
  colData(sce_prediction_assigned)$cell_type1 
  
  MetaData_Integrated = data.frame(Barcode = Barcode_Donor_All_subset$Barcode, CellType = colData(sce_prediction_assigned)$cell_type1, Donor = Barcode_Donor_All_subset$Donor)
  
  library(dplyr)
  MetaData_Integrated %>%
    group_by(Donor, yan) %>%
    summarize(count=n())
  
  MetaData_Integrated$Donor_subtype[MetaData_Integrated$yan == 'Vulnerable'] = paste0(MetaData_Integrated$Donor[MetaData_Integrated$yan == 'Vulnerable'], 'V')
  MetaData_Integrated$Donor_subtype[MetaData_Integrated$yan == 'Resistant'] = paste0(MetaData_Integrated$Donor[MetaData_Integrated$yan == 'Resistant'], 'R')
  
  #---------------------------------- Save results
  setwd(OutputDir)
  fwrite(MetaData_Integrated, paste0(OutputDir, PI , '_Neurons_Subset_MetaData.txt'), quote = F, row.names = F, sep = '\t')
  fwrite(data.frame(GeneSymbol=row.names(normcounts(sce_prediction_assigned)), normcounts(sce_prediction_assigned)), paste0(OutputDir, PI , '_Neurons_Subset_GeneCount.txt'), quote = F, row.names = F, sep = '\t')
  #----------------------------------
  
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
