
#------------------------- single cell -------------------------

#------------------------- parameters
PI='GSE184950'
outputDir = paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis/', PI, '/')
dir.create(outputDir, recursive = T)

inputReadCount = paste0('/Users/isarnassiri/Documents/Parkinson/Data/',PI,'_RAW/')

#------------------------- MetaData
library(data.table)
MetaData = fread(paste0(inputReadCount, '/MetaData.txt'), stringsAsFactors = F, header = T)
MetaData = MetaData[which(MetaData$Region == 'SN' & MetaData$Diagnosis %in% c("PD")),]
MetaData = as.data.frame(MetaData)

table(MetaData$Diagnosis)

write.table(MetaData$Sample_ID, paste0(inputReadCount, 'Samples.txt'), quote = F, row.names = F, sep = '\t', col.names = F)

#------------------------- CellType Annotations
CellTypeEnrichment = fread(paste0('/Users/isarnassiri/Documents/Parkinson/Analysis/',PI,'/',PI,'_CellType_Enrichment/assignmentTable_dedup.txt'), stringsAsFactors = F, header = T)
CellTypeEnrichment = as.data.frame(CellTypeEnrichment)

CellTypeEnrichment_subset = CellTypeEnrichment[which(CellTypeEnrichment$Sample_ID %in% MetaData$Sample_ID),]
head(CellTypeEnrichment_subset)

dim(CellTypeEnrichment_subset[which(CellTypeEnrichment_subset$Diagnosis == 'normal'),])
dim(CellTypeEnrichment_subset[which(CellTypeEnrichment_subset$Diagnosis == 'PD'),])

table(is.na(CellTypeEnrichment_subset$cell_orig))

#------------------------- read the scRNAseq profile
setwd(inputReadCount)

library(data.table)
library(Seurat)

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

TP_profile = TP_profile[apply(TP_profile, 1, function(x) !all(x==0)),]
dim(TP_profile)
TP_profile = TP_profile[apply(TP_profile, 2, function(x) !all(x==0)),]
dim(TP_profile)

# TP_profile_normal = TP_profile_normal[apply(TP_profile_normal, 1, function(x) !all(x==0)),]
# dim(TP_profile_normal)
# TP_profile_normal = TP_profile_normal[apply(TP_profile_normal, 2, function(x) !all(x==0)),]
# dim(TP_profile_normal)

colnames(TP_profile)[1:10]
TP_profile[1:5, 1:5]

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

table(is.na(CellTypeEnrichment_subset_PD$cell_orig))
table(is.na(colnames(TP_profile)))

setdiff(gsub('\\-.*','-1',CellTypeEnrichment_subset_PD$cell), colnames(TP_profile) )
setdiff(colnames(TP_profile), gsub('\\-.*','-1',CellTypeEnrichment_subset_PD$cell) )

CellTypeEnrichment_subset_PD = CellTypeEnrichment_subset_PD[match(colnames(TP_profile), CellTypeEnrichment_subset_PD$cell_orig),]
identical(gsub('\\-.*','-1',CellTypeEnrichment_subset_PD$cell), colnames(TP_profile) )
dim(CellTypeEnrichment_subset_PD)
dim(TP_profile)

setwd(outputDir)
fwrite(cbind(GENES = row.names(TP_profile), TP_profile), paste0(PI, "_PD.txt"), quote = F, row.names = F, sep = '\t')
fwrite(cbind(SpotID = gsub('\\-.*','-1',CellTypeEnrichment_subset_PD$cell), CellType = CellTypeEnrichment_subset_PD$geneSet), paste0(PI, '_CellTypeLabel_PD.txt'), quote = F, row.names = F, sep = '\t')

table(CellTypeEnrichment_subset_PD$geneSet)

#--- save normal ---
# CellTypeEnrichment_subset_normal = CellTypeEnrichment_subset[which(CellTypeEnrichment_subset$cell_orig %in% colnames(TP_profile_normal)),]
# CellTypeEnrichment_subset_normal = CellTypeEnrichment_subset_normal[!duplicated(CellTypeEnrichment_subset$cell_orig),]
# CellTypeEnrichment_subset_normal = CellTypeEnrichment_subset_normal[!is.na(CellTypeEnrichment_subset_normal$cell_orig),]
# 
# dim(TP_profile_normal)
# if(length(setdiff(CellTypeEnrichment_subset$cell_orig, colnames(TP_profile_normal) )) > 0)
# {
#   TP_profile_normal = TP_profile_normal[,-which(colnames(TP_profile_normal) %in% setdiff(colnames(TP_profile_normal), CellTypeEnrichment_subset_normal$cell_orig) )]
# }
# 
# table(is.na(CellTypeEnrichment_subset_normal$cell_orig))
# table(is.na(colnames(TP_profile_normal)))
# 
# setdiff(gsub('\\-.*','-1',CellTypeEnrichment_subset_normal$cell), colnames(TP_profile_normal) )
# setdiff(colnames(TP_profile_normal), gsub('\\-.*','-1',CellTypeEnrichment_subset_normal$cell) )
# 
# CellTypeEnrichment_subset_normal = CellTypeEnrichment_subset_normal[match(colnames(TP_profile_normal), CellTypeEnrichment_subset_normal$cell_orig),]
# identical(gsub('\\-.*','-1',CellTypeEnrichment_subset_normal$cell), colnames(TP_profile_normal) )
# dim(CellTypeEnrichment_subset_normal)
# dim(TP_profile_normal)
# 
# setwd(outputDir)
# fwrite(cbind(GENES = row.names(TP_profile_normal), TP_profile_normal), paste0(PI, "_normal.txt"), quote = F, row.names = F, sep = '\t')
# fwrite(cbind(SpotID = gsub('\\-.*','-1',CellTypeEnrichment_subset_normal$cell), CellType = CellTypeEnrichment_subset_normal$geneSet), paste0(PI, '_CellTypeLabel_normal.txt'), quote = F, row.names = F, sep = '\t')
# 
# table(CellTypeEnrichment_subset_normal$geneSet)



#------------------------- GeoMx -------------------------
library(data.table)
# the following profiles were generated using Step0_GeoMx_Analysis.R
datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/")
MetaData = fread(paste0(datadir, 'MetaData.txt'), header = T, stringsAsFactors = F)
MetaData = as.data.frame(MetaData)

GeoMx_TP = fread(paste0(datadir, 'GeoMx_readCount.txt'), header = T, stringsAsFactors = F)
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
sum(is.na(GeoMx_TP_vulnerable))

#---------------- scmap [projecting cells from a scRNA-seq experiment on to the cell-types or cells identified in a different experiment.]
library(SingleCellExperiment)
library(scmap)

prediction_all = data.frame(row.names = c(colnames(GeoMx_TP_vulnerable), colnames(GeoMx_TP_resistant)), cell_type1=c(rep('Vulnerable', length(colnames(GeoMx_TP_vulnerable))), rep('Resistant', length(colnames(GeoMx_TP_resistant)))))

sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(cbind(GeoMx_TP_vulnerable, GeoMx_TP_resistant))), colData = prediction_all)
sce <- sce[!duplicated(rownames(sce)), ]
counts(sce) = normcounts(sce)
sce

# use gene names as feature symbols
rowData(sce)$feature_symbol <- rownames(sce)

library(standR)
sce_tmm <- geomxNorm(sce, method = "TMM")
# logcounts(sce) <- log2(normcounts(sce) + 1)

#select the most informative features (genes) from our input dataset:
sce <- selectFeatures(sce_tmm, suppress_plot = FALSE)
table(rowData(sce)$scmap_features)

# index of a reference dataset 
sce <- indexCluster(sce)  #writes the scmap_cluster_index item of the metadata slot of the reference dataset.
head(metadata(sce)$scmap_cluster_index)
dim(metadata(sce)$scmap_cluster_index)

getwd()
setwd(outputDir)
pdf('Selected_genes_training.pdf', width = 10, height = 10, useDingbats = F)
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))
dev.off()

# pdf('sce_prediction_clustering_method.pdf', width = 10, height = 10, useDingbats = F)
# 
# plot(
#   getSankey(
# colData(sce)$cell_type1, 
# scmapCluster_results$scmap_cluster_labs[,'yan'],
# plot_height = 400
#   )
# )
# dev.off()


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

grep('CERT1',rownames(sce_prediction))

table(scmapCluster_results$scmap_cluster_labs)
dim(sce_prediction@colData)

sce_prediction@colData$cell_type1 = scmapCluster_results$scmap_cluster_labs
sce_prediction_assigned <- sce_prediction[, sce_prediction$cell_type1 %in% c("Resistant", "Vulnerable")]

# pdf('sce_prediction_clustering_method.pdf', width = 10, height = 10, useDingbats = F)
# 
# plot(
#   getSankey(
# colData(sce)$cell_type1, 
# scmapCluster_results$scmap_cluster_labs[,'yan'],
# plot_height = 400
#   )
# )
# dev.off()
 
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
setwd(outputDir)
fwrite(MetaData_Integrated, paste0(outputDir, PI , '_Neurons_Subset_MetaData.txt'), quote = F, row.names = F, sep = '\t')
fwrite(data.frame(GeneSymbol=row.names(normcounts(sce_prediction_assigned)), normcounts(sce_prediction_assigned)), paste0(outputDir, PI , '_Neurons_Subset_GeneCount.txt'), quote = F, row.names = F, sep = '\t')
#----------------------------------

## ----
library(Seurat)
target_Data_Seurat <- CreateSeuratObject(counts = normcounts(sce_prediction_assigned), project = "TS_SC_Integration",  meta.data = MetaData_Integrated, assay = "TS_SC")

library(scater)
target_Data_sce <- as.SingleCellExperiment(target_Data_Seurat)

library(aggregateBioVar)
aggregate_counts = aggregateBioVar(scExp = target_Data_sce, cellVar = "yan", subjectVar = "Donor_subtype")

State = as.character(colData(aggregate_counts$AllCells)$Donor_subtype)
State[grep('V', State)] = 'Vulnerable'
State[grep('R', State)] = 'Resistant'
colData(aggregate_counts$AllCells)$State = State

## ----DESeq2Aggregate----------------------------------------------------------
library(DESeq2)
subj_dds_dataset <-
  DESeqDataSetFromMatrix(
countData = assay(aggregate_counts$AllCells, "counts"),
colData = colData(aggregate_counts$AllCells),
design = ~ State
  )

subj_dds <- DESeq(subj_dds_dataset)

resultsNames(subj_dds)

# using other methods does not change the results
# library("apeglm")
# res.ape <- lfcShrink(dds=subj_dds, coef=2, type="apeglm")
# dim(res.ape[which(res.ape$padj < 0.01),])
# 
# library(ashr)
# res.ash <- lfcShrink(dds=subj_dds, coef=2, type="ashr")
# dim(res.ash[which(res.ash$padj < 0.01),])
# 
# res.norm <- lfcShrink(dds=subj_dds, coef=2, type="normal")
# intersect(row.names(res.norm[which(res.norm$padj < 0.01),]), row.names(subj_dds_transf[which(subj_dds_transf$padj < 0.01),]))

subj_dds_results <-
  results( subj_dds, contrast = c("State", "Vulnerable", "Resistant") ) # State PD [Vulnerable] vs Control

subj_dds_transf <- as.data.frame(subj_dds_results) %>%
  bind_cols(feature = rownames(subj_dds_results)) %>%
  mutate(log_padj = -log(subj_dds_results$padj, base = 10))

table(subj_dds_transf$pvalue < 0.05)
table(subj_dds_transf$log2FoldChange[which(subj_dds_transf$pvalue < 0.05)] > 0) # up and down regulation

# Count_aggregated = as.data.frame(aggregate_counts$CellType@assays@data$counts)
# dim(Count_aggregated)
# Count_aggregated = data.frame(GeneSymbol = row.names(Count_aggregated), Count_aggregated)

#-------------------------------- pathway Enrichment

#--- gene enrichment analysis
library('org.Hs.eg.db')

# library("DOSE")
# edo <- enrichDGN(de[!is.na(de)])
# edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
# 
# library(enrichplot)
# barplot(edo, showCategory=20)
# dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
# cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

BD = 'c5.all.v2023.1.Hs.entrez.gmt'

library(qusage)
library(clusterProfiler)
gmtfile <- paste0("/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Scripts_USED_visium/msigdb_v2023_Hs/", BD)
c5 <- read.gmt(gmtfile)

#-------------------------------- pathway Enrichment

#--- gene enrichment analysis
library('org.Hs.eg.db')

# library("DOSE")
# edo <- enrichDGN(de[!is.na(de)])
# edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
# 
# library(enrichplot)
# barplot(edo, showCategory=20)
# dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
# cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

BD = 'c5.all.v2023.1.Hs.entrez.gmt'

library(qusage)
library(clusterProfiler)
gmtfile <- paste0("/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Scripts_USED_visium/msigdb_v2023_Hs/", BD)
c5 <- read.gmt(gmtfile)

# use mapIds method to obtain Entrez IDs
input = subj_dds_transf$feature[which(subj_dds_transf$pvalue < 0.05)]
de = as.character(mapIds(org.Hs.eg.db, input, 'ENTREZID', 'SYMBOL'))

egmt <- enricher(de[!is.na(de)], TERM2GENE=c5)
 
library(cowplot)
library(ggplot2)
p2 <- dotplot(egmt, showCategory=5) + ggtitle("")
p2

paste0(outputDir, Profile, '_DE.txt')
pdf(paste0('PathwayEnrichment.pdf') ,width = 10, height = 10,useDingbats = FALSE)
print(plot_grid(p2, ncol=1))
dev.off()

#---- pathview - visualization of hsa00190
library('org.Hs.eg.db')
keytypes(org.Hs.eg.db)
names(input) = as.character(mapIds(org.Hs.eg.db, input, 'ENTREZID', 'SYMBOL'))

getwd()
library("pathview")
###################################################
### code chunk number 15: kegg.native
###################################################
pv.out <- pathview(gene.data = names(input), pathway.id = '05012',
   species = "hsa", kegg.native = T,
   same.layer = F, high=list(gene="gold"), out.suffix = 'VR')
head(pv.out$plot.data.gene)
 
intersect(as.character(input), pv.out$plot.data.gene$labels)

library(dplyr)
library(annotables)
devtools::install_github("stephenturner/annotables")
grch38_tx2gene
grch38
# convert to the tibble
grch38 = tbl_df(grch38)

Annotated_Genes = grch38 %>% 
dplyr::filter(symbol %in% as.character(subj_dds_transf$feature[which(subj_dds_transf$pvalue < 0.05)])) %>% 
dplyr::select(symbol, description)

subj_dds_transf = subj_dds_transf[order(subj_dds_transf$pvalue, decreasing = F),]
subj_dds_transf$feature[which(subj_dds_transf$pvalue < 0.05)]

Annotated_Genes = Annotated_Genes[match(subj_dds_transf$feature[which(subj_dds_transf$pvalue < 0.05)], Annotated_Genes$symbol),]






#-------------------------------- visualization [test]
library(hgnc)
library(data.table)
library(ComplexHeatmap)

hgnc_tbl <- import_hgnc_dataset()

Controls = c("C1", "F1")
PD = c("E4", "D2")

outputDir = '/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis/'
CellTypes = c("Astrocytes", "Microglia", "Neurons", "Oligodendrocyte_Progenitor_Cells", "Oligodendrocytes", "Pericytes")
i=1

for(i in 1:length(CellTypes))
{
print(CellTypes[i])

setwd('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis')
Files = list.files(pattern = CellTypes[i])

if(length(grep('_used', Files))!=0)
{
  Files = Files[-grep('_used', Files)]
}

options(datatable.fread.datatable=FALSE)
Dorsal_Count = fread(Files[grep(paste0('^Dorsal_',CellTypes[i],'_Count'), Files)], stringsAsFactors = F, header = T)
Dorsal_DE = fread(Files[grep(paste0('^Dorsal_',CellTypes[i],'_DE'), Files)], stringsAsFactors = F, header = T)

D=it=1
for(D in 1:dim(Dorsal_DE)[1])
{
  DE = Dorsal_DE[D,]
  Expression = Dorsal_Count[which(Dorsal_Count$GeneSymbol == DE$feature),]
  
  if(DE$log2FoldChange > 0 & !is.na(DE$log2FoldChange))
  {
if(Expression$D2 > Expression$C1 & Expression$D2 > Expression$F1 & Expression$E4 > Expression$C1 & Expression$E4 > Expression$F1)
{
  if(it==1){Dorsal_DE_used = DE; it=it+1}else{Dorsal_DE_used = rbind(Dorsal_DE_used, DE)}
}
  }
  
  if(DE$log2FoldChange < 0 & !is.na(DE$log2FoldChange))
  {
if(Expression$D2 < Expression$C1 & Expression$D2 < Expression$F1 & Expression$E4 < Expression$C1 & Expression$E4 < Expression$F1)
{
  if(it==1){Dorsal_DE_used = DE; it=it+1}else{Dorsal_DE_used = rbind(Dorsal_DE_used, DE)}
}
  }
  
}

Ventral_Count = fread(Files[grep(paste0('^Ventral_',CellTypes[i],'_Count'), Files)], stringsAsFactors = F, header = T)
Ventral_DE = fread(Files[grep(paste0('^Ventral_',CellTypes[i],'_DE'), Files)], stringsAsFactors = F, header = T)

D=it=1
for(D in 1:dim(Ventral_DE)[1])
{
  DE = Ventral_DE[D,]
  Expression = Ventral_Count[which(Ventral_Count$GeneSymbol == DE$feature),]
  
  if(DE$log2FoldChange > 0 & !is.na(DE$log2FoldChange))
  {
if(Expression$D2 > Expression$C1 & Expression$D2 > Expression$F1 & Expression$E4 > Expression$C1 & Expression$E4 > Expression$F1)
{
  if(it==1){Ventral_DE_used = DE; it=it+1}else{Ventral_DE_used = rbind(Ventral_DE_used, DE)}
}
  }
  
  if(DE$log2FoldChange < 0 & !is.na(DE$log2FoldChange))
  {
if(Expression$D2 < Expression$C1 & Expression$D2 < Expression$F1 & Expression$E4 < Expression$C1 & Expression$E4 < Expression$F1)
{
  if(it==1){Ventral_DE_used = DE; it=it+1}else{Ventral_DE_used = rbind(Ventral_DE_used, DE)}
}
  }
}

Dorsal_Ventral_Count = fread(Files[grep(paste0('Dorsal_Ventral_',CellTypes[i],'_Count'), Files)], stringsAsFactors = F, header = T)
Dorsal_Ventral_DE = fread(Files[grep(paste0('Dorsal_Ventral_',CellTypes[i],'_DE'), Files)], stringsAsFactors = F, header = T)

D=it=1
for(D in 1:dim(Dorsal_Ventral_DE)[1])
{
  DE = Dorsal_Ventral_DE[D,]
  Expression = Dorsal_Ventral_Count[which(Dorsal_Ventral_Count$GeneSymbol == DE$feature),]
  
  if(DE$log2FoldChange > 0 & !is.na(DE$log2FoldChange))
  {
if(Expression$D2 > Expression$C1 & Expression$D2 > Expression$F1 & Expression$E4 > Expression$C1 & Expression$E4 > Expression$F1)
{
  if(it==1){Dorsal_Ventral_DE_used = DE; it=it+1}else{Dorsal_Ventral_DE_used = rbind(Dorsal_Ventral_DE_used, DE)}
}
  }
  
  if(DE$log2FoldChange < 0 & !is.na(DE$log2FoldChange))
  {
if(Expression$D2 < Expression$C1 & Expression$D2 < Expression$F1 & Expression$E4 < Expression$C1 & Expression$E4 < Expression$F1)
{
  if(it==1){Dorsal_Ventral_DE_used = DE; it=it+1}else{Dorsal_Ventral_DE_used = rbind(Dorsal_Ventral_DE_used, DE)}
}
  }
}

#------ DE
DE <- list(Dorsal = Dorsal_DE_used$feature[which(Dorsal_DE_used$padj < 1e-2)],
   Ventral = Ventral_DE_used$feature[which(Ventral_DE_used$padj < 1e-2)],
   Dorsal_Ventral = Dorsal_Ventral_DE_used$feature[which(Dorsal_Ventral_DE_used$padj < 1e-4)])
setdiffDE = lapply(1:length(DE), function(n) setdiff(DE[[n]], unlist(DE[-n])))
str(setdiffDE)

names(setdiffDE) = c('Dorsal', 'Ventral', 'Dorsal_Ventral')

Dorsal_DE_used = Dorsal_DE_used[order(Dorsal_DE_used$padj, decreasing = F),]
Dorsal_DE_used_subset = Dorsal_DE_used[which(Dorsal_DE_used$feature %in% setdiffDE[['Dorsal']]), ]

Ventral_DE_used = Ventral_DE_used[order(Ventral_DE_used$padj, decreasing = F),]
Ventral_DE_used_subset = Ventral_DE_used[which(Ventral_DE_used$feature %in% setdiffDE[['Ventral']]), ]

Dorsal_Ventral_DE_used = Dorsal_Ventral_DE_used[order(Dorsal_Ventral_DE_used$padj, decreasing = F),]
Dorsal_Ventral_DE_used_subset = Dorsal_Ventral_DE_used[which(Dorsal_Ventral_DE_used$feature %in% setdiffDE[['Dorsal_Ventral']]), ]

#------ Count
Dorsal_Count_subet = Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% setdiffDE[['Dorsal']]),]
Ventral_Count_subet = Ventral_Count[which(Ventral_Count$GeneSymbol %in% setdiffDE[['Ventral']]),]
Dorsal_Ventral_Count_subet = Dorsal_Ventral_Count[which(Dorsal_Ventral_Count$GeneSymbol %in% setdiffDE[['Dorsal_Ventral']]),]

#------ Count subset
CountMatrix = rbindlist(list(Dorsal_Count_subet[match(Dorsal_DE_used_subset$feature[1:10], Dorsal_Count_subet$GeneSymbol),], Ventral_Count_subet[match(Ventral_DE_used_subset$feature[1:10], Ventral_Count_subet$GeneSymbol),], Dorsal_Ventral_Count_subet[match(Dorsal_Ventral_DE_used_subset$feature[1:10], Dorsal_Ventral_Count_subet$GeneSymbol),]))
CountMatrix = as.data.frame(CountMatrix)
rownames(CountMatrix) = CountMatrix$GeneSymbol
CountMatrix = CountMatrix[,-1]
dim(CountMatrix)

# --The differentiation of other genes into colors is inhibited by genes that are highly expressed.
CountMatrixlog10 = log10(CountMatrix+1)
#------ 

#------ annotations
State = as.character(colnames(CountMatrix))
State[which(State %in% Controls)] = 'Control'
State[which(State %in% PD)] = 'PD'

annotation_col = data.frame(
  State = factor(State, levels = c('Control', 'PD'))
)
rownames(annotation_col) = as.character(colnames(CountMatrix))
str(annotation_col)

annotation_row = data.frame(
  GeneClass = factor(rep(c('Dorsal', 'Ventral', 'Dorsal_Ventral'), c(10, 10, 10)))
)
rownames(annotation_row) = c(Dorsal_DE_used_subset$feature[1:10], Ventral_DE_used_subset$feature[1:10], Dorsal_Ventral_DE_used_subset$feature[1:10])
identical(row.names(CountMatrix), row.names(annotation_row))

ann_colors = list(
  State = c(Control = "#1B9E77", PD = "#D95F02"),
  GeneClass = c(Dorsal = "#7570B3", Ventral = "#E7298A", Dorsal_Ventral = "#66A61E")
)

p <- pheatmap(CountMatrix, annotation_col = annotation_col, annotation_row = annotation_row, 
  annotation_colors = ann_colors, 
  row_split = annotation_row$GeneClass,
  column_split = annotation_col$State, name = "Expression", scale='row',
  display_numbers = as.matrix(CountMatrix), fontsize_number=13  )
# ref: https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
p

pdf(paste0('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Figures/', CellTypes[i], '_pheatmap.pdf'), width = 8, height = 10,useDingbats = FALSE)
print(p)
dev.off()

Dorsal_DE_used_Annotated = data.frame(symbol = Dorsal_DE_used_subset$feature[1:10])
Dorsal_Ventral_DE_used_Annotated = data.frame(symbol =Dorsal_Ventral_DE_used_subset$feature[1:10])
Ventral_DE_used_Annotated = data.frame(symbol =Ventral_DE_used_subset$feature[1:10])

Dorsal_DE_used_Annotated = merge(Dorsal_DE_used_Annotated, hgnc_tbl[,c('symbol', 'name', 'gene_group')])
Dorsal_Ventral_DE_used_Annotated = merge(Dorsal_Ventral_DE_used_Annotated, hgnc_tbl[,c('symbol', 'name', 'gene_group')])
Ventral_DE_used_Annotated = merge(Ventral_DE_used_Annotated, hgnc_tbl[,c('symbol', 'name', 'gene_group')])

DE_Annotated = rbindlist(list(Dorsal_DE_used_Annotated,Dorsal_Ventral_DE_used_Annotated,Ventral_DE_used_Annotated))

fwrite(Dorsal_Ventral_DE_used, paste0(outputDir, Profile, '_',CellType, '_DE_used.txt'), quote = F, row.names = F, sep = '\t')
fwrite(Dorsal_DE_used, paste0(outputDir, Profile, '_',CellType, '_DE.txt'), quote = F, row.names = F, sep = '\t')
fwrite(Ventral_DE_used, paste0(outputDir, Profile, '_',CellType, '_DE.txt'), quote = F, row.names = F, sep = '\t')

fwrite(DE_Annotated, paste0('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Figures/DE_Annotated_', CellTypes[i], '.txt'), quote = F, row.names = F, sep = '\t')
}

#---------------- compare with published articles

outputDir = '/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis/'
setwd(outputDir)
Files = list.files(pattern = "Neurons")

library(data.table)
Dorsal_Ventral_Count = fread(Files[grep(paste0('Dorsal_Ventral_Neurons_Count'), Files)], stringsAsFactors = F, header = T)
Dorsal_Ventral_DE = fread(Files[grep(paste0('Dorsal_Ventral_Neurons_DE'), Files)], stringsAsFactors = F, header = T)

Input_Comparison = fread('Input_Comparison.txt', stringsAsFactors = F, header = F)
colnames(Input_Comparison) = c('Rank', 'GeneSymbol')

length(which(Input_Comparison$V2 %in% Dorsal_Ventral_Count$GeneSymbol))

merged = merge(Input_Comparison, Dorsal_Ventral_Count, by = 'GeneSymbol')

#---------------- compare with published articles

outputDir = '/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis/'
setwd(outputDir)
Files = list.files(pattern = "Neurons")

library(data.table)
Dorsal_Ventral_Count = fread(Files[grep(paste0('Dorsal_Ventral_Neurons_Count'), Files)], stringsAsFactors = F, header = T)

Dorsal_Count = fread('Dorsal_Neurons_Count_aggregated.txt', stringsAsFactors = F, header = T)
Ventral_Count = fread('Ventral_Neurons_Count_aggregated.txt', stringsAsFactors = F, header = T)

c("GAPDH",
  "B2M",
  "TH",
  "DAT",
  "VMAT2",
  "GFAP",
  "PMP22",
  "GAD1",
  "PCP4",
  "RAB3B",
  "HCN1",
  "SPARC",
  "SNX8",
  "MT1G",
  "ANXA1",
  "ATP13A4",
  "LYPD1")


df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('C1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('C1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("PCP4", "HCN1", "RAB3B")

df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("FJX1", "LYPD1", "CALB1")),c('F1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("FJX1", "LYPD1", "CALB1")),c('F1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("FJX1", "LYPD1", "CALB1")

df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("SNCA", "LIX1", "GRIK1")),c('C1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("SNCA", "LIX1", "GRIK1")),c('C1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("SNCA", "LIX1", "GRIK1")


df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('F1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('F1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("PCP4", "HCN1", "RAB3B")



# use mapIds method to obtain Entrez IDs
input = subj_dds_transf$feature[which(subj_dds_transf$padj < 0.05)]
de = as.character(mapIds(org.Hs.eg.db, input, 'ENTREZID', 'SYMBOL'))

egmt <- enricher(de[!is.na(de)], TERM2GENE=c5)

egmt@result[grep('GLIOGENESIS', egmt@result$Description),'p.adjust']

egmt@result$Description = gsub('^[^_]+','',egmt@result$Description)

library(cowplot)
library(ggplot2)
p2 <- dotplot(egmt, showCategory=20) + ggtitle("")
p2

paste0(outputDir, Profile, '_DE.txt')
pdf(paste0('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Figures/RESULTS_IntegrationAnalysis_', Profile,'_', CellType, '.pdf') ,width = 10, height = 10,useDingbats = FALSE)
print(plot_grid(p2, ncol=1))
dev.off()

#---- pathview - visualization of hsa00190
library('org.Hs.eg.db')
keytypes(org.Hs.eg.db)
names(input) = as.character(mapIds(org.Hs.eg.db, input, 'ENTREZID', 'SYMBOL'))

setwd('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Figures/')
library("pathview")
###################################################
### code chunk number 15: kegg.native
###################################################
pv.out <- pathview(gene.data = names(input), pathway.id = '05012',
   species = "hsa", kegg.native = T,
   same.layer = F, high=list(gene="gold"), out.suffix = paste0(Profile, '_', CellType))
head(pv.out$plot.data.gene)
}
}

#-------------------------------- visualization
library(hgnc)
library(data.table)
library(ComplexHeatmap)

hgnc_tbl <- import_hgnc_dataset()

Controls = c("C1", "F1")
PD = c("E4", "D2")

outputDir = '/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis/'
CellTypes = c("Astrocytes", "Microglia", "Neurons", "Oligodendrocyte_Progenitor_Cells", "Oligodendrocytes", "Pericytes")
i=1

for(i in 1:length(CellTypes))
{
  print(CellTypes[i])
  
  setwd('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis')
  Files = list.files(pattern = CellTypes[i])
  
  if(length(grep('_used', Files))!=0)
  {
Files = Files[-grep('_used', Files)]
  }
  
  options(datatable.fread.datatable=FALSE)
  Dorsal_Count = fread(Files[grep(paste0('^Dorsal_',CellTypes[i],'_Count'), Files)], stringsAsFactors = F, header = T)
  Dorsal_DE = fread(Files[grep(paste0('^Dorsal_',CellTypes[i],'_DE'), Files)], stringsAsFactors = F, header = T)
  
  D=it=1
  for(D in 1:dim(Dorsal_DE)[1])
  {
DE = Dorsal_DE[D,]
Expression = Dorsal_Count[which(Dorsal_Count$GeneSymbol == DE$feature),]

if(DE$log2FoldChange > 0 & !is.na(DE$log2FoldChange))
{
  if(Expression$D2 > Expression$C1 & Expression$D2 > Expression$F1 & Expression$E4 > Expression$C1 & Expression$E4 > Expression$F1)
  {
if(it==1){Dorsal_DE_used = DE; it=it+1}else{Dorsal_DE_used = rbind(Dorsal_DE_used, DE)}
  }
}

if(DE$log2FoldChange < 0 & !is.na(DE$log2FoldChange))
{
  if(Expression$D2 < Expression$C1 & Expression$D2 < Expression$F1 & Expression$E4 < Expression$C1 & Expression$E4 < Expression$F1)
  {
if(it==1){Dorsal_DE_used = DE; it=it+1}else{Dorsal_DE_used = rbind(Dorsal_DE_used, DE)}
  }
}

  }
  
  Ventral_Count = fread(Files[grep(paste0('^Ventral_',CellTypes[i],'_Count'), Files)], stringsAsFactors = F, header = T)
  Ventral_DE = fread(Files[grep(paste0('^Ventral_',CellTypes[i],'_DE'), Files)], stringsAsFactors = F, header = T)
  
  D=it=1
  for(D in 1:dim(Ventral_DE)[1])
  {
DE = Ventral_DE[D,]
Expression = Ventral_Count[which(Ventral_Count$GeneSymbol == DE$feature),]

if(DE$log2FoldChange > 0 & !is.na(DE$log2FoldChange))
{
  if(Expression$D2 > Expression$C1 & Expression$D2 > Expression$F1 & Expression$E4 > Expression$C1 & Expression$E4 > Expression$F1)
  {
if(it==1){Ventral_DE_used = DE; it=it+1}else{Ventral_DE_used = rbind(Ventral_DE_used, DE)}
  }
}

if(DE$log2FoldChange < 0 & !is.na(DE$log2FoldChange))
{
  if(Expression$D2 < Expression$C1 & Expression$D2 < Expression$F1 & Expression$E4 < Expression$C1 & Expression$E4 < Expression$F1)
  {
if(it==1){Ventral_DE_used = DE; it=it+1}else{Ventral_DE_used = rbind(Ventral_DE_used, DE)}
  }
}
  }
  
  Dorsal_Ventral_Count = fread(Files[grep(paste0('Dorsal_Ventral_',CellTypes[i],'_Count'), Files)], stringsAsFactors = F, header = T)
  Dorsal_Ventral_DE = fread(Files[grep(paste0('Dorsal_Ventral_',CellTypes[i],'_DE'), Files)], stringsAsFactors = F, header = T)
  
  D=it=1
  for(D in 1:dim(Dorsal_Ventral_DE)[1])
  {
DE = Dorsal_Ventral_DE[D,]
Expression = Dorsal_Ventral_Count[which(Dorsal_Ventral_Count$GeneSymbol == DE$feature),]

if(DE$log2FoldChange > 0 & !is.na(DE$log2FoldChange))
{
  if(Expression$D2 > Expression$C1 & Expression$D2 > Expression$F1 & Expression$E4 > Expression$C1 & Expression$E4 > Expression$F1)
  {
if(it==1){Dorsal_Ventral_DE_used = DE; it=it+1}else{Dorsal_Ventral_DE_used = rbind(Dorsal_Ventral_DE_used, DE)}
  }
}

if(DE$log2FoldChange < 0 & !is.na(DE$log2FoldChange))
{
  if(Expression$D2 < Expression$C1 & Expression$D2 < Expression$F1 & Expression$E4 < Expression$C1 & Expression$E4 < Expression$F1)
  {
if(it==1){Dorsal_Ventral_DE_used = DE; it=it+1}else{Dorsal_Ventral_DE_used = rbind(Dorsal_Ventral_DE_used, DE)}
  }
}
  }
  
  #------ DE
  DE <- list(Dorsal = Dorsal_DE_used$feature[which(Dorsal_DE_used$padj < 1e-2)],
 Ventral = Ventral_DE_used$feature[which(Ventral_DE_used$padj < 1e-2)],
 Dorsal_Ventral = Dorsal_Ventral_DE_used$feature[which(Dorsal_Ventral_DE_used$padj < 1e-4)])
  setdiffDE = lapply(1:length(DE), function(n) setdiff(DE[[n]], unlist(DE[-n])))
  str(setdiffDE)
  
  names(setdiffDE) = c('Dorsal', 'Ventral', 'Dorsal_Ventral')
  
  Dorsal_DE_used = Dorsal_DE_used[order(Dorsal_DE_used$padj, decreasing = F),]
  Dorsal_DE_used_subset = Dorsal_DE_used[which(Dorsal_DE_used$feature %in% setdiffDE[['Dorsal']]), ]
  
  Ventral_DE_used = Ventral_DE_used[order(Ventral_DE_used$padj, decreasing = F),]
  Ventral_DE_used_subset = Ventral_DE_used[which(Ventral_DE_used$feature %in% setdiffDE[['Ventral']]), ]
  
  Dorsal_Ventral_DE_used = Dorsal_Ventral_DE_used[order(Dorsal_Ventral_DE_used$padj, decreasing = F),]
  Dorsal_Ventral_DE_used_subset = Dorsal_Ventral_DE_used[which(Dorsal_Ventral_DE_used$feature %in% setdiffDE[['Dorsal_Ventral']]), ]
  
  #------ Count
  Dorsal_Count_subet = Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% setdiffDE[['Dorsal']]),]
  Ventral_Count_subet = Ventral_Count[which(Ventral_Count$GeneSymbol %in% setdiffDE[['Ventral']]),]
  Dorsal_Ventral_Count_subet = Dorsal_Ventral_Count[which(Dorsal_Ventral_Count$GeneSymbol %in% setdiffDE[['Dorsal_Ventral']]),]
  
  #------ Count subset
  CountMatrix = rbindlist(list(Dorsal_Count_subet[match(Dorsal_DE_used_subset$feature[1:10], Dorsal_Count_subet$GeneSymbol),], Ventral_Count_subet[match(Ventral_DE_used_subset$feature[1:10], Ventral_Count_subet$GeneSymbol),], Dorsal_Ventral_Count_subet[match(Dorsal_Ventral_DE_used_subset$feature[1:10], Dorsal_Ventral_Count_subet$GeneSymbol),]))
  CountMatrix = as.data.frame(CountMatrix)
  rownames(CountMatrix) = CountMatrix$GeneSymbol
  CountMatrix = CountMatrix[,-1]
  dim(CountMatrix)
  
  # --The differentiation of other genes into colors is inhibited by genes that are highly expressed.
  CountMatrixlog10 = log10(CountMatrix+1)
  #------ 
  
  #------ annotations
  State = as.character(colnames(CountMatrix))
  State[which(State %in% Controls)] = 'Control'
  State[which(State %in% PD)] = 'PD'
  
  annotation_col = data.frame(
State = factor(State, levels = c('Control', 'PD'))
  )
  rownames(annotation_col) = as.character(colnames(CountMatrix))
  str(annotation_col)
  
  annotation_row = data.frame(
GeneClass = factor(rep(c('Dorsal', 'Ventral', 'Dorsal_Ventral'), c(10, 10, 10)))
  )
  rownames(annotation_row) = c(Dorsal_DE_used_subset$feature[1:10], Ventral_DE_used_subset$feature[1:10], Dorsal_Ventral_DE_used_subset$feature[1:10])
  identical(row.names(CountMatrix), row.names(annotation_row))
  
  ann_colors = list(
State = c(Control = "#1B9E77", PD = "#D95F02"),
GeneClass = c(Dorsal = "#7570B3", Ventral = "#E7298A", Dorsal_Ventral = "#66A61E")
  )
  
  p <- pheatmap(CountMatrix, annotation_col = annotation_col, annotation_row = annotation_row, 
annotation_colors = ann_colors, 
row_split = annotation_row$GeneClass,
column_split = annotation_col$State, name = "Expression", scale='row',
display_numbers = as.matrix(CountMatrix), fontsize_number=13  )
  # ref: https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
  p
  
  pdf(paste0('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Figures/', CellTypes[i], '_pheatmap.pdf'), width = 8, height = 10,useDingbats = FALSE)
  print(p)
  dev.off()
  
  Dorsal_DE_used_Annotated = data.frame(symbol = Dorsal_DE_used_subset$feature[1:10])
  Dorsal_Ventral_DE_used_Annotated = data.frame(symbol =Dorsal_Ventral_DE_used_subset$feature[1:10])
  Ventral_DE_used_Annotated = data.frame(symbol =Ventral_DE_used_subset$feature[1:10])
  
  Dorsal_DE_used_Annotated = merge(Dorsal_DE_used_Annotated, hgnc_tbl[,c('symbol', 'name', 'gene_group')])
  Dorsal_Ventral_DE_used_Annotated = merge(Dorsal_Ventral_DE_used_Annotated, hgnc_tbl[,c('symbol', 'name', 'gene_group')])
  Ventral_DE_used_Annotated = merge(Ventral_DE_used_Annotated, hgnc_tbl[,c('symbol', 'name', 'gene_group')])
  
  DE_Annotated = rbindlist(list(Dorsal_DE_used_Annotated,Dorsal_Ventral_DE_used_Annotated,Ventral_DE_used_Annotated))
  
  fwrite(Dorsal_Ventral_DE_used, paste0(outputDir, Profile, '_',CellType, '_DE_used.txt'), quote = F, row.names = F, sep = '\t')
  fwrite(Dorsal_DE_used, paste0(outputDir, Profile, '_',CellType, '_DE.txt'), quote = F, row.names = F, sep = '\t')
  fwrite(Ventral_DE_used, paste0(outputDir, Profile, '_',CellType, '_DE.txt'), quote = F, row.names = F, sep = '\t')
  
  fwrite(DE_Annotated, paste0('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Figures/DE_Annotated_', CellTypes[i], '.txt'), quote = F, row.names = F, sep = '\t')
}

#---------------- compare with published articles

outputDir = '/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis/'
setwd(outputDir)
Files = list.files(pattern = "Neurons")

library(data.table)
Dorsal_Ventral_Count = fread(Files[grep(paste0('Dorsal_Ventral_Neurons_Count'), Files)], stringsAsFactors = F, header = T)
Dorsal_Ventral_DE = fread(Files[grep(paste0('Dorsal_Ventral_Neurons_DE'), Files)], stringsAsFactors = F, header = T)

Input_Comparison = fread('Input_Comparison.txt', stringsAsFactors = F, header = F)
colnames(Input_Comparison) = c('Rank', 'GeneSymbol')

length(which(Input_Comparison$V2 %in% Dorsal_Ventral_Count$GeneSymbol))

merged = merge(Input_Comparison, Dorsal_Ventral_Count, by = 'GeneSymbol')

#---------------- compare with published articles

outputDir = '/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Analysis/SVG/Inputs/INPUTS_for_SVG_Analysis/'
setwd(outputDir)
Files = list.files(pattern = "Neurons")

library(data.table)
Dorsal_Ventral_Count = fread(Files[grep(paste0('Dorsal_Ventral_Neurons_Count'), Files)], stringsAsFactors = F, header = T)

Dorsal_Count = fread('Dorsal_Neurons_Count_aggregated.txt', stringsAsFactors = F, header = T)
Ventral_Count = fread('Ventral_Neurons_Count_aggregated.txt', stringsAsFactors = F, header = T)

c("GAPDH",
  "B2M",
  "TH",
  "DAT",
  "VMAT2",
  "GFAP",
  "PMP22",
  "GAD1",
  "PCP4",
  "RAB3B",
  "HCN1",
  "SPARC",
  "SNX8",
  "MT1G",
  "ANXA1",
  "ATP13A4",
  "LYPD1")


df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('C1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('C1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("PCP4", "HCN1", "RAB3B")

df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("FJX1", "LYPD1", "CALB1")),c('F1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("FJX1", "LYPD1", "CALB1")),c('F1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("FJX1", "LYPD1", "CALB1")

df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("SNCA", "LIX1", "GRIK1")),c('C1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("SNCA", "LIX1", "GRIK1")),c('C1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("SNCA", "LIX1", "GRIK1")


df = data.frame(Dorsal_Count[which(Dorsal_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('F1')], Ventral_Count[which(Ventral_Count$GeneSymbol %in% c("PCP4", "HCN1", "RAB3B")),c('F1')])

colnames(df) = c('Dorsal', 'Ventral')
row.names(df) = c("PCP4", "HCN1", "RAB3B")



