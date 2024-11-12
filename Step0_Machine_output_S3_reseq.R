
#============================================ Repeat the analysis using new parameters

#============================================
library(GeomxTools)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(tidyverse)

## Introduction
# This summarises the QC of the GeoMx pilot run of 10th October, utilising the 'GeoMxTools' pipeline.

#============================================ Collate files
datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/")
setwd(datadir)

DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/pkcs"), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)

# I added NucleiCount manually to the following file from Annie Initial Dataset.xlsx file in the 'DSPDA Files' folder
SampleAnnotationFile <-
  dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)

library(readxl)

#============================================ load data
Data <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = excel_sheets(path = SampleAnnotationFile),
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi","roi","state", "NucleiCount"),
                         experimentDataColNames = c("panel"))

#============================================

library(knitr)
pkcs <- annotation(Data)
modules <- gsub(".pkc", "", pkcs)

# Add one to all counts in an expression matrix
Data <- shiftCountsOne(Data, useDALogic = TRUE)

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 10000,     # for WGT is 10K. Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)

Data <- setSegmentQCFlags(Data, qcCutoffs = QC_params)        

#============================================ Collate QC Results
QCResults <- protocolData(Data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))

QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

col_by <- "segment"

#============================================ Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

colnames(sData(Data))

QC_histogram(sData(Data), "NucleiCount", col_by, 20)

QC_histogram(sData(Data), "Trimmed (%)", col_by, 80)

# Stitching: These fragments are computationally merged into a single, longer read using specialized algorithms. This is essential because short reads alone might not uniquely map to the target gene sequences.
QC_histogram(sData(Data), "Stitched (%)", col_by, 80)

QC_histogram(sData(Data), "Aligned (%)", col_by, 75)

QC_histogram(sData(Data), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")

QC_histogram(sData(Data), "area", col_by, 1000, scale_trans = "log10")

## Generate negative probe count
# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(Data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(Data)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData & plot
negCols <- paste0("NegGeoMean_", modules)
pData(Data)[, negCols] <- sData(Data)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(Data), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(Data) <- pData(Data)[, !colnames(pData(Data)) %in% negCols]

#============================================ QC as tables
# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(Data)$NTC), col.names = c("NTC Count", "# of Segments"))

## ----QCSummaryTable, results = "asis"-----------------------------------------
kable(QC_Summary, caption = "QC Summary Table for each Segment")
#============================================

SamplesMetadata = read_excel('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/annotation/Annotation_File.xlsx')

#============================================

sData_df <- sData(Data)[, dimnames(sData(Data))[[2]]]

library(tidyr)
new_sData_df <- unnest(sData_df, QCFlags)
new_sData_df <- unnest(new_sData_df, colnames(new_sData_df)[grep('%',colnames(new_sData_df))])
new_sData_df <- unnest(new_sData_df, NegGeoMean)

colnames(new_sData_df)

library(stringr)
string = str_split_fixed(SamplesMetadata$Sample_ID, "-", 4)
SamplesMetadata$batch = string[,3]

# LowSaturation
kable(table(State = SamplesMetadata$segment[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowSaturation])]))
kable(table(Batch = SamplesMetadata$batch[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowSaturation])]))

# LowAligned
kable(table(State = SamplesMetadata$segment[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowAligned])]))
kable(table(Batch = SamplesMetadata$batch[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowAligned])]))

# LowReads
kable(table(State = SamplesMetadata$segment[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowReads])]))
kable(table(Batch = SamplesMetadata$batch[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowReads])]))

# LowTrimmed
kable(table(State = SamplesMetadata$segment[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowTrimmed])]))
kable(table(Batch = SamplesMetadata$batch[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[new_sData_df$LowTrimmed])]))


#============================================ summary of samples
# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols

df_pData = pData(Data)
df_pData$Donor = gsub(' .*', '', df_pData$`slide name`)

table(df_pData$Donor)
table(df_pData$segment)

#============================================

#============================================ save QC results

sData_df <- sData(Data)[, dimnames(sData(Data))[[2]]]
dim(sData_df)

colnames(sData_df)
summary(sData_df$NucleiCount[sData_df$segment == 'control'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-with-LBP'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-without-LBP'])


library(tidyr)
new_sData_df <- unnest(sData_df, QCFlags)
new_sData_df <- unnest(new_sData_df, colnames(new_sData_df)[grep('%',colnames(new_sData_df))])
new_sData_df <- unnest(new_sData_df, NegGeoMean)

#new_sData_df = new_sData_df[-grep('NegGeoMean',colnames(new_sData_df)),]
library(data.table)
fwrite(new_sData_df, paste0(datadir, 'QC_All.txt'), quote = F, row.names = F, sep = '\t')

dim(Data[, QCResults$QCStatus == "PASS"])
#============================================

#============================================ removeQCSampleProbe, eval = TRUE-----------------------------------------
Data_backup = Data
QCResults_backup = QCResults

Data <- Data[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(Data)

#============================================ Set Probe QC Flags
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
Data <- setBioProbeQCFlags(Data, 
                           qcCutoffs = list(minProbeRatio = 0.1,
                                            percentFailGrubbs = 20), 
                           removeLocalOutliers = TRUE)

ProbeQCResults <- sData(Data)[["QCFlags"]]

# I do not have GlobalGrubbsOutlier
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0))

kable(qc_df, caption = "Probes flagged or passed as outliers")

#============================================  Exclude Outlier Probes
Data <- 
  subset(Data, 
         fData(Data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(Data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(Data)


#============================================ Create Gene-level Count Data
# The count for any gene with multiple probes per segment is calculated as the geometric mean of those probes.
# https://bioconductor.riken.jp/packages/3.15/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
length(unique(featureData(Data)[["TargetName"]]))
 
# collapse to targets
target_Data <- aggregateCounts(Data)
dim(target_Data)
# exprs(target_Data)[1:5, 1:2]

exprs(Data)[1:5, 1:2]

#============================================ number of detected genes per segment
d <- data.frame(exprs(target_Data));
d.tst <- ifelse(d > 1, 1, 0) %>% as.data.frame

library(plyr)
colwise(sum)(d.tst)

#============================================ Limit of Quantification
# The LOQ is calculated based on the distribution of negative control probes and is intended to approximate the quantifiable limit of gene expression per segment.
# Please note that this process is more stable in larger segments.

# Define LOQ SD threshold and minimum value
cutoff <- 1
minLOQ <- 2 # We also recommend that a minimum LOQ of 2 be used if the LOQ calculated in a segment is below this threshold.

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_Data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_Data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_Data)[, vars[1]] * 
             pData(target_Data)[, vars[2]] ^ cutoff)
  }
}
pData(target_Data)$LOQ <- LOQ

#============================================ Filtering
#  After determining the limit of quantification (LOQ) per segment, we recommend filtering out either segments and/or genes with abnormally low signal. Filtering is an important step to focus on the true biological data of interest.

LOQ_Mat <- c()
module = modules
Mat_i <- t(esApply(target_Data, MARGIN = 1,
                   FUN = function(x) {
                     x > LOQ[, module]
                   }))
LOQ_Mat <- Mat_i

d <- data.frame(Mat_i);
d.tst <- ifelse(d == 'TRUE', 1, 0) %>% as.data.frame
library(plyr)
colwise(sum)(d.tst)

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_Data)$Module == module
  Mat_i <- t(esApply(target_Data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_Data)$TargetName, ]

#============================================ Segment Gene Detection
# We first filter out segments with exceptionally low signal. These segments will have a small fraction of panel genes detected above the LOQ relative to the other segments in the study.

# Save detection rate information to pheno data
pData(target_Data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_Data)$GeneDetectionRate <-
  pData(target_Data)$GenesDetected / nrow(target_Data)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_Data)$DetectionThreshold <- 
  cut(pData(target_Data)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_Data),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

#============================================ save QC incluidng gene detection rate
dimnames(pData(target_Data))[[2]]
target_Data_df <- pData(target_Data)[, c("GenesDetected", "GeneDetectionRate", "DetectionThreshold")]

target_Data_df$SampleID = gsub('.dcc', '', row.names(target_Data_df))
QC_Complete = merge(target_Data_df, new_sData_df, by = 'SampleID')
dim(QC_Complete)

getwd()

QC_Complete$GeneDetectionFail = rep('NA', dim(QC_Complete)[1])
QC_Complete$GeneDetectionFail[which(QC_Complete$DetectionThreshold %in% c('<1%', '1-5%'))] = 'TRUE'
QC_Complete$GeneDetectionFail[-which(QC_Complete$DetectionThreshold %in% c('<1%', '1-5%'))] = 'FALSE'
table(QC_Complete$GeneDetectionFail)

fwrite(QC_Complete, paste0(datadir, 'QC_All_inc_GeneDetection.txt'), quote = F, row.names = F, sep = '\t')

#============================================ visualize the intersection
# library(data.table)
# QC_Complete = fread(paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/QC_All_inc_GeneDetection.txt'), stringsAsFactors = F)
# table(QC_Complete$GeneDetectionPass)
# QC_Complete = as.data.frame(QC_Complete)

# QC_Completesubset = QC_Complete[, which(colnames(QC_Complete) %in% c('GeneDetectionFail','LowReads','LowTrimmed','LowStitched','LowAligned','LowSaturation','HighNTC'))]
# 
# colnames(QC_Completesubset)
# 
# for(i in 1:dim(QC_Completesubset)[1])
# {
#   QC_Completesubset$SatLA[i] = all(QC_Completesubset[i,c('LowStitched','LowAligned')])
#   QC_Completesubset$SatLR[i] = all(QC_Completesubset[i,c('LowSaturation','LowReads')])
#   QC_Completesubset$SatLT[i] = all(QC_Completesubset[i,c('LowSaturation','LowTrimmed')])
#   QC_Completesubset$SatLS[i] = all(QC_Completesubset[i,c('LowSaturation','LowStitched')])
# }
# 
# gsub(" ", "','", capture.output(cat(colnames(QC_Completesubset))))
# 
# cTable <- data.frame(matrix(ncol = 12, nrow = 12))
# row.names(cTable) = c('LowReads','LowTrimmed','LowStitched','LowAligned','LowSaturation','HighNTC','GeneDetectionFail','LL','SatNTC','GeneDNTC','GeneDSat','GeneDSatNTC')
# colnames(cTable) = c('LowReads','LowTrimmed','LowStitched','LowAligned','LowSaturation','HighNTC','GeneDetectionFail','LL','SatNTC','GeneDNTC','GeneDSat','GeneDSatNTC')
# 
# for (r in 1:nrow(cTable)) {
#   for (c in 1:ncol(cTable)) {
#     
#     print(table(QC_Completesubset[,row.names(cTable)[r]], QC_Completesubset[,colnames(cTable)[c]]))
#     
#     cTable[r, c] = table(QC_Completesubset[,row.names(cTable)[r]], QC_Completesubset[,colnames(cTable)[c]])[2,2]
#     
#   }
# }
# 
# View(cTable[c('LowReads','LowTrimmed','LowStitched','LowAligned','LowSaturation','HighNTC','GeneDetectionFail'),c('LowReads','LowTrimmed','LowStitched','LowAligned','LowSaturation','HighNTC','GeneDetectionFail')])
# View(cTable[c('LL','SatNTC','GeneDNTC','GeneDSat','GeneDSatNTC'),c('LL','SatNTC','GeneDNTC','GeneDSat','GeneDSatNTC')])

#============================================ gene detection as table

# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_Data)$DetectionThreshold,
            pData(target_Data)$segment))

#============================================
toMatch <- gsub('.dcc','',row.names(pData(target_Data))[which(pData(target_Data)$DetectionThreshold %in% c('<1%', '1-5%'))])

table(pData(target_Data)$DetectionThreshold %in% c('<1%', '1-5%'))
kable(table(Batch = SamplesMetadata$batch[which(SamplesMetadata$Sample_ID %in% toMatch )]))

#============================================ filter Segments based on gene detection rate
target_Data <-
  target_Data[, -which(pData(target_Data)$DetectionThreshold %in% c('<1%', '1-5%'))]
dim(target_Data)

#============================================ summary of samples
#--- some test for gene counts
summary(pData(target_Data)$GenesDetected)

d <- data.frame(Mat_i);
d.tst <- ifelse(d == 'TRUE', 1, 0) %>% as.data.frame

library(plyr)
colwise(sum)(d.tst)

d <- data.frame(exprs(target_Data));
d.tst <- ifelse(d == 1, 1, 0) %>% as.data.frame
library(plyr)
colwise(sum)(d.tst)

pData(target_Data)$GenesDetected <- as.numeric(colwise(sum)(d.tst))
pData(target_Data)$GeneDetectionRate <-
  pData(target_Data)$GenesDetected / nrow(target_Data)

dim(target_Data)
table(pData(target_Data)$segment)

df_pData = pData(target_Data)
df_pData$Donor = gsub(' .*', '', df_pData$`slide name`)

table(df_pData$Donor)
table(df_pData$segment)

dim(pData(target_Data))
sData_df <- sData(target_Data)[, dimnames(sData(target_Data))[[2]]]
dim(sData_df)

summary(sData_df$NucleiCount[sData_df$segment == 'control'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-with-LBP'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-without-LBP'])

#------------ Convert NanoStringGeoMxSet to SpatialExperiment
# The three errors that can occur when trying to coerce to Seurat are:
#   
# object must be on the target level
# object should be normalized, if you want raw data you can set forceRaw to TRUE
# normalized count matrix name must be valid


#============================================
target_Data <- normalize(target_Data, norm_method="quant", desiredQuantile = .75, toElt = "q_norm")
target_Data <- normalize(target_Data, norm_method = "neg", fromElt = "exprs", toElt = "neg_norm")
exprs(target_Data)[1:5, 1:2]

data.frame(assayData(target_Data)[["exprs"]][seq_len(3), seq_len(3)])
featureType(target_Data)
data.frame(assayData(target_Data)[["q_norm"]][seq_len(3), seq_len(3)])

library(stringr)
string = str_split_fixed(gsub('.dcc','',colnames(exprs(target_Data))), "-", 4)

# boxplot(exprs(target_Data)[,1:50],
#         col = "#9EDAE5", main = "Raw Counts",
#         log = "y", names = gsub('.dcc|DSP-1001660016912-A-','',string[1:50,4]), xlab = "Segment",
#         ylab = "Counts, Raw")
# 
# boxplot(assayDataElement(target_Data, elt = "q_norm"),
#         col = "#2CA02C", main = "Q3 Norm Counts",
#         log = "y", names = gsub('.dcc|DSP-1001660016912-A-','',colnames(exprs(target_Data))), xlab = "Segment",
#         ylab = "Counts, Q3 Normalized")
# 
# boxplot(assayDataElement(target_Data, elt = "neg_norm"),
#         col = "#FF7F0E", main = "Neg Norm Counts",
#         log = "y", names = gsub('.dcc|DSP-1001660016912-A-','',colnames(exprs(target_Data))), xlab = "Segment",
#         ylab = "Counts, Neg. Normalized")

#============================================ using Seurat and scater for further analysis
options(Seurat.object.assay.version = "v3")
target_Data_Seurat <- as.SpatialExperiment(target_Data, normData = "q_norm")
counts <- SummarizedExperiment::assay(target_Data_Seurat)
meta <- SummarizedExperiment::colData(target_Data_Seurat) %>% as.data.frame

library(Seurat)
target_Data_Seurat <- CreateSeuratObject(counts = counts, project = "Seurat",  meta.data = meta, assay = "GeoMx")

head(target_Data_Seurat, 3) 
target_Data_Seurat$segment

library(scater)
sce <- as.SingleCellExperiment(target_Data_Seurat)
# sce_mat <- assay(pbmc.sce) 
# sce_mat = as.matrix(sce_mat)
# sce_mat[1:5,1:5]

sce <- logNormCounts(sce)

library(dittoSeq)
dittoPlot(sce, "PEG10", group.by = "segment")

#============================================
library(standR)
library(SpatialExperiment)

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

#============================================ generate SpatialExperiment
spe <- as.SpatialExperiment(target_Data, normData = "neg_norm", forceRaw = T)
spe

assayNames(spe)

duplicated(names(assay(spe, assayNames(spe))[1,]))

# -- add raw counts to the object
Temp = target_Data@assayData$exprs
row.names(Temp) = target_Data@featureData@data$TargetName

assay(spe, 'counts') = Temp
assayNames(spe)[which(assayNames(spe) == 'GeoMx')] = 'logcounts'
assay(spe, 'logcounts') = Temp

assayNames(spe)
assay(spe)[1:5,1:5]

#============================================ Sample metadata is stored in the colData of the object.
colData(spe)[1:5,1:5]

which(names(colData(spe)) == 'segment')

# colData(spe)[which(as.data.frame(colData(spe)[4])[,1] == 'BP+TH+'), 4] = "LB509+TH+"
# colData(spe)[which(as.data.frame(colData(spe)[4])[,1] == 'BP+'), 4] = "TH+"

MD = data.frame((colData(spe)))
fwrite(MD, paste0(datadir, 'MetaData.txt'), quote = F, row.names = F, sep = '\t')

metadata(spe)$NegProbes = negativeGeoMeans

#============================================ Import from DGEList object
# Alternatively, standR provides a function to generate a spatial experiment object from a DGEList object, which would be useful for users who used edgeR package and have existing analyses and implementations using DGEList objects to port across to the standR workflow.
dge <- edgeR::SE2DGEList(spe)
spe2 <- readGeoMxFromDGE(dge)
spe2

#============================================ Check QCFlags
# In the meta data file generated by Nanostring, there is a column called “QCFlags”, which indicates bad quality tissue samples in their preliminary QC step. If the data is fine from their QC, you will see NA/empty cells in this column.
colData(spe2)$QCFlags 
protocolData(target_Data)[["QCFlags"]]

#============================================ Sample level QC
# To visualise sample metadata
library(ggplot2)
library(ggalluvial)
names(colData(spe2))
spe2$`slide name` = 'SN'
#plotSampleInfo(spe2, column2plot = c("slide name", "segment", "roi"))

#============================================ Gene level QC
length(SummarizedExperiment::assayNames(spe2))

# sample_fraction	Double. Genes with low count in more than this threshold of the samples will be removed. Default is 0.9
# min_count Integer. Minimum read count to calculate count threshold. Default is 5.
spe2_subset <- addPerROIQC(spe2, rm_genes = TRUE, min_count = 4)
metadata(spe2_subset) |> names()

# table(pData(target_Data)$segment)[1]/sum(table(pData(target_Data)$segment))

colData(spe2_subset)$Names = gsub('.dcc|DSP-1001660016912-A-','', row.names(colData(spe2_subset)))
#plotGeneQC(spe2_subset, ordannots = "Names", col = Names, point_size = 2)

library(stringr)
string = str_split_fixed(gsub('.dcc','',colData(spe2_subset)$Names), "-", 4)
colData(spe2_subset)$Names = string[,4]

colData(spe2_subset)$state = factor(colData(spe2_subset)$state, levels =  c("PD-with-LBP", "PD-without-LBP", "control"))

#============================================ Relative log expression distribution
plotRLExpr(spe2_subset)

plotRLExpr(spe2_subset, ordannots = "state", assay = 2, color = state)
plotRLExpr(spe2_subset, ordannots = "state", assay = 2, color = area)

drawPCA(spe2_subset, assay = 2, color = segment)

set.seed(100)
spe2_subset <- scater::runPCA(spe2_subset)
pca_results <- reducedDim(spe2_subset, "PCA")
drawPCA(spe2_subset, precomputed = pca_results, col = segment)
plotScreePCA(spe2_subset, precomputed = pca_results)

standR::plotMDS(spe2_subset, assay = 2, color = segment)

#============================================ Normalization
spe_tmm <- geomxNorm(spe2_subset, method = "TMM")
plotRLExpr(spe_tmm, assay = 2, color = segment) + ggtitle("TMM")
drawPCA(spe_tmm, assay = 2, color = segment)

#============================================ Batch correction
spe_tmm <- findNCGs(spe_tmm, batch_name = "segment", top_n = 300)
metadata(spe_tmm) |> names()

for(i in seq(10)){
  spe_ruv <- geomxBatchCorrection(spe_tmm, factors = "segment", 
                                  NCGs = metadata(spe_tmm)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = segment, title = paste0("k = ", i)))
  
}

spe_ruv <- geomxBatchCorrection(spe_tmm, factors = "segment", NCGs = metadata(spe_tmm)$NCGs, k = 5)
print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = segment, title = paste0("k = ", i)))

plotRLExpr(spe_ruv, assay = 2, color = segment) + ggtitle("TMM + RUV4")

#============================================ Differential expression analysis with limma-voom pipeline
colData(spe_ruv)[,seq(ncol(colData(spe_ruv))-1, ncol(colData(spe_ruv)))] |> head()

library(edgeR)
library(limma)

dge <- SE2DGEList(spe_ruv)

# 
design <- model.matrix(~0 + segment + ruv_W1 + ruv_W2 , data = colData(spe_ruv))
colnames(design)

colnames(design) <- gsub("^segment","",colnames(design))
colnames(design) <- gsub("\\+","",colnames(design))
#colnames(design) = c("TH", "LB509TH",  "ruv_W1",   "ruv_W2")

colnames(design) = c("control", "PDwLBP", "PDwoLBP",  "ruv_W1",   "ruv_W2")

# LB509TH - TH (TH is reference)
contr.matrix <- makeContrasts(BvT = PDwLBP - control, levels = colnames(design))

keep <- filterByExpr(dge, design)
table(keep)
rownames(dge)[!keep]

dge_all <- dge[keep, ]

dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)

plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))
bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)

highbcv <- bcv_df$BCV > 0.8
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)

v <- voom(dge_all, design, plot = TRUE) 

fit <- lmFit(v)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit <- decideTests(efit, p.value = 0.05)
summary_efit <- summary(results_efit)
summary_efit

library(ggrepel)
library(tidyverse)
de_results_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)
#options(ggrepel.max.overlaps = 10)

#---------------------------- calculating sample size ----------------------------
#- effect size (Log fold change (LFC) estimates)
de_genes_toptable_BvT_001 <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.01)
dim(de_genes_toptable_BvT_001)
summary(de_genes_toptable_BvT_001$logFC) # The fact that limma gives log2 results is mentioned many times in the documentation.
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.7048 -0.8966 -0.6015 -0.1493  0.9815  2.0953 
# A log2 fold change of 1 means the gene is expressed twice as high in the test condition.
# A log2 fold change of 2 means the gene is expressed four times higher in the test condition.
# A log2 fold change of 3 means the gene is expressed eight times higher in the test condition, and so on.

# Cohen suggests that d values of 0.2, 0.5, and 0.8 represent small, medium, and large effect sizes respectively.
# I select 0.5

results <- pwr::pwr.t.test(n = NULL,
                           sig.level = 0.05, 
                           type = "two.sample", 
                           alternative = "two.sided", 
                           power = 0.80, 
                           d = 0.5)

plot(results) +
  ggplot2::theme_minimal(base_size = 14) +
  labs(title = 'Optimizing Sample Size for 2-Sided t test',
       subtitle = "")

#------------------------------------------------------------------------------------
table(de_results_BvT$adj.P.Val < 0.05 & de_results_BvT$logFC > 0)
table(de_results_BvT$adj.P.Val < 0.05 & de_results_BvT$logFC < 0)
dim(de_results_BvT)

de_results_BvT %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                     ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
  ggplot(aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = de_genes_toptable_BvT %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("LB509+TH+ vs LB509-TH+ as reference (limma-voom)") +
  scale_color_manual(values = c("dodgerblue","lightblue","orange2")) +
  theme(text = element_text(size=15))

dev.off()

library(DT)
updn_cols <- c(RColorBrewer::brewer.pal(6, 'Greens')[2], RColorBrewer::brewer.pal(6, 'Purples')[2])

library(data.table)
fwrite(de_genes_toptable_BvT, paste0(datadir, 'GeoMx_DEA.txt'), quote = F, row.names = F, sep = '\t')

de_genes_toptable_BvT %>% 
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  DT::datatable(caption = 'B cell zone vs. T cell zone in Lymph node (limma-voom)') %>%
  DT::formatStyle('logFC',
                  valueColumns = 'logFC',
                  backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
  DT::formatSignif(1:4, digits = 4)

library(msigdb)
library(GSEABase)

msigdb_hs <- getMsigdb(version = '7.2')
msigdb_hs <- appendKEGG(msigdb_hs)
# load PPI from the msigdb package
# ppi = getIMEX('hs', inferred = TRUE)

sc <- listSubCollections(msigdb_hs)

gsc <- c(subsetCollection(msigdb_hs, c('h')),
         subsetCollection(msigdb_hs, 'c2', sc[grepl("^CP:",sc)]),
         subsetCollection(msigdb_hs, 'c5', sc[grepl("^GO:",sc)])) %>%
  GeneSetCollection()

#------- Enrichment analysis
# Preprocessing is conducted on these genesets, filtering out genesets with less than 5 genes and creating indices vector list for formatting on the results before applying fry.

fry_indices <- ids2indices(lapply(gsc, geneIds), rownames(v), remove.empty = FALSE)
names(fry_indices) <- sapply(gsc, setName)

gsc_category <- sapply(gsc, function(x) bcCategory(collectionType(x)))
gsc_category <- gsc_category[sapply(fry_indices, length) > 5]

gsc_subcategory <- sapply(gsc, function(x) bcSubCategory(collectionType(x)))
gsc_subcategory <- gsc_subcategory[sapply(fry_indices, length) > 5]

fry_indices <- fry_indices[sapply(fry_indices, length) > 5]

names(gsc_category) = names(gsc_subcategory) = names(fry_indices)

# Now we run fry with all the gene sets we filtered.

fry_indices_cat <- split(fry_indices, gsc_category[names(fry_indices)])
fry_res_out <- lapply(fry_indices_cat, function (x) {
  limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)
})

post_fry_format <- function(fry_output, gsc_category, gsc_subcategory){
  names(fry_output) <- NULL
  fry_output <- do.call(rbind, fry_output)
  fry_output$GenesetName <- rownames(fry_output)
  fry_output$GenesetCat <- gsc_category[rownames(fry_output)]
  fry_output$GenesetSubCat <- gsc_subcategory[rownames(fry_output)]
  return(fry_output)
}

fry_res_sig <- post_fry_format(fry_res_out, gsc_category, gsc_subcategory) %>%
  as.data.frame() %>%
  filter(FDR < 0.05) 

# The output is a data.frame object. We can either output the whole table, or inspect the top N gene sets in a bar plot.
# We can see many immune-related gene sets are significantly enriched, B cell-related gene-sets are enriched in up-regulated genes while T-cell related gene-sets are enriched in down-regulated genes.

fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Up") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "orange2") +
  theme_bw() +
  coord_flip() +
  ggtitle("Up-regulated")
 
fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Down") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  theme_bw() +
  coord_flip() +
  ggtitle("Down-regulated")

# Visualization
# An alternative way to summarise the GSEA output is to visualise common gene sets as a group.
# We can use the igraph and vissE package to perform clustering on the enriched gene sets and visualise the gene sets using word cloud-based algorithm and network-based visualisation. For more information about vissE, check out here.
library(vissE)
library(igraph)

fry_out = fry_res_sig
de_table = de_genes_toptable_BvT
topN = 6
title = ""
specific_clusters = NA

dovissE <- function(fry_out, de_table, topN = 6, title = "", specific_clusters = NA){
  
  n_row = min(1000, nrow(fry_out))
  gs_sig_name <- fry_out %>% 
    filter(FDR < 0.05) %>%
    arrange(FDR) %>% 
    .[1:n_row,] %>% 
    rownames()
  gsc_sig <- gsc[as.numeric(gs_sig_name),]
  
  gs_ovlap <- computeMsigOverlap(gsc_sig, thresh = 0.15)
  gs_ovnet <- computeMsigNetwork(gs_ovlap, gsc)
  
  max(fry_out[as.numeric(gs_sig_name),]$FDR)
  
  gs_stats <- -log10(fry_out[as.numeric(gs_sig_name),]$FDR)
  names(gs_stats) <- as.numeric(gs_sig_name)
  
  #identify clusters
  grps = cluster_walktrap(gs_ovnet)
  #extract clustering results
  grps = groups(grps)
  #sort by cluster size
  grps = grps[order(sapply(grps, length), decreasing = TRUE)]
  
  # write output
  output_clusters <- list()
  for(i in seq(length(grps))){
    output_clusters[[i]] <- data.frame(geneset = grps[[i]], cluster = paste0("cluster",names(grps)[i]))
  }
  output_clusters <<- output_clusters %>% bind_rows()
  
  if(is.na(specific_clusters)){
    grps <- grps[1:topN]
  } else {
    grps <- grps[specific_clusters %>% as.character()]
  }
  
  #plot the top 12 clusters
  set.seed(36) #set seed for reproducible layout
  p1 <<- plotMsigNetwork(gs_ovnet, markGroups = grps, 
                         genesetStat = gs_stats, rmUnmarkedGroups = TRUE) +
    scico::scale_fill_scico(name = "-log10(FDR)")
  
  p2 <<- plotMsigWordcloud(gsc, grps, type = 'Name')
  
  genes <- unique(unlist(geneIds(gsc_sig)))
  genes_logfc <- de_table %>% rownames_to_column() %>% filter(rowname %in% genes) %>% .$logFC
  names(genes_logfc) <- de_table %>% rownames_to_column() %>% filter(rowname %in% genes) %>% .$rowname
  p3 <<- plotGeneStats(genes_logfc, gsc, grps) +
    geom_hline(yintercept = 0, colour = 2, lty = 2) +
    ylab("logFC")
  
  # p4 <- plotMsigPPI(ppi, gsc, grps[1:topN], geneStat = genes_logfc) +
  #  guides(col=guide_legend(title="logFC"))
  
  print(p2 + p1 + p3 + patchwork::plot_layout(ncol = 3) +
          patchwork::plot_annotation(title = title))  
  
}

dovissE(fry_res_sig, de_genes_toptable_BvT, topN = 3, title = "B cell zone vs. T cell zone in Lymph node." )

#-------- consider the SAVER approach for data normalization
# setwd('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output')
# saver = fread('SAVER.txt', stringsAsFactors = F, header = T)
# saver = as.data.frame(saver)
# row.names(saver) = saver$GeneID
# saver = saver[,-1]
# dim(saver)
# 
# boxplot(saver[which(rownames(saver) %in% rownames(spe_ruv)),])
# 
# dim(colData(spe_ruv))
# names(metadata(spe))
# 
# assay(spe2_subset, 'logcounts')[1:10,1]
# assay(spe_ruv, 'logcounts')[1:10,1:10]
# 
# # logcounts is not log, and saves the results of the trimmed mean of M values (TMM) + RUV4
# boxplot(assay(spe2_subset, 'logcounts'))
# boxplot(assay(spe_ruv, 'logcounts'))
# 
# # this is log of TMM + RUV4
# plotRLExpr(spe_ruv, assay = 2, color = segment) + ggtitle("TMM + RUV4")
# 
# assay(spe_ruv, 'logcounts')[1:5,1:5]
# count_mat = data.frame(GeneSymbol = rownames(assay(spe_ruv, 'logcounts')), assay(spe_ruv, 'logcounts'))
# 
# library(data.table)
# fwrite(count_mat, paste0(datadir, 'GeoMx_readCount_TMMRUV4.txt'), quote = F, row.names = F, sep = '\t')

#---------------------- visualization (pairplot)

library(data.table)
datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/")
count_mat2 = fread(paste0(datadir, 'GeoMx_readCount.txt'), stringsAsFactors = F, header = T)
count_mat2 = as.data.frame(count_mat2)
row.names(count_mat2) = make.names(count_mat2$GeneSymbol,unique=T) 
count_mat2 = count_mat2[,-1]
colnames(count_mat2) = gsub('\\.','-',gsub('.dcc','',colnames(count_mat2)))

count_mat3 = fread(paste0(datadir, 'GeoMx_readCount_TMMRUV4.txt'), stringsAsFactors = F, header = T)
count_mat3 = as.data.frame(count_mat3)
row.names(count_mat3) = make.names(count_mat3$GeneSymbol,unique=T) 
count_mat3 = count_mat3[,-1]
colnames(count_mat3) = gsub('\\.','-',gsub('.dcc','',colnames(count_mat3)))

count_mat2 = count_mat2[which(row.names(count_mat2) %in% row.names(count_mat3)), which(colnames(count_mat2) %in% colnames(count_mat3))]
dim(count_mat2)

MD = fread(paste0(datadir, 'MetaData.txt'), stringsAsFactors = F, header = T)
MD = as.data.frame(MD)
colnames(MD)

Reads <- apply(count_mat2[,which(colnames(count_mat2) %in% MD$SampleID)], 2, function(x) sum(x))
#-- nGenesPerCell
nGenesPerCell <- apply(count_mat2, 2, function(x) sum(x>0))

identical(names(Reads), MD$SampleID)
identical(names(nGenesPerCell), MD$SampleID)

plot(as.numeric(nGenesPerCell), MD$GenesDetected)

dataset = data.frame(Nuclei_Count = MD$NucleiCount, Genes_Detected = nGenesPerCell,
                     Area_Size = MD$area,
                     Reads_Count = as.numeric(Reads) )

colnames(dataset) = c("Nuclei Count", "Genes Detected", "Area Size", "Reads Count")
variables <- c("Nuclei Count", "Genes Detected", "Area Size", "Reads Count")
str(dataset)

# ---------------------------------- 

.plotMarginalCor <- function(variable, cexYlab = 1.3, lwd = 2, rugs = FALSE) {
  
  # histogram with density estimator
  
  variable <- variable[!is.na(variable)]
  
  density <- density(variable)
  h <- hist(variable, plot = FALSE)
  jitVar <- jitter(variable)
  yhigh <- max(max(h$density), max(density$y))
  ylow <- 0
  xticks <- pretty(c(variable, h$breaks), min.n = 3)
  plot(range(xticks), c(ylow, yhigh), type = "n", axes = FALSE, ylab = "", 
       xlab = "")
  h <- hist(variable, freq = FALSE, main = "", ylim = c(ylow, yhigh), xlab = "", 
            ylab = " ", axes = FALSE, col = "grey", add = TRUE, nbreaks = round(length(variable)/5))
  ax1 <- axis(1, line = 0.3, at = xticks, lab = xticks)
  par(las = 0)
  ax2 <- axis(2, at = c(0, max(max(h$density), max(density$y))/2, max(max(h$density), 
                                                                      max(density$y))), labels = c("", "Density", ""), lwd.ticks = 0, 
              pos = range(ax1) - 0.08 * diff(range(ax1)), cex.axis = 2, mgp = c(3, 0.7, 0))
  
  if (rugs) 
    rug(jitVar)
  
  lines(density$x[density$x >= min(ax1) & density$x <= max(ax1)], density$y[density$x >= 
                                                                              min(ax1) & density$x <= max(ax1)], lwd = lwd)
}

.poly.pred <- function(fit, line = FALSE, xMin, xMax, lwd) {
  
  # predictions of fitted model
  
  # create function formula
  f <- vector("character", 0)
  
  for (i in seq_along(coef(fit))) {
    
    if (i == 1) {
      
      temp <- paste(coef(fit)[[i]])
      f <- paste(f, temp, sep = "")
    }
    
    if (i > 1) {
      
      temp <- paste("(", coef(fit)[[i]], ")*", "x^", i - 1, sep = "")
      f <- paste(f, temp, sep = "+")
    }
  }
  
  x <- seq(xMin, xMax, length.out = 100)
  predY <- eval(parse(text = f))
  
  if (line == FALSE) {
    return(predY)
  }
  
  if (line) {
    lines(x, predY, lwd = lwd)
  }
}


.plotScatter <- function(xVar, yVar, cexPoints = 1.3, cexXAxis = 1.3, 
                         cexYAxis = 1.3, lwd = 2) {
  
  # displays scatterplot
  
  d <- data.frame(xx = xVar, yy = yVar)
  d <- na.omit(d)
  xVar <- d$xx
  yVar <- d$yy
  
  fit <- lm(yy ~ poly(xx, 1, raw = TRUE), d)
  
  xlow <- min((min(xVar) - 0.1 * min(xVar)), min(pretty(xVar)))
  xhigh <- max((max(xVar) + 0.1 * max(xVar)), max(pretty(xVar)))
  xticks <- pretty(c(xlow, xhigh))
  ylow <- min((min(yVar) - 0.1 * min(yVar)), min(pretty(yVar)), min(.poly.pred(fit, 
                                                                               line = FALSE, xMin = xticks[1], xMax = xticks[length(xticks)], 
                                                                               lwd = lwd)))
  yhigh <- max((max(yVar) + 0.1 * max(yVar)), max(pretty(yVar)), max(.poly.pred(fit, 
                                                                                line = FALSE, xMin = xticks[1], xMax = xticks[length(xticks)], 
                                                                                lwd = lwd)))
  
  
  yticks <- pretty(c(ylow, yhigh))
  
  yLabs <- vector("character", length(yticks))
  
  for (i in seq_along(yticks)) {
    
    if (yticks[i] < 10^6) {
      
      yLabs[i] <- format(yticks[i], digits = 3, scientific = FALSE)
      
    } else {
      
      yLabs[i] <- format(yticks[i], digits = 3, scientific = TRUE)
    }
  }
  
  plot(xVar, yVar, col = "black", pch = 21, bg = "grey", ylab = "", 
       xlab = "", axes = FALSE, ylim = range(yticks), xlim = range(xticks), 
       cex = cexPoints)
  .poly.pred(fit, line = TRUE, xMin = xticks[1], xMax = xticks[length(xticks)], 
             lwd = lwd)
  
  par(las = 1)
  
  axis(1, line = 0.4, labels = xticks, at = xticks, cex.axis = cexXAxis)
  axis(2, line = 0.2, labels = yLabs, at = yticks, cex.axis = cexYAxis)
  
  invisible(max(nchar(yLabs)))
}

.plotCorValue <- function(xVar, yVar, cexText = 1.7, cexCI = 1.7, hypothesis = "correlated", 
                          pearson = TRUE, kendallsTauB = FALSE, spearman = FALSE, confidenceInterval = 0.95) {
  
  # displays correlation value
  
  CIPossible <- TRUE
  
  tests <- c()
  
  if (pearson) 
    tests <- c(tests, "pearson")
  
  if (spearman) 
    tests <- c(tests, "spearman")
  
  if (kendallsTauB) 
    tests <- c(tests, "kendall")
  
  plot(1, 1, type = "n", axes = FALSE, ylab = "", xlab = "")
  
  lab <- vector("list")
  
  for (i in seq_along(tests)) {
    
    if (round(cor.test(xVar, yVar, method = tests[i])$estimate, 8) == 
        1) {
      
      CIPossible <- FALSE
      
      if (tests[i] == "pearson") {
        lab[[i]] <- bquote(italic(r) == "1.000")
      }
      
      if (tests[i] == "spearman") {
        lab[[i]] <- bquote(italic(rho) == "1.000")
      }
      
      if (tests[i] == "kendall") {
        lab[[i]] <- bquote(italic(tau) == "1.000")
      }
      
    } else if (round(cor.test(xVar, yVar, method = tests[i])$estimate, 
                     8) == -1) {
      
      CIPossible <- FALSE
      
      if (tests[i] == "pearson") {
        lab[[i]] <- bquote(italic(r) == "-1.000")
      }
      
      if (tests[i] == "spearman") {
        lab[[i]] <- bquote(italic(rho) == "-1.000")
      }
      
      if (tests[i] == "kendall") {
        lab[[i]] <- bquote(italic(tau) == "-1.000")
      }
      
    } else {
      
      if (tests[i] == "pearson") {
        lab[[i]] <- bquote(italic(r) == .(formatC(round(cor.test(xVar, 
                                                                 yVar, method = tests[i])$estimate, 3), format = "f", 
                                                  digits = 3)))
      }
      
      if (tests[i] == "spearman") {
        lab[[i]] <- bquote(rho == .(formatC(round(cor.test(xVar, 
                                                           yVar, method = tests[i])$estimate, 3), format = "f", 
                                            digits = 3)))
      }
      
      if (tests[i] == "kendall") {
        lab[[i]] <- bquote(tau == .(formatC(round(cor.test(xVar, 
                                                           yVar, method = tests[i])$estimate, 3), format = "f", 
                                            digits = 3)))
      }
    }
  }
  
  if (length(tests) == 1) {
    ypos <- 1
  }
  
  if (length(tests) == 2) {
    ypos <- c(1.1, 0.9)
  }
  
  if (length(tests) == 3) {
    ypos <- c(1.2, 1, 0.8)
  }
  
  
  for (i in seq_along(tests)) {
    
    text(1, ypos[i], labels = lab[[i]], cex = cexText)
  }
  
  
  if (hypothesis == "correlated" & length(tests) == 1 & any(tests == 
                                                            "pearson")) {
    
    alternative <- "two.sided"
    ctest <- cor.test(xVar, yVar, method = tests, conf.level = confidenceInterval)
  }
  
  if (hypothesis != "correlated" & length(tests) == 1 & any(tests == 
                                                            "pearson")) {
    
    if (hypothesis == "correlatedPositively") {
      
      ctest <- cor.test(xVar, yVar, method = tests, alternative = "greater", 
                        conf.level = confidenceInterval)
      
    } else if (hypothesis == "correlatedNegatively") {
      
      ctest <- cor.test(xVar, yVar, method = tests, alternative = "less", 
                        conf.level = confidenceInterval)
    }
    
  }
  
  if (any(tests == "pearson") & length(tests) == 1 && CIPossible) {
    
    CIlow <- formatC(round(ctest$conf.int[1], 3), format = "f", digits = 3)
    CIhigh <- formatC(round(ctest$conf.int[2], 3), format = "f", digits = 3)
    
    text(1, 0.7, labels = paste(100 * confidenceInterval, "% CI: [", 
                                CIlow, ", ", CIhigh, "]", sep = ""), cex = cexCI)
  }
  
}

### matrix plot ###
l <- length(variables)

par(mfrow = c(l, l), cex.axis = 1.3, mar = c(3, 4, 2, 1.5) + 0.1, oma = c(0, 2.2, 2, 0))

for (row in seq_len(l)) {
  
  for (col in seq_len(l)) {
    
    if (row == col) {
      .plotMarginalCor(dataset[[variables[row]]])  # plot marginal (histogram with density estimator)
    }
    if (col > row) {
      .plotScatter(dataset[[variables[col]]], dataset[[variables[row]]])  # plot scatterplot
    }
    if (col < row) {
      if (l < 7) {
        .plotCorValue(dataset[[variables[col]]], dataset[[variables[row]]], 
                      cexCI = 1.2)  # plot r= ...
      }
      if (l >= 7) {
        .plotCorValue(dataset[[variables[col]]], dataset[[variables[row]]], 
                      cexCI = 1.2)
      }
    }
  }
}

textpos <- seq(1/(l * 2), (l * 2 - 1)/(l * 2), 2/(l * 2))
for (t in seq_along(textpos)) {
  mtext(text = variables[t], side = 3, outer = TRUE, at = textpos[t], 
        cex = 0.9, line = -0.8)
  mtext(text = variables[t], side = 2, outer = TRUE, at = rev(textpos)[t], 
        cex = 0.9, line = -0.1)
}


#------------ comparison of GeoMx and snRNA-seq DEGs after integration

GeoMx_DEGs = fread('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/GeoMx_DEA.txt', stringsAsFactors = F)
GeoMx_DEGs = as.data.frame(GeoMx_DEGs)
dim(GeoMx_DEGs)

GeoMx_DEGs_sig = GeoMx_DEGs[which(GeoMx_DEGs$adj.P.Val < 0.05),]
dim(GeoMx_DEGs_sig)

GeoMx_DEGs_up = GeoMx_DEGs[which(GeoMx_DEGs$logFC > 0),]
dim(GeoMx_DEGs_up)

GeoMx_DEGs_down = GeoMx_DEGs[which(GeoMx_DEGs$logFC < 0),]
dim(GeoMx_DEGs_down)


snRNAseq_DEGs = fread('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/Annotated_DEGenes.txt', stringsAsFactors = F)
snRNAseq_DEGs = as.data.frame(snRNAseq_DEGs)
dim(snRNAseq_DEGs)

snRNAseq_DEGs_sig = snRNAseq_DEGs[which(snRNAseq_DEGs$padj < 0.05),]
dim(snRNAseq_DEGs_sig)

snRNAseq_DEGs_up = snRNAseq_DEGs_sig[which(snRNAseq_DEGs_sig$log2FoldChange > 0),]
dim(snRNAseq_DEGs_up)

snRNAseq_DEGs_down = snRNAseq_DEGs_sig[which(snRNAseq_DEGs_sig$log2FoldChange < 0),]
dim(snRNAseq_DEGs_down)


length(intersect(GeoMx_DEGs_sig$TargetName, snRNAseq_DEGs_sig$feature))
length(intersect(GeoMx_DEGs_up$TargetName, snRNAseq_DEGs_up$feature))
length(intersect(GeoMx_DEGs_down$TargetName, snRNAseq_DEGs_down$feature))




#============================================ test ======================================================================================== 
library(GeomxTools)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(tidyverse)

## Introduction
# This summarises the QC of the GeoMx pilot run of 10th October, utilising the 'GeoMxTools' pipeline.

#============================================ Collate files
datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/")
setwd(datadir)

DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/pkcs"), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)

# I added NucleiCount manually to the following file from Annie Initial Dataset.xlsx file in the 'DSPDA Files' folder
SampleAnnotationFile <-
  dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)

library(readxl)

#============================================ load data
Data <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = excel_sheets(path = SampleAnnotationFile),
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi","roi","state", "NucleiCount"),
                         experimentDataColNames = c("panel"))

#============================================ summary of empty samples
library(data.table)
excluded = fread('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/excluded.txt', stringsAsFactors = F, header = F)
colnames(excluded) = 'Sample_ID'
excluded$Sample_ID = gsub('.dcc','',excluded$Sample_ID)

SamplesMetadata = read_excel('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/annotation/Annotation_File.xlsx')

intersect(SamplesMetadata$Sample_ID, excluded$Sample_ID)

excluded_SamplesMetadata = merge(SamplesMetadata, excluded, by = 'Sample_ID')
#============================================

library(knitr)
pkcs <- annotation(Data)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

#============================================ Sample Overview
# Serial sections of a single brain section were analysed, with 7 ROI from BP+TH+ cell bodies and 7 ROI from BP+ cell bodies.

library(dplyr)
library(ggforce)
count_mat = data.frame(GeneSymbol = Data@featureData@data$TargetName, Data@assayData$exprs)

fwrite(count_mat, paste0(datadir, 'GeoMx_readCount.txt'), quote = F, row.names = F, sep = '\t')
dim(count_mat)
# summary(count_mat[,-1])

## Preprocessing and Segment QC
# Before we begin, we will shift any expression counts with a value of 0 to 1 to enable in downstream transformations.
# We then assess sequencing quality and adequate tissue sampling for every ROI/AOI segment.

# Every ROI/AOI segment will be tested for:
#   - Raw sequencing reads: segments with >1000 raw reads are removed.
# - % Aligned,% Trimmed, or % Stitched sequencing reads: segments below ~80% for one or more of these QC parameters are removed.
# - % Sequencing saturation ([1-deduplicated reads/aligned reads]%): segments below ~50% require additional sequencing to capture full sample diversity and are not typically analyzed until improved.
# - Negative Count: this is the geometric mean of the several unique negative probes in the GeoMx panel that do not target mRNA and establish the background count level per segment; segments with low negative counts (1-10) are not necessarily removed but may be studied closer for low endogenous gene signal and/or insufficient tissue sampling.
# - No Template Control (NTC) count: values >1,000 could indicate contamination for the segments associated with this NTC; however, in cases where the NTC count is between 1,000- 10,000, the segments may be used if the NTC data is uniformly low (e.g. 0-2 counts for all probes).
# - Nuclei: >100 nuclei per segment is generally recommended; however, this cutoff is highly study/tissue dependent and may need to be reduced; what is most important is consistency in the nuclei distribution for segments within the study.
# - Area: generally correlates with nuclei; a strict cutoff is not generally applied based on area.
# 
# In this case the nuclei QC metric was ignored.

Data <- shiftCountsOne(Data, useDALogic = TRUE)

summary(Data@assayData$exprs)

# QC_params <-
#   list(minSegmentReads = 1000, # Minimum number of reads (1000)
#        percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
#        percentStitched = 80,   # Minimum % of reads stitched (80%)
#        percentAligned = 75,    # Minimum % of reads aligned (80%)
#        percentSaturation = 50, # Minimum sequencing saturation (50%)
#        minNegativeCount = 1,   # Minimum negative control counts (10)
#        maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
#        minNuclei = 20,         # Minimum # of nuclei estimated (100)
#        minArea = 1000)         # Minimum segment area (5000)

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 10,   # Minimum negative control counts (10)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 100,         # Minimum # of nuclei estimated (100)
       minArea = 5000)         # Minimum segment area (5000)

Data <- setSegmentQCFlags(Data, qcCutoffs = QC_params)        

QC_Table = data.frame(Aligned = sData(Data)$Aligned, DeduplicatedReads = sData(Data)$DeduplicatedReads)
# fwrite(QC_Table, paste0(datadir, 'QC_Results.txt'), quote = F, row.names = F, sep = '\t')

#============================================ Collate QC Results
QCResults <- protocolData(Data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))

QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

col_by <- "segment"

#============================================ Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(Data), "NucleiCount", col_by, 20)

QC_histogram(sData(Data), "Trimmed (%)", col_by, 80)

# Stitching: These fragments are computationally merged into a single, longer read using specialized algorithms. This is essential because short reads alone might not uniquely map to the target gene sequences.
QC_histogram(sData(Data), "Stitched (%)", col_by, 80)

QC_histogram(sData(Data), "Aligned (%)", col_by, 75)

QC_histogram(sData(Data), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")

QC_histogram(sData(Data), "area", col_by, 1000, scale_trans = "log10")

## Generate negative probe count
# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(Data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(Data)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData & plot
negCols <- paste0("NegGeoMean_", modules)
pData(Data)[, negCols] <- sData(Data)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(Data), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(Data) <- pData(Data)[, !colnames(pData(Data)) %in% negCols]

#============================================ QC as tables

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(Data)$NTC), col.names = c("NTC Count", "# of Segments"))

## ----QCSummaryTable, results = "asis"-----------------------------------------
kable(QC_Summary, caption = "QC Summary Table for each Segment")
#============================================

#============================================ explore NTC  
highNTC = gsub('.dcc','',row.names(QCResults[QCResults$HighNTC,]))

dccFiles <- dir('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/dccs/', pattern=".dcc$", full.names=TRUE)
dccData <- sapply(dccFiles, readDccFile, simplify = FALSE)

names(dccData) = gsub('.*\\//|.dcc','',names(dccData))

colnames(SamplesMetadata)
Sample_ID_NTC = SamplesMetadata[which(SamplesMetadata$`slide name` == 'No Template Control'),"Sample_ID"]

string = gsub('A01|DSP-', '', Sample_ID_NTC$Sample_ID)
number=3
string = substr(string, nchar(string) - number + 1, nchar(string))

for(n in string)
{
  # print(grep(n, names(dccData)))
  
  selected = grep(n, names(dccData))
  
  print(names(dccData)[which(names(dccData)[selected] %in% Sample_ID_NTC$Sample_ID)])
  
  t=1
  for(i in selected)
  {
    # print(dim(dccData[[i]][[4]]))
    temp = dccData[[i]][[4]]
    
    #temp = temp[1:50,]
    
    colnames(temp)[1] = gsub(paste0('.*', n), '', names(dccData)[i])
    colnames(temp)[2] = paste0('Count_', gsub(paste0('.*', n), '', names(dccData)[i]) )
    
    if(names(dccData)[i] %in% c(highNTC, Sample_ID_NTC$Sample_ID)) # if its NTC has a problem
    {
      #print(names(dccData)[i])
      #if(t==1){RESULT = temp; t=t+1}else{RESULT = cbind(RESULT, temp)}
      
      if(t==1){RESULT = temp; t=t+1}else{
        
        temp2 = merge(RESULT, temp, by = 'row.names')
        
        row.names(temp2) = temp2$Row.names
        
        toMatch <- c('Count_')
        matches <- unique(grep(paste(toMatch, collapse="|"), colnames(temp2), value=F))
        temp2 = temp2[,matches] 
        temp2 = temp2[order(temp2[,1], decreasing = T),]
        
        fwrite(temp2, paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/NTC/NTC_', n, '_', colnames(temp2)[2], '.txt'), quote = F, sep = '\t', row.names = T)
        
        print(dim(temp2))
      }
    }
  }
}

#============================================


