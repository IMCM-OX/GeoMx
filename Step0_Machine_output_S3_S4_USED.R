# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html#1_Introduction
# https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html

# line 248 : Data <- Data[, QCResults$QCStatus == "PASS"]

#============================================ prepare annotation files
# datadir0 <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/")
# setwd(datadir0)
# 
# library(data.table)
# PD_FULL = fread(paste0('PD_FULL_Annotation.txt'), stringsAsFactors = F)
# CONTROL_FULL = fread(paste0('CONTROL_FULL_Annotation.txt'), stringsAsFactors = F)
# SELECTED_SAMPLES = fread(paste0('PD_CONTROL_SELECTED_SAMPLES.txt'), stringsAsFactors = F)
# 
# length(intersect(SELECTED_SAMPLES$Case, PD_FULL$Case))
# length(intersect(SELECTED_SAMPLES$Case, CONTROL_FULL$Case))
# 
# PD_FULL_subset = PD_FULL[which(PD_FULL$Case %in% SELECTED_SAMPLES$Case),]
# 
# fwrite(PD_FULL_subset, paste0(datadir0, 'PD_FULL_GeoMx_16n.txt'), quote = F, row.names = F, sep = '\t')
# 
# PD_FULL_subset = PD_FULL_subset[!is.na(PD_FULL_subset$clusters),]
# 
# PD_FULL_subset[which(PD_FULL_subset$`LB Braak stage` == "5 or 6"),"LB Braak stage"] = "5"
# 
# PD_FULL_subset$`LB Braak stage` = paste0('Braak', PD_FULL_subset$`LB Braak stage`)
# PD_FULL_subset$Dementia = paste0('Dementia_', PD_FULL_subset$Dementia)
# 
# PD_FULL_subset = as.data.frame(PD_FULL_subset)
# class(PD_FULL_subset)
# 
# # xtabs(~clusters+Gender,data=PD_FULL_subset), 
# INPUTCCA = cbind(xtabs(~clusters+Dementia,data=PD_FULL_subset), xtabs(~clusters+`LB Braak stage`,data=PD_FULL_subset))
# 
# rowsums = as.matrix(apply(INPUTCCA, 1, sum))
# rowsums
# 
# colsums = as.matrix(apply(INPUTCCA, 2, sum))
# t(colsums)
# 
# HCexp = rowsums %*%t (colsums) / sum(colsums)
# 
# library("vcd")
# mosaicplot(INPUTCCA+1, shade=TRUE, las=1, type="pearson", cex.axis=0.7, main="")
# 
# HC = as.data.frame.matrix(INPUTCCA+1)
# coaHC = dudi.coa(HC,scannf=FALSE,nf=2)
# round(coaHC$eig[1:3]/sum(coaHC$eig)*100)
# 
# fviz_ca_biplot(coaHC, repel=TRUE, col.col="brown", col.row="purple") +
#   ggtitle("") + ylim(c(-0.5,0.5))


#============================================  
library(GeomxTools)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(tidyverse)

## Introduction
# This summarises the QC of the GeoMx pilot run of 10th October, utilising the 'GeoMxTools' pipeline.

#============================================ Collate files
datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/")
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

# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

#============================================ this is before aggregation and it is not useable
# library(dplyr)
# library(ggforce)
# count_mat = data.frame(GeneSymbol = Data@featureData@data$TargetName, Data@assayData$exprs)
# 
# library(data.table)
# fwrite(count_mat, paste0(datadir, 'GeoMx_readCount.txt'), quote = F, row.names = F, sep = '\t')

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
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 20000,     # for WGT is 10K. Maximum counts observed in NTC well (1000)
       minNuclei = 5,         # Minimum # of nuclei estimated (100)
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
#============================================ read in metadata for all samples

SamplesMetadata = read_excel(paste0(datadir, 'annotation/Annotation_File.xlsx'))
dim(SamplesMetadata)
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

#============================================
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

# NTC
kable(table(State = SamplesMetadata$segment[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[sData_df$QCFlags$HighNTC] )]))
kable(table(Batch = SamplesMetadata$batch[which(SamplesMetadata$Sample_ID %in% new_sData_df$SampleID[sData_df$QCFlags$HighNTC] )]))

#============================================ summary of samples
# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols

df_pData = pData(Data)
df_pData$Donor = gsub(' .*', '', df_pData$`slide name`)

table(df_pData$Donor)
table(df_pData$segment)

#============================================ save QC results

sData_df <- sData(Data)[, dimnames(sData(Data))[[2]]]
dim(sData_df)

colnames(sData_df)
summary(sData_df$NucleiCount[sData_df$segment == 'control'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-with-LB'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-without-LB'])

library(tidyr)
new_sData_df <- unnest(sData_df, QCFlags)
new_sData_df <- unnest(new_sData_df, colnames(new_sData_df)[grep('%',colnames(new_sData_df))])
new_sData_df <- unnest(new_sData_df, NegGeoMean)

new_sData_df$RUN = gsub('.*Run ', '', new_sData_df$`scan name`)
table(new_sData_df$RUN)
table(new_sData_df$segment)
table(new_sData_df$segment[which(new_sData_df$RUN == 1)])

#new_sData_df = new_sData_df[-grep('NegGeoMean',colnames(new_sData_df)),]
library(data.table)

new_sData_df$SampleID
QCResults$SampleID = gsub('.dcc','',row.names(QCResults))

new_sData_df = merge(new_sData_df, QCResults[,c('SampleID', 'QCStatus')], by = 'SampleID')

fwrite(new_sData_df, paste0(datadir, 'QC_All.txt'), quote = F, row.names = F, sep = '\t')
dim(new_sData_df)

#============================================ summary of metadata
SamplesMetadata = SamplesMetadata[which(SamplesMetadata$Sample_ID %in% sData_df$SampleID),]

table(SamplesMetadata$segment)
length(unique(SamplesMetadata$Sample_ID))
length(unique(gsub(" .*", "",SamplesMetadata$`slide name`)))-1 # -1 is for NO
length(grep('PD', unique(gsub(" .*", "",SamplesMetadata$`slide name`)))) # -1 is for NO

table(SamplesMetadata$batch, SamplesMetadata$segment)

dim(Data[, QCResults$QCStatus == "PASS"])

#============================================ removeQCSampleProbe, eval = TRUE-----------------------------------------
Data <- Data[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(Data)

#============================================

#---- metadata for QC passed samples
sData_df <- sData(Data)[, dimnames(sData(Data))[[2]]]
dim(sData_df)

sData_df$RUN = gsub('.*Run ', '', sData_df$`scan name`)

table(sData_df$RUN)
table(sData_df$segment)
table(sData_df$segment[which(sData_df$RUN == 1)])
table(sData_df$segment[which(sData_df$RUN == 2)])

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
dim(LOQ_Mat)

#============================================ Segment Gene Detection
# We first filter out segments with exceptionally low signal. These segments will have a small fraction of panel genes detected above the LOQ relative to the other segments in the study.

#---- test the GenesDetected
#View(target_Data@assayData$exprs)
nGenesPerCell <- apply(target_Data@assayData$exprs, 2, function(x) sum(x>1))

identical(names(nGenesPerCell), row.names(pData(target_Data)))

#View(data.frame(original = pData(target_Data)$GenesDetected, extimate = nGenesPerCell))

library(dplyr)
library(ggforce)
count_mat = data.frame(GeneSymbol = target_Data@featureData@data$TargetName, target_Data@assayData$exprs)
dim(count_mat)
library(data.table)
fwrite(count_mat, paste0(datadir, 'GeoMx_readCount.txt'), quote = F, row.names = F, sep = '\t')

#---- 

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

#============================================
# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_Data)$DetectionThreshold, pData(target_Data)$segment))
#View(pData(target_Data))

#============================================ save QC incluidng gene detection rate [filters out segments]
dimnames(pData(target_Data))[[2]]
target_Data_df <- pData(target_Data)[, c("GenesDetected", "GeneDetectionRate", "DetectionThreshold")]

target_Data_df$SampleID = gsub('.dcc', '', row.names(target_Data_df))
QC_Complete = merge(target_Data_df, new_sData_df, by = 'SampleID')
dim(QC_Complete)

QC_Complete$GeneDetectionFail = rep('NA', dim(QC_Complete)[1])
# QC_Complete$GeneDetectionFail[which(QC_Complete$DetectionThreshold %in% c('<1%'))] = 'TRUE'
# QC_Complete$GeneDetectionFail[-which(QC_Complete$DetectionThreshold %in% c('<1%'))] = 'FALSE'
QC_Complete$GeneDetectionFail[which(QC_Complete$DetectionThreshold %in% c('<1%', '1-5%'))] = 'TRUE'
QC_Complete$GeneDetectionFail[-which(QC_Complete$DetectionThreshold %in% c('<1%', '1-5%'))] = 'FALSE'
table(QC_Complete$GeneDetectionFail)
table(QC_Complete$RUN)

table(QC_Complete$segment)
table(QC_Complete$segment[which(QC_Complete$RUN == 1)])
table(QC_Complete$segment[which(QC_Complete$RUN == 2)])

QC_Complete = merge(QC_Complete, new_sData_df[,c('SampleID', 'QCStatus')], by = 'SampleID')
kable(table(QC_Complete$`slide name`[which(QC_Complete$GeneDetectionFail != "TRUE")]))

getwd()
fwrite(QC_Complete, paste0(datadir, 'QC_PASS_inc_GeneDetection.txt'), quote = F, row.names = F, sep = '\t')

dim(QC_Complete)
kable(table(QC_Complete$`slide name`[which(QC_Complete$GeneDetectionFail != "TRUE" & QC_Complete$NTC < 8000)]))

#--- keep QC passed segments
#rm(list = c('Data_backup', 'Data', 'QCResults_backup', 'QCResults'))
# target_Data
# target_Data = target_Data[which(QC_Complete$GeneDetectionFail != "TRUE" & QC_Complete$QCStatus == "PASS")]
#---

#============================================ Gene Detection Rate
library(scales) # for percent

# Next, we determine the detection rate for genes across the study. To illustrate this idea, we create a small gene list (goi) to review.
# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_Data)]
fData(target_Data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_Data)$DetectionRate <-
  fData(target_Data)$DetectedSegments / nrow(pData(target_Data))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "SNCA")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_Data)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_Data)[goi, "DetectionRate"]))
goi_df

#============================================ Gene Filtering [filters out genes]
# We will graph the total number of genes detected in different percentages of segments. Based on the visualization below, we can better understand global gene detection in our study and select how many low detected genes to filter out of the dataset. Gene filtering increases performance of downstream statistical tests and improves interpretation of true biological signal.
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))

plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_Data)$DetectionRate >= x)}))

plot_detect$Rate <- plot_detect$Number / nrow(fData(target_Data))

rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_Data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_Data <- 
  target_Data[fData(target_Data)$DetectionRate >= 0.1 |
                    fData(target_Data)$TargetName %in% neg_probes, ]
dim(target_Data)
# Features  Samples 
# 12927      274 

# retain only detected genes of interest
# goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
#          "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
# goi <- goi[goi %in% rownames(target_Data)]

#============================================ 
# library(reshape2)  # for melt
# library(cowplot)   # for plot_grid
# 
# # Graph Q3 value vs negGeoMean of Negatives
# ann_of_interest <- "segment"
# Stat_data <- 
#   data.frame(row.names = colnames(exprs(target_Data)),
#              Segment = colnames(exprs(target_Data)),
#              Annotation = pData(target_Data)[, ann_of_interest],
#              Q3 = unlist(apply(exprs(target_Data), 2,
#                                quantile, 0.75, na.rm = TRUE)),
#              NegProbe = exprs(target_Data)[neg_probes, ])
# 
# Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
#                     variable.name = "Statistic", value.name = "Value")
# 
# plt1 <- ggplot(Stat_data_m,
#                aes(x = Value, fill = Statistic)) +
#   geom_histogram(bins = 40) + theme_bw() +
#   scale_x_continuous(trans = "log2") +
#   facet_wrap(~Annotation, nrow = 1) + 
#   scale_fill_brewer(palette = 3, type = "qual") +
#   labs(x = "Counts", y = "Segments, #")
# 
# plt2 <- ggplot(Stat_data,
#                aes(x = NegProbe, y = Q3, color = Annotation)) +
#   geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
#   geom_point() + guides(color = "none") + theme_bw() +
#   scale_x_continuous(trans = "log2") + 
#   scale_y_continuous(trans = "log2") +
#   theme(aspect.ratio = 1) +
#   labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")
# 
# plt3 <- ggplot(Stat_data,
#                aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
#   geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
#   geom_point() + theme_bw() +
#   scale_x_continuous(trans = "log2") + 
#   scale_y_continuous(trans = "log2") +
#   theme(aspect.ratio = 1) +
#   labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
# 
# btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
#                      rel_widths = c(0.43,0.57))
# plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

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
kable(table(pData(target_Data)$DetectionThreshold, pData(target_Data)$segment))

#============================================
toMatch <- gsub('.dcc','',row.names(pData(target_Data))[which(pData(target_Data)$DetectionThreshold %in% c('<1%', '1-5%'))])

table(pData(target_Data)$DetectionThreshold %in% c('<1%', '1-5%')) # , '1-5%'
kable(table(Batch = SamplesMetadata$batch[which(SamplesMetadata$Sample_ID %in% toMatch )]))

#============================================ filter Segments based on gene detection rate [total number of genes is around 12k]
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

table(QC_Complete$segment)
# control    PD-with-LB PD-without-LB 
# 254           190           184 

df_pData = pData(target_Data)
df_pData$Run = gsub('.*Run ', '', df_pData$`scan name`)
df_pData$Donor = gsub(' .*', '', df_pData$`slide name`)

table(df_pData$Donor)
table(df_pData$segment)
table(df_pData$segment[which(df_pData$Run == 1)])

dim(pData(target_Data))

sData_df <- sData(target_Data)[, dimnames(sData(target_Data))[[2]]]
dim(sData_df)

summary(sData_df$NucleiCount[sData_df$segment == 'control'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-with-LB'])
summary(sData_df$NucleiCount[sData_df$segment == 'PD-without-LB'])

#------------ Convert NanoStringGeoMxSet to SpatialExperiment
# The three errors that can occur when trying to coerce to Seurat are:
#   
# object must be on the target level
# object should be normalized, if you want raw data you can set forceRaw to TRUE
# normalized count matrix name must be valid


#============================================
# target_Data <- normalize(target_Data, norm_method="quant", desiredQuantile = .75, toElt = "q_norm")
# target_Data <- normalize(target_Data, norm_method = "neg", fromElt = "exprs", toElt = "neg_norm")
# exprs(target_Data)[1:5, 1:2]
# 
# data.frame(assayData(target_Data)[["exprs"]][seq_len(3), seq_len(3)])
# featureType(target_Data)
# data.frame(assayData(target_Data)[["q_norm"]][seq_len(3), seq_len(3)])

# library(stringr)
# string = str_split_fixed(gsub('.dcc','',colnames(exprs(target_Data))), "-", 4)

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

#============================================ check the expression of gene to interest 
# to run the dittoPlot I need a SingleCellExperiment object

# target_Data <- normalize(target_Data, norm_method="quant", desiredQuantile = .75, toElt = "q_norm")
# target_Data <- normalize(target_Data, norm_method = "neg", fromElt = "exprs", toElt = "neg_norm")
# exprs(target_Data)[1:5, 1:2]
# 
# data.frame(assayData(target_Data)[["exprs"]][seq_len(3), seq_len(3)])
# featureType(target_Data)
# data.frame(assayData(target_Data)[["q_norm"]][seq_len(3), seq_len(3)])

# options(Seurat.object.assay.version = "v3")
# target_Data_Seurat <- as.SpatialExperiment(target_Data, normData = "exprs", forceRaw = TRUE)
# 
# counts <- SummarizedExperiment::assay(target_Data_Seurat)
# meta <- SummarizedExperiment::colData(target_Data_Seurat) %>% as.data.frame
# 
# library(Seurat)
# target_Data_Seurat <- CreateSeuratObject(counts = counts, project = "Seurat",  meta.data = meta, assay = "GeoMx")
# 
# head(target_Data_Seurat, 3)
# target_Data_Seurat$segment
# 
# library(scater)
# sce <- as.SingleCellExperiment(target_Data_Seurat)
# # sce_mat <- assay(pbmc.sce)
# # sce_mat = as.matrix(sce_mat)
# # sce_mat[1:5,1:5]
# 
# sce <- logNormCounts(sce)
# 
# library(dittoSeq)
# dittoPlot(sce, "PEG10", group.by = "segment")

#============================================
library(standR)
library(SpatialExperiment)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

#============================================ generate SpatialExperiment
spe <- as.SpatialExperiment(target_Data, normData = "exprs", forceRaw = T)
spe

assayNames(spe)

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
dim(MD)

metadata(spe)$NegProbes = negativeGeoMeans

#============================================ Check QCFlags
# In the meta data file generated by Nanostring, there is a column called “QCFlags”, which indicates bad quality tissue samples in their preliminary QC step. If the data is fine from their QC, you will see NA/empty cells in this column.
# colData(spe)$QCFlags 
# protocolData(target_Data)[["QCFlags"]]

#============================================ Sample level QC
# To visualise sample metadata
library(ggplot2)
library(ggalluvial)
# names(colData(spe2))
# spe$`slide name` = 'SN'
# #plotSampleInfo(spe, column2plot = c("slide name", "segment", "roi"))
# 
# length(SummarizedExperiment::assayNames(spe))

#============================================ Gene level QC
# sample_fraction	Double. Genes with low count in more than this threshold of the samples will be removed. Default is 0.9
# min_count Integer. Minimum read count to calculate count threshold. Default is 5.

# addPerROIQC uses 'logcounts' - 
# The sample_fraction parameter takes a value between 0 and 1, representing the fraction of sequencing reads to be used for QC calculations. For example, sample_fraction = 0.8 means that 80% of the reads in each ROI will be considered.
# A higher sample_fraction retains more information but may be more sensitive to outliers.

spe2_subset <- addPerROIQC(spe, rm_genes = TRUE, min_count = 3, sample_fraction = 0.95)
dim(spe2_subset)
metadata(spe2_subset) |> names()

# table(pData(target_Data)$segment)[1]/sum(table(pData(target_Data)$segment))

colData(spe2_subset)$Names = gsub('.dcc','', row.names(colData(spe2_subset)))
plotGeneQC(spe2_subset, ordannots = "Names", col = Names, point_size = 2)

# ref: https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html

library(stringr)
string = str_split_fixed(gsub('.dcc','',colData(spe2_subset)$Names), "-", 4)
colData(spe2_subset)$Names = string[,4]

colData(spe2_subset)$state = factor(colData(spe2_subset)$state, levels =  c("PD-with-LB", "PD-without-LB", "control"))

#============================================ Relative log expression distribution

# plotRLExpr(spe2_subset, ordannots = "state", assay = 2, color = state)
# plotRLExpr(spe2_subset, ordannots = "state", assay = 2, color = area)
# 
# drawPCA(spe2_subset, assay = 2, color = segment)
# 
# set.seed(100)
# spe2_subset <- scater::runPCA(spe2_subset)
# pca_results <- reducedDim(spe2_subset, "PCA")
# drawPCA(spe2_subset, precomputed = pca_results, col = segment)
# plotScreePCA(spe2_subset, precomputed = pca_results)
# 
# standR::plotMDS(spe2_subset, assay = 2, color = segment)

#============================================ check integrity of metadata
# View(as.data.frame(colData(spe_ruv)))
# colData(spe_ruv)$slide.name
# colData(spe_ruv)$scan.name
# colData(spe_ruv)$segment
# colnames(colData(spe_ruv))
# 
table(colData(spe2_subset)$segment)
# control    PD-with-LB PD-without-LB 
# 235           119           166 

samples = unique(gsub(' Run.*|_Run.*| Control.*|_.*','', colData(spe2_subset)$`slide name`))
samples = unique(gsub('PD ', 'PD', samples))
length(samples)

length(unique(colData(spe2_subset)$`slide name`))
length(unique(gsub(' Run.*|_Run.*| Control.*|_.*','', colData(spe2_subset)$`slide name`)))

length(unique(gsub(' Run.*|_Run.*| Control.*|_.*','', SamplesMetadata$`slide name`)))
setdiff(unique(gsub(' Run.*|_Run.*| Control.*|_.*','', SamplesMetadata$`slide name`)), samples)

SamplesMetadata = read_excel(paste0(datadir, 'annotation/Annotation_File.xlsx'))
dim(SamplesMetadata)

table(SamplesMetadata$`slide name` == SamplesMetadata$`scan name`)

table(SamplesMetadata$segment[grep('wo', SamplesMetadata$roi)])
table(SamplesMetadata$segment[grep('w LB', SamplesMetadata$roi)])

length(unique(gsub(' Run.*|_Run.*| Control.*|_.*','', SamplesMetadata$`slide name`)))



#============================================ Normalization

#------------------ I should select the Run ------------------
Run = 'Run 1'
spe2_subset <- spe2_subset[, grep(Run, colData(spe2_subset)$`scan name`)]
dim(spe2_subset)

#============================================ I keep all genes

# spe2_subset_Run <- spe2[, grep(Run, colData(spe2)$scan.name)]
# 
# library(dplyr)
# library(ggforce)
# count_mat = data.frame(GeneSymbol = rownames(spe2_subset_Run), counts(spe2_subset_Run))
# fwrite(count_mat, paste0(datadir, 'GeoMx_readCount_',gsub(' ','',Run),'.txt'), quote = F, row.names = F, sep = '\t')
# dim(count_mat)
#============================================

#-------------------------------------- summary of segments
table(colData(spe2_subset)$segment)
# Run 1
# control    PD-with-LB PD-without-LB 
# 112            52            86 
# Run 2
# control    PD-with-LB PD-without-LB 
# 123            67            80 

samples = unique(gsub(' Run.*|_Run.*| Control.*|_.*','', colData(spe2_subset)$`slide name`))
samples = unique(gsub('PD ', 'PD', samples))
length(samples)

SamplesMetadata = read_excel(paste0(datadir, 'annotation/Annotation_File.xlsx'))
dim(SamplesMetadata)
SamplesMetadata$Run = gsub('.*Run ', '', SamplesMetadata$`scan name`)

SamplesMetadata = SamplesMetadata[!is.na(SamplesMetadata$Run),]
dim(SamplesMetadata)

table(SamplesMetadata$Run)
table(SamplesMetadata$segment)
table(SamplesMetadata$segment[which(SamplesMetadata$segment == "PD-with-LB")])
table(SamplesMetadata$segment[which(SamplesMetadata$Run == "1")])
table(SamplesMetadata$segment[which(SamplesMetadata$Run == "2")])

#--------------------------------------
spe2_subset@assays
spe2_subset@assays@data@listData$logcounts = NULL

spe_tmm1 <- geomxNorm(spe2_subset, method = "TMM")
# I test to be sure that I save TMM in logcounts
spe_tmm1@assays # logcounts [assay = 2] contain the TMM

#============================================ Import from DGEList object [this is essential for running drawPCA]
# Alternatively, standR provides a function to generate a spatial experiment object from a DGEList object, which would be useful for users who used edgeR package and have existing analyses and implementations using DGEList objects to port across to the standR workflow.
dge <- edgeR::SE2DGEList(spe2_subset)

spe_tmm <- readGeoMxFromDGE(dge)
spe_tmm

# the functin revises the fromat but changes the content of logcounts, so I need to revise it.
assay(spe_tmm, 'logcounts') = assay(spe_tmm1, 'logcounts')

# assay(spe_tmm1, 'logcounts')[1:5,1:5]
# assay(spe_tmm, 'logcounts')[1:5,1:5]

#============================================ 
plotRLExpr(spe_tmm, assay = 2, color = segment) + ggtitle("TMM")
drawPCA(spe_tmm, assay = 2, color = segment)

#============================================ Batch correction

S3Metadata = read_excel('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3/annotation/Annotation_File.xlsx')
dim(S3Metadata)

#--- grep multiple pattern [removed accoroding heatmap]
# toMatch <- c("PD1121", "PD1141", "PD1151")
# # OR
# matches <- unique(grep(paste(toMatch,collapse="|"), S3Metadata$`slide name`, value=TRUE))

S4Metadata = read_excel('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S4/annotation/Annotation_File.xlsx')
dim(S4Metadata)

row.names(colData(spe_tmm)) = gsub('.dcc', '', row.names(colData(spe_tmm)))

colData(spe_tmm)$batch = 1:dim(colData(spe_tmm))[1]
colData(spe_tmm)$batch[which(row.names(colData(spe_tmm)) %in% S3Metadata$Sample_ID)] = 'S3'
colData(spe_tmm)$batch[which(row.names(colData(spe_tmm)) %in% S4Metadata$Sample_ID)] = 'S4'

table(colData(spe_tmm)$batch)

#================================================== combined batch
slide_dateMetadata = read_excel('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k//annotation2//Metatdata_S3S4_with_slide_date.xlsx')
dim(slide_dateMetadata)

slide_dateMetadata = as.data.frame(slide_dateMetadata)
slide_dateMetadata = slide_dateMetadata[which(slide_dateMetadata$Sample_ID %in% row.names(colData(spe_tmm))),]
dim(slide_dateMetadata)
length(row.names(colData(spe_tmm)))

slide_dateMetadata = slide_dateMetadata[match(row.names(colData(spe_tmm)),slide_dateMetadata$Sample_ID),]
identical(row.names(colData(spe_tmm)),slide_dateMetadata$Sample_ID)

table(slide_dateMetadata$Batch_slides)
colData(spe_tmm)$batchslide = slide_dateMetadata$Batch_slides
table(colData(spe_tmm)$batchslide)

colData(spe_tmm)$batch_combined = paste0(colData(spe_tmm)$batch, '_', colData(spe_tmm)$batchslide)

library(openxlsx)
write.xlsx(as.data.frame(colData(spe_tmm)), '/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/annotation2/Metatdata_S3S4_with_batch_combined.xlsx')


#==================================================

# Always combine your relevant batch columns into a single CombinedBatch factor before running findNCGs. This provides the most comprehensive and accurate definition of what constitutes a "batch" for the purpose of identifying negative control genes, which are crucial for subsequent normalization and batch correction steps in standR.

spe_tmm <- findNCGs(spe_tmm, batch_name = "batch", top_n = 1000)
# findNCGs: This function aims to identify network community genes (NCGs), which are genes that are highly influential within their respective network communities.
# top_n: This parameter specifies the number of top-ranked NCGs to be returned for each network community.
metadata(spe_tmm) |> names()

for(i in seq(5)){
  spe_ruv <- geomxBatchCorrection(spe_tmm, factors = "segment", 
                                  NCGs = metadata(spe_tmm)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = segment, title = paste0("k = ", i)))
  
}

spe_ruv <- geomxBatchCorrection(spe_tmm, factors = "segment", NCGs = metadata(spe_tmm)$NCGs, k = 5)

print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = segment, title = paste0("k = ", 5)))
plotRLExpr(spe_ruv, assay = 2, color = segment) + ggtitle("TMM + RUV4")

# spe_ruv@assays
# spe_ruv@assays@data@listData$logcounts

count_mat = data.frame(GeneSymbol = row.names(spe_ruv@assays@data@listData$logcounts), spe_ruv@assays@data@listData$logcounts)
fwrite(count_mat, paste0(datadir, 'GeoMx_readCount_TMMRUV_R1.txt'), quote = F, row.names = F, sep = '\t')

format(object.size(spe2_subset), units = "Mb")
#rm(list = c('spe_tmm', 'spe2_subset'))

#============================================ compare different types of normalization
# library(edgeR)
# library(limma)
# 
# dge <- SE2DGEList(spe_ruv)
# design <- model.matrix(~0 + segment + ruv_W1 + ruv_W2, data = colData(spe_ruv)) # Adding ruv_W1 + ruv_W2 reduces the number of DEGs by up to 3.
# colnames(design)
# 
# v <- list()
# # 1. Normalizing only to library size:
# v$libsize <- voom(dge, design)
# # 2. TMM normalization:
# dge <- calcNormFactors(dge)
# v$tmm <- voom(dge, design)
# # 3. Cyclic loess normalization:
# v$cycloess <- voom(dge, design, normalize="cyclicloess")
# # 4. Quantile normalization:
# v$quantile <- voom(dge, design, normalize="quantile")
# 
# names(spe_ruv@assays)
# spe_ruv@assays@data@listData$LimmaTMM = v$tmm$E
# spe_ruv@assays@data@listData$Limmacyclicloess = v$cycloess$E
# spe_ruv@assays@data@listData$LimmaQuantile = v$quantile$E
# spe_ruv@assays@data@listData$libsize = v$libsize$E

# print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = segment, title = "TMM + RUV4"))
# print(plotPairPCA(spe_ruv, assay = 3, n_dimension = 4, color = segment, title = "Limma TMM"))
# print(plotPairPCA(spe_ruv, assay = 4, n_dimension = 4, color = segment, title = "Cyclic Loess"))
# print(plotPairPCA(spe_ruv, assay = 5, n_dimension = 4, color = segment, title = "Quantile"))
# print(plotPairPCA(spe_ruv, assay = 6, n_dimension = 4, color = segment, title = "Library size"))
# 
# plotRLExpr(spe_ruv, assay = 2, color = segment) + ggtitle("TMM + RUV4")
# plotRLExpr(spe_ruv, assay = 3, color = segment) + ggtitle("Limma TMM")
# plotRLExpr(spe_ruv, assay = 4, color = segment) + ggtitle("Cyclic Loess")
# plotRLExpr(spe_ruv, assay = 5, color = segment) + ggtitle("Quantile")
# plotRLExpr(spe_ruv, assay = 6, color = segment) + ggtitle("Library size")

## -----------------------------------------------------------------------------

# names = data.frame(names1 = names(spe_ruv@assays),
#                    names2 = c('Counts','TMM + RUV4','Limma TMM','Cyclic Loess','Quantile','Library size'))
# 
# i=1
# for(i in 1:length(names(spe_ruv@assays)) )
# {
#   name = names(spe_ruv@assays)[i]
#   print(names$names2[names$names1 == name])
#   setwd(datadir)
#   
#   pdf(paste0(name, "_", Run, "_plotPairPCA.pdf"), width = 10, height = 10)
#   
#   library(factoextra)
#   library(tidyverse)
#   identical(colnames(spe_ruv@assays@data@listData$logcounts), row.names(colData(spe_ruv)))
#   res.pca <- prcomp(t(spe_ruv@assays@data@listData[[name]]),  scale = TRUE)
#   p = fviz_pca_biplot(res.pca,
#                       select.var = list(cos2 = 5),
#                       # Individuals
#                       geom.ind = "point",
#                       fill.ind = colData(spe_ruv)$segment, col.ind = "black",
#                       pointshape = 21, pointsize = 4,
#                       palette = "jco",
#                       addEllipses = TRUE,
#                       # Variables
#                       alpha.var ="contrib", col.var = "contrib",
#                       gradient.cols = "RdYlBu",
#                       legend.title = list(fill = "Species", color = "Contrib",
#                                           alpha = "Contrib")
#   ) + ggtitle(names$names2[names$names1 == name]) + coord_fixed() + theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"), legend.text=element_text(size=20), legend.title=element_text(size=20), plot.title = element_text(size = 30, face="bold"))
#   
#   print(p)
#   dev.off()
# }

#============================================ UMAP  

library(ggspavis)
library(Seurat)
## Set seed
set.seed(987)

## Compute UMAP on top 50 PCs
spe_ruv <- scater::runPCA(spe_ruv)
spe_ruv <- scater::runUMAP(spe_ruv, dimred = "PCA")

reducedDimNames(spe_ruv)
dim(reducedDim(spe_ruv, "UMAP"))

## Update column names for easier plotting
colnames(reducedDim(spe_ruv, "UMAP")) <- paste0("UMAP", 1:2)

PCA = ggplot(data = as.data.frame(spe_ruv@int_colData@listData$reducedDims$PCA),
       aes(x = PC1, y = PC2, colour = spe_ruv@colData$segment)) + 
  geom_point(size = 5) + 
  scale_colour_brewer(type = "qual") + 
  labs(title = "Reduced dimensions: PCA",
       x = "PC1",
       y = "PC2",
       colour = "Layers") +
  theme_classic()

getwd()
pdf(paste0(gsub(' ','',Run), "_TMM_RUV_PCA.pdf"), width = 8, height = 8)
PCA
dev.off()
PCA

UMAP = ggplot(data = as.data.frame(spe_ruv@int_colData@listData$reducedDims$UMAP),
       aes(x = UMAP1, y = UMAP2, colour = spe_ruv@colData$segment)) + 
  geom_point(size = 5) + 
  scale_colour_brewer(type = "qual") + 
  labs(title = "Reduced dimensions: UMAP",
       x = "UMAP1",
       y = "UMAP2",
       colour = "Layers") +
  theme_classic()
UMAP
getwd()
pdf(paste0(gsub(' ','',Run), "_TMM_RUV_UMAP.pdf"), width = 8, height = 8)
UMAP
dev.off()

#---- visualization of expression level per gene

# meta.data = data.frame(slide.name = colData(spe_ruv)$slide.name, state = colData(spe_ruv)$state)
# row.names(meta.data) = row.names(colData(spe_ruv))
# 
# meta.data$slide.name[which(meta.data$state == 'PD-with-LB')] = paste0('V', meta.data$slide.name[which(meta.data$state == 'PD-with-LB')])
# meta.data$slide.name[which(meta.data$state == 'PD-without-LB')] = paste0('R', meta.data$slide.name[which(meta.data$state == 'PD-without-LB')])
# 
# library(Seurat)
# target_Data_Seurat <- CreateSeuratObject(counts = logcounts(spe_ruv), project = "Seurat",  meta.data = meta.data, assay = "GeoMx")
# 
# library(scater)
# sce <- as.SingleCellExperiment(target_Data_Seurat)
# # assay(sce) # uses the normalized values
# 
# query = 'KCNJ6'
# query = 'CALB1'
# query = "SNCA"
# # ENO1
# 
# library(dittoSeq)
# dittoPlot(sce, query, group.by = "state")
# 
# dittoPlot(sce, query, group.by = "state",
#           plots = c("vlnplot", "jitter", "boxplot"),
#           # change the color and size of jitter points
#           jitter.color = "blue", jitter.size = 0.7,
#           # change the outline color and width, and remove the fill of boxplots
#           boxplot.color = "white", boxplot.width = 0.1,
#           boxplot.fill = FALSE,
#           # change how the violinplot widths are normalized across groups
#           vlnplot.scaling = "count"
# )

#============================================ UMAP _ tSNE [does not use TMM + RUV4]
# library(umap)
# library(Rtsne)
# 
# # update defaults for umap to contain a stable random_state (seed)
# custom_umap <- umap::umap.defaults
# custom_umap$random_state <- 42
# # run UMAP
# umap_out <-
#   umap(t(log2(assayDataElement(target_Data , elt = "q_norm"))),  
#        config = custom_umap)
# #> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# #> Also defined by 'spam'
# pData(target_Data)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
# 
# colnames(pData(target_Data))
# 
# kable(table(pData(target_Data)$DetectionThreshold, pData(target_Data)))
# 
# ggplot(pData(target_Data),
#        aes(x = UMAP1, y = UMAP2, color = area)) +
#   geom_point(size = 3) +
#   theme_bw()
# 
# # run tSNE
# set.seed(42) # set the seed for tSNE as well
# tsne_out <-
#   Rtsne(t(log2(assayDataElement(target_Data , elt = "q_norm"))),
#         perplexity = ncol(target_Data)*.15)
# pData(target_Data)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
# ggplot(pData(target_Data),
#        aes(x = tSNE1, y = tSNE2, color = segment)) +
#   geom_point(size = 3) +
#   theme_bw()


#============================================ celltype deconvolution

datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/")
setwd(datadir)

#---------- countFile
library(data.table)
count_matrix = fread(paste0(datadir, 'GeoMx_readCount_TMMRUV_R1.txt'), stringsAsFactors = F, header = T)
count_matrix = as.data.frame(count_matrix)
colnames(count_matrix)[colnames(count_matrix) == 'GeneSymbol'] = 'TargetName'

#- add negativeProbef
negativeProbefData <- subset(fData(target_Data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

count_neg_probes = data.frame(assayData(target_Data)$exprs[grep(neg_probes, row.names(assayData(target_Data)$exprs)),])
count_neg_probes = t(count_neg_probes)
colnames(count_neg_probes) = gsub('.dcc','',colnames(count_neg_probes))
colnames(count_neg_probes) = gsub('-','.',colnames(count_neg_probes))
count_neg_probes = cbind(TargetName = neg_probes, count_neg_probes)
row.names(count_neg_probes) = NULL
count_neg_probes = count_neg_probes[,which(colnames(count_neg_probes) %in% colnames(count_matrix))]
length(count_neg_probes)
dim(count_matrix)

count_neg_probes = count_neg_probes[match(colnames(count_matrix), names(count_neg_probes))]
identical(names(count_neg_probes), colnames(count_matrix))

# add Negative probes to the count matrix
count_matrix = rbind(count_neg_probes, count_matrix)
count_matrix[,-1] <- as.data.frame(sapply(count_matrix[,-1], as.numeric)) #<- sapply is here

#---------- featureAnnoFile
featureAnnotationFile = data.frame(TargetName = count_matrix$TargetName)

#---------- sampleAnnoFile
sampleAnnotationFile = fread(paste0(datadir, 'QC_All.txt'), stringsAsFactors = F, header = T)
sampleAnnotationFile = as.data.frame(sampleAnnotationFile)
# "slide name" "SlideName"
# "scan name"  "ScanLabel"
# "segment" "SegmentLabel"
# "SampleID" "SegmentDisplayName"
# "segment" "Type"

sampleAnnotationFile$SlideName = sampleAnnotationFile$`slide name`
sampleAnnotationFile$ScanLabel = sampleAnnotationFile$`scan name`
sampleAnnotationFile$ROILabel = sampleAnnotationFile$Well
sampleAnnotationFile$SegmentLabel = sampleAnnotationFile$segment
sampleAnnotationFile$RawReads = sampleAnnotationFile$Raw
sampleAnnotationFile$Type = sampleAnnotationFile$segment
sampleAnnotationFile$SegmentDisplayName = gsub('-','.', sampleAnnotationFile$SampleID)
sampleAnnotationFile$ROICoordinateX = 1:dim(sampleAnnotationFile)[1]
sampleAnnotationFile$ROICoordinateY = 1:dim(sampleAnnotationFile)[1]

# sampleAnnotationFile = sampleAnnotationFile[grep('Run 1', sampleAnnotationFile$`slide name`),]
# count_matrix = count_matrix[,which(gsub('\\.', '-', colnames(count_matrix)) %in% sampleAnnotationFile$SampleID)]
# dim(count_matrix)
# sampleAnnotationFile = sampleAnnotationFile[which(sampleAnnotationFile$SampleID %in% gsub('\\.', '-', colnames(count_matrix))),]
# dim(sampleAnnotationFile)

#-------------------------
library(SpatialDecon)
library(standR)
spe2 <- readGeoMx(count_matrix, sampleAnnotationFile, featureAnnotationFile, rmNegProbe = FALSE, NegProbeName = "NegProbe-WTX",
                  colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
                  coord.colnames = c("ROICoordinateX", "ROICoordinateY"))

spd <- prepareSpatialDecon(spe2)
names(spd)

#------------------------- scRANseq Cell-Types -------------------------

PI='GSE243639'
library(data.table)
#--- MetaData
MetaData = fread('/Users/isarnassiri/Documents/Parkinson/Data/GSE243639_RAW/GSE243639_Clinical_data.csv')
MetaData = as.data.frame(MetaData)
table(MetaData$`Clinical diagnosis`)
#MetaData = MetaData[which(MetaData$`Clinical diagnosis` != "Parkinson's"),]
length(MetaData$`Sample ID`)

#--- CellType Annotations
CellTypeEnrichment = fread('/Users/isarnassiri/Documents/Parkinson/Data/GSE243639_RAW/AllCells_Annotation.txt')
CellTypeEnrichment = as.data.frame(CellTypeEnrichment)

CellTypeEnrichment$Samples = gsub('_.*', '', CellTypeEnrichment$CELL_ID)

CellTypeEnrichment_subset = CellTypeEnrichment[which(CellTypeEnrichment$Samples %in% MetaData$`Sample ID`),]
dim(CellTypeEnrichment)
dim(CellTypeEnrichment_subset)

CellTypeEnrichment_subset$CELL_ID = gsub('\\.', '-', CellTypeEnrichment_subset$CELL_ID)
CellTypeEnrichment_subset = CellTypeEnrichment_subset[-which(CellTypeEnrichment_subset$IDENT %in% c("OPC", 'VC', 'T cells')),]

CellTypeEnrichment_subset = CellTypeEnrichment_subset[!is.na(CellTypeEnrichment_subset$CELL_ID),]
dim(CellTypeEnrichment_subset)

#--- Read Count
TP_profile = fread('/Users/isarnassiri/Documents/Parkinson/Data/GSE243639_RAW/GSE243639_Filtered_count_table.csv')
TP_profile = as.data.frame(TP_profile)
dim(TP_profile)
row.names(TP_profile) = TP_profile$V1
TP_profile = TP_profile[,-1]
colnames(TP_profile) = gsub('\\.', '-', colnames(TP_profile))

#------
files = list.files('/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis', pattern = 'sc_brain_region_')
setwd('/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis')

require(data.table) ## 1.9.2 or 1.9.3
genesets = rbindlist(lapply(files, fread), idcol="ID")

genesets$ID = gsub('.*rna_|.tsv', '', files)[genesets$ID]
table(genesets$ID)
table(duplicated(genesets$Gene))

TP_profile_subsetCellType = TP_profile[which(row.names(TP_profile) %in% genesets$Gene),]
dim(TP_profile_subsetCellType)

TP_profile_subsetCellType = TP_profile_subsetCellType[, which(colnames(TP_profile_subsetCellType) %in% CellTypeEnrichment_subset$CELL_ID),]
dim(TP_profile_subsetCellType)

TP_profile_subsetCellType = TP_profile_subsetCellType[apply(TP_profile_subsetCellType, 1, function(x) !all(x==0)),]
dim(TP_profile_subsetCellType)

genesets = genesets[which(genesets$Gene %in% row.names(TP_profile_subsetCellType)),]

kable(table(Cell_Type = genesets$ID[which(genesets$Gene %in% row.names(TP_profile_subsetCellType))]))
#------

CellTypeEnrichment_subset = CellTypeEnrichment_subset[match(colnames(TP_profile_subsetCellType), CellTypeEnrichment_subset$CELL_ID),]
dim(CellTypeEnrichment_subset)
identical(colnames(TP_profile_subsetCellType), CellTypeEnrichment_subset$CELL_ID)

CellTypeEnrichment_subset$Type = rep(1, dim(CellTypeEnrichment_subset)[1])

colnames(TP_profile_subsetCellType) = CellTypeEnrichment_subset$IDENT
table(colnames(TP_profile_subsetCellType))
class(TP_profile_subsetCellType)
dim(TP_profile_subsetCellType)

patterns <- names(table(colnames(TP_profile_subsetCellType)))
scRNAseq <- sapply(patterns, function(xx) rowSums(TP_profile_subsetCellType[,grep(xx, names(TP_profile_subsetCellType)), drop=FALSE]))  # loop through

heatmap(sweep(scRNAseq, 1, apply(scRNAseq, 1, max), "/"), labRow = NA, margins = c(10, 5))

#---------- prune cells using Multiple pairwise-comparison [the proposed cell should be dominant one]

genesets = genesets[match(row.names(scRNAseq),genesets$Gene),]
identical(row.names(scRNAseq),genesets$Gene)

scRNAseq = as.data.frame(scRNAseq)
scRNAseq$sigDiff = scRNAseq$CellType = 1:length(scRNAseq$Astro)

# Using options(scipen = ...)
options(scipen = 999) # Prevent scientific notation

for(it in 1:dim(scRNAseq)[1])
{
  print(it)
  print(genesets$ID[it])
  
  if(genesets$ID[it] == "Microglia")
  {
    results = pairwise.wilcox.test(as.numeric(TP_profile_subsetCellType[it,]), factor(CellTypeEnrichment_subset$IDENT, levels = c("Micro", "Astro",  "Neurons", "Oligo")), p.adjust.method = "BH")
  }
  
  if(genesets$ID[it] == "Astrocyte")
  {
    results = pairwise.wilcox.test(as.numeric(TP_profile_subsetCellType[it,]), factor(CellTypeEnrichment_subset$IDENT, levels = c("Astro", "Micro", "Neurons", "Oligo")), p.adjust.method = "BH")
  }
  
  if(genesets$ID[it] == "Neurons")
  {
    results = pairwise.wilcox.test(as.numeric(TP_profile_subsetCellType[it,]), factor(CellTypeEnrichment_subset$IDENT, levels = c("Neurons", "Micro", "Astro",  "Oligo")), p.adjust.method = "BH")
  }
  
  if(genesets$ID[it] == "Oligodendrocyte")
  {
    results = pairwise.wilcox.test(as.numeric(TP_profile_subsetCellType[it,]), factor(CellTypeEnrichment_subset$IDENT, levels = c("Oligo", "Micro", "Astro",  "Neurons")), p.adjust.method = "BH")
  }
  
  print(results$p.value[,1])
  scRNAseq$sigDiff[it] = all(results$p.value[,1] < 0.001) # TRUE is 1
  scRNAseq$CellType[it] = genesets$ID[it]
 
}

options(scipen = 0) # Reset to default

#--- Select genes whose expression is compatible with the proposed marker gene per cell type.
#-- select cell type with max expression
scRNAseq$max <- apply(scRNAseq[,c(1:4)], 1, function(x) colnames(scRNAseq[,c(1:4)])[which.max(x)])
scRNAseq$GeneNames <- row.names(scRNAseq)

#-- rename
df = data.frame(max = c('Astro', 'Micro', 'Neurons', 'Oligo'), CellType = c('Astrocyte', 'Microglia', 'Neurons', 'Oligodendrocyte'))
scRNAseq = merge(scRNAseq, df, by = 'max')

#-- subset [keep marker genes that annotation and dominant cell type are consistence]
scRNAseq_subset = scRNAseq[which(scRNAseq$CellType.x == scRNAseq$CellType.y & scRNAseq$sigDiff == 1),]
kable(table(scRNAseq$CellType.x))
kable(table(scRNAseq_subset$CellType.x))

scRNAseq_subset = scRNAseq_subset[,-which(colnames(scRNAseq_subset) %in% c("max", "CellType.y"))]

colnames(scRNAseq_subset)[which(colnames(scRNAseq_subset) == "CellType.x")] = "CellType"
colnames(scRNAseq_subset)

table(scRNAseq$sigDiff)
kable(table(scRNAseq$max))

# scRNAseqsubset = scRNAseq[apply(scRNAseq[,c(1:4)], 1, function(x) (sum(x)>50)),]
# scRNAseqsubset = scRNAseqsubset[which(scRNAseqsubset$sigDiff == 1),]
# dim(scRNAseqsubset)
kable(table(scRNAseq_subset$CellType))

#----------

patterns <- names(table(colnames(TP_profile_subsetCellType)))
scRNAseq2 <- sapply(patterns, function(xx) rowSums(TP_profile_subsetCellType[,grep(xx, names(TP_profile_subsetCellType)), drop=FALSE]))  # loop through

scRNAseq2 = scRNAseq2[which(row.names(scRNAseq2) %in% scRNAseq_subset$GeneNames), ]
dim(scRNAseq2)

scRNAseq2 = scRNAseq2[apply(scRNAseq2[,c(1:4)], 1, function(x) (sum(x)>100)),]

kable(table(Cell_Type = genesets$ID[which(genesets$Gene %in% row.names(scRNAseq2))]))

colnames(scRNAseq2) = df$CellType[match(colnames(scRNAseq2), df$max)]

heatmap(sweep(scRNAseq2, 1, apply(scRNAseq2, 1, max), "/"), labRow = NA, margins = c(12, 5))
#dev.off()

#------------------------- 
res <- spatialdecon(norm = spd$normCount,
                    bg = spd$backGround,
                    X = scRNAseq2,
                    align_genes = TRUE)

# If I proceed with a min value higher than 5 in the function addPerROIQC (820), this function does not work and I get the following error message.
# Only 60 genes are shared between norm and X, which may not be enough to support accurate deconvolution.

#-------------------------
str(res$beta)
colnames(res$prop_of_all)
subset_prop <- res$prop_of_all

colnames(colData(spe_ruv))

#------ input for visualziation
sampleAnnotationFile_Run1 = sampleAnnotationFile[which(sampleAnnotationFile$SampleID %in% colData(spe_ruv)$SampleID ),]
dim(sampleAnnotationFile_Run1)

subset_prop_Run1 = subset_prop[,which(gsub('\\.','-',colnames(subset_prop)) %in% sampleAnnotationFile_Run1$SampleID)]

row.names(subset_prop_Run1)
class(subset_prop_Run1)
subset_prop_Run1 = subset_prop_Run1[,order(as.numeric(subset_prop_Run1[which(row.names(subset_prop_Run1) == "Neurons"),]), decreasing = T)]
subset_prop_Run1 = subset_prop_Run1[,!is.na(subset_prop_Run1[3,])]

sampleAnnotationFile_Run1 = sampleAnnotationFile_Run1[which( sampleAnnotationFile_Run1$SampleID %in%  gsub('\\.','-',colnames(subset_prop_Run1)) ),]
dim(subset_prop_Run1)
dim(sampleAnnotationFile_Run1)

length(setdiff(gsub('\\.','-',colnames(subset_prop_Run1)), sampleAnnotationFile_Run1$SampleID))

sampleAnnotationFile_Run1 = sampleAnnotationFile_Run1[match(gsub('\\.','-',colnames(subset_prop_Run1)), sampleAnnotationFile_Run1$SampleID),]
identical(gsub('\\.','-',colnames(subset_prop_Run1)), sampleAnnotationFile_Run1$SampleID)

colnames(subset_prop_Run1) = make.names(sampleAnnotationFile_Run1$segment, unique = T)
subset_prop_Run1 = as.data.frame(subset_prop_Run1)

cbind(subset_prop_Run1[,1:10], subset_prop_Run1[,(dim(subset_prop_Run1)[2]-10):(dim(subset_prop_Run1)[2])]) %>%
  as.data.frame() %>%
  rownames_to_column("CellTypes") %>%
  gather(samples, prop, -CellTypes) %>%
  ggplot(aes(samples, prop, fill = CellTypes)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

subset_prop_Run1 %>%
  as.data.frame() %>%
  rownames_to_column("CellTypes") %>%
  gather(samples, prop, -CellTypes) %>%
  ggplot(aes(samples, prop, fill = CellTypes)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")


library(dplyr)
summary_input = t(subset_prop_Run1)
summary_input = as.data.frame(summary_input)
summary_input$Type = sampleAnnotationFile_Run1$segment

dim(summary_input)

summary_result = summary_input %>%
  group_by(Type) %>%
  summarise_all(list(median))

#--- table of cell types
library(pander)
panderOptions('digits', 3)
pander(summary_result)

#============================================ Differential expression analysis ============================================

#--- add celltype proportions as covariance
subset_prop_perRun = subset_prop[,which(gsub('\\.','-', colnames(subset_prop)) %in% colData(spe_ruv)$SampleID)]
subset_prop_perRun = t(subset_prop_perRun)
subset_prop_perRun = as.data.frame(subset_prop_perRun)
row.names(subset_prop_perRun) = gsub('\\.', '-', row.names(subset_prop_perRun))

identical(row.names(subset_prop_perRun), colData(spe_ruv)$SampleID)

class(subset_prop_perRun)
colnames((subset_prop_perRun))

table(is.na(subset_prop_perRun$Astro))

subset_prop_perRun$Astr0o[is.na(subset_prop_perRun$Astro)] = 0
subset_prop_perRun$Micro[is.na(subset_prop_perRun$Micro)] = 0
subset_prop_perRun$Oligo[is.na(subset_prop_perRun$Oligo)] = 0

View(as.data.frame(colData(spe_ruv)))
colData(spe_ruv)$Astro = subset_prop_perRun$Astro
colData(spe_ruv)$Micro = subset_prop_perRun$Micro
colData(spe_ruv)$Oligo = subset_prop_perRun$Oligo

colData(spe_ruv)$Donors = gsub(' .*|\\/', '', colData(spe_ruv)$slide.name)
colData(spe_ruv)$Donors[-grep('PD', colData(spe_ruv)$Donors)] = paste0('C', colData(spe_ruv)$Donors[-grep('PD', colData(spe_ruv)$Donors)])

length(unique(colData(spe_ruv)$Donors))

#============================================ PCA

setwd(datadir)

library(factoextra)
library(tidyverse)

names(spe_ruv@assays@data@listData)

identical(colnames(spe_ruv@assays@data@listData$logcounts), row.names(colData(spe_ruv)))

res.pca <- prcomp(t(spe_ruv@assays@data@listData[["logcounts"]]), scale = TRUE)

p = fviz_pca_biplot(res.pca,
                    axes = c(1, 2),
                    # this stands for PC1 vs PC2
                    select.var = list(cos2 = 5),
                    # Individuals
                    geom.ind = "point",
                    labelsize = 6,
                    fill.ind = colData(spe_ruv)$segment, col.ind = "black",
                    pointshape = 21, pointsize = 4,
                    palette = "jco",
                    addEllipses = TRUE,
                    # Variables
                    alpha.var ="contrib", col.var = "contrib",
                    gradient.cols = "RdYlBu",
                    legend.title = list(fill = "Neurons", color = "Contrib",
                                        alpha = "Contrib")
) + ggtitle('TMM + RUV4') + coord_fixed() + theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"), legend.text=element_text(size=20), legend.title=element_text(size=20), plot.title = element_text(size = 30, face="bold"))

print(p)

getwd()

pdf(paste0("plotPairPCA_before_regressout.pdf"), width = 12, height = 14)
print(p)
dev.off()


#========= regress out

Expression = spe_ruv@assays@data@listData[["logcounts"]]
class(Expression)
dim(Expression)

CellTypeProportions = as.data.frame(res$prop_of_all[-which(row.names(res$prop_of_all) == 'Neurons'),])
dim(CellTypeProportions)

# replace na with 0
CellTypeProportions[1,][is.na(CellTypeProportions[1,])] = 0
CellTypeProportions[2,][is.na(CellTypeProportions[2,])] = 0
CellTypeProportions[3,][is.na(CellTypeProportions[3,])] = 0

NeuronalCellType = colData(spe_ruv)$segment
class(NeuronalCellType)

colnames(CellTypeProportions) = gsub('\\.','-',colnames(CellTypeProportions))
CellTypeProportions = CellTypeProportions[,match(colnames(Expression), colnames(CellTypeProportions))]
identical(colnames(Expression), colnames(CellTypeProportions))
dim(CellTypeProportions)

message("Regress out cell types")

length(NeuronalCellType)
dim(t(as.matrix(CellTypeProportions)))
dim((Expression))

# NeuronalCellType_M = NeuronalCellType
# NeuronalCellType_M[NeuronalCellType_M == "control"] = 1
# NeuronalCellType_M[NeuronalCellType_M == "PD-with-LB"] = 2
# NeuronalCellType_M[NeuronalCellType_M == "PD-without-LB"] = 3
# NeuronalCellType_M = as.numeric(NeuronalCellType_M)


#-- per gene I regress out the effect of Glial cells
j=1
for(j in 1:dim(Expression)[1])
{

  # row.names(Expression[j,])
  # colnames(Expression)[1]
  # as.numeric(Expression[j,1])
  # NeuronalCellType[1]
  # t(CellTypeProportions)[1,]
  
  
  pcFit <- lm(as.numeric(Expression[j,]) ~ NeuronalCellType + t(as.matrix(CellTypeProportions)) )
  
  summary(pcFit)
  
  test.cov = t(as.matrix(CellTypeProportions))
  remodelled <- as.numeric(Expression[j,] - rowSums(sapply(1:3, function(i) pcFit$coefficients[i+3]*test.cov[,i])))
  print(summary(remodelled))
  
  Expression[j,] = remodelled
  print(j)
}

dim(Expression)
dim(t(spe_ruv@assays@data@listData[["logcounts"]]))
    
res.pca <- prcomp(t(Expression),  scale = TRUE)

p = fviz_pca_biplot(res.pca,
                    select.var = list(cos2 = 5),
                    # Individuals
                    geom.ind = "point",
                    labelsize = 6,
                    fill.ind = colData(spe_ruv)$segment, col.ind = "black",
                    pointshape = 21, pointsize = 4,
                    palette = "jco",
                    addEllipses = TRUE,
                    # Variables
                    alpha.var ="contrib", col.var = "contrib",
                    gradient.cols = "RdYlBu",
                    legend.title = list(fill = "Neurons", color = "Contrib",
                                        alpha = "Contrib")
) + ggtitle('TMM + RUV4') + coord_fixed() + theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"), legend.text=element_text(size=20), legend.title=element_text(size=20), plot.title = element_text(size = 30, face="bold"))
print(p)

pdf(paste0("plotPairPCA_after_regressout.pdf"), width = 12, height = 14)
print(p)
dev.off()

library(dplyr)
library(ggforce)
count_normalized_mat = data.frame(GeneSymbol = row.names(Expression), Expression)
dim(count_normalized_mat)

library(data.table)
fwrite(count_normalized_mat, paste0(datadir, 'GeoMx_TMM_RUV4_CellTypeCorrection_readCount.txt'), quote = F, row.names = F, sep = '\t')


#================================================== further exploration of data before DEGA ==================================================

# There is a striking pattern that calls for explanation. This horseshoe or arch structure in the points is often an indicator of a sequential latent ordering or gradient in the data (Diaconis, Goel, and Holmes 2008). 
# Diaconis, Persi, Sharad Goel, and Susan Holmes. 2008. “Horseshoes in Multidimensional Scaling and Kernel Methods.” Annals of Applied Statistics 2: 777. 

## ----------------------------------------------------------------------------- ranked PCA

library("ade4")
library("factoextra")
library("sva")

identical(row.names(colData(spe_ruv)), colnames(assay(spe_ruv, 'counts')))

#--- ranked PCA for counts
Input = assay(spe_ruv, 'counts')
colnames(Input) = colData(spe_ruv)$Donors

# Instead of using the continuous, somehow normalized data, we use a robust analysis replacing the values by their ranks. The lower values are considered ties encoded as a threshold chosen to reflect the number of expected taxa thought to be present:
rankthreshPCA = function(x, threshold = 3000) {
  ranksM = apply(x, 2, rank)
  ranksM[ranksM < threshold] = threshold
  ranksM = threshold - ranksM
  dudi.pca(t(ranksM), scannf = FALSE, nf = 2)
}

pcaGeoMx = rankthreshPCA(Input)

pdf(paste0("RankedPCAPlot_batch_counts.pdf"), width = 8, height = 6)

fviz(pcaGeoMx, element = "ind", axes = c(1, 2), geom = c("point", "text"),
     habillage = colData(spe_ruv)$batch, repel = TRUE, palette = "Dark2",
     addEllipses = TRUE, ellipse.type = "convex", label = "none") + 
  ggtitle("") + coord_fixed()

dev.off()

pdf(paste0("RankedPCAPlot_state_counts.pdf"), width = 8, height = 6)

fviz(pcaGeoMx, element = "ind", axes = c(1, 2), geom = c("point", "text"),
     habillage = colData(spe_ruv)$state, repel = TRUE, palette = "Dark2",
     addEllipses = TRUE, ellipse.type = "convex", label = "none") + 
  ggtitle("") + coord_fixed()

dev.off()


#--- ranked PCA for logcounts
identical(row.names(colData(spe_ruv)), colnames(assay(spe_ruv, 'logcounts')))

Input = assay(spe_ruv, 'logcounts')
colnames(Input) = colData(spe_ruv)$Donors

pcaGeoMx = rankthreshPCA(Input)

pdf(paste0("RankedPCAPlot_batch_logcounts.pdf"), width = 8, height = 6)

fviz(pcaGeoMx, element = "ind", axes = c(1, 2), geom = c("point", "text"),
     habillage = colData(spe_ruv)$batch, repel = TRUE, palette = "Dark2",
     addEllipses = TRUE, ellipse.type = "convex", label = "none") + 
  ggtitle("") + coord_fixed()

dev.off()

pdf(paste0("RankedPCAPlot_state_logcounts.pdf"), width = 8, height = 6)

fviz(pcaGeoMx, element = "ind", axes = c(1, 2), geom = c("point", "text"),
     habillage = colData(spe_ruv)$state, repel = TRUE, palette = "Dark2",
     addEllipses = TRUE, ellipse.type = "convex", label = "none") + 
  ggtitle("") + coord_fixed()

dev.off()

#--- ranked PCA for log-counts + regress-out glial cells
Input = Expression
identical(row.names(colData(spe_ruv)), colnames(Input))
colnames(Input) = colData(spe_ruv)$Donors

pcaGeoMx = rankthreshPCA(Input)

pdf(paste0("RankedPCAPlot_batch_logcounts_regress-out.pdf"), width = 8, height = 6)

fviz(pcaGeoMx, element = "ind", axes = c(1, 2), geom = c("point", "text"),
     habillage = colData(spe_ruv)$batch, repel = TRUE, palette = "Dark2",
     addEllipses = TRUE, ellipse.type = "convex", label = "none") + 
  ggtitle("") + coord_fixed()

dev.off()

pdf(paste0("RankedPCAPlot_state_logcounts_regress-out.pdf"), width = 8, height = 6)

fviz(pcaGeoMx, element = "ind", axes = c(1, 2), geom = c("point", "text"),
     habillage = colData(spe_ruv)$state, repel = TRUE, palette = "Dark2",
     addEllipses = TRUE, ellipse.type = "convex", label = "none") + 
  ggtitle("") + coord_fixed()

dev.off()

## --------- batch correction using ComBat; This is not effective  
# model0 = model.matrix(~1, colData(spe_ruv)$batch)
# 
# combatIBD = ComBat(dat = Input, batch = colData(spe_ruv)$batch, mod = model0)
# 
# pcaDayBatRM = rankthreshPCA(combatIBD)
# 
# fviz(pcaDayBatRM, element = "ind", geom = c("point", "text"),
#      habillage = colData(spe_ruv)$batch, repel=TRUE, palette = "Dark2", addEllipses = TRUE,
#      ellipse.type = "convex", axes =c(1,2), label = "none") + coord_fixed() + ggtitle("")
# 
# fviz(pcaDayBatRM, element = "ind", axes = c(1, 2), geom = c("point", "text"),
#      habillage = colData(spe_ruv)$state, repel = TRUE, palette = "Dark2",
#      addEllipses = TRUE, ellipse.type = "convex", label = "none") + 
#   ggtitle("") + coord_fixed()
# 
# 
# fviz_eig(pcaDayBatRM, bar_width = 0.6) + ggtitle("")
# 
# # While it's theoretically possible to apply RUV4 and ComBat sequentially for batch correction, it's generally not recommended and could lead to unintended consequences. 
# # These methods address unwanted variation in different ways, and applying them one after the other might over-correct the data or remove biologically relevant signals.


#---- tSNE and MDS plots

library(data.table)
# datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/")
count_mat3 = Expression # fread(paste0(datadir, 'GeoMx_TMM_RUV4_CellTypeCorrection_readCount.txt'), stringsAsFactors = F, header = T)
count_mat3 = as.data.frame(count_mat3)
# row.names(count_mat3) = make.names(count_mat3$GeneSymbol,unique=T) 
# count_mat3 = count_mat3[,-1]
# colnames(count_mat3) = gsub('\\.','-',gsub('.dcc','',colnames(count_mat3)))
dim(count_mat3)

MD = fread(paste0(datadir, 'MetaData.txt'), stringsAsFactors = F, header = T)
MD = as.data.frame(MD)
colnames(MD)

MD = MD[which(MD$SampleID %in% colnames(count_mat3)),]
MD = MD[match(colnames(count_mat3), MD$SampleID),]

identical(colnames(count_mat3), MD$SampleID)

library(SummarizedExperiment)

cellt = MD$state

palette.colors(palette = "Okabe-Ito")
# [1] "#000000" "#E69F00" "#56B4E9" "#009E73"
# [5] "#F0E442" "#0072B2" "#D55E00" "#CC79A7"
# [9] "#999999"

colsn = c("#0072B2", "#D55E00", "#009E73")

dim(count_mat3)
length(cellt)

# compute the distances between the rows of a data matrix.
# distance "manhattan"

dist1n.l1 = dist(t(count_mat3), "manhattan")
dist2n.euclid = dist(t(count_mat3))

# The Manhattan distance, also known as the L1 distance
ce1Mds = cmdscale(dist1n.l1, k = 10, eig = TRUE)

# L2 distance is also known as the Euclidean distance. 
ce2Mds = cmdscale(dist2n.euclid, k = 10, eig = TRUE)

# ----------------------------------------- Toy example of euclidean distance
# Sample Gene Expression Matrix (genes as rows, samples as columns)
expression_matrix <- matrix(
  c(10, 12, 5, 8,   # Gene 1 expression across 4 samples
    15, 11, 9, 13,  # Gene 2 expression across 4 samples
    6,  8, 12, 7),  # Gene 3 expression across 4 samples
  nrow = 3,
  byrow = TRUE
)
colnames(expression_matrix) <- paste0("Sample-", 1:4)
rownames(expression_matrix) <- paste0("Gene-", 1:3)

print("Gene Expression Matrix:")
kable(expression_matrix)

#-------------------------------------------------------------------------------
# 1. Euclidean Distance
#-------------------------------------------------------------------------------

# We need to transpose the matrix so that samples are rows for dist()
transposed_matrix <- t(expression_matrix)

# Calculate Euclidean distance
kable(as.matrix(dist(transposed_matrix, method = "euclidean", diag = TRUE)))
# euclidean_distances <- dist(transposed_matrix, method = "euclidean", diag = TRUE)

# Formula Explanation for Euclidean Distance between Sample_1 and Sample_2:
# Sample_1 = (10, 15, 6)
# Sample_2 = (12, 11, 8)
# d_E(Sample_1, Sample_2) = sqrt((10-12)^2 + (15-11)^2 + (6-8)^2)
#                      = sqrt((-2)^2 + (4)^2 + (-2)^2)
#                      = sqrt(4 + 16 + 4)
#                      = sqrt(24)
#                      ≈ 4.899

# ----------------------------------------- End of the toy example

# L1 distance is also known as the manhattan distance.  
perc1  = round(100*sum(ce1Mds$eig[1:2])/sum(ce1Mds$eig))
perc2  = round(100*sum(ce2Mds$eig[1:2])/sum(ce2Mds$eig))

library("dplyr")
library("ggplot2")
plotscree = function(x, m = length(x$eig)) {
  ggplot(tibble(eig = x$eig[seq_len(m)], k = seq(along = eig)),
         aes(x = k, y = eig)) + theme_minimal() +
    scale_x_discrete("k", limits = as.factor(seq_len(m))) + 
    geom_bar(stat = "identity", width = 0.5, fill = "#ffd700", col = "#0057b7")
}

plotscree(ce1Mds, m = 4)
plotscree(ce2Mds, m = 4)

# The Manhattan distance, also known as the L1 distance
c1mds = ce1Mds$points[, 1:2] |>
  `colnames<-`(paste0("L1_PCo", 1:2)) |>
  as_tibble()

pdf(paste0("PCo_manhattan_logcounts_regress-out.pdf"), width = 8, height = 6)

ggplot(c1mds, aes(x = L1_PCo1, y = L1_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
  scale_colour_manual(values = colsn) + guides(color = "none") + theme_linedraw() + geom_point(size = 4)

dev.off()

# L2 distance is also known as the Euclidean distance. 
c2mds = ce2Mds$points[, 1:2] |>
  `colnames<-`(paste0("L2_PCo", 1:2)) |>
  as_tibble()

c2mds = as.data.frame(c2mds)
row.names(c2mds) = colnames(count_mat3)
class(c2mds)

pdf(paste0("PCoPlot_euclidean_logcounts_regress-out.pdf"), width = 8, height = 6)

ggplot(c2mds, aes(x = L2_PCo1, y = L2_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
  scale_colour_manual(values = colsn) + guides(color = "none") + theme_linedraw() + geom_point(size = 4)

dev.off()

# L2 distance is also known as the Euclidean distance. [without DSP-1001660034021-C-A10]

pdf(paste0("PCoPlot_euclidean_logcounts_regress-out_withoutDSP-1001660034021-C-A10.pdf"), width = 8, height = 6)

ggplot(c2mds[-which(row.names(c2mds) == 'DSP-1001660034021-C-A10'),], aes(x = L2_PCo1, y = L2_PCo2, color = cellt[-which(row.names(c2mds) == 'DSP-1001660034021-C-A10')])) +
  geom_point(aes(color = cellt[-which(row.names(c2mds) == 'DSP-1001660034021-C-A10')]), alpha = 0.6) +
  scale_colour_manual(values = colsn) + guides(color = "none") + theme_linedraw() + geom_point(size = 4)

dev.off()



# legend
ggpcolor = ggplot(c1mds,aes(x=L1_PCo1,y=L1_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
  scale_colour_manual(values=colsn, name = "cell type")
g_legend = function(a) {
  gt = ggplot_gtable(ggplot_build(a))
  leg = which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  gt$grobs[[leg]]
}

library(grid)
grid.draw(g_legend(ggpcolor))

# "t-SNE"
library("Rtsne")
restsne = Rtsne(t(count_mat3), dims = 2, perplexity = 30, verbose = FALSE,
                max_iter = 900)

dftsne = restsne$Y[, 1:2] |>
  `colnames<-`(paste0("tSNE", 1:2)) |>
  as_tibble()

pdf(paste0("tsnePlot_logcounts_regress-out.pdf"), width = 8, height = 6)

ggplot(dftsne,aes(x = tSNE1, y = tSNE2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) + geom_point(size = 4) +
  scale_color_manual(values = colsn) + guides(color = "none") + theme_linedraw()

dev.off()

#-- explore the outlier in the tSNE plot
#View(cbind(MD[, c('SampleID', 'state')], dftsne))

QC_PASS_inc_GeneDetection_830 = fread('QC_PASS_inc_GeneDetection_830.txt', stringsAsFactors = F, header = T)
QC_PASS_inc_GeneDetection_830 = as.data.frame(QC_PASS_inc_GeneDetection_830)
# View(QC_PASS_inc_GeneDetection_830[which(QC_PASS_inc_GeneDetection_830$SampleID== 'DSP-1001660034021-C-A10'),])

# DSP-1001660034021-C-A10 - min estimated nuclei 20
# PD-without-LB
# 5.6253434
# -2.2945440
# DSP-1001660034021-C-A10 - min estimated nuclei 5
# -10.0889329
# 11.72131035


#---- visualization of expression level per gene using revised Expression profile
meta.data = data.frame(slide.name = colData(spe_ruv)$slide.name, state = colData(spe_ruv)$state)
row.names(meta.data) = row.names(colData(spe_ruv))

meta.data$slide.name[which(meta.data$state == 'PD-with-LB')] = paste0('V', meta.data$slide.name[which(meta.data$state == 'PD-with-LB')])
meta.data$slide.name[which(meta.data$state == 'PD-without-LB')] = paste0('R', meta.data$slide.name[which(meta.data$state == 'PD-without-LB')])

identical(colnames(Expression), row.names(meta.data))

library(Seurat)
target_Data_Seurat <- CreateSeuratObject(counts = Expression, project = "Seurat",  meta.data = meta.data, assay = "GeoMx")

library(scater)
sce <- as.SingleCellExperiment(target_Data_Seurat)
# assay(sce) # uses the normalized values

query = 'KCNJ6'
query = 'CALB1'
query = "SNCA"
# ENO1

library(dittoSeq)
# dittoPlot(sce, query, group.by = "state")

dittoPlot(sce, query, group.by = "state",
          plots = c("vlnplot", "jitter", "boxplot"),
          # change the color and size of jitter points
          jitter.color = "blue", jitter.size = 0.7,
          # change the outline color and width, and remove the fill of boxplots
          boxplot.color = "white", boxplot.width = 0.1,
          boxplot.fill = FALSE,
          # change how the violinplot widths are normalized across groups
          vlnplot.scaling = "count"
)

#-----------------
# print(fviz_eig(res.pca, addlabels = TRUE))
# Graph of the variables
# fviz_pca_var(res.pca, col.var = "black")

# fviz_cos2(res.pca, choice = "var", axes = 1:2)
# fviz_pca_var(res.pca, col.var = "cos2",
#              gradient.cols = c("black", "orange", "green"),
#              repel = TRUE)

#--- variable selection using PCA before and after cell type correction
name = 'logcounts'

if(name == 'logcounts')
{
  #before
  before.pca <- prcomp(t(spe_ruv@assays@data@listData[["logcounts"]]),  scale = TRUE)
  
  #after
  after.pca <- prcomp(t(Expression),  scale = TRUE)
  
  # attr(spe_ruv@int_colData@listData$reducedDims$PCA, "varExplained")
  # attr(spe_ruv@int_colData@listData$reducedDims$PCA, "rotation")
  
  #print(name)
  
  # Loadings before
  loadings <- before.pca$rotation
  
  number_of_PC_to_keep = 2;
  pc_importance_short <- abs(loadings[, 1:number_of_PC_to_keep])
  variable_importance_short <- rowSums(pc_importance_short)
  selected_variables_short <- names(sort(variable_importance_short, decreasing = TRUE))
  length(selected_variables_short)
  
  #Example using percentile.
  threshold_percentile = 0.95;
  threshold_value = quantile(variable_importance_short, threshold_percentile);
  selected_variables_percentile_before = variable_importance_short[which(variable_importance_short > threshold_value)]
  length(selected_variables_percentile_before)
  
  
  # Loadings after
  loadings <- after.pca$rotation
  
  number_of_PC_to_keep = 2;
  pc_importance_short <- abs(loadings[, 1:number_of_PC_to_keep])
  variable_importance_short <- rowSums(pc_importance_short)
  selected_variables_short <- names(sort(variable_importance_short, decreasing = TRUE))
  length(selected_variables_short)
  
  #Example using percentile.
  threshold_percentile = 0.95;
  threshold_value = quantile(variable_importance_short, threshold_percentile);
  selected_variables_percentile_after = variable_importance_short[which(variable_importance_short > threshold_value)]
  length(selected_variables_percentile_after)
  
  length(unique(names(selected_variables_percentile_before)))
  length(intersect(names(selected_variables_percentile_before), names(selected_variables_percentile_after)))
  setdiff(names(selected_variables_percentile_before), names(selected_variables_percentile_after))
  
  x <- list()
  
  x[['All Cell Types']] = names(selected_variables_percentile_before)
  x[['Neurons']] = names(selected_variables_percentile_after)
  
  pdf(paste0("PCA-associated-genes_regress_out.pdf"), width = 10, height = 7)
  
  names(x)
  
  library(ggvenn)
  print(ggvenn(
    x, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 8, text_size = 8
  ))
  
  dev.off()
  
  
  getwd()
  write.table(names(selected_variables_percentile_after), paste0("After_regressout_cellTypes_selected_variables_95percentile.txt"), sep = '\t', col.names = F)
}

#================================================== co-expression analysis

colData(spe_ruv)$segment
identical(colnames(Expression), row.names(colData(spe_ruv)))

sample_annot = data.frame(SampleName = row.names(colData(spe_ruv)), Class = colData(spe_ruv)$segment)

#- NOTE ---- remove outlier samples
Expression = Expression[, -which(colnames(Expression) == 'DSP-1001660034021-C-A10')]
sample_annot = sample_annot[-which(sample_annot$SampleName == 'DSP-1001660034021-C-A10'),]
identical(colnames(Expression), sample_annot$SampleName)

spe_ruv = spe_ruv[,-which(colnames(assay(spe_ruv, 'logcounts')) == 'DSP-1001660034021-C-A10')]
assay(spe_ruv, 'logcounts') = Expression
identical(colnames(assay(spe_ruv, 'logcounts')), colnames(Expression))

#-----

#sample_annot$Class[-which(sample_annot$Class == 'control')] = 'PD'

#BiocManager::install("CEMiTool")
library("CEMiTool")
# run cemitool with sample annotation
cem <- cemitool(Expression, sample_annot, apply_vst = FALSE, cor_method = "spearman") # cor_method = "spearman", apply_vst Logical. If TRUE, will apply Variance Stabilizing Transform

sample_annotation(cem,
                  sample_name_column="SampleName",
                  class_column="Class") <- sample_annot

# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")

# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[2]

# read GMT file
gmt_fname <- '/Users/isarnassiri/Documents/GeoMx_Analysis/reference_datasets/msigdb_v2024.1.Hs_GMTs/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt'
# gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
length(unique(gmt_in$term))
length(unique(gmt_in$gene))

# perform over representation analysis
cem <- mod_ora(cem, gmt_in)

# plot ora results
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[1]

# read interactions
library(data.table)
int_df = fread('/Users/isarnassiri/Documents/GeoMx_Analysis/reference_datasets/BIOGRID-PROJECT-alzheimers_disease_project-4.4.245/BIOGRID-PROJECT-alzheimers_disease_project-INTERACTIONS-4.4.245.tab3.txt', stringsAsFactors = F, header = T)
int_df = as.data.frame(int_df)
int_df = int_df[which(int_df$`Organism Name Interactor A` == "Homo sapiens"),]
table(int_df$`Organism Name Interactor B`)

int_df = data.frame(gene1symbol = int_df$`Official Symbol Interactor A`, gene2symbol = int_df$`Official Symbol Interactor B`)

# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]

# save all plots
save_plots(cem, "all", force=TRUE, directory="/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/cemitool/Plots")

# create report as html document
generate_report(cem, directory="/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/cemitool/Report" , force=TRUE  )

# write analysis results into files
write_files(cem, directory="/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/cemitool/Tables", force=TRUE  )

table(cem@module$modules)

M1G = cem@module$genes[cem@module$modules == "M1"]
write.table(M1G, "/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/GeoMx_M1.txt", quote = F, row.names = F, sep = '\t', col.names = F)

M2G = cem@module$genes[cem@module$modules == "M2"]
write.table(M2G, "/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/GeoMx_M2.txt", quote = F, row.names = F, sep = '\t', col.names = F)

#================================================== compare co-expression modules and gene sets
#------ selected brain regions

files = list.files('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/USED/', pattern = 'Region_')
datadir2 = '/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/USED/'
setwd(datadir2)

require(data.table) ## 1.9.2 or 1.9.3
genesets_brain_region = rbindlist(lapply(files, fread), idcol="ID")

genesets_brain_region$ID = gsub('.*rna_|.tsv', '', files)[genesets_brain_region$ID]
table(genesets_brain_region$ID)
table(duplicated(genesets_brain_region$Gene))

toMatch <- c("Cancer", "Mouse", "mouse", "pig", "blood", "Blood", "Evidence", "HPA evidence", "Cancer", "cancer", "cell line", "RNA single cell", "RNA tissue", "Tissue expression cluster", "Cell line expression cluster", "Single cell expression cluster", "UniProt evidence", "NeXtProt evidence", "CCD Protein", "CCD Transcript" , "Interactions"  )
matches <- unique(grep(paste(toMatch,collapse="|"), colnames(genesets_brain_region), value=TRUE))

genesets_brain_region = as.data.frame(genesets_brain_region)
class(genesets_brain_region)
genesets_brain_region = genesets_brain_region[,-which(colnames(genesets_brain_region) %in% matches)]
colnames(genesets_brain_region)

table((genesets_brain_region$`RNA single nuclei brain distribution`))
genesets_brain_region = genesets_brain_region[-which(genesets_brain_region$`RNA single nuclei brain distribution` %in% c('', 'Not detected')),]
dim(genesets_brain_region)

table(genesets_brain_region$ID)

#-------------------- comparison with GeoMx co-expression module

library(data.table)
M1G = fread("/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/GeoMx_M1.txt", stringsAsFactors = F, header = F)

files = list.files('/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis', pattern = 'sc_brain_region_')
setwd('/Users/isarnassiri/Documents/GeoMx_Analysis/Analysis')
genesets = rbindlist(lapply(files, fread), idcol="ID")
genesets$ID = gsub('.*rna_|.tsv', '', files)[genesets$ID]
dim(genesets)

library(knitr)
kable(table(genesets$ID))
kable(table(genesets_brain_region$ID))

library(dplyr)
# I have some duplicate genes per group
# View(genesets_brain_region %>% group_by(ID, Gene) %>% filter(n()>1) )

x <- list(
  `Amygdala` = unique(genesets_brain_region$Gene[genesets_brain_region$ID == 'Region_amygdala']),
  `Hippocampus` = unique(genesets_brain_region$Gene[genesets_brain_region$ID == 'Region_hippocampal_formation']),
  `Midbrain` = unique(genesets_brain_region$Gene[genesets_brain_region$ID == 'Region_midbrain'])
)

intersection_selected_brain_regions <- Reduce(intersect, x)
length(intersection_selected_brain_regions)

celltypes = c("Neurons", "Microglia", "Astrocyte", "Oligodendrocyte")

totalGenesCellType = 0
for(CellTypeName in celltypes)
{
  eval(call("<-", as.name(paste0(CellTypeName, '_input')), genesets$Gene[which(genesets$ID == CellTypeName & genesets$Gene %in% intersection_selected_brain_regions)]))
  totalGenesCellType = totalGenesCellType + length(genesets$Gene[which(genesets$ID == CellTypeName & genesets$Gene %in% intersection_selected_brain_regions)])
}
totalGenesCellType

GWAS_input = fread('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/USED/GWAS/GWAS_hits_L32.txt', stringsAsFactors = F, header = T)
GWAS_input = GWAS_input$SYMBOL[which(GWAS_input$SYMBOL %in% intersection_selected_brain_regions)]

PD_input = fread('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/USED/PD.txt', stringsAsFactors = F, header = F)
PD_input = PD_input$V1[which(PD_input$V1 %in% intersection_selected_brain_regions)]

Drug_input = fread('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/USED/Drugable.txt', stringsAsFactors = F, header = T)
Drug_input = Drug_input$`Gene Name`[which(Drug_input$`Gene Name` %in% intersection_selected_brain_regions)]

ND_input = fread('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/USED/ND.txt', stringsAsFactors = F, header = F)
ND_input = ND_input$V1[which(ND_input$V1 %in% intersection_selected_brain_regions)]
class(ND_input)

# Hallmark gene sets
gmt_fname <- '/Users/isarnassiri/Documents/GeoMx_Analysis/reference_datasets/msigdb_v2024.1.Hs_GMTs/h.all.v2024.1.Hs.symbols.gmt'
gmt_in <- read_gmt(gmt_fname)
length(unique(gmt_in$term))
length(unique(gmt_in$gene))

Hypoxia_input = gmt_in[which(gmt_in$term %in% "HALLMARK_HYPOXIA"),]
Hypoxia_input = Hypoxia_input$gene[which(Hypoxia_input$gene %in% intersection_selected_brain_regions)]
length(Hypoxia_input)

Inflammation_input = gmt_in[which(gmt_in$term %in% "HALLMARK_INFLAMMATORY_RESPONSE"),]
Inflammation_input = Inflammation_input$gene[which(Inflammation_input$gene %in% intersection_selected_brain_regions)]
length(Inflammation_input)

GeoMx_input = fread('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/GeoMx_M1.txt', stringsAsFactors = F, header = F)
GeoMx_input = GeoMx_input$V1

ND_input = unique(ND_input, PD_input)
length(ND_input)

all_objects <- ls()
selected_objects = all_objects[grep('_input', all_objects)]
selected_objects = selected_objects[-grep("summary_input", selected_objects)]

setdiff(GeoMx_input, intersection_selected_brain_regions)

for(CellTypeName in celltypes)
{
  eval(call("<-", as.name(paste0(CellTypeName, '_input')), genesets$Gene[which(genesets$ID == CellTypeName & genesets$Gene %in% intersection_selected_brain_regions)]))
  totalGenesCellType = totalGenesCellType + length(genesets$Gene[which(genesets$ID == CellTypeName & genesets$Gene %in% intersection_selected_brain_regions)])
}

number_of_rows = length(unique(as.character(unlist(mget(selected_objects)))))
# number_of_rows = length(cem@module$genes[cem@module$modules == "M1"])
number_of_columns = length(selected_objects)

upsetPlotInput = matrix(0, nrow = number_of_rows, ncol = number_of_columns)
colnames(upsetPlotInput) = selected_objects # gsub("_input","",selected_objects)
row.names(upsetPlotInput) = unique(as.character(unlist(mget(selected_objects)))) # cem@module$genes[cem@module$modules == "M1"] #

for (n in selected_objects) {
  print(n)
  print(length(upsetPlotInput[which(row.names(upsetPlotInput) %in% get(n)), which(colnames(upsetPlotInput) == n)]))

  upsetPlotInput[which(row.names(upsetPlotInput) %in% get(n)),which(colnames(upsetPlotInput) == n)] = 1

  # if(length(upsetPlotInput[which(row.ns(upsetPlotInput) %in% get(n)),which(colns(upsetPlotInput) == n)]) > 0)
  # {
  #   upsetPlotInput[which(row.ns(upsetPlotInput) %in% get(n)),which(colns(upsetPlotInput) == n)] = 1
  # }else{
  #   upsetPlotInput = upsetPlotInput[,-which(colns(upsetPlotInput) == n)]
  # }
}

class(upsetPlotInput)
colnames(upsetPlotInput) = gsub("_input","", colnames(upsetPlotInput))

library(UpSetR)
my_upset_plot <- upset(as.data.frame(upsetPlotInput),
                       nintersects = NA,
                       nsets = 20,
                       order.by = "freq",
                       decreasing = T,
                       mb.ratio = c(0.6, 0.4),
                       number.angles = 0,
                       text.scale = 1.1,
                       point.size = 2.8,
                       line.size = 1,
                       set_size.show = TRUE
)
print(my_upset_plot)
# Save as PDF (vector graphics, good for print)
pdf("/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/upset_plot.pdf", width = 10, height = 6) # Adjust width and height as needed
print(my_upset_plot) # You need to explicitly print the plot object
dev.off()

for(n in selected_objects)
{
  print(length(mget(n)))
}

merged_vector <- unlist(mget(selected_objects[which(selected_objects %in% c("GWAS_input", "Astrocyte_input", "Microglia_input", "Neurons_input", "Oligodendrocyte_input"))]))
length(merged_vector)

length(unique(merged_vector))
length(unique(c(intersect(ND_input, GeoMx_input))))
length(unique(c(intersect(Hypoxia_input, GeoMx_input))))
length(unique(c(intersect(Inflammation_input, GeoMx_input))))
length(unique(c(intersect(Drug_input, GeoMx_input))))

Selected_Genes = unique(c(unique(merged_vector),
unique(c(intersect(ND_input, GeoMx_input))),
unique(c(intersect(Hypoxia_input, GeoMx_input))),
unique(c(intersect(Inflammation_input, GeoMx_input))),
unique(c(intersect(ND_input, Drug_input))),
unique(c(intersect(Drug_input, GeoMx_input)))))

Selected_Genes = Selected_Genes[Selected_Genes != "ENSG00000265118"]

#-- selected genes for Xenium panel design
getwd()
write.table(Selected_Genes, '/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/USED/Selected_Genes.txt', quote = F, row.names = F, sep = '\t', col.names = F)

#================================================== end: further exploration of data before DEGA ==================================================

setwd(datadir)

# in the case that you want to see the impact of normalization method on number of DEGs
names = data.frame(names1 = names(spe_ruv@assays),
                   names2 = c('Counts','TMM + RUV4','Limma TMM','Cyclic Loess','Quantile','Library size'))
name = "logcounts"

#------------------ select a reference ------------------

LMM = FALSE # for with Lewy body vs. without Lewy body
EFA_usage = FALSE

# LB509TH - TH (TH is reference)
# contr.matrix <- makeContrasts(BvT = PDwLBP - PDwoLBP, levels = colnames(design))

states = c("PDwLBP - control", "PDwoLBP - control", "PDwLBP - PDwoLBP")
state = "PDwLBP - control"
 
for(state in states)
{

#============================================ 
library(stringr)
state1 = str_split_fixed(state, " - ", 2)[1]
state2 = str_split_fixed(state, " - ", 2)[2]

if(state1 == 'PDwLBP'){state1 = "PD-with-LB"}
if(state1 == 'PDwoLBP'){state1 = "PD-without-LB"}
if(state2 == 'PDwLBP'){state2 = "PD-with-LB"}
if(state2 == 'PDwoLBP'){state2 = "PD-without-LB"}

spe_ruv3 <- spe_ruv[, grep(Run, colData(spe_ruv)$scan.name)]
spe_ruv3 <- spe_ruv3[, which(colData(spe_ruv3)$segment %in% c(state1, state2))]

#--- test possibilities
# colnames(as.data.frame(colData(spe_ruv3)))
# spe_ruv3 <- spe_ruv3[, which(colData(spe_ruv3)$NucleiCount >= 10)]
# as.data.frame(assay(spe_ruv3, name))[,which(colnames(as.data.frame(assay(spe_ruv3, name))) %in% c('DSP-1001660036859-D-B06', 'DSP-1001660034021-C-B12'))]
#---

colData(spe_ruv3)$slide.name[grep('PD ', colData(spe_ruv3)$slide.name)] = gsub('PD ', 'PD', colData(spe_ruv3)$slide.name[grep('PD ', colData(spe_ruv3)$slide.name)])
colData(spe_ruv3)$slide.name[-grep('PD', colData(spe_ruv3)$slide.name)] = paste0('C', colData(spe_ruv3)$slide.name[-grep('PD', colData(spe_ruv3)$slide.name)])

colData(spe_ruv3)$slide.name[which(colData(spe_ruv3)$state == 'PD-with-LB')] = paste0('V', colData(spe_ruv3)$slide.name[which(colData(spe_ruv3)$state == 'PD-with-LB')])
colData(spe_ruv3)$slide.name[which(colData(spe_ruv3)$state == 'PD-without-LB')] = paste0('R', colData(spe_ruv3)$slide.name[which(colData(spe_ruv3)$state == 'PD-without-LB')])

colData(spe_ruv3)$slide.name = gsub(' .*|_.*','', colData(spe_ruv3)$slide.name)
colData(spe_ruv3)$slide.name = gsub('\\/', '.', colData(spe_ruv3)$slide.name)

table(colData(spe_ruv3)$state)
# PD-with-LB PD-without-LB 
# 54            95 

# if(EFA_usage == TRUE & state == "PDwLBP - PDwoLBP")
# {
#   tier_PDwithoutLBcompacta = fread("PD-without-LBcompacta.txt", header = T, stringsAsFactors = F)
#   tier_PDwithoutLBcompacta = as.data.frame(tier_PDwithoutLBcompacta)
# 
#   tier_PDwithLBcompacta = fread("PD-with-LBcompacta.txt", header = T, stringsAsFactors = F)
#   tier_PDwithLBcompacta = as.data.frame(tier_PDwithLBcompacta)
# 
#   #-------
#   spe_ruv3 = spe_ruv3[, which(row.names(colData(spe_ruv3)) %in% gsub('\\.','-', c(tier_PDwithLBcompacta$x, tier_PDwithoutLBcompacta$x)))]
# }

#============================================ keep pair samples for PDwLBP vs PDwoLBP

# if(state == "PDwLBP - PDwoLBP")
# {
# 
#   #--- more than one segment per donor
#   count_segments_per_donor = as.data.frame(table(colnames(input_heatmap)))
#   count_segments_per_donor_selected = count_segments_per_donor[which(count_segments_per_donor$Freq > 1),]
#   
#   CountMatrix = CountMatrix[,which(colnames(CountMatrix) %in% as.character(count_segments_per_donor_selected$Var1))]
#   #---
#   
#   toMatch <- gsub('^.', '', colnames(CountMatrix))[duplicated(gsub('^.', '', colnames(CountMatrix) )) ]
#   
#   # OR
#   matches <- unique(grep(paste(toMatch, collapse="|"), colnames(CountMatrix), value=TRUE))
#   
#   CountMatrix_subset = CountMatrix[,which(colnames(CountMatrix) %in% matches)]
#   dim(CountMatrix_subset)
# 
#   #-------
#   spe_ruv3 = spe_ruv3[, which(colData(spe_ruv3)$slide.name %in% matches)]
# }

#============================================ pseudobulk
input_heatmap = assay(spe_ruv3, 'logcounts')
CN = colData(spe_ruv3)$scan.name

identical(colnames(input_heatmap), row.names(colData(spe_ruv3)))
colnames(input_heatmap) = colData(spe_ruv3)$slide.name
table(colData(spe_ruv3)$slide.name)

patterns <- names(table(colnames(input_heatmap)))
CountMatrix <- sapply(patterns, function(xx) rowSums(input_heatmap[,grep(xx, names(input_heatmap)), drop=FALSE]))  # loop through
dim(CountMatrix)

#============================================  

#--- summary of segments used for DEGs analysis
library(edgeR)
library(limma)
# batch is a blocking factor and using in as a covarinace make the model two-factor analysis (lower level of noise and more parameters that beed to be estimated (the fit has fewer degreees of freedom))
# set up a paired analysis by adding donors as a blocking factor.

rm(design)
# design <- model.matrix(~0 + segment + Donors + Astro + Micro + Oligo, data = colData(spe_ruv3)) # Adding ruv_W1 + ruv_W2 reduces the number of DEGs by up to 3.
# View( as.data.frame(colData(spe_ruv3)))

# I tried batch_combined and batchslide as covariance. batchslide generates more DEGs and I donot get PD as first enrich disease. I decided to proceed with batch_combined.

dim(colData(spe_ruv3))
colnames(colData(spe_ruv3))
design <- model.matrix(~0 + segment + batch_combined, data = colData(spe_ruv3)) # Adding ruv_W1 + ruv_W2 reduces the number of DEGs by up to 3.
colnames(design)
dim(design)

# Use interaction terms: If you suspect that the biological effect of interest might differ across batches, consider including interaction terms between the batch and the condition in your model.
# Account for Batch Effects in Your Statistical Model in Addition to Correction: This is often recommended over aggressive batch removal.
# + batch:state - does not work in makeContrasts Donors + 

# if(state == "PDwLBP - PDwoLBP")
# {
#   rm(design)
#   # design <- model.matrix(~0 + segment + Donors + Astro + Micro + Oligo, data = colData(spe_ruv3)) # Adding ruv_W1 + ruv_W2 reduces the number of DEGs by up to 3.
#   # View( as.data.frame(colData(spe_ruv3)))
#   design <- model.matrix(~0 + segment + ruv_W1 + ruv_W2 + Donors + Astro + Micro + Oligo, data = colData(spe_ruv3)) # Adding ruv_W1 + ruv_W2 reduces the number of DEGs by up to 3.
#   colnames(design)
# }

colnames(design) <- gsub("^segment","", colnames(design))
colnames(design) <- gsub("\\+","",colnames(design)) # Double check this part
colnames(design)[which(colnames(design) == "PD-with-LB")] = "PDwLBP"
colnames(design)[which(colnames(design) == "PD-without-LB")] = "PDwoLBP"

contr.matrix <- makeContrasts(BvT = state, levels = colnames(design))

query = gsub(' - .*','',state)
control = gsub('.* - ','',state)

dim(design)
dim(colData(spe_ruv3))

#============================================ 
dim(design)
dim(assay(spe_ruv3, name))  
  
print(name)

if(state == "PDwLBP - PDwoLBP")
{
  # Skew towards 1 (or a spike at 1): This can happen if:
  # 1. The statistical test is not appropriate for your data (e.g., assumptions are violated).
  fit <- lmFit(assay(spe_ruv3, name), design = design, method = 'robust', maxiter = 25, psi = 'psi.huber') # 
  
  # 2. Genes with very low or zero counts are being inappropriately analyzed, leading to p-values of 1.
  # If I filter genes I get a perfect p-valule plot; but I leave it as it is because the number of up-reg is too low.
  # spe_ruv4 <- addPerROIQC(spe_ruv3, rm_genes = TRUE, min_count = 3, sample_fraction = 0.95)
  # dim(spe_ruv4)
  # fit <- lmFit(assay(spe_ruv4, name), design = design, method = 'robust', maxiter = 25, psi = 'psi.huber')

}else{
  fit <- lmFit(assay(spe_ruv3, name), design = design)
}

# "ls" (Least Squares): This is the default method for lmFit. It uses ordinary least squares (OLS) regression to fit a linear model for each gene independently.
# "robust" (Robust Regression): When to use: Robust regression is less sensitive to outliers in the data. If your data contains extreme expression values for some genes that might unduly influence the least squares fit, robust regression can provide more stable and reliable estimates of coefficients and standard errors.
# What it does: Instead of minimizing the sum of squared residuals, it minimizes a robust loss function (e.g., Huber or Tukey's bisquare), giving less weight to large residuals (outliers). This makes the fit less susceptible to individual data points that deviate significantly from the overall trend.

fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

# lmFit(method="robust"): Deals with individual expression value outliers (observation-based). arrayWeights(): Deals with outlier arrays (sample-based).
# eBayes(robust=TRUE): Deals with outlier (hypervariable) genes in the empirical Bayes moderation step.

# results_efit <- decideTests(efit, p.value = 0.025)
# summary_efit <- summary(results_efit)
# print(summary_efit)

de_results_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
dim(de_results_BvT)

print(state)
print(dim(de_results_BvT[which(de_results_BvT$adj.P.Val < 0.025 & (de_results_BvT$logFC) > 0.1),]))
print(dim(de_results_BvT[which(de_results_BvT$adj.P.Val < 0.025 & (de_results_BvT$logFC) < -0.1),]))
print(dim(de_results_BvT[which(de_results_BvT$adj.P.Val < 0.025 & abs(de_results_BvT$logFC) > 0.1),]))
print(dim(de_results_BvT))


#----- intersection with KEGG PD
ND_input = fread('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/ND.txt', stringsAsFactors = F, header = F)
dim(ND_input)
print(length(intersect(ND_input$V2, row.names(de_results_BvT[which(de_results_BvT$adj.P.Val < 0.025 & abs(de_results_BvT$logFC) > 0.1),]))))
#-----

#}

library(data.table)
de_results_BvT$Gene = row.names(de_results_BvT)
fwrite(de_results_BvT, paste0(datadir, Run, "_", query, "_", control, '_GeoMx_DEA.txt'), quote = F, row.names = F, sep = '\t')

## --- "Visual estimation of the FDR with the p-value histogram."
alpha = binw = 0.025
pi0 = 2 * mean(de_results_BvT$adj.P.Val > 0.5)
  
pi0 * alpha / mean(de_results_BvT$adj.P.Val <= alpha)
(pi0 * binw * nrow(de_results_BvT))/table(de_results_BvT$adj.P.Val < alpha)[1]
  
#----
library(ggplot2)
  
sum(summary_efit[1] + summary_efit[2] + summary_efit[3])
  
p = ggplot(as(de_results_BvT, "data.frame"), aes(x = adj.P.Val)) +
    geom_histogram(binwidth = binw, fill = "Royalblue", boundary = 0, alpha=0.6) +
    ggtitle(paste0(names$names2[names$names1 == name], ' (DEGs: ',sum(summary_efit[1] + summary_efit[3]),', Up-regulated: ',summary_efit[3],', Down-regulated: ',summary_efit[1],')')) +
    geom_hline(yintercept = pi0 * binw * nrow(de_results_BvT), color="black", linetype="dashed", size=1) +
    geom_vline(xintercept = alpha, col = "red") +
    annotate("text", x = 0.5, y = 100, label = paste("FP =", round((pi0 * binw * nrow(de_results_BvT)), 0), " FDRe:", round(pi0 * alpha / mean(de_results_BvT$adj.P.Val <= alpha), 4), " FDR:", round((pi0 * binw * nrow(de_results_BvT))/table(de_results_BvT$adj.P.Val < alpha)[1], 4)), color = "black", size = 10) +
    theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"), legend.text=element_text(size=20), legend.title=element_text(size=20), plot.title = element_text(size = 19, face="bold"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "#EEEEEE"),
          panel.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA)) +
    xlab("P-value") + ylab("Count")
  
print(p)

setwd(datadir)
pdf(paste0(state, "_", name, "_", Run, "_PVALUES.pdf"), width = 10, height = 7)
print(p)
dev.off()


## ----------------------------------------------------------------------------- DEG analysis using LMM
if(state == "PDwLBP - PDwoLBP" & LMM)
{
  #============================================ Within Slide Analysis - DEGs
  # https://bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
 
  # Subset by genes and ROIs
  target_Data_subset <- target_Data[fData(target_Data)$TargetName %in% row.names(assay(spe_ruv3, "logcounts")), gsub('.dcc','',row.names(pData(target_Data))) %in% colnames(assay(spe_ruv3, "logcounts"))]
  dim(target_Data_subset)
  table(pData(target_Data_subset)$segment)
  
  assayDataElement(object = target_Data_subset, elt = "exprs2") <- assayDataApply(target_Data_subset, 2, FUN = log, base = 2, elt = "exprs")
  
  # I need the object in NanoStringGeoMxSet but spe_ruv3 is in SpatialExperiment format so I need to repalce the expression values
  for(i in 1:dim(assayDataElement(object = target_Data_subset, elt = "exprs"))[2])
  {
    print(i)
    assayDataElement(object = target_Data_subset, elt = "exprs")[,i] = as.numeric(assay(spe_ruv3, "logcounts")[,i])
  }
  
  identical(as.numeric(assayDataElement(object = target_Data_subset, elt = "exprs")[,i]), as.numeric(assay(spe_ruv3, "logcounts")[,i]))
  
  table(pData(target_Data_subset)$segment)
  target_Data_subset = target_Data_subset[,which(pData(target_Data_subset)$segment != "control")]
  table(pData(target_Data_subset)$segment)
  
  # convert test variables to factors
  pData(target_Data_subset)$testRegion <- factor(pData(target_Data_subset)$segment, c("PD-with-LB", "PD-without-LB"))
  pData(target_Data_subset)[["slide"]] <-  factor(pData(target_Data_subset)[["slide name"]])
  
  # add batch_combined
  pData(target_Data_subset)[["batchslide"]] <-  factor(colData(spe_ruv3)$batchslide)
  pData(target_Data_subset)[["batch_combined"]] <-  factor(colData(spe_ruv3)$batch_combined)
  
  colnames(colData(spe_ruv3))
  # ~ FixedEffect1 (condition) + FixedEffect2 + ... + (1 | BatchFactor1)
  # (1 | random intercept)
  # Crossed vs. Nested Random Effects:
  # (1 | MainBatch:SubBatch)
  # (1 | random slope:random intercept)
  # for independent technical batches like SequencingRun and ProcessingDate, crossed effects are usually appropriate.
  
  # run LMM:
  # formula follows conventions defined by the lme4 package
  mixedOutmc <-
    mixedModelDE(target_Data_subset,
                 elt = "exprs",
                 modelFormula = ~ testRegion + (1 + testRegion | batch_combined),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = TRUE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <-
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  
  # Benjamini-Hochberg (FDR control)
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate",
                       "Pr(>|t|)", "FDR")]

  results <- c()
  results <- rbind(results, r_test)
  table(results$FDR < 0.025 & results$Estimate < 0)
  table(results$FDR < 0.025 & results$Estimate > 0)
  
  fwrite(results, paste0(datadir, 'Within_Slide_Analysis_DEGs_wol_vs_wl.txt'), quote = F, row.names = F, sep = '\t')
  
  colnames(results)
  colnames(de_results_BvT)
  
  de_results_BvT = data.frame(logFC = results$Estimate, adj.P.Val =  results$FDR, row.names = results$Gene)
  
  #====== pvalue plot
  
  ## --- "Visual estimation of the FDR with the p-value histogram."
  alpha = binw = 0.025
  pi0 = 2 * mean(de_results_BvT$adj.P.Val > 0.5)
  
  pi0 * alpha / mean(de_results_BvT$adj.P.Val <= alpha)
  (pi0 * binw * nrow(de_results_BvT))/table(de_results_BvT$adj.P.Val < alpha)[1]
  
  #---- 
  library(ggplot2)
  
  sum(summary_efit[1] + summary_efit[2] + summary_efit[3])
  
  p = ggplot(as(de_results_BvT, "data.frame"), aes(x = adj.P.Val)) +
    geom_histogram(binwidth = binw, fill = "Royalblue", boundary = 0, alpha=0.6) +
    ggtitle(paste0(names$names2[names$names1 == name], ' (DEGs: ',sum(summary_efit[1] + summary_efit[3]),', Up-regulated: ',summary_efit[3],', Down-regulated: ',summary_efit[1],')')) +
    geom_hline(yintercept = pi0 * binw * nrow(de_results_BvT), color="black", linetype="dashed", size=1) +
    geom_vline(xintercept = alpha, col = "red") +
    annotate("text", x = 0.5, y = 100, label = paste("FP =", round((pi0 * binw * nrow(de_results_BvT)), 0), " FDRe:", round(pi0 * alpha / mean(de_results_BvT$adj.P.Val <= alpha), 4), " FDR:", round((pi0 * binw * nrow(de_results_BvT))/table(de_results_BvT$adj.P.Val < alpha)[1], 4)), color = "black", size = 10) +
    theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"), legend.text=element_text(size=20), legend.title=element_text(size=20), plot.title = element_text(size = 19, face="bold"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "#EEEEEE"),
          panel.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA)) +
    xlab("P-value") + ylab("Count")
  
  print(p)
  
  setwd(datadir)
  pdf(paste0(state, "_", name, "_", Run, "_LMM_PVALUES.pdf"), width = 10, height = 7)
  print(p)
  dev.off()
}

#====== volcono plot
library(ggrepel)
# Categorize de_results_BvT based on P-value & FDR for plotting
de_results_BvT$Color <- "NS or FC < 0.5"
de_results_BvT$Color[which(de_results_BvT$adj.P.Val < 0.025)] <- "FDR < 0.025"
de_results_BvT$Color[which(de_results_BvT$adj.P.Val < 0.01 & de_results_BvT$logFC > 0.1)] <- "FDR < 0.01 & up-regulated"
de_results_BvT$Color[which(de_results_BvT$adj.P.Val < 0.01 & de_results_BvT$logFC < -0.1)] <- "FDR < 0.01 & down-regulated"
de_results_BvT$Color[which(abs(de_results_BvT$logFC) < 0.1)] <- "NS or FC < 0.1"
de_results_BvT$Color <- factor(de_results_BvT$Color,
                               levels = c("NS or FC < 0.1", "FDR < 0.025", "FDR < 0.01 & up-regulated", "FDR < 0.01 & down-regulated"))
table(de_results_BvT$Color)
  
# pick top TargetNames for either side of volcano to label
# order TargetNames for convenience:
de_results_BvT$invert_P <- (-log10(de_results_BvT$adj.P.Val)) * sign(de_results_BvT$logFC)
top_g <- c()

colnames(de_results_BvT)
de_results_BvT$TargetName = row.names(de_results_BvT)

top_g <- c(top_g,
           de_results_BvT[, 'TargetName'][
             order(de_results_BvT[, 'invert_P'], decreasing = TRUE)[1:30]],
           de_results_BvT[, 'TargetName'][
             order(de_results_BvT[, 'invert_P'], decreasing = FALSE)[1:30]])

top_g <- unique(top_g)
de_results_BvT <- de_results_BvT[, -which(colnames(de_results_BvT) == 'invert_P')] # remove invert_P from matrix
colnames(de_results_BvT)
table(de_results_BvT$Color)

library(dplyr)
# Graph de_results_BvT
DEG_Plot = ggplot(de_results_BvT,
                  aes(x = logFC, y = -log10(adj.P.Val),
                      color = Color, label = TargetName)) +
  geom_point() +
  geom_vline(xintercept = c(0.1, -0.1), lty = "dashed") +
  geom_hline(yintercept = -log10(0.025), lty = "dashed") +
  ggtitle(name) + 
  labs(x = expression(paste('Log'['2'],' fold change')),
       y = expression(paste('Log'['10'],'P')) ,
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.01 & down-regulated` = "#56B4E9",
                                `FDR < 0.01 & up-regulated` = "#D55E00",
                                `FDR < 0.025` = "lightblue",
                                `P < 0.01` = "orange2",
                                `NS or FC < 0.1` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(de_results_BvT, TargetName %in% top_g & adj.P.Val < 0.025),
                  size = 5, point.padding = 0.2, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom")
print(DEG_Plot)

setwd(datadir)
pdf(paste0(state, "_", name, "_", Run, "_Volcanoplot.pdf"), width = 10, height = 7)
print(DEG_Plot)
dev.off()

#====== heatmap

de_results_BvT_sig = de_results_BvT[order(de_results_BvT$adj.P.Val, decreasing = F),]

upregulated = row.names(de_results_BvT_sig)[which(de_results_BvT_sig$adj.P.Val < 0.01 & de_results_BvT_sig$logFC > 0)]
downregulated = row.names(de_results_BvT_sig)[which(de_results_BvT_sig$adj.P.Val < 0.01 & de_results_BvT_sig$logFC < 0)]

if(state != "PDwLBP - PDwoLBP")
{
  CountMatrix_subset = CountMatrix
}

length(which(row.names(CountMatrix_subset) %in% c(upregulated[1:10], downregulated[1:10]) ))
CountMatrix_subset = CountMatrix_subset[which(row.names(CountMatrix_subset) %in% c(upregulated[1:10], downregulated[1:10]) ),]
dim(CountMatrix_subset)

#------ heatmap
State = as.character(colnames(CountMatrix_subset))
State[grep('V', State)] = 'Vulnerable Neurons'
State[grep('R', State)] = 'Resistant Neurons'
State[grep('C', State)] = 'Healthy Neurons'

annotation_col = data.frame(
  State = factor(State, levels = unique(State))
)
rownames(annotation_col) = as.character(colnames(CountMatrix_subset))
str(annotation_col)

annotation_row = data.frame(
  GeneClass = factor(rep(c('Up-regulated', 'Down-regulated'), c(10, 10)), levels = c('Up-regulated', 'Down-regulated'))
)

rownames(annotation_row) = c(upregulated[1:10], downregulated[1:10])
CountMatrix_subset = CountMatrix_subset[match(rownames(annotation_row), row.names(CountMatrix_subset)),]
identical(row.names(CountMatrix_subset), row.names(annotation_row))

ann_colors = list(
  State = c(S1 = "#1B9E77", S2 = "#D95F02"),
  GeneClass = c(`Up-regulated` = "#7570B3", `Down-regulated` = "#E7298A")
)

names(ann_colors$State)[1] = unique(State)[1]
names(ann_colors$State)[2] = unique(State)[2]

#-- revise colnames
colnames(CountMatrix_subset)[grep('V', colnames(CountMatrix_subset))] = rep(paste0('Donor-', 1:length(grep('V', colnames(CountMatrix_subset))), '-V' ), each = 1)
colnames(CountMatrix_subset)[grep('R', colnames(CountMatrix_subset))] = rep(paste0('Donor-', 1:length(grep('R', colnames(CountMatrix_subset))), '-R' ), each = 1)
colnames(CountMatrix_subset)[grep('C', colnames(CountMatrix_subset))] = rep(paste0('Donor-', 1:length(grep('C', colnames(CountMatrix_subset))), '-C' ), each = 1)
#--

library(hgnc)
library(data.table)
library(ComplexHeatmap)
p <- pheatmap(CountMatrix_subset, annotation_col = annotation_col, annotation_row = annotation_row, 
              annotation_colors = ann_colors, 
              row_split = annotation_row$GeneClass,
              column_split = annotation_col$State, name = "Expression", scale='row',
              fontsize_number=8, display_numbers = round(as.matrix(CountMatrix_subset), digits = 1), cluster_cols=FALSE  )
p

sampleAnnotationFile = fread(paste0(datadir, 'QC_All.txt'), stringsAsFactors = F, header = T)
sampleAnnotationFile = as.data.frame(sampleAnnotationFile)

sampleAnnotationFile[which(sampleAnnotationFile$`slide name` %in% c('12/046 Run 1')), ]
sampleAnnotationFile[which(sampleAnnotationFile$`slide name` %in% c('13/021 Run 1')), ]
sampleAnnotationFile[which(sampleAnnotationFile$`slide name` %in% c('20/009 Run 1')), ]
sampleAnnotationFile[which(sampleAnnotationFile$`slide name` %in% c('PD1101 Run 1')), ]
sampleAnnotationFile[which(sampleAnnotationFile$`slide name` %in% c('PD1310 Run 1')), ]

unique(gsub('.*\\-', '', gsub('-A01.dcc', '', sampleAnnotationFile$NTC_ID))[which(sampleAnnotationFile$NTC > 8000)])

## ----------------------------------------------------------------------------- visualization of DEGs as table
# de_genes_toptable_BvT %>% 
#   dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
#   DT::datatable(caption = '(limma-voom)') %>%
#   DT::formatStyle('logFC',
#                   valueColumns = 'logFC',
#                   backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
#   DT::formatSignif(1:4, digits = 4)

## ----------------------------------------------------------------------------- Enrichment analysis
# library(msigdb)
# library(GSEABase)
# 
# #--- you can inactive this part as it takes long time
# 
# msigdb_hs <- getMsigdb(version = '7.2')
# msigdb_hs <- appendKEGG(msigdb_hs)
# 
# # load PPI from the msigdb package
# # ppi = getIMEX('hs', inferred = TRUE)
# 
# sc <- listSubCollections(msigdb_hs)
# 
# gsc <- c(subsetCollection(msigdb_hs, c('h')),
#          subsetCollection(msigdb_hs, 'c2', sc[grepl("^CP:",sc)]),
#          subsetCollection(msigdb_hs, 'c5', sc[grepl("^GO:",sc)])) %>%
#   GeneSetCollection()
# 
# # Preprocessing is conducted on these genesets, filtering out genesets with less than 5 genes and creating indices vector list for formatting on the results before applying fry.
# GeneList = assay(spe_ruv3, "logcounts")
# 
# GeneList_subset = row.names(GeneList[which(row.names(GeneList) %in% row.names(de_genes_toptable_BvT[which(de_genes_toptable_BvT$adj.P.Val < 1e-5), ])),])
# 
# length(GeneList_subset)
# dim(GeneList)
# 
# fry_indices <- ids2indices(lapply(gsc, geneIds), , remove.empty = FALSE)
# names(fry_indices) <- sapply(gsc, setName)
# 
# gsc_category <- sapply(gsc, function(x) bcCategory(collectionType(x)))
# gsc_category <- gsc_category[sapply(fry_indices, length) > 5]
# 
# gsc_subcategory <- sapply(gsc, function(x) bcSubCategory(collectionType(x)))
# gsc_subcategory <- gsc_subcategory[sapply(fry_indices, length) > 5]
# 
# fry_indices <- fry_indices[sapply(fry_indices, length) > 5]
# 
# names(gsc_category) = names(gsc_subcategory) = names(fry_indices)
# 
# # Now we run fry with all the gene sets we filtered.
# 
# fry_indices_cat <- split(fry_indices, gsc_category[names(fry_indices)])
# 
# fry_res_out <- lapply(fry_indices_cat, function (x) {
#   limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)
# })
# 
# post_fry_format <- function(fry_output, gsc_category, gsc_subcategory){
#   names(fry_output) <- NULL
#   fry_output <- do.call(rbind, fry_output)
#   fry_output$GenesetName <- rownames(fry_output)
#   fry_output$GenesetCat <- gsc_category[rownames(fry_output)]
#   fry_output$GenesetSubCat <- gsc_subcategory[rownames(fry_output)]
#   return(fry_output)
# }
# 
# fry_res_sig <- post_fry_format(fry_res_out, gsc_category, gsc_subcategory) %>%
#   as.data.frame() %>%
#   filter(FDR < 0.05) 
# 
# # The output is a data.frame object. We can either output the whole table, or inspect the top N gene sets in a bar plot.
# # We can see many immune-related gene sets are significantly enriched, B cell-related gene-sets are enriched in up-regulated genes while T-cell related gene-sets are enriched in down-regulated genes.
# 
# #============================================ visualization of enrichment analysis results
# fry_res_sig %>%
#   arrange(FDR) %>%
#   filter(Direction == "Up") %>%
#   .[seq(20),] %>%
#   mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
#   ggplot(aes(GenesetName, -log(FDR))) +
#   geom_bar(stat = "identity", fill = "orange2") +
#   theme_bw() +
#   coord_flip() +
#   ggtitle("Up-regulated")
# 
# fry_res_sig %>%
#   arrange(FDR) %>%
#   filter(Direction == "Down") %>%
#   .[seq(20),] %>%
#   mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
#   ggplot(aes(GenesetName, -log(FDR))) +
#   geom_bar(stat = "identity", fill = "dodgerblue") +
#   theme_bw() +
#   coord_flip() +
#   ggtitle("Down-regulated")
# 
# # Visualization
# # An alternative way to summarise the GSEA output is to visualise common gene sets as a group.
# # We can use the igraph and vissE package to perform clustering on the enriched gene sets and visualise the gene sets using word cloud-based algorithm and network-based visualisation. For more information about vissE, check out here.
# library(vissE)
# library(igraph)
# 
# fry_out = fry_res_sig
# de_table = de_genes_toptable_BvT
# topN = 6
# title = ""
# specific_clusters = NA
# 
# dovissE <- function(fry_out, de_table, topN = 6, title = "", specific_clusters = NA){
#   
#   n_row = min(1000, nrow(fry_out))
#   gs_sig_name <- fry_out %>% 
#     filter(FDR < 0.05) %>%
#     arrange(FDR) %>% 
#     .[1:n_row,] %>% 
#     rownames()
#   gsc_sig <- gsc[as.numeric(gs_sig_name),]
#   
#   gs_ovlap <- computeMsigOverlap(gsc_sig, thresh = 0.15)
#   gs_ovnet <- computeMsigNetwork(gs_ovlap, gsc)
#   
#   max(fry_out[as.numeric(gs_sig_name),]$FDR)
#   
#   gs_stats <- -log10(fry_out[as.numeric(gs_sig_name),]$FDR)
#   names(gs_stats) <- as.numeric(gs_sig_name)
#   
#   #identify clusters
#   grps = cluster_walktrap(gs_ovnet)
#   #extract clustering results
#   grps = groups(grps)
#   #sort by cluster size
#   grps = grps[order(sapply(grps, length), decreasing = TRUE)]
#   
#   # write output
#   output_clusters <- list()
#   for(i in seq(length(grps))){
#     output_clusters[[i]] <- data.frame(geneset = grps[[i]], cluster = paste0("cluster",names(grps)[i]))
#   }
#   output_clusters <<- output_clusters %>% bind_rows()
#   
#   if(is.na(specific_clusters)){
#     grps <- grps[1:topN]
#   } else {
#     grps <- grps[specific_clusters %>% as.character()]
#   }
#   
#   #plot the top 12 clusters
#   set.seed(36) #set seed for reproducible layout
#   p1 <<- plotMsigNetwork(gs_ovnet, markGroups = grps, 
#                          genesetStat = gs_stats, rmUnmarkedGroups = TRUE) +
#     scico::scale_fill_scico(name = "-log10(FDR)")
#   
#   p2 <<- plotMsigWordcloud(gsc, grps, type = 'Name')
#   
#   genes <- unique(unlist(geneIds(gsc_sig)))
#   genes_logfc <- de_table %>% rownames_to_column() %>% filter(rowname %in% genes) %>% .$logFC
#   names(genes_logfc) <- de_table %>% rownames_to_column() %>% filter(rowname %in% genes) %>% .$rowname
#   p3 <<- plotGeneStats(genes_logfc, gsc, grps) +
#     geom_hline(yintercept = 0, colour = 2, lty = 2) +
#     ylab("logFC")
#   
#   # p4 <- plotMsigPPI(ppi, gsc, grps[1:topN], geneStat = genes_logfc) +
#   #  guides(col=guide_legend(title="logFC"))
#   
#   print(p2 + p1 + p3 + patchwork::plot_layout(ncol = 3) +
#           patchwork::plot_annotation(title = title))  
#   
# }
# 
# dovissE(fry_res_sig, de_genes_toptable_BvT, topN = 3, title = "" )


## ----------------------------------------------------------------------------- pathway Enrichment v1

#---- pathview - visualization 

library('org.Hs.eg.db')
keytypes(org.Hs.eg.db)

input = de_results_BvT[which(de_results_BvT$adj.P.Val < 0.05 & abs(de_results_BvT$logFC) > 0.1),] # 
input$TargetName = row.names(input)

input = input[-which(is.na(as.character(mapIds(org.Hs.eg.db, input$TargetName, 'ENTREZID', 'SYMBOL')))),]
row.names(input) = as.character(mapIds(org.Hs.eg.db, input$TargetName, 'ENTREZID', 'SYMBOL'))
dim(input)

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

input_used = input$logFC
names(input_used) = row.names(input)

# pv.out <- pathview(gene.data = input_used, pathway.id = '05012',
#                    species = "hsa", kegg.native = T,
#                    same.layer = F, high=list(gene="gold"), out.suffix = paste0(gsub(' ','',Run), "_", query, "_", control, '_VR'))

list_pathways = c('05010', '05012', '05014', '05016', '05017', '05020', '05022', 'hsa04210', 'hsa04217', 'hsa04140', 'hsa04150', 'hsa00190', 'hsa04137', 'hsa04620', 'hsa04621', 'hsa04060', 'hsa04610', 'hsa04120', 'hsa04141', 'hsa04728', 'hsa04724', 'hsa04727', 'hsa04310', 'hsa04340', 'hsa04722', 'hsa04151', 'hsa04010', 'hsa04020')

for(p in 1:length(list_pathways))
{
  pv.out <- pathview(gene.data = input_used, pathway.id = list_pathways[p],
                     species = "hsa", kegg.native = T,
                     same.layer = F, high = "darkred", low = "dodgerblue", out.suffix = paste0(gsub(' ','',Run), "_", query, "_", control, '_VR'))
}

## ----------------------------------------------------------------------------- pathway Enrichment v2

#--- gene enrichment analysis
library('org.Hs.eg.db')
BD = 'h.all.v2024.1.Hs.entrez.gmt'

library(qusage)
library(clusterProfiler)
gmtfile <- paste0("/Users/isarnassiri/Documents/GeoMx_Analysis/reference_datasets/msigdb_v2024.1.Hs_GMTs/", BD)
c5 <- read.gmt(gmtfile)

# use mapIds method to obtain Entrez IDs
# input = de_results_BvT[de_results_BvT$adj.P.Val < 0.05,]
de = as.character(mapIds(org.Hs.eg.db, input$TargetName, 'ENTREZID', 'SYMBOL'))

table(c5$term)
c5$term = sub("^[^_]*_", "", c5$term)

egmt <- enricher(de[!is.na(de)], TERM2GENE=c5)

if(as.numeric(table(egmt@result$p.adjust < 0.05)['TRUE'])>5)
{
  library(cowplot)
  library(ggplot2)
  p2 <- dotplot(egmt, showCategory=5) + ggtitle("")
  p2
  
  pdf(paste0(gsub(' ','',Run), "_", query, "_", control, '_PathwayEnrichment.pdf'), width = 10, height = 10, useDingbats = FALSE)
  print(plot_grid(p2, ncol=1))
  dev.off()
  
  fwrite(egmt@result, paste0(gsub(' ','',Run), "_", query, "_", control, '_PathwayEnrichment.txt') , quote = F, sep = '\t', row.names = F)
}

#--- gene enrichment analysis
library("DOSE")
edo <- enrichDGN(de[!is.na(de)])

#--- grep multiple pattern
toMatch <- c("Disease", "Disorders", "syndrome", "Parkin")

# OR
matches <- unique(grep(paste(toMatch,collapse="|"), edo@result$Description, value=TRUE))

#-- keep terms related to diseases
edo@result = edo@result[which(edo@result$Description %in% matches),]

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

# View(edo@result)

library(enrichplot)
# barplot(edo, showCategory=10)
# dotplot(edo, showCategory=15) + ggtitle(print(paste0(  S )))
 
if(as.numeric(table(edo@result$p.adjust < 0.05)['TRUE'])>5)
{
  pdf(paste0(gsub(' ','',Run), "_", query, "_", control, '_DiseaseEnrichment_dotplot.pdf'),width = 8, height = 9, useDingbats = FALSE)
  print(dotplot(edo, showCategory=15))
  dev.off()
  
  pdf(paste0(gsub(' ','',Run), "_", query, "_", control, '_DiseaseEnrichment_cnetplot.pdf'),width = 8, height = 9, useDingbats = FALSE)
  print(cnetplot(edox, circular = TRUE, colorEdge = TRUE))
  dev.off()
  
  fwrite(edo@result, paste0(gsub(' ','',Run), "_", query, "_", control, '_DiseaseEnrichment.txt') , quote = F, sep = '\t', row.names = F)
}
}




## ----------------------------------------------------------------------------- downstream analysis

## ----------------------------------------------------------------------------- intersection of DEGs per state
setwd('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4/')

listFiles = list.files(pattern = '_GeoMx_DEA.txt')

#--- gene enrichment analysis
library(data.table)
library('org.Hs.eg.db')
library(qusage)
library(clusterProfiler)
library(cowplot)
library(ggplot2)
library("DOSE")

#--
BD = 'c5.all.v2023.1.Hs.entrez.gmt'
gmtfile <- paste0("/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Scripts_USED_visium/msigdb_v2023_Hs/", BD)
c5 <- read.gmt(gmtfile)
c5$term = sub("^[^_]*_", "", c5$term)
#--

for(filename in listFiles)
{
  print(gsub('Run 1_|_GeoMx_DEA.txt','',filename))
  Temp = data.frame(fread(filename, stringsAsFactors = F, header = T))
  Temp = Temp[which(Temp$adj.P.Val < 0.025),]
  print(dim(Temp))
  
  assign(gsub('Run 1_|_GeoMx_DEA.txt','',filename), Temp$TargetName)
  
  input_used = Temp$logFC
  names(input_used) = Temp$TargetName
  assign(paste0(gsub('Run 1_|_GeoMx_DEA.txt','',filename), '_KEGG'), input_used)
  
  #---------------------- Pathway enrichment
  
  # use mapIds method to obtain Entrez IDs
  input = Temp$TargetName
  de = as.character(mapIds(org.Hs.eg.db, input, 'ENTREZID', 'SYMBOL'))
  
  egmt <- enricher(de[!is.na(de)], TERM2GENE=c5)
  
  if(length(unique(egmt@result$p.adjust < 0.05))>1)
  {
    
    #-- all
    p2 <- dotplot(egmt, showCategory=5) + ggtitle(gsub('Run 1_|_GeoMx_DEA.txt','',filename))
    
    pdf(paste0(gsub('Run 1_|_GeoMx_DEA.txt','',filename), '_setDiff_PathwayEnrichment.pdf'), width = 8, height = 10, useDingbats = FALSE)
    print(plot_grid(p2, ncol=1))
    dev.off()
    
    fwrite(egmt@result, paste0(gsub('Run 1_|_GeoMx_DEA.txt','',filename), '_setDiff_PathwayEnrichment.txt'), quote = F, sep = '\t', row.names = F)
    
    #-- disease
    edo <- enrichDGN(de[!is.na(de)])
    edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
    
    pdf(paste0(gsub('Run 1_|_GeoMx_DEA.txt','',filename), '_setDiff_DiseaseEnrichment.pdf'), width = 8, height = 10, useDingbats = FALSE)
    print(dotplot(edo, showCategory=10))
    dev.off()
    
    fwrite(edo@result, paste0(gsub('Run 1_|_GeoMx_DEA.txt','',filename), '_setDiff_DiseaseEnrichment.txt'), quote = F, sep = '\t', row.names = F)
  }
  
  # x2 <- list(
  #   `PD-without-LB vs Control` = get('PDwoLBP_control_KEGG'), 
  #   `PD-with-LB vs Control` = get('PDwLBP_control_KEGG'), 
  #   `PD-with-LB vs PD-without-LB` = get('PDwLBP_PDwoLBP_KEGG')
  # )
  # 
  # list_pathways = c('05010', '05012', '05014', '05016', '05017', '05020', '05022')
  # 
  # for(p in 1:length(list_pathways))
  # {
  #   pv.out <- pathview(gene.data = get('PDwLBP_PDwoLBP_KEGG'), pathway.id = list_pathways[p],
  #                      species = "hsa", kegg.native = T,
  #                      same.layer = F, high=list(gene="gold"), out.suffix = 'PDwoLBP_control_KEGG')
  # }
}

#---------------- Venn Diagram
x <- list(
  `PD-without-LB vs Control` = get('PDwoLBP_control'), 
  `PD-with-LB vs Control` = get('PDwLBP_control'), 
  `PD-with-LB vs PD-without-LB` = get('PDwLBP_PDwoLBP')
)

str(x)

# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)
pdf(paste0("VennDiagram_DEGs.pdf"), width = 10, height = 7)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

dev.off()

#devtools::install_github("vqf/nVennR")
# library(nVennR)
# plotVenn(x, outFile="a.svg")

#----------------


#--------------------------------- compare with Variance-defining genes in PCA (n 311)


x <- list()

#--- PD_KEGG
setwd('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script')
library(data.table)
PD_KEGG = fread("PD.txt", stringsAsFactors = F, header = T)
PD_KEGG = data.frame(PD_KEGG)
colnames(PD_KEGG)[2] = 'Gene.ID'

x[['Genes associated with PD']] = PD_KEGG$Gene.ID

#--- 5k
setwd('/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script')
Panel5k = fread("Panel5k.txt", stringsAsFactors = F, header = T)
Panel5k = data.frame(Panel5k)
colnames(Panel5k)[1] = 'Gene.ID'

x[['Genes in Xenium 5k panel']] = Panel5k$Gene.ID

#--- Drug targets
Drugable = fread("Drugable.txt", stringsAsFactors = F, header = T)
Drugable = data.frame(Drugable)
colnames(Drugable)[1] = 'Gene.ID'

x[['Druggable Genes']] = Drugable$Gene.ID

#--- Potential Drug targets
Drugable_Potential = fread("Drugable_Potential.txt", stringsAsFactors = F, header = T)
Drugable_Potential = data.frame(Drugable_Potential)
colnames(Drugable_Potential)[1] = 'Gene.ID'

x[['Potential Druggable Genes']] = Drugable_Potential$Gene.ID


setwd("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4")
selected_variables_percentile = fread("logcounts_selected_variables_percentile.txt", stringsAsFactors = F, header = F)

x[['PCA-associated genes']] = selected_variables_percentile$V1

print(length(intersect(Panel5k$Gene.ID, selected_variables_percentile$V1)))
print(length(intersect(Panel5k$Gene.ID, selected_variables_percentile$V1 ))/length(selected_variables_percentile$V1))

#--- GeoMx
setwd('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4/')
library(data.table)

# library(readxl)
# GeoMx_metadata = read_excel(paste0(datadir, 'annotation/Annotation_File.xlsx'))
# dim(GeoMx_metadata)

filenames =  list.files('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4/', pattern = '_GeoMx_DEA.txt')
require(data.table) ## 1.9.2 or 1.9.3
GeoMx_DEG = rbindlist(lapply(filenames, fread))
GeoMx_DEG = GeoMx_DEG[which(GeoMx_DEG$adj.P.Val < 0.01), ]
GeoMx_DEG = GeoMx_DEG[order(GeoMx_DEG$adj.P.Val, decreasing = F),]
GeoMx_DEG = GeoMx_DEG[!duplicated(GeoMx_DEG$TargetName),]
dim(GeoMx_DEG)

x[['GeoMx DEGs']] = GeoMx_DEG$TargetName

pdf(paste0("/Users/isarnassiri/Documents/Xenium/Custom_Panel_and_Script/PCA-associated-genes.pdf"), width = 10, height = 7)

names(x)

library(ggvenn)
print(ggvenn(
  x[c(1:2,5:6)], 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
))

dev.off()

for(n in names(x) )
{
  print(n)
  
  print(length(intersect(x[["PCA-associated genes"]], x[[n]])))
  print(length(intersect(x[["PCA-associated genes"]], x[[n]]))/length(x[["PCA-associated genes"]]))
  
}

print(length(intersect(Panel5k$Gene.ID, selected_variables_percentile$V1)))
print(length(intersect(Panel5k$Gene.ID, selected_variables_percentile$V1 ))/length(selected_variables_percentile$V1))













## ----------------------------------------------------------------------------- psuedotime analysis

colData(spe2_subset)$state
spe_ruv@assays
spe_ruv@assays@data@listData$logcounts
spe_ruv@colData$segment
table(colData(spe_ruv)$segment)
names(spe_ruv@assays)
molecules(spe_ruv)

library(data.table)
fwrite(de_results_BvT, paste0(datadir, Run, "_", query, "_", control, '_GeoMx_DEA.txt'), quote = F, row.names = F, sep = '\t')



#------------------------------------------ consider the SAVER approach for data normalization
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

# library(data.table)
# datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/")
# count_mat2 = fread(paste0(datadir, 'GeoMx_readCount.txt'), stringsAsFactors = F, header = T)
# count_mat2 = as.data.frame(count_mat2)
# row.names(count_mat2) = make.names(count_mat2$GeneSymbol,unique=T) 
# count_mat2 = count_mat2[,-1]
# colnames(count_mat2) = gsub('\\.','-',gsub('.dcc','',colnames(count_mat2)))
# 
# count_mat3 = fread(paste0(datadir, 'GeoMx_readCount_TMMRUV4.txt'), stringsAsFactors = F, header = T)
# count_mat3 = as.data.frame(count_mat3)
# row.names(count_mat3) = make.names(count_mat3$GeneSymbol,unique=T) 
# count_mat3 = count_mat3[,-1]
# colnames(count_mat3) = gsub('\\.','-',gsub('.dcc','',colnames(count_mat3)))
# 
# count_mat2 = count_mat2[which(row.names(count_mat2) %in% row.names(count_mat3)), which(colnames(count_mat2) %in% colnames(count_mat3))]
# dim(count_mat2)
# 
# MD = fread(paste0(datadir, 'MetaData.txt'), stringsAsFactors = F, header = T)
# MD = as.data.frame(MD)
# colnames(MD)
# 
# Reads <- apply(count_mat2[,which(colnames(count_mat2) %in% MD$SampleID)], 2, function(x) sum(x))
# #-- nGenesPerCell
# nGenesPerCell <- apply(count_mat2, 2, function(x) sum(x>0))
# 
# identical(names(Reads), MD$SampleID)
# identical(names(nGenesPerCell), MD$SampleID)
# 
# plot(as.numeric(nGenesPerCell), MD$GenesDetected)
# 
# dataset = data.frame(Nuclei_Count = MD$NucleiCount, Genes_Detected = nGenesPerCell,
#                      Area_Size = MD$area,
#                      Reads_Count = as.numeric(Reads) )
# 
# colnames(dataset) = c("Nuclei Count", "Genes Detected", "Area Size", "Reads Count")
# variables <- c("Nuclei Count", "Genes Detected", "Area Size", "Reads Count")
# str(dataset)
# 
# # ---------------------------------- 
# 
# .plotMarginalCor <- function(variable, cexYlab = 1.3, lwd = 2, rugs = FALSE) {
#   
#   # histogram with density estimator
#   
#   variable <- variable[!is.na(variable)]
#   
#   density <- density(variable)
#   h <- hist(variable, plot = FALSE)
#   jitVar <- jitter(variable)
#   yhigh <- max(max(h$density), max(density$y))
#   ylow <- 0
#   xticks <- pretty(c(variable, h$breaks), min.n = 3)
#   plot(range(xticks), c(ylow, yhigh), type = "n", axes = FALSE, ylab = "", 
#        xlab = "")
#   h <- hist(variable, freq = FALSE, main = "", ylim = c(ylow, yhigh), xlab = "", 
#             ylab = " ", axes = FALSE, col = "grey", add = TRUE, nbreaks = round(length(variable)/5))
#   ax1 <- axis(1, line = 0.3, at = xticks, lab = xticks)
#   par(las = 0)
#   ax2 <- axis(2, at = c(0, max(max(h$density), max(density$y))/2, max(max(h$density), 
#                                                                       max(density$y))), labels = c("", "Density", ""), lwd.ticks = 0, 
#               pos = range(ax1) - 0.08 * diff(range(ax1)), cex.axis = 2, mgp = c(3, 0.7, 0))
#   
#   if (rugs) 
#     rug(jitVar)
#   
#   lines(density$x[density$x >= min(ax1) & density$x <= max(ax1)], density$y[density$x >= 
#                                                                               min(ax1) & density$x <= max(ax1)], lwd = lwd)
# }
# 
# .poly.pred <- function(fit, line = FALSE, xMin, xMax, lwd) {
#   
#   # predictions of fitted model
#   
#   # create function formula
#   f <- vector("character", 0)
#   
#   for (i in seq_along(coef(fit))) {
#     
#     if (i == 1) {
#       
#       temp <- paste(coef(fit)[[i]])
#       f <- paste(f, temp, sep = "")
#     }
#     
#     if (i > 1) {
#       
#       temp <- paste("(", coef(fit)[[i]], ")*", "x^", i - 1, sep = "")
#       f <- paste(f, temp, sep = "+")
#     }
#   }
#   
#   x <- seq(xMin, xMax, length.out = 100)
#   predY <- eval(parse(text = f))
#   
#   if (line == FALSE) {
#     return(predY)
#   }
#   
#   if (line) {
#     lines(x, predY, lwd = lwd)
#   }
# }
# 
# 
# .plotScatter <- function(xVar, yVar, cexPoints = 1.3, cexXAxis = 1.3, 
#                          cexYAxis = 1.3, lwd = 2) {
#   
#   # displays scatterplot
#   
#   d <- data.frame(xx = xVar, yy = yVar)
#   d <- na.omit(d)
#   xVar <- d$xx
#   yVar <- d$yy
#   
#   fit <- lm(yy ~ poly(xx, 1, raw = TRUE), d)
#   
#   xlow <- min((min(xVar) - 0.1 * min(xVar)), min(pretty(xVar)))
#   xhigh <- max((max(xVar) + 0.1 * max(xVar)), max(pretty(xVar)))
#   xticks <- pretty(c(xlow, xhigh))
#   ylow <- min((min(yVar) - 0.1 * min(yVar)), min(pretty(yVar)), min(.poly.pred(fit, 
#                                                                                line = FALSE, xMin = xticks[1], xMax = xticks[length(xticks)], 
#                                                                                lwd = lwd)))
#   yhigh <- max((max(yVar) + 0.1 * max(yVar)), max(pretty(yVar)), max(.poly.pred(fit, 
#                                                                                 line = FALSE, xMin = xticks[1], xMax = xticks[length(xticks)], 
#                                                                                 lwd = lwd)))
#   
#   
#   yticks <- pretty(c(ylow, yhigh))
#   
#   yLabs <- vector("character", length(yticks))
#   
#   for (i in seq_along(yticks)) {
#     
#     if (yticks[i] < 10^6) {
#       
#       yLabs[i] <- format(yticks[i], digits = 3, scientific = FALSE)
#       
#     } else {
#       
#       yLabs[i] <- format(yticks[i], digits = 3, scientific = TRUE)
#     }
#   }
#   
#   plot(xVar, yVar, col = "black", pch = 21, bg = "grey", ylab = "", 
#        xlab = "", axes = FALSE, ylim = range(yticks), xlim = range(xticks), 
#        cex = cexPoints)
#   .poly.pred(fit, line = TRUE, xMin = xticks[1], xMax = xticks[length(xticks)], 
#              lwd = lwd)
#   
#   par(las = 1)
#   
#   axis(1, line = 0.4, labels = xticks, at = xticks, cex.axis = cexXAxis)
#   axis(2, line = 0.2, labels = yLabs, at = yticks, cex.axis = cexYAxis)
#   
#   invisible(max(nchar(yLabs)))
# }
# 
# .plotCorValue <- function(xVar, yVar, cexText = 1.7, cexCI = 1.7, hypothesis = "correlated", 
#                           pearson = TRUE, kendallsTauB = FALSE, spearman = FALSE, confidenceInterval = 0.95) {
#   
#   # displays correlation value
#   
#   CIPossible <- TRUE
#   
#   tests <- c()
#   
#   if (pearson) 
#     tests <- c(tests, "pearson")
#   
#   if (spearman) 
#     tests <- c(tests, "spearman")
#   
#   if (kendallsTauB) 
#     tests <- c(tests, "kendall")
#   
#   plot(1, 1, type = "n", axes = FALSE, ylab = "", xlab = "")
#   
#   lab <- vector("list")
#   
#   for (i in seq_along(tests)) {
#     
#     if (round(cor.test(xVar, yVar, method = tests[i])$estimate, 8) == 
#         1) {
#       
#       CIPossible <- FALSE
#       
#       if (tests[i] == "pearson") {
#         lab[[i]] <- bquote(italic(r) == "1.000")
#       }
#       
#       if (tests[i] == "spearman") {
#         lab[[i]] <- bquote(italic(rho) == "1.000")
#       }
#       
#       if (tests[i] == "kendall") {
#         lab[[i]] <- bquote(italic(tau) == "1.000")
#       }
#       
#     } else if (round(cor.test(xVar, yVar, method = tests[i])$estimate, 
#                      8) == -1) {
#       
#       CIPossible <- FALSE
#       
#       if (tests[i] == "pearson") {
#         lab[[i]] <- bquote(italic(r) == "-1.000")
#       }
#       
#       if (tests[i] == "spearman") {
#         lab[[i]] <- bquote(italic(rho) == "-1.000")
#       }
#       
#       if (tests[i] == "kendall") {
#         lab[[i]] <- bquote(italic(tau) == "-1.000")
#       }
#       
#     } else {
#       
#       if (tests[i] == "pearson") {
#         lab[[i]] <- bquote(italic(r) == .(formatC(round(cor.test(xVar, 
#                                                                  yVar, method = tests[i])$estimate, 3), format = "f", 
#                                                   digits = 3)))
#       }
#       
#       if (tests[i] == "spearman") {
#         lab[[i]] <- bquote(rho == .(formatC(round(cor.test(xVar, 
#                                                            yVar, method = tests[i])$estimate, 3), format = "f", 
#                                             digits = 3)))
#       }
#       
#       if (tests[i] == "kendall") {
#         lab[[i]] <- bquote(tau == .(formatC(round(cor.test(xVar, 
#                                                            yVar, method = tests[i])$estimate, 3), format = "f", 
#                                             digits = 3)))
#       }
#     }
#   }
#   
#   if (length(tests) == 1) {
#     ypos <- 1
#   }
#   
#   if (length(tests) == 2) {
#     ypos <- c(1.1, 0.9)
#   }
#   
#   if (length(tests) == 3) {
#     ypos <- c(1.2, 1, 0.8)
#   }
#   
#   
#   for (i in seq_along(tests)) {
#     
#     text(1, ypos[i], labels = lab[[i]], cex = cexText)
#   }
#   
#   
#   if (hypothesis == "correlated" & length(tests) == 1 & any(tests == 
#                                                             "pearson")) {
#     
#     alternative <- "two.sided"
#     ctest <- cor.test(xVar, yVar, method = tests, conf.level = confidenceInterval)
#   }
#   
#   if (hypothesis != "correlated" & length(tests) == 1 & any(tests == 
#                                                             "pearson")) {
#     
#     if (hypothesis == "correlatedPositively") {
#       
#       ctest <- cor.test(xVar, yVar, method = tests, alternative = "greater", 
#                         conf.level = confidenceInterval)
#       
#     } else if (hypothesis == "correlatedNegatively") {
#       
#       ctest <- cor.test(xVar, yVar, method = tests, alternative = "less", 
#                         conf.level = confidenceInterval)
#     }
#     
#   }
#   
#   if (any(tests == "pearson") & length(tests) == 1 && CIPossible) {
#     
#     CIlow <- formatC(round(ctest$conf.int[1], 3), format = "f", digits = 3)
#     CIhigh <- formatC(round(ctest$conf.int[2], 3), format = "f", digits = 3)
#     
#     text(1, 0.7, labels = paste(100 * confidenceInterval, "% CI: [", 
#                                 CIlow, ", ", CIhigh, "]", sep = ""), cex = cexCI)
#   }
#   
# }
# 
# ### matrix plot ###
# l <- length(variables)
# 
# par(mfrow = c(l, l), cex.axis = 1.3, mar = c(3, 4, 2, 1.5) + 0.1, oma = c(0, 2.2, 2, 0))
# 
# for (row in seq_len(l)) {
#   
#   for (col in seq_len(l)) {
#     
#     if (row == col) {
#       .plotMarginalCor(dataset[[variables[row]]])  # plot marginal (histogram with density estimator)
#     }
#     if (col > row) {
#       .plotScatter(dataset[[variables[col]]], dataset[[variables[row]]])  # plot scatterplot
#     }
#     if (col < row) {
#       if (l < 7) {
#         .plotCorValue(dataset[[variables[col]]], dataset[[variables[row]]], 
#                       cexCI = 1.2)  # plot r= ...
#       }
#       if (l >= 7) {
#         .plotCorValue(dataset[[variables[col]]], dataset[[variables[row]]], 
#                       cexCI = 1.2)
#       }
#     }
#   }
# }
# 
# textpos <- seq(1/(l * 2), (l * 2 - 1)/(l * 2), 2/(l * 2))
# for (t in seq_along(textpos)) {
#   mtext(text = variables[t], side = 3, outer = TRUE, at = textpos[t], 
#         cex = 0.9, line = -0.8)
#   mtext(text = variables[t], side = 2, outer = TRUE, at = rev(textpos)[t], 
#         cex = 0.9, line = -0.1)
# }
# 
# 
# #------------ comparison of GeoMx and snRNA-seq DEGs after integration
# 
# GeoMx_DEGs = fread('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/GeoMx_DEA.txt', stringsAsFactors = F)
# GeoMx_DEGs = as.data.frame(GeoMx_DEGs)
# dim(GeoMx_DEGs)
# 
# GeoMx_DEGs_sig = GeoMx_DEGs[which(GeoMx_DEGs$adj.P.Val < 0.05),]
# dim(GeoMx_DEGs_sig)
# 
# GeoMx_DEGs_up = GeoMx_DEGs[which(GeoMx_DEGs$logFC > 0),]
# dim(GeoMx_DEGs_up)
# 
# GeoMx_DEGs_down = GeoMx_DEGs[which(GeoMx_DEGs$logFC < 0),]
# dim(GeoMx_DEGs_down)
# 
# 
# snRNAseq_DEGs = fread('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/Annotated_DEGenes.txt', stringsAsFactors = F)
# snRNAseq_DEGs = as.data.frame(snRNAseq_DEGs)
# dim(snRNAseq_DEGs)
# 
# snRNAseq_DEGs_sig = snRNAseq_DEGs[which(snRNAseq_DEGs$padj < 0.05),]
# dim(snRNAseq_DEGs_sig)
# 
# snRNAseq_DEGs_up = snRNAseq_DEGs_sig[which(snRNAseq_DEGs_sig$log2FoldChange > 0),]
# dim(snRNAseq_DEGs_up)
# 
# snRNAseq_DEGs_down = snRNAseq_DEGs_sig[which(snRNAseq_DEGs_sig$log2FoldChange < 0),]
# dim(snRNAseq_DEGs_down)
# 
# 
# length(intersect(GeoMx_DEGs_sig$TargetName, snRNAseq_DEGs_sig$feature))
# length(intersect(GeoMx_DEGs_up$TargetName, snRNAseq_DEGs_up$feature))
# length(intersect(GeoMx_DEGs_down$TargetName, snRNAseq_DEGs_down$feature))
# 
# 
# 
# 
# #============================================ test ======================================================================================== 
# library(GeomxTools)
# library(GeoMxWorkflows)
# library(NanoStringNCTools)
# library(tidyverse)
# 
# ## Introduction
# # This summarises the QC of the GeoMx pilot run of 10th October, utilising the 'GeoMxTools' pipeline.
# 
# #============================================ Collate files
# datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/")
# setwd(datadir)
# 
# DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
#                 full.names = TRUE, recursive = TRUE)
# PKCFiles <- dir(file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/pkcs"), pattern = ".pkc$",
#                 full.names = TRUE, recursive = TRUE)
# 
# # I added NucleiCount manually to the following file from Annie Initial Dataset.xlsx file in the 'DSPDA Files' folder
# SampleAnnotationFile <-
#   dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
#       full.names = TRUE, recursive = TRUE)
# 
# library(readxl)
# 
# #============================================ load data
# Data <-
#   readNanoStringGeoMxSet(dccFiles = DCCFiles,
#                          pkcFiles = PKCFiles,
#                          phenoDataFile = SampleAnnotationFile,
#                          phenoDataSheet = excel_sheets(path = SampleAnnotationFile),
#                          phenoDataDccColName = "Sample_ID",
#                          protocolDataColNames = c("aoi","roi","state", "NucleiCount"),
#                          experimentDataColNames = c("panel"))
# 
# #============================================ summary of empty samples
# library(data.table)
# excluded = fread('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/excluded.txt', stringsAsFactors = F, header = F)
# colnames(excluded) = 'Sample_ID'
# excluded$Sample_ID = gsub('.dcc','',excluded$Sample_ID)
# 
# SamplesMetadata = read_excel('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/annotation/Annotation_File.xlsx')
# 
# intersect(SamplesMetadata$Sample_ID, excluded$Sample_ID)
# 
# excluded_SamplesMetadata = merge(SamplesMetadata, excluded, by = 'Sample_ID')
# #============================================
# 
# library(knitr)
# pkcs <- annotation(Data)
# modules <- gsub(".pkc", "", pkcs)
# kable(data.frame(PKCs = pkcs, modules = modules))
# 
# #============================================ Sample Overview
# # Serial sections of a single brain section were analysed, with 7 ROI from BP+TH+ cell bodies and 7 ROI from BP+ cell bodies.
# 
# library(dplyr)
# library(ggforce)
# count_mat = data.frame(GeneSymbol = Data@featureData@data$TargetName, Data@assayData$exprs)
# 
# fwrite(count_mat, paste0(datadir, 'GeoMx_readCount.txt'), quote = F, row.names = F, sep = '\t')
# dim(count_mat)
# # summary(count_mat[,-1])
# 
# ## Preprocessing and Segment QC
# # Before we begin, we will shift any expression counts with a value of 0 to 1 to enable in downstream transformations.
# # We then assess sequencing quality and adequate tissue sampling for every ROI/AOI segment.
# 
# # Every ROI/AOI segment will be tested for:
# #   - Raw sequencing reads: segments with >1000 raw reads are removed.
# # - % Aligned,% Trimmed, or % Stitched sequencing reads: segments below ~80% for one or more of these QC parameters are removed.
# # - % Sequencing saturation ([1-deduplicated reads/aligned reads]%): segments below ~50% require additional sequencing to capture full sample diversity and are not typically analyzed until improved.
# # - Negative Count: this is the geometric mean of the several unique negative probes in the GeoMx panel that do not target mRNA and establish the background count level per segment; segments with low negative counts (1-10) are not necessarily removed but may be studied closer for low endogenous gene signal and/or insufficient tissue sampling.
# # - No Template Control (NTC) count: values >1,000 could indicate contamination for the segments associated with this NTC; however, in cases where the NTC count is between 1,000- 10,000, the segments may be used if the NTC data is uniformly low (e.g. 0-2 counts for all probes).
# # - Nuclei: >100 nuclei per segment is generally recommended; however, this cutoff is highly study/tissue dependent and may need to be reduced; what is most important is consistency in the nuclei distribution for segments within the study.
# # - Area: generally correlates with nuclei; a strict cutoff is not generally applied based on area.
# # 
# # In this case the nuclei QC metric was ignored.
# 
# Data <- shiftCountsOne(Data, useDALogic = TRUE)
# 
# summary(Data@assayData$exprs)
# 
# # QC_params <-
# #   list(minSegmentReads = 1000, # Minimum number of reads (1000)
# #        percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
# #        percentStitched = 80,   # Minimum % of reads stitched (80%)
# #        percentAligned = 75,    # Minimum % of reads aligned (80%)
# #        percentSaturation = 50, # Minimum sequencing saturation (50%)
# #        minNegativeCount = 1,   # Minimum negative control counts (10)
# #        maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
# #        minNuclei = 20,         # Minimum # of nuclei estimated (100)
# #        minArea = 1000)         # Minimum segment area (5000)
# 
# QC_params <-
#   list(minSegmentReads = 1000, # Minimum number of reads (1000)
#        percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
#        percentStitched = 80,   # Minimum % of reads stitched (80%)
#        percentAligned = 80,    # Minimum % of reads aligned (80%)
#        percentSaturation = 50, # Minimum sequencing saturation (50%)
#        minNegativeCount = 10,   # Minimum negative control counts (10)
#        maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
#        minNuclei = 100,         # Minimum # of nuclei estimated (100)
#        minArea = 5000)         # Minimum segment area (5000)
# 
# Data <- setSegmentQCFlags(Data, qcCutoffs = QC_params)        
# 
# QC_Table = data.frame(Aligned = sData(Data)$Aligned, DeduplicatedReads = sData(Data)$DeduplicatedReads)
# # fwrite(QC_Table, paste0(datadir, 'QC_Results.txt'), quote = F, row.names = F, sep = '\t')
# 
# #============================================ Collate QC Results
# QCResults <- protocolData(Data)[["QCFlags"]]
# flag_columns <- colnames(QCResults)
# QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
#                          Warning = colSums(QCResults[, flag_columns]))
# 
# QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
#   ifelse(sum(x) == 0L, "PASS", "WARNING")
# })
# 
# QC_Summary["TOTAL FLAGS", ] <-
#   c(sum(QCResults[, "QCStatus"] == "PASS"),
#     sum(QCResults[, "QCStatus"] == "WARNING"))
# 
# col_by <- "segment"
# 
# #============================================ Graphical summaries of QC statistics plot function
# QC_histogram <- function(assay_data = NULL,
#                          annotation = NULL,
#                          fill_by = NULL,
#                          thr = NULL,
#                          scale_trans = NULL) {
#   plt <- ggplot(assay_data,
#                 aes_string(x = paste0("unlist(`", annotation, "`)"),
#                            fill = fill_by)) +
#     geom_histogram(bins = 50) +
#     geom_vline(xintercept = thr, lty = "dashed", color = "black") +
#     theme_bw() + guides(fill = "none") +
#     facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
#     labs(x = annotation, y = "Segments, #", title = annotation)
#   if(!is.null(scale_trans)) {
#     plt <- plt +
#       scale_x_continuous(trans = scale_trans)
#   }
#   plt
# }
# 
# QC_histogram(sData(Data), "NucleiCount", col_by, 20)
# 
# QC_histogram(sData(Data), "Trimmed (%)", col_by, 80)
# 
# # Stitching: These fragments are computationally merged into a single, longer read using specialized algorithms. This is essential because short reads alone might not uniquely map to the target gene sequences.
# QC_histogram(sData(Data), "Stitched (%)", col_by, 80)
# 
# QC_histogram(sData(Data), "Aligned (%)", col_by, 75)
# 
# QC_histogram(sData(Data), "Saturated (%)", col_by, 50) +
#   labs(title = "Sequencing Saturation (%)",
#        x = "Sequencing Saturation (%)")
# 
# QC_histogram(sData(Data), "area", col_by, 1000, scale_trans = "log10")
# 
# ## Generate negative probe count
# # calculate the negative geometric means for each module
# negativeGeoMeans <- 
#   esBy(negativeControlSubset(Data), 
#        GROUP = "Module", 
#        FUN = function(x) { 
#          assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
#        }) 
# protocolData(Data)[["NegGeoMean"]] <- negativeGeoMeans
# 
# # explicitly copy the Negative geoMeans from sData to pData & plot
# negCols <- paste0("NegGeoMean_", modules)
# pData(Data)[, negCols] <- sData(Data)[["NegGeoMean"]]
# for(ann in negCols) {
#   plt <- QC_histogram(pData(Data), ann, col_by, 2, scale_trans = "log10")
#   print(plt)
# }
# 
# # detatch neg_geomean columns ahead of aggregateCounts call
# pData(Data) <- pData(Data)[, !colnames(pData(Data)) %in% negCols]
# 
# #============================================ QC as tables
# 
# # show all NTC values, Freq = # of Segments with a given NTC count:
# kable(table(NTC_Count = sData(Data)$NTC), col.names = c("NTC Count", "# of Segments"))
# 
# ## ----QCSummaryTable, results = "asis"-----------------------------------------
# kable(QC_Summary, caption = "QC Summary Table for each Segment")
# #============================================
# 
# #============================================ explore NTC  
# highNTC = gsub('.dcc','',row.names(QCResults[QCResults$HighNTC,]))
# 
# dccFiles <- dir('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/dccs/', pattern=".dcc$", full.names=TRUE)
# dccData <- sapply(dccFiles, readDccFile, simplify = FALSE)
# 
# names(dccData) = gsub('.*\\//|.dcc','',names(dccData))
# 
# colnames(SamplesMetadata)
# Sample_ID_NTC = SamplesMetadata[which(SamplesMetadata$`slide name` == 'No Template Control'),"Sample_ID"]
# 
# string = gsub('A01|DSP-', '', Sample_ID_NTC$Sample_ID)
# number=3
# string = substr(string, nchar(string) - number + 1, nchar(string))
# 
# for(n in string)
# {
#   # print(grep(n, names(dccData)))
#   
#   selected = grep(n, names(dccData))
#   
#   print(names(dccData)[which(names(dccData)[selected] %in% Sample_ID_NTC$Sample_ID)])
#   
#   t=1
#   for(i in selected)
#   {
#     # print(dim(dccData[[i]][[4]]))
#     temp = dccData[[i]][[4]]
#     
#     #temp = temp[1:50,]
#     
#     colnames(temp)[1] = gsub(paste0('.*', n), '', names(dccData)[i])
#     colnames(temp)[2] = paste0('Count_', gsub(paste0('.*', n), '', names(dccData)[i]) )
#     
#     if(names(dccData)[i] %in% c(highNTC, Sample_ID_NTC$Sample_ID)) # if its NTC has a problem
#     {
#       #print(names(dccData)[i])
#       #if(t==1){RESULT = temp; t=t+1}else{RESULT = cbind(RESULT, temp)}
#       
#       if(t==1){RESULT = temp; t=t+1}else{
#         
#         temp2 = merge(RESULT, temp, by = 'row.names')
#         
#         row.names(temp2) = temp2$Row.names
#         
#         toMatch <- c('Count_')
#         matches <- unique(grep(paste(toMatch, collapse="|"), colnames(temp2), value=F))
#         temp2 = temp2[,matches] 
#         temp2 = temp2[order(temp2[,1], decreasing = T),]
#         
#         fwrite(temp2, paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_reseq_used_for_analysis/NTC/NTC_', n, '_', colnames(temp2)[2], '.txt'), quote = F, sep = '\t', row.names = T)
#         
#         print(dim(temp2))
#       }
#     }
#   }
# }

#============================================ test

#---- test
# data <- de_results_BvT$adj.P.Val[which(de_results_BvT$adj.P.Val > 0.5)]  # Generate sample data
# length(data)
# hist_result <- hist(data, breaks = 100, plot = F)  # Create histogram, but don't plot
# print(mean(hist_result$counts))
# Calculate average height of bins in the range of 0.5 to 1
# high_p_value_bins <- which(hist_result$breaks >= 0.5)
# background_level <- mean(!is.na(hist_result$counts[high_p_value_bins])) # Exclude the first break
# 
# print(paste("Estimated background level:", background_level))
# abline(h = background_level, col = "red")
#----
# library(ggplot2)
# 
# sum(summary_efit[1] + summary_efit[2] + summary_efit[3])
# 
# p = ggplot(as(de_results_BvT, "data.frame"), aes(x = adj.P.Val)) +
#   geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 1, alpha=0.6) + ggtitle(paste0(names$names2[names$names1 == name], ' (DEGs: ',sum(summary_efit[1] + summary_efit[3]),', Up-regulated: ',summary_efit[3],', Down-regulated: ',summary_efit[1],')')) + geom_hline(aes(yintercept=mean(hist_result$counts)), color="black", linetype="dashed", size=1) +
#   annotate("text", x = 0.5, y = 100, label = paste("Mean =", round(mean(hist_result$counts), 2), " FDR:", round(mean(hist_result$counts)/length(de_results_BvT$adj.P.Val[which(de_results_BvT$adj.P.Val < 0.01)]), 4)), color = "black", size = 10) +
#   theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"), legend.text=element_text(size=20), legend.title=element_text(size=20), plot.title = element_text(size = 19, face="bold"), 
#         panel.border = element_rect(color = "black", fill = NA),
#         panel.grid = element_line(color = "#EEEEEE"),
#         panel.background = element_rect(fill = NA),
#         legend.key = element_rect(fill = NA)) + 
#   xlab("P-value") + ylab("Count") 
#   
# print(p)
# 
# setwd(datadir)
# pdf(paste0(name, "_PVALUES.pdf"), width = 10, height = 7)
# print(p)
# dev.off()


#---------------------------- select samples for transferring to GSK
# setwd('/Users/isarnassiri/Documents/GeoMx_Analysis')
# library(data.table)
# List_excluding_samples_GSK = fread(paste0('List_excluding_samples_GSK.txt'), stringsAsFactors = F)
# List_excluding_samples_GSK$SampleIDs = gsub('-','.',List_excluding_samples_GSK$SampleIDs)
# 
# GeoMx_TMM_RUV4_CellTypeCorrection_readCount = fread(paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/GeoMx_TMM_RUV4_CellTypeCorrection_readCount.txt'), stringsAsFactors = F)
# GeoMx_TMM_RUV4_CellTypeCorrection_readCount = as.data.frame(GeoMx_TMM_RUV4_CellTypeCorrection_readCount)
# 
# GeoMx_readCount = fread(paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4_limma_20k/GeoMx_readCount.txt'), stringsAsFactors = F)
# GeoMx_readCount = as.data.frame(GeoMx_readCount)
# colnames(GeoMx_readCount)[-1] = gsub('.dcc', '', colnames(GeoMx_readCount)[-1])
# 
# GeoMx_TMM_RUV4_CellTypeCorrection_readCount = GeoMx_TMM_RUV4_CellTypeCorrection_readCount[,-which(colnames(GeoMx_TMM_RUV4_CellTypeCorrection_readCount) %in% c(List_excluding_samples_GSK$SampleIDs))]
# GeoMx_readCount = GeoMx_readCount[,-which(colnames(GeoMx_readCount) %in% c(List_excluding_samples_GSK$SampleIDs))]
# 
# fwrite(GeoMx_TMM_RUV4_CellTypeCorrection_readCount, paste0('GeoMx_TMM_RUV4_CellTypeCorrection_readCount.txt'), quote = F, sep = '\t', row.names = T)
# fwrite(GeoMx_readCount, paste0('GeoMx_readCount.txt'), quote = F, sep = '\t', row.names = T)






