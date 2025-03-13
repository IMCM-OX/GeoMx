
#============================================ Repeat the analysis using new parameters

library(GeomxTools)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(tidyverse)

## Introduction
# This summarises the QC of the GeoMx pilot run of 10th October, utilising the 'GeoMxTools' pipeline.

#============================================ Collate files
datadir <- file.path("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S3_S4/")
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
library(dplyr)
library(ggforce)
count_mat = data.frame(GeneSymbol = Data@featureData@data$TargetName, Data@assayData$exprs)

library(data.table)
fwrite(count_mat, paste0(datadir, 'GeoMx_readCount.txt'), quote = F, row.names = F, sep = '\t')

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
       maxNTCCount = 20000,     # for WGT is 10K. Maximum counts observed in NTC well (1000)
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
fwrite(new_sData_df, paste0(datadir, 'QC_All.txt'), quote = F, row.names = F, sep = '\t')

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
kable(table(pData(target_Data)$DetectionThreshold,
            pData(target_Data)$segment))

#============================================ save QC incluidng gene detection rate
dimnames(pData(target_Data))[[2]]
target_Data_df <- pData(target_Data)[, c("GenesDetected", "GeneDetectionRate", "DetectionThreshold")]

target_Data_df$SampleID = gsub('.dcc', '', row.names(target_Data_df))
QC_Complete = merge(target_Data_df, new_sData_df, by = 'SampleID')
dim(QC_Complete)

getwd()

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

fwrite(QC_Complete, paste0(datadir, 'QC_All_inc_GeneDetection.txt'), quote = F, row.names = F, sep = '\t')

#---
rm(list = c('Data_backup', 'Data', 'QCResults_backup', 'QCResults'))

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

#============================================ Gene Filtering
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

target_Data <- normalize(target_Data, norm_method="quant", desiredQuantile = .75, toElt = "q_norm")
target_Data <- normalize(target_Data, norm_method = "neg", fromElt = "exprs", toElt = "neg_norm")
exprs(target_Data)[1:5, 1:2]

data.frame(assayData(target_Data)[["exprs"]][seq_len(3), seq_len(3)])
featureType(target_Data)
data.frame(assayData(target_Data)[["q_norm"]][seq_len(3), seq_len(3)])

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

#============================================ Normalization

# View(as.data.frame(colData(spe_ruv)))
# colData(spe_ruv)$slide.name
# colData(spe_ruv)$scan.name
# colData(spe_ruv)$segment
# colnames(colData(spe_ruv))
# 
table(colData(spe2_subset)$segment)
# control    PD-with-LB PD-without-LB 
# 235           119           166 

#------------------ I should select the Run ------------------
Run = 'Run 1'
spe2_subset <- spe2_subset[, grep(Run, colData(spe2_subset)$scan.name)]

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

SamplesMetadata = read_excel(paste0(datadir, 'annotation/Annotation_File.xlsx'))
dim(SamplesMetadata)
SamplesMetadata$Run = gsub('.*Run ', '', SamplesMetadata$`scan name`)

SamplesMetadata = SamplesMetadata[!is.na(SamplesMetadata$Run),]
dim(SamplesMetadata)

table(SamplesMetadata$segment)
table(SamplesMetadata$segment[which(SamplesMetadata$segment == "PD-with-LB")])
table(SamplesMetadata$segment[which(SamplesMetadata$Run == "1")])
table(SamplesMetadata$segment[which(SamplesMetadata$Run == "2")])
#--------------------------------------

spe2_subset@assays
spe2_subset@assays@data@listData$logcounts = NULL
spe_tmm <- geomxNorm(spe2_subset, method = "TMM")
# I test to be sure that I save TMM in logcounts
spe_tmm@assays # logcounts [assay = 2] contain the TMM

plotRLExpr(spe_tmm, assay = 2, color = segment) + ggtitle("TMM")
drawPCA(spe_tmm, assay = 2, color = segment)

#============================================ Batch correction
spe_tmm <- findNCGs(spe_tmm, batch_name = "segment", top_n = 300)
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


format(object.size(spe2_subset), units = "Gb")
rm(list = c('spe_tmm', 'spe2_subset'))

#============================================ UMAP  

library(ggspavis)
## Set seed
set.seed(987)

## Compute UMAP on top 50 PCs
spe_ruv <- scater::runPCA(spe_ruv)
spe_ruv <- runUMAP(spe_ruv, dimred = "PCA")
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

UMAP = ggplot(data = as.data.frame(spe_ruv@int_colData@listData$reducedDims$UMAP),
       aes(x = UMAP1, y = UMAP2, colour = spe_ruv@colData$segment)) + 
  geom_point(size = 5) + 
  scale_colour_brewer(type = "qual") + 
  labs(title = "Reduced dimensions: UMAP",
       x = "UMAP1",
       y = "UMAP2",
       colour = "Layers") +
  theme_classic()

getwd()
pdf(paste0(gsub(' ','',Run), "_TMM_RUV_UMAP.pdf"), width = 8, height = 8)
UMAP
dev.off()

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

#============================================ Differential expression analysis ============================================

#--- summary of segments used for DEGs analysis
table(colData(spe_ruv)$segment)

library(edgeR)
library(limma)

dge <- SE2DGEList(spe_ruv)

design <- model.matrix(~0 + segment + ruv_W1 + ruv_W2 , data = colData(spe_ruv))
colnames(design)

colnames(design) <- gsub("^segment","", colnames(design))
gsub("\\+","",colnames(design)) 
colnames(design) <- gsub("\\+","",colnames(design)) # Double check this part
colnames(design) = c("control", "PDwLBP", "PDwoLBP",  "ruv_W1",   "ruv_W2")

#------------------ select a reference ------------------

# LB509TH - TH (TH is reference)
# contr.matrix <- makeContrasts(BvT = PDwLBP - PDwoLBP, levels = colnames(design))

states = c("PDwLBP - control", "PDwLBP - PDwoLBP", "PDwoLBP - control")
state = "PDwLBP - PDwoLBP"

for(state in states)
{

contr.matrix <- makeContrasts(BvT = state, levels = colnames(design))
  
query = gsub(' - .*','',state)
control = gsub('.* - ','',state)

#============================================ Limma normalization

v <- list()
# 1. Normalizing only to library size:
v$libsize <- voom(dge, design)
# 2. TMM normalization:
dge <- calcNormFactors(dge)
v$tmm <- voom(dge, design)
# 3. Cyclic loess normalization:
v$cycloess <- voom(dge, design, normalize="cyclicloess")
# 4. Quantile normalization:
v$quantile <- voom(dge, design, normalize="quantile")

names(spe_ruv@assays)
spe_ruv@assays@data@listData$LimmaTMM = v$tmm$E
spe_ruv@assays@data@listData$Limmacyclicloess = v$cycloess$E
spe_ruv@assays@data@listData$LimmaQuantile = v$quantile$E
spe_ruv@assays@data@listData$libsize = v$libsize$E
# 
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

## ----------------------------------------------------------------------------- PCA plots

names = data.frame(names1 = names(spe_ruv@assays),
                   names2 = c('Counts','TMM + RUV4','Limma TMM','Cyclic Loess','Quantile','Library size'))

for(i in 1:length(names(spe_ruv@assays)) )
{
  name = names(spe_ruv@assays)[i]
  # print(names$names2[names$names1 == name])
  setwd(datadir)

  pdf(paste0(name, "_plotPairPCA.pdf"), width = 10, height = 10)

  library(factoextra)
  library(tidyverse)
  identical(colnames(spe_ruv@assays@data@listData$logcounts), row.names(colData(spe_ruv)))
  res.pca <- prcomp(t(spe_ruv@assays@data@listData[[name]]),  scale = TRUE)
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
                  legend.title = list(fill = "Species", color = "Contrib",
                                      alpha = "Contrib")
  ) + ggtitle(names$names2[names$names1 == name]) + coord_fixed() + theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"), legend.text=element_text(size=20), legend.title=element_text(size=20), plot.title = element_text(size = 30, face="bold"))

  print(p)
  dev.off()
  
  
  #-----------------
  print(fviz_eig(res.pca, addlabels = TRUE))
  # Graph of the variables
  # fviz_pca_var(res.pca, col.var = "black")
  
  fviz_cos2(res.pca, choice = "var", axes = 1:2)
  # fviz_pca_var(res.pca, col.var = "cos2",
  #              gradient.cols = c("black", "orange", "green"),
  #              repel = TRUE)
  
  #--- variable selection using PCA for comparison with 5k
 
  if(name == 'logcounts')
  {
  # attr(spe_ruv@int_colData@listData$reducedDims$PCA, "varExplained")
  # attr(spe_ruv@int_colData@listData$reducedDims$PCA, "rotation")
  
    #print(name)
    
  # Loadings
  loadings <- res.pca$rotation
  
  number_of_PC_to_keep = 2;
  pc_importance_short <- abs(loadings[, 1:number_of_PC_to_keep])
  variable_importance_short <- rowSums(pc_importance_short)
  selected_variables_short <- names(sort(variable_importance_short, decreasing = TRUE))
  length(selected_variables_short)
  
  #Example using percentile.
  threshold_percentile = 0.9;
  threshold_value = quantile(variable_importance_short, threshold_percentile);
  selected_variables_percentile = variable_importance_short[which(variable_importance_short > threshold_value)]
  length(selected_variables_percentile)
  
  getwd()
  write.table(selected_variables_percentile, paste0(name, "_selected_variables_percentile.txt"), sep = '\t', col.names = F)
  }
  #-----------------
  
}

#------------------

# As I use this potentially can omit the GeoMx pipeline gene filtering
keep <- filterByExpr(dge, design) # , min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7
table(keep)

dge_all <- dge[keep, ]
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)

plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))
bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)

## ----------------------------------------------------------------------------- histogram of differentially expressed genes

# highbcv <- bcv_df$BCV > 0.8
# highbcv_df <- bcv_df[highbcv, ]
# points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
# text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)
# v <- voom(dge_all, design, plot = F, normalize="quantile")

names(spe_ruv@assays)
gsub(" ", "','",capture.output(cat(names(spe_ruv@assays))))

names = data.frame(names1 = names(spe_ruv@assays),
names2 = c('counts','TMM+RUV4','Limma TMM','Cyclic Loess','Quantile','Library size'))

getwd()

for(name in names(spe_ruv@assays))
{
print(name)
fit <- lmFit(assay(spe_ruv, name), design = design)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit <- decideTests(efit, p.value = 0.025)
summary_efit <- summary(results_efit)
print(summary_efit)

de_results_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
dim(de_results_BvT)

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

}

## ----------------------------------------------------------------------------- DEG analysis

fit <- lmFit(assay(spe_ruv, "logcounts"), design = design)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit <- decideTests(efit, p.value = 0.025)
summary_efit <- summary(results_efit)
print(summary_efit)

library(ggrepel)
library(tidyverse)
de_results_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
dim(de_results_BvT)
dim(assay(spe_ruv, "logcounts"))
de_results_BvT$TargetName = row.names(de_results_BvT)

de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)
dim(de_genes_toptable_BvT)
# options(ggrepel.max.overlaps = 10)

#============================================ calculating sample size 

#- effect size (Log fold change (LFC) estimates)
# de_genes_toptable_BvT_001 <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.01)
# dim(de_genes_toptable_BvT_001)
# summary(de_genes_toptable_BvT_001$logFC) # The fact that limma gives log2 results is mentioned many times in the documentation.
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # -1.7048 -0.8966 -0.6015 -0.1493  0.9815  2.0953 
# # A log2 fold change of 1 means the gene is expressed twice as high in the test condition.
# # A log2 fold change of 2 means the gene is expressed four times higher in the test condition.
# # A log2 fold change of 3 means the gene is expressed eight times higher in the test condition, and so on.
# 
# # Cohen suggests that d values of 0.2, 0.5, and 0.8 represent small, medium, and large effect sizes respectively.
# # I select 0.5
# 
# results <- pwr::pwr.t.test(n = NULL,
#                            sig.level = 0.05, 
#                            type = "two.sample", 
#                            alternative = "two.sided", 
#                            power = 0.80, 
#                            d = 0.5)
# 
# plot(results) +
#   ggplot2::theme_minimal(base_size = 14) +
#   labs(title = 'Optimizing Sample Size for 2-Sided t test',
#        subtitle = "")

## ----------------------------------------------------------------------------- visualization of DEGs
summary(de_results_BvT$logFC)
table(de_results_BvT$adj.P.Val < 0.05 & de_results_BvT$logFC > 0) # using logFC > 0.1 gives me same results
table(de_results_BvT$adj.P.Val < 0.05 & de_results_BvT$logFC < 0)
dim(de_results_BvT)

# de_results_BvT %>% 
#   mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
#                      ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
#   ggplot(aes(AveExpr, logFC, col = DE)) + 
#   geom_point(shape = 1, size = 1) + 
#   geom_text_repel(data = de_genes_toptable_BvT %>% 
#                     mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
#                                        ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
#                     rownames_to_column(), aes(label = rowname)) +
#   theme_bw() +
#   xlab("Average log-expression") +
#   ylab("Log-fold-change") +
#   ggtitle("LB509+TH+ vs LB509-TH+ as reference (limma-voom)") +
#   scale_color_manual(values = c("dodgerblue","lightblue","orange2")) +
#   theme(text = element_text(size=15))
# 
# dev.off()

library(DT)
updn_cols <- c(RColorBrewer::brewer.pal(6, 'Greens')[2], RColorBrewer::brewer.pal(6, 'Purples')[2])

library(data.table)
fwrite(de_results_BvT, paste0(datadir, Run, "_", query, "_", control, '_GeoMx_DEA_20KNTC.txt'), quote = F, row.names = F, sep = '\t')

#======
# temp = fread(paste0(datadir, 'GeoMx_DEA.txt'), stringsAsFactors = F, header = T)
# 
# dim(temp)
# dim(de_genes_toptable_BvT)
# length(which(de_genes_toptable_BvT$TargetName %in% temp$TargetName))

#======
de_results_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
class(de_results_BvT)

library(ggrepel)
# Categorize de_results_BvT based on P-value & FDR for plotting
de_results_BvT$Color <- "NS or FC < 0.5"
de_results_BvT$Color[which(de_results_BvT$adj.P.Val < 0.05)] <- "FDR < 0.05"
de_results_BvT$Color[which(de_results_BvT$adj.P.Val < 0.01 & de_results_BvT$logFC > 0)] <- "FDR < 0.01 & up-regulated"
de_results_BvT$Color[which(de_results_BvT$adj.P.Val < 0.01 & de_results_BvT$logFC < 0)] <- "FDR < 0.01 & down-regulated"
de_results_BvT$Color[which(abs(de_results_BvT$logFC) < 0.1)] <- "NS or FC < 0.1"
de_results_BvT$Color <- factor(de_results_BvT$Color,
                               levels = c("NS or FC < 0.1", "FDR < 0.05", "FDR < 0.01 & up-regulated", "FDR < 0.01 & down-regulated"))
table(de_results_BvT$Color)

# pick top TargetNames for either side of volcano to label
# order TargetNames for convenience:
de_results_BvT$invert_P <- (-log10(de_results_BvT$P.Value)) * sign(de_results_BvT$logFC)
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
       aes(x = logFC, y = -log10(P.Value),
           color = Color, label = TargetName)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = expression(paste('Log'['10'],'P')),
       y = expression(paste('Log'['2'],' fold change')),
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.01 & down-regulated` = "dodgerblue",
                                `FDR < 0.01 & up-regulated` = "darkred",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.01` = "orange2",
                                `NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(de_results_BvT, TargetName %in% top_g & adj.P.Val < 0.001),
                  size = 5, point.padding = 0.2, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom")

getwd()
pdf(paste0(gsub(' ','',Run), "_", query, "_", control, "_DEG_Plot.pdf"), width = 10.5, height = 8)
print(DEG_Plot)
dev.off()


## ----------------------------------------------------------------------------- visualization of DEGs as table
# de_genes_toptable_BvT %>% 
#   dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
#   DT::datatable(caption = '(limma-voom)') %>%
#   DT::formatStyle('logFC',
#                   valueColumns = 'logFC',
#                   backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
#   DT::formatSignif(1:4, digits = 4)

#============================================ Within Slide Analysis - DEGs
# https://bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

# table(pData(target_Data)$segment)
# target_Data_subset = target_Data[,which(pData(target_Data)$segment != "control")]
# 
# # convert test variables to factors
# pData(target_Data_subset)$testRegion <- factor(pData(target_Data_subset)$segment, c("control", "PD-with-LB", "PD-without-LB"))
# pData(target_Data_subset)[["slide"]] <-  factor(pData(target_Data_subset)[["slide name"]])
# assayDataElement(object = target_Data_subset, elt = "log_q") <- assayDataApply(target_Data_subset, 2, FUN = log, base = 2, elt = "q_norm")
# 
# # run LMM:
# # formula follows conventions defined by the lme4 package
# results <- c()
# 
# mixedOutmc <-
#   mixedModelDE(target_Data_subset,
#                elt = "log_q",
#                modelFormula = ~ testRegion + (1 + testRegion | slide),
#                groupVar = "testRegion",
#                nCores = parallel::detectCores(),
#                multiCore = TRUE)
# 
# # format results as data.frame
# r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
# tests <- rownames(r_test)
# r_test <- as.data.frame(r_test)
# r_test$Contrast <- tests
# 
# # use lapply in case you have multiple levels of your test factor to
# # correctly associate gene name with it's row in the results table
# r_test$Gene <-
#   unlist(lapply(colnames(mixedOutmc),
#                 rep, nrow(mixedOutmc["lsmeans", ][[1]])))
# 
# r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
# r_test <- r_test[, c("Gene", "Contrast", "Estimate",
#                      "Pr(>|t|)", "FDR")]
# results <- rbind(results, r_test)
# 
# fwrite(results, paste0(datadir, 'Within_Slide_Analysis_DEGs_wol_vs_wl.txt'), quote = F, row.names = F, sep = '\t')


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
# GeneList = assay(spe_ruv, "logcounts")
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

input = de_results_BvT[de_results_BvT$adj.P.Val < 0.01,]
input = input[-which(is.na(as.character(mapIds(org.Hs.eg.db, input$TargetName, 'ENTREZID', 'SYMBOL')))),]
row.names(input) = as.character(mapIds(org.Hs.eg.db, input$TargetName, 'ENTREZID', 'SYMBOL'))

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

list_pathways = c('05010', '05012', '05014', '05016', '05017', '05020', '05022')

for(p in 1:length(list_pathways))
{
  pv.out <- pathview(gene.data = input_used, pathway.id = list_pathways[p],
                     species = "hsa", kegg.native = T,
                     same.layer = F, high=list(gene="gold"), out.suffix = paste0(gsub(' ','',Run), "_", query, "_", control, '_VR'))
}

## ----------------------------------------------------------------------------- pathway Enrichment v2

#--- gene enrichment analysis
library('org.Hs.eg.db')
BD = 'c5.all.v2023.1.Hs.entrez.gmt'

library(qusage)
library(clusterProfiler)
gmtfile <- paste0("/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Scripts_USED_visium/msigdb_v2023_Hs/", BD)
c5 <- read.gmt(gmtfile)

# use mapIds method to obtain Entrez IDs
input = de_results_BvT[de_results_BvT$adj.P.Val < 0.01,]
de = as.character(mapIds(org.Hs.eg.db, input$TargetName, 'ENTREZID', 'SYMBOL'))

table(c5$term)
c5$term = sub("^[^_]*_", "", c5$term)

egmt <- enricher(de[!is.na(de)], TERM2GENE=c5)

if(length(unique(egmt@result$p.adjust < 0.05))>1)
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
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

library(enrichplot)
# barplot(edo, showCategory=10)
# dotplot(edo, showCategory=15) + ggtitle(print(paste0(  S )))
# cnetplot(edox, circular = TRUE, colorEdge = TRUE)

if(length(unique(edo@result$p.adjust < 0.05))>1)
{
  pdf(paste0(gsub(' ','',Run), "_", query, "_", control, '_DiseaseEnrichment.pdf'),width = 8, height = 9, useDingbats = FALSE)
  print(dotplot(edo, showCategory=15))
  dev.off()

  fwrite(edo@result, paste0(gsub(' ','',Run), "_", query, "_", control, '_DiseaseEnrichment.txt') , quote = F, sep = '\t', row.names = F)
}
}
