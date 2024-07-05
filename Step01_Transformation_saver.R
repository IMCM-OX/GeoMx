#---------------
output = '/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/GeoMx/'

setwd(output)

#---------------
library(data.table)
library(SAVER)
library(Seurat)
#========================================== transform assigned cells per selected pairs of individuals

#========= read gene expression profile
read_count_S1S2 <- fread('GeoMx_readCount.txt', stringsAsFactors = F, header = T)
read_count_S1S2 = as.data.frame(read_count_S1S2)
dim(read_count_S1S2)

rownames(read_count_S1S2) = make.names(read_count_S1S2$GeneSymbol,unique=T)
read_count_S1S2 = read_count_S1S2[,-which(colnames(read_count_S1S2) == 'GeneSymbol')]
colnames(read_count_S1S2)[1:10]

#========= subset gene expression profile

read_count_S1S2 = read_count_S1S2[apply(read_count_S1S2, 1, function(x) !all(x==0)),] # I did this in the QC step
dim(read_count_S1S2)

MetaData <- fread('MetaData.txt', stringsAsFactors = F, header = T)
MetaData = as.data.frame(MetaData)
dim(MetaData)

colnames(read_count_S1S2) = gsub('\\.','-', gsub('.dcc', '', colnames(read_count_S1S2)))

read_count_S1S2 = read_count_S1S2[, which(colnames(read_count_S1S2) %in% MetaData$SampleID)] # I did this in the QC step
dim(read_count_S1S2)

#========= run saver for assigned cells per selected pairs of individuals
set.seed(123)# if you do not use seed, you get different values per run
estimate_S1S2 <- saver(read_count_S1S2, ncores = 6, estimates.only = TRUE)
estimate_S1S2 = data.frame(GeneID=row.names(estimate_S1S2), estimate_S1S2)
estimate_S1S2[1:5,1:5]
#---------------

dir.create(paste0(output, '/SAVER/'))
setwd(paste0(output, '/SAVER/'))

fwrite(estimate_S1S2, 'SAVER.txt', row.names = F, quote = F, sep = '\t')
