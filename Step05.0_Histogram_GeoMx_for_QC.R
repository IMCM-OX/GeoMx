
#-------------- plotting a cumulative frequency curve on a histogram in R

palette.colors(palette = "Okabe-Ito")
#------- first run
library(data.table)
setwd("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output/")

QC_Results = fread('QC_Results.txt', stringsAsFactors = F)

#---
x = (- QC_Results$DeduplicatedReads + QC_Results$Aligned)/QC_Results$Aligned
h <- hist(x, plot = FALSE)

pdf(paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Figures/QC_Results_S1.pdf'), width = 10, height = 10)
plot(h, col = "#0072B2", main='Proportion of Aligned Reads after Deduplication S1', ylim = c(0, sum(h$counts)))
lines(h$mids, cumsum(h$counts))
axis(
  side = 4,
  at = cumsum(h$counts),
  labels = round(cumsum(h$density), digits = 2)
)
dev.off()
#---

S1_raw = fread('GeoMx_readCount.txt', stringsAsFactors = F)
S1_raw = as.data.frame(S1_raw)
rownames(S1_raw) = make.names(S1_raw$GeneSymbol,unique=T)
S1_raw = S1_raw[,-1]
colnames(S1_raw) = gsub('\\.','-',gsub('.dcc','',colnames(S1_raw)))

MetaData = fread('MetaData.txt', stringsAsFactors = F)
MetaData = as.data.frame(MetaData)

S1_raw_Control = S1_raw[,which(colnames(S1_raw) %in% MetaData$SampleID[which(MetaData$segment == "LB509+TH+")])]
S1_raw_Case = S1_raw[,-which(colnames(S1_raw) %in% MetaData$SampleID[which(MetaData$segment == "LB509+TH+")])]

S1_raw_Control$RowSum = rowSums(S1_raw_Control)
S1_raw_Case$RowSum = rowSums(S1_raw_Case)

#-----------------------------
S1_N = fread('GeoMx_readCount_TMMRUV4.txt', stringsAsFactors = F)
S1_N = as.data.frame(S1_N)
rownames(S1_N) = make.names(S1_N$GeneSymbol,unique=T)
S1_N = S1_N[,-1]
colnames(S1_N) = gsub('.dcc','',colnames(S1_N))
colnames(S1_N) = gsub('\\.','-',gsub('.dcc','',colnames(S1_N)))

S1_N_Control = S1_N[,-which(colnames(S1_N) %in% MetaData$SampleID[which(MetaData$segment == "LB509+TH+")])]
S1_N_Case = S1_N[,which(colnames(S1_N) %in% MetaData$SampleID[which(MetaData$segment == "LB509+TH+")])]

S1_N_Control$RowSum = rowSums(S1_N_Control)
S1_N_Case$RowSum = rowSums(S1_N_Case)

cases = c('S1_N_Case', 'S1_N_Control', 'S1_raw_Case', 'S1_raw_Control')

for(c in 1:length(cases))
{

  assign('temp', cases[c])
  temp = get(temp)
  
  x <- temp$RowSum[which(temp$RowSum < summary(temp$RowSum)[[5]])]
  
  pdf(paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Figures/', cases[c],'.pdf'), width = 10, height = 10)
  h <- hist(x, plot = FALSE)
  plot(h, col = "#0072B2", main=cases[c], ylim = c(0, sum(h$counts)))
  lines(h$mids, cumsum(h$counts))
  axis(
    side = 4,
    at = cumsum(h$counts),
    labels = round(cumsum(h$density), digits = 2)
  )
dev.off()
}
getwd()

#------- second run -------
library(data.table)
setwd("/Users/isarnassiri/Documents/GeoMx_Analysis/Machine_output_S2/")

QC_Results = fread('QC_Results.txt', stringsAsFactors = F)

#---
x = (- QC_Results$DeduplicatedReads + QC_Results$Aligned)/QC_Results$Aligned
h <- hist(x, plot = FALSE)

pdf(paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Figures/QC_Results_S2.pdf'), width = 10, height = 10)
plot(h, col = "#D55E00", main='Proportion of Aligned Reads after Deduplication S2', ylim = c(0, sum(h$counts)))
lines(h$mids, cumsum(h$counts))
axis(
  side = 4,
  at = cumsum(h$counts),
  labels = round(cumsum(h$density), digits = 2)
)
dev.off()
#---

S2_raw = fread('GeoMx_readCount.txt', stringsAsFactors = F)
S2_raw = as.data.frame(S2_raw)
rownames(S2_raw) = make.names(S2_raw$GeneSymbol,unique=T)
S2_raw = S2_raw[,-1]
colnames(S2_raw) = gsub('\\.','-',gsub('.dcc','',colnames(S2_raw)))

MetaData = fread('MetaData.txt', stringsAsFactors = F)
MetaData = as.data.frame(MetaData)

S2_raw_Control = S2_raw[,-which(colnames(S2_raw) %in% MetaData$SampleID[which(MetaData$segment == "LB509+")])]
S2_raw_Case = S2_raw[,which(colnames(S2_raw) %in% MetaData$SampleID[which(MetaData$segment == "LB509+")])]

S2_raw_Control$RowSum = rowSums(S2_raw_Control)
S2_raw_Case$RowSum = rowSums(S2_raw_Case)

#-----------------------------
S2_N = fread('GeoMx_readCount_TMMRUV4.txt', stringsAsFactors = F)
S2_N = as.data.frame(S2_N)
rownames(S2_N) = make.names(S2_N$GeneSymbol,unique=T)
S2_N = S2_N[,-1]

colnames(S2_N) = gsub('.dcc','',colnames(S2_N))
colnames(S2_N) = gsub('\\.','-',gsub('.dcc','',colnames(S2_N)))

S2_N_Control = S2_N[,-which(colnames(S2_N) %in% MetaData$SampleID[which(MetaData$segment == "LB509+")])]
S2_N_Case = S2_N[,which(colnames(S2_N) %in% MetaData$SampleID[which(MetaData$segment == "LB509+")])]

S2_N_Control$RowSum = rowSums(S2_N_Control)
S2_N_Case$RowSum = rowSums(S2_N_Case)

cases = c('S2_N_Case', 'S2_N_Control', 'S2_raw_Case', 'S2_raw_Control')

for(c in 1:length(cases))
{
  
  assign('temp', cases[c])
  temp = get(temp)
  
  x <- temp$RowSum[which(temp$RowSum < summary(temp$RowSum)[[5]])]
  
  pdf(paste0('/Users/isarnassiri/Documents/GeoMx_Analysis/Figures/', cases[c],'.pdf'), width = 10, height = 10)
  h <- hist(x, plot = FALSE)
  plot(h, col = "#D55E00", main=cases[c], ylim = c(0, sum(h$counts)))
  lines(h$mids, cumsum(h$counts))
  axis(
    side = 4,
    at = cumsum(h$counts),
    labels = round(cumsum(h$density), digits = 2)
  )
  dev.off()
}
getwd()




