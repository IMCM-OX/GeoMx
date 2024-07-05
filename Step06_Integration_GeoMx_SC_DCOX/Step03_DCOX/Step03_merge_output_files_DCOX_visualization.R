#--------------- 
SAVER='/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/DCOA/SAVER/'
input = '/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/DCOA/DCOXA/'
parent = '/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/DCOA/'

setwd(parent)
dir.create('DCOX_Plots/', recursive = T)
#---------------
library(data.table)
listfiles = list.files(input, pattern = '_DCOG.txt')

listfiles_Pair = listfiles

setwd(input)
merged_Pair = rbindlist(lapply(listfiles_Pair, fread))

merged_Pair = as.data.frame(merged_Pair)

merged_Pair = merged_Pair[order(-as.numeric(merged_Pair$TPrate), merged_Pair$pValDiff_adj),]

# View(merged_Pair[which(merged_Pair$Gene2 == 'IFI6'),])
# IFI6_SYT1-_RealCluster_correlationplot.pdf - IFI6 is gene 2

dim(merged_Pair[which(merged_Pair$pValDiff_adj < 1e-8),])
merged_Pair = merged_Pair[which(merged_Pair$pValDiff_adj < 1e-5),]

# fwrite(merged_Pair, paste0('DCOX_Plots/', gsub('.txt','',Pairs[i]),'_merged.txt'), quote = F, row.names = F, sep = '\t')


#-------------------------------- test
temp = merged_Pair#[which(merged_Pair$pValDiff_adj < 1e-6 & merged_Pair$TPrate >= 0.5),] # 
dim(temp)
length(unique(temp$Gene1))+length(unique(temp$Gene2))

Gene1 = 'MT.CO2'
Gene2 = 'HSP90AB1'

Gene1 = 'MT.CO2'
Gene2 = 'KAZN'

Gene1 = 'TUBB2A'
Gene2 = 'HSP90AB1'
View(temp[which(temp$Gene2 == Gene1 & temp$Gene1 == Gene2),])

DEGs = fread('/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/Figures/Annotated_DEGenes.txt', stringsAsFactors = F, header = T)
DEGs = DEGs[which(DEGs$padj < 0.01),]

intersect(DEGs$feature, c(unique(temp$Gene1), unique(temp$Gene2)))[which(intersect(DEGs$feature, c(unique(temp$Gene1), unique(temp$Gene2))) %in% c('TUBB2A', 'HSP90AB1', 'MT.CO2','KAZN' ))]

length(intersect(DEGs$feature, c(unique(temp$Gene1), unique(temp$Gene2))))
length(DEGs$feature)

# annotation
library(dplyr)
library(annotables)
grch38 = tbl_df(grch38)
View(grch38[which(grch38$symbol %in% c('TUBB2A', 'HSP90AB1', 'MT-CO2','KAZN' )),])

View(temp[which(temp$Gene2 %in% c('TUBB2A', 'HSP90AB1', 'MT-CO2','KAZN' ) | temp$Gene1 %in% c('TUBB2A', 'HSP90AB1', 'MT-CO2','KAZN' )),])

#--------------------------------

#========= read inputs

estimate = fread(paste0(SAVER, 'SAVER.txt'), stringsAsFactors = F, header = T)
estimate = as.data.frame(estimate)

row.names(estimate) = make.names(estimate[,1], unique = T)
estimate = estimate[,-c(1)]
estimate[1:5,1:5]
dim(estimate)

#=========  donors ID
MetaData <- fread(paste0(SAVER, 'Metadata_SAVER.txt'), stringsAsFactors = F, header = T)
MetaData = as.data.frame(MetaData)
dim(MetaData)

colnames(estimate) = gsub('\\.','-', colnames(estimate))

S1_BARCODE = colnames(estimate)[which(colnames(estimate) %in% MetaData$Barcode[MetaData$yan == 'Vulnerable'])]
S2_BARCODE = colnames(estimate)[which(colnames(estimate) %in% MetaData$Barcode[MetaData$yan == 'Resistant'])]

length(S1_BARCODE)
length(S2_BARCODE)

print(length(merged_Pair$Gene1))

#========================================== visualization
t=1
for (t in 1:length(merged_Pair$Gene1)) # length(merged_Pair$Gene1)
{

# Gene1 = 'VDAC1' # in RESULT_merged, Gene2 is the original Gene1
# Gene2 = 'PTBP1'
  
Gene1=merged_Pair$Gene2[t]# in RESULT_merged, Gene2 is the original Gene1
Gene2=merged_Pair$Gene1[t]

estimate_Gene1 = estimate[which(row.names(estimate) == Gene1),]
estimate_Gene2 = estimate[which(row.names(estimate) == Gene2),]
#----------
input_mixture_model = t(rbind(as.matrix(estimate_Gene1), as.matrix(estimate_Gene2)))
input_mixture_model = as.data.frame(input_mixture_model, stringsAsFactors = F)

# library(flexmix)
# library(mvtnorm)
# mixed_model = flexmix(cbind(as.numeric(input_mixture_model[,1]),as.numeric(input_mixture_model[,2]))~1,
# k=2, data=input_mixture_model,
# model = FLXMCmvnorm(diag = T), #diag = T to get fix results FLXMCmvpois(), FLXMCmvnorm
# control = list(tolerance = 1e-15, iter.max = 1000))

#--------------------------- assignments comparison
input_mixture_model$real_clusters = colnames(estimate)
input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S1_BARCODE)]=1
input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S2_BARCODE)]=2

# input_mixture_model$predicted_clusters = clusters(mixed_model)
# 
# if(input_mixture_model$predicted_clusters[1] == 2)
# {
#   input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 1)]=0
#   input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 2)]=1
#   input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 0)]=2
# }

#---------------------------- visualization (real clusters)

S1=unique(input_mixture_model$real_clusters)[1]
S2=unique(input_mixture_model$real_clusters)[2]

expr0=data.frame(Gene1=as.numeric(input_mixture_model[which(input_mixture_model[,'real_clusters']==S1),1]),Gene2=as.numeric(input_mixture_model[which(input_mixture_model[,'real_clusters']==S1),2]) )
expr1=data.frame(Gene1=as.numeric(input_mixture_model[which(input_mixture_model[,'real_clusters']==S2),1]),Gene2=as.numeric(input_mixture_model[which(input_mixture_model[,'real_clusters']==S2),2]) )

expr0$Sample = rep(paste0('Sample-0'), dim(expr0)[1])
expr1$Sample = rep(paste0('Sample-1'), dim(expr1)[1])

input.plot = do.call(rbind,list(expr0,expr1))
input.plot$Sample = factor(input.plot$Sample, levels = unique(input.plot$Sample))


#---------- remove outliers
#' Detect outliers using IQR method
#'
#' @param x A numeric vector
#' @param na.rm Whether to exclude NAs when computing quantiles
#'
is_outlier <- function(x, na.rm = FALSE) {
  qs = quantile(x, probs = c(0.05, 0.95), na.rm = na.rm)

  lowerq <- qs[1]
  upperq <- qs[2]
  iqr = upperq - lowerq

  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)

  # Return logical vector
  x > extreme.threshold.upper | x < extreme.threshold.lower
}

#' Remove rows with outliers in given columns
#'
#' Any row with at least 1 outlier will be removed
#'
#' @param df A data.frame
#' @param cols Names of the columns of interest. Defaults to all columns.
#'
#'
remove_outliers <- function(df, cols = names(df)) {
  for (col in cols) {
    cat("Removing outliers in column: ", col, " \n")
    df <- df[!is_outlier(df[[col]]),]
  }
  df
}

vars_of_interest <- c("Gene1", "Gene2")
input.plot <- remove_outliers(input.plot, vars_of_interest)
#----------

library(grDevices)
library(grid)
library(ggpubr)
library(ggplot2)
theme_set(theme_bw())
g=ggplot(input.plot, aes( Gene1,Gene2))+labs(subtitle="Sample-specific Correlation for Expression Values",title="", x = Gene1, y = Gene2)
p=g+geom_jitter(aes(col=Sample,size=Gene1))+geom_smooth(aes(col=Sample),method="lm",se=F)+scale_color_manual(values=c("#999999","#E69F00","#56B4E9"))+geom_rug(aes(color=Sample))+theme(axis.text=element_text(size=30),axis.title=element_text(size=30),plot.subtitle=element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20))

pdf(file=paste0(parent, '/DCOX_Plots/', Gene1, '_', Gene2,'-', "_RealCluster_correlationplot.pdf"), width = 15, height = 10, useDingbats = F)
print(p)
dev.off()
}


#-- Note
# If pair[i] is F1_F2
# Sample-0 is F1 and Sample-1 is F2


