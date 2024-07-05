#--------------- 
SAVER='/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/GeoMx/SAVER/'

setwd('/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/GeoMx/DCOXA/')
dir.create('DCOX_Plots/', recursive = T)
#---------------
library(data.table)
listfiles = list.files(pattern = '_DCOG.txt')
Pairs = unique(gsub('-.*', '', listfiles))

length(Pairs)
i=1
for (i in 1:length(Pairs))
{

Pairs[i]

listfiles_Pair = list.files(pattern = Pairs[i])
merged_Pair = rbindlist(lapply(listfiles_Pair, fread))

merged_Pair = as.data.frame(merged_Pair)

merged_Pair = merged_Pair[order(-as.numeric(merged_Pair$TPrate), merged_Pair$pValDiff_adj),]
print(i)

merged_Pair = merged_Pair[order(-as.numeric(merged_Pair$TPrate), merged_Pair$pValDiff_adj),]

# fwrite(merged_Pair, paste0('DCOX_Plots/', gsub('.txt','',Pairs[i]),'_merged.txt'), quote = F, row.names = F, sep = '\t')

#========= read inputs

estimate = fread(paste0(SAVER, 'SAVER.txt'), stringsAsFactors = F, header = T)

estimate = as.data.frame(estimate)

row.names(estimate) = make.names(estimate[,1], unique = T)
estimate = estimate[,-c(1)]
estimate[1:5,1:5]
dim(estimate)

#=========  donors ID

MetaData <- fread('/well/singlecell/projects_Isar/scRNAseq_Parkinson/RAW_DATA/GeoMx/MetaData.txt', stringsAsFactors = F, header = T)
MetaData = as.data.frame(MetaData)
dim(MetaData)

colnames(estimate) = gsub('\\.','-', gsub('.dcc', '', colnames(estimate)))

S1_BARCODE = colnames(estimate)[which(colnames(estimate) %in% MetaData$SampleID[MetaData$diseased == 'diseased'])]
S2_BARCODE = colnames(estimate)[which(colnames(estimate) %in% MetaData$SampleID[MetaData$diseased == 'healthy'])]

length(S1_BARCODE)
length(S2_BARCODE)

print(length(merged_Pair$Gene1))

#========================================== visualization
t=1
for (t in 1:length(merged_Pair$Gene1)) # length(merged_Pair$Gene1)
{

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

library(grDevices)
library(grid)
library(ggpubr)
library(ggplot2)
theme_set(theme_bw())
g=ggplot(input.plot, aes( Gene1,Gene2))+labs(subtitle="Sample-specific Correlation for Expression Values",title="", x = Gene1, y = Gene2)
p=g+geom_jitter(aes(col=Sample,size=Gene1))+geom_smooth(aes(col=Sample),method="lm",se=F)+scale_color_manual(values=c("#999999","#E69F00","#56B4E9"))+geom_rug(aes(color=Sample))+theme(axis.text=element_text(size=30),axis.title=element_text(size=30),plot.subtitle=element_text(size=30),legend.title=element_text(size=30),legend.text=element_text(size=20))

pdf(file=paste0('DCOX_Plots/', Gene1, '_', Gene2,'-', gsub('.txt','',Pairs[i]), "_RealCluster_correlationplot.pdf"), width = 15, height = 10, useDingbats = F)
print(p)
dev.off()
}
}

#-- Note
# If pair[i] is F1_F2
# Sample-0 is F1 and Sample-1 is F2


