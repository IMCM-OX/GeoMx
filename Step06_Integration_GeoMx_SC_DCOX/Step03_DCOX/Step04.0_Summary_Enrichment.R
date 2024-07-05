
setwd('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/outputs_spaceranger/unzip/')
FC_all = gsub('./', '', list.dirs(path = ".", full.names = TRUE, recursive = FALSE))
FC_all
length(FC_all)

#---------------
library(data.table)
library(Seurat)
#========================================== barplots

for(i in 1:length(FC_all))
{
  print(i)
#========= read gene expression profile
S1 <- Read10X(paste(FC_all[i],'/outs/filtered_feature_bc_matrix', sep=''))
S1 = as.data.frame(S1)
# dim(S1)

S1 = S1[apply(S1, 1, function(x) !all(x==0)),]

S1_T = S1
S1_T[S1_T > 0] <- 1 

# summary(as.numeric(colSums(S1)))
# summary(as.numeric(rowSums(S1)))

if(i==1){RESULTS = data.frame(Sample = FC_all[i], Number_Genes = dim(S1)[1], Number_Spots = dim(S1)[2], Mean_Gene_per_Spot = summary(colSums(S1_T))[[4]], Max_Gene_per_Spot = summary(colSums(S1_T))[[6]] )}else{RESULTS = rbind(RESULTS, data.frame(Sample = FC_all[i], Number_Genes = dim(S1)[1], Number_Spots = dim(S1)[2], Mean_Gene_per_Spot = summary(colSums(S1_T))[[4]], Max_Gene_per_Spot = summary(colSums(S1_T))[[6]] ))}
}

#=========
library(ggplot2)

RESULT = within(RESULTS, Name <- factor(RESULTS$Sample, labels = unique(RESULTS$Sample)))

Samples_Details = fread('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Samples_Details.txt', stringsAsFactors = F, header = T)
Samples_Details = as.data.frame(Samples_Details)
colnames(Samples_Details)

RESULT = merge(RESULT, Samples_Details, by = 'Sample')
RESULTS$Name = paste(RESULT$Donor_ID, RESULT$Region, sep = '-')
table(RESULTS$Name)

RESULTS = RESULTS[order(RESULTS$Name),]

pdf(file = "/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Number_Genes.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches

p = barplot(height=RESULTS$Number_Genes, names=RESULTS$Name, xlab="", ylab = 'Number of Genes', las=2, mgp=c(4,0.5,0), cex.lab=1.2)
text(p, RESULTS$Number_Genes-5500, paste(RESULTS$Number_Genes, sep=""), xpd = T, srt=90,  adj = c(0,0)) 

dev.off()

pdf(file = "/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Number_Spots.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches

p = barplot(height=RESULTS$Number_Spots, names=RESULTS$Name, xlab="", ylab = 'Number of Spots', las=2, mgp=c(4,0.5,0), cex.lab=1.2, col="#69b3a2")
text(p, RESULTS$Number_Spots-1500, paste(RESULTS$Number_Spots, sep=""), xpd = T, srt=90,  adj = c(0,0)) 

dev.off()

RESULTS = RESULTS[order(RESULTS$Max_Gene_per_Spot, decreasing = T),]

pdf(file = "/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Max_Gene_per_Spot.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches

p = barplot(height=RESULTS$Max_Gene_per_Spot, names=RESULTS$Name, xlab="", ylab = 'Max Gene per Spot', las=2, mgp=c(4,0.5,0), cex.lab=1.2, col="#6989b3")
text(p, RESULTS$Max_Gene_per_Spot-50, paste(RESULTS$Max_Gene_per_Spot, sep=""), xpd = T,  adj = c(0.5,-0.5)) 

dev.off()

RESULTS = RESULTS[order(RESULTS$Mean_Gene_per_Spot, decreasing = T),]
RESULTS$Mean_Gene_per_Spot = round(RESULTS$Mean_Gene_per_Spot, digits = 0)
pdf(file = "/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Mean_Gene_per_Spot.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches

p = barplot(height=RESULTS$Mean_Gene_per_Spot, names=RESULTS$Name, xlab="", ylab = 'Mean Gene per Spot', las=2, mgp=c(4,0.5,0), cex.lab=1.2, col="#b38c69")
text(p, RESULTS$Mean_Gene_per_Spot-20, paste(RESULTS$Mean_Gene_per_Spot, sep=""), xpd = T,  adj = c(0.5,-0.5)) 

dev.off()

#========================================== summary of DCOX gens

setwd('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/DCOX_Plots/')
DCOX_all = gsub('./', '', list.files(path = ".", full.names = TRUE, recursive = FALSE, pattern = '.txt'))
length(DCOX_all)
DCOX_all

require(data.table) ## 1.9.2 or 1.9.3
DCOX_all_merged = rbindlist(lapply(DCOX_all, fread), use.names=TRUE, fill=TRUE, idcol="ID")

DCOX_all_merged$ID = gsub('_merged.txt', '', DCOX_all[DCOX_all_merged$ID])
DCOX_all_merged$GenePair = paste(DCOX_all_merged$Gene2, DCOX_all_merged$Gene1, sep = '_')

length(unique(DCOX_all_merged$ID))

table(DCOX_all_merged$ID)
length(unique(c(DCOX_all_merged$Gene1, DCOX_all_merged$Gene2)))

library(stringr)
DCOX_all_merged$Sample1 = str_split_fixed(DCOX_all_merged$ID, "_", 2 )[,1]
DCOX_all_merged$Sample2 = str_split_fixed(DCOX_all_merged$ID, "_", 2 )[,2]

Sample_Type = fread('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Sample_Type.txt', stringsAsFactors = F, header = F)
Sample_Type = as.data.frame(Sample_Type)
colnames(Sample_Type) = c('Sample1', 'State_Sample1') 

DCOX_all_merged = merge(DCOX_all_merged, Sample_Type, by = 'Sample1')
colnames(Sample_Type) = c('Sample2', 'State_Sample2') 
DCOX_all_merged = merge(DCOX_all_merged, Sample_Type, by = 'Sample2')

DCOX_all_merged$State_Pair = paste(DCOX_all_merged$State_Sample1, DCOX_all_merged$State_Sample2, sep = '_')
DCOX_all_merged$State_Pair[which(DCOX_all_merged$State_Pair == 'PD_Control')] = 'Control_PD'

DCOX_all_merged = as.data.frame(DCOX_all_merged)

# to enrich parkinson disease
DCOX_all_merged_subset = DCOX_all_merged[which(DCOX_all_merged$TPrate >= 0.99 & DCOX_all_merged$pValDiff_adj < 1e-10),]

#- remove extra col
table(c(DCOX_all_merged_subset$Sample1, DCOX_all_merged_subset$Sample2))
dim(table(c(DCOX_all_merged_subset$Sample1, DCOX_all_merged_subset$Sample2)))
setdiff(names(table(c(DCOX_all_merged$Sample1, DCOX_all_merged$Sample2))), names(table(c(DCOX_all_merged_subset$Sample1, DCOX_all_merged_subset$Sample2))))
DCOX_all_merged_subset = DCOX_all_merged_subset[,-which(colnames(DCOX_all_merged_subset) %in% c('Sample2', 'Sample1'))]

# to cover CAB39_HSPA5
# DCOX_all_merged_subset = DCOX_all_merged[which(DCOX_all_merged$TPrate >= 0.97 & DCOX_all_merged$pValDiff_adj < 1e-8),]
# test = DCOX_all_merged_subset[which(DCOX_all_merged_subset$GenePair == 'CAB39_HSPA5'),]

table(DCOX_all_merged$State_Pair)

length(unique(unlist(c(DCOX_all_merged[which(DCOX_all_merged$State_Pair == 'Control_Control'), c('Gene1', 'Gene2')]))))
length(unique(unlist(c(DCOX_all_merged[which(DCOX_all_merged$State_Pair == 'Control_PD'), c('Gene1', 'Gene2')]))))
length(unique(unlist(c(DCOX_all_merged[which(DCOX_all_merged$State_Pair == 'PD_PD'), c('Gene1', 'Gene2')]))))

table(DCOX_all_merged_subset$State_Pair)

gsub(" ", "','",capture.output(cat(c(length(unique(unlist(c(DCOX_all_merged_subset[which(DCOX_all_merged_subset$State_Pair == 'Control_Control'), c('Gene1', 'Gene2')])))),
                                     length(unique(unlist(c(DCOX_all_merged_subset[which(DCOX_all_merged_subset$State_Pair == 'PD_PD'), c('Gene1', 'Gene2')])))),
                                     length(unique(unlist(c(DCOX_all_merged_subset[which(DCOX_all_merged_subset$State_Pair == 'Control_PD'), c('Gene1', 'Gene2')]))))))))

df = data.frame(name = c('Control', 'PD', 'Control_PD' ), value = c(1929, 1656, 2972))

p = ggplot(df, aes(x=name, y=value)) +
  geom_segment( aes(xend=name, yend=0)) +
  geom_point( size=4, color="orange") +
  theme_bw()+ xlab("")  +
  geom_label(aes(name, value , label = signif(value)), 
             colour = "darkred", nudge_x = 0.35, size = 4)


fwrite(DCOX_all_merged, '/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/DCOX.txt', quote = F, row.names = F)

DCOX_all_merged[which(DCOX_all_merged$GenePair == 'CAB39_HSPA5' & DCOX_all_merged$ID == 'E2_E3'),]
DCOX_all_merged[which(DCOX_all_merged$GenePair == 'FBXO28_UCHL1' & DCOX_all_merged$ID == 'C3_F3'),]
DCOX_all_merged[which(DCOX_all_merged$GenePair == 'ARPC4_UCHL1' & DCOX_all_merged$ID == 'C1_F1'),]
DCOX_all_merged[which(DCOX_all_merged$GenePair == 'GPR161_OLFM1' & DCOX_all_merged$ID == 'C4_F1'),]
DCOX_all_merged[which(DCOX_all_merged$ID == 'C1_E3'),]
# UCHL1_ACTB-C3_E2_RealCluster_correlationplot

 

#========================================== Enrichment
#input = unique(unlist(c(DCOX_all_merged_subset[which(DCOX_all_merged_subset$State_Pair == 'Control_PD'), c('Gene1', 'Gene2')])))
input = unique(unlist(c(DCOX_all_merged_subset[, c('Gene1', 'Gene2')])))

#--- gene enrichment analysis

library('org.Hs.eg.db')

# use mapIds method to obtain Entrez IDs
de = as.character(mapIds(org.Hs.eg.db, input, 'ENTREZID', 'SYMBOL'))

# library("DOSE")
# edo <- enrichDGN(de[!is.na(de)])
# edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
# 
# library(enrichplot)
# barplot(edo, showCategory=20)
# dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
#cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

BD = 'c5.all.v2023.1.Hs.entrez.gmt'

library(qusage)
library(clusterProfiler)
gmtfile <- paste0("/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Scripts_USED_visium/msigdb_v2023_Hs/", BD)
c5 <- read.gmt(gmtfile)
egmt <- enricher(de[!is.na(de)], TERM2GENE=c5)

library(cowplot)
library(ggplot2)
p2 <- dotplot(egmt, showCategory=15) + ggtitle("")
p2

pdf(paste0("/Users/isar/Downloads/RESULTS_", BD, ".pdf") ,width = 10, height = 10,useDingbats = FALSE)
print(plot_grid(p2, ncol=1))
dev.off()

#========================================== visualization


#---------------
library(data.table)
 
#========= read inputs
Pairs = 'D3_E4'
Gene2 = 'CTNNA2'
merged_Pair = DCOX_all_merged[which(DCOX_all_merged$ID == Pairs),]

setwd('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen')
estimate = fread(paste0('SAVER/', Pairs, '_SAVER.txt'), stringsAsFactors = F, header = T)
estimate = as.data.frame(estimate)

row.names(estimate) = make.names(estimate$GeneID, unique = T)
estimate = estimate[,-c(1)]
estimate[1:5,1:5]
dim(estimate)

BARCODES = strsplit(Pairs, "\\_")

S1_BARCODE = colnames(estimate)[grep(BARCODES[[1]][1], colnames(estimate))]
S2_BARCODE = colnames(estimate)[grep(BARCODES[[1]][2], colnames(estimate))]

length(S1_BARCODE)
length(S2_BARCODE)

candidate_genes = merged_Pair[which(merged_Pair$Gene2 == Gene2),]

t=1
for (t in 1:dim(candidate_genes)[1]) # length(merged_Pair$Gene1)
{  
  Gene1=candidate_genes$Gene2[t]# in RESULT_merged, Gene2 is the original Gene1
  Gene2=candidate_genes$Gene1[t]
  
  estimate_Gene1 = estimate[which(row.names(estimate) == Gene1),]
  estimate_Gene2 = estimate[which(row.names(estimate) == Gene2),]
  #----------
  input_mixture_model = t(rbind(as.matrix(estimate_Gene1), as.matrix(estimate_Gene2)))
  input_mixture_model = as.data.frame(input_mixture_model, stringsAsFactors = F)
  
  #--------------------------- assignments comparison
  input_mixture_model$real_clusters = colnames(estimate)
  input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S1_BARCODE)]=1
  input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S2_BARCODE)]=2
  
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
  
  pdf(file=paste0('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/test/', Gene1, '_', Gene2,'-', Pairs, "_RealCluster_correlationplot.pdf"), width = 15, height = 10, useDingbats = F)
  print(p)
  dev.off()
}

#========================================== pathview - visualization of hsa00190

INPUT_pathview = as.vector(unique(unlist(c(DCOX_all_merged_subset[which(DCOX_all_merged_subset$State_Pair == 'Control_PD'), c('Gene1', 'Gene2')]))))

library('org.Hs.eg.db')
keytypes(org.Hs.eg.db)
names(INPUT_pathview) = as.character(mapIds(org.Hs.eg.db, INPUT_pathview, 'ENTREZID', 'SYMBOL'))

setwd('/Users/isarnassiri/Documents/Parkinson/scRNAseq-Parkkinen/Figures/')
library("pathview")
###################################################
### code chunk number 15: kegg.native
###################################################
i <- 1
pv.out <- pathview(gene.data = names(INPUT_pathview), pathway.id = '05012',
                   species = "hsa", out.suffix = "PD", kegg.native = T)
list.files(pattern="hsa05012", full.names=T)
str(pv.out)
head(pv.out$plot.data.gene)


###################################################
### code chunk number 16: kegg.native_2layer
###################################################
pv.out <- pathview(gene.data = names(INPUT_pathview), pathway.id = '05012',
                   species = "hsa", out.suffix = "PD.2layer", kegg.native = T,
                   same.layer = F, high=list(gene="gold"))

