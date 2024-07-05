# module purge
# module add R/4.1.0-foss-2021a
# module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0
#---------------
Pair=1
path = '/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/DCOA/SAVER/'
OutPutPath = '/Users/isarnassiri/Documents/Parkinson/GeoMx_SC_Integration/DCOA/DCOXA/'
dir.create(OutPutPath, recursive = T)
#---------------

setwd(path)

library(stringr)
library(data.table)

library(DGCA) #install.packages('DGCA')
library(grDevices)
library(grid)
library(flexmix)
library(mvtnorm)
library(glmnet)

#=========
listfiles = list.files( pattern = '^SAVER.txt$')
length(listfiles)

# list files in folder matching an exact filename
# ^ is the regular expression denoting the beginning of the string,
# \\ escapes the . to make it a literal .,
# $ is the regular expression denoting the end of the string.

MetaData <- fread(paste0(path, 'Metadata_SAVER.txt'), stringsAsFactors = F, header = T)
MetaData = as.data.frame(MetaData)
dim(MetaData)

#========= read inputs

estimate = fread(listfiles[Pair], stringsAsFactors = F, header = T)
estimate = as.data.frame(estimate)

row.names(estimate) = make.names(estimate$GeneID, unique = T)
estimate = estimate[,-c(1)]
estimate[1:5,1:5]
dim(estimate)

colnames(estimate) = gsub('\\.','-', colnames(estimate))

S1_BARCODE = colnames(estimate)[which(colnames(estimate) %in% MetaData$Barcode[MetaData$yan == 'Vulnerable'])]
S2_BARCODE = colnames(estimate)[which(colnames(estimate) %in% MetaData$Barcode[MetaData$yan == 'Resistant'])]

length(S1_BARCODE)
length(S2_BARCODE)

#========= consider all genes one by one to find the DCGs

considered_GenePairs = vector()
i=1
for (i in 1:dim(estimate)[1]) {

tryCatch({

SelectedGene=row.names(estimate)[i]
AllGenes=row.names(estimate)[-which(row.names(estimate) == SelectedGene)]

#========= variable selection (select relevant genes to the indicated gene)

#-- Lasso
x = t(estimate[-which(row.names(estimate) == SelectedGene),])
dim(x)

y = as.numeric(estimate[which(row.names(estimate) == SelectedGene),])

set.seed(123)
train=sample(seq(dim(x)[1]),(dim(x)[1]/2), replace=FALSE)

lasso.tr=glmnet(x[train,],y[train], standardize=F)
pred=predict(lasso.tr,x[-train,])
rmse= sqrt(apply((y[-train]-pred)^2,2,mean))
lam.best=lasso.tr$lambda[order(rmse)[1]]
lam.best

#------
lasso.tr=glmnet(x ,y, standardize=F)
#------

Lasso_coefficient <- coef(lasso.tr, s=lam.best)
inds<-which(Lasso_coefficient!=0)
variables<-row.names(Lasso_coefficient)[inds]
variables<-variables[!(variables %in% '(Intercept)')];

length(which(0 != Lasso_coefficient[,1]))
results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
colnames(results) <- "coef"
results <- as.data.frame(results[-1,,FALSE])# FALSE is about inactivation of drop paremters
dim(results)[1]

AllGenes = row.names(results)

print(length(AllGenes))
print(i)
# }, error=function(e){})
#   
# }
#==========================================

RESULT_merged = data.frame()

if(length(AllGenes) > 0)
{
for(j in 1:length(AllGenes))
{
Gene1=SelectedGene;
Gene2=AllGenes[j];

if( length(considered_GenePairs[which(considered_GenePairs %in% c(paste(Gene1, Gene2, sep = '_'), paste(Gene2, Gene1, sep = '_')))]) > 0)
{
  print('already is considered.')
}else{

#---------- save considered GenePairs
considered_GenePairs[length(considered_GenePairs)+1] = paste(Gene1, Gene2, sep = '_')
considered_GenePairs[length(considered_GenePairs)+1] = paste(Gene2, Gene1, sep = '_')
#----------

estimate_Gene1 = estimate[which(row.names(estimate) == Gene1),]
estimate_Gene2 = estimate[which(row.names(estimate) == Gene2),]
#----------

input_mixture_model = t(rbind(as.matrix(estimate_Gene1), as.matrix(estimate_Gene2)))
input_mixture_model = as.data.frame(input_mixture_model, stringsAsFactors = F)

identical(row.names(input_mixture_model), colnames(estimate))

#------------------------------------------------------------ remove outliers
#' Detect outliers using IQR method
#' 
#' @param x A numeric vector
#' @param na.rm Whether to exclude NAs when computing quantiles
#' 
is_outlier <- function(x, na.rm = FALSE) {
  qs = quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  
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
dim(input_mixture_model)
vars_of_interest <- c("Gene1", "Gene2")
colnames(input_mixture_model)[c(1,2)] <- c("Gene1", "Gene2")
input_mixture_model <- remove_outliers(input_mixture_model, vars_of_interest)
#------------------------------------------------------------

mixed_model = flexmix(cbind(as.numeric(input_mixture_model[,1]),as.numeric(input_mixture_model[,2]))~1, 
k=2, data=input_mixture_model,
model = FLXMCmvnorm(diag = T), #diag = T to get fix results FLXMCmvpois(), FLXMCmvnorm
control = list(tolerance = 1e-15, iter.max = 1000))

# print(mixed_model)
# posterior(mixed_model)
# parameters(mixed_model)

# print(length(unique(clusters(mixed_model))))
# print(j)
# }

if(length(unique(clusters(mixed_model)))==2)
{
#========================================== DGCA
#------ real clusters
# design_matrix = cbind(input_mixture_model[,'predicted_clusters'], input_mixture_model[,'predicted_clusters'])
design_matrix = cbind(colnames(estimate), colnames(estimate))
design_matrix[-which(design_matrix[,1]%in%S1_BARCODE),1]=0
design_matrix[which(design_matrix[,1]%in%S1_BARCODE),1]=1
design_matrix[-which(design_matrix[,2]%in%S2_BARCODE),2]=0
design_matrix[which(design_matrix[,2]%in%S2_BARCODE),2]=1
colnames(design_matrix) <- c("G0", "G1")

design_matrix = as.data.frame(design_matrix)
str(design_matrix)
design_matrix[,1] = as.numeric(design_matrix[,1])
design_matrix[,2] = as.numeric(design_matrix[,2])
design_matrix = as.matrix(design_matrix)
ddcor_G1_G0 = ddcorAll(inputMat = estimate[which(row.names(estimate) %in% c(Gene1,Gene2)),], design = design_matrix,compare = c("G1","G0"),adjust = "fdr", nPerm = 0, corrType = "spearman", splitSet = Gene1, sigThresh = 0.05, sortBy = "pValDiff_adj", verbose = TRUE)
head(ddcor_G1_G0) 

#--------------- intersection between real and predicted clusters
input_mixture_model$real_clusters = row.names(input_mixture_model)
input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S1_BARCODE)]=1
input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S2_BARCODE)]=2

input_mixture_model$predicted_clusters = clusters(mixed_model)

if(input_mixture_model$predicted_clusters[1] == 2)
{
  input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 1)]=0
  input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 2)]=1
  input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 0)]=2
}

table(input_mixture_model$predicted_clusters)
table(input_mixture_model$real_clusters)

table = table(input_mixture_model$real_clusters,input_mixture_model$predicted_clusters)
TPrate = (table[1,1] + table[2,2])/length(clusters(mixed_model))
print(TPrate)

SM = summary(mixed_model)

#AIC: akaike information criterion. Measure model quality, wherein low values are better (if it is negative, more negative is better).
#We use the BIC to select optimal number of factor.Lower BIC always prefer.
#---------------
RESULT = data.frame()
RESULT = ddcor_G1_G0#[which(ddcor_G1_G0$Gene1 == Gene2), ]
RESULT$AIC = SM@AIC
RESULT$TPrate = TPrate
print(RESULT)
RESULT = RESULT[which(RESULT$pValDiff_adj < 1e-4 ),]# & RESULT$G0_pVal != 0.0
# RESULT = RESULT[which(RESULT$Classes %in% c('-/+', '+/-')),]

if(dim(RESULT)[1]>0)
{
if(j==1){RESULT_merged = RESULT}
if(j!=1 & dim(RESULT_merged)[1] != 0){RESULT_merged = rbind(RESULT_merged, RESULT)}
if(j!=1 & dim(RESULT_merged)[1] == 0){RESULT_merged = RESULT}
rm(list = c('RESULT', 'TPrate', 'ddcor_G1_G0', 'SM'))
}
print(RESULT_merged)

}# if cluster > 2
} # considered_Pairs
}#for j

#========================================== 
if(dim(RESULT_merged)[1] != 0 )
{
#========================================== save results
  write.table(RESULT_merged, paste0(OutPutPath, '/', Gene1, '_DCOG.txt'), quote = F, row.names = F, sep = '\t')
} 

}

}, error=function(e){})

}# for i
 
