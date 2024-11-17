# create directories ------------------------------------------------------

#create output directory
dir.create('downstream_analysis/differential_expression')

outdir<-'downstream_analysis/differential_expression/'
indir<-'counts/'

# load libraries ----------------------------------------------------------
library(DESeq2)
library(limma)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(plyr)
library(clusterProfiler)
library(ComplexHeatmap)


# read in data ------------------------------------------------------------

RNAseq_counts<-read.csv('all_counts_WT_SA_ASCL1_filtered.csv')
sample_info<-read.csv('sample_information.csv')

#factorise sample_info
sample_info$condition<-factor(sample_info$condition,levels=c('ESC','WTASCL105','SAASCL101','SAASCL105'))
sample_info$replicate<-factor(sample_info$replicate)
str(sample_info)

# generate correlation matrix and PCA plot ---------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=counts_data,colData=sample_info,design=~replicate+condition)
dds <- dds[ rowSums(counts(dds)) > 15, ] 

# normalisation and preprocessing
dds <- DESeq(dds)

# VST normalisation
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd))) 

# Correlation plot of all samples
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

tiff(paste0(outdir1,'vst_corrPlot_all.tiff'), height = 2000, width = 2200,res=300)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
dev.off()


# plot PCAs for different condition combinations
generate_PCA <- function(vsdx,groups,colourset,setname,outfolder) {
  vsd.sub <- vsdx[ , vsdx$condition %in% groups]
  pca_data <- DESeq2::plotPCA(vsd.sub, intgroup = c( "replicate", "condition"), returnData = TRUE) 
  percent_var <- round(100 * attr(pca_data, "percentVar")) 
  
  f1<-ggplot(pca_data, aes(x = PC1, y = PC2, fill = factor(condition),col = factor(condition), shape = factor(replicate))) + geom_point(size =3) + 
    scale_shape_manual(values=c(21,22,23,24)) + scale_fill_manual(values=colourset) + scale_colour_manual(values=colourset) + 
    xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) + theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text=element_text(size=12))
  ggsave(plot=f1,file=paste0(outfolder,setname,"_PCA.tiff"),width=3.75,height=2.5,dpi=300)
  
}

setnames<-c('all','ASCL1')
groupsets<-list(c('ESC','WTASCL105','SAASCL101','SAASCL105'),c('WTASCL105','SAASCL101','SAASCL105'))
coloursets<-list(c('grey58','royalblue3','darkorange','mediumvioletred'),c('royalblue3','darkorange','mediumvioletred'))

for(i in 1:length(groupsets)){
  generate_PCA(vsd,groupsets[[i]],coloursets[[i]],setnames[i],outdir1)
}


# conduct differential expression analysis ----------------------------------------

counts_matrix<-filtered_counts 
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

#factorise sample_info
sample_info$condition<-factor(sample_info$condition,levels=c('ESC','WTASCL105','SAASCL101','SAASCL105'))
sample_info$replicate<-factor(sample_info$replicate)
str(sample_info)

colnames(counts_matrix)==sample_info$sample

# do the analysis
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = sample_info,design= ~ condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 

#WT ASCL1 05 v control
resLFC <- lfcShrink(DESeq2output, coef="condition_WTASCL105_vs_ESC", type="apeglm") 

#SA ASCL1 01 v control
resLFC2 <- lfcShrink(DESeq2output, coef="condition_SAASCL101_vs_ESC", type="apeglm") 

#SA ASCL1 05 v control
resLFC3 <- lfcShrink(DESeq2output, coef="condition_SAASCL105_vs_ESC", type="apeglm") 

rm(DESeq2data,DESeq2output)


# refactorise to run different comparison
sample_info$condition<-factor(sample_info$condition,levels=c('WTASCL105','ESC','SAASCL101','SAASCL105'))
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = sample_info,design = ~ condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 


#SA ASCL1 01 v WT ASCL1 05 
resLFC4 <- lfcShrink(DESeq2output, coef="condition_SAASCL101_vs_WTASCL105", type="apeglm") 

#SA ASCL1 05 v WT ASCL1 05 
resLFC5 <- lfcShrink(DESeq2output, coef="condition_SAASCL105_vs_WTASCL105", type="apeglm") 

rm(DESeq2data,DESeq2output)


# refactorise to run different comparison
sample_info$condition<-factor(sample_info$condition,levels=c('SAASCL101','ESC','WTASCL105','SAASCL105'))
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = sample_info,design = ~ condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 


#SA ASCL1 05 v WT ASCL1 01 
resLFC6 <- lfcShrink(DESeq2output, coef="condition_SAASCL105_vs_SAASCL101", type="apeglm") 

rm(DESeq2data,DESeq2output)


# sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

resOrdered2 <- resLFC2[order(resLFC2$pvalue),]
summary(resLFC2)
sum(resLFC2$padj < 0.05, na.rm=TRUE)

resOrdered3 <- resLFC3[order(resLFC3$pvalue),]
summary(resLFC3)
sum(resLFC3$padj < 0.05, na.rm=TRUE)

resOrdered4 <- resLFC4[order(resLFC4$pvalue),]
summary(resLFC4)
sum(resLFC4$padj < 0.05, na.rm=TRUE)

resOrdered5 <- resLFC5[order(resLFC5$pvalue),]
summary(resLFC5)
sum(resLFC5$padj < 0.05, na.rm=TRUE)

resOrdered6 <- resLFC6[order(resLFC6$pvalue),]
summary(resLFC6)
sum(resLFC6$padj < 0.05, na.rm=TRUE)



# quick look at results
png(paste0(outdir1,'MAplot_WTASCL105vESC.png'),width=1400,height=1200,res=300)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png(paste0(outdir1,'MAplot_SAASCL101vESC.png'),width=1400,height=1200,res=300)
DESeq2::plotMA(resLFC2, ylim=c(-4,4))
dev.off()

png(paste0(outdir1,'MAplot_SAASCL105vESC.png'),width=1400,height=1200,res=300)
DESeq2::plotMA(resLFC3, ylim=c(-4,4))
dev.off()

png(paste0(outdir1,'MAplot_SAASCL101vWTASCL105.png'),width=1400,height=1200,res=300)
DESeq2::plotMA(resLFC4, ylim=c(-4,4))
dev.off()

png(paste0(outdir1,'MAplot_SAASCL105vWTASCL105.png'),width=1400,height=1200,res=300)
DESeq2::plotMA(resLFC5, ylim=c(-4,4))
dev.off()

png(paste0(outdir1,'MAplot_SAASCL105vSAASCL101.png'),width=1400,height=1200,res=300)
DESeq2::plotMA(resLFC6, ylim=c(-4,4))
dev.off()



# results table formatting
results1<-as.data.frame(resOrdered)
nrow(results1)
results1$gene<-rownames(results1)
results1<-results1[,c(6,1:5)]

results2<-as.data.frame(resOrdered2)
nrow(results2)
results2$gene<-rownames(results2)
results2<-results2[,c(6,1:5)]

results3<-as.data.frame(resOrdered3)
nrow(results3)
results3$gene<-rownames(results3)
results3<-results3[,c(6,1:5)]

results4<-as.data.frame(resOrdered4)
nrow(results4)
results4$gene<-rownames(results4)
results4<-results4[,c(6,1:5)]

results5<-as.data.frame(resOrdered5)
nrow(results5)
results5$gene<-rownames(results5)
results5<-results5[,c(6,1:5)]

results6<-as.data.frame(resOrdered6)
nrow(results6)
results6$gene<-rownames(results6)
results6<-results6[,c(6,1:5)]


# write results
write.csv(results1, file=paste0(outdir1,"WTASCL105vESC_DESeq2output.csv"),row.names=FALSE,quote=FALSE)
write.csv(results2, file=paste0(outdir1,"SAASCL101vESC_DESeq2output.csv"),row.names=FALSE,quote=FALSE)
write.csv(results3, file=paste0(outdir1,"SAASCL105vESC_DESeq2output.csv"),row.names=FALSE,quote=FALSE)
write.csv(results4, file=paste0(outdir1,"SAASCL101vWTASCL105_DESeq2output.csv"),row.names=FALSE,quote=FALSE)
write.csv(results5, file=paste0(outdir1,"SAASCL105vWTASCL105_DESeq2output.csv"),row.names=FALSE,quote=FALSE)
write.csv(results6, file=paste0(outdir1,"SAASCL105vSAASCL101_DESeq2output.csv"),row.names=FALSE,quote=FALSE)



# collate results ---------------------------------------------------------

#combine differential expression tables from the different comparisons conducted
collate_results<-function(results_list,results_names){
  for(i in 1:length(results_list)){
    if(i==1){
      all_results<-results_list[[i]][,c('gene','log2FoldChange','padj')]
      names(all_results)[2:3]<-paste0(gsub('_DESeq2output.csv','',results_names[i]),'_',names(all_results)[2:3])
    }else{
      resultsx<-results_list[[i]][,c('gene','log2FoldChange','padj')]
      names(resultsx)[2:3]<-paste0(gsub('_nDESeq2output.csv','',results_names[i]),'_',names(resultsx)[2:3])
      all_results<-join(all_results,resultsx,by='gene',type='full')
      rm(resultsx)
    }
  }
  return(all_results)
} 


filenames <- dir(outdir, pattern = "*_DESeq2output.csv")
length(filenames)
filenames

allinputs<-lapply(as.list(filenames),function(x) read.delim(paste0(data_path,x),sep=',',stringsAsFactors = FALSE,header=TRUE))
length(allinputs)

all_DE<-collate_results(allinputs, filenames)
write.csv(all_DE,'all_DE_collated.csv',quote=FALSE,row.names=FALSE,sep=',')


