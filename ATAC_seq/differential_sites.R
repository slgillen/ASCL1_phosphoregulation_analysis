#this script conducts DiffBind analysis to identify differentially accessible regions
#using R version: 4.1.2

# load libraries --------------------
library(DiffBind) #version 3.4
library(BiocParallel)

datadir<-'DiffBind_outputs/'


# read in data ------------------------------------------------------------
samples <- read.csv("ATAC_DiffBind_SampleSheet.csv",stringsAsFactors=FALSE) 
print(names(samples))

ASCL1_ATAC <- dba(sampleSheet=samples,minOverlap=1)


# get full consensus peakset across conditions --------------------
ATAC_consensus <-  dba.peakset(ASCL1_ATAC, consensus = DBA_CONDITION, minOverlap=2)
ExptConsensus <-  dba(ATAC_consensus, mask=ATAC_consensus$masks$Consensus,minOverlap=1) 
ExptConsensus
ConsensusPeaks <- dba.peakset(ExptConsensus, bRetrieve=TRUE,DataType=DBA_DATA_FRAME,minOverlap=1)
ConsensusPeaksx<-ConsensusPeaks[,1:3]
nrow(ConsensusPeaksx) 
write.table(ConsensusPeaksx,paste0(datadir,'ATAC_ConsensusPeaks_allin2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")



# get peakset counts for consensus regions defined --------------------
#sort config 
ASCL1_ATAC$config$singleEnd <- FALSE #using paired-end data
ASCL1_ATAC$config$cores <- 20
ASCL1_ATAC$config$doBlacklist <- FALSE
ASCL1_ATAC$config$doGreylist <- FALSE
ASCL1_ATAC$config

#get the counts
ConsensusPeaks_v2 <- dba.peakset(ExptConsensus, bRetrieve=TRUE)
ASCL1_ATAC_count <- dba.count(ASCL1_ATAC,peaks=ConsensusPeaks_v2,bParallel=TRUE,summits=TRUE,filter=1,minCount=0,bUseSummarizeOverlaps=TRUE) 
ASCL1_ATAC_count
reads <- dba.peakset(ASCL1_ATAC_count, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
apply(reads[,4:23],2,function(x) sum(x))

write.csv(reads, file=paste0(datadir,'ATAC_inconsensuspeaks_reads.csv'))


# PCA plots with consensus data ----------------------------------------------
info <- dba.show(ASCL1_ATAC_count)
write.table(data.frame(info,stringsAsFactors=FALSE),paste0(datadir,'info_ATAC.txt'),col.names=TRUE,row.names=FALSE,quote=FALSE)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
print(libsizes)

tiff(paste0(datadir,'ATAC_PCA_normalised_Consensus.tiff'),width=1400,height=1400,res=300)
dba.plotPCA(ASCL1_ATAC_count,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('black','grey48','royalblue3','darkorange','mediumvioletred'),score=DBA_SCORE_NORMALIZED)
dev.off()


  
# data normalisation ----------------------------------------------
ATAC_normalize<-dba.normalize(ASCL1_ATAC_count,method=DBA_DESEQ2)
print(ATAC_normalize$norm)
normlibs <- cbind(FullLibSize=ATAC_normalize$norm$DESeq2$lib.sizes, NormFacs=ATAC_normalize$norm$DESeq2$norm.facs,NormLibSize=round(ATAC_normalize$norm$DESeq2$lib.sizes/ATAC_normalize$norm$DESeq2$norm.facs))
rownames(normlibs) <- info$ID
print(normlibs)


# conduct differential analysis ----------------------------------------------
ATAC_contrast <- dba.contrast(ATAC_normalize,design="~Replicate+Condition")

ATAC_analyze <- dba.analyze(ATAC_contrast,method = DBA_DESEQ2,bGreylist=FALSE,bBlacklist=FALSE) 
dba.show(ATAC_analyze, bContrasts=TRUE)
write.table(dba.show(ATAC_analyze, bContrasts=TRUE),file=paste0(datadir,'ATAC_analyze_info.txt'),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)



# write outputs for the different comparisons ----------------------------------------------
report <- dba.report(ATAC_analyze, contrast = 7, th = 1, bFlip = T)
write.table(report, paste0(datadir,"d3ESCvd2ESC_ConsensusSummits.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(ATAC_analyze, contrast = 11, th = 1, bFlip = T)
write.table(report, paste0(datadir,"d3WT05vd3ESC_ConsensusSummits.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(ATAC_analyze, contrast = 12, th = 1, bFlip = T) 
write.table(report, paste0(datadir,"d3SA01vd3ESC_ConsensusSummits.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(ATAC_analyze, contrast = 13, th = 1, bFlip = T) 
write.table(report, paste0(datadir,"d3SA05vd3ESC_ConsensusSummits.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(ATAC_analyze, contrast = 14, th = 1, bFlip = T) 
write.table(report, paste0(datadir,"d3SA01vd3WT05_ConsensusSummits.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(ATAC_analyze, contrast = 15, th = 1, bFlip = T) 
write.table(report, paste0(datadir,"d3SA05vd3WT05_ConsensusSummits.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(ATAC_analyze, contrast = 16, th = 1, bFlip = F) 
write.table(report, paste0(datadir,"d3SA05vd3SA01_ConsensusSummits.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

 
