
## script for biost 578 final project

## Contents:
# 1. Get the data from SRA database
# 2. Align using bowtie
# 3. Read in unaligned seqs; make Figure S1
# 4. Convert SAM to BAM
# 5. Plot aligned reads to check results
# 6. Filter aligned reads
# 7. Try qAlign
# 8. Read in counts results
# 9. Make data.frame of counts
# 10. Make some plots, do FQ normalization
# 11. Read in website data, make exprSet object
# 12. Make plots, do w/in lane normalization
# 13. Do pseudo-datasets, make plots


library(EDASeq)
library(affy); library(edgeR)


### 
# 1. Get the data from SRA database

library(SRAdb) # allows us to query the SRA database

sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)

dbListTables(sra_con)
dbListFields(sra_con,"study"); dbListFields(sra_con,"sra")
rs <- dbGetQuery(sra_con,"select * from study limit 3")
rs <- dbGetQuery(sra_con,"select * from sra limit 3")
rs <- dbGetQuery(sra_con,"select * from submission limit 3")

# yeast data: SRA048710
# MAQC data: SRA010153
rs = listSRAfile( c("SRA048710","SRA010153"), sra_con, fileType = 'sra' )
dim(rs) # 56 6
fqinfo <- getFASTQinfo( c("SRA048710","SRA010153"), srcType = 'fasp' ) # there are 56 files
getSRAfile( c("SRA048710","SRA010153"), sra_con, fileType = 'fastq' ) # took forever!!
#Files are saved to: 
#  '/Users/c8linmch'

# 26 and 32 cant be unzipped, so need to download them again
getSRAfile(c("SRR390926"),sra_con,fileType='fastq')
getSRAfile(c("SRR390932"),sra_con,fileType='fastq')


###
# 2. Align using bowtie

# align using bowtie, w unique mapping and up to 2 mismatches
# read-count per gene is defined as # reads w 5-end falling w/in region
# genes w avg read count <10 for each of 3 growth conditions are filtered out
# retained 5690/6575 genes in YEAST

# did this from the tower computer via batch.

biocLite("Rbowtie")
library(Rbowtie)

#refFiles <- dir(system.file(package="Rbowtie", "samples", "refs"), full=TRUE)
refFiles <- "/Users/c8linmch/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa"
indexDir <- file.path("/Users/c8linmch/S288C_reference_genome_R64-1-1_20110203/refsIndex")
tmp <- bowtie_build(references=refFiles, outdir=indexDir, prefix="index", force=TRUE) 
head(tmp)
#[1] "Settings:"
#[2] " Output files: \"/tmp/RtmpbsK0To/refsIndex/index.*.ebwt\""
#[3] " Line rate: 6 (line is 64 bytes)"
#[4] " Lines per side: 1 (side is 64 bytes)"
#[5] " Offset rate: 5 (one in 32)"
#[6] " FTable chars: 10"

#biocLite("BSgenome.Scerevisiae.UCSC.sacCer1")
#library("BSgenome.Scerevisiae.UCSC.sacCer1")
# its called Scerevisiae
# downloaded ref sequence from: 
# http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tar

#from command line at the tower:
#wget "http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz" 

# readsFiles is the name of a file containing short reads to be aligned with bowtie, 
# and samFiles specifies the name of the output file with the generated alignments
## need to do this one by one!!

# called from pearson at the tower:
# qsub -N bowtie_align -e mchughc@uw.edu bowtie_align.sh


###
# 3. Read in unaligned seqs; make Figure S1

files <- list.files("/Users/c8linmch/yeast",full.names=TRUE)
# exclude the alignments.sam file
files <- files[-1]
names(files) <- fqinfo$library.name

fqinfo <- getFASTQinfo( c("SRA048710"), srcType = 'fasp' ) # there are 56 files
met <- DataFrame(conditions=fqinfo$library.name,row.names=files)
fastq <- FastqFileList(files)
elementMetadata(fastq) <- met
fastq

pdf("/Users/c8linmch/biost_578/fig_S1a.pdf")
barplot(fastq,col=c("blue",rep("red",8),"blue","blue",rep("red",3)))
legend("topright",c("Protocol 1","Protocol 2"),fill=c("red","blue"))
dev.off()


###
# 4. Convert SAM to BAM

# called from pearson at the tower:
# qsub -N sam_bam -e mchughc@uw.edu sam_bam.sh


###
# 5. Plot aligned reads to check results

biocLite("EDASeq")
library(EDASeq); library(Rsamtools)
library(SRAdb) # allows us to query the SRA database

# make a data.frame of the metadata, save it so dont have to keep
# querying the database. takes too long!
#sqlfile <- getSRAdbFile()
#sra_con <- dbConnect(SQLite(),sqlfile)
#fqinfo <- getFASTQinfo( c("SRA048710"), srcType = 'fasp' ) 
#save(fqinfo,file="/home/staff/mchughc/biost578/wq2014/yeast_info.RData")

fqinfo <- get(load("/home/staff/mchughc/biost578/wq2014/yeast_info.RData"))
met <- DataFrame(conditions=fqinfo$library.name,row.names=fqinfo$run)


files <- list.files("/home/staff/mchughc/biost578/wq2014/yeast",full.names=TRUE)
bais <- grep(".bai",files)
files <- files[-bais]
bams <- files[grep(".bam",files)]

bfs <- BamFileList(bams)
elementMetadata(bfs) <- met
bfs
pdf("/home/staff/mchughc/biost578/wq2014/fig_S1b.pdf")
barplot(bfs,las=2,col=c("blue",rep("red",8),"blue","blue",rep("red",3)))
legend("topright",c("Protocol 1","Protocol 2"),fill=c("red","blue"))
dev.off()


###
# 6. Filter aligned reads

# read count is num of reads with 5prime end falling w/in gene region
# genes w avg read ct <10 for each of 3 growth conditions were filtered out
# ie gene j was filtered if max k in{ypd,del,gly} y bar_j,k <10, where y bar_j,k is avg
# read count for gene j in condition k

library(Rsamtools)

?summarizeOverlaps
# how do we get the genes w 5prime end falling in region? thats what we want to count

#ScanBamParam()
#countBam(bamFile,index=bamFile,param=ScanBamParam())

filterBam

# sacCer3=yeastgenome.org download of current version
#biocLite("BSgenome.Scerevisiae.UCSC.sacCer3")
library("BSgenome.Scerevisiae.UCSC.sacCer3")
Scerevisiae
# its called Scerevisiae

biocLite("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

exbygene <- exonsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, "gene")

# we are now ready to count the number of sequences that fall within each gene
bamFile <- "/home/staff/mchughc/biost578/wq2014/yeast/SRR390921.bam.bam"
se <- summarizeOverlaps(exbygene, bamFile, mode="IntersectionNotEmpty")
se 
# object of class `SummarizedExperiment`.

# extract the counts
head(table(assays(se)$counts)) # 0 for everything
rowData(se)

file1 <- readGAlignments("/home/staff/mchughc/biost578/wq2014/yeast/SRR390921.bam.bam")
gr1 <- as(file1,"GRanges")
cts <- countOverlaps(query=exbygene,subject=gr1,type="end") 
# returns the overlap hit count for each range in query

bamFiles <- list.files("/home/staff/mchughc/biost578/wq2014/yeast", "bam$", full=TRUE)
names(bamFiles) <- sub("\\..*","",basename(bamFiles))

aln <- readGAlignments(bamFiles[7])
seqlevels(aln)
seqlevels(exbygene)
newLevels <- seqlevels(aln)
seqlevels(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
# ok, great! we just need to rename the ref genome w our chr names so that we can get counts!

refGenes <- renameSeqlevels(exbygene,c(chrI=newLevels[1]))
refGenes <- renameSeqlevels(refGenes,c(chrII=newLevels[2]))
refGenes <- renameSeqlevels(refGenes,c(chrIII=newLevels[3]))
refGenes <- renameSeqlevels(refGenes,c(chrIV=newLevels[4]))
refGenes <- renameSeqlevels(refGenes,c(chrV=newLevels[5]))
refGenes <- renameSeqlevels(refGenes,c(chrVI=newLevels[6]))
refGenes <- renameSeqlevels(refGenes,c(chrVII=newLevels[7]))
refGenes <- renameSeqlevels(refGenes,c(chrVIII=newLevels[8]))
refGenes <- renameSeqlevels(refGenes,c(chrIX=newLevels[9]))
refGenes <- renameSeqlevels(refGenes,c(chrX=newLevels[10]))
refGenes <- renameSeqlevels(refGenes,c(chrXI=newLevels[11]))
refGenes <- renameSeqlevels(refGenes,c(chrXII=newLevels[12]))
refGenes <- renameSeqlevels(refGenes,c(chrXIII=newLevels[13]))
refGenes <- renameSeqlevels(refGenes,c(chrXIV=newLevels[14]))
refGenes <- renameSeqlevels(refGenes,c(chrXV=newLevels[15]))
refGenes <- renameSeqlevels(refGenes,c(chrXVI=newLevels[16]))
refGenes <- renameSeqlevels(refGenes,c(chrM=newLevels[17]))
seqlevels(refGenes); seqlevels(aln)
# great, everything matches now! fewf.

se <- summarizeOverlaps(refGenes, bamFiles[7], mode="IntersectionNotEmpty")
se 
# object of class `SummarizedExperiment`.

# extract the counts
head(table(assays(se)$counts))
#  0   1   2   3   4   5 
#716 180 108  91  90  57
rowData(se)

# great!!!! run in batch for all, save the se objects.

rm(list=ls())


###
# 7. Try qAlign

#qsub -N qAlign -e mchughc@uw.edu qAlign.sh
# didn't work. oh well.


###
# 8. Read in counts results

counts <- get(load("/projects/geneva/gcc-fs2/biost578_wq2014/all_counts.RData"))
dim(counts); head(counts) # 6534 14

counts <- data.frame(counts)

# get mean by growth condition
fqinfo <- get(load("/projects/geneva/gcc-fs2/biost578_wq2014/yeast_info.RData"))
fqinfo[,c("run","library.name")]

ct_g <- apply(counts[,1:3],1,mean)
ct_y <- apply(counts[,4:11],1,mean)
ct_d <- apply(counts[,12:14],1,mean)

counts$g_avg <- ct_g
counts$y_avg <- ct_y
counts$d_avg <- ct_d

counts$max <- apply(counts[,15:17],1,max)
sum(counts$max<10)
sum(counts$max>=10) # 6205

# redo counts w stricter version -- "Union" option
counts <- get(load("/projects/geneva/gcc-fs2/biost578_wq2014/all_countsUnion.RData"))
dim(counts); head(counts) # 6534 14

counts <- data.frame(counts)

ct_g <- apply(counts[,1:3],1,mean)
ct_y <- apply(counts[,4:11],1,mean)
ct_d <- apply(counts[,12:14],1,mean)

counts$g_avg <- ct_g
counts$y_avg <- ct_y
counts$d_avg <- ct_d

counts$max <- apply(counts[,15:17],1,max)
sum(counts$max<10) # 329
sum(counts$max>=10) # 6205
sum(counts$max>10) # 6198

counts$filter <- counts$max>10
table(counts$filter) # 6198 true

# redo counts w stricter version -- "IntersectionStrict" option
counts <- get(load("/projects/geneva/gcc-fs2/biost578_wq2014/all_countsIntersectStrict.RData"))
dim(counts); head(counts) # 6534 14

counts <- data.frame(counts)

ct_g <- apply(counts[,1:3],1,mean)
ct_y <- apply(counts[,4:11],1,mean)
ct_d <- apply(counts[,12:14],1,mean)

counts$g_avg <- ct_g
counts$y_avg <- ct_y
counts$d_avg <- ct_d

counts$max <- apply(counts[,15:17],1,max)
sum(counts$max<10) # 414
sum(counts$max>=10) # 6120
sum(counts$max>10) # 6111

counts$filter <- counts$max>10
table(counts$filter) # 6111 true


# redo counts w stricter version -- "IntersectionStrict" option
counts <- get(load("/projects/geneva/gcc-fs2/biost578_wq2014/all_countsIntStrictGene.RData"))
dim(counts); head(counts) # 6534 14

counts <- data.frame(counts)

ct_g <- apply(counts[,1:3],1,mean)
ct_y <- apply(counts[,4:11],1,mean)
ct_d <- apply(counts[,12:14],1,mean)

counts$g_avg <- ct_g
counts$y_avg <- ct_y
counts$d_avg <- ct_d

counts$max <- apply(counts[,15:17],1,max)
sum(counts$max<10) # 414
sum(counts$max>=10) # 6120
sum(counts$max>10) # 6111

counts$filter <- counts$max>10
table(counts$filter) # 6111 true

# well, this is what we have to go on, i guess.


###
# 9. Make data.frame of counts

library(EDASeq)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
gns <- genes(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

counts <- get(load("/projects/geneva/gcc-fs2/biost578_wq2014/all_countsIntStrictGene.RData"))
rownames(counts) <- names(gns)

data(yeastGC); data(yeastLength)
feature <- data.frame(gc=yeastGC,length=yeastLength)

mnG <- rowMeans(counts[,1:3])
mnY <- rowMeans(counts[,4:11])
mnD <- rowMeans(counts[,12:14])

mns <- data.frame(cbind(mnG,mnY,mnD))
mns$mx <- apply(mns,1,max)

filter <- mns$mx>10
table(filter)
#filter
#FALSE  TRUE
#423  6111

common <- intersect(names(yeastGC),rownames(counts[filter,]))
length(common) # 6111

data <- newSeqExpressionSet(exprs=as.matrix(counts[common,]),
                            featureData=feature[common,],
                            phenoData=data.frame(
                              conditions=c(rep("Gly",3),rep("YPD",8),rep("Del",3)),
                              library.prep=c(2,rep(1,8),2,2,1,1,1),
                              flow.cells=c("6247L","620AY","620AY",rep(c("428R1","4328B"),3),
                                rep("61MKN",2),rep("428R1",3)),
                              culture.prep=c(paste("G",1:3,sep=""),"Y1","Y1","Y2","Y2","Y7","Y7","Y4","Y4",
                                                                 "D1","D2","D7"),
                              row.names=colnames(counts)))
save(data,file="/projects/geneva/gcc-fs2/biost578_wq2014/geneExprSet_data.RData")

data2 <- newSeqExpressionSet(exprs=as.matrix(counts[common,]+1),
                             featureData=feature[common,],
                             phenoData=data.frame(
                               conditions=c(rep("Gly",3),rep("YPD",8),rep("Del",3)),
                               library.prep=c(2,rep(1,8),2,2,1,1,1),
                               flow.cells=c("6247L","620AY","620AY",rep(c("428R1","4328B"),3),
                                 rep("61MKN",2),rep("428R1",3)),
                               culture.prep=c(paste("G",1:3,sep=""),"Y1","Y1","Y2","Y2","Y7","Y7","Y4","Y4",
                                 "D1","D2","D7"),
                               row.names=colnames(counts)))
save(data2,file="/projects/geneva/gcc-fs2/biost578_wq2014/geneExprSet_data_plusOne.RData")

pdf("/projects/geneva/gcc-fs2/biost578_wq2014/boxplot_laneCounts.pdf")
boxplot(data,xlab=pData(data)$culture.prep)
dev.off()

pdf("/projects/geneva/gcc-fs2/biost578_wq2014/boxplot_laneCounts_plus1.pdf")
boxplot(data2)
dev.off()

rm(list=ls())


###
# 10. Make some plots, do FQ normalization

data <- get(load("/Users/c8linmch//Dropbox/biost578_finalProject/geneExprSet_data.RData"))
data2 <- get(load("/Users/c8linmch//Dropbox/biost578_finalProject/geneExprSet_data_plusOne.RData"))

# FQ between-lane normalization
fq_bw <- betweenLaneNormalization(data2,which="full")

pdf("boxplot_bwLane_FQ_plusOne.pdf")
colCond <- rep("blue",14)
colCond[pData(fq_bw)$conditions=="YPD"] <- "red"
colCond[pData(fq_bw)$conditions=="Del"] <- "green"
boxplot(fq_bw,axes=F,col=colCond,main="FQ between-lane normalization")
axis(2)
axis(1,lab=as.character(pData(fq_bw)$culture.prep),at=1:14)
legend("topleft",c("YPD","Del","Gly"),fill=c("red","green","blue"),bg="white")
dev.off()

pdf("boxplot_plusOne.pdf")
colCond <- rep("blue",14)
colCond[pData(data2)$conditions=="YPD"] <- "red"
colCond[pData(data2)$conditions=="Del"] <- "green"
boxplot(data2,axes=F,col=colCond,main="Unnormalized")
axis(2)
axis(1,lab=as.character(pData(data2)$culture.prep),at=1:14)
legend("topleft",c("YPD","Del","Gly"),fill=c("red","green","blue"),bg="white")
dev.off()

data_exclProt2 <- data2[,pData(data2)$library.prep==1]
dim(exprs(data_exclProt2)) # 6111 11
pData(data_exclProt2)

fq_bw_exclProt2 <- betweenLaneNormalization(data_exclProt2,which="full")
pdf("boxplot_bwLane_FQ_exclProt2_plusOne.pdf")
colCond <- rep("blue",11)
colCond[pData(fq_bw_exclProt2)$conditions=="YPD"] <- "red"
colCond[pData(fq_bw_exclProt2)$conditions=="Del"] <- "green"
boxplot(fq_bw_exclProt2,axes=F,col=colCond,main="FQ between-lane normalization\nExcl Protocol 2 Samples")
axis(2)
axis(1,lab=as.character(pData(fq_bw_exclProt2)$culture.prep),at=1:11)
legend("topleft",c("YPD","Del","Gly"),fill=c("red","green","blue"),bg="white")
dev.off()



ypds <- fq_bw[,pData(fq_bw)$conditions=="YPD"]

biasPlot(ypds[,pData(ypds)$library.prep!=2],"gc",log=T,ylim=c(0,6),xlab="GC-content",ylab="log(count+1)",
         col=c(rep("black",2),rep("red",2),rep("green",2)))
legend("topleft",c("Y1","Y2","Y4"),fill=c("black","red","green"))

# stratified boxplot of count log-ratio vs gc-content
# YPD Y1 flow cell 428R1 vs 4328B
pdf("boxplot_fig2a.pdf")
lfc <- log(exprs(fq_bw)[,4]+0.5)-log(exprs(fq_bw)[,5]+0.5)
biasBoxplot(lfc,fData(fq_bw)$gc,ylab="log-fold-change")
dev.off()

# YPD Y1 vs Y2 lane 428R1
pdf("boxplot_fig2b.pdf")
lfc <- log(exprs(fq_bw)[,6]+0.5)-log(exprs(fq_bw)[,4]+0.5)
biasBoxplot(lfc,fData(fq_bw)$gc,ylab="log-fold-change")
dev.off()

fq_bw <- betweenLaneNormalization(data,"full")
lfc <- log(exprs(fq_bw)[,4]+1)
biasBoxplot(lfc,fData(fq_bw)$gc,ylab=c("log(count+1)"))

dataWithin <- withinLaneNormalization(data,"gc",which="full")
dataNorm <- betweenLaneNormalization(dataWithin,which="full")
lfc <- log(exprs(dataNorm)[,4]+1)
biasBoxplot(lfc,fData(dataNorm)$gc)

rm(list=ls())


###
# 11. Read in website data, make exprSet object

# data downloaded from http://cgrlucb.wikispaces.com/RNASeqSpring2012
# load in, filter genes based on avg read count
geneLevelCounts <- read.table("yeastRNASeqRisso2011/geneLevelCounts.txt",header=TRUE,as.is=TRUE)
dim(geneLevelCounts) # 6575 14

laneInfo <- read.table("yeastRNASeqRisso2011/laneInfo.txt",header=TRUE,as.is=TRUE)
laneInfo

library(ShortRead)
fa <- readFasta("yeastRNASeqRisso2011/Scer.fasta")
abc <- alphabetFrequency(sread(fa), baseOnly=TRUE)
rownames(abc) <- sapply(strsplit(as.character(id(fa))," "),function(x) x[1])
alphabet <- abc[,1:4]
gc <- rowSums(alphabet[,2:3])/rowSums(alphabet)
length <- width(sread(fa))
head(gc)
head(length)

geneInfo <- data.frame(length=length, gc=gc)
head(geneInfo)

meansY <- rowMeans(geneLevelCounts[,laneInfo$conditions=="YPD"])
meansG <- rowMeans(geneLevelCounts[,laneInfo$conditions=="Gly"])
meansD <- rowMeans(geneLevelCounts[,laneInfo$conditions=="Del"])
mns <- data.frame("Y"=meansY,"G"=meansG,"D"=meansD)
filter <- apply(mns,1,function(x){max(x)>=10})
table(filter) # 5690 true; great!
geneLevelCounts <- geneLevelCounts[filter,]

library(EDASeq)
data <- newSeqExpressionSet(exprs = as.matrix(geneLevelCounts),
                            featureData = geneInfo[rownames(geneLevelCounts), ],
                            phenoData = laneInfo)
data
head(exprs(data))
pData(data)
head(fData(data))

save(data,file="geneExprSet_data_fromWebsite.RData")

rm(list=ls())


###
# 12. Make plots, do w/in lane normalization

data <- get(load("geneExprSet_data_fromWebsite.RData"))
exprs(data) <- exprs(data)+1

colors <- rep("red",14)
colors[pData(data)$conditions=="Del"] <- "green"
colors[pData(data)$conditions=="Gly"] <- "blue"

png("figure_s2_onlineData.png",width=960)
par(mfrow=c(1,2))
boxplot(data, col=colors,ylab="log(count+1)",main="Unnormalized",axes=F)
axis(2)
axis(1,lab=pData(data)$lib_prep,at=1:14)
legend("topleft",c("YPD","Del","Gly"),fill=c("red","green","blue"))

dat_norm <- betweenLaneNormalization(data,which="full")
boxplot(dat_norm, col=colors,ylab="log(count+1)",main="After FQ between-lane normalization",axes=F)
axis(2)
axis(1,lab=pData(data)$lib_prep,at=1:14)
dev.off()

ypds <- dat_norm[,pData(dat_norm)$conditions=="YPD"]
dim(ypds) # 5690 8; as expected
colors <- rep("black",8)
colors[pData(ypds)$lib_prep=="Y2"] <- "red"
colors[pData(ypds)$lib_prep=="Y4"] <- "green"
colors[pData(ypds)$lib_prep=="Y7"] <- "blue"
pdf("figure_1.pdf")
biasPlot(ypds,"gc",log=T,ylim=c(0,6),xlab="GC-content",ylab="log(count+1)",
         col=colors)
legend("topleft",c("Y1","Y2","Y4","Y7"),fill=c("black","red","green","blue"))
dev.off()


png("figure_2.png",width=960)
par(mfrow=c(1,2))
lfc <- log(exprs(dat_norm)[,2]+0.5)-log(exprs(dat_norm)[,1]+0.5)
biasBoxplot(lfc,fData(dat_norm)$gc,ylab="log-fold-change",
            main="Same Library Prep")

# YPD Y1 vs Y2 lane 428R1
lfc <- log(exprs(dat_norm)[,3]+0.5)-log(exprs(dat_norm)[,1]+0.5)
biasBoxplot(lfc,fData(dat_norm)$gc,ylab="log-fold-change",
            main="Different Library Preps")
dev.off()


dataWithin <- withinLaneNormalization(data,"gc",which="full")
dataNorm <- betweenLaneNormalization(dataWithin,which="full")
lfc <- log(exprs(dataNorm)[,3]+0.5)-log(exprs(dataNorm)[,1]+0.5)
pdf("figure_3c.pdf")
biasBoxplot(lfc,fData(dataNorm)$gc,ylab="log-fold-change",ylim=c(-2,2),
            main="FQ within-lane normalization\nFollowed by FQ between-lane normalization")
dev.off()

dataWithin <- withinLaneNormalization(data,"gc",which="loess")
dataNorm <- betweenLaneNormalization(dataWithin,which="full")
lfc <- log(exprs(dataNorm)[,3]+0.5)-log(exprs(dataNorm)[,1]+0.5)
pdf("figure_3a.pdf")
biasBoxplot(lfc,fData(dataNorm)$gc,ylab="log-fold-change",ylim=c(-2,2),
            main="Loess within-lane normalization\nFollowed by FQ between-lane normalization")
dev.off()

rm(list=ls())


###
# 13. Do pseudo-datasets, make plots

data <- get(load("geneExprSet_data_fromWebsite.RData"))
exprs(data) <- exprs(data)+1

dataWithin <- withinLaneNormalization(data,"gc",which="full",num.bins=50)
dataNorm <- betweenLaneNormalization(dataWithin,which="full")

YPD_lanes = row.names(pData(data)[which(pData(data)$conditions=="YPD"),])
groups = combn(8,4)[,1:35]

YPD_lanes
groups

grp0 = grp1 = matrix(0,nrow=dim(exprs(data))[1],ncol=35)
for (i in 1:35){
  grp0[,i] = apply(exprs(data)[,YPD_lanes[groups[,i]]],1,mean)
  grp1[,i] = apply(exprs(data)[,YPD_lanes[-groups[,i]]],1,mean)
  #grp0[,i] = exprs(dataii)[,YPD_lanes[groups[,i]]]
  #grp1[,i] = exprs(dataii)[,YPD_lanes[-groups[,i]]]
}

lfc <- log(grp1+1)-log(grp0+1)
pdf("figure_s11a.pdf")
boxplot(lfc,ylim=c(-2,2),outline=FALSE,xlab="pseudo-dataset",ylab="log-fold-change",
        main="Unnormalized")
abline(h=0,lty=2)
dev.off()

databw <- betweenLaneNormalization(data,which="full")
grp0bw = grp1bw = matrix(0,nrow=dim(exprs(data))[1],ncol=35)
for (i in 1:35){
  grp0bw[,i] = apply(exprs(databw)[,YPD_lanes[groups[,i]]],1,mean)
  grp1bw[,i] = apply(exprs(databw)[,YPD_lanes[-groups[,i]]],1,mean)
}

lfc = log2(grp1bw+0.5)-log2(grp0bw+0.5)
pdf("figure_s11b.pdf")
boxplot(lfc,ylim=c(-0.2,0.2),outline=FALSE,xlab="pseudo-dataset",ylab="log-fold-change",
        main="Only between-lane FQ normalization")
abline(h=0,lty=2)
dev.off()


grp0norm = grp1norm = matrix(0,nrow=dim(exprs(data))[1],ncol=35)
for (i in 1:35){
  grp0norm[,i] = apply(exprs(dataNorm)[,YPD_lanes[groups[,i]]],1,mean)
  grp1norm[,i] = apply(exprs(dataNorm)[,YPD_lanes[-groups[,i]]],1,mean)
}
lfc = log(grp1norm+0.5)-log(grp0norm+0.5)
pdf("figure_s11e.pdf")
boxplot(lfc,ylim=c(-0.2,0.2),outline=FALSE,xlab="pseudo-dataset",ylab="log-fold-change",
        main="FQ Normalization")
abline(h=0,lty=2)
dev.off()


lfc <- log(grp1norm+1)-log(grp0norm+1)

bias = apply(lfc,1,mean)
lfc2 = lfc^2
mse = apply(lfc2,1,mean)
boxplot(bias,outline=FALSE,ylim=c(-0.2,0.4),col="blue",xlab="FQ",ylab="bias")
abline(h=0,lty=2)

save(bias,file="FQ_bias.RData")
save(mse,file="FQ_mse.RData")

boxplot(mse,ylim=c(0,0.12),outline=FALSE,col="blue",xlab="FQ",ylab="MSE")

library(edgeR)
design <- model.matrix(~conditions, data=pData(dataNorm))
disp <- estimateGLMCommonDisp(exprs(dataNorm),design)#, offset=-offst(dataOffset))
fit <- glmFit(exprs(dataNorm), design, disp)#, offset=-offst(dataOffset))
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

nominalt1 <- seq(from=0,to=1,by=0.001)
actualt1 <- nominalt1
for(i in 1:length(nominalt1)){
  actualt1[i] <- sum(lrt$table$PValue<nominalt1[i])/nrow(lrt)
}


ps <- data.frame("pseudo"=c(rep("grp1",4),rep("grp0",4)))
design <- model.matrix(~pseudo,data=ps)

nominalt1 <- seq(from=0,to=1,by=0.01)
actualt1 <- matrix(0,nrow=length(nominalt1),ncol=35)

for (i in 1:35){
  dat <- cbind(exprs(dataNorm)[,YPD_lanes[groups[,i]]],exprs(dataNorm)[,YPD_lanes[-groups[,i]]])
  disp <- estimateGLMCommonDisp(dat,design)
  fit <- glmFit(dat,design,disp)
  lrt <- glmLRT(fit,coef=2)
  for(j in 1:length(nominalt1)){
    actualt1[j,i] <- sum(lrt$table$PValue<nominalt1[j])/nrow(lrt)
  }
  #plot(nominalt1,actualt1[,i]-nominalt1,type="l")
}

pdf("figure_s12d.pdf")
plot(nominalt1,actualt1[,1]-nominalt1,type="l",ylim=c(-0.4,0.4),
     xlab="nominal Type I error rate",ylab="actual Type I error rate - nominal Type I error rate",
     main="FQ")
for(i in 1:35){
  points(nominalt1,actualt1[,i]-nominalt1,type="l")
}
diffEr <- actualt1-nominalt1
meds <- apply(diffEr,1,median)
abline(h=0,lty=2)
points(nominalt1,meds,col="red",type="l",lwd=1.5)
dev.off()

rm(list=ls())


## Now for Loess Normalization pseudo-groups, plots

data <- get(load("geneExprSet_data_fromWebsite.RData"))
exprs(data) <- exprs(data)+1

# Loess Normalization for Fig S11c
dataWithin <- withinLaneNormalization(data,"gc",which="loess")
dataLoess <- betweenLaneNormalization(dataWithin,which="full")

YPD_lanes = row.names(pData(data)[which(pData(data)$conditions=="YPD"),])
groups = combn(8,4)[,1:35]

grp0loess = grp1loess = matrix(0,nrow=dim(exprs(data))[1],ncol=35)
for (i in 1:35){
  grp0loess[,i] = apply(exprs(dataLoess)[,YPD_lanes[groups[,i]]],1,mean)
  grp1loess[,i] = apply(exprs(dataLoess)[,YPD_lanes[-groups[,i]]],1,mean)
}
lfc = log(grp1loess+0.5)-log(grp0loess+0.5)
pdf("figure_s11c.pdf")
boxplot(lfc,ylim=c(-0.2,0.2),outline=FALSE,xlab="pseudo-dataset",ylab="log-fold-change",
        main="Loess Normalization")
abline(h=0,lty=2)
dev.off()


lfc <- log(grp1loess+1)-log(grp0loess+1)

bias = apply(lfc,1,mean)
lfc2 = lfc^2
mse = apply(lfc2,1,mean)
boxplot(bias,outline=FALSE,ylim=c(-0.2,0.4),col="red",xlab="Loess",ylab="bias")
abline(h=0,lty=2)

save(bias,file="loess_bias.RData") #for Fig S9
save(mse,file="loess_mse.RData") #for Fig S10

boxplot(mse,ylim=c(0,0.12),outline=FALSE,col="red",xlab="Loess",ylab="MSE")

library(edgeR)

design <- model.matrix(~conditions, data=pData(dataLoess))
disp <- estimateGLMCommonDisp(exprs(dataLoess),design)#, offset=-offst(dataOffset))
fit <- glmFit(exprs(dataLoess), design, disp)#, offset=-offst(dataOffset))
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

nominalt1 <- seq(from=0,to=1,by=0.001)
actualt1 <- nominalt1
for(i in 1:length(nominalt1)){
  actualt1[i] <- sum(lrt$table$PValue<nominalt1[i])/nrow(lrt)
}


ps <- data.frame("pseudo"=c(rep("grp1",4),rep("grp0",4)))
design <- model.matrix(~pseudo,data=ps)

nominalt1 <- seq(from=0,to=1,by=0.01)
actualt1 <- matrix(0,nrow=length(nominalt1),ncol=35)

for (i in 1:35){
  dat <- cbind(exprs(dataLoess)[,YPD_lanes[groups[,i]]],exprs(dataLoess)[,YPD_lanes[-groups[,i]]])
  disp <- estimateGLMCommonDisp(dat,design)
  fit <- glmFit(dat,design,disp)
  lrt <- glmLRT(fit,coef=2)
  for(j in 1:length(nominalt1)){
    actualt1[j,i] <- sum(lrt$table$PValue<nominalt1[j])/nrow(lrt)
  }
  #plot(nominalt1,actualt1[,i]-nominalt1,type="l")
}

pdf("figure_s12b.pdf")
plot(nominalt1,actualt1[,1]-nominalt1,type="l",ylim=c(-0.4,0.4),
     xlab="nominal Type I error rate",ylab="actual Type I error rate - nominal Type I error rate",
     main="Loess")
for(i in 1:35){
  points(nominalt1,actualt1[,i]-nominalt1,type="l")
}
diffEr <- actualt1-nominalt1
meds <- apply(diffEr,1,median)
abline(h=0,lty=2)
points(nominalt1,meds,col="red",type="l",lwd=1.5)
dev.off()

rm(list=ls())

# End Loess Normalization