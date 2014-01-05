
## ----setup, include=FALSE, cache=FALSE-----------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold')
options(replace.assign=FALSE,width=90)


## ----autoinstall, eval=FALSE---------------------------------------------
## setRepositories(ind=1:6)
## install.packages("CAGExploreR")


## ----loading-------------------------------------------------------------
library(CAGExploreR)


## ----getinternalEnsemblAnot,eval=TRUE------------------------------------
data(EnsemblAnot_chr22) #loads an object EnsemblAnot
EnsemblAnot = prepEnsemblAnot(EnsemblAnot)
#You now skip to the next section


## ----getEnsemblAnot,eval=FALSE-------------------------------------------
## EnsemblAnot = getEnsemblAnot()


## ----prepEnsemblAnot,eval=FALSE------------------------------------------
## EnsemblAnot = prepEnsemblAnot(EnsemblAnot)


## ----saveEnsemblAnot,eval=FALSE------------------------------------------
## save(EnsemblAnot,file="EnsemblAnot.RData")
## #can append today's date to the file for version tracking:
## #save(EnsemblAnot,file=paste0("EnsemblAnot",Sys.Date(),".RData"))


## ----loadEnsemblAnot,eval=FALSE------------------------------------------
## load("EnsemblAnot.RData")
## #replace name by the one you chose


## ----files, eval=FALSE---------------------------------------------------
## setwd("D:/bam") #set the working directory to where the bam and bai files are located
## files = dir()[grep("bam",dir())]
## my.bai.files = dir()[grep("bam.bai",dir())]
## my.bam.files = setdiff(files, my.bai.files)
## my.ids = c("a549.1","a549.2","mcf7.1","mcf7.2")


## ----filesB, eval=TRUE---------------------------------------------------
my.bai.files = c("wgEncodeRikenCageA549CellPapAlnRep1.bam.bai"

,"wgEncodeRikenCageA549CellPapAlnRep2.bam.bai"

,"wgEncodeRikenCageMcf7CellPapAlnRep1.bam.bai"

,"wgEncodeRikenCageMcf7CellPapAlnRep2.bam.bai")

my.bam.files = c("wgEncodeRikenCageA549CellPapAlnRep1.bam"

,"wgEncodeRikenCageA549CellPapAlnRep2.bam"

,"wgEncodeRikenCageMcf7CellPapAlnRep1.bam"

,"wgEncodeRikenCageMcf7CellPapAlnRep2.bam")

my.ids = c("a549.1","a549.2","mcf7.1","mcf7.2")


## ----promdefs, eval=TRUE, cache=TRUE,echo=1:2,comment=NA-----------------
data(F5.hg19.promoters) #this object contains both promoter and gene regions
head(F5.hg19.promoters)
unique.genes = unique(sapply(strsplit(as.character(F5.hg19.promoters$gene),"@"),function(x) x[1]))
luq = length(unique.genes)
n.full.gene.regions = length(grep("@",as.character(F5.hg19.promoters$gene)))
prom.only = F5.hg19.promoters[-grep("@",as.character(F5.hg19.promoters$gene)),]
n.promoters = nrow(prom.only)
av.width = mean(prom.only$end - prom.only$start)
sd.width = sd(prom.only$end - prom.only$start)


## ----lookpromdefs2, eval=TRUE, echo=FALSE,cache=TRUE---------------------
data(mpromdb.hg19.promoters)
unique.genes2 = unique(sapply(strsplit(as.character(mpromdb.hg19.promoters$gene),"@"),function(x) x[1]))
luq2 = length(unique.genes2)
n.full.gene.regions2 = length(grep("@",as.character(mpromdb.hg19.promoters$gene)))
prom.only2 = mpromdb.hg19.promoters[-grep("@",as.character(mpromdb.hg19.promoters$gene)),]
n.promoters2 = nrow(prom.only2)
av.width2 = mean(prom.only2$end - prom.only2$start)
sd.width2 = sd(prom.only2$end - prom.only2$start)


## ----counttags,eval=FALSE------------------------------------------------
## my.promoters = definePromoters(F5.hg19.promoters)
## my.genes = defineGenes(F5.hg19.promoters)
## 
## mcf7a549.raw.counts.F5 = countTags(my.bam.files, my.bai.files, my.ids, my.genes, my.promoters)
## #This takes 2-10 minutes per BAM file


## ----loaddata,comment=NA,cache=TRUE--------------------------------------
data(mcf7a549.raw.counts.F5)
mcf7a549.raw.counts.F5 #use of the head() optional here


## ----subsetdata22,comment=NA---------------------------------------------
get = which(osc2info(mcf7a549.raw.counts.F5$counts$region)$chr=="chr22")
mcf7a549.raw.counts.F5$counts = mcf7a549.raw.counts.F5$counts[get]


## ----dgelist-------------------------------------------------------------
#osc2info converts genomic regions in the format "chr:start..end,strand" to an R list
annotations = osc2info(mcf7a549.raw.counts.F5$counts$region)
my.ids = c("a549.1","a549.2","mcf7.1","mcf7.2")

my.data = DGEList(

	counts = mcf7a549.raw.counts.F5$counts[,-c(1,2),with=FALSE], #remove columns 1&2
	lib.size = mcf7a549.raw.counts.F5$depth,
	group = my.ids,                          #will not pool replicates
	#group = c("a549","a549","mcf7","mcf7"), #pools replicates
	#group = Select(my.ids,".",1),           #pools replicates (alternative)
	
	genes = data.frame(
		chr = annotations$chr,
		strand = annotations$strand,
		start = annotations$start,
		end = annotations$end,
		gene = mcf7a549.raw.counts.F5$counts$gene
		),
		
	remove.zeros = FALSE
)


## ----dgelistview,comment=NA----------------------------------------------
my.data


## ----edgeroptional,eval=FALSE--------------------------------------------
## #To normalize samples
## my.data = calcNormFactors(my.data)
## #To estimate negative binomial dispersion (common or tagwise)
## my.data = estimateCommonDisp(my.data,verbose=T)
## my.data = estimateTagwiseDisp(my.data)


## ----pool,comment=NA-----------------------------------------------------
#This step is REQUIRED even if not pooling replicates
data.not.pooled = pool(my.data)
data.not.pooled


## ----results,comment=NA--------------------------------------------------
results = diffcomp(data.not.pooled)
head(results)


## ----subsetting,comment=NA-----------------------------------------------
significant = subset(results,fdr<0.001 & coverage >= 0.1 & coverage <= 1 & RepAgree==1)
head(significant)
#This is how many significant, high-coverage genes we are left with:
nrow(significant)
#This is the number of genes in which the dominant promoter switches between conditions:
sum(significant$dominant.promoter.switch != "")
#This is the number of genes that have differential promoter composition but very little to no differential gene expression:
sum(significant$geneHetero<0.01)
#Differences in gene expression and differences in promoter composition do not appear to be correlated:
cor(significant$entropy.Reduction,significant$geneHetero)


## ----results.detailed,comment=NA-----------------------------------------
my.data2 = DGEList(

	counts = mcf7a549.raw.counts.F5$counts[,-c(1,2),with=FALSE], #remove columns 1&2
	lib.size = mcf7a549.raw.counts.F5$depth,
	#group = my.ids,                          #will not pool replicates
	group = c("a549","a549","mcf7","mcf7"), #pools replicates
	#group = Select(my.ids,".",1),           #pools replicates (alternative)
	
	genes = data.frame(
		chr = annotations$chr,
		strand = annotations$strand,
		start = annotations$start,
		end = annotations$end,
		gene = mcf7a549.raw.counts.F5$counts$gene
		),
		
	remove.zeros = FALSE
)

data.pooled = pool(my.data2)

results.detailed = diffcomp(data.pooled,detailed=TRUE)
head(results.detailed,15)


## ----subsetting2,comment=NA----------------------------------------------
significant.genes = rownames(subset(results,fdr<0.001 & coverage >= 0.5 & coverage <= 1  & RepAgree==1 & dominant.promoter.switch != ""))
sig.results.detailed = subset(results.detailed,gene %in% significant.genes)

#We are down to just this number of pair-wise promoter comparisons:
nrow(sig.results.detailed)
head(sig.results.detailed)


## ----plot.detailed,echo=1,comment=NA,fig.cap="Differential promoter composition for GAL3ST1 gene with 13 promoters"----
plotcomp(data.not.pooled, "GAL3ST1",EnsemblAnot)
dev.off()


## ----plot.html,eval=FALSE------------------------------------------------
## html.report(data.not.pooled,rownames(significant)[1:30],EnsemblAnot)


## ----plot.html2,eval=FALSE-----------------------------------------------
## html.report(data.not.pooled,my.genes,EnsemblAnot,fig.dir="myFigures",report.name="my.report")


## ----genomewideplot,fig.cap="Genome-wide volcano plot for comparing differential promoter composition between 2 cell lines",out.width="4in",fig.show="asis",fig.pos="!hb"----
plot(y=-log(results.detailed$fdr,base=10),x=results.detailed$log2OR,
xlim=c(-10,10),ylim=c(0,50),col=rgb(0,100,0,60,maxColorValue=255),

pch=20,ylab="-log10(q-value)",xlab="log2(odds ratio)",
main="MCF7 v.s. A549 cell lines")

abline(v=0,h=0)


## ----sessioninfo---------------------------------------------------------
sessionInfo()
date()


