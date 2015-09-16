require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
source("~/compbio/projects/yuping/bin/bedEnrich_single.R")
library(foreach)
library(doMC)
registerDoMC(cores=10)

Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]

key="GraspGWAS"
#annotationBedFiles=Sys.glob("/home/unix/zzhang/hptmp/yp/CellSpecificEnh/*gz")
#
querybed="test/ceu.bed"


querybed=Args[2]
key=Args[3]
annotationBedFiles=Sys.glob(sprintf("/home/unix/zzhang/hptmp/yp/%s/*gz",key))

outdir="reports"
if(length(Args)>3){
        outdir=Args[4]
}


ChIPseeker_Enrich<-function(querybed,annotationBedFiles){
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tmpdir="~/htpmp/yp/tmp/"
system(sprintf("mkdirhier %s",tmpdir))
topPeakbed=sprintf("%s/top.%s",tmpdir,basename(querybed))
query.gr=readPeakFile(querybed)

topN=min(length(query.gr)*0.05,10000)
o=order(-query.gr$V4)
write.table(as.data.frame(query.gr[o[1:topN]]),file=topPeakbed,sep="\t",quote=F,row.names=F,col.names=F)

results=enrichPeakOverlap(queryPeak = topPeakbed,
targetPeak = unlist(annotationBedFiles),
TxDb = txdb,
pAdjustMethod = "BH",
nShuffle = 1000,
chainFile = NULL,
verbose = FALSE)

annotNames=results[,2]
frac=results[,5]/results[,4]
p=results[,6]

ret=cbind(results[,c("qLen","tLen","N_OL")],results[,7],p)
rownames(ret)=annotNames
colnames(ret)=c("QueryCount","AnnotationCount","OverlapAnnotCount","genomewide.shuffle.fdr","genomewide.shuffle.P")
ret
}

FlankingShuffle_Enrich<-function(querybed,annotationBedFiles){
	query.gr=readPeakFile(querybed)
	topN=min(length(query.gr)*0.05,10000)
	o=order(-query.gr$V4)
	top.gr=query.gr[o[1:topN]]
	results=foreach(fl=annotationBedFiles,.combine=rbind)%dopar%{
		annot.gr=readPeakFile(fl)	
		FlankingShuffleTest(top.gr,annot.gr,nperm=100)
	}	
	rownames(results)=unlist(lapply(annotationBedFiles,basename))
	colnames(results)=c("OverlapQueryCount","FlankingShuffleTest.P")
	results
}

KSTest_Enrich<-function(querybed,annotationBedFiles){
        query.gr=readPeakFile(querybed)
        results=foreach(fl=annotationBedFiles,.combine=rbind)%dopar%{
                annot.gr=readPeakFile(fl)
                KSTest(query.gr,annot.gr)
        }
        rownames(results)=unlist(lapply(annotationBedFiles,basename))
        colnames(results)=c("KSTest.stat","KSTest.P")
        results
}



PermutationTest_Enrich<-function(querybed,annotationBedFiles,confoundTable=NULL){
        query.gr=readPeakFile(querybed)
        results=foreach(fl=annotationBedFiles,.combine=rbind)%dopar%{
                annot.gr=readPeakFile(fl)
                PermutationTest_confoundTable(query.gr,annot.gr,confoundTable)
        }
        rownames(results)=unlist(lapply(annotationBedFiles,basename))
	colnames(results)=c("PermuteTest.stat","PermuteTest.P")
        results
}


AnnovaTest_Enrich<-function(querybed,annotationBedFiles,confoundTable=NULL){
        query.gr=readPeakFile(querybed)
        results=foreach(fl=annotationBedFiles,.combine=rbind)%dopar%{
                annot.gr=readPeakFile(fl)
                AnnovaTest_confoundTable(query.gr,annot.gr,confoundTable)
        }
        rownames(results)=unlist(lapply(annotationBedFiles,basename))
 	colnames(results)=c("AnnovaTest.stat","AnnovaTest.P")
        results
}



print(basename(querybed))
r1=ChIPseeker_Enrich(querybed,annotationBedFiles)
r5=FlankingShuffle_Enrich(querybed,annotationBedFiles)
r2=KSTest_Enrich(querybed,annotationBedFiles)
r3=PermutationTest_Enrich(querybed,annotationBedFiles)
r4=AnnovaTest_Enrich(querybed,annotationBedFiles)

o=match(rownames(r2),rownames(r1))
combine=cbind(r1[o,],r5,r2,r3,r4)
write.table(combine,file=sprintf("%s/%s_%s.report",outdir,basename(querybed),key),sep="\t",quote=F)



