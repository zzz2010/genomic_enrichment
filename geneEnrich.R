require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(goseq)
library(topGO)
library(chipenrich)
require(ChIPseeker)
source("~/compbio/common_functions.R")


Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]

#querybed="test/ceu.bed"

querybed=Args[2]

outdir="reports"
if(length(Args)>2){
	outdir=Args[3]
}

tmpdir="~/htpmp/yp/tmp/"
system(sprintf("mkdirhier %s",tmpdir))
topPeakbed=sprintf("%s/top.%s",tmpdir,basename(querybed))
query.gr=readPeakFile(querybed)

topN=min(length(query.gr)*0.05,10000)
o=order(-query.gr$V4)
write.table(as.data.frame(query.gr[o[1:topN]]),file=topPeakbed,sep="\t",quote=F,row.names=F,col.names=F)


broadenrich<-function(topPeakbed){
results = chipenrich(peaks = topPeakbed, genesets =c("GOBP","GOCC","GOMF","kegg_pathway","ehmn_pathway_gene","gene_expression","biocarta_pathway","metabolite"), locusdef = "nearest_tss", max_geneset_size = 500, qc_plots = F,genome="hg19",use_mappability=F,method='broadenrich',n_cores=12)

topResults=results$results[results$results$FDR<0.05,]
topResults=topResults[order(topResults$FDR),]
topResults
}


peakAnno <- as.data.frame(annotatePeak(querybed, level="gene",tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db"))

genescores=aggregate(peakAnno$V4,by=list(peakAnno$ENSEMBL),FUN=max)


goseq_enrich<-function(genescores){
genes=genescores[,2]>(mean(genescores[,2])+2*sd(genescores[,2]))
names(genes)=genescores[,1]

pwf=nullp(genes,"hg19","ensGene" )
GO.wall=goseq(pwf,"hg19","ensGene",test.cats=c("GO:CC", "GO:BP", "GO:MF"))

#en2eg=as.list(org.Hs.egENSEMBL2EG)
#eg2kegg=as.list(org.Hs.egPATH)
#grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
#kegg=lapply(en2eg,grepKEGG,eg2kegg)
#KEGG=goseq(pwf,gene2cat=kegg)

topResults=GO.wall[GO.wall$over_represented_pvalue<0.01,]
topResults

}


topGO_enrich<-function(genescores){
FDR=1-genescores[,2]
names(FDR)<-genescores[,1]
p=mean(FDR)-2*sd(FDR)
BP=topGO(FDR,ontology="BP",p=p)
CC=topGO(FDR,ontology="CC",p=p)
MF=topGO(FDR,ontology="MF",p=p)
topResults=rbind(BP,CC,MF)
topResults[topResults$elimKS<0.01,]

}

outfn=sprintf("%s/%s.GO",outdir,basename(querybed))
t1=broadenrich(topPeakbed)
write.table(t1,file=outfn,sep="\t",quote=F,append = F)

t2=goseq_enrich(genescores)
write.table(t2,file=outfn,sep="\t",quote=F,append = T)

t3=topGO_enrich(genescores)
write.table(t3,file=outfn,sep="\t",quote=F,append = T)



