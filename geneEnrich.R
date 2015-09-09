library(topGO)
library(goseq)
library(chipenrich)
require(ChIPseeker)
require(clusterProfiler)


Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]

querybed="test/ceu.bed"

#querybed=Args[2]

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tmpdir="~/htpmp/yp/tmp/"
system(sprintf("mkdirhier %s",tmpdir))
topPeakbed=sprintf("%s/top.%s",tmpdir,basename(querybed))
query.gr=readPeakFile(querybed)

topN=min(length(query.gr)*0.05,10000)
o=order(-query.gr$V4)
write.table(as.data.frame(query.gr[o[1:topN]]),file=topPeakbed,sep="\t",quote=F,row.names=F,col.names=F)


results = chipenrich(peaks = topPeakbed, genesets =c("GOBP","GOCC","GOMF","kegg_pathway","ehmn_pathway_gene","gene_expression","biocarta_pathway","metabolite"), locusdef = "nearest_tss", max_geneset_size = 100, qc_plots = F,genome="hg19",out_name =sprintf("%s.chipenrich",basename(Args[2])),read_length=500,n_cores=12)

