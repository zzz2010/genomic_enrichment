library(foreach)
require(ChIPseeker)
require(GenomicRanges)
#library(digest)
library(doMC)
registerDoMC(cores=10)


PermutationTest_confoundTable<-function(query.gr,annot.gr,confoundTable=NULL){
topN=min(length(query.gr)*0.05,10000)
o=order(-query.gr$V4)
topIndex=o[1:topN]

#query.gr$V4[unique(findOverlaps(query.gr,annot.gr)@queryHits)]=2000  ##debug purpose
observed.count=length(subsetByOverlaps(query.gr[topIndex],annot.gr))

###sth we need to control during the resampling###
confoundMat=data.frame(length=log(query.gr@ranges@width))

if(!is.null(confoundTable)){
confoundMat=cbind(confoundMat,confoundTable)
}

confoundMat.quantile=apply(round(apply(confoundMat,2,rank)/nrow(confoundMat),1),1,function(x) sum(x*10^(1:length(x))) )

lookupTable=aggregate(1:length(confoundMat.quantile),by=list(confoundMat.quantile),FUN=c)

topKey.count=table(confoundMat.quantile[topIndex])

nperm=1000
permCountV= rep(0,nperm);
for(ii in 1:length(topKey.count)){
	key=names(topKey.count)[ii]
	groupid=match(key,lookupTable[[1]])[1]
	sel.gr=query.gr[unlist(lookupTable[groupid,2])]
	tot.count=length(sel.gr)
	annot.count=length(subsetByOverlaps(sel.gr,annot.gr))
	if(annot.count>0){
	permCountV=permCountV+rhyper(nperm,annot.count,tot.count-annot.count,topKey.count[ii])
	}
}

p=sum(permCountV>=observed.count)/nperm
z=(observed.count-mean(permCountV))/(sd(permCountV)+1E-4)
c(z,p)
}


KSTest<-function(query.gr,annot.gr ){
	ol=findOverlaps(query.gr,annot.gr)	
	hits=unique(ol@queryHits)
	annot.scores=query.gr$V4[hits]
	other.scores=query.gr$V4[-hits]
	a=ks.test(annot.scores,other.scores,alternative="less") ##seesm less is the right direction
	c(a$statistic,a$p.value)
}

FlankingShuffleTest<-function(query.gr,annot.gr,nperm=1000,expand=4){
	ol=findOverlaps(query.gr,annot.gr)
        hits=unique(ol@queryHits)
	observed.count=length(hits)
	perm_count=foreach(i=1:nperm)%dopar%{
		delta=as.integer(expand*(runif(length(query.gr))-0.5)*query.gr@ranges@width)
		rand.gr=GRanges(seqnames=query.gr@seqnames,ranges=IRanges(query.gr@ranges@start+delta+1,query.gr@ranges@start+delta+query.gr@ranges@width ) )
		ol=findOverlaps(rand.gr,annot.gr)
	        hits=unique(ol@queryHits)
		length(hits)
	}
	pvalue=mean(observed.count<perm_count)
	c(observed.count,pvalue)
}

AnnovaTest_confoundTable<-function(query.gr,annot.gr,confoundTable=NULL){
ol=findOverlaps(query.gr,annot.gr)
hits=unique(ol@queryHits)

x=log(table(c(1:length(query.gr),ol@queryHits)))
###sth we need to control during the resampling###
confoundMat=data.frame(length=log(query.gr@ranges@width))

if(!is.null(confoundTable)){
confoundMat=cbind(confoundMat,confoundTable)
}
	covar=as.matrix(confoundMat)
	#x=(1:length(query.gr))%in%hits
	y=query.gr$V4
	fit0=lm(y~1+covar)
	fit1 <- lm(y~ covar+1 + x)
	a=anova(fit0, fit1)
	p=a[2,6]
	F=a[2,5]
	##one sided p
	if(fit1$coefficients["x"]>0){
		p=p/2
	}else{
		p=1-p/2
	}	
	c(F,p)
}

test<-function(){
querybed="test/e124.dhs.bed"
annotBed="/home/unix/zzhang/hptmp/yp/CellSpecificEnh/E124_25_imputed12marks_dense.bed.gz"

annot.gr=readPeakFile(annotBed)

query.gr=readPeakFile(querybed)

confoundTable=NULL #matrix(runif(length(query.gr)*4),nrow=length(query.gr))
FlankingShuffleTest(query.gr,annot.gr )
AnnovaTest_confoundTable(query.gr,annot.gr,confoundTable)
KSTest(query.gr,annot.gr)
PermutationTest_confoundTable(query.gr,annot.gr,confoundTable )

}


