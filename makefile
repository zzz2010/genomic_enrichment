~/hptmp/yp/totalPeaks/%.bed:~/compbio/projects/yuping/data/totalPeaks/%.H12_peaks.txt.2
	awk '{OFS="\t";print $$1,$$3,$$4,$$5}' $< > $@

~/hptmp/yp/totalPeaks_10kbp/%.bed:~/compbio/projects/yuping/data/totalPeaks/%.H12_peaks.txt.2
	awk '{OFS="\t";m=int(($$3+$$4)/2);print $$1,m-5000,m+5000,$$5}' $< > $@

all.totalPeaks_10kbp:
	mkdirhier ~/hptmp/yp/totalPeaks_10kbp
	ls ~/compbio/projects/yuping/data/totalPeaks/total.*.H12_peaks.txt.2|sed 's/.H12_peaks.txt.2//'|xargs -n 1 basename |awk -v d="~/hptmp/yp/totalPeaks_10kbp/" '{print d"/"$$1".bed"}'|xargs -n1  make 

all.totalPeaks:
	mkdirhier ~/hptmp/yp/totalPeaks
	ls ~/compbio/projects/yuping/data/totalPeaks/total.*.H12_peaks.txt.2|sed 's/.H12_peaks.txt.2//'|xargs -n 1 basename |awk -v d="~/hptmp/yp/totalPeaks/" '{print d"/"$$1".bed"}'|xargs -n1  make 

all.report.%:bedEnrich.R
	ls ~/hptmp/yp/totalPeaks/*.bed| xargs -n 1 -I {} Rscript $< {} $*

all.metaPop.report:
	ls /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/output/*.bed|xargs -n 1 -I {} Rscript bedEnrich.R  {} CellSpecificEnh /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/reports &
	ls /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/output/*.bed|xargs -n 1 -I {} Rscript bedEnrich.R  {} GraspGWAS /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/reports &
	ls /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/output/*.bed|xargs -n 1 -I {} Rscript bedEnrich.R  {} H3K27ac /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/reports &
	ls /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/output/*.bed|xargs -n 1 -I {} Rscript geneEnrich.R {} /broad/compbio/zhizhuo/projects/yuping/result/2015-09-10/region_matrix/reports 

all.GO:geneEnrich.R
	ls ~/hptmp/yp/totalPeaks/*.bed| xargs -n 1 -I {} Rscript $< {} reports

all.GO.10k:geneEnrich.R
	ls ~/hptmp/yp/totalPeaks/*.bed| xargs -n 1 -I {} Rscript $< {} reports_10kbp
	

all.10kbp.report.%:bedEnrich.R
	ls ~/hptmp/yp/totalPeaks_10kbp/*.bed| xargs -n 1 -I {} Rscript $< {} $*  reports_10kbp



