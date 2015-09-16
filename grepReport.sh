echo  "Population	GWAS	QueryCount      AnnotationCount OverlapAnnotCount       genomewide.shuffle.fdr  genomewide.shuffle.P    OverlapQueryCount       FlankingShuffleTest.P   KSTest.stat     KSTest.P        PermuteTest.stat        PermuteTest.P   AnnovaTest.stat AnnovaTest.P"
grep $1 reports/*GWAS.report |sort -k8|sed 's/.bed_GraspGWAS.report//'|sed 's/:/\t/'|sed 's/.rsid.bed.gz//'
grep $1 reports_10kbp/*GWAS.report |sort -k8|sed 's/.bed_GraspGWAS.report//'|sed 's/:/\t/'|sed 's/.rsid.bed.gz//'
