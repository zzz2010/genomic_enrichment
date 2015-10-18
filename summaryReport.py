import os,sys
import glob
import numpy as np


Fls=glob.glob(sys.argv[1]+"/*.report")

def getPvalue(line):
	startcol=5  #0-based
	comps=line.strip().split()
	nscore=(len(comps)-startcol)/2
	scorelist=list()
	for i in range(nscore):
		p=float(comps[i*2+startcol]) ##need to change 3 to 4 later
		scorelist.append(p)
	return float(sorted(scorelist)[nscore-1])

def getPopulation(bname):
	return bname.split(".")[1]

def getTerm(line):
	return line.split()[0].split(".")[0]
term_pops=dict()
for fl in Fls:
	bname=os.path.basename(fl)
	pop=getPopulation(bname)
	lines=open(fl).readlines()
	for line in lines[1:]:
		p=getPvalue(line)
		if p<0.1:
		#	print p,bname,line
			t=getTerm(line)
			if t not in term_pops:
				term_pops[t]=list()
			term_pops[t].append(pop)

for t in term_pops:
	print t, ",".join(term_pops[t])


Fls2=glob.glob(sys.argv[1]+"/*.GO")
go_pops=dict()
for fl in Fls2:
	bname=os.path.basename(fl)
	pop=getPopulation(bname)
	mode=0
	for line in open(fl):
		if line.startswith("Geneset.Type"):
			continue
		if line.startswith("category"):
			mode+=1
			continue
		if line.startswith("GO.ID"):
			mode+=1
			continue
		comps=line.strip().split("\t")
		if mode==0:
			term=comps[3]
		if mode==1:
			term=comps[6]
			continue
		if mode==2:
			term=comps[2]
			continue
		if term not in go_pops:
			go_pops[term]=set()
		go_pops[term].add(pop)

for term in go_pops:
	pops=go_pops[term]
	if len(pops)<2:
		continue
	print "GO|"+term+"|"+",".join(pops)
