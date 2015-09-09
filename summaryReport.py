import os,sys
import glob
import numpy as np


Fls=glob.glob("reports/*.report")

def getPvalue(line):
	comps=line.strip().split()
	nscore=(len(comps)-1)/2
	scorelist=list()
	for i in range(nscore):
		p=float(comps[i*2+2])
		scorelist.append(p)
	return float(sorted(scorelist)[1])

for fl in Fls:
	bname=os.path.basename(fl)
	lines=open(fl).readlines()
	for line in lines[1:]:
		p=getPvalue(line)
		if p<0.05/200:
			print p,bname,line
