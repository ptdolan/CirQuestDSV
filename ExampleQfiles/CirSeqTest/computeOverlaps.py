import os
from operator import itemgetter
import numpy as np
#Read in TEs
#Read in Contigs

def faiDict(faiDir):
	print("Building sequence dictionary...")
	faList=[]
	for root,dirs,files in os.walk(faiDir):
		for F in files:
			if ".fa.fai" in F:
				faList.append(root+"/"+F)
	print (faList)
	for F in faList:
		with open(F,'r') as IF:
			faiFile= [line.strip().split("\t") for line in IF]
		lengths=[int(line[1]) for line in faiFile]
		convertedList=[[ line[j] if j==0 else int(line[j]) for j in range(len(line))] for line in faiFile]
		convertedList.sort(key=itemgetter(1),reverse=True)
		sortedLengths=np.array([int(line[1]) for line in convertedList])
		LenSum=np.cumsum(sortedLengths)
		adjusted=sortedLengths+LenSum
		print(adjusted)
	return(adjusted)

def bedRead(bed):
	print("read Bed: "+bed)
	outDict=dict()
	with open(bed,'r') as IF:
		bedFile= [line.strip().split("\t") for line in IF]
		intList=[list(range(int(line[1]),int(line[2]),1)) for line in bedFile]
		contigList=[line[0] for line in bedFile]
		dictInput=zip(contigList,intList)
		for line in dictInput:
			if line[0] in outDict.keys():
				currList=outDict[line[0]]
				currList.append(line[1])
				outDict.update({line[0]:currList})
			else:
				outDict.update({line[0]:line[1]})
	return(outDict)

def eveRead(eveFile):
	print("Building sequence dictionary...")
	faList=[]
	for root,dirs,files in os.walk(faiDir):
		for F in files:
			#print (F)
			if ".fa.fai" in F:
				faList.append(root+"/"+F)
		for F in faList:
			with open(F,'r') as IF:
				faiFile= [line.strip().split("\t") for line in IF]
				lengths=[int(line[1]) for line in faiFile]
				convertedList=[[ line[j] if j==0 else int(line[j]) for j in range(len(line))] for line in faiFile]
				convertedList.sort(key=itemgetter(1),reverse=True)
	return(convertedList)
def ComputeIntList(input):
	[range(start,end,1)]


faiDir= "/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/"
teFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs_TEs.bed"

eveFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs_EVEs_sorted.bed_withTaxonomy.txt"

positions=faiDict(faiDir)
positions=bedRead(teFile)
print(positions["000001F"])
