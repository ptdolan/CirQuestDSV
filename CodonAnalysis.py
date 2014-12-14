def readNtFreqs():
    infile = "/Users/ptdolan/Desktop/AndinoLab/DengueData/SetA/Human/p1/Q20threshold.txt"
    freqTable = open(infile).readlines()
    freqTable=[line.strip("\n").split("\t") for line in freqTable]

    return freqTable

def translate(ntFreqTable, start, end): ##Converts nt counts to codon counts (3x4 numpy array), readDepth (3xN list), and sequence (3xN list).
    freqs=[pos[2:] for pos in ntFreqTable]
    ntSeq=[pos[1] for pos in ntFreqTable]
    codonSeq   = [ntSeq[x:x+3] for x in range(start, end, 3)]
    codonCounts = [freqs[x:x+3] for x in range(start, end, 3)] 
    x = np.array(codonCounts, dtype='|S4')
    codonCounts = x.astype(np.int)
    readDepth=[]
    for codonPos in range(0,len(codonCounts)):
        ntDepth=[sum(codonCounts[codonPos][pos]) for pos in range(0,3) ]
        readDepth.append(ntDepth)
    #print readDepth
    return codonSeq,codonCounts,readDepth
    
def codonSubCounts(codonSeq, codonCounts, readDepth): # Determine counts for specific codon substitution. 
    subCounts=np.zeros([64,64])
    subReadDepth=np.ones([64,64])
    nucList=["A","C","G","T"]
    codonList=[[i,j,k] for i in nucList for j in nucList for k in nucList ] #make alphabetical Codon Table
    asCodonIDs = [ codonList.index(codon) for codon in codonSeq ] # WT sequence as codon ID numbers
    for pos in range(0, len(codonCounts)): # Sum up all counts and Depths 
        wt=codonList.index(codonSeq[pos])
        scalar = {0:16,1:4,2:1}
        for site in range(0,3):
            for nuc in range(0,4):
                mutation = (scalar[site]*nuc-scalar[site]*nucList.index(codonSeq[pos][site]))
                mutID = wt+mutation
                subCounts[wt][mutID]+=codonCounts[pos][site][nuc]
                subReadDepth[wt][mutID] += readDepth[pos][site]
    return subCounts, subReadDepth

def bootstrapFreqs(counts,depth):
	pass #TO DO.
	
def meanSubFreqs(counts,depth):
    meanFreqs=counts/depth
    return meanFreqs

def output(meanFreqs):
    pass


start=96 
end=10272 
import sys, csv 
import rpy2
import numpy as np 

## MAIN 

ntFreqs=readNtFreqs() 
codonSeq,codonCounts,readDepth=translate(ntFreqs, start, end) 
subCountTable, subReadDepth=codonSubCounts(codonSeq, codonCounts, readDepth) 
bsFreqs = bootstrapFreqs(subCountTable, subReadDepth) 
meanFreqs = meanSubFreqs(subCountTable, subReadDepth)
print meanFreqs
