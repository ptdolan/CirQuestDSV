# # # CodonAnalysis.py
#
# Rewritten python version of CirSeq Analysis Code. Converts Q-file to Codon Counts (Assumes point mutations). Will include bin-test, F.E.T., Bootstrapping, HMM(?).
# Author: P.T. Dolan (ptdolan@stanford.edu)
# Usage: python CodonAnalysis.py
#
# # #

def makeCodonList(nucList):
    codonList = [[i,j,k] for i in nucList for j in nucList for k in nucList ] # make alphabetical Codon Table
    return codonList

def readNtFreqs():
    infile = "/Users/ptdolan/Research/CirSeq/Dengue/Shuhei_ExperimentalData/SetA/Human/p1/Q20threshold.txt"
    freqTable = open(infile).readlines()
    freqTable=[line.strip("\n").split("\t") for line in freqTable]
    #print freqTable
    return freqTable

def translate(ntFreqTable, start, end): ##Converts nt counts to codon counts (3x4 numpy array), readDepth (3xN list), and sequence (3xN list). DONE
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
        #print codonSeq
        #print readDepth
    return codonSeq,codonCounts,readDepth

def codonSubCounts(codonSeq, codonCounts, readDepth, nucList, codonList): # Determine counts for specific codon substitution.  DONE
    subCounts=np.zeros([64,64])
    #depthList=[]
    #countList=[]
    subReadDepth=np.ones([64,64])                         # !!! pseudoCounts
    asCodonIDs = [ codonList.index(codon) for codon in codonSeq ] # WT sequence as codon ID numbers
    for pos in range(0, len(codonCounts)):                # Sum up all counts and Depths... For each position
        wt=codonList.index(codonSeq[pos])                 #For each site in codon.
        scalar = {0:16,1:4,2:1}                           #Adjusts movement on 'subCounts' table. 16 for first site, 4 for second, 1 for third.
        for site in range(0,3):
        
            for nuc in range(0,4):                        #for each Mutation
                mutation = (scalar[site]*nuc-scalar[site]*nucList.index(codonSeq[pos][site]))
                mutID = wt+mutation
                subCounts[wt][mutID] += codonCounts[pos][site][nuc]         #Add reads to cumulative count for each transition
                subReadDepth[wt][mutID] += readDepth[pos][site]             #Add reads to read depth for each transition
    return subCounts, subReadDepth

def meanSubFreqs(counts,depth): # Simple Calculation of mean Substitution Frequencies. DONE, BUT MAY BE CONVERTED TO HMM...
    meanFreqs = counts / depth
    return meanFreqs

def bootstrapFreqs(codonSeq, codonCounts, readDepth, codonList):  # bootstrap counts or deviation from MLE? --- NEED TO THINK ABOUT BEST APPROACH HERE
    bootedFreqs=[]
    bootedReads=[]
    bootMeanCounts=[]
    bootDepthCounts=[]
    BSmeans =[]
    scalar = {0:16,1:4,2:1}
    print "Bootstrapping Frequencies"
    for codon in codonList:
        countList = []
        depthList = []
        BSdepth = []
        BSmeans = []
        posList = [] #for each Codon
        posDList = []
        for pos in range(0, len(codonCounts)):    # for each position
            siteList = []                         # make new list for sites counts
            siteDList = []                        # make new list for sites readDepth
            if codonSeq[pos] == codon:            # collect reads for each codon
                for site in range(0,3):           # for each site
                    nucList =  [codonCounts[pos][site][nuc] for nuc in range(0,4)]          # for each Mutation makes list of 4 elements
                    nucDList = [readDepth[pos][site] for nuc in range (0,4)]                # and list of depths (all the same for a given site)
                    siteList.append(nucList)                                                #
                    siteDList.append(nucDList)
            if codonSeq[pos] == codon:
                posList.append(siteList)
                posDList.append(siteDList)
            countList.append(posList)
            depthList.append(posDList)
        print codon
        print len(countList), countList[0][0]
    
        
        #bootstrapping
#    
#    nboot = 1                                 # Define number of boots
#    
#    for nuc in range(0,4):
#        for count in range( 0, len(countList)):
#            bootMeanCounts=[countList[ random.randint(0, len(countList[0]) -1)] for repl in range(0,nboot)]
#            print count
#        BSmeans.append( np.mean (bootMeanCounts) )
#        BSdepth.append( np.mean (bootDepthCounts) )
#    bootedFreqs.append(BSmeans)
#    bootedReads.append(BSdepth)
#
    return bootedFreqs

def binTest(count,depth,mean): #If binomial is preferred approach for characterizing counts, include this. --- NEEDS WORK.
    np.binom(20,100,0.5)

# # # # # #// Globals and Imports //# # # # # #

start=96            #Start of Coding Sequence
end=10272           #End of Coding Sequence

outfile = "/Users/ptdolan/GitHub/CirQuestDSV/codonCounts.txt"
nucList=["A","C","G","T"]

import sys, csv, random
#import rpy2
import numpy as np

# # # # # #//       MAIN      //# # # # # #

codonList = makeCodonList
ntFreqs   = readNtFreqs()
codonSeq,codonCounts,readDepth = translate(ntFreqs, start, end)

subCountTable, subReadDepth = codonSubCounts(codonSeq, codonCounts, readDepth, nucList, codonList)
print subCountTable
print subReadDepth

#bsFreqs = bootstrapFreqs(codonSeq, codonCounts, readDepth, codonList)
#print "BootedMeans: ", bsFreqs

meanFreqs = meanSubFreqs(subCountTable, subReadDepth)

np.savetxt(outfile, meanFreqs, delimiter="\t")