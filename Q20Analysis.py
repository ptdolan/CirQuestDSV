#Q20Analysis.py
#	Purpose: convert all Q20 files (all files ending in "Q20.txt") in an entire directory tree to a codon annotated version. 
#	Also annotates amino acids with 18 residue properties.
#	Usage: python Q20Analysis.py <q20file> <translationStart> <translationStop>


import csv, os, sys, argparse, numpy as np

bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

ntDict={
'A':0,
'C':1,
'G':2,
'T':3}

resDict={
'Ala':0,  'A':0,
'Arg':1,  'R':1,
'Asn':2,  'N':2, 
'Asp':3,  'D':3,
'Cys':4,  'C':4, 
'Glu':5,  'E':5, 
'Gln':6,  'Q':6,
'Gly':7,  'G':7, 
'His':8,  'H':8, 
'Ile':9,  'I':9 , 
'Leu':10, 'L':10,
'Lys':11, 'K':11,
'Met':12, 'M':12, 
'Phe':13, 'F':13, 
'Pro':14, 'P':14 , 
'Ser':15, 'S':15, 
'Thr':16, 'T':16, 
'Trp':17, 'W':17, 
'Tyr':18, 'Y':18, 
'Val':19, 'V':19,
'Stop':20,'*':20,
"UTR":21, "U":21}

termDict={
0:'acidic',
1:'acyclic',
2:'aliphatic',
3:'aromatic',
4:'basic',
5:'buried',
6:'charged',
7:'cyclic',
8:'hydrophobic',
9:'large',
10:'medium',
11:'negative',
12:'neutral',
13:'polar',
14:'positive',
15:'small',
16:'surface',
17:'stop'}

resAnnotMatrix=[##Binary Vector -- could be mproved to be continuous data, see K. Chou's stuff.
#A  R  N  D  C  E  Q  G  H  I  L  K  M  F  P  S  T  W  Y  V  *
[0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],# acidic
[1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0],# acyclic
[1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],# aliphatic
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0],# aromatic
[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],# basic
[1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0],# buried
[0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],# charged
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],# cyclic
[1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],# hydrophobic
[0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0],# large
[0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0],# medium
[0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],# negative
[1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0],# neutral
[0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0],# polar
[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],# positive
[1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],# small
[0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0],# surface
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]# stop

pechmann=[["S","NS","C","C","NS","C","NS","NS","NS","NS","NS","NS","NC","NS","NS","C","C","NC","NS","NS","X","U"],
["NS","S","NS","NS","NC","C","NS","NS","NS","NS","NS","NS","NS","NS","C","C","NS","NS","C","C","X","U"],
["C","NS","S","C","NS","C","C","NS","NS","NS","NS","C","NS","NS","NS","NS","NS","C","NS","C","X","U"],
["C","NS","C","S","NS","C","NS","NS","C","NS","NS","NS","NS","C","NS","NS","NS","C","NS","NS","X","U"],
["NS","NC","NS","NS","S","NS","NS","C","NS","C","NS","NS","NS","NS","NS","NC","NS","NC","NS","C","X","U"],
["C","C","C","C","NS","S","NS","NS","NS","NS","NS","NS","NS","NS","C","NC","NS","NC","NC","NS","X","U"],
["NS","NS","C","NS","NS","NS","S","NS","NS","C","NS","C","NS","NS","C","C","C","NS","NS","C","X","U"],
["NS","NS","NS","NS","C","NS","NS","S","NC","C","C","C","NS","NS","NC","NC","NC","C","NS","NS","X","U"],
["NS","NS","NS","C","NS","NS","NS","NC","S","NS","C","C","NS","C","C","NS","C","NS","NS","NS","X","U"],
["NS","NS","NS","NS","C","NS","C","C","NS","S","C","NS","C","C","C","NC","NS","C","C","NS","X","U"],
["NS","NS","NS","NS","NS","NS","NS","C","C","C","S","NS","NS","NS","C","NS","C","C","NS","NS","X","U"],
["NS","NS","C","NS","NS","NS","C","C","C","NS","NS","S","NS","NS","NS","C","C","NS","NS","C","X","U"],
["NC","NS","NS","NS","NS","NS","C","NS","NS","C","NS","NS","S","C","C","C","NC","NS","NS","NS","X","U"],
["NS","NS","NS","C","NS","NS","C","NS","C","C","NS","NS","C","S","C","NS","NS","NS","NS","NS","X","U"],
["NS","C","NS","NS","NS","C","C","NC","C","C","C","NS","C","C","S","NC","C","NS","C","NS","X","U"],
["C","C","NS","NS","NC","NC","NS","NC","NS","NC","NS","C","C","NS","NC","S","C","NS","NS","C","X","U"],
["C","NS","NS","NS","NS","NS","NS","NC","C","NS","C","C","NC","NS","C","C","S","NS","NS","NS","X","U"],
["NC","NS","C","C","NC","NC","NS","C","NS","C","C","NS","NS","NS","NS","NS","NS","S","NS","NS","X","U"],
["NS","C","NS","NS","NS","NC","NS","NS","NS","C","NS","NS","NS","NS","C","C","NS","NS","S","NS","X","U"],
["NS","C","C","NS","C","NS","C","NS","NS","NS","NS","C","NS","NS","NS","C","NS","NS","NS","S","X","U"],
["X", "X","X","X", "X","X", "X","X", "X", "X", "X", "X","X", "X", "X", "X","X", "X", "X", "X","S","U"],
["U", "U","U","U", "U","U", "U","U", "U", "U", "U", "U","U", "U", "U", "U","U","U", "U", "U", "U","U"]]


#########FUNCTIONS#############

##isSyn: translates and annotates codon substitution
def isSyn(q20rows, p, subNT):
	wt = [entry[1] for entry in q20rows[:]]
	wtCodon = "".join(wt)
	
	wtRes = codon_table[wtCodon]
	mut = wt[:]
	mut[p] = subNT
	mutCodon = "".join(mut)
	mutRes = codon_table["".join(mutCodon)]
	if wtCodon == mutCodon:
		SNS = "WT"
	else:
		if wtRes == mutRes:
			SNS = "S"
		else:
			SNS = "NS"
	#print (wtCodon+"->"+mutCodon+":"+SNS)
	return [SNS, wtCodon, mutCodon, wtRes, mutRes]

## annotate: mainfunction for formatting q20 files
def annotate(root,Q20file):
	print ("Annotating:", Q20file)
	with open(root+"/"+Q20file,'r') as IF:
		q20 = [element for element in [line.strip().split("\t") for line in IF]]
		print ("Length of input reference: "+str(len(q20)))
	#define annotation columns and tables	
	ntTable = []
	codTable = []
	coverage = []
	res = []

	#begin nt annotation
	posVector =   [row[0] for row in q20 for sub in ("A","C","G","T")]
	wtNtVector  = [row[1] for row in q20 for sub in ("A","C","G","T")]
	mutNtVector = [sub    for row in q20 for sub in ("A","C","G","T")]
	coverage    = [ (int(row[2])+int(row[3])+int(row[4])+int(row[5])) for row in q20 for sub in ("A","C","G","T") ]	
	countVector = [ float(row[ntDict[sub]+2]) for row in q20 for sub in ("A","C","G","T") ]
	freqVector =  np.divide(countVector,coverage)

	SNS = []
	wtC = []
	muC = []
	wtR = []
	muR = []
	wtAnnot = []
	muAnnot = []
	mutClass = []
	posStack=[]
	wtStack=[]
	resStack=[]
	muStack=[]
	covStack=[]
	countStack=[]
	freqStack=[]
	orfN = []
	orfcounter = -1

	for [start,stop] in intervals:
		resPosVector = [(int(row[0])-start) for row in q20 for sub in ("A","C","G","T")]
		orfcounter+=1
		start=int(start)
		stop=int(stop)

		posStack.extend(posVector)
		wtStack.extend(wtNtVector)
		muStack.extend(mutNtVector)
		covStack.extend(coverage)
		countStack.extend(countVector)
		freqStack.extend(freqVector)

		print ("Start:"+str(start)+"  Stop:"+str(stop+1)+"  ORF Length:"+str((stop+1-start)/3))
		for pos in range(0,start-1):				#5'UTR annotation
			for m in range(4):
				resPosVector[pos*4+m] = 0
				synInfo = ["U","U","U","U","U"]		
				SNS.append(synInfo[0])
				wtC.append(synInfo[1])
				muC.append(synInfo[2])
				wtR.append(synInfo[3])
				muR.append(synInfo[4])
				wtAnnot.append([0 for P in range(18)])
				muAnnot.append([0 for P in range(18)])
				orfN.append(str(orfcounter))
				mutClass.append(pechmann[resDict[synInfo[3]]][resDict[synInfo[4]]])				

		for pos in range(start-1,stop,3):		# for codon 	### Coding Sequence
			for p in range(3): 										# for codon nt position
				for sub in ("A","C","G","T"):	
					synInfo = isSyn(q20[pos:pos+3], p, sub)					# for each sub
					SNS.append(synInfo[0])
					wtC.append(synInfo[1])
					muC.append(synInfo[2])
					wtR.append(synInfo[3])
					muR.append(synInfo[4])
					orfN.append(str(orfcounter))
					wtAnnot.append([resAnnotMatrix[P][resDict[synInfo[3]]] for P in range(18)])
					muAnnot.append([resAnnotMatrix[P][resDict[synInfo[4]]] for P in range(18)])
					mutClass.append(pechmann[resDict[synInfo[3]]][resDict[synInfo[4]]])

		for pos in range((stop),len(q20)):			#3'UTR annotation
			for m in range(4) :
				resPosVector[pos*4+m] = 0
				synInfo = ["U","U","U","U","U"]	
				SNS.append(synInfo[0])
				wtC.append(synInfo[1])
				muC.append(synInfo[2])
				wtR.append(synInfo[3])
				muR.append(synInfo[4])
				orfN.append(str(orfcounter))
				wtAnnot.append([0 for P in range(18)])
				muAnnot.append([0 for P in range(18)])
				mutClass.append(pechmann[resDict[synInfo[3]]][resDict[synInfo[4]]])	

		resStack.extend(resPosVector)
		wtVectorStr = ["\t".join([str(i) for i in row]) for row in wtAnnot]
		muVectorStr = ["\t".join([str(i) for i in row]) for row in muAnnot]
	q20annot = zip(posStack, resStack, wtStack, muStack, orfN,countStack, covStack, freqStack, SNS, mutClass, wtC, muC, wtR, muR, wtVectorStr, muVectorStr)
	return q20annot

def outputFormat(root,file,annotQ20):
	filename = root+"/"+file.split(".txt")[0]+"_annot.txt"
	with open(filename,'w') as OF:
		OF.write('ntpos\tresPos\twtNT\tmutNT\tORF\tcount\tcoverage\tfreq\tsynNonsyn\tmutClass\twtCodon\tmuCodon\twtRes\tmuRes\tacidic_wt\tacyclic_wt\taliphatic_wt\taromatic_wt\tbasic_wt\tburied_wt\tcharged_wt\tcyclic_wt\thydrophobic_wt\tlarge_wt\tmedium_wt\tnegative_wt\tneutral_wt\tpolar_wt\tpositive_wt\tsmall_wt\tsurface_wt\tstop_wt\tacidic_mu\tacyclic_mu\taliphatic_mu\taromatic_mu\tbasic_mu\tburied_mu\tcharged_mu\tcyclic_mu\thydrophobic_mu\tlarge_mu\tmedium_mu\tnegative_mu\tneutral_mu\tpolar_mu\tpositive_mu\tsmall_mu\tsurface_mu\tstop_mu\n')
		for row in annotQ20:
			count = 0
			for element in row:
				count+=1
				OF.write(str(element))
				if count < len(row):
					OF.write("\t")
			OF.write("\n")
	return filename

def combineQ20s(inputDir):
	print ("Combining Q20s...")
	fileList = []
	for root,dirs,files in os.walk(inputDir):
		for file in files: 
			if file [-9:] == 'annot.txt':
				fileList.append(file)
	uniqueFiles = set(fileList)

	for filename in uniqueFiles:
		print ("..."+filename)
		first = 0
		for root,dirs,files in os.walk(inputDir):
			for F in files:
				if filename==F:
					if first == 0: 
						first = 1
						print (filename)
						with open(root+"/"+F,'r') as IF:
							masterQ20 = [element for element in [line.strip().split("\t") for line in IF]]
					elif first==1:
						with open(root+"/"+F,'r') as IF:
							currentQ20 = [element for element in [line.strip().split("\t") for line in IF]]
						for line in range(1,len(currentQ20)):
							for i in [4,5]:
								masterQ20[line][i] = masterQ20[line][i]+currentQ20[line][i]
							masterQ20[line][6] = masterQ20[line][4]/masterQ20[line][5]
		with open(inputDir+"/master_"+filename,'w') as OF:
			for row in masterQ20:
				count = 0
				for element in row:
					count+=1
					OF.write(str(element))
					if count < len(row):
						OF.write("\t")
				OF.write("\n")
			print ("wrote "+inputDir+"/master_"+filename)


############ MAIN ###############

inputDir = sys.argv[1]

translationbreaks=list(sys.argv[2:]) #97 for DENV #10272 for DENVprint(translationbreaks)
intervals=[[int(breaks) for breaks in translationbreaks[i:i+2]] for i in range(0,len(translationbreaks),2)]

print("ORF Coordinates:")
print (intervals)
print("\n")

for root,dirs,files in os.walk(inputDir):
	for file in files:
		if "Q20.txt"==file[-7:]:				#check if it is a Q20file
			annotQ20=annotate(root,file)
			print(root+file)
			OF = outputFormat(root,file,annotQ20)

#combineQ20s(inputDir) #uncomment if combining q20counts !!!!! unstable !!!!