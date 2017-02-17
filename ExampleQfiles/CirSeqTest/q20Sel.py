#q20Sel.py
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as pp
from scipy import stats

mufile = "" 
qfile = "/GitHub/CirQuestDSV/exampleQ20_annot.txt"

mutDict={"AC":0,
"AG":1,
"AT":2,
"CA":3,
"CG":4,
"CT":5,
"GA":6,
"GC":7,
"GT":8,
"TA":9,
"TC":10,
"TG":11}

with open(mufile,'r') as mu:
	names = [line.strip().split("\t")[0] for line in mu]

with open(mufile,'r') as mu:
	muTable = [line.strip().split("\t")[1] for line in mu]
muTable = np.array(muTable[1:]).astype(float)

with open(qfile,'r') as q:
	qTable = [line.strip().split("\t") for line in q]
qTable=qTable[1:]				#strip labels

start=97
end=10272
proteinLength=3391

#print (muTable)
selecVec=np.array([])
cEntVec=np.array([])
dNdSVec=np.array([])
for resPos in range(1,proteinLength):				#Very slow implementation....
	codon   = [line for line in qTable if line[1]==str(resPos)]
	NSyn    = [line for line in codon if line[7]=="NS"]
	Syn     = [line for line in codon if line[7]=="S"]
	coverage = sum([float(line[5]) for line in codon])/len(codon)

	selec = sum( [float(line[4])/(muTable[ mutDict[line[2]+line[3]]]*coverage) for line in codon if line[2]!=line[3]])

	pNS  =    sum( [float(line[4])/coverage for line in NSyn])/(len(NSyn))
	#print(pNS)
	
	if len(Syn)>0:
		pS  = sum( [float(line[4])/coverage for line in Syn]) /(len( Syn))
	else: pS = 0.0
	#print(pS)
	print (pNS, pS)
	
	dn= -.75* (np.log(1-((4.0*pNS)/3.0)))
	ds= -.75*(np.log(1-((4.0*pS)/3.0)))
	dnds = dn/ds
	print(dnds)
	
	nstates=coverage*len(codon)

	cEnt = -sum( [(float(line[4])/nstates)*np.log2((float(line[4])/nstates)) for line in codon if float(line[4])>0.0] )
	
	selecVec=np.append(selecVec,selec)
	cEntVec=np.append(cEntVec,cEnt)
	dNdSVec=np.append(dNdSVec,dnds)


np.savetxt("huCEnt.txt", cEntVec)
np.savetxt("huSelec.txt", selecVec)
np.savetxt("hudNdS.txt", dNdSVec)
#print(selecVec)

pp.plot(selecVec,range(len(selecVec)))
pp.show()

pp.scatter(cEntVec,selecVec)
pp.show()


