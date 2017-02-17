#CMHsel

FILTER=0.005 # change for labeling

library(limma)
library(ggplot2)
library(reshape2)
#rootDir<-"~/GitHub/CirQuestDSV/DengueData/" #put input directory here. DENGUE DATA

rootDir<-"~/Research/CirSeq/Hong/Processed_Data/Renamed/" #put input directory here.

load(paste(rootDir,"allq20.Rdata",sep="/"))

#Functions
cmhO<-function(amu,awt,bmu,bwt){
  #print(amu)
  n = sum(amu,awt,bmu,bwt)
  OR<-log2(((amu*bwt)/n)/((awt*bmu)/n))
  return(sum(OR))
}

windowCMHO<-function(wts,mus,a,b,window=c(0,9,25,50)){
  print(paste("Windows:", window))
  lim=max(mus$Pos)
  ORmatrix<-1:lim
  for(win in window){
    print(paste("Current Window:", win))
    output<-lapply(1:lim,
                   FUN = function(row){
                     OR=0
                     #print(row/lim)
                     if(any( mus$Pos %in% c(row:(row-win)) )){
                       OR<-cmhO(
                         mus[mus$Pos%in%c(row:(row-win)),a],
                         wts[mus$Pos%in%c(row:(row-win)),a],
                         mus[mus$Pos%in%c(row:(row-win)),b],
                         wts[mus$Pos%in%c(row:(row-win)),b])
                     }
                     return(c(row,OR))})
    CMHvec<-matrix(unlist(output),ncol=2,byrow=T)[,2]
    ORmatrix<-cbind(ORmatrix,CMHvec)
  }
  print("returning... Matrix")
  DF<-as.data.frame(ORmatrix)
  colnames(DF)<- c("Pos",window)
  melted<-melt(DF,id.vars = 1,variable.name = "Window",measure.vars = 2:5)
  return(melted)}

generateCounts<-function(q20cast){
  NTcounter=0
  print("There's...")
  q20cast$wtCount<- 0
  q20cast$muCount<- 0
  nucs<-as.character(sort(unique(as.character(q20cast$wtNT))))
  for(NT in nucs){
    print(paste(NT,"'s",sep=""))
    NTcounter<-NTcounter+1
    wtVec<-q20cast[as.character(q20cast$wtNT)==NT,5+NTcounter]
    sumVec<-rowSums(q20cast[as.character(q20cast$wtNT)==NT,6:9])
    
    q20cast$wtCount[q20cast$wtNT==NT]<-wtVec
    #print(wtVec)
    
    q20cast$muCount[q20cast$wtNT==NT]<-(sumVec-wtVec)
  }
  return(q20cast)
}


p<-dcast(allq20,rep+passage+genotype+ntpos+wtNT~mutNT,value.var = "count",fill = 0)

pcounts<-generateCounts(p)
head(pcounts)
#pcounts$wtCount

wts<-dcast(pcounts,ntpos+wtNT~genotype+passage,value.var = "wtCount",fun.aggregate = sum)
mus<-dcast(pcounts,ntpos+wtNT~genotype+passage,value.var = "muCount",fun.aggregate = sum)

FET_Pg<-unlist(lapply(1:length(wts$WT_7),function(X){ fisher.test(alternative = 'g',matrix(c(wts$G64S_7[X],mus$G64S_7[X],wts$WT_7[X],mus$WT_7[X]),byrow = T,nrow = 2))$p.value }))
FET_Pl<-unlist(lapply(1:length(wts$WT_7),function(X){ fisher.test(alternative = 'l',matrix(c(wts$G64S_7[X],mus$G64S_7[X],wts$WT_7[X],mus$WT_7[X]),byrow = T,nrow = 2))$p.value }))

FET_OR<-unlist(lapply(1:length(wts$WT_7),function(X){ fisher.test(alternative = 'g',matrix(c(wts$G64S_7[X],mus$G64S_7[X],wts$WT_7[X],mus$WT_7[X]),byrow = T,nrow = 2))$estimate}))
plot(log10(FET_OR),type ='h')
