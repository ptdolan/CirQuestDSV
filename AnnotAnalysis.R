

inputDir<-"~/GitHub/CirQuestDSV/" #put input directory here.

FILTER=0.005 # change for labeling

library(limma)
library(ggplot2)
library(reshape2)

placeNames<-function(Limits){
  namePos<-NULL
  N=1
  while(N < length(Limits)){
    namePos<-append(namePos,mean(Limits[N],Limits[N+1]))
    N=N+1}
  return(namePos)
}
nameList<-NULL

########## plot frequencies ##########


compq20<-NULL

for (q20file in list.files(inputDir,pattern = 'annot.txt',full.names = T)){ #Read annotated CirSeq output file
  
  par (mfrow=c(1,1), pch = 16, mar= c(5,5,5,5),col = rgb(0,0,0,0.5), las = 1 )
  name<-strsplit2(strsplit2(q20file,"_annot.txt")[1],"//")[2]
  q20<-read.delim(q20file,head=T)
  q20$name<-as.factor(name)
  q20$HsN<-rep(unlist(by(q20$count,q20$ntpos,function(X){sum(X/sum(X)*log2(X/sum(X)))}[[1]])),each=4)
  filter = q20[q20$wtNT!=q20$mutNT,] #filter out WT
  compq20<-rbind(compq20,filter)
  
  assign(name,read.delim(q20file,sep="\t"))  #?
  q20<-get(name)
  pdf(paste(inputDir,name,"_ManhattanPlots.pdf"),width=7,height=5)
  plot(log10(filter$freq) ~ filter$ntpos, main = name,ylim = c(-5.5,0.5),ylab = "Allele Frequency",xlab = "Genome Position",
       pch = 16,
       cex = ifelse(filter$count>4,(4.3+log10(filter$freq))/5+0.2,0.05),
       col =  ifelse(filter$synNonsyn=="U", rgb(0,0,0,0.3), 
              ifelse(filter$synNonsyn=="S", rgb(0,0,.8,0.5),
              ifelse(filter$synNonsyn=="NS", rgb(0.8,0,0,0.6),'green'))))
  abline(h=log10(FILTER))
  biggies<-filter[filter$freq>FILTER,] #greater than 5e-3 (1 in 200)
  text(cex=0.4,log10(biggies$freq)~biggies$ntpos, labels=paste(biggies$wtRes,biggies$resPos,biggies$muRes,sep=""),pos=2)
  
  legend(x=0,y=0,legend = c("S","NS","UTR"),col = c(rgb(0,0,0.8,0.5),rgb(0.8,0,0.6),rgb(0,0,0,0.3)), pch = 16,cex = 0.7)
  
  write.csv(file=paste(inputDir,"/",name,"_0.005CO.csv",sep=""),filter[order(filter$freq,decreasing = T),])
  
  print(name)
  dev.off()
}
