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
  
  #pdf(paste(inputDir,name,"_ManhattanPlot.pdf"),width=7,height=5)
  Qf<-ggplot(filter)
  
  PP<-Qf+geom_point(aes(ntpos,freq,cex = freq,color=synNonsyn),alpha=.6)+
    ylim(1e-6,1)+geom_text(cex=2,position = position_nudge(y=.15),data=filter[filter$freq>FILTER,],aes(ntpos,freq,label=paste(wtRes,resPos,muRes)))+
    scale_y_log10()+ggtitle(name)+scale_size(trans = "sqrt",range=c(.5,2.5))+
    theme_classic()+ylab("Frequency")+xlab("Genome Position")+scale_color_brewer(palette = "Set1")+facet_grid(ORF~.)
  
  ggsave(filename = paste(inputDir,name,"_ManhattanPlot.pdf"),
         device = 'pdf',PP)
  
  write.csv(file=paste(inputDir,"/",name,"_0.005CO.csv",sep=""),filter[order(filter$freq,decreasing = T),])
  
  print(name)
}
