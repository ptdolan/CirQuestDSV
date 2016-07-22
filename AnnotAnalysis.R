inputDir<-"~/GitHub/CirQuestDSV/" #put input directory here.

FILTER=0.005 # change for labeling

library(limma)
library(ggplot2)
library(reshape2)

########## plot frequencies ##########
allq20<-NULL
for (q20file in list.files(inputDir,pattern = 'annot.txt',full.names = T)){ #Read annotated CirSeq output file
  par (mfrow=c(1,1), pch = 16, mar= c(5,5,5,5),col = rgb(0,0,0,0.5), las = 1 )
  name<-strsplit2(strsplit2(q20file,"_annot.txt")[1],"/")[length(strsplit2(strsplit2(q20file,"_annot.txt")[1],"/"))]
  print(paste("Processing:",name))
  
  passage<-as.integer(strsplit2(name,split = "-Q20")[1])
  
  q20<-read.delim(q20file,head=T)
  q20$name<-as.factor(name)
  q20$ntID<-as.factor(with(q20,paste( wtNT,ntpos,mutNT,sep="" )))
  q20$resID<-as.factor(with(q20,make.unique(paste( wtRes,resPos,muRes,sep="" ))))
  #NT entropy
  q20$HsN<-rep(unlist(by(q20$count,q20$ntpos,function(X){sum(X/sum(X)*log2(X/sum(X)))}[[1]])),each=4)
  q20$HsN[is.na(q20$HsN)]<-0.0
  
  #filt
  if(is.na(passage)){
    print("No passage")
    filter = q20[(q20$wtNT!=q20$mutNT)&q20$count>3,] #filter out WT and low counts(by binom)
    }
  else{print("not NaN")
    q20$passage<-as.factor(passage)
    filter = q20[(q20$wtNT!=q20$mutNT)&q20$count>3,] #filter out WT and low counts(by binom)
    allq20<-rbind(allq20,filter)}

  
  #pdf(paste(inputDir,name,"_ManhattanPlot.pdf"),width=7,height=5)
  
##### GGPLOT2 ######
  
  Qf<-ggplot(filter)
  PP<-Qf+geom_hline(aes(yintercept = FILTER),lty=3,color="grey")+geom_point(aes(ntpos,freq,cex = freq,color=synNonsyn),alpha=0.4)+
    geom_text(cex=2,data=filter[filter$freq>FILTER,],aes(ntpos,freq,label=ifelse(wtRes!="U",paste(wtRes,resPos,muRes),paste(wtNT,ntpos,mutNT))))+
    ggtitle(name)+scale_size(trans = "sqrt",range=c(.5,2.5))+
    theme_classic()+ylab("Frequency")+xlab("Genome Position")+
    scale_color_brewer(palette = "Set1")+facet_grid(ORF~.)
  
  ggsave(filename = paste(inputDir,name,"_log10_ManhattanPlot.pdf"),
         device = 'pdf',PP+scale_y_log10(limits=c(10E-7,1)))##Edit here ("limits=c(10E-7,1)) for different fixed depth or delete for autoscaling

  ggsave(filename = paste(inputDir,name,"_sqrt_ManhattanPlot.pdf"),
         device = 'pdf',PP+scale_y_sqrt(limits=c(0,1),breaks=c(0,.001,.01,.05,.1,.2,.4,.6,.8,1)))
  
  ##### Write output table##########
  write.csv(file=paste(inputDir,"/",name,"_",FILTER,"CO.csv",sep=""),filter[order(filter$freq,decreasing = T),])
}
Qa<-ggplot(allq20)
TT<-Qa+geom_path( aes(passage,freq, group=ntID,colour=synNonsyn),alpha=0.4)+
  ggtitle(name)+
  theme_classic()+ylab("Frequency")+xlab("Passage")+
  scale_color_brewer(palette = "Set1")

ggsave(filename = paste(inputDir,"dir_sqrt_TrajPlot.pdf"),
       device = 'pdf',TT+scale_y_sqrt(breaks=c(0,.001,.01,.05,.1,.2,.4,.6,.8,1))+facet_grid(ORF~.))

ggsave(filename = paste(inputDir,"dir_log_TrajPlot.pdf"),
       device = 'pdf',TT+scale_y_log10(breaks=c(.000001,.00001,.0001,.001,.01,.1,1))+facet_grid(ORF~.))
