# # # #
#
# AnnotAnalysis.R: For analysis and plotting of annotated q20file output from q20analysis.py . 
#
# # # # 

FILTER=0.005 # change for labeling

library(limma)
library(ggplot2)
library(reshape2)
rootDir<-"~/Research/CirSeq/Hong/Processed_Data/NewWTData_Nov9/" #put input directory here.
########## plot frequencies ##########

allq20<-NULL
fileList<-list.files(recursive = T,paste(rootDir,sep = ""),pattern = 'annot.txt',full.names = T)
if(length(fileList>0)){
  figureDir<-"DSVOutput_Figures/"
  outputDir<-"DSVOutput_Tables/"
  entropyDir<-"DSVOutput_Figures/Entropy/"
  dir.create(paste(rootDir,figureDir,sep=""))
  dir.create(paste(rootDir,outputDir,sep=""))
  dir.create(paste(rootDir,entropyDir,sep=""))
  
  for (q20file in fileList){ #Read annotated CirSeq output file
    print(q20file)
    name<-strsplit2(strsplit2(q20file,"_annot.txt")[1],"/")
    
    #Read File
    name
    q20<-read.delim(q20file,head=T)
    q20$genotype<-as.factor("WT")
    q20$rep<-as.factor(name[10])
    passage<-strsplit2(name[11],"p")[2]
    q20$passage<-as.numeric(passage)
    #Name Interpretation
    headerInfo<-paste(name[10:11],collapse = "-")
    headerVec<-strsplit2(headerInfo,"-")
    headLen<-length(headerVec)
    
    if(is.na(passage)){ #No passage info.
      print(paste("No passage info:", name))
      passage<-0
      print(passage)
      header<-paste(headerVec,collapse = "")
    }else{if(headLen>1){
      HeadLabel<-paste(headerVec[1:(headLen-1)],collapse = "")
      print(paste("Header:", HeadLabel, "Passage:", passage))
      header<-paste(headerVec[1:(headLen-1)],collapse = "")
    }else{
      header="Sample"
      print(paste("No Header Info:",headerInfo,"Passage:", passage))
    } 
      
    } 
    
    #labels
    q20$name<-as.factor(headerInfo)
    q20$ntID<-as.factor(with(q20,paste( wtNT,ntpos,mutNT,sep="" )))
    q20$resID<-as.factor(with(q20,make.unique(paste( wtRes,resPos,muRes,sep="" ))))
    q20$header<-as.factor(header)
    q20$HsN<-rep(unlist(by(q20$count,q20$ntpos,function(X){-sum(na.rm = T, (X/sum(X)) * log((X/sum(X)),base = 4))}[[1]]) ),each=4)
    q20$HsN[is.na(q20$HsN)]<-0.0
    
    filter<- q20[(q20$wtNT!=q20$mutNT)&q20$count>3,] #filter out WT and low counts(by binom)
    allq20<- rbind(allq20,q20)
    save(list = "allq20",file =  paste(rootDir,"allq20.Rdata",sep = "/"))
#     Q<-ggplot(q20[q20$wtNT!=q20$mutNT,])
    
    ### Entropy
#     Qh<-Q+geom_line(aes(ntpos,HsN),position = 'dodge',stat="identity")+
#       geom_smooth(span=10,aes(ntpos,HsN),method = )+
#       ylab("Entropy")+xlab("Genome Position")
#     
#     ggsave(filename = paste(rootDir,figureDir,"Entropy/",headerInfo,"_Entropy.pdf",sep=""),
#            device = 'pdf',Qh)
#     
#     PU<-Q+geom_hline(aes(yintercept = FILTER),lty=3,color="grey")+
#       geom_point(aes(ntpos,freq,cex = freq,color=synNonsyn,alpha=count<3))+
#       scale_alpha_discrete(range=c(.5,.1))+
#       geom_text(cex=2,
#                 data=q20[q20$wtNT!=q20$mutNT&q20$freq>FILTER,],
#                 aes(ntpos,freq,label=ifelse(wtRes!="U",as.character(resID),as.character(ntID))))+
#       theme_classic()+ylab("Frequency")+xlab("Genome Position")+
#       scale_color_brewer(palette = "Set1")+ggtitle(headerInfo)+facet_grid(ORF~.)
#     PU+scale_y_log10(limits=c(10E-7,1))
#     
#     ggsave(filename = paste(rootDir,figureDir,headerInfo,"_log10_ManhattanPlot.pdf",sep=""),
#            device = 'pdf',PU+scale_y_log10(limits=c(10E-7,1))) #Edit here ("limits=c(10E-7,1)) for different fixed depth or delete for autoscaling
#     
#     ggsave(filename = paste(rootDir,figureDir,headerInfo,"_sqrt_ManhattanPlot.pdf",sep=""),
#            device = 'pdf',PU+scale_y_sqrt(limits=c(0,1),breaks=c(0,.001,.01,.05,.1,.2,.4,.6,.8,1)))
#     
    ##### Write output table##########
    write.csv(file=paste(rootDir,outputDir,headerInfo,"_",FILTER,"CO.csv",sep=""),filter[order(filter$freq,decreasing = T),])
  }
  
  #WholeDirectoryOutput
  output<-dcast(data = allq20,formula = wtNT+ntpos+mutNT+wtRes+resPos+muRes~header+passage,value.var = "freq")
  outputCounts<-dcast(data = allq20,formula = wtNT+ntpos+mutNT+wtRes+resPos+muRes~header+passage,value.var = "count")
  outputCover<-dcast(data = allq20,formula = wtNT+ntpos+mutNT+wtRes+resPos+muRes~header+passage,value.var = "coverage")
  write.csv(file=paste(rootDir,outputDir,"DirFreqTable.csv",sep=""),output)
  #maskedHeat<- d3heatmap::d3heatmap( output[ rowSums(output[,7:ncol(output)])>2, 7:ncol(output) ], labRow=
  #                                 matrix(apply(output[ rowSums(output[,7:ncol(output)])>2, 1:6 ],MARGIN = 1,paste,collapse=".")))
  
  masked<-output[ which(rowSums(output[,7:ncol(output)])<2),]
  
  
  
  #Princomp
  
  
  PCallq20<-prcomp(masked[,7:ncol(masked)])

  Dist<-dist(t(masked[,7:ncol(masked)]),method = "man")
  distDF<-as.data.frame(cmdscale(Dist))
  distDF<-data.frame(distDF,genotype=strsplit2(rownames(distDF),split = "_")[,3])
  distDF$passage<-strsplit2(rownames(distDF),split="_")[,4]
  PCDF<-data.frame(PCallq20$rotation,genotype=strsplit2(rownames(distDF),split = "Rep")[,1],passage=strsplit2(rownames(distDF),split="_")[,2])
  ggplot(PCDF)+geom_point(aes(PC1,PC2,color=genotype))#+geom_text(aes(PC1,PC2,label=passage),nudge_x = -0.01)
  ggplot(PCDF)+geom_point(aes(PC1,PC3,color=genotype))+geom_text(aes(PC1,PC3,label=))
  mdsPlot<-ggplot(distDF)+geom_point(aes(V1,V2,color=genotype),alpha=0.5)+scale_color_brewer(palette="Paired")+geom_text(aes(V1,V2,label=passage),nudge_x = -0.1)
  
  ggsave(width=7, height=5,paste(rootDir,figureDir,"geneticDistPlot.pdf",sep=""),mdsPlot)
  
  maskDF<-data.frame(masked)
  
  ggplot(maskDF)+geom_point(aes(ntpos,WTRep1_7),col="red")+geom_point(aes(ntpos,G64SD79HRep1_7))+scale_y_sqrt()+theme_classic()
  ggplot(maskDF)+geom_point(aes(G64SD79HRep1_7,WTRep1_7),alpha=.8)+scale_x_sqrt(limits=c(0,.25))+scale_y_sqrt(limits=c(0,.25))
  
  #Trajectories
  Qa<-ggplot(allq20[allq20$passage!=0,])
  
  ggsave(filename = paste(rootDir,figureDir,headerInfo,"Coverage.tiff",sep=""),
         Qa+
           geom_point(aes(ntpos,coverage,color = factor(passage)),cex = 0.1)+
           geom_smooth(se = T,aes(ntpos,coverage,color=factor(passage)),color="black")+
           theme_classic()+ylab("Reads per position")+xlab("Genome Position")+
           scale_color_brewer(palette = "Set1")+scale_color_brewer("seq",palette = "RdBu")+
           scale_y_log10(limits=c(1,3e6))+facet_wrap(~header))
  
  TT<-Qa+geom_path( aes(passage,freq, group=ntID,colour=synNonsyn),alpha=0.4)+
    theme_classic()+ylab("Frequency")+xlab("Passage")+scale_x_discrete()+
    scale_color_brewer(palette = "Set1")
  
  ggsave(filename = paste(rootDir,figureDir,"dir_sqrt_TrajPlot.pdf",sep=""),
         device = 'pdf',TT+scale_y_sqrt(breaks=c(0.00001,0.0001,.001,.01,.05,.1,.2,.4,.6,.8,1))+facet_wrap(ORF~header))
  
  ggsave(filename = paste(rootDir,figureDir,"dir_log_TrajPlot.pdf",sep=""),
         device = 'pdf',TT+scale_y_log10(breaks=c(10**c(-6:0)))+facet_wrap(ORF~header))
  
  #Qa+geom_bar(aes(ntpos,HsN),position = "dodge",stat = "identity")+facet_grid(passage~header)
  ggsave(filename = paste(rootDir,figureDir,"allFreqPlots.tiff",sep=""),dpi = 300,width = 7.2,height = 9.8,
         Qa+geom_point(aes(ntpos,freq,color=synNonsyn,alpha=count>3),cex=0.3)+
           scale_alpha_discrete(range=c(0.1,0.5))+
           scale_y_log10()+
           facet_grid(passage~header)
  )
  
  ggsave(filename = paste(rootDir,figureDir,"allFreqPlots.eps",sep=""),width = 7.2,height = 9.8,
         Qa+geom_point(aes(ntpos,freq,color=synNonsyn,alpha=count>3),cex=0.3)+
           scale_alpha_discrete(range=c(0.1,0.5))+
           scale_y_log10()+
           facet_grid(passage~header)
  )
  #Qa+geom_bar(aes(ntpos,HsN),position = "dodge",stat = "identity")+facet_grid(passage~header)
  
}else{print("Error: No Files Found.")}

