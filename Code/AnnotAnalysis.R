# # # #
#
# AnnotAnalysis.R: For analysis and plotting of annotated q20file output from q20analysis.py . 
# # Usage: commandline Rscript  /path/to/AnnotAnalysis.R /path/to/annotQ20DIR/ ##NOTE: not file. 
#
# # # # 
print("Generating R plots and text outputs...")

FILTER=0.005 # change for labeling

require(limma)
require(ggplot2)
require(reshape2)
args = commandArgs(trailingOnly=TRUE)
rootDir<-args
########## plot frequencies ##########
rootDir="~/Research/CirSeq/Hong/newcopiedset/"#Can change input dir here if not using command line arg
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
    L=length(name)
    headString<-strsplit2(name[L],"-")[1] ## "..._1"
    headString
    splitHead=strsplit2(headString,"_")
    if(length(splitHead)>1){
      headerInfo=paste(splitHead[1:(length(splitHead)-1)],collapse = " ")
      passage=splitHead[length(splitHead)]
    }
    else{passage<-headerInfo<-splitHead}
    #Read File
    q20<-read.delim(q20file,head=T)
    #print(name)
    q20$passage<-as.numeric(passage)
    q20$ORF<-as.factor(q20$ORF)
    #Name Interpretation
    headLen=length(splitHead)
    if(is.na(passage)){ #No passage info.
      print(paste("No passage info:", name))
      passage<-0
      header<-headerInfo
    }else{if(headLen>1){
      print(paste("Header:", headerInfo, "Passage:", passage))
      header<-headerInfo
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
    q20$Conf<-q20$count<4
    allq20<- rbind(allq20,q20)
    
    Q<-ggplot(q20[q20$wtNT!=q20$mutNT,])

    ### Entropy
    Qh<-Q+geom_line(aes(ntpos,HsN),position = 'dodge',stat="identity")+
      ylab("Entropy")+xlab("Genome Position")

    ggsave(width=6,height=4, filename = paste(rootDir,figureDir,"Entropy/",headString,"_Entropy.pdf",sep=""),
           device = 'pdf',Qh)

    PU<-Q+geom_hline(aes(yintercept = FILTER),lty=3,color="grey")+
      geom_point(aes(ntpos,freq,cex = freq,color=synNonsyn,alpha=Conf))+
      scale_alpha_discrete(range=c(.5,.1))+
      scale_size(range = c(0.2,3))+
      geom_text(cex=2,
                data=q20[q20$wtNT!=q20$mutNT&q20$freq>FILTER,],
                aes(ntpos,freq,label=ifelse(wtRes!="U",as.character(resID),as.character(ntID))))+
      theme_classic()+ylab("Frequency")+xlab("Genome Position")+
      scale_color_brewer(palette = "Set1")+ggtitle(headerInfo)+facet_grid(ORF~.)
    
    PU+scale_y_log10(limits=c(10E-7,1))

    ggsave(width=6,height=4, filename = paste(rootDir,figureDir,headString,"_log10_ManhattanPlot.pdf",sep=""),
           device = 'pdf',PU+scale_y_log10(limits=c(10E-7,1))) #Edit here ("limits=c(10E-7,1)) for different fixed depth or delete for autoscaling

    ggsave(width=6,height=4, filename = paste(rootDir,figureDir,headString,"_sqrt_ManhattanPlot.pdf",sep=""),
           device = 'pdf',PU+scale_y_sqrt(limits=c(0,1),breaks=c(0,.001,.01,.05,.1,.2,.4,.6,.8,1)))

    ##### Write output table##########
    write.csv(file=paste(rootDir,outputDir,headString,"_",FILTER,"CO.csv",sep=""),filter[order(filter$freq,decreasing = T),])
  }
  save(list = "allq20",file =  paste(rootDir,"allq20.Rdata",sep = "/"))
  
  #WholeDirectoryOutputs
  output<-dcast(data = allq20,formula = wtNT+ntpos+mutNT+wtRes+resPos+muRes+Conf~header+passage,value.var = "freq")
  outputCounts<-dcast(data = allq20,formula = wtNT+ntpos+mutNT+wtRes+resPos+muRes+Conf~header+passage,value.var = "count")
  outputCover<-dcast(data = allq20,formula = wtNT+ntpos+mutNT+wtRes+resPos+muRes+Conf~header+passage,value.var = "coverage")
  write.csv(file=paste(rootDir,outputDir,"DirFreqTable.csv",sep=""),output)
  
  if(ncol(output)>8){
    masked<-data.frame(na.omit(output[output$wtNT!=output$mutNT&output$wtRes!="U",]))
    #Princomp
    
#masked=#Some all integer DF 
#MDS
    #Euclidian ("Crow's Flight") distance
    DistE<-dist(t(masked[,8:ncol(masked)]),method = "euclidian")
    distDF<-as.data.frame(cmdscale(DistE))
    distDF$header<-strsplit2(rownames(distDF),split="_")[,1]
    distDF$passage<-strsplit2(rownames(distDF),split="_")[,2]
    mdsPlot<-ggplot(distDF)+geom_point(aes(V1,V2,color=header),alpha=0.5)+ggtitle("Euclidian Genetic Distance MDS")+xlab("Genetic Distance - Dim. 1")+ylab("Genetic Distance - Dim. 2")+scale_color_brewer(palette="Paired")+geom_text(aes(V1,V2,label=passage))
    
    #Manhattan ("Taxi Cab") distance
    DistM<-dist(t(masked[,8:ncol(masked)]),method = "manhattan")
    distDF<-as.data.frame(cmdscale(DistM))
    distDF$header<-strsplit2(rownames(distDF),split="_")[,1]
    distDF$passage<-strsplit2(rownames(distDF),split="_")[,2]
    mdsPlot<-ggplot(distDF)+geom_point(aes(V1,V2,color=header),alpha=0.5)+ggtitle("Manhattan Genetic Distance MDS")+xlab("Genetic Distance - Dim. 1")+ylab("Genetic Distance - Dim. 2")+scale_color_brewer(palette="Paired")+geom_text(aes(V1,V2,label=passage))
    
#PCA    
    PCallq20<-prcomp(masked[,8:ncol(masked)])
    PCDF<-data.frame(PCallq20$rotation,header=strsplit2(rownames(distDF),split = "_")[,1],passage=strsplit2(rownames(distDF),split="_")[,2])
    ggsave(width=6,height=4, filename = paste(rootDir,figureDir,"AllPop_PCA.pdf",sep=""),
           ggplot(PCDF)+geom_point(aes(PC1,PC2,color=header))+
           geom_text(aes(PC1,PC2,label=passage),cex=2.5)+
           scale_color_brewer(palette="Paired"))
    
    ggsave(width=6, height=4,paste(rootDir,figureDir,"geneticDistPlot_MDS.pdf",sep=""),mdsPlot)
    
    #Trajectories
    Qa<-ggplot(allq20[allq20$passage!=0,])
    
    ggsave(width=6,height=4, filename = paste(rootDir,figureDir,headerInfo,"Coverage.tiff",sep=""),
           Qa+
             geom_point(aes(ntpos,coverage,color = factor(passage)),cex = 0.1)+
             geom_smooth(se = T,aes(ntpos,coverage,color=factor(passage)),color="black")+
             theme_classic()+ylab("Reads per position")+xlab("Genome Position")+
             scale_color_brewer(palette = "Set1")+scale_color_brewer("seq",palette = "Spectral")+
             scale_y_log10(limits=c(1,3e6))+facet_wrap(~header))
    
    TT<-Qa+geom_path( aes(passage,freq, group=ntID,colour=synNonsyn),alpha=0.4)+
      theme_classic()+ylab("Frequency")+xlab("Passage")+scale_x_discrete()+
      scale_color_brewer(palette = "Set1")
    
    ggsave(width=6,height=4, filename = paste(rootDir,figureDir,"dir_sqrt_TrajPlot.pdf",sep=""),
           device = 'pdf',TT+scale_y_sqrt(breaks=c(0.00001,0.0001,.001,.01,.05,.1,.2,.4,.6,.8,1))+facet_wrap(ORF~header))
    
    ggsave(width=6,height=4, filename = paste(rootDir,figureDir,"dir_log_TrajPlot.pdf",sep=""),
           device = 'pdf',TT+scale_y_log10(breaks=c(10**c(-6:0)))+facet_wrap(ORF~header))
    
    #Qa+geom_bar(aes(ntpos,HsN),position = "dodge",stat = "identity")+facet_grid(passage~header)
    ggsave(filename = paste(rootDir,figureDir,"allFreqPlots.tiff",sep=""),dpi = 300,width = 7.2,height = 9.8,
           Qa+geom_point(aes(ntpos,freq,color=synNonsyn,alpha=count>3),cex=0.3)+
             scale_alpha_discrete(range=c(0.1,0.5))+
             scale_y_log10()+
             facet_grid(passage~header)
    )
    
    # ggsave(filename = paste(rootDir,figureDir,"allFreqPlots.",sep=""),width = 7.2,height = 9.8,
    #        Qa+geom_point(aes(ntpos,freq,color=synNonsyn,alpha=count>3),cex=0.3)+
    #          scale_alpha_discrete(range=c(0.1,0.5))+
    #          scale_y_log10()+
    #          facet_grid(passage~header)
    # )
    #Qa+geom_bar(aes(ntpos,HsN),position = "dodge",stat = "identity")+facet_grid(passage~header)
    print("Completed Successfully!")
    }
}else{print("Error: No Files Found. Add as Arguments?")}
