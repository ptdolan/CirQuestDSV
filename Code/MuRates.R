library(ggplot2)
library(data.table)
library(reshape2)

#MutationRatesAnalysis
procHdir<-function(hdir="",patt){
  hdir<-ifelse(hdir=="",getwd(),hdir)
  list<-list.files(hdir,pattern = patt)
  all<-data.frame()
  for (i in list){
    curr<-fread(paste(hdir,i,sep = ""),header = T,sep = "\t")
    #print(head(curr))
    curr$name<-(limma::strsplit2(i,"_")[2])
    curr$passage<-(limma::strsplit2(i,"_")[3])
    all<-rbind(all,curr)
  }
  return(all)
}

getRef<-function(hdir="",patt){
  hdir<-ifelse(hdir=="",getwd(),hdir)
  list<-list.files(hdir,pattern = patt)
  refInput<-read.delim(paste(hdir,list[1],sep = ""),header = T,sep = "\t")
  #print(refInput)
  REF<-as.character(refInput$wtNT)
  #print(REF)
  consensus<-by(refInput$wtNT,refInput$ntpos,function(X){as.character(X[1])})
  region<-lapply(1:length(consensus),FUN = function(X){
        paste(consensus[max(X-11,1):min(length(consensus),(X+11))],collapse="")
      })
  return(unlist(region))
}

region<-getRef(hdir,"phred")

ggplot2::ggplot(all)+
  geom_pointrange(fatten = 0.8,aes(name,y = MLE,ymin=MLE-SE,ymax=MLE+SE),alpha=0.5,cex=0.5)+
  scale_y_log10()+facet_wrap(~Type,)

hdir="~/Research/Zika/CirSeq_InDelRun/71017/MLE_Mu_Estimates/"
allMus<-procHdir(hdir,patt="MLE")

hdir="~/Research/Zika/CirSeq_InDelRun/71017/MM_071017_phred20/"
all20s<-procHdir(hdir,patt="annot.txt")

all20s<-all20s[all20s$wtNT!=all20s$mutNT,]
all20s$mut=paste(all20s$wtNT,all20s$mutNT,sep="")

all20s$region<-rep(region,each=3)
features<-c(1,107+3*c(6,128,217,292,797,1158,1370,1521,2125,2270,2774),10378)
featurenames<-c("utr5","C","pr","M","E","NS1","NS2A","NS2B","NS3","NS4A","NS4B","NS5","utr3")
featuretable<-data.frame(features,featurenames)
str(featuretable)
all20s$feature<-"utr3"
for(row in 2:nrow(featuretable)){
  feat<-as.character(featuretable$featurenames[row-1])
  all20s$feature[all20s$ntpos>=featuretable$features[row-1]&all20s$ntpos<featuretable$features[row]]=feat
}
all20s$feature<-as.factor(all20s$feature)


#####

filt20s<-all20s

ggplot(filt20s)+geom_vline(data=allMus,aes(xintercept = MLE,alpha = 0.3,color=Type))+geom_density(data=filt20s,aes(freq,fill=mut))+
  scale_x_log10()+facet_grid(wtNT~mutNT)

newall<-data.table()
by(filt20s,filt20s$mut,function(f){
  new<-(f[order(f$freq),])
  print(tail(new))
  l<-max(new$ntpos)
  new<-new[new$ntpos>12&(new$ntpos<l-12),]
  print(tail(new))
  l<-nrow(new)
  bottom<-new$region[1:(l*.20)]
  top<-new$region[(l*.80):l]
  all<-new$region
  newall<-rbind(newall,new)
  write.table(file = paste("~/Research/CirSeq/Hong/ZikaEditing_8-17/Seqs_top_",f$mut[1],".txt",sep=""),top,quote = F,row.names = F,col.names = F)
  write.table(file = paste("~/Research/CirSeq/Hong/ZikaEditing_8-17/Seqs_bottom_",f$mut[1],".txt",sep=""),bottom,quote = F,row.names = F,col.names = F)
  write.table(file = paste("~/Research/CirSeq/Hong/ZikaEditing_8-17/Seqs_all_",f$mut[1],".txt",sep=""),all,quote = F,row.names = F,col.names = F)
})



ggplot(filt20s[mut=="ag"&filt20s$passage=="p1",])+geom_point(aes(ntpos,freq,color=feature),alpha=0.4,stat='identity')+
  geom_vline(data=featuretable, xintercept = features)+scale_y_log10()

out<-dcast(filt20s,formula = ntpos+passage+name~mut,na.rm = T,value.var = "freq",fill=0)
melted<-(melt(out,variable.name = 'mut',id.vars = c("ntpos","passage","name")))
rolled<-data.table(apply(out[,-c(1:3)],MARGIN = 2,function(X){zoo::rollmedian(X,k=200,fill=0)}))
merged<-data.table(out,rolled=rolled)
melted<-melt(merged,variable.name = 'mut',id.vars = c("ntpos","passage","name"))

ggplot(melted)+
  geom_line(aes(ntpos,value,color=mut))+
  scale_y_log10()+
  facet_grid(passage~name)

muttable<-data.table(melted)
