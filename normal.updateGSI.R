Coloncancerdatapre<-function(){
  rlt<-list()
  setwd("/home/shg047/oasis/monod/hapinfo/June")
  saminfo<-read.table("/home/shg047/oasis/monod/saminfo.txt",sep="\t")
  #data<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/June/monod.mhl.june25.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
  #save(data,file="monod.mhl.june25.RData")
  load("monod.mhl.june25.RData")
  data<-data[,-grep("CTT-|PC-P|PC-T",colnames(data))]
  colnames(data)[grep("STL",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("STL",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("WB",colnames(data))]<-"WBC"
  colnames(data)[grep("N37",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("N37",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("methylC",colnames(data))]<-"H1"
  colnames(data)[grep("6-T|Colon_Tumor_Primary",colnames(data))]<-"CCT"
  colnames(data)[grep("7-T-",colnames(data))]<-"LCT"
  colnames(data)[grep("6-P-",colnames(data))]<-"CCP"
  colnames(data)[grep("7-P-",colnames(data))]<-"LCP"
  colnames(data)[grep("NC-P-",colnames(data))]<-"NCP"
  data1<-data[,grep("Brain|Stomach|Lung|Heart|Colon|WBC|Intestine|Liver|Esophagus|Kidney",colnames(data))]
  data1_ref<-data1[,-(grep("LCP",colnames(data1)))]
  base=unlist(lapply(strsplit(colnames(data1_ref),"[.]"),function(x) x[[1]]))
  colnames(data1_ref)<-base
  gsi<-data.frame(GSIMethDecon(data1_ref))
  group<-as.character(unique(gsi$group))
  signatures<-AverageWithGroup(data1_ref)  # QR
  rlt$GSI<-gsi
  rlt$data<-data1
  rlt$ref<-data1_ref
  rlt$signature<-signatures
  return(rlt)
}

GSIMethDecon<-function(data){
  data<-data.matrix(data)
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-gmaxgroup<-avebase<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(data[i,][which(index==group[j])],na.rm=T))/10^(mean(data[i,][which(index==gmax)],,na.rm=T)))/(length(group)-1)
      gsit<-gsit+tmp
    }
    ave<-tapply(data[i,], index, function(x) mean(x,na.rm=T))
    gmaxgroup<-c(gmaxgroup,gmax)
    gsi<-c(gsi,gsit)
    avebase<-rbind(avebase,ave)
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi,avebase)
  return(rlt)
}


AverageWithGroup<-function(data){
  base=unlist(lapply(strsplit(colnames(data),"[.]"),function(x) x[[1]]))
  matrix=apply(data,1,function(x) tapply(x,base,function(x) mean(x,na.rm=T)))
  matrix<-t(matrix)
  rownames(matrix)=rownames(data)
  matrix<-matrix[!rowSums(!is.finite(matrix)),]
  return(matrix)
}

RandomSamplingMean<-function(data,number=round(ncol(data)/2)){
  rlt<-c()
  for(i in 1:ncol(data)){
    VirtualSample<-rowMeans(data[,sample(1:ncol(data),number)],na.rm=T)  
    rlt<-cbind(rlt,VirtualSample)
  }
  rlt
}

TopGSIByCategory<-function(gsi,top=2,thresHigh=0.3,thresLow=0.1,plotfigure=F,figprefix="tissuespecific"){
  GSIRlt<-c()
  group<-as.character(unique(gsi$group))
  rank<-c(rep(top,length(group)))
  otf1<-paste(figprefix,"boxplot.pdf",sep="")
  otf2<-paste(figprefix,"heatmap.pdf",sep="")
  parnum<-ceiling(sqrt(length(group)))
  pdf(otf1)
  par(mfrow=c(parnum,parnum),oma = c(2,2,2,2) + 0.1,mar = c(2,2,2,2) + 0.1,cex.axis=0.75, cex.lab=0.75)
  for (i in 1:length(group)){
    subset=gsi[which(gsi$group==group[i]),]
    rexclhigh<-which(apply(subset,1,function(x) x[grep(group[i],colnames(gsi))]<0.2))
    xx<-subset[,-grep(group[i],colnames(gsi))]
    rexcllow<-which(apply(xx,1,function(x) any(as.numeric(x[4:length(x)])>0.1)))
    rexcl<-c(rexclhigh,rexcllow)
    subset=subset[-rexcl,]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
    if(plotfigure==TRUE){
      zz=subset[which(subset$group==group[i]),]
      boxplot(zz[,4:ncol(zz)],horizontal=T,las=2,col="red")
    }
  }
  dev.off()
  
  if(plotfigure==TRUE){
    HeatMap(data=data.matrix(subset[,4:ncol(subset)]),phen=gsub("AVE.","",colnames(subset)[4:ncol(subset)]),figure=otf2)
  }
  return(GSIRlt)
}

HeatMap<-function(data,phen,figure="heatmap.pdf",cexRow = 0.01,cexCol = 1.2,Colv=T,Rowv=T){
  library("gplots")
  colnames(data)=phen
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  colors <-greenred(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  ColSideColors=sidecol(phen)
  pdf(figure)
  heatmap.2(data,trace="none",cexRow = cexRow,cexCol = cexCol, ColSideColors=ColSideColors,density.info="none",col=colors,Colv=Colv,Rowv=Rowv,keysize=0.9, margins = c(5, 10))
  dev.off()
}

FigurePrepareSimulation<-function(){
  library("ggplot2")
  Fig<-data.frame(seq(0,1,by=0.05),Rlt$out.all[,grep("CCT",colnames(Rlt$out.all))],Rlt$out.all[,grep("WB",colnames(Rlt$out.all))])
  colnames(Fig)<-c("Simulated","CCT","WB")
  Fig<-100*Fig
  c <- ggplot(Fig, aes(Simulated,CCT))
  c <- c + xlab("Simulated contribution (%)")+ylab("Predicted contribution (%)")
  c <- c + stat_smooth(se = TRUE,n=10,size=1.5) + geom_point(size=3)
  c <- c + stat_smooth(se = TRUE,n=10,size=1.5) + geom_point(size=3)
  c <- c + theme_bw()
  c <- c + theme(axis.text=element_text(size=10), axis.title=element_text(size=14,face="bold"))
  ggsave(c,file="coloncancer-deconvolution.simultaion-2.pdf")
  save.image(file="Colon-Deconvolution.RData")
}


ci95<-function(x){
  x<-data.frame(x)
  out<-c()
  for(i in 1:ncol(x)){
    error <- qt(0.975,df=length(x[,i])-1)*sd(x[,i])/sqrt(length(x[,i]))
    m<-round(mean(x[,i]),4)
    d<-round(mean(x[,i])-error,4)
    u<-round(mean(x[,i])+error,4)
    rlt<-paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
    out<-cbind(out,rlt)
  }
  return(out)
}

outputprefix="Normal.plasma.0.05.top100"
library("DeconRNASeq")
library("ggplot2")
Colon<-Coloncancerdatapre()
save(Colon,file=paste(outputprefix,".input.RData",sep=""))
data1=Colon$data
data1_ref=Colon$ref
signatures<-AverageWithGroup(data1_ref)  
topgsi<-TopGSIByCategory(Colon$GSI,top=100,thresHigh=0.3,thresLow=0.1,plotfigure=TRUE,figprefix=outputprefix)
VirtualMatrix<-data.frame(RandomSamplingMean(data1[,grep("NCP",colnames(data1))]))
VirtualMatrix<-VirtualMatrix[-which(! is.finite(rowSums(VirtualMatrix))),]
head(VirtualMatrix)
# Deconvolution-QR
DeconData<-data.frame(VirtualMatrix,signatures[match(rownames(VirtualMatrix),rownames(signatures)),])
DeconData.tmp<-DeconData[na.omit(match(topgsi[,1],rownames(DeconData))),]
VirtualMatrix=data.frame(DeconData.tmp[,grep("VirtualSample",colnames(DeconData))])
Signatures=data.frame(DeconData.tmp[,-grep("VirtualSample",colnames(DeconData))])
library(car)
NewVirtualMatrix=logit(VirtualMatrix)
NewSignatures=logit(Signatures)
head(NewVirtualMatrix)
head(NewSignatures)
Rlt<-DeconRNASeq(NewVirtualMatrix,NewSignatures,checksig=FALSE,known.prop = F, use.scale = TRUE, fig = TRUE)
output=rbind(Rlt$out.all,colMeans(Rlt$out.all),ci95(Rlt$out.all))
write.table(output,file=paste(outputprefix,".deconvolution.result.txt",sep=""),sep="\t",quote=F,col.names=NA,row.names=T)
save.image(file=paste(outputprefix,"workspace.image.RData",sep="."))

